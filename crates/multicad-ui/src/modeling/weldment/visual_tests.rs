//! CAD-395 weldment schematics: the drawn section is the chosen stock's own
//! section placed by the member anchor and rotation; callouts edit the
//! member's staged fields and named sections open on demand or validation.
use super::*;
use crate::modeling::diagram::harness::{self, WIDTHS};
use multicad_core::{
    commands::Channel,
    document::Entity,
    materials::{ProfileShape, Stock, prepare_library},
    model::{ProfileSection, SectionLoop},
};

fn profile(name: &str, w: f64, h: f64, wall: Option<f64>, shape: ProfileShape) -> Entity {
    LibraryDefinition::Stock(Stock {
        material: EntityId::new(),
        specification: StockKind::Profile {
            designation: name.into(),
            length_mm: 6000.,
            width_mm: w,
            height_mm: h,
            wall_mm: wall,
            shape,
            section: None,
            corner_radius_mm: None,
            web_mm: None,
        },
    })
    .entity(EntityId::new(), name)
    .unwrap()
}

/// SHS 30 × 2, L 40 × 30 × 4 and CHS 30 × 2 on one steel material.
fn session() -> (Session, [EntityId; 3]) {
    let mut s = Session::new(Document::new("Weldment diagrams").unwrap()).unwrap();
    let steel = Entity::new(
        "Steel",
        EntityData::Material {
            density_kg_m3: None,
            color_rgba: [120, 124, 130, 255],
            section_hatch: None,
        },
    );
    let mut stocks = [
        profile("SHS 30", 30., 30., Some(2.), ProfileShape::Rectangular),
        profile("L 40x30", 40., 30., Some(4.), ProfileShape::Angle),
        profile("CHS 30", 30., 30., Some(2.), ProfileShape::Round),
    ];
    for stock in &mut stocks {
        let Some(LibraryDefinition::Stock(mut definition)) = materials::read(stock).unwrap() else {
            panic!("stock")
        };
        definition.material = steel.id;
        *stock = LibraryDefinition::Stock(definition)
            .entity(stock.id, &stock.name)
            .unwrap();
    }
    let ids = [stocks[0].id, stocks[1].id, stocks[2].id];
    let mut entities = vec![steel];
    entities.extend(stocks);
    s.commit(prepare_library(&s, 0, Actor::local(Channel::Ui), &entities).unwrap())
        .unwrap();
    (s, ids)
}

fn panel(s: &Session, stock: EntityId) -> WeldmentPanel {
    let mut p = WeldmentPanel::new(s.document(), "Rail".into());
    p.stock = Some(stock);
    p.path_kind = PathKind::Points;
    p
}

fn show<'a>(p: &'a mut WeldmentPanel, d: &'a Document) -> impl FnMut(&mut egui::Ui) + 'a {
    move |ui| {
        let _ = p.show(ui, d, true, PreviewState::Idle);
    }
}

fn member(p: &WeldmentPanel, d: &Document) -> MemberDefinition {
    let WeldmentCommand::Members { requests, .. } = p.request(d, true).unwrap() else {
        panic!("members")
    };
    requests[0].member.clone()
}

#[test]
fn the_drawn_section_is_the_stocks_own_section_placed_by_anchor_and_rotation() {
    let (s, [shs, angle_stock, chs]) = session();
    let d = s.document();
    for stock in [shs, angle_stock, chs] {
        let p = panel(&s, stock);
        let section = p.section(d).expect("stock section");
        let expected =
            ProfileSection::from_stock(&match materials::read(d.get(stock).unwrap()).unwrap() {
                Some(LibraryDefinition::Stock(s)) => s,
                _ => panic!("stock"),
            })
            .unwrap();
        assert_eq!(section, expected);
        let loops = visual::placed_loops(&section, [AxisAnchor::Min; 2], 0.);
        assert_eq!(loops.len(), 1 + section.holes().len(), "outer plus holes");
        if let SectionLoop::Polygon(outer) = section.outer(0.) {
            let drawn: Vec<[f32; 2]> = outer.iter().map(|p| [p[0] as f32, p[1] as f32]).collect();
            assert_eq!(
                loops[0], drawn,
                "a Min/Min anchor draws the section as stored"
            );
        }
        // A centred anchor moves the section around the path; rotation turns it.
        let centred = visual::placed_loops(&section, [AxisAnchor::Center; 2], 0.);
        let turned = visual::placed_loops(&section, [AxisAnchor::Center; 2], 90.);
        assert_ne!(loops, centred);
        assert_ne!(centred, turned);
        let schematic = visual::diagram(&p.diagram_values(d));
        let fields: Vec<_> = p.linked_fields().into_iter().map(|(k, _)| k).collect();
        for key in schematic.keys() {
            assert!(
                key == visual::ANCHOR || fields.contains(&key),
                "{key} has a staged field"
            );
        }
    }
    let p = panel(&s, angle_stock);
    let first = visual::placed_loops(&p.section(d).unwrap(), [AxisAnchor::Min; 2], 0.);
    assert_eq!(first[0].len(), 6, "an angle is its L outline");
    let p = panel(&s, chs);
    assert_eq!(
        visual::placed_loops(&p.section(d).unwrap(), [AxisAnchor::Min; 2], 0.).len(),
        2,
        "a round tube has its bore"
    );
}

#[test]
fn clicking_the_rotation_callout_edits_the_member_rotation_at_both_widths() {
    let (s, [shs, ..]) = session();
    let d = s.document();
    for width in WIDTHS {
        let ctx = egui::Context::default();
        let mut p = panel(&s, shs);
        let painted =
            harness::edit_through_callout(&ctx, width, "∠ Rotation 0°", "45", &mut show(&mut p, d));
        assert!(painted.has("∠ Rotation 45°"));
        // Every callout keeps its own target at this width.
        let schematic = visual::diagram(&p.diagram_values(d));
        harness::assert_callouts_apart(&schematic, width, 260., &|k| p.callout_value(k));
        let m = member(&p, d);
        assert_eq!(m.rotation, Some(angle(45.)));
        assert_eq!(m.extend_start, None, "only the rotation changed");
    }
}

#[test]
fn extensions_and_points_are_collapsed_retained_and_revealed() {
    let (s, [shs, ..]) = session();
    let d = s.document();
    let ctx = egui::Context::default();
    let mut p = panel(&s, shs);
    p.extend_end.text = "25".into();
    let painted = harness::settle(&ctx, 360., &mut show(&mut p, d));
    assert!(
        !painted.has("Extend start (mm)"),
        "extensions start collapsed"
    );
    assert!(painted.has_prefix("End extensions · 0 mm / 25 mm"));
    assert!(!painted.has("End X (mm)"), "path points start collapsed");
    assert_eq!(member(&p, d).extend_end, Some(length(25.)));
    p.extend_start.text = "=".into();
    let painted = harness::settle(&ctx, 360., &mut show(&mut p, d));
    assert_eq!(p.link.invalid, Some(visual::EXTEND_START));
    assert!(painted.has("Extend start (mm)"), "validation reveals it");
    p.extend_start.text = "0".into();
    let end = painted.find("End 1000, 0, 0 mm").expect("end callout");
    harness::frame(
        &ctx,
        360.,
        harness::click(end.center()),
        &mut show(&mut p, d),
    );
    let painted = harness::settle(&ctx, 360., &mut show(&mut p, d));
    assert!(
        painted.has("End X (mm)"),
        "the End callout opens the points"
    );
    assert_eq!(p.link.active(), Some(visual::END[0]));
}
