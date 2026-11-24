//! CAD-395 base-sheet schematic: stock thickness, inside radius and
//! K-factor are explained on a representative bend; callouts edit the
//! form's staged size and override fields.
use super::*;
use crate::modeling::diagram::{
    Key, Mark,
    harness::{self, WIDTHS},
};
use multicad_core::{
    commands::Channel,
    document::Entity,
    materials::{Stock, prepare_library},
    model::length,
};

/// DC01 2.0: 2 mm sheet, 3 mm inside radius, K 0.4.
fn session() -> Session {
    let mut s = Session::new(Document::new("Sheet diagrams").unwrap()).unwrap();
    let steel = Entity::new(
        "DC01",
        EntityData::Material {
            density_kg_m3: None,
            color_rgba: [180, 184, 190, 255],
            section_hatch: None,
        },
    );
    let stock = LibraryDefinition::Stock(Stock {
        material: steel.id,
        specification: StockKind::Sheet {
            thickness_mm: 2.,
            length_mm: 2000.,
            width_mm: 1000.,
            inside_bend_radius_mm: 3.,
            k_factor: 0.4,
        },
    })
    .entity(EntityId::new(), "DC01 2.0")
    .unwrap();
    s.commit(prepare_library(&s, 0, Actor::local(Channel::Ui), &[steel, stock]).unwrap())
        .unwrap();
    s
}

fn show<'a>(p: &'a mut SheetPanel, d: &'a Document) -> impl FnMut(&mut egui::Ui) + 'a {
    move |ui| {
        let _ = p.show(ui, d, true, PreviewState::Idle);
    }
}

#[test]
fn the_bend_section_reads_the_stock_and_the_overrides_it_explains() {
    let s = session();
    let d = s.document();
    let mut p = SheetPanel::new(d, "Cover".into());
    let stock = visual::diagram(&p.diagram_values(d));
    assert!(
        stock
            .marks
            .iter()
            .any(|m| matches!(m, Mark::Note { text, .. } if text == "t 2 mm · stock")),
        "thickness is the stock's, as a note"
    );
    let keys = stock.keys();
    for key in [
        visual::RADIUS_OVERRIDE,
        visual::K_OVERRIDE,
        visual::LENGTH,
        visual::WIDTH,
    ] {
        assert!(keys.contains(&key), "{key}");
    }
    assert_eq!(p.callout_value(d, visual::RADIUS_OVERRIDE), "3 mm");
    assert_eq!(p.callout_value(d, visual::K_OVERRIDE), "0.4");
    p.radius_override = true;
    p.k_override = true;
    p.k_factor.text = "0.33".into();
    let overridden = visual::diagram(&p.diagram_values(d));
    let fields: Vec<Key> = p.linked_fields().into_iter().map(|(k, _)| k).collect();
    for key in overridden.keys() {
        assert!(fields.contains(&key), "{key} edits a staged field");
    }
    assert_ne!(overridden, stock, "the neutral axis follows the K-factor");
    p.reverse = true;
    assert_ne!(visual::diagram(&p.diagram_values(d)), overridden);
}

#[test]
fn clicking_the_length_callout_edits_the_base_length_at_both_widths() {
    let s = session();
    let d = s.document();
    for width in WIDTHS {
        let ctx = egui::Context::default();
        let mut p = SheetPanel::new(d, "Cover".into());
        let painted = harness::edit_through_callout(
            &ctx,
            width,
            "Length (u) 200 mm",
            "250",
            &mut show(&mut p, d),
        );
        assert!(painted.has("Length (u) 250 mm"));
        let sheet = p.definition(d).unwrap();
        let SheetBase::Rectangle {
            length: l,
            width: w,
            ..
        } = sheet.base
        else {
            panic!("rectangle")
        };
        assert_eq!((l, w), (length(250.), length(100.)));
        for radius_override in [false, true] {
            p.radius_override = radius_override;
            p.k_override = radius_override;
            let schematic = visual::diagram(&p.diagram_values(d));
            harness::assert_callouts_apart(&schematic, width, 220., &|k| p.callout_value(d, k));
        }
    }
}

#[test]
fn bend_overrides_and_position_are_collapsed_retained_and_revealed() {
    let s = session();
    let d = s.document();
    let ctx = egui::Context::default();
    let mut p = SheetPanel::new(d, "Cover".into());
    let painted = harness::settle(&ctx, 360., &mut show(&mut p, d));
    assert!(
        !painted.has("Override inside radius"),
        "overrides start collapsed"
    );
    assert!(painted.has_prefix("Bend overrides · R stock 3 mm · K stock 0.4"));
    assert!(!painted.has("Origin X (mm)"), "position starts collapsed");
    // The stock radius callout opens the overrides at its checkbox.
    let callout = painted.find("R Stock radius 3 mm").expect("radius callout");
    harness::frame(
        &ctx,
        360.,
        harness::click(callout.center()),
        &mut show(&mut p, d),
    );
    let painted = harness::settle(&ctx, 360., &mut show(&mut p, d));
    assert!(painted.has("Override inside radius"));
    assert_eq!(p.link.active(), Some(visual::RADIUS_OVERRIDE));

    // A collapsed override is still in the definition; validation reveals it.
    let ctx = egui::Context::default();
    let mut p = SheetPanel::new(d, "Cover".into());
    p.k_override = true;
    p.k_factor.text = "0.45".into();
    p.origin[1].text = "20".into();
    let painted = harness::settle(&ctx, 360., &mut show(&mut p, d));
    assert!(!painted.has("K-factor (ratio)"));
    assert!(painted.has_prefix("Base position · 0, 20, 0 mm"));
    let sheet = p.definition(d).unwrap();
    assert_eq!(sheet.k_factor, p.k_factor.resolve(d).ok());
    assert!(sheet.k_factor.is_some());
    p.k_factor.text = "=".into();
    let painted = harness::settle(&ctx, 360., &mut show(&mut p, d));
    assert_eq!(p.link.invalid, Some(visual::K_FACTOR));
    assert!(
        painted.has("K-factor (ratio)"),
        "validation reveals the override"
    );
    p.k_factor.text = "0.45".into();
    p.origin[2].text = "=".into();
    let painted = harness::settle(&ctx, 360., &mut show(&mut p, d));
    assert_eq!(p.link.invalid, Some(visual::ORIGIN[2]));
    assert!(painted.has("Origin Z (mm)"));
}
