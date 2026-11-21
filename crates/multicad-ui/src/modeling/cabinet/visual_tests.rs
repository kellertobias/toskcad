//! CAD-389 observable behaviour of the visual cabinet assistant.
use super::*;
use crate::modeling::diagram::Mark;
use multicad_core::{
    commands::{Actor, Channel},
    document::{Entity, EntityData},
    materials::{Grain, Stock, prepare_library},
};
use multicad_model::Compiler;

fn session() -> Session {
    let mut s = Session::new(Document::new("Visual cabinet").unwrap()).unwrap();
    let birch = Entity::new(
        "Birch",
        EntityData::Material {
            density_kg_m3: None,
            color_rgba: [210, 180, 140, 255],
            section_hatch: None,
        },
    );
    let stock = |t: f64, name: &str| {
        LibraryDefinition::Stock(Stock {
            material: birch.id,
            specification: StockKind::Panel {
                thickness_mm: t,
                length_mm: 2500.,
                width_mm: 1250.,
                grain: Grain::AlongLength,
            },
        })
        .entity(EntityId::new(), name)
        .unwrap()
    };
    let entities = [birch.clone(), stock(18., "Ply 18"), stock(6., "Ply 6")];
    s.commit(prepare_library(&s, 0, Actor::local(Channel::Ui), &entities).unwrap())
        .unwrap();
    s
}

struct Painted {
    texts: Vec<(String, egui::Rect)>,
    event: Option<ModelEvent>,
}
fn frame(
    ctx: &egui::Context,
    form: &mut CabinetPanel,
    s: &Session,
    width: f32,
    events: Vec<egui::Event>,
) -> Painted {
    let mut event = None;
    let output = ctx.run_ui(
        egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(
                egui::Pos2::ZERO,
                egui::vec2(width, 4000.),
            )),
            focused: true,
            events,
            ..Default::default()
        },
        |ui| {
            event = form.show(ui, s.document(), true, PreviewState::Idle);
        },
    );
    let texts = output
        .shapes
        .iter()
        .filter_map(|shape| match &shape.shape {
            egui::Shape::Text(t) => Some((
                t.galley.text().to_owned(),
                egui::Rect::from_min_size(t.pos, t.galley.size()),
            )),
            _ => None,
        })
        .collect();
    output.drop_without_applying_deltas();
    Painted { texts, event }
}
fn settle(ctx: &egui::Context, form: &mut CabinetPanel, s: &Session, width: f32) -> Painted {
    for _ in 0..10 {
        frame(ctx, form, s, width, vec![]);
    }
    frame(ctx, form, s, width, vec![])
}
fn find(p: &Painted, prefix: &str) -> Option<egui::Rect> {
    p.texts
        .iter()
        .find(|(t, _)| t.starts_with(prefix))
        .map(|(_, r)| *r)
}
fn click(pos: egui::Pos2) -> Vec<egui::Event> {
    [true, false]
        .into_iter()
        .flat_map(|pressed| {
            [
                egui::Event::PointerMoved(pos),
                egui::Event::PointerButton {
                    pos,
                    button: egui::PointerButton::Primary,
                    pressed,
                    modifiers: Default::default(),
                },
            ]
        })
        .collect()
}
fn replace_text(text: &str) -> Vec<egui::Event> {
    vec![
        egui::Event::Key {
            key: egui::Key::A,
            physical_key: None,
            pressed: true,
            repeat: false,
            modifiers: egui::Modifiers::COMMAND,
        },
        egui::Event::Text(text.into()),
    ]
}
fn carcass_width(form: &CabinetPanel, s: &Session) -> Dimension {
    let CabinetCommand::Carcass { request, .. } = form.request(s.document(), true).unwrap() else {
        unreachable!()
    };
    request.definition.width
}

#[test]
fn construction_and_back_choices_draw_what_they_produce() {
    let s = session();
    let mut form = CabinetPanel::cabinet(s.document(), "Choices".into());
    let stocks = panel_stocks(s.document());
    let widest = |d: &diagram::Diagram| {
        d.marks
            .iter()
            .filter_map(|m| match m {
                Mark::Path {
                    points, fill: true, ..
                } => {
                    let xs = points.iter().map(|p| p[0]);
                    let ys = points.iter().map(|p| p[1]);
                    Some((
                        xs.clone().fold(f32::MIN, f32::max) - xs.fold(f32::MAX, f32::min),
                        ys.clone().fold(f32::MIN, f32::max) - ys.fold(f32::MAX, f32::min),
                    ))
                }
                _ => None,
            })
            .collect::<Vec<_>>()
    };
    form.construction = Construction::SidesFullHeight;
    let sides = form.schematic(s.document(), &stocks);
    form.construction = Construction::TopBottomFullWidth;
    let tops = form.schematic(s.document(), &stocks);
    assert_ne!(sides, tops);
    // Full-width top: some panel spans the whole 600 mm width.
    assert!(widest(&tops).iter().any(|(w, _)| (*w - 600.).abs() < 1e-3));
    assert!(!widest(&sides).iter().any(|(w, _)| (*w - 600.).abs() < 1e-3));
    assert_ne!(
        visual::construction_art(Construction::SidesFullHeight),
        visual::construction_art(Construction::TopBottomFullWidth)
    );
    let keys = |form: &CabinetPanel| form.schematic(s.document(), &stocks).keys();
    form.back = BackKind::Inset;
    assert!(keys(&form).contains(&visual::RECESS));
    form.back = BackKind::Overlay;
    assert!(!keys(&form).contains(&visual::RECESS));
    let overlay = form.schematic(s.document(), &stocks);
    form.back = BackKind::None;
    assert_ne!(form.schematic(s.document(), &stocks), overlay);
    let arts = [BackKind::None, BackKind::Inset, BackKind::Overlay].map(visual::back_art);
    assert!(arts[0] != arts[1] && arts[1] != arts[2] && arts[0] != arts[2]);
    // Every shelf row is drawn and dimensioned by its own retained slot.
    let shelf = visual::shelf_keys(0).0;
    assert!(keys(&form).contains(&shelf));
    for key in keys(&form) {
        assert!(
            form.linked_fields().iter().any(|(k, _)| *k == key),
            "{key} has a field"
        );
    }
}

#[test]
fn essentials_first_advanced_groups_collapsed_with_retained_values() {
    let s = session();
    for width in [200., 320.] {
        let ctx = egui::Context::default();
        let mut form = CabinetPanel::cabinet(s.document(), "Essentials".into());
        let request = form.request(s.document(), true).unwrap();
        let painted = settle(&ctx, &mut form, &s, width);
        for essential in [
            "Width (mm)",
            "Floor-to-top height (mm)",
            "Depth (mm)",
            "Toe-kick clearance (mm)",
            "Full-height sides",
            "Inset back",
            "Add shelf",
        ] {
            let rect = find(&painted, essential)
                .unwrap_or_else(|| panic!("{essential} visible at {width}"));
            assert!(
                rect.left() >= 0. && rect.right() <= width + 0.5,
                "{essential} unclipped at {width}"
            );
        }
        for advanced in ["Edge strip thickness (mm)", "Origin X (mm)", "Shelf setback (mm)"] {
            assert!(
                find(&painted, advanced).is_none(),
                "{advanced} starts collapsed"
            );
        }
        assert!(find(&painted, "Edge banding · front edges 2 mm").is_some());
        assert!(find(&painted, "Shelf setback and end clearance · 20 mm · 1 mm").is_some());
        for (text, rect) in &painted.texts {
            if text.starts_with("Width ")
                || text.starts_with("Depth ")
                || text.starts_with("Floor to top")
            {
                assert!(
                    rect.left() >= 0. && rect.right() <= width + 0.5,
                    "callout {text} unclipped at {width}"
                );
            }
        }
        let CabinetCommand::Carcass { request: after, .. } =
            form.request(s.document(), true).unwrap()
        else {
            unreachable!()
        };
        let CabinetCommand::Carcass {
            request: before, ..
        } = request
        else {
            unreachable!()
        };
        assert_eq!(
            after.definition, before.definition,
            "painting never changes the definition"
        );
    }
}

#[test]
fn callouts_edit_the_same_fields_and_reveal_collapsed_groups() {
    let s = session();
    let ctx = egui::Context::default();
    let mut form = CabinetPanel::cabinet(s.document(), "Callouts".into());
    let painted = settle(&ctx, &mut form, &s, 320.);
    let width = find(&painted, "Width 600 mm").expect("width callout");
    frame(&ctx, &mut form, &s, 320., click(width.center()));
    frame(&ctx, &mut form, &s, 320., vec![]);
    let typed = frame(&ctx, &mut form, &s, 320., replace_text("800"));
    assert!(
        matches!(typed.event, Some(ModelEvent::Changed)),
        "edits restart the automatic preview"
    );
    assert_eq!(carcass_width(&form, &s), length(800.));
    let painted = settle(&ctx, &mut form, &s, 320.);
    assert!(
        find(&painted, "Width 800 mm").is_some(),
        "diagram follows the staged value"
    );
    // The setback lives in a collapsed group; its callout reveals and focuses it.
    let setback = find(&painted, "Setback 20 mm").expect("setback callout");
    frame(&ctx, &mut form, &s, 320., click(setback.center()));
    let painted = settle(&ctx, &mut form, &s, 320.);
    assert!(
        find(&painted, "Shelf setback (mm)").is_some(),
        "callout reveals its group"
    );
    frame(&ctx, &mut form, &s, 320., replace_text("35"));
    let CabinetCommand::Carcass { request, .. } = form.request(s.document(), true).unwrap() else {
        unreachable!()
    };
    assert_eq!(request.definition.shelves.setback, length(35.));
    // A shelf height callout edits that shelf's own row.
    let painted = settle(&ctx, &mut form, &s, 320.);
    let shelf = find(&painted, "Shelf 1 underside 360 mm").expect("shelf callout");
    frame(&ctx, &mut form, &s, 320., click(shelf.center()));
    frame(&ctx, &mut form, &s, 320., vec![]);
    frame(&ctx, &mut form, &s, 320., replace_text("300"));
    let options = form.options.as_ref().unwrap().spec(s.document()).unwrap();
    assert_eq!(options.shelves[0].underside, length(300.));
}

#[test]
fn validation_identifies_and_reveals_the_affected_group() {
    let s = session();
    let ctx = egui::Context::default();
    let mut form = CabinetPanel::cabinet(s.document(), "Invalid".into());
    form.band.text = "=".into();
    let painted = settle(&ctx, &mut form, &s, 320.);
    assert_eq!(form.link.invalid, Some(visual::BAND[0]));
    assert!(find(&painted, "Edge strip thickness (mm)").is_some());
    assert!(form.request(s.document(), true).is_err());
    form.band.text = "2".into();
    form.origin[1].text = "=".into();
    let painted = settle(&ctx, &mut form, &s, 320.);
    assert_eq!(form.link.invalid, Some(visual::PLACEMENT[1]));
    assert!(find(&painted, "Origin Y (mm)").is_some());
    form.origin[1].text = "0".into();
    assert!(form.request(s.document(), true).is_ok());
    settle(&ctx, &mut form, &s, 320.);
    assert_eq!(form.link.invalid, None);
}

#[test]
fn editing_a_representative_cabinet_reopens_with_its_diagram_and_applies_once() {
    let mut s = session();
    let form = CabinetPanel::cabinet(s.document(), "Edited".into());
    let CabinetCommand::Carcass { request, .. } = form.request(s.document(), true).unwrap() else {
        unreachable!()
    };
    let prepared =
        model::prepare_carcass(&s, s.revision(), Actor::local(Channel::Ui), &request).unwrap();
    Compiler::default()
        .compile(prepared.model.prepared.document(), 1)
        .unwrap();
    s.commit(prepared.model.prepared).unwrap();
    let group = model::read_carcass(s.document(), request.carcass)
        .unwrap()
        .group
        .unwrap();
    let target = cabinet_edit::carcass_target(s.document(), group).unwrap();
    let mut form = CabinetPanel::existing(s.document(), target.feature).unwrap();
    for width in [200., 320.] {
        let ctx = egui::Context::default();
        let painted = settle(&ctx, &mut form, &s, width);
        assert!(
            find(&painted, "Width 600 mm").is_some(),
            "reopened diagram at {width}"
        );
        assert!(find(&painted, "Inset back").is_some());
    }
    let ctx = egui::Context::default();
    let painted = settle(&ctx, &mut form, &s, 320.);
    let overlay = find(&painted, "Overlay back").unwrap();
    frame(&ctx, &mut form, &s, 320., click(overlay.center()));
    assert_eq!(form.back, BackKind::Overlay);
    let prepared = form
        .request(s.document(), true)
        .unwrap()
        .prepare(&s, Actor::local(Channel::Ui))
        .unwrap();
    Compiler::default()
        .compile(prepared.prepared.document(), s.revision() + 1)
        .unwrap();
    let before = s.revision();
    s.commit(prepared.prepared).unwrap();
    assert_eq!(s.revision(), before + 1);
    assert!(matches!(
        model::read_carcass(s.document(), request.carcass)
            .unwrap()
            .back,
        Back::Overlay { .. }
    ));
}
