//! CAD-430: the panel inspector edits position (X/Y/Z with units) and
//! orientation (face direction preset and rotation). Each edit runs through
//! the full application with real pointer and key events and applies as
//! one canonical transaction.
use super::ui_harness::Harness;
use super::*;
use multicad_core::{
    materials::{Grain, LibraryDefinition, Stock, StockKind, prepare_library},
    model::{self, PanelFacing},
};
use std::time::{Duration, Instant};

fn ui() -> Actor {
    Actor::local(Channel::Ui)
}
/// One 600 × 400 panel of 18 mm stock on XY at the origin.
fn panel_document() -> (Document, EntityId, EntityId) {
    let mut s = Session::new(Document::new("Panel placement").unwrap()).unwrap();
    let birch = Entity::new(
        "Birch",
        EntityData::Material {
            density_kg_m3: None,
            color_rgba: [210, 180, 140, 255],
            section_hatch: None,
        },
    );
    let stock = LibraryDefinition::Stock(Stock {
        material: birch.id,
        specification: StockKind::Panel {
            thickness_mm: 18.,
            length_mm: 2500.,
            width_mm: 1250.,
            grain: Grain::AlongLength,
        },
    })
    .entity(EntityId::new(), "Ply 18")
    .unwrap();
    let stock_id = stock.id;
    s.commit(prepare_library(&s, 0, ui(), &[birch, stock]).unwrap())
        .unwrap();
    let request = model::PanelRequest::new(
        "Panel 1",
        model::PanelDefinition {
            output: model::BodyOutput {
                body: EntityId::new(),
            },
            source: model::PanelSource::Rectangle {
                stock: stock_id,
                plane: model::PanelPlane::Xy,
                origin: [0., 0., 0.].map(model::length),
                length: model::length(600.),
                width: model::length(400.),
            },
            grain: model::GrainAxis::U,
            edges: Default::default(),
            oversize: None,
        },
    );
    let p = model::prepare_panel(&s, s.revision(), ui(), &request).unwrap();
    s.commit(p.prepared).unwrap();
    (s.document().clone(), request.feature, request.body)
}
fn open(h: &mut Harness, feature: EntityId) {
    h.until("the current model", |app| {
        let p = app.current();
        p.scene
            .is_current(p.session.document().id(), p.session.revision())
    });
    let p = h.app.current_mut();
    p.sketch
        .begin_feature_edit(&p.session, &p.resources, &p.scene, feature, Channel::Ui)
        .unwrap();
    let deadline = Instant::now() + Duration::from_secs(20);
    while !h.texts().iter().any(|t| t == "Position X") {
        assert!(
            Instant::now() < deadline,
            "timed out waiting for the panel inspector"
        );
    }
}
fn placement(h: &Harness, part: EntityId) -> model::PanelPlacement {
    let p = h.app.current();
    model::panel_placement(p.session.document(), p.session.parameters(), part).unwrap()
}
fn history(h: &Harness, part: EntityId) -> usize {
    match &h.app.current().session.document().get(part).unwrap().data {
        EntityData::Part { features, .. } => features.len(),
        _ => panic!("Part expected"),
    }
}
fn applied(h: &mut Harness, what: &str, revision: u64) {
    h.until(what, |app| app.current().session.revision() == revision);
    assert!(!h.modeling(), "{what} closes the inspector edit");
}

#[test]
fn the_panel_inspector_edits_position_face_direction_and_rotation_one_step_each() {
    let (document, feature, part) = panel_document();
    let mut h = Harness::new(document);
    open(&mut h, feature);
    let texts = h.texts();
    for shown in [
        "Position X",
        "Position Y",
        "Position Z",
        "Face direction",
        "+Z · face up",
        "Rotation about face",
    ] {
        assert!(texts.iter().any(|t| t == shown), "{shown} in {texts:?}");
    }

    // Position: typed value, Enter applies one transaction in the feature.
    let start = h.revision();
    h.type_into_field_and_enter("Position X", "150");
    applied(&mut h, "the position edit", start + 1);
    let read = placement(&h, part);
    assert_eq!(read.origin_mm, [150., 0., 0.]);
    assert_eq!(read.facing, Some(PanelFacing::PosZ));
    assert_eq!(history(&h, part), 1);

    // Orientation preset: open the combo, pick face down, Enter applies.
    open(&mut h, feature);
    h.click_label("+Z · face up");
    let deadline = Instant::now() + Duration::from_secs(20);
    while !h.texts().iter().any(|t| t == "−Z · face down") {
        assert!(
            Instant::now() < deadline,
            "timed out waiting for the presets"
        );
    }
    h.click_label("−Z · face down");
    h.until("the face-down preview", |app| {
        let p = app.current();
        p.sketch.model_can_apply(&p.session, &p.resources)
    });
    h.key(egui::Key::Enter);
    applied(&mut h, "the preset edit", start + 2);
    let read = placement(&h, part);
    assert_eq!(read.facing, Some(PanelFacing::NegZ));
    assert_eq!(read.origin_mm, [150., 0., 0.]);
    // Carried by one trailing Move / rotate.
    assert_eq!(history(&h, part), 2);

    // Rotation about the face: the same Move / rotate changes in place.
    open(&mut h, feature);
    h.type_into_field_and_enter("Rotation about face", "90");
    applied(&mut h, "the rotation edit", start + 3);
    let read = placement(&h, part);
    assert_eq!(read.facing, Some(PanelFacing::NegZ));
    assert!((read.rotation_deg - 90.).abs() < 1e-9, "{read:?}");
    assert_eq!(read.origin_mm, [150., 0., 0.]);
    assert_eq!(history(&h, part), 2);

    // One Undo reverts only the rotation.
    let p = h.app.current_mut();
    let revision = p.session.revision();
    p.session.undo(revision, ui()).unwrap();
    let read = placement(&h, part);
    assert_eq!(read.facing, Some(PanelFacing::NegZ));
    assert!(read.rotation_deg.abs() < 1e-9, "{read:?}");

    // Reopened, the inspector shows the placement made so far.
    open(&mut h, feature);
    let texts = h.texts();
    assert!(texts.iter().any(|t| t == "−Z · face down"), "{texts:?}");
    assert!(texts.iter().any(|t| t == "150"), "{texts:?}");
}
