//! CAD-410: Enter applies the active modelling tool exactly as Apply does,
//! Escape cancels it, and an invalid value in a focused field is not applied.
//! Every step runs through the full application key handling.
use super::ui_harness::Harness;
use super::*;
use crate::modeling::ModelingTool;
use multicad_core::{
    constraint_sketch::source::Source,
    model::{self, PrimitiveParameters, PrimitiveRequest},
    references::{GeometryIntent, TopologyKind},
    sketch::Datum,
};
use std::{
    collections::BTreeSet,
    time::{Duration, Instant},
};

fn box_document() -> (Document, model::BodyInput) {
    let mut session = Session::new(Document::new("Enter applies Fillet").unwrap()).unwrap();
    let request = PrimitiveRequest::new(
        "Stock",
        PrimitiveParameters::Box {
            size: [model::length(30.), model::length(20.), model::length(10.)],
            frame: Default::default(),
        },
    );
    session
        .commit(
            model::prepare_primitive(&session, 0, Actor::local(Channel::Ui), &request)
                .unwrap()
                .prepared,
        )
        .unwrap();
    let body = model::BodyInput {
        body: request.body,
        producer: request.feature,
    };
    (session.document().clone(), body)
}
fn rectangle_sketch() -> (Document, EntityId) {
    let (sketch, edits) = Source::create(Datum::Xy, "Rectangle").unwrap();
    let document = Document::new("Enter applies Extrude")
        .unwrap()
        .edited(&edits)
        .unwrap();
    let (mut source, _, _) = Source::read(&document, sketch).unwrap();
    source.add_rectangle([0., 0.], [30., 20.]).unwrap();
    let document = document
        .edited(&source.replace(&document, sketch).unwrap())
        .unwrap();
    (document, sketch)
}
fn last_extrusion(document: &Document) -> model::ExtrudeDefinition {
    document
        .entities()
        .filter_map(|entity| model::read_extrude(document, entity.id).ok())
        .last()
        .expect("an applied extrusion")
}
fn distance_mm(definition: &model::ExtrudeDefinition) -> f64 {
    let model::Extent::OneSided {
        end: model::Termination::Distance { distance },
    } = &definition.parameters.extent
    else {
        panic!("one-sided distance extrusion expected")
    };
    match distance {
        multicad_core::constraint_sketch::Dimension::Literal { value, unit } if unit == "mm" => {
            *value
        }
        other => panic!("literal millimetres expected, got {other:?}"),
    }
}

#[test]
fn enter_applies_fillet_with_picked_edges_like_apply_and_escape_cancels() {
    let (document, body) = box_document();
    let mut h = Harness::new(document);
    let open = |h: &mut Harness| {
        h.app
            .current_mut()
            .sketch
            .set_tree_selection_for_tests(BTreeSet::from([body.body]));
        h.app
            .invoke(
                Command::Modeling(ModelingTool::EdgeTreatment),
                &[],
                Channel::Ui,
            )
            .unwrap();
        h.frame(vec![]);
        assert!(
            h.app
                .current()
                .sketch
                .model_for_tests()
                .edge_treatment_active()
        );
    };
    open(&mut h);
    // Picking is live: there is no separate control to finish it.
    assert!(
        h.app
            .current()
            .sketch
            .model_for_tests()
            .edge_treatment_picking()
    );
    let texts = h.texts();
    assert!(
        !texts
            .iter()
            .any(|t| t.contains("Finish picking") || t.contains("reselect edges")),
        "{texts:?}"
    );
    // Enter with nothing picked does not apply.
    h.key(egui::Key::Enter);
    assert_eq!(h.revision(), 0);
    assert!(!h.modeling());
    open(&mut h);
    // Pick one edge; the preview builds while picking stays active.
    let p = h.app.current_mut();
    let scene = p.scene.visible().unwrap().clone();
    let topology = scene.geometry_snapshot(body.producer).unwrap();
    let edge = topology
        .elements
        .iter()
        .find(|e| e.key.kind == TopologyKind::Edge)
        .unwrap();
    let capture = GeometryIntent::capture(&topology, edge.key).unwrap();
    p.sketch
        .model_for_tests_mut()
        .set_edge_treatment_edges(&p.session, &p.resources, &scene, vec![capture])
        .unwrap();
    h.until("the fillet preview", |app| {
        let p = app.current();
        p.sketch.model_can_apply(&p.session, &p.resources)
    });
    assert!(
        h.app
            .current()
            .sketch
            .model_for_tests()
            .edge_treatment_picking()
    );
    // Enter commits exactly what Apply would.
    h.key(egui::Key::Enter);
    h.until("Enter to apply the fillet", |app| {
        app.current().session.revision() == 1
    });
    assert!(!h.modeling());
    let document = h.app.current().session.document().clone();
    let EntityData::Part { features, .. } = &document.get(body.body).unwrap().data else {
        panic!("Part expected")
    };
    let (saved, _) = model::read_edge_treatment(&document, *features.last().unwrap()).unwrap();
    assert_eq!(
        saved.parameters.profile,
        model::EdgeTreatmentProfile::Fillet
    );
    // Escape cancels a reopened tool without touching the document.
    h.until("the applied model", |app| {
        let p = app.current();
        p.scene
            .is_current(p.session.document().id(), p.session.revision())
    });
    open(&mut h);
    h.key(egui::Key::Escape);
    assert!(!h.modeling());
    assert_eq!(h.app.current().session.document(), &document);
    assert_eq!(h.revision(), 1);
}

#[test]
fn enter_applies_extrude_and_an_invalid_focused_value_shows_its_error_instead() {
    let (document, sketch) = rectangle_sketch();
    let mut h = Harness::new(document);
    let open = |h: &mut Harness| {
        h.app
            .current_mut()
            .sketch
            .set_tree_selection_for_tests(BTreeSet::from([sketch]));
        h.app
            .invoke(Command::Modeling(ModelingTool::Extrude), &[], Channel::Ui)
            .unwrap();
        h.until("the extrude preview", |app| {
            let p = app.current();
            p.sketch.model_can_apply(&p.session, &p.resources)
        });
    };
    open(&mut h);
    // Invalid text in the focused Distance field: Enter shows the error and
    // nothing is applied, however long the preview is given.
    h.type_into_field_and_enter("Distance", "abc");
    let deadline = Instant::now() + Duration::from_millis(600);
    while Instant::now() < deadline {
        h.frame(vec![]);
        std::thread::sleep(Duration::from_millis(5));
    }
    assert_eq!(h.revision(), 0);
    assert!(h.modeling());
    let texts = h.texts();
    assert!(
        texts.iter().any(|t| t.starts_with("Enter a number")),
        "the invalid value's error is shown: {texts:?}"
    );
    // A valid value typed into the focused field: Enter commits it and applies.
    h.type_into_field_and_enter("Distance", "12");
    h.until("Enter to apply the extrusion", |app| {
        app.current().session.revision() == 1
    });
    assert!(!h.modeling());
    let applied = h.app.current().session.document().clone();
    assert!((distance_mm(&last_extrusion(&applied)) - 12.).abs() < 1e-9);
    // Enter with no text field focused applies the reopened tool too.
    h.until("the applied model", |app| {
        let p = app.current();
        p.scene
            .is_current(p.session.document().id(), p.session.revision())
    });
    open(&mut h);
    h.key(egui::Key::Enter);
    h.until("Enter to apply the second extrusion", |app| {
        app.current().session.revision() == 2
    });
    assert!(!h.modeling());
}
