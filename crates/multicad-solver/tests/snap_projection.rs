//! CAD-383 snap-time projection: the face references a gesture snapped to,
//! captured against one base and merged, and the drawn geometry are one
//! checked, undoable transaction.
use multicad_core::{
    commands::{Actor, Channel, Prepared, Request, Session},
    constraint_sketch::{ConstraintKind, CurveKind, Dimension, source::Source},
    document::{Definition, Document, Entity, EntityData, EntityId},
    references::{
        self, AnalyticGeometry, GeometryIntent, Origin, TopologyElement, TopologyKey, TopologyKind,
        TopologySnapshot,
    },
    sketch::Datum,
    sketch_draw::{DrawInput, DrawRequest, PointInput},
    sketch_external::{self, Behavior, CaptureInput, CaptureRequest, Method},
};
use multicad_precision::PortableMath as _;
use multicad_solver::{
    ExternalEdit,
    draw::{prepare_draw_on, prepare_draw_on_with},
    prepare_external_reference,
};
use std::collections::{BTreeMap, BTreeSet};

struct Fixture {
    session: Session,
    sketch: EntityId,
    owner: EntityId,
    map: BTreeMap<EntityId, TopologySnapshot>,
}
fn fixture() -> Fixture {
    let (sketch, edits) = Source::create(Datum::Xy, "Face sketch").unwrap();
    let mut document = Document::new("Snap-time projection")
        .unwrap()
        .edited(&edits)
        .unwrap();
    let source = Entity::new(
        "Analytic face corners",
        EntityData::Feature {
            definition: Definition {
                type_id: "test.exact-descriptor".into(),
                schema_version: 1,
                data: serde_json::json!({}),
            },
            inputs: vec![],
            suppressed: false,
        },
    );
    let owner = source.id;
    document.insert(source).unwrap();
    let corner = |ordinal: u32, position_mm: [f64; 3]| TopologyElement {
        key: TopologyKey {
            kind: TopologyKind::Vertex,
            ordinal,
        },
        origins: BTreeSet::from([Origin {
            source: owner,
            role: format!("face.corner.{ordinal}"),
        }]),
        boundary: BTreeSet::new(),
        frames: vec![],
        geometry: Some(AnalyticGeometry::Point { position_mm }),
    };
    let map = BTreeMap::from([(
        owner,
        TopologySnapshot {
            owner,
            revision: 0,
            elements: vec![corner(0, [0., 0., 0.]), corner(1, [30., 20., 0.])].into(),
        },
    )]);
    Fixture {
        session: Session::new(document).unwrap(),
        sketch,
        owner,
        map,
    }
}
fn actor() -> Actor {
    Actor::local(Channel::Ui)
}
fn capture(f: &Fixture, ordinal: u32) -> (ExternalEdit, references::ReferenceReport, EntityId) {
    let key = TopologyKey {
        kind: TopologyKind::Vertex,
        ordinal,
    };
    let request = CaptureRequest::new(
        f.session.document(),
        f.session.revision(),
        CaptureInput {
            sketch: f.sketch,
            source: f.owner,
            intent: GeometryIntent::capture(&f.map[&f.owner], key).unwrap(),
            method: Method::Coplanar,
            behavior: Behavior::ReferenceOnly,
        },
        &f.map,
    )
    .unwrap();
    let point = request.output.points()[0];
    let edit =
        prepare_external_reference(&f.session, f.session.revision(), actor(), &request, &f.map)
            .unwrap();
    let report = references::inspect(
        edit.edit.preview.as_ref().unwrap().document(),
        f.session.revision(),
        &f.map,
    );
    (edit, report, point)
}
/// The captures merged against their common base into one candidate, and
/// the reference report that resolves it.
fn merged(f: &Fixture, captures: &[&ExternalEdit]) -> (Prepared, references::ReferenceReport) {
    let previews: Vec<_> = captures
        .iter()
        .map(|c| c.edit.preview.as_ref().unwrap().document())
        .collect();
    let edits =
        sketch_external::merge_capture_previews(f.session.document(), f.sketch, &previews).unwrap();
    let base = f
        .session
        .prepare(Request {
            expected_revision: f.session.revision(),
            actor: actor(),
            label: "Project face references".into(),
            edits,
        })
        .unwrap();
    let report = references::inspect(base.document(), f.session.revision(), &f.map);
    (base, report)
}
fn source(session: &Session, sketch: EntityId) -> Source {
    let EntityData::Sketch { definition, .. } = &session.document().get(sketch).unwrap().data
    else {
        panic!()
    };
    Source::from_definition(definition).unwrap()
}
fn rectangle(a: EntityId, b: EntityId) -> DrawRequest {
    DrawRequest {
        input: DrawInput::Rectangle {
            a: PointInput::Existing(a),
            b: PointInput::Existing(b),
        },
        construction: false,
    }
}

#[test]
fn rectangle_between_two_projected_face_corners_is_one_undoable_transaction() {
    let mut f = fixture();
    let (first, _, a) = capture(&f, 0);
    let (second, _, b) = capture(&f, 1);
    let before = f.session.document().clone();
    let (base, report) = merged(&f, &[&first, &second]);
    let drawn = prepare_draw_on(
        &f.session,
        f.session.revision(),
        actor(),
        f.sketch,
        &base,
        Some(&report),
        &rectangle(a, b),
    )
    .unwrap();
    assert_eq!(drawn.created_curves.len(), 4);
    assert_eq!(f.session.document(), &before, "preparing changes nothing");
    let solved = drawn.edit.evaluation.candidate.clone().unwrap();
    let corners: BTreeSet<_> = solved
        .curves
        .iter()
        .filter(|c| drawn.created_curves.contains(&c.id))
        .flat_map(|c| match c.kind {
            CurveKind::Line { start, end } => [start, end],
            _ => panic!("rectangle sides are lines"),
        })
        .collect();
    assert!(
        corners.contains(&a) && corners.contains(&b),
        "rectangle reuses both projected corner points by identity"
    );
    let at = |id| {
        solved
            .points
            .iter()
            .find(|p| p.id == id)
            .unwrap()
            .position_mm
    };
    assert_eq!(at(a), [0., 0.]);
    assert_eq!(at(b), [30., 20.]);
    f.session
        .commit(drawn.edit.preview.unwrap())
        .unwrap();
    let references: BTreeSet<_> = sketch_external::bindings(f.session.document(), f.sketch)
        .unwrap()
        .into_iter()
        .map(|binding| binding.reference)
        .collect();
    assert_eq!(
        references,
        BTreeSet::from([first.reference, second.reference])
    );
    assert_eq!(source(&f.session, f.sketch).sketch.curves.len(), 4);
    f.session.undo(f.session.revision(), actor()).unwrap();
    assert_eq!(
        f.session.document(),
        &before,
        "one Undo removes both references and the rectangle"
    );
    f.session.redo(f.session.revision(), actor()).unwrap();
    assert_eq!(
        sketch_external::bindings(f.session.document(), f.sketch)
            .unwrap()
            .len(),
        2
    );
}

#[test]
fn finish_hook_constraints_stale_captures_and_unsolvable_sizes_create_nothing() {
    let mut f = fixture();
    let (first, _, a) = capture(&f, 0);
    let (base, report) = merged(&f, &[&first]);
    let line = DrawRequest {
        input: DrawInput::Line {
            start: PointInput::Existing(a),
            end: PointInput::Position([10., 0.]),
        },
        construction: false,
    };
    let draft = |f: &Fixture, size: f64| {
        prepare_draw_on_with(
            &f.session,
            f.session.revision(),
            actor(),
            f.sketch,
            &base,
            Some(&report),
            &line,
            |source, created| {
                source.add_constraint(ConstraintKind::Length {
                    line: created[0],
                    dimension: Dimension::Literal {
                        value: size,
                        unit: "mm".into(),
                    },
                });
                Ok(())
            },
        )
    };
    let sized = draft(&f, 12.5).unwrap();
    let solved = sized.edit.evaluation.candidate.unwrap();
    let end = solved
        .curves
        .iter()
        .find_map(|c| match c.kind {
            CurveKind::Line { start, end } if start == a => Some(end),
            _ => None,
        })
        .unwrap();
    let end = solved.points.iter().find(|p| p.id == end).unwrap();
    assert!((end.position_mm[0].portable_hypot(end.position_mm[1]) - 12.5).abs() < 1e-7);
    // A hook error refuses the whole edit.
    let refused = prepare_draw_on_with(
        &f.session,
        f.session.revision(),
        actor(),
        f.sketch,
        &base,
        Some(&report),
        &line,
        |_, _| Err("entered sizes conflict".into()),
    );
    assert!(refused.is_err());
    // After any other commit, the held capture is stale and composes nothing.
    let (other, ..) = capture(&f, 1);
    f.session.commit(other.edit.preview.unwrap()).unwrap();
    let before = f.session.document().clone();
    assert!(draft(&f, 12.5).is_err());
    assert_eq!(f.session.document(), &before);
}
