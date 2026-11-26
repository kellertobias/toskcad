//! CAD-404: Extrude picks Join, Cut or New body by default from the
//! extrusion's real contact with existing bodies. Every case runs the
//! production controller against the native model worker, applies on the
//! original Session and reads the canonical operation record back.
use super::*;
use crate::modeling::native_fixture::NativeBackend;
use multicad_core::{
    commands::Request,
    constraint_sketch::{Curve, CurveKind, source::Source},
    document::Document,
    model::{
        self, AutomaticOperation, BodyInput, ExtrudeGroupKind, Operation, PrimitiveFrame,
        PrimitiveParameters, PrimitiveRequest, length,
    },
    project_file::{self, ProjectFile},
    references::TopologyKind,
    sketch::{self, CreationPlane, Datum, PlaneOptions},
    sketch_regions,
};
use multicad_kernel::NativeRegionQueries;
use multicad_model::{Compiler, Model};
use std::time::{Duration, Instant};

fn actor() -> Actor {
    Actor::local(Channel::Ui)
}
fn compile(session: &Session) -> Model {
    Compiler::default()
        .compile(session.document(), session.revision())
        .unwrap()
}
fn box_part(session: &mut Session, name: &str, origin: [f64; 3], size: [f64; 3]) -> BodyInput {
    let request = PrimitiveRequest::new(
        name,
        PrimitiveParameters::Box {
            size: size.map(length),
            frame: PrimitiveFrame {
                origin: origin.map(length),
                ..Default::default()
            },
        },
    );
    let prepared =
        model::prepare_primitive(session, session.revision(), actor(), &request).unwrap();
    session.commit(prepared.prepared).unwrap();
    BodyInput {
        body: request.body,
        producer: request.feature,
    }
}

struct Harness {
    session: Session,
    controller: ModelToolController,
    stock: BodyInput,
}
impl Harness {
    /// A 20 × 20 × 10 box named "Box" with `draw` sketched on its top face.
    fn on_box_face(
        draw: impl FnOnce(&mut Source, sketch::Frame),
        pick: fn(f64) -> bool,
        others: &[(&str, [f64; 3], [f64; 3])],
    ) -> Self {
        let mut session =
            Session::new(Document::new("CAD-404 default operation").unwrap()).unwrap();
        let stock = box_part(&mut session, "Box", [0.; 3], [20., 20., 10.]);
        let compiled = compile(&session);
        let face = compiled
            .solid(stock.producer)
            .unwrap()
            .faces()
            .unwrap()
            .into_iter()
            .find(|f| {
                f.plane
                    .is_some_and(|(p, n)| (p[2] - 10.).abs() < 1e-7 && n[2] > 0.9)
            })
            .unwrap();
        let snapshot = &compiled.snapshot().shapes[&stock.producer].topology;
        let element = snapshot
            .elements
            .iter()
            .find(|e| e.key.kind == TopologyKind::Face && e.key.ordinal == face.index)
            .unwrap();
        let anchor = element.frames.first().unwrap().origin.clone();
        let created = sketch::prepare_new(
            &session,
            session.revision(),
            actor(),
            "Top face sketch",
            CreationPlane::Face {
                snapshot,
                key: element.key,
                anchor: Some(&anchor),
            },
            PlaneOptions::default(),
        )
        .unwrap();
        let (sketch, frame) = (created.sketch, created.frame);
        session.commit(created.prepared).unwrap();
        let mut source = Source::read_for_evaluation(
            session.document(),
            sketch,
            session.revision(),
            Some(&compile(&session).snapshot().references),
        )
        .unwrap()
        .0;
        draw(&mut source, frame);
        session
            .apply(Request {
                expected_revision: session.revision(),
                actor: actor(),
                label: "Draw on the box face".into(),
                edits: source.replace(session.document(), sketch).unwrap(),
            })
            .unwrap();
        let compiled = compile(&session);
        let region = multicad_solver::evaluate_document_sketch(
            session.document(),
            session.revision(),
            sketch,
            Some(&compiled.snapshot().references),
        )
        .unwrap()
        .regions_with_queries(&NativeRegionQueries::analytic(Default::default()).unwrap())
        .unwrap()
        .regions
        .iter()
        .find(|region| pick(region.area_mm2))
        .expect("picked face region")
        .intent();
        for (name, origin, size) in others {
            box_part(&mut session, name, *origin, *size);
        }
        Self::begin(session, stock, sketch, vec![region])
    }
    /// A 10 × 10 rectangle on the top face of the box.
    fn rectangle() -> Self {
        Self::rectangle_with(&[])
    }
    fn rectangle_with(others: &[(&str, [f64; 3], [f64; 3])]) -> Self {
        Self::on_box_face(
            |source, frame| {
                source
                    .add_rectangle(
                        frame.coordinates([5., 5., 10.]),
                        frame.coordinates([15., 15., 10.]),
                    )
                    .unwrap();
            },
            |area| (area - 100.).abs() < 1e-6,
            others,
        )
    }
    /// The product owner's model: a circle around the face corner. The face
    /// clips it to the partial disc that lies on the box.
    fn corner_circle() -> Self {
        Self::on_box_face(
            |source, frame| {
                let center = source.add_point(frame.coordinates([18., 18., 10.]));
                source.sketch.curves.push(Curve {
                    id: EntityId::new(),
                    construction: false,
                    kind: CurveKind::Circle {
                        center,
                        radius_mm: 6.,
                    },
                });
            },
            |area| area < 100.,
            &[],
        )
    }
    /// A datum XY sketch far from the box.
    fn separate() -> Self {
        let mut session = Session::new(Document::new("CAD-404 separate").unwrap()).unwrap();
        let stock = box_part(&mut session, "Box", [100., 100., 0.], [20., 20., 10.]);
        let (sketch, edits) = Source::create(Datum::Xy, "Free sketch").unwrap();
        session
            .apply(Request {
                expected_revision: session.revision(),
                actor: actor(),
                label: "Free sketch".into(),
                edits,
            })
            .unwrap();
        let mut source = Source::read(session.document(), sketch).unwrap().0;
        source.add_rectangle([0., 0.], [10., 10.]).unwrap();
        session
            .apply(Request {
                expected_revision: session.revision(),
                actor: actor(),
                label: "Draw".into(),
                edits: source.replace(session.document(), sketch).unwrap(),
            })
            .unwrap();
        let profiles = sketch_regions::extract(&source)
            .regions
            .iter()
            .map(|region| region.intent())
            .collect();
        Self::begin(session, stock, sketch, profiles)
    }
    fn begin(
        session: Session,
        stock: BodyInput,
        sketch: EntityId,
        profiles: Vec<sketch_regions::RegionIntent>,
    ) -> Self {
        let mut controller = ModelToolController::default();
        controller.set_backend(Box::new(NativeBackend::new()));
        controller
            .begin(&session, &Default::default(), sketch, profiles, Channel::Ui)
            .unwrap();
        Self {
            session,
            controller,
            stock,
        }
    }
    fn controls(&self) -> &OperationControls {
        self.controller.operations.as_ref().unwrap()
    }
    /// Edit the signed distance as the arrow or typed value does, preview, and
    /// wait until the final native candidate is ready to apply.
    fn distance(&mut self, millimeters: f64) {
        let Some(ModelForm::Extrude(panel)) = &mut self.controller.form else {
            panic!("extrude form")
        };
        panel.set_gizmo_distance(millimeters);
        self.controller
            .preview(&self.session, &Default::default(), Channel::Ui)
            .unwrap();
        let until = Instant::now() + Duration::from_secs(30);
        while !self
            .controller
            .can_apply(&self.session, &Default::default())
        {
            self.controller.poll(&self.session, &Default::default());
            if let Some(error) = self.controller.error() {
                panic!("{error}");
            }
            assert!(Instant::now() < until, "native preview did not settle");
            std::thread::sleep(Duration::from_millis(2));
        }
    }
    fn summary(&self) -> Option<String> {
        self.controls().automatic_summary()
    }
    /// Apply and return the one recorded extrusion operation.
    fn apply(&mut self) -> (EntityId, model::ExtrudeGroupState) {
        self.controller
            .apply(&mut self.session, &Default::default())
            .unwrap();
        group(self.session.document())
    }
}
fn group(document: &Document) -> (EntityId, model::ExtrudeGroupState) {
    let mut groups: Vec<_> = document
        .entities()
        .filter_map(|entity| {
            model::inspect_extrude_group(document, entity.id)
                .ok()
                .map(|state| (entity.id, state))
        })
        .collect();
    assert_eq!(groups.len(), 1);
    groups.pop().unwrap()
}
fn targets(state: &model::ExtrudeGroupState) -> Vec<EntityId> {
    state.definition.results.iter().map(|r| r.body).collect()
}
/// Save to project bytes and reopen, as File › Save and Open do.
fn reopen(document: &Document) -> Document {
    let mut bytes = vec![];
    project_file::write(
        &ProjectFile::new(document.clone()).unwrap(),
        &mut bytes,
        Default::default(),
    )
    .unwrap();
    project_file::read(
        &mut &bytes[..],
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap()
    .project
    .document
}

#[test]
fn outward_from_a_box_face_defaults_to_join_with_the_box_and_survives_save_and_reopen() {
    let mut h = Harness::rectangle();
    assert!(h.controls().automatic());
    h.distance(5.);
    assert_eq!(h.controls().operation, Operation::Join);
    assert_eq!(
        h.controls().decision(),
        Some(&AutomaticOperation::Join(h.stock.clone()))
    );
    assert_eq!(h.summary().as_deref(), Some("Join · touches Box"));
    let (id, state) = h.apply();
    assert_eq!(state.definition.kind, ExtrudeGroupKind::Boolean);
    assert_eq!(state.operation().unwrap(), Operation::Join);
    assert_eq!(targets(&state), vec![h.stock.body]);
    assert!(state.definition.automatic);
    // The canonical record survives save/reopen and reopens with its note.
    let reopened = reopen(h.session.document());
    let saved = model::inspect_extrude_group(&reopened, id).unwrap();
    assert_eq!(saved.definition, state.definition);
    assert_eq!(saved.schema_version, model::EXTRUDE_GROUP_V2_SCHEMA);
    let panel = crate::modeling::group_edit::GroupEditPanel::new(&reopened, saved).unwrap();
    assert_eq!(
        panel.automatic_note().as_deref(),
        Some("Automatic · Join · touches Box")
    );
}

#[test]
fn inward_into_a_box_face_defaults_to_cut_with_the_box() {
    let mut h = Harness::rectangle();
    h.distance(-5.);
    assert_eq!(h.controls().operation, Operation::Cut);
    assert_eq!(h.summary().as_deref(), Some("Cut · enters Box"));
    let (_, state) = h.apply();
    assert_eq!(state.operation().unwrap(), Operation::Cut);
    assert_eq!(targets(&state), vec![h.stock.body]);
    assert!(state.definition.automatic);
}

#[test]
fn an_extrusion_away_from_every_body_defaults_to_new_body() {
    let mut h = Harness::separate();
    h.distance(5.);
    assert_eq!(h.controls().operation, Operation::NewBody);
    assert_eq!(h.controls().decision(), Some(&AutomaticOperation::NewBody));
    assert_eq!(h.summary().as_deref(), Some("New body · touches no body"));
    let (_, state) = h.apply();
    assert_eq!(state.definition.kind, ExtrudeGroupKind::NewBody);
    assert!(state.definition.automatic);
}

#[test]
fn the_default_follows_distance_and_direction_edits_until_the_author_chooses() {
    let mut h = Harness::rectangle();
    h.distance(5.);
    assert_eq!(h.controls().operation, Operation::Join);
    h.distance(-5.);
    assert_eq!(h.controls().operation, Operation::Cut);
    h.distance(8.);
    assert_eq!(h.controls().operation, Operation::Join);
}

#[test]
fn an_explicit_operation_is_never_overridden_by_distance_or_direction_edits() {
    let mut h = Harness::rectangle();
    h.distance(5.);
    assert_eq!(h.controls().operation, Operation::Join);
    // The author picks New body although the extrusion touches the box.
    h.controller
        .operations
        .as_mut()
        .unwrap()
        .choose(Operation::NewBody);
    h.controller.operation_changed();
    assert!(!h.controls().automatic());
    assert_eq!(h.summary(), None);
    for distance in [-5., 7., -3.] {
        h.distance(distance);
        assert_eq!(h.controls().operation, Operation::NewBody, "{distance}");
        assert_eq!(h.controls().decision(), None);
    }
    let (_, state) = h.apply();
    assert_eq!(state.definition.kind, ExtrudeGroupKind::NewBody);
    assert!(!state.definition.automatic);
    // Choosing Automatic again hands the decision back to the next preview.
    let mut h = Harness::rectangle();
    let controls = h.controller.operations.as_mut().unwrap();
    controls.choose(Operation::Cut);
    controls.choose_automatic();
    h.distance(5.);
    assert_eq!(h.controls().operation, Operation::Join);
}

#[test]
fn several_touching_bodies_are_listed_and_fall_back_to_new_body() {
    // A lid resting above the face covers the extrusion's far end.
    let mut h = Harness::rectangle_with(&[("Lid", [0., 0., 15.], [20., 20., 5.])]);
    let lid = h
        .session
        .document()
        .entities()
        .find(|entity| entity.name == "Lid")
        .unwrap()
        .id;
    h.distance(5.);
    assert_eq!(h.controls().operation, Operation::NewBody);
    let Some(AutomaticOperation::Several(several)) = h.controls().decision() else {
        panic!("{:?}", h.controls().decision())
    };
    let mut bodies: Vec<_> = several.iter().map(|t| t.body).collect();
    bodies.sort();
    let mut expected = vec![h.stock.body, lid];
    expected.sort();
    assert_eq!(bodies, expected);
    let summary = h.summary().unwrap();
    assert!(
        summary.starts_with("New body · several bodies touch: "),
        "{summary}"
    );
    assert!(
        summary.contains("Box") && summary.contains("Lid"),
        "{summary}"
    );
    let (_, state) = h.apply();
    assert_eq!(state.definition.kind, ExtrudeGroupKind::NewBody);
}

#[test]
fn reported_corner_circle_joins_outward_and_cuts_inward_on_the_box() {
    for (distance, operation, summary) in [
        (5., Operation::Join, "Join · touches Box"),
        (-5., Operation::Cut, "Cut · enters Box"),
    ] {
        let mut h = Harness::corner_circle();
        let before = compile(&h.session)
            .solid(h.stock.producer)
            .unwrap()
            .volume_mm3()
            .unwrap();
        h.distance(distance);
        assert_eq!(h.controls().operation, operation, "{distance}");
        assert_eq!(h.summary().as_deref(), Some(summary));
        let (_, state) = h.apply();
        assert_eq!(state.operation().unwrap(), operation);
        assert_eq!(targets(&state), vec![h.stock.body]);
        // The applied Boolean changes the box's material by the disc prism.
        let compiled = compile(&h.session);
        let result = state.definition.results[0].feature;
        let after = compiled.solid(result).unwrap().volume_mm3().unwrap();
        let disc = 55.822171121641276 * 5.;
        let expected = if operation == Operation::Join {
            before + disc
        } else {
            before - disc
        };
        assert!((after - expected).abs() < 1e-4, "{after} != {expected}");
    }
}
