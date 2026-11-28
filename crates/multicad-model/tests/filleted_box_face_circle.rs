//! CAD-424: a circle sketched on a face of a box whose twelve edges are all
//! filleted extrudes outward as Join and inward as Cut through the model
//! commands the Extrude panel uses.
use multicad_core::{
    commands::{Actor, Channel, Request, Session},
    constraint_sketch::source::Source,
    document::{Document, EntityId},
    model::{
        self, BodyInput, CutTargetQuery, EdgeTreatmentParameters, EdgeTreatmentProfile,
        EdgeTreatmentRequest, Extent, ExtrudeBooleanRequest, ExtrudeParameters, ExtrudeRequest,
        Operation, PrimitiveFrame, PrimitiveParameters, PrimitiveRequest, Termination, length,
    },
    references::{AnalyticGeometry, GeometryIntent, TopologyKind},
    sketch::{self, CreationPlane, PlaneOptions},
    sketch_regions::RegionIntent,
};
use multicad_kernel::NativeRegionQueries;
use multicad_model::{Compiler, GeometryBackend, Model};
use std::f64::consts::PI;

const SIZE: [f64; 3] = [40., 30., 20.];
const RADIUS: f64 = 3.;

fn actor() -> Actor {
    Actor::local(Channel::Ui)
}

fn compile(document: &Document, revision: u64) -> Model {
    Compiler::default()
        .compile_with_backend(document, revision, GeometryBackend::Rust, None)
        .unwrap()
}

struct Fixture {
    session: Session,
    body: BodyInput,
    sketch: EntityId,
    volume: f64,
}

impl Fixture {
    fn new(center: [f64; 2], radius: f64) -> Self {
        let mut session = Session::new(Document::new("Rounded box face circle").unwrap()).unwrap();
        let stock = PrimitiveRequest::new(
            "Box",
            PrimitiveParameters::Box {
                size: SIZE.map(length),
                frame: PrimitiveFrame::default(),
            },
        );
        let prepared =
            model::prepare_primitive(&session, session.revision(), actor(), &stock).unwrap();
        session.commit(prepared.prepared).unwrap();
        let compiled = compile(session.document(), session.revision());
        let topology = &compiled.snapshot().shapes[&stock.feature].topology;
        let edges: Vec<_> = topology
            .elements
            .iter()
            .filter(|e| {
                e.key.kind == TopologyKind::Edge
                    && matches!(e.geometry, Some(AnalyticGeometry::Line { .. }))
            })
            .map(|e| GeometryIntent::capture(topology, e.key).unwrap())
            .collect();
        assert_eq!(edges.len(), 12);
        let input = BodyInput {
            body: stock.body,
            producer: stock.feature,
        };
        let fillet = EdgeTreatmentRequest::new(
            input,
            EdgeTreatmentParameters {
                profile: EdgeTreatmentProfile::Fillet,
                size: length(RADIUS),
            },
            edges,
        );
        let prepared =
            model::prepare_edge_treatment(&session, session.revision(), actor(), &fillet).unwrap();
        session.commit(prepared.prepared).unwrap();
        let compiled = compile(session.document(), session.revision());
        let body = BodyInput {
            body: stock.body,
            producer: fillet.feature,
        };
        let shape = &compiled.snapshot().shapes[&fillet.feature];
        let volume = shape.volume_mm3;
        let snapshot = &shape.topology;
        let face = compiled
            .solid(fillet.feature)
            .unwrap()
            .faces()
            .unwrap()
            .into_iter()
            .find(|f| {
                f.plane
                    .is_some_and(|(p, n)| (p[2] - SIZE[2]).abs() < 1e-7 && n[2] > 0.9)
            })
            .expect("planar top face");
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
            "Face circle",
            CreationPlane::Face {
                snapshot,
                key: element.key,
                anchor: Some(&anchor),
            },
            PlaneOptions::default(),
        )
        .unwrap();
        let sketch = created.sketch;
        let frame = created.frame;
        session.commit(created.prepared).unwrap();
        let mut source = Source::read_for_evaluation(
            session.document(),
            sketch,
            session.revision(),
            Some(
                &compile(session.document(), session.revision())
                    .snapshot()
                    .references,
            ),
        )
        .unwrap()
        .0;
        let c = source.add_point(frame.coordinates([center[0], center[1], SIZE[2]]));
        source.add_circle(c, radius);
        session
            .apply(Request {
                expected_revision: session.revision(),
                actor: actor(),
                label: "Draw a face circle".into(),
                edits: source.replace(session.document(), sketch).unwrap(),
            })
            .unwrap();
        Self {
            session,
            body,
            sketch,
            volume,
        }
    }

    fn circle_region(&self) -> RegionIntent {
        let compiled = compile(self.session.document(), self.session.revision());
        let regions = multicad_solver::evaluate_document_sketch(
            self.session.document(),
            self.session.revision(),
            self.sketch,
            Some(&compiled.snapshot().references),
        )
        .unwrap()
        .regions_with_queries(&NativeRegionQueries::analytic(Default::default()).unwrap())
        .unwrap();
        regions
            .regions
            .iter()
            .min_by(|a, b| a.area_mm2.total_cmp(&b.area_mm2))
            .unwrap()
            .intent()
    }

    fn extrusion(&self, distance: f64) -> ExtrudeRequest {
        ExtrudeRequest::new(
            self.sketch,
            vec![self.circle_region()],
            ExtrudeParameters {
                extent: Extent::OneSided {
                    end: Termination::Distance {
                        distance: length(distance),
                    },
                },
                ..Default::default()
            },
        )
    }
}

fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-6 * b.abs().max(1.), "{a} != {b}");
}

/// Closed-form volume of the box with every edge rounded: the inner box
/// swept by a ball.
fn rounded_volume() -> f64 {
    let i = SIZE.map(|s| s - 2. * RADIUS);
    i[0] * i[1] * i[2]
        + 2. * RADIUS * (i[0] * i[1] + i[1] * i[2] + i[0] * i[2])
        + PI * RADIUS * RADIUS * (i[0] + i[1] + i[2])
        + 4. / 3. * PI * RADIUS.powi(3)
}

#[test]
fn a_circle_on_the_flat_face_of_a_rounded_box_joins_outward() {
    let mut f = Fixture::new([20., 15.], 6.);
    near(f.volume, rounded_volume());
    let request = ExtrudeBooleanRequest::new(
        f.extrusion(5.),
        Operation::Join,
        vec![f.body.clone()],
        false,
    );
    let prepared =
        model::prepare_extrude_boolean(&f.session, f.session.revision(), actor(), &request)
            .unwrap();
    // The preview compiles the candidate without touching the document.
    let preview = compile(prepared.prepared.document(), f.session.revision() + 1);
    let shape = &preview.snapshot().shapes[&request.booleans[0].feature];
    assert_eq!(shape.solid_count, 1);
    near(shape.volume_mm3, rounded_volume() + 36. * PI * 5.);
    f.session.commit(prepared.prepared).unwrap();
    let model = compile(f.session.document(), f.session.revision());
    let head = model.snapshot().bodies[&f.body.body];
    let shape = &model.snapshot().shapes[&head];
    assert_eq!(shape.solid_count, 1);
    near(shape.volume_mm3, rounded_volume() + 36. * PI * 5.);
    // Twelve edge rounds, eight corner blends and the boss wall.
    let joined = model.solid(head).unwrap();
    assert_eq!(
        joined
            .faces()
            .unwrap()
            .iter()
            .filter(|face| face.plane.is_none())
            .count(),
        21
    );
    let tool = request.extrusion.outputs[0].feature;
    assert!(shape.topology.elements.iter().any(
        |e| e.key.kind == TopologyKind::Face && e.origins.iter().any(|o| o.source == tool)
    ));
}

#[test]
fn a_circle_on_the_flat_face_of_a_rounded_box_cuts_inward() {
    for depth in [5., SIZE[2]] {
        let mut f = Fixture::new([20., 15.], 6.);
        let extrusion = f.extrusion(-depth);
        let tools =
            model::prepare_extrude(&f.session, f.session.revision(), actor(), &extrusion).unwrap();
        let document = tools.prepared.document();
        let revision = f.session.revision() + 1;
        let candidate = compile(document, revision);
        let selection = candidate
            .cut_targets(
                document,
                revision,
                &CutTargetQuery {
                    tools: extrusion.outputs.iter().map(|o| o.feature).collect(),
                    component: None,
                    include: vec![],
                    exclude: vec![],
                    contact: false,
                },
                None,
            )
            .unwrap();
        assert_eq!(selection.targets, vec![f.body.clone()]);
        let request =
            ExtrudeBooleanRequest::new(extrusion, Operation::Cut, selection.targets, false);
        let prepared =
            model::prepare_extrude_boolean(&f.session, f.session.revision(), actor(), &request)
                .unwrap();
        let preview = compile(prepared.prepared.document(), f.session.revision() + 1);
        let shape = &preview.snapshot().shapes[&request.booleans[0].feature];
        assert_eq!(shape.solid_count, 1);
        near(shape.volume_mm3, rounded_volume() - 36. * PI * depth);
        f.session.commit(prepared.prepared).unwrap();
        let model = compile(f.session.document(), f.session.revision());
        let head = model.snapshot().bodies[&f.body.body];
        near(
            model.snapshot().shapes[&head].volume_mm3,
            rounded_volume() - 36. * PI * depth,
        );
    }
}

#[test]
fn a_circle_reaching_the_fillet_is_refused_in_plain_words() {
    // The flat top face ends at y = 27; this circle reaches y = 30.
    let f = Fixture::new([20., 24.], 6.);
    let request = ExtrudeBooleanRequest::new(
        f.extrusion(5.),
        Operation::Join,
        vec![f.body.clone()],
        false,
    );
    let prepared =
        model::prepare_extrude_boolean(&f.session, f.session.revision(), actor(), &request)
            .unwrap();
    let error = Compiler::default()
        .compile_with_backend(
            prepared.prepared.document(),
            f.session.revision() + 1,
            GeometryBackend::Rust,
            None,
        )
        .err()
        .expect("the fillet zone is refused");
    assert!(error.message.contains("fillets"), "{error:?}");
}
