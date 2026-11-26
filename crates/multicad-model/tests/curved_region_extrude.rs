//! CAD-406 / CAD-407: a circle sketched on a box face is clipped by the face
//! and its region extrudes as Join or Cut through the model commands the
//! Extrude panel uses. The Cut's target discovery and preview evaluate the
//! exact Rust result; excluding the box leaves no target.
use multicad_core::{
    commands::{Actor, Channel, Request, Session},
    constraint_sketch::source::Source,
    document::{Document, EntityId},
    model::{
        self, BodyInput, CutTargetQuery, Extent, ExtrudeBooleanRequest, ExtrudeParameters,
        ExtrudeRequest, Operation, PrimitiveFrame, PrimitiveParameters, PrimitiveRequest,
        Termination, length,
    },
    references::TopologyKind,
    sketch::{self, CreationPlane, PlaneOptions},
    sketch_regions::RegionIntent,
};
use multicad_kernel::NativeRegionQueries;
use multicad_model::{Compiler, GeometryBackend, Model};
use multicad_precision::PortableMath as _;
use std::f64::consts::PI;

const SIZE: [f64; 3] = [40., 30., 20.];
const BOX: f64 = SIZE[0] * SIZE[1] * SIZE[2];

fn actor() -> Actor {
    Actor::local(Channel::Ui)
}

fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-6 * b.abs().max(1.), "{a} != {b}");
}

fn compile(document: &Document, revision: u64) -> Model {
    Compiler::default()
        .compile_with_backend(document, revision, GeometryBackend::Rust, None)
        .unwrap()
}

/// Area of a disk of radius `r` beyond a chord at distance `d` from its centre.
fn segment_area(r: f64, d: f64) -> f64 {
    r * r * (d / r).portable_acos() - d * (r * r - d * d).sqrt()
}

/// Area of the part of a disk with x > a and y > b relative to its centre.
fn corner_area(r: f64, a: f64, b: f64) -> f64 {
    let f = |x: f64| 0.5 * (x * (r * r - x * x).sqrt() + r * r * (x / r).portable_asin()) - b * x;
    f((r * r - b * b).sqrt()) - f(a)
}

/// Circles on the top face: centre, radius and the area the face keeps.
fn circles() -> Vec<(&'static str, [f64; 2], f64, f64)> {
    vec![
        ("corner centre", [40., 30.], 10., PI * 25.),
        (
            "around corner",
            [37., 27.],
            8.,
            PI * 64. - 2. * segment_area(8., 3.) + corner_area(8., 3., 3.),
        ),
        ("one edge", [20., 27.], 8., PI * 64. - segment_area(8., 3.)),
        ("inside", [20., 15.], 6., PI * 36.),
    ]
}

struct Fixture {
    session: Session,
    stock: PrimitiveRequest,
    sketch: EntityId,
}

impl Fixture {
    fn new(center: [f64; 2], radius: f64) -> Self {
        let mut session = Session::new(Document::new("Curved face region").unwrap()).unwrap();
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
        let snapshot = &compiled.snapshot().shapes[&stock.feature].topology;
        let face = compiled
            .solid(stock.feature)
            .unwrap()
            .faces()
            .unwrap()
            .into_iter()
            .find(|f| {
                f.plane
                    .is_some_and(|(p, n)| (p[2] - SIZE[2]).abs() < 1e-7 && n[2] > 0.9)
            })
            .unwrap();
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
            stock,
            sketch,
        }
    }

    /// The region the face clips from the circle: the smaller cell.
    fn circle_region(&self, area: f64) -> RegionIntent {
        let compiled = compile(self.session.document(), self.session.revision());
        let report = &compiled.snapshot().references;
        let regions = multicad_solver::evaluate_document_sketch(
            self.session.document(),
            self.session.revision(),
            self.sketch,
            Some(report),
        )
        .unwrap()
        .regions_with_queries(&NativeRegionQueries::analytic(Default::default()).unwrap())
        .unwrap();
        let region = regions
            .regions
            .iter()
            .min_by(|a, b| a.area_mm2.total_cmp(&b.area_mm2))
            .unwrap();
        near(region.area_mm2, area);
        region.intent()
    }

    fn extrusion(&self, region: RegionIntent, distance: f64) -> ExtrudeRequest {
        ExtrudeRequest::new(
            self.sketch,
            vec![region],
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

    fn target(&self) -> BodyInput {
        BodyInput {
            body: self.stock.body,
            producer: self.stock.feature,
        }
    }

    /// Prepare the extrusion with its Boolean, compile it as the preview
    /// does, and return the previewed result volume.
    fn preview(&self, request: &ExtrudeBooleanRequest) -> (model::ModelPrepared, f64) {
        let prepared = model::prepare_extrude_boolean(
            &self.session,
            self.session.revision(),
            actor(),
            request,
        )
        .unwrap();
        let before = self.session.document().clone();
        let preview = compile(prepared.prepared.document(), self.session.revision() + 1);
        assert_eq!(&before, self.session.document());
        let shape = &preview.snapshot().shapes[&request.booleans[0].feature];
        assert_eq!(shape.solid_count, 1);
        let volume = shape.volume_mm3;
        (prepared, volume)
    }
}

#[test]
fn a_face_clipped_circle_joins_the_box_exactly() {
    for (name, center, radius, area) in circles() {
        let mut f = Fixture::new(center, radius);
        let region = f.circle_region(area);
        let request = ExtrudeBooleanRequest::new(
            f.extrusion(region, 5.),
            Operation::Join,
            vec![f.target()],
            false,
        );
        let (prepared, previewed) = f.preview(&request);
        near(previewed, BOX + area * 5.);
        f.session.commit(prepared.prepared).unwrap();
        let model = compile(f.session.document(), f.session.revision());
        let head = model.snapshot().bodies[&f.stock.body];
        let shape = &model.snapshot().shapes[&head];
        near(shape.volume_mm3, BOX + area * 5.);
        assert_eq!(shape.solid_count, 1, "{name}");
        // The joined box keeps faces named after both operands.
        let joined = model.solid(head).unwrap();
        assert_eq!(
            joined
                .faces()
                .unwrap()
                .iter()
                .filter(|f| f.plane.is_none())
                .count(),
            1,
            "{name}"
        );
        let tool = request.extrusion.outputs[0].feature;
        assert!(shape.topology.elements.iter().any(
            |e| e.key.kind == TopologyKind::Face && e.origins.iter().any(|o| o.source == tool)
        ));
    }
}

#[test]
fn a_face_clipped_circle_cut_discovers_the_box_previews_and_commits() {
    for (name, center, radius, area) in circles() {
        for depth in [5., SIZE[2]] {
            let mut f = Fixture::new(center, radius);
            let region = f.circle_region(area);
            let extrusion = f.extrusion(region, -depth);
            // Target discovery runs on the candidate with the tool body.
            let tools =
                model::prepare_extrude(&f.session, f.session.revision(), actor(), &extrusion)
                    .unwrap();
            let document = tools.prepared.document();
            let revision = f.session.revision() + 1;
            let candidate = compile(document, revision);
            let query = |exclude: Vec<EntityId>| CutTargetQuery {
                tools: extrusion.outputs.iter().map(|o| o.feature).collect(),
                component: None,
                include: vec![],
                exclude,
            };
            let selection = candidate
                .cut_targets(document, revision, &query(vec![]), None)
                .unwrap_or_else(|e| panic!("{name}: {e}"));
            assert_eq!(selection.targets, vec![f.target()], "{name}");
            // Excluding the box leaves no target for the preview to cut.
            let excluded = candidate
                .cut_targets(document, revision, &query(vec![f.stock.body]), None)
                .unwrap();
            assert!(excluded.targets.is_empty(), "{name}");
            assert_eq!(excluded.defaults, vec![f.target()], "{name}");
            // The removed material is the region prism inside the box.
            let removed = candidate
                .solid(f.stock.feature)
                .unwrap()
                .rust_occupied_overlap_volume(
                    candidate.solid(extrusion.outputs[0].feature).unwrap(),
                    Default::default(),
                )
                .unwrap();
            near(removed, area * depth);
            let request =
                ExtrudeBooleanRequest::new(extrusion, Operation::Cut, selection.targets, false);
            let (prepared, previewed) = f.preview(&request);
            near(previewed, BOX - area * depth);
            f.session.commit(prepared.prepared).unwrap();
            let model = compile(f.session.document(), f.session.revision());
            let head = model.snapshot().bodies[&f.stock.body];
            near(
                model.snapshot().shapes[&head].volume_mm3,
                BOX - area * depth,
            );
        }
    }
}
