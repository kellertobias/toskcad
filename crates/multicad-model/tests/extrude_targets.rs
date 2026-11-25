//! CAD-405: automatic Join/Cut targets of an extrusion that starts on a
//! body's own face. Every hit is decided on genuine native solids.
use multicad_core::{
    commands::{Actor, Channel, Request, Session},
    constraint_sketch::{Curve, CurveKind, source::Source},
    document::{Document, EntityId},
    model::{
        self, AutomaticTargetError, BodyInput, CutTargetQuery, CutTargetSelection, Extent,
        ExtrudeBooleanRequest, ExtrudeParameters, ExtrudeRequest, Operation, PrimitiveFrame,
        PrimitiveParameters, PrimitiveRequest, Termination, length,
    },
    references::TopologyKind,
    sketch::{self, CreationPlane, PlaneOptions},
    sketch_regions::RegionIntent,
};
use multicad_kernel::NativeRegionQueries;
use multicad_model::{Compiler, Model};

fn actor() -> Actor {
    Actor::local(Channel::Mcp)
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-5, "{a} != {b}");
}
fn compile(s: &Session) -> Model {
    Compiler::default()
        .compile(s.document(), s.revision())
        .unwrap()
}
fn box_part(s: &mut Session, name: &str, origin: [f64; 3], size: [f64; 3]) -> BodyInput {
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
    let prepared = model::prepare_primitive(s, s.revision(), actor(), &request).unwrap();
    s.commit(prepared.prepared).unwrap();
    BodyInput {
        body: request.body,
        producer: request.feature,
    }
}
fn distance(value: f64) -> ExtrudeParameters {
    ExtrudeParameters {
        extent: Extent::OneSided {
            end: Termination::Distance {
                distance: length(value),
            },
        },
        ..Default::default()
    }
}

/// The reported model: a 20 × 20 × 10 box with a sketch on its top face.
struct Fixture {
    session: Session,
    stock: BodyInput,
    sketch: EntityId,
    frame: sketch::Frame,
}
impl Fixture {
    fn new() -> Self {
        let mut session = Session::new(Document::new("CAD-405 face extrusion").unwrap()).unwrap();
        let stock = box_part(&mut session, "Box", [0.; 3], [20., 20., 10.]);
        let compiled = compile(&session);
        let snapshot = &compiled.snapshot().shapes[&stock.producer].topology;
        let element = top_face_element(&compiled, &stock);
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
        Self {
            session,
            stock,
            sketch,
            frame,
        }
    }
    fn draw(&mut self, edit: impl FnOnce(&mut Source, sketch::Frame)) {
        let mut source = Source::read_for_evaluation(
            self.session.document(),
            self.sketch,
            self.session.revision(),
            Some(&compile(&self.session).snapshot().references),
        )
        .unwrap()
        .0;
        edit(&mut source, self.frame);
        self.session
            .apply(Request {
                expected_revision: self.session.revision(),
                actor: actor(),
                label: "Draw on the box face".into(),
                edits: source
                    .replace(self.session.document(), self.sketch)
                    .unwrap(),
            })
            .unwrap();
    }
    /// A circle centred near the face corner; the face clips it to one cell.
    fn corner_circle() -> (Self, RegionIntent) {
        let mut f = Self::new();
        f.draw(|source, frame| {
            let center = source.add_point(frame.coordinates([18., 18., 10.]));
            source.sketch.curves.push(Curve {
                id: EntityId::new(),
                construction: false,
                kind: CurveKind::Circle {
                    center,
                    radius_mm: 6.,
                },
            });
        });
        let regions = f.regions();
        // Two face cells: the partial disc and the face around it.
        assert_eq!(regions.len(), 2);
        let disc = regions
            .into_iter()
            .find(|(_, area)| *area < 100.)
            .expect("partial circle region");
        near(disc.1, 55.822171121641276);
        (f, disc.0)
    }
    fn rectangle() -> (Self, RegionIntent) {
        let mut f = Self::new();
        f.draw(|source, frame| {
            source
                .add_rectangle(
                    frame.coordinates([5., 5., 10.]),
                    frame.coordinates([15., 15., 10.]),
                )
                .unwrap();
        });
        let rectangle = f
            .regions()
            .into_iter()
            .find(|(_, area)| (area - 100.).abs() < 1e-6)
            .expect("rectangle region");
        (f, rectangle.0)
    }
    fn regions(&self) -> Vec<(RegionIntent, f64)> {
        let compiled = compile(&self.session);
        let report = &compiled.snapshot().references;
        multicad_solver::evaluate_document_sketch(
            self.session.document(),
            self.session.revision(),
            self.sketch,
            Some(report),
        )
        .unwrap()
        .regions_with_queries(&NativeRegionQueries::analytic(Default::default()).unwrap())
        .unwrap()
        .regions
        .iter()
        .map(|region| (region.intent(), region.area_mm2))
        .collect()
    }
    fn extrusion(&self, region: &RegionIntent, value: f64) -> ExtrudeRequest {
        ExtrudeRequest::new(self.sketch, vec![region.clone()], distance(value))
    }
    fn discover(&self, extrusion: &ExtrudeRequest, contact: bool) -> CutTargetSelection {
        discover(&self.session, extrusion, contact).unwrap()
    }
    fn add_box(&mut self, name: &str, origin: [f64; 3], size: [f64; 3]) -> BodyInput {
        box_part(&mut self.session, name, origin, size)
    }
}
fn top_face_element<'a>(
    compiled: &'a Model,
    body: &BodyInput,
) -> &'a multicad_core::references::TopologyElement {
    let face = compiled
        .solid(body.producer)
        .unwrap()
        .faces()
        .unwrap()
        .into_iter()
        .find(|f| {
            f.plane
                .is_some_and(|(p, n)| (p[2] - 10.).abs() < 1e-7 && n[2] > 0.9)
        })
        .unwrap();
    compiled.snapshot().shapes[&body.producer]
        .topology
        .elements
        .iter()
        .find(|e| e.key.kind == TopologyKind::Face && e.key.ordinal == face.index)
        .unwrap()
}
/// Discovery exactly as the window and MCP run it: on the tools-only
/// candidate, at the candidate revision, in the document-root scope.
fn discover(
    session: &Session,
    extrusion: &ExtrudeRequest,
    contact: bool,
) -> Result<CutTargetSelection, multicad_model::ModelError> {
    let tools = model::prepare_extrude(session, session.revision(), actor(), extrusion).unwrap();
    let document = tools.prepared.document();
    let revision = session.revision() + 1;
    Compiler::default()
        .compile(document, revision)
        .unwrap()
        .cut_targets(
            document,
            revision,
            &CutTargetQuery {
                tools: extrusion.outputs.iter().map(|o| o.feature).collect(),
                component: None,
                include: vec![],
                exclude: vec![],
                contact,
            },
            None,
        )
}
fn bodies(selection: &CutTargetSelection) -> Vec<EntityId> {
    selection.targets.iter().map(|t| t.body).collect()
}
fn apply_boolean(
    session: &Session,
    extrusion: ExtrudeRequest,
    operation: Operation,
    targets: Vec<BodyInput>,
) -> (Model, EntityId) {
    let request = ExtrudeBooleanRequest::new(extrusion, operation, targets, false);
    let feature = request.booleans[0].feature;
    let prepared =
        model::prepare_extrude_boolean(session, session.revision(), actor(), &request).unwrap();
    let model = Compiler::default()
        .compile(prepared.prepared.document(), session.revision() + 1)
        .unwrap();
    (model, feature)
}

#[test]
fn corner_circle_join_discovers_the_box_that_owns_the_sketch_face() {
    let (mut f, disc) = Fixture::corner_circle();
    let far = f.add_box("Disjoint", [100., 0., 0.], [10., 10., 10.]);
    let outward = f.extrusion(&disc, 5.);
    let join = f.discover(&outward, true);
    assert_eq!(join.targets, vec![f.stock.clone()]);
    assert_eq!(join.single_target(), Ok(&f.stock));
    assert!(!bodies(&join).contains(&far.body));
    // Contact alone removes no material, so Cut finds nothing to cut.
    assert!(f.discover(&outward, false).targets.is_empty());
}

#[test]
fn cut_discovers_the_body_an_inward_corner_extrusion_enters() {
    let (mut f, disc) = Fixture::corner_circle();
    f.add_box("Disjoint", [100., 0., 0.], [10., 10., 10.]);
    let inward = f.extrusion(&disc, -5.);
    assert_eq!(f.discover(&inward, false).targets, vec![f.stock.clone()]);
    assert_eq!(f.discover(&inward, true).targets, vec![f.stock.clone()]);
}

#[test]
fn planar_join_and_cut_from_the_face_use_the_discovered_box() {
    let (f, rectangle) = Fixture::rectangle();
    let outward = f.extrusion(&rectangle, 5.);
    let join = f.discover(&outward, true);
    let target = join.single_target().unwrap().clone();
    assert_eq!(target, f.stock);
    let (model, feature) = apply_boolean(&f.session, outward, Operation::Join, vec![target]);
    let shape = &model.snapshot().shapes[&feature];
    near(shape.volume_mm3, 20. * 20. * 10. + 10. * 10. * 5.);
    assert_eq!(shape.solid_count, 1);
    near(shape.bounds.max[2], 15.);

    let inward = f.extrusion(&rectangle, -5.);
    let cut = f.discover(&inward, false);
    assert_eq!(cut.targets, vec![f.stock.clone()]);
    let (model, feature) = apply_boolean(&f.session, inward, Operation::Cut, cut.targets);
    let shape = &model.snapshot().shapes[&feature];
    near(shape.volume_mm3, 20. * 20. * 10. - 10. * 10. * 5.);
    assert_eq!(shape.solid_count, 1);
}

#[test]
fn touching_only_contact_counts_for_join_but_edge_contact_and_disjoint_bodies_do_not() {
    let (mut f, rectangle) = Fixture::rectangle();
    // The tool spans x, y 5..15 and z 10..15. This block meets it along one
    // vertical edge only; the other is far away.
    let edge = f.add_box("Edge contact", [15., 15., 10.], [10., 10., 5.]);
    let far = f.add_box("Disjoint", [100., 0., 0.], [10., 10., 10.]);
    let outward = f.extrusion(&rectangle, 5.);
    let join = f.discover(&outward, true);
    assert_eq!(join.targets, vec![f.stock.clone()]);
    for other in [&edge, &far] {
        assert!(!bodies(&join).contains(&other.body));
    }
    // A lid resting exactly on the extrusion's end cap: shared face, zero
    // common volume. Both bodies are now plausible Join targets.
    let lid = f.add_box("Lid", [0., 0., 15.], [20., 20., 4.]);
    let join = f.discover(&outward, true);
    let mut expected = vec![f.stock.body, lid.body];
    expected.sort();
    let mut found = bodies(&join);
    found.sort();
    assert_eq!(found, expected);
    let Err(AutomaticTargetError::Several(listed)) = join.single_target() else {
        panic!("several touching bodies must be listed, not chosen");
    };
    assert_eq!(listed, join.targets);
    // Cut keeps volumetric semantics: neither contact removes material.
    assert!(f.discover(&outward, false).targets.is_empty());
}

#[test]
fn a_face_profile_joins_its_own_body_and_cuts_into_it() {
    let mut session = Session::new(Document::new("CAD-405 face profile").unwrap()).unwrap();
    let stock = box_part(&mut session, "Box", [0.; 3], [20., 20., 10.]);
    let other = box_part(&mut session, "Disjoint", [100., 0., 0.], [10., 10., 10.]);
    let compiled = compile(&session);
    let element = top_face_element(&compiled, &stock);
    let face = model::capture_profile_face(
        &compiled.snapshot().shapes[&stock.producer].topology,
        element.key,
        None,
    )
    .unwrap();
    let outward = ExtrudeRequest::from_faces(vec![face.clone()], distance(5.));
    let join = discover(&session, &outward, true).unwrap();
    assert_eq!(join.targets, vec![stock.clone()]);
    assert!(!bodies(&join).contains(&other.body));
    let (model, feature) = apply_boolean(&session, outward, Operation::Join, join.targets);
    near(model.snapshot().shapes[&feature].volume_mm3, 6000.);
    let inward = ExtrudeRequest::from_faces(vec![face], distance(-5.));
    assert_eq!(
        discover(&session, &inward, false).unwrap().targets,
        vec![stock]
    );
}

#[test]
fn an_offset_sketch_plane_is_not_treated_as_lying_on_the_face() {
    let mut session = Session::new(Document::new("CAD-405 offset").unwrap()).unwrap();
    let stock = box_part(&mut session, "Box", [0.; 3], [20., 20., 10.]);
    let compiled = compile(&session);
    let snapshot = &compiled.snapshot().shapes[&stock.producer].topology;
    let element = top_face_element(&compiled, &stock);
    let anchor = element.frames.first().unwrap().origin.clone();
    let created = sketch::prepare_new(
        &session,
        session.revision(),
        actor(),
        "Offset sketch",
        CreationPlane::Face {
            snapshot,
            key: element.key,
            anchor: Some(&anchor),
        },
        PlaneOptions {
            normal_offset_mm: 2.,
            angle_rad: 0.,
        },
    )
    .unwrap();
    let mut f = Fixture {
        session,
        stock,
        sketch: created.sketch,
        frame: created.frame,
    };
    f.session.commit(created.prepared).unwrap();
    f.draw(|source, frame| {
        source
            .add_rectangle(
                frame.coordinates([5., 5., 12.]),
                frame.coordinates([15., 15., 12.]),
            )
            .unwrap();
    });
    let rectangle = f
        .regions()
        .into_iter()
        .find(|(_, area)| (area - 100.).abs() < 1e-6)
        .expect("rectangle region")
        .0;
    // Floating 2 mm above the face, the tool touches nothing.
    assert!(
        f.discover(&f.extrusion(&rectangle, 5.), true)
            .targets
            .is_empty()
    );
}
