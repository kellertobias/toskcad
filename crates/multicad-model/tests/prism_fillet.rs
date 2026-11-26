//! CAD-411: the reported body, built through the model commands that the UI
//! and MCP share. A sketch on the top face of a box holds a circle around a
//! face corner; the face clips it to a region with two straight edges and an
//! arc, which is extruded as a new body. Its top cap edges are then filleted
//! or chamfered.
use multicad_core::{
    commands::{Actor, Channel, Request, Session},
    constraint_sketch::{Curve, CurveKind, source::Source},
    document::{Document, EntityId},
    model::{
        self, BodyInput, EdgeTreatmentParameters, EdgeTreatmentProfile, EdgeTreatmentRequest,
        Extent, ExtrudeParameters, ExtrudeRequest, PrimitiveFrame, PrimitiveParameters,
        PrimitiveRequest, Termination, length,
    },
    references::{AnalyticGeometry, GeometryIntent, Origin, TopologyKind},
    sketch::{self, CreationPlane, PlaneOptions},
};
use multicad_kernel::{
    CurvedEdge2d, CurvedProfile, EdgeTreatment, NativeRegionQueries, PlaneFrame, Solid,
    TreatmentEdge,
};
use multicad_model::{Compiler, GeometryBackend, Model};

fn actor() -> Actor {
    Actor::local(Channel::Ui)
}

fn compile(document: &Document, revision: u64) -> Result<Model, String> {
    Compiler::default()
        .compile_with_backend(document, revision, GeometryBackend::Rust, None)
        .map_err(|error| error.to_string())
}

/// The reported session and the new body.
fn reported() -> (Session, BodyInput) {
    let mut session = Session::new(Document::new("Corner segment").unwrap()).unwrap();
    let stock = PrimitiveRequest::new(
        "Box",
        PrimitiveParameters::Box {
            size: [length(20.), length(20.), length(10.)],
            frame: PrimitiveFrame::default(),
        },
    );
    let prepared = model::prepare_primitive(&session, session.revision(), actor(), &stock).unwrap();
    session.commit(prepared.prepared).unwrap();
    let compiled = compile(session.document(), session.revision()).unwrap();
    let topology = &compiled.snapshot().shapes[&stock.feature].topology;
    let top = topology
        .elements
        .iter()
        .find(|e| {
            e.key.kind == TopologyKind::Face
                && matches!(e.geometry, Some(AnalyticGeometry::Plane { frame })
                    if frame.normal()[2] > 0.9 && (frame.origin_mm[2] - 10.).abs() < 1e-9)
        })
        .unwrap();
    let anchor = top.frames.first().unwrap().origin.clone();
    let created = sketch::prepare_new(
        &session,
        session.revision(),
        actor(),
        "Corner circle",
        CreationPlane::Face {
            snapshot: topology,
            key: top.key,
            anchor: Some(&anchor),
        },
        PlaneOptions::default(),
    )
    .unwrap();
    let (sketch, frame) = (created.sketch, created.frame);
    session.commit(created.prepared).unwrap();
    let report = compile(session.document(), session.revision())
        .unwrap()
        .snapshot()
        .references
        .clone();
    let mut source = Source::read_for_evaluation(
        session.document(),
        sketch,
        session.revision(),
        Some(&report),
    )
    .unwrap()
    .0;
    let center = source.add_point(frame.coordinates([2., 3., 10.]));
    source.sketch.curves.push(Curve {
        id: EntityId::new(),
        construction: false,
        kind: CurveKind::Circle {
            center,
            radius_mm: 8.,
        },
    });
    session
        .apply(Request {
            expected_revision: session.revision(),
            actor: actor(),
            label: "Circle around the face corner".into(),
            edits: source.replace(session.document(), sketch).unwrap(),
        })
        .unwrap();
    let report = compile(session.document(), session.revision())
        .unwrap()
        .snapshot()
        .references
        .clone();
    let queries = NativeRegionQueries::analytic(Default::default()).unwrap();
    let semantic = multicad_solver::evaluate_document_sketch(
        session.document(),
        session.revision(),
        sketch,
        Some(&report),
    )
    .unwrap();
    let inside = semantic.frame.coordinates([1., 1., 10.]);
    let region = semantic
        .regions_with_queries(&queries)
        .unwrap()
        .regions
        .into_iter()
        .find(|r| r.contains_with_queries(inside, &queries).unwrap())
        .expect("the clipped circle region");
    let request = ExtrudeRequest::new(
        sketch,
        vec![region.intent()],
        ExtrudeParameters {
            extent: Extent::OneSided {
                end: Termination::Distance {
                    distance: length(5.),
                },
            },
            ..Default::default()
        },
    );
    let output = BodyInput {
        body: request.outputs[0].body,
        producer: request.outputs[0].feature,
    };
    let prepared = model::prepare_extrude(&session, session.revision(), actor(), &request).unwrap();
    session.commit(prepared.prepared).unwrap();
    (session, output)
}

/// The top cap edges of the new body: the arc and its two straight edges.
fn top_edges(session: &Session, body: &BodyInput) -> Vec<GeometryIntent> {
    let compiled = compile(session.document(), session.revision()).unwrap();
    let topology = &compiled.snapshot().shapes[&body.producer].topology;
    let top = |g: AnalyticGeometry| match g {
        AnalyticGeometry::Line { a_mm, b_mm } => {
            (a_mm[2] - 15.).abs() < 1e-9 && (b_mm[2] - 15.).abs() < 1e-9
        }
        AnalyticGeometry::CircularArc { frame, .. } => (frame.origin_mm[2] - 15.).abs() < 1e-9,
        _ => false,
    };
    let edges: Vec<_> = topology
        .elements
        .iter()
        .filter(|e| e.key.kind == TopologyKind::Edge && e.geometry.is_some_and(top))
        .map(|e| GeometryIntent::capture(topology, e.key).unwrap())
        .collect();
    assert_eq!(edges.len(), 3, "arc and two straight top edges");
    edges
}

fn request(
    body: BodyInput,
    profile: EdgeTreatmentProfile,
    size: f64,
    edges: Vec<GeometryIntent>,
) -> EdgeTreatmentRequest {
    EdgeTreatmentRequest::new(
        body,
        EdgeTreatmentParameters {
            profile,
            size: length(size),
        },
        edges,
    )
}

/// The same body and treatment straight from the kernel, as a cross-check
/// that the model compiles the reported selection rather than another one.
fn direct_volume(treatment: EdgeTreatment) -> f64 {
    let id = EntityId::new();
    let (c, r) = ([2_f64, 3.], 8_f64);
    let x = c[0] + (r * r - c[1] * c[1]).sqrt();
    let y = c[1] + (r * r - c[0] * c[0]).sqrt();
    let outer = vec![
        CurvedEdge2d::Line {
            a: [0., 0.],
            b: [x, 0.],
        },
        CurvedEdge2d::Arc {
            center: c,
            a: [x, 0.],
            b: [0., y],
            clockwise: false,
        },
        CurvedEdge2d::Line {
            a: [0., y],
            b: [0., 0.],
        },
    ];
    let origins: Vec<_> = (0..3)
        .map(|i| Origin {
            source: id,
            role: format!("edge/{i}"),
        })
        .collect();
    let solid = Solid::rust_extrude_normal_profile_named(
        &CurvedProfile {
            frame: PlaneFrame {
                origin: [0., 0., 10.],
                ..Default::default()
            },
            outer,
            holes: vec![],
        },
        5.,
        id,
        &origins,
        Default::default(),
    )
    .unwrap();
    let snapshot = solid.topology_snapshot(id, 0).unwrap();
    let feature = EntityId::new();
    let edges: Vec<_> = snapshot
        .elements
        .iter()
        .filter(|e| {
            e.key.kind == TopologyKind::Edge
                && match e.geometry {
                    Some(AnalyticGeometry::Line { a_mm, b_mm }) => a_mm[2] == 15. && b_mm[2] == 15.,
                    Some(AnalyticGeometry::CircularArc { frame, .. }) => frame.origin_mm[2] == 15.,
                    _ => false,
                }
        })
        .enumerate()
        .map(|(i, e)| TreatmentEdge {
            key: e.key,
            origin: Origin {
                source: feature,
                role: format!("treated/{i}"),
            },
        })
        .collect();
    assert_eq!(edges.len(), 3);
    solid
        .treat_edges_named(&edges, treatment, feature, Default::default())
        .unwrap()
        .volume_mm3()
        .unwrap()
}

fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-7 * b.abs(), "{a} != {b}");
}

#[test]
fn reported_corner_segment_top_edges_fillet_and_chamfer_through_model_commands() {
    let (session, body) = reported();
    let revision = session.revision();
    let base = compile(session.document(), session.revision())
        .unwrap()
        .snapshot()
        .shapes[&body.producer]
        .volume_mm3;
    for (profile, treatment) in [
        (
            EdgeTreatmentProfile::Fillet,
            EdgeTreatment::Fillet { radius_mm: 1. },
        ),
        (
            EdgeTreatmentProfile::Chamfer,
            EdgeTreatment::Chamfer { distance_mm: 1. },
        ),
    ] {
        let edges = top_edges(&session, &body);
        let request = request(body.clone(), profile, 1., edges);
        let prepared =
            model::prepare_edge_treatment(&session, session.revision(), actor(), &request).unwrap();
        let compiled = compile(prepared.prepared.document(), session.revision() + 1)
            .unwrap_or_else(|error| panic!("{profile:?}: {error}"));
        let shape = &compiled.snapshot().shapes[&request.feature];
        assert!(shape.volume_mm3 < base);
        near(shape.volume_mm3, direct_volume(treatment));
        assert!(
            compiled
                .solid(request.feature)
                .unwrap()
                .supports_step_export()
        );
        // Preparing and compiling the candidate leaves the session untouched.
        assert_eq!(session.revision(), revision);
    }
}

#[test]
fn reported_selection_refusal_is_plain_language_and_keeps_the_request() {
    let (session, body) = reported();
    let edges = top_edges(&session, &body);
    // Taller than the 5 mm body.
    let request = request(body, EdgeTreatmentProfile::Fillet, 6., edges.clone());
    let prepared =
        model::prepare_edge_treatment(&session, session.revision(), actor(), &request).unwrap();
    let error = compile(prepared.prepared.document(), session.revision() + 1)
        .err()
        .expect("a 6 mm fillet cannot fit a 5 mm body");
    assert!(
        error.contains("too large for the height of the body"),
        "{error}"
    );
    assert!(!error.contains("rectangular box"), "{error}");
    assert_eq!(
        request.edges, edges,
        "the edge selection is kept for repair"
    );
}
