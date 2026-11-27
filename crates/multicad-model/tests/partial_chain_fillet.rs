//! CAD-421: fillets and chamfers on part of a tangent edge chain, through
//! the model commands that the UI and MCP share. The reported bodies are an
//! L whose inner corner is rounded in the sketch and a plate whose pocket
//! has rounded corners. One straight edge of the smooth inner chain is
//! treated, and the blend ends in flat faces at the rounded corners.
use multicad_core::{
    commands::{Actor, Channel, Session},
    constraint_sketch::{Curve, CurveKind, source::Source},
    document::{Document, EntityId},
    model::{
        self, BodyInput, EdgeTreatmentParameters, EdgeTreatmentProfile, EdgeTreatmentRequest,
        Extent, ExtrudeParameters, ExtrudeRequest, Termination, length,
    },
    references::{AnalyticGeometry, GeometryIntent, Origin, TopologyKind},
    sketch::Datum,
    sketch_regions,
};
use multicad_kernel::{
    CurvedEdge2d, CurvedProfile, EdgeTreatment, PlaneFrame, Solid, TreatmentEdge,
};
use multicad_model::{Compiler, GeometryBackend, Model};

/// Corner points, and for each edge from point `i` to `i + 1` an optional
/// arc center with its direction.
type Loop = (Vec<[f64; 2]>, Vec<Option<([f64; 2], bool)>>);

fn actor() -> Actor {
    Actor::local(Channel::Mcp)
}

fn compile(document: &Document, revision: u64) -> Result<Model, String> {
    Compiler::default()
        .compile_with_backend(document, revision, GeometryBackend::Rust, None)
        .map_err(|error| error.to_string())
}

fn rounded_l() -> Vec<Loop> {
    vec![(
        vec![
            [0., 0.],
            [20., 0.],
            [20., 8.],
            [11., 8.],
            [8., 11.],
            [8., 20.],
            [0., 20.],
        ],
        vec![None, None, None, Some(([11., 11.], true)), None, None, None],
    )]
}

fn pocket_plate() -> Vec<Loop> {
    vec![
        (
            vec![[0., 0.], [30., 0.], [30., 20.], [0., 20.]],
            vec![None; 4],
        ),
        (
            vec![
                [12., 6.],
                [18., 6.],
                [20., 8.],
                [20., 12.],
                [18., 14.],
                [12., 14.],
                [10., 12.],
                [10., 8.],
            ],
            vec![
                None,
                Some(([18., 8.], false)),
                None,
                Some(([18., 12.], false)),
                None,
                Some(([12., 12.], false)),
                None,
                Some(([12., 8.], false)),
            ],
        ),
    ]
}

const HEIGHT: f64 = 5.;

/// A top-plane sketch of `loops` extruded as a new body.
fn session(loops: &[Loop]) -> (Session, BodyInput) {
    let (sketch, edits) = Source::create(Datum::Xy, "Profile").unwrap();
    let document = Document::new("Partial chain")
        .unwrap()
        .edited(&edits)
        .unwrap();
    let mut profile = Source::read(&document, sketch).unwrap().0;
    for (points, arcs) in loops {
        let ids: Vec<_> = points.iter().map(|&p| profile.add_point(p)).collect();
        for (i, arc) in arcs.iter().enumerate() {
            let (a, b) = (ids[i], ids[(i + 1) % ids.len()]);
            match *arc {
                None => {
                    profile.add_line(a, b);
                }
                Some((center, clockwise)) => {
                    let center = profile.add_point(center);
                    profile.sketch.curves.push(Curve {
                        id: EntityId::new(),
                        construction: false,
                        kind: CurveKind::Arc {
                            center,
                            start: a,
                            end: b,
                            clockwise,
                        },
                    });
                }
            }
        }
    }
    let region = sketch_regions::extract(&profile)
        .regions
        .into_iter()
        .find(|r| r.holes.len() == loops.len() - 1)
        .unwrap()
        .intent();
    let document = document
        .edited(&profile.replace(&document, sketch).unwrap())
        .unwrap();
    let mut session = Session::new(document).unwrap();
    let request = ExtrudeRequest::new(
        sketch,
        vec![region],
        ExtrudeParameters {
            extent: Extent::OneSided {
                end: Termination::Distance {
                    distance: length(HEIGHT),
                },
            },
            ..Default::default()
        },
    );
    let body = BodyInput {
        body: request.outputs[0].body,
        producer: request.outputs[0].feature,
    };
    let prepared = model::prepare_extrude(&session, session.revision(), actor(), &request).unwrap();
    session.commit(prepared.prepared).unwrap();
    (session, body)
}

/// The line edge from `a` to `b` at height `z`, in either direction.
fn is_line(geometry: Option<AnalyticGeometry>, z: f64, a: [f64; 2], b: [f64; 2]) -> bool {
    let at = |p: [f64; 3], q: [f64; 2]| {
        (p[0] - q[0]).abs() < 1e-9 && (p[1] - q[1]).abs() < 1e-9 && (p[2] - z).abs() < 1e-9
    };
    matches!(geometry, Some(AnalyticGeometry::Line { a_mm, b_mm })
        if (at(a_mm, a) && at(b_mm, b)) || (at(a_mm, b) && at(b_mm, a)))
}

fn edge(session: &Session, body: &BodyInput, z: f64, a: [f64; 2], b: [f64; 2]) -> GeometryIntent {
    let compiled = compile(session.document(), session.revision()).unwrap();
    let topology = &compiled.snapshot().shapes[&body.producer].topology;
    let found = topology
        .elements
        .iter()
        .find(|e| e.key.kind == TopologyKind::Edge && is_line(e.geometry, z, a, b))
        .expect("profile line edge");
    GeometryIntent::capture(topology, found.key).unwrap()
}

/// The same body and treatment straight from the kernel.
fn direct_volume(
    loops: &[Loop],
    z: f64,
    a: [f64; 2],
    b: [f64; 2],
    treatment: EdgeTreatment,
) -> f64 {
    let id = EntityId::new();
    let wires: Vec<Vec<CurvedEdge2d>> = loops
        .iter()
        .map(|(points, arcs)| {
            (0..points.len())
                .map(|i| {
                    let (a, b) = (points[i], points[(i + 1) % points.len()]);
                    match arcs[i] {
                        None => CurvedEdge2d::Line { a, b },
                        Some((center, clockwise)) => CurvedEdge2d::Arc {
                            center,
                            a,
                            b,
                            clockwise,
                        },
                    }
                })
                .collect()
        })
        .collect();
    let count = wires.iter().map(Vec::len).sum::<usize>();
    let origins: Vec<_> = (0..count)
        .map(|i| Origin {
            source: id,
            role: format!("edge/{i}"),
        })
        .collect();
    let solid = Solid::rust_extrude_normal_profile_named(
        &CurvedProfile {
            frame: PlaneFrame::default(),
            outer: wires[0].clone(),
            holes: wires[1..].to_vec(),
        },
        HEIGHT,
        id,
        &origins,
        Default::default(),
    )
    .unwrap();
    let snapshot = solid.topology_snapshot(id, 0).unwrap();
    let feature = EntityId::new();
    let key = snapshot
        .elements
        .iter()
        .find(|e| e.key.kind == TopologyKind::Edge && is_line(e.geometry, z, a, b))
        .unwrap()
        .key;
    let edges = [TreatmentEdge {
        key,
        origin: Origin {
            source: feature,
            role: "treated/0".into(),
        },
    }];
    solid
        .treat_edges_named(&edges, treatment, feature, Default::default())
        .unwrap()
        .volume_mm3()
        .unwrap()
}

fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-7 * b.abs(), "{a} != {b}");
}

fn treats_one_chain_edge(loops: &[Loop], z: f64, a: [f64; 2], b: [f64; 2]) {
    let (session, body) = session(loops);
    let revision = session.revision();
    let base = compile(session.document(), revision)
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
        let request = EdgeTreatmentRequest::new(
            body.clone(),
            EdgeTreatmentParameters {
                profile,
                size: length(1.),
            },
            vec![edge(&session, &body, z, a, b)],
        );
        let prepared =
            model::prepare_edge_treatment(&session, session.revision(), actor(), &request).unwrap();
        let compiled = compile(prepared.prepared.document(), revision + 1)
            .unwrap_or_else(|error| panic!("{profile:?}: {error}"));
        let shape = &compiled.snapshot().shapes[&request.feature];
        assert!(shape.volume_mm3 < base);
        near(shape.volume_mm3, direct_volume(loops, z, a, b, treatment));
        assert!(
            compiled
                .solid(request.feature)
                .unwrap()
                .supports_step_export()
        );
        assert_eq!(session.revision(), revision);
    }
}

#[test]
fn rounded_l_single_inner_edge_fillet_and_chamfer_on_both_caps_through_model_commands() {
    treats_one_chain_edge(&rounded_l(), HEIGHT, [20., 8.], [11., 8.]);
    treats_one_chain_edge(&rounded_l(), 0., [8., 11.], [8., 20.]);
}

#[test]
fn pocket_single_inner_side_fillet_and_chamfer_on_both_caps_through_model_commands() {
    treats_one_chain_edge(&pocket_plate(), HEIGHT, [12., 6.], [18., 6.]);
    treats_one_chain_edge(&pocket_plate(), 0., [10., 8.], [10., 12.]);
}

#[test]
fn partial_chain_refusal_is_plain_language_and_keeps_the_request() {
    let (session, body) = session(&rounded_l());
    let edges = vec![edge(&session, &body, HEIGHT, [20., 8.], [11., 8.])];
    let request = EdgeTreatmentRequest::new(
        body,
        EdgeTreatmentParameters {
            profile: EdgeTreatmentProfile::Fillet,
            size: length(6.),
        },
        edges.clone(),
    );
    let prepared =
        model::prepare_edge_treatment(&session, session.revision(), actor(), &request).unwrap();
    let error = compile(prepared.prepared.document(), session.revision() + 1)
        .err()
        .expect("a 6 mm fillet cannot fit a 5 mm body");
    assert!(
        error.contains("too large for the height of the body"),
        "{error}"
    );
    assert!(!error.contains("tangent chain"), "{error}");
    assert_eq!(
        request.edges, edges,
        "the edge selection is kept for repair"
    );
}
