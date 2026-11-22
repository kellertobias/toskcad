use multicad_core::{
    document::EntityId,
    drawings::{VectorCurve, VectorRole},
    project_file,
    references::{AnalyticGeometry, TopologyKind},
};
use multicad_kernel::{
    DrawingProjection, DrawingSolid, GeometryPolicy, MeshOptions, PlaneFrame, ProjectedVector,
    Solid, project_drawing,
};
use multicad_model::{Compiler, GeometryBackend};
use std::f64::consts::FRAC_PI_2;

#[path = "support/saved_prism.rs"]
mod saved_prism;

fn toolcase() -> multicad_core::document::Document {
    project_file::read(
        &mut include_bytes!("../../../examples/toolcase.tcd").as_slice(),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap()
    .project
    .document
}

/// One expected projected profile edge: its 2D ends and, for the quarter
/// arc, its center and radius.
struct ExpectedEdge {
    ends: [[f64; 2]; 2],
    arc: Option<([f64; 2], f64)>,
}

fn matches_expected(actual: &ProjectedVector, expected: &ExpectedEdge) -> bool {
    let near = |a: f64, b: f64| (a - b).abs() < 1e-6;
    let same_point = |a: [f64; 2], b: [f64; 2]| near(a[0], b[0]) && near(a[1], b[1]);
    let [a, b] = [actual.preview_mm[0], *actual.preview_mm.last().unwrap()];
    let [c, d] = expected.ends;
    if !(same_point(a, c) && same_point(b, d) || same_point(a, d) && same_point(b, c)) {
        return false;
    }
    match (&actual.curve, expected.arc) {
        (
            VectorCurve::Ellipse {
                center,
                major_mm,
                minor_mm,
                start_rad,
                end_rad,
                ..
            },
            Some((expected_center, radius)),
        ) => {
            same_point(*center, expected_center)
                && near(*major_mm, radius)
                && near(*minor_mm, radius)
                && near(end_rad - start_rad, FRAC_PI_2)
        }
        (VectorCurve::Line { .. }, None) => true,
        _ => false,
    }
}

/// The saved profile's edges projected onto a drawing plane with axes
/// `(u, v)`: the near and far caps of a straight sweep along the view
/// direction both project onto these.
fn projected_profile(prism: &saved_prism::Prism, u: [f64; 3], v: [f64; 3]) -> Vec<ExpectedEdge> {
    let dot = |a: [f64; 3], b: [f64; 3]| (0..3).map(|i| a[i] * b[i]).sum::<f64>();
    let project = |p: [f64; 3]| [dot(p, u), dot(p, v)];
    let n = prism.loop_uv.len();
    (0..n)
        .map(|i| {
            let [a, b] = [prism.loop_mm[i], prism.loop_mm[(i + 1) % n]];
            ExpectedEdge {
                ends: [project(a), project(b)],
                arc: prism.arcs[i].map(|center| {
                    let c = prism.world(center);
                    let radius = (0..3).map(|k| (a[k] - c[k]).powi(2)).sum::<f64>().sqrt();
                    (project(c), radius)
                }),
            }
        })
        .collect()
}

#[test]
fn saved_toolcase_floor_profiles_rebuild_and_round_trip_with_analytic_quarter_arcs() {
    let document = toolcase();
    let features = [
        "015a44ba-56d1-477e-aaa0-c9016233401a",
        "16cc5d02-e3cc-472f-a2bd-9a1334340d67",
    ]
    .map(|id| EntityId::parse(id).unwrap());
    let compile = |backend| {
        Compiler::default()
            .compile_sources_with_backend(&document, 0, &features, backend, None)
            .unwrap_or_else(|cause| panic!("{backend:?}: {cause}"))
    };
    let rust = compile(GeometryBackend::Rust);
    let dir = tempfile::tempdir().unwrap();
    let near = |actual: f64, expected: f64| {
        assert!((actual - expected).abs() < 1e-6, "{actual} != {expected}");
    };
    for (index, feature) in features.into_iter().enumerate() {
        let actual = rust.solid(feature).unwrap();
        // Independent expectation from the saved sketch plane, profile
        // points, quarter arc and extrusion distance.
        let prism = saved_prism::extruded_prism(&document, feature);
        assert!(actual.uses_rust_geometry());
        assert_eq!(actual.solid_count().unwrap(), 1);
        let expected_volume = prism.volume_mm3();
        near(actual.volume_mm3().unwrap(), expected_volume);
        let a = actual.bounds().unwrap();
        let [min, max] = prism.bounds();
        for axis in 0..3 {
            near(a.min[axis], min[axis]);
            near(a.max[axis], max[axis]);
        }
        let topology = actual.topology_snapshot(feature, 0).unwrap();
        let arcs = topology
            .elements
            .iter()
            .filter_map(|element| {
                if let Some(AnalyticGeometry::CircularArc {
                    radius_mm,
                    sweep_rad,
                    ..
                }) = element.geometry
                {
                    assert_eq!(element.key.kind, TopologyKind::Edge);
                    Some((radius_mm, sweep_rad.abs()))
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();
        assert_eq!(arcs.len(), 2);
        for (radius, sweep) in arcs {
            near(radius, [4., 2.][index]);
            near(sweep, FRAC_PI_2);
        }
        assert_eq!(actual.faces().unwrap().len(), 7);
        for sign in [-1., 1.] {
            let options = DrawingProjection {
                frame: PlaneFrame {
                    v: [0., 0., sign],
                    ..PlaneFrame::default()
                },
                hidden_lines: true,
                ..DrawingProjection::default()
            };
            let project = |solid, options| {
                project_drawing(&[DrawingSolid { solid, hatch: None }], options).unwrap()
            };
            let vectors = project(actual, options);
            // Looking along the extrusion, the near and far profile caps
            // project onto the same five analytic edges.
            let expected_layer = projected_profile(&prism, options.frame.u, options.frame.v);
            for role in [VectorRole::Visible, VectorRole::Hidden] {
                let actual_layer: Vec<_> = vectors.iter().filter(|v| v.role == role).collect();
                assert_eq!(actual_layer.len(), expected_layer.len());
                for edge in &expected_layer {
                    assert_eq!(
                        actual_layer
                            .iter()
                            .filter(|target| matches_expected(target, edge))
                            .count(),
                        1,
                        "{role:?} analytic edge {:?} / {:?} must match one Rust edge",
                        edge.ends,
                        edge.arc
                    );
                }
            }
            let visible_only = project(
                actual,
                DrawingProjection {
                    hidden_lines: false,
                    ..options
                },
            );
            assert_eq!(
                visible_only.iter().collect::<Vec<_>>(),
                vectors
                    .iter()
                    .filter(|v| v.role == VectorRole::Visible)
                    .collect::<Vec<_>>()
            );
            let ellipse = visible_only
                .iter()
                .find_map(|vector| match vector.curve {
                    VectorCurve::Ellipse {
                        major_mm,
                        minor_mm,
                        start_rad,
                        end_rad,
                        ..
                    } => Some((major_mm, minor_mm, end_rad - start_rad)),
                    _ => None,
                })
                .unwrap();
            near(ellipse.0, [4., 2.][index]);
            near(ellipse.1, ellipse.0);
            near(ellipse.2, FRAC_PI_2);
        }
        let receipt = actual
            .manufacturing_mesh_with_receipt(MeshOptions::default())
            .unwrap();
        assert!(receipt.edge_pairs > 0 && receipt.vertex_links > 0);
        let path = dir.path().join(format!("floor-{index}.step"));
        actual.write_step(&path).unwrap();
        let policy = GeometryPolicy::default();
        for readback in [
            Solid::read_step(&path, policy).unwrap(),
            Solid::import_rust_step(&std::fs::read_to_string(&path).unwrap(), feature, policy)
                .unwrap(),
        ] {
            assert_eq!(readback.solid_count().unwrap(), 1);
            near(readback.volume_mm3().unwrap(), expected_volume);
        }
    }
}

#[test]
fn saved_toolcase_ribs_and_nut_pockets_rebuild_with_bounded_straight_sweeps() {
    let document = toolcase();
    let features = document
        .entities()
        .filter(|entity| {
            matches!(
                entity.name.as_str(),
                "Small drawer sliding rib extrusion" | "M4 captive nut pocket · 7.4 AF extrusion"
            )
        })
        .map(|entity| entity.id)
        .collect::<Vec<_>>();
    assert!(!features.is_empty());
    let rust = Compiler::default()
        .compile_sources_with_backend(&document, 0, &features, GeometryBackend::Rust, None)
        .unwrap_or_else(|cause| panic!("Rust: {cause}"));
    let [mut ribs, mut pockets] = [0, 0];
    for feature in features {
        let a = rust.solid(feature).unwrap();
        assert!(a.uses_rust_geometry());
        assert_eq!(a.solid_count().unwrap(), 1);
        // Independent expectation from the saved sketch plane, the saved
        // profile points and the authored extrusion distance.
        let prism = saved_prism::extruded_prism(&document, feature);
        if document.get(feature).unwrap().name == "Small drawer sliding rib extrusion" {
            assert_eq!(prism.loop_uv.len(), 4);
            ribs += 1;
        } else {
            assert_eq!(prism.loop_uv.len(), 6);
            pockets += 1;
        }
        assert!((a.volume_mm3().unwrap() - prism.volume_mm3()).abs() < 1e-6);
        let ([min, max], b) = (prism.bounds(), a.bounds().unwrap());
        for axis in 0..3 {
            assert!((b.min[axis] - min[axis]).abs() < 1e-6);
            assert!((b.max[axis] - max[axis]).abs() < 1e-6);
        }
    }
    assert_eq!([ribs, pockets], [16, 8]);
}
