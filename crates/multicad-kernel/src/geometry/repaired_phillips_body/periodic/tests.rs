use super::*;
use multicad_rust_kernel::{
    head_profile::{Head, Profile},
    helical_thread::{Spec as ThreadSpec, Thread},
    phillips_drive::{Cutter, Domain, Inner, Spec},
    profiled_screw::SourceSpec,
};
fn source(left: bool, plain: f64) -> Rc<OwnedPeriodicMaterial> {
    let head = Head::Crowned {
        radius_mm: 2.,
        height_mm: 1.6,
        crown_radius_mm: 3.2,
        top_land_radius_mm: 1.235645633015642,
        neck_radius_mm: 0.1,
    };
    let thread = Thread::new(ThreadSpec {
        major_radius_mm: 1.,
        pitch_mm: 0.4,
        length_mm: 11.5,
        depth_mm: 0.21650635094610965,
        crest_width_mm: 0.05,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: left,
        root_radius_mm: None,
    })
    .unwrap();
    let profile = Profile::new_for_thread_union(head, 1., plain).unwrap();
    let cutter = Cutter::new(
        Spec {
            mouth_diameter_mm: 1.9,
            mouth_radius_mm: 0.3,
            b_width_mm: 0.61,
            g_diameter_mm: 0.81,
            tilt_angle_deg: 7.,
            outer_flank_angle_deg: 26.5,
            floor_angle_deg: 28.,
            outer_corner_angle_deg: 92.,
            inner: Inner::Rounded { radius_mm: 0.25 },
        },
        Domain {
            corner_start_depth_mm: -3.8,
            corner_end_depth_mm: 3.208443210824065,
            corner_reach_mm: 16.,
        },
    )
    .unwrap();
    Rc::new(
        OwnedPeriodicMaterial::from_authored_parts(
            &thread,
            &profile,
            &cutter,
            SourceSpec {
                head: Some(head),
                shank_length_mm: plain,
                minimum_wall_mm: 1e-5,
            },
            Default::default(),
            None,
        )
        .unwrap(),
    )
}
fn feature() -> EntityId {
    EntityId::from_uuid_bits(0x00ca_d352_0000_0000_0000_0000_0000_0021)
}
#[test]
fn actual_long_m2_periodic_display_and_fresh_manufacturing_have_literal_embedding() {
    for left in [false, true] {
        let start = std::time::Instant::now();
        let source = source(left, 0.5);
        let part = Solid::rust_periodic_repaired_phillips_material(
            Rc::clone(&source),
            feature(),
            PlaneFrame::default(),
            Default::default(),
            None,
        )
        .unwrap();
        let create = start.elapsed();
        let start = std::time::Instant::now();
        let q = MeshOptions {
            linear_deflection_mm: 0.05,
            angular_deflection_rad: 0.3,
            maximum_triangles: 1_000_000,
        };
        let (first, receipt) = part
            .rust_periodic_repaired_phillips_display_mesh_with_receipt(q, Default::default(), None)
            .unwrap();
        assert!(std::ptr::eq(receipt.source(), source.as_ref()));
        assert_eq!(receipt.occurrence(), PlaneFrame::default());
        assert_eq!(receipt.quality(), q);
        let (generous, generous_receipt) = part
            .rust_periodic_repaired_phillips_display_mesh_with_receipt(
                MeshOptions {
                    maximum_triangles: 2_000_000,
                    ..q
                },
                Default::default(),
                None,
            )
            .unwrap();
        assert_eq!(generous_receipt.quality(), q);
        assert!(generous_receipt.reused_display_receipt());
        assert_eq!(generous.vertices, first.vertices);
        assert_eq!(generous.indices, first.indices);
        assert!(generous.indices.len() / 3 <= 1_000_000);
        assert!(receipt.embedding_pairs() > 0 && !receipt.reused_display_receipt());
        assert!(
            receipt.local_source_surface_error_mm() > 0.
                && receipt.local_source_surface_error_mm() <= q.linear_deflection_mm
        );
        assert_eq!(
            receipt.new_source_volume_mm3(),
            part.rust_repaired_phillips_volume_bounds().unwrap()
        );
        assert!(receipt.literal_signed_volume_mm3() > 0.);
        let literal = receipt.literal_volume_enclosure_mm3();
        assert!(literal[0] > 0. && literal[1] >= literal[0]);
        assert!(
            receipt.source_volume_error_bound_mm3() > 0.
                && receipt.source_volume_error_bound_mm3() < 5.
        );
        assert_eq!(receipt.query_work_used() + receipt.remaining_work(), WORK);
        assert!(
            receipt.source_work_used() > 2_000_000
                && receipt.source_work_used() < receipt.query_work_used()
        );
        assert!(
            receipt.combined_live_upper_bound_bytes()
                >= receipt.consumer_accounted_peak_bytes()
                    + receipt.producer_combined_upper_bound_bytes()
                    + receipt.validator_workspace_upper_bound_bytes()
        );
        let cold = start.elapsed();
        let start = std::time::Instant::now();
        let (repeat, warm_receipt) = part
            .rust_periodic_repaired_phillips_display_mesh_with_receipt(q, Default::default(), None)
            .unwrap();
        assert_eq!(repeat, first);
        assert!(warm_receipt.reused_display_receipt());
        let warm = start.elapsed();
        let start = std::time::Instant::now();
        let manufacturing = part
            .manufacturing_mesh_with_receipt(MeshOptions {
                maximum_triangles: 2_000_000,
                ..q
            })
            .unwrap();
        assert_eq!(manufacturing.mesh, first);
        assert!(manufacturing.edge_pairs > 0 && manufacturing.vertex_links > 0);
        assert!(first.indices.len() > 1000);
        assert!(
            first
                .face_indices
                .iter()
                .all(|&i| (i as usize) < source.face_count())
        );
        eprintln!(
            "actual longM2 periodic left{left}: create{create:?} cold{cold:?} warm{warm:?} manufacturing{:?} triangles{}",
            start.elapsed(),
            first.indices.len() / 3
        );
    }
}
#[test]
fn actual_periodic_native_mesh_limits_cancellation_and_pose_preserve_source_and_cache() {
    let source = source(false, 0.03);
    let part = Solid::rust_periodic_repaired_phillips_material(
        Rc::clone(&source),
        feature(),
        PlaneFrame::default(),
        Default::default(),
        None,
    )
    .unwrap();
    let q = MeshOptions {
        linear_deflection_mm: 0.1,
        angular_deflection_rad: 0.4,
        maximum_triangles: 100_000,
    };
    for budget in [
        RepairedPhillipsMeshBudget {
            maximum_work: 0,
            ..Default::default()
        },
        RepairedPhillipsMeshBudget {
            maximum_work: WORK + 1,
            ..Default::default()
        },
        RepairedPhillipsMeshBudget {
            maximum_bytes: 0,
            ..Default::default()
        },
        RepairedPhillipsMeshBudget {
            maximum_bytes: BYTES + 1,
            ..Default::default()
        },
    ] {
        assert!(
            part.rust_periodic_repaired_phillips_mesh_with_receipt(q, budget, None)
                .is_err()
        );
    }
    for bad in [
        MeshOptions {
            linear_deflection_mm: 0.,
            ..q
        },
        MeshOptions {
            angular_deflection_rad: 0.,
            ..q
        },
        MeshOptions {
            maximum_triangles: 0,
            ..q
        },
    ] {
        assert_eq!(
            part.rust_periodic_repaired_phillips_mesh_with_receipt(bad, Default::default(), None)
                .err()
                .unwrap()
                .kind,
            GeometryErrorKind::Input
        );
    }
    let (first, _) = part
        .rust_periodic_repaired_phillips_display_mesh_with_receipt(q, Default::default(), None)
        .unwrap();
    for at in [0, 1, 10, 100, 1000] {
        let calls = std::cell::Cell::new(0);
        let cancel = || {
            let n = calls.get();
            calls.set(n + 1);
            n == at
        };
        assert_eq!(
            part.rust_periodic_repaired_phillips_mesh_with_receipt(
                q,
                Default::default(),
                Some(&cancel)
            )
            .err()
            .unwrap()
            .kind,
            GeometryErrorKind::Cancelled
        );
    }
    let (same, receipt) = part
        .rust_periodic_repaired_phillips_display_mesh_with_receipt(q, Default::default(), None)
        .unwrap();
    assert!(receipt.reused_display_receipt());
    assert_eq!(same, first);
    let (_, fresh) = part
        .rust_periodic_repaired_phillips_mesh_with_receipt(q, Default::default(), None)
        .unwrap();
    let used = fresh.query_work_used();
    assert!(
        part.rust_periodic_repaired_phillips_mesh_with_receipt(
            q,
            RepairedPhillipsMeshBudget {
                maximum_work: used,
                ..Default::default()
            },
            None
        )
        .is_ok()
    );
    assert_eq!(
        part.rust_periodic_repaired_phillips_mesh_with_receipt(
            q,
            RepairedPhillipsMeshBudget {
                maximum_work: used - 1,
                ..Default::default()
            },
            None
        )
        .err()
        .unwrap()
        .kind,
        GeometryErrorKind::Limit
    );
    let cube = Solid::box_solid([1.; 3], Default::default()).unwrap();
    assert!(
        cube.rust_periodic_repaired_phillips_mesh_with_receipt(q, Default::default(), None)
            .is_err()
    );
    assert!(
        part.rust_repaired_phillips_mesh_with_receipt(q, Default::default(), None)
            .is_err()
    );
    let frame = PlaneFrame {
        origin: [10., 20., 30.],
        u: [0.8, 0.6, 0.],
        v: [0.36, -0.48, -0.8],
    };
    let posed = Solid::rust_periodic_repaired_phillips_material(
        source,
        feature(),
        frame,
        Default::default(),
        None,
    )
    .unwrap();
    let (mesh, receipt) = posed
        .rust_periodic_repaired_phillips_mesh_with_receipt(q, Default::default(), None)
        .unwrap();
    assert_eq!(receipt.occurrence(), frame);
    assert!(!mesh.indices.is_empty() && receipt.embedding_pairs() > 0);
    assert!(receipt.source_volume_error_bound_mm3() > 0.);
}
fn world(source: &OwnedPeriodicMaterial, frame: PlaneFrame, p: Point) -> [[f64; 2]; 3] {
    let c = source.chart();
    let z = I::point(p[2]) - I::point(c.outer_top_z);
    let built: [I; 3] = std::array::from_fn(|i| {
        I::point(c.columns[0][i])
            + I::point(c.columns[1][i]) * I::point(p[0])
            + I::point(c.columns[2][i]) * I::point(p[1])
            + I::point(c.columns[3][i]) * z
    });
    let n = frame.normal();
    std::array::from_fn(|i| {
        let v = I::point(frame.origin[i])
            + I::point(frame.u[i]) * built[0]
            + I::point(frame.v[i]) * built[1]
            + I::point(n[i]) * built[2];
        [v.lo, v.hi]
    })
}
#[test]
fn actual_authored_periodic_native_ownership_mass_and_whole_box_queries() {
    for (left, plain, frame) in [
        (false, 0.5, PlaneFrame::default()),
        (
            true,
            0.03,
            PlaneFrame {
                origin: [10., 20., 30.],
                u: [0.8, 0.6, 0.],
                v: [0.36, -0.48, -0.8],
            },
        ),
    ] {
        let source = source(left, plain);
        let part = Solid::rust_periodic_repaired_phillips_material(
            Rc::clone(&source),
            feature(),
            frame,
            Default::default(),
            None,
        )
        .unwrap();
        assert!(part.rust_repaired_phillips_source().is_none());
        let (actual, occurrence) = part.rust_periodic_repaired_phillips_source().unwrap();
        assert!(std::ptr::eq(actual, source.as_ref()));
        assert_eq!(occurrence, frame);
        let mass = part.rust_repaired_phillips_volume_bounds().unwrap();
        assert!(mass[0] > 0. && mass[1] >= mass[0]);
        for (p, expected) in [
            ([0., 0., 5.], Classification::Inside),
            ([3., 0., 5.], Classification::Outside),
            ([0., 0., -1.6 - plain + 0.2], Classification::Outside),
            ([1.3, 0., -0.5 - plain], Classification::Inside),
        ] {
            let bounds = world(&source, frame, p);
            let query = part
                .rust_periodic_repaired_phillips_classify_world_box(
                    bounds,
                    Default::default(),
                    None,
                )
                .unwrap();
            assert_eq!(query.location(), expected);
            assert!(std::ptr::eq(query.source(), source.as_ref()));
            assert_eq!(query.world_bounds_mm(), bounds);
            assert_eq!(query.occurrence(), frame);
            for (b, x) in query.local_bounds_mm().into_iter().zip(p) {
                assert!(b[0] <= x && x <= b[1]);
            }
            assert_eq!(query.work_used() + query.remaining_work(), MAX);
            let point = bounds.map(|v| v[0].midpoint(v[1]));
            assert_eq!(
                part.rust_thread_contains(point).unwrap(),
                expected == Classification::Inside
            );
        }
    }
}

#[test]
fn actual_periodic_source_topology_retains_two_sided_finite_incidence() {
    use multicad_core::references::TopologyKind;
    for left in [false, true] {
        let source = source(left, 0.5);
        let part = Solid::rust_periodic_repaired_phillips_material(
            Rc::clone(&source),
            feature(),
            PlaneFrame::default(),
            Default::default(),
            None,
        )
        .unwrap();
        let snapshot = part.topology_snapshot(feature(), 7).unwrap();
        snapshot.validate().unwrap();
        assert_eq!(snapshot.revision, 7);
        let faces: Vec<_> = snapshot
            .elements
            .iter()
            .filter(|e| e.key.kind == TopologyKind::Face)
            .collect();
        assert_eq!(faces.len(), source.face_count());
        let edges: Vec<_> = snapshot
            .elements
            .iter()
            .filter(|e| e.key.kind == TopologyKind::Edge)
            .collect();
        for edge in &edges {
            assert_eq!(edge.boundary.len(), 2);
            for origin in &edge.origins {
                assert_eq!(
                    faces.iter().filter(|f| f.boundary.contains(origin)).count(),
                    2
                );
            }
        }
        assert!(
            edges
                .iter()
                .any(|e| e.origins.iter().any(|o| o.role == "thread/cap/0/rim/3")
                    && e.boundary.iter().any(|o| o.role.starts_with("screw.head/")))
        );
        assert!(!edges.iter().any(|e| {
            e.origins
                .iter()
                .any(|o| o.role == "screw.head/thread.junction/3")
        }));
        for vertex in snapshot
            .elements
            .iter()
            .filter(|e| e.key.kind == TopologyKind::Vertex)
        {
            assert!(!vertex.boundary.is_empty());
            assert!(vertex.geometry.is_some());
        }
    }
}

#[test]
fn actual_source_planar_faces_keep_bounded_frames_through_occurrence() {
    use multicad_core::references::{AnalyticGeometry, TopologyKind};
    let source = source(true, 0.5);
    let occurrence = PlaneFrame {
        origin: [10., -3., 4.],
        u: [0., 1., 0.],
        v: [0., 0., 1.],
    };
    let part = Solid::rust_periodic_repaired_phillips_material(
        Rc::clone(&source),
        feature(),
        occurrence,
        Default::default(),
        None,
    )
    .unwrap();
    let topology = part.topology_snapshot(feature(), 1).unwrap();
    let tip = topology
        .elements
        .iter()
        .find(|e| {
            e.key.kind == TopologyKind::Face
                && e.origins.iter().any(|o| o.role == "screw.thread/end")
        })
        .unwrap();
    let Some(AnalyticGeometry::Plane { frame }) = tip.geometry else {
        panic!("actual tip plane")
    };
    frame.validate().unwrap();
    assert!((frame.origin_mm[0] - 21.5).abs() < 1e-9);
    assert!((frame.origin_mm[1] + 3.).abs() < 1e-9);
    assert!((frame.origin_mm[2] - 4.).abs() < 1e-9);
    assert_eq!(tip.frames.len(), 1);
    for face in topology
        .elements
        .iter()
        .filter(|e| e.key.kind == TopologyKind::Face)
    {
        let curved = face
            .origins
            .iter()
            .any(|o| o.role.starts_with("thread/start/") || o.role == "screw.head/crown");
        if curved {
            assert!(face.geometry.is_none());
        }
    }
}

#[test]
fn actual_periodic_faces_are_source_measured_lazy_and_cached() {
    for left in [false, true] {
        let source = source(left, 0.5);
        let part = Solid::rust_periodic_repaired_phillips_material(
            Rc::clone(&source),
            feature(),
            PlaneFrame::default(),
            Default::default(),
            None,
        )
        .unwrap();
        let (body, _) = part.repaired_phillips_body().unwrap();
        assert!(body.periodic_properties.borrow().is_none());
        let before = body.resident_bytes();
        let now = std::time::Instant::now();
        let faces = part.faces().unwrap();
        let cold = now.elapsed();
        assert_eq!(faces.len(), source.face_count());
        assert!(
            faces
                .iter()
                .all(|f| f.area_mm2 > 0. && f.center_mm.into_iter().all(f64::is_finite))
        );
        assert!(body.resident_bytes() > before);
        let now = std::time::Instant::now();
        assert_eq!(faces, part.faces().unwrap());
        eprintln!(
            "Native actual periodic faces LH{left}: cold{cold:?} warm{:?}",
            now.elapsed()
        );
    }
}

#[test]
fn actual_periodic_face_centroid_transport_encloses_nonconstant_affine_weights() {
    let source = source(false, 0.5);
    let occurrence = PlaneFrame {
        origin: [10., -3., 4.],
        u: [1. + 1e-12, 0., 0.],
        v: [0., 1., 0.],
    };
    let part = Solid::rust_periodic_repaired_phillips_material(
        Rc::clone(&source),
        feature(),
        occurrence,
        Default::default(),
        None,
    )
    .unwrap();
    let faces = part.faces().unwrap();
    let tip = &faces[source.head().view().segment_count() + 4];
    assert!(tip.plane.is_some());
    assert!((tip.center_mm[2] - (4. + 11.5 * (1. + 1e-12))).abs() < 1e-8);
    assert!(faces.iter().all(|f| f.area_mm2 > 0.));
}

#[test]
fn actual_periodic_edge_display_selects_every_finite_owner_and_sews_original_vertices() {
    use multicad_core::references::{AnalyticGeometry, TopologyKind};
    for left in [false, true] {
        let source = source(left, 0.5);
        let occurrence = PlaneFrame {
            origin: [10., -3., 4.],
            u: [0., 1., 0.],
            v: [0., 0., 1.],
        };
        let part = Solid::rust_periodic_repaired_phillips_material(
            source,
            feature(),
            occurrence,
            Default::default(),
            None,
        )
        .unwrap();
        let topology = part.topology_snapshot(feature(), 1).unwrap();
        let now = std::time::Instant::now();
        let edges = part.display_edges(0.05).unwrap();
        // The Model compiler uses this combined path; direct edge display
        // alone must not mask a missing body dispatch here.
        let combined = part
            .topology_with_display(EntityId::new(), 13, Some(0.05))
            .unwrap();
        assert_eq!(combined.edges, edges);
        assert_eq!(combined.vertices, part.display_vertices().unwrap());

        eprintln!(
            "PERIODIC source edges LH{left}: {}curves {}points {:?}",
            edges.len(),
            edges.iter().map(|e| e.points_mm.len()).sum::<usize>(),
            now.elapsed()
        );
        let source_edges: Vec<_> = topology
            .elements
            .iter()
            .filter(|e| e.key.kind == TopologyKind::Edge)
            .collect();
        assert_eq!(edges.len(), source_edges.len());
        // Independent dense source queries check the observable polyline
        // accuracy, rather than reusing the producer's chord-bound calculation.
        let (owner, _) = part.rust_periodic_repaired_phillips_source().unwrap();
        let mut exact = owner.source_curve_query(16_000_000, None).unwrap();
        let distance = |p: Point, polyline: &[Point]| {
            polyline
                .windows(2)
                .map(|v| {
                    let d: Point = std::array::from_fn(|i| v[1][i] - v[0][i]);
                    let norm = d.into_iter().map(|x| x * x).sum::<f64>();
                    let t = if norm == 0. {
                        0.
                    } else {
                        (0..3).map(|i| (p[i] - v[0][i]) * d[i]).sum::<f64>() / norm
                    }
                    .clamp(0., 1.);
                    (0..3)
                        .map(|i| (p[i] - v[0][i] - t * d[i]).powi(2))
                        .sum::<f64>()
                        .sqrt()
                })
                .fold(f64::INFINITY, f64::min)
        };
        for band in 0..4 {
            for kind in 0..3 {
                let role = if kind == 0 {
                    format!("thread/rail/{band}")
                } else {
                    format!("thread/cap/{}/rim/{band}", kind - 1)
                };
                let edge = source_edges
                    .iter()
                    .find(|e| e.origins.iter().any(|o| o.role == role))
                    .unwrap();
                let polyline = &edges[edge.key.ordinal as usize].points_mm;
                for i in 0..=64 {
                    let t = i as f64 / 64.;
                    let p = if kind == 0 {
                        exact.thread_rail_point_mm(band, t)
                    } else {
                        exact.thread_rim_point_mm(kind - 1, band, t)
                    }
                    .unwrap();
                    let p = p.map(|v| v[0].midpoint(v[1]));
                    let world = std::array::from_fn(|i| {
                        occurrence.origin[i]
                            + occurrence.u[i] * p[0]
                            + occurrence.v[i] * p[1]
                            + occurrence.normal()[i] * p[2]
                    });
                    assert!(
                        distance(world, polyline) <= 0.05 + 1e-10,
                        "source curve {role} parameter{t}"
                    );
                }
            }
        }

        for (i, edge) in edges.iter().enumerate() {
            assert_eq!(edge.index, i as u32);
            assert!(
                edge.points_mm.len() >= 2 && edge.points_mm.iter().flatten().all(|p| p.is_finite())
            );
        }
        for vertex in topology
            .elements
            .iter()
            .filter(|e| e.key.kind == TopologyKind::Vertex)
        {
            let Some(AnalyticGeometry::Point { position_mm }) = vertex.geometry else {
                panic!("actual source vertex")
            };
            for origin in &vertex.boundary {
                let edge = source_edges
                    .iter()
                    .find(|e| e.origins.contains(origin))
                    .unwrap();
                let rendered = &edges[edge.key.ordinal as usize].points_mm;
                assert!(
                    rendered.first() == Some(&position_mm) || rendered.last() == Some(&position_mm),
                    "original endpoint {origin:?}"
                );
            }
        }
        for edge in source_edges.iter().filter(|e| {
            e.origins
                .iter()
                .any(|o| o.role.starts_with("screw.head/") && o.role.ends_with("/start"))
        }) {
            let p = &edges[edge.key.ordinal as usize].points_mm;
            assert_eq!(p.first(), p.last());
        }
        assert_eq!(
            part.display_vertices().unwrap().len(),
            topology
                .elements
                .iter()
                .filter(|e| e.key.kind == TopologyKind::Vertex)
                .count()
        );
        assert!(!part.supports_step_export());
    }
}
