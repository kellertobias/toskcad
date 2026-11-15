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
        let manufacturing = part.manufacturing_mesh_with_receipt(q).unwrap();
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
