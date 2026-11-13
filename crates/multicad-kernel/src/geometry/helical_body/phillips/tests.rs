use super::*;
use multicad_rust_kernel::{
    head_profile::{Head, Profile},
    helical_thread,
    profiled_screw::SourceSpec,
};
pub(super) fn fixture(index: usize, plain: f64, left: bool) -> PhillipsScrew {
    let (m, r, b, g, inner, beta, rs, pitch, head) = match index {
        0 => (
            1.9,
            0.3,
            0.61,
            0.81,
            source::Inner::Rounded { radius_mm: 0.25 },
            7.,
            1.,
            0.4,
            Head::Crowned {
                radius_mm: 2.,
                height_mm: 1.6,
                crown_radius_mm: 3.2,
                top_land_radius_mm: 1.25,
                neck_radius_mm: 0.1,
            },
        ),
        1 => (
            3.,
            0.5,
            0.955,
            1.27,
            source::Inner::ShallowV {
                angle_deg: 138.,
                e_width_mm: 0.435,
            },
            7.,
            1.5,
            0.5,
            Head::Crowned {
                radius_mm: 2.8,
                height_mm: 2.4,
                crown_radius_mm: 3.,
                top_land_radius_mm: 2.1,
                neck_radius_mm: 0.1,
            },
        ),
        2 => (
            3.2,
            0.5,
            0.955,
            1.27,
            source::Inner::ShallowV {
                angle_deg: 138.,
                e_width_mm: 0.435,
            },
            7.,
            1.5,
            0.5,
            Head::Countersunk {
                radius_mm: 2.75,
                height_mm: 1.65,
                included_angle_deg: 90.,
                edge_chamfer_mm: 0.05,
                neck_radius_mm: 0.1,
            },
        ),
        3 => (
            4.4,
            0.6,
            1.455,
            2.29,
            source::Inner::ShallowV {
                angle_deg: 140.,
                e_width_mm: 0.815,
            },
            5.75,
            2.,
            0.7,
            Head::Crowned {
                radius_mm: 4.,
                height_mm: 3.1,
                crown_radius_mm: 6.5,
                top_land_radius_mm: 3.,
                neck_radius_mm: 0.2,
            },
        ),
        4 => (
            6.9,
            0.8,
            2.395,
            3.81,
            source::Inner::ShallowV {
                angle_deg: 146.,
                e_width_mm: 2.005,
            },
            5.75,
            3.,
            1.,
            Head::Crowned {
                radius_mm: 6.,
                height_mm: 4.6,
                crown_radius_mm: 10.,
                top_land_radius_mm: 4.2,
                neck_radius_mm: 0.25,
            },
        ),
        _ => (
            9.,
            1.,
            3.465,
            5.08,
            source::Inner::ShallowV {
                angle_deg: 153.,
                e_width_mm: 2.415,
            },
            7.,
            4.,
            1.25,
            Head::Crowned {
                radius_mm: 8.,
                height_mm: 6.,
                crown_radius_mm: 13.,
                top_land_radius_mm: 5.4,
                neck_radius_mm: 0.4,
            },
        ),
    };
    let spec = source::Spec {
        mouth_diameter_mm: m,
        mouth_radius_mm: r,
        b_width_mm: b,
        g_diameter_mm: g,
        tilt_angle_deg: beta,
        outer_flank_angle_deg: 26.5,
        floor_angle_deg: 28.,
        outer_corner_angle_deg: 92.,
        inner,
    };
    let tip = (m - g) / (2. * 26.5_f64.to_radians().tan()) + g / 2. * 28_f64.to_radians().tan();
    let cavity = source::Cavity::new(
        source::Cutter::new(
            spec,
            source::Domain {
                corner_start_depth_mm: -2. * m,
                corner_end_depth_mm: tip + m,
                corner_reach_mm: 100.,
            },
        )
        .unwrap(),
        1e-5,
        100_000,
    )
    .unwrap();
    let thread = Thread::new(helical_thread::Spec {
        major_radius_mm: rs - 0.03,
        pitch_mm: pitch,
        length_mm: 8.,
        depth_mm: 0.6 * pitch,
        crest_width_mm: 0.125 * pitch,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: left,
        root_radius_mm: None,
    })
    .unwrap();
    PhillipsScrew::new(
        thread,
        Profile::new_for_thread_union(head, rs, plain).unwrap(),
        SourceSpec {
            head: Some(head),
            shank_length_mm: plain,
            minimum_wall_mm: 1e-5,
        },
        cavity,
        100_000,
        &|| false,
    )
    .unwrap()
}
fn feature() -> EntityId {
    EntityId::from_uuid_bits(0xa129352)
}
#[test]
fn typed_source_properties_and_original_geometry_survive_native_cardinal_and_proper_poses() {
    for (index, plain, frame) in [
        (0, 0., PlaneFrame::default()),
        (
            1,
            1.,
            PlaneFrame {
                origin: [10., 20., 30.],
                u: [0.8, 0.6, 0.],
                v: [-0.36, 0.48, 0.8],
            },
        ),
        (2, 0., PlaneFrame::default()),
        (5, 1., PlaneFrame::default()),
    ] {
        let source = fixture(index, plain, index == 1);
        let spec = source.thread().spec();
        let profile = source.profile().segments().collect::<Vec<_>>();
        let cavity = source.cavity().source().spec();
        let s = Solid::rust_phillips_source_for_test(
            source,
            feature(),
            frame,
            GeometryPolicy::default(),
            &|| false,
        )
        .unwrap();
        let (retained, actual_frame) = s.rust_phillips_screw_geometry().unwrap();
        assert_eq!(actual_frame, frame);
        assert_eq!(retained.thread().spec(), spec);
        assert_eq!(retained.profile().segments().collect::<Vec<_>>(), profile);
        assert_eq!(retained.cavity().source().spec(), cavity);
        let faces = s.faces().unwrap();
        assert_eq!(faces.len(), retained.faces().len());
        assert!(
            faces
                .iter()
                .all(|f| f.area_mm2 > 0. && f.center_mm.iter().all(|v| v.is_finite()))
        );
        let topology = s.topology_snapshot(feature(), 1).unwrap();
        topology.validate().unwrap();
        assert_eq!(
            topology
                .elements
                .iter()
                .filter(|e| e.key.kind == TopologyKind::Face)
                .count(),
            faces.len()
        );
        assert!(!s.display_edges(0.1).unwrap().is_empty());
        assert!(!s.display_vertices().unwrap().is_empty());
        assert!(s.volume_mm3().unwrap() > 0.);
        assert!(
            s.rust_thread_contains(position(frame, [0., 0., 4.]))
                .unwrap()
        );
        assert!(
            !s.rust_thread_contains(position(frame, [20., 0., 0.]))
                .unwrap()
        );
    }
}
#[test]
fn actual_native_meshes_keep_source_volume_and_literal_placement_receipts_without_step_admission() {
    for frame in [
        PlaneFrame::default(),
        PlaneFrame {
            origin: [10., 20., 30.],
            u: [0.8, 0.6, 0.],
            v: [-0.36, 0.48, 0.8],
        },
    ] {
        let source = fixture(1, 0., false);
        let original = source.source_volume_bounds().unwrap();
        let s = Solid::rust_phillips_source_for_test(
            source,
            feature(),
            frame,
            GeometryPolicy::default(),
            &|| false,
        )
        .unwrap();
        assert!(s.rust_phillips_mesh_volume_pair().is_none());
        let receipt = s
            .manufacturing_mesh_with_receipt(Default::default())
            .unwrap();
        assert_eq!(receipt.edge_pairs as usize * 2, receipt.mesh.indices.len());
        let pair = s.rust_phillips_mesh_volume_pair().unwrap();
        assert!(pair.source_world_mm3[0] > 0. && pair.literal_world_mm3[0] > 0.);
        assert!(
            (pair.literal_world_mm3[0] - pair.source_world_mm3[1]).abs()
                <= pair.source_to_literal_bound_mm3 + original.volume_mm3()[1]
                    - original.volume_mm3()[0]
        );
        assert!(pair.literal_vertex_transfer_bound_mm3 >= 0.);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("not-qualified.step");
        std::fs::write(&path, b"keep destination").unwrap();
        assert!(s.write_step(&path).is_err());
        assert_eq!(std::fs::read(path).unwrap(), b"keep destination");
    }
}

#[test]
fn failed_source_and_mesh_proofs_publish_no_partial_state_or_stale_pose_receipt() {
    use std::cell::Cell;
    let source = fixture(1, 1., false);
    let spec = source.thread().spec();
    let policy = GeometryPolicy::default();
    let fail = Solid::rust_phillips_source_for_test(
        source.clone(),
        feature(),
        PlaneFrame::default(),
        policy,
        &|| true,
    )
    .err()
    .unwrap();
    assert_eq!(fail.kind, GeometryErrorKind::Cancelled);
    let calls = Cell::new(0);
    let fail = Solid::rust_phillips_source_for_test(
        source.clone(),
        feature(),
        PlaneFrame::default(),
        policy,
        &|| {
            calls.set(calls.get() + 1);
            calls.get() > 20
        },
    )
    .err()
    .unwrap();
    assert_eq!(fail.kind, GeometryErrorKind::Cancelled);
    for policy in [
        GeometryPolicy {
            linear_tolerance_mm: 0.,
            ..policy
        },
        GeometryPolicy {
            minimum_feature_mm: 1.,
            ..policy
        },
        GeometryPolicy {
            linear_tolerance_mm: 1e-14,
            ..policy
        },
    ] {
        assert!(
            Solid::rust_phillips_source_for_test(
                source.clone(),
                feature(),
                PlaneFrame::default(),
                policy,
                &|| false
            )
            .is_err()
        );
    }
    assert!(
        Solid::rust_phillips_source_for_test(
            source.clone(),
            feature(),
            PlaneFrame {
                u: [1.01, 0., 0.],
                ..Default::default()
            },
            policy,
            &|| false
        )
        .is_err()
    );
    let solid = Solid::rust_phillips_source_for_test(
        source.clone(),
        feature(),
        PlaneFrame::default(),
        policy,
        &|| false,
    )
    .unwrap();
    assert!(
        solid
            .mesh(MeshOptions {
                maximum_triangles: 1,
                ..Default::default()
            })
            .is_err()
    );
    assert!(solid.rust_phillips_mesh_volume_pair().is_none());
    solid.mesh(Default::default()).unwrap();
    let pair = solid.rust_phillips_mesh_volume_pair().unwrap();
    assert_eq!(pair.linear_deflection_mm, 0.1);
    assert_eq!(pair.angular_deflection_rad, 0.3);
    assert_eq!(pair.triangle_limit, 1_000_000);
    assert!(
        solid
            .mesh(MeshOptions {
                maximum_triangles: 1,
                ..Default::default()
            })
            .is_err()
    );
    assert_eq!(
        solid
            .rust_phillips_mesh_volume_pair()
            .unwrap()
            .literal_world_mm3,
        pair.literal_world_mm3
    );
    let moved = solid
        .transformed(RigidTransform {
            translation_mm: [1., 2., 3.],
            ..Default::default()
        })
        .unwrap();
    assert!(moved.rust_phillips_mesh_volume_pair().is_none());
    assert!(moved.faces().is_ok());
    let far = solid
        .transformed(RigidTransform {
            translation_mm: [1e6, -2e6, 3e6],
            ..Default::default()
        })
        .unwrap();
    assert!(far.faces().is_ok());
    assert!(far.volume_mm3().unwrap() > 0.);
    let (body, frame) = solid.helical_body().unwrap();
    assert!(body.volume_mm3() > 0.);
    assert!(
        body.volume_error_bound(Solid::rust_mesh_options(Default::default()))
            .is_infinite()
    );
    let actual = solid.rust_phillips_screw_geometry().unwrap().0;
    for point in [
        [0., 0., 4.],
        [0., 0., actual.head_top_z() + 0.2],
        [1.4, 0., 4.],
    ] {
        assert_eq!(body.contains(point, 0.), actual.contains(point, 0.));
    }
    assert_eq!(frame, PlaneFrame::default());
    assert_eq!(source.thread().spec(), spec);
}

#[test]
fn admitted_near_proper_columns_keep_actual_determinant_and_source_geometry() {
    let source = fixture(0, 1., false);
    let volume = source.source_volume_bounds().unwrap();
    let frame = PlaneFrame {
        u: [1. + 2e-11, 0., 0.],
        ..Default::default()
    };
    let solid = Solid::rust_phillips_source_for_test(
        source,
        feature(),
        frame,
        GeometryPolicy::default(),
        &|| false,
    )
    .unwrap();
    assert!(solid.faces().is_ok());
    assert!(solid.volume_mm3().unwrap() > volume.volume_mm3()[0].midpoint(volume.volume_mm3()[1]));
    solid
        .manufacturing_mesh_with_receipt(Default::default())
        .unwrap();
    let pair = solid.rust_phillips_mesh_volume_pair().unwrap();
    assert!(pair.source_world_mm3[0] > volume.volume_mm3()[0]);
    assert_eq!(pair.frame, frame);
    assert_eq!(solid.rust_phillips_screw_geometry().unwrap().1, frame);
}

#[test]
fn every_metric_cavity_keeps_default_quality_in_cardinal_and_proper_native_meshes() {
    for index in 0..6 {
        for frame in [
            PlaneFrame::default(),
            PlaneFrame {
                origin: [10., 20., 30.],
                u: [0.8, 0.6, 0.],
                v: [-0.36, 0.48, 0.8],
            },
        ] {
            let source = fixture(index, 1., index % 2 == 1);
            let solid = Solid::rust_phillips_source_for_test(
                source,
                feature(),
                frame,
                GeometryPolicy::default(),
                &|| false,
            )
            .unwrap_or_else(|e| panic!("H{index} source {frame:?}: {e:?}"));
            let receipt = solid
                .manufacturing_mesh_with_receipt(Default::default())
                .unwrap_or_else(|e| panic!("H{index} default MFG {frame:?}: {e:?}"));
            assert_eq!(receipt.edge_pairs as usize * 2, receipt.mesh.indices.len());
            let source = solid.rust_phillips_screw_geometry().unwrap().0;
            assert!(
                receipt
                    .mesh
                    .face_indices
                    .iter()
                    .all(|&i| (i as usize) < source.faces().len())
            );
            let pair = solid.rust_phillips_mesh_volume_pair().unwrap();
            assert_eq!(pair.angular_deflection_rad, 0.3);
            assert_eq!(pair.linear_deflection_mm, 0.1);
            assert_eq!(pair.frame, frame);
        }
    }
}

#[test]
fn transported_sources_and_wrong_native_getters_remain_explicitly_separate() {
    let source = fixture(1, 1., false);
    let segments = source.profile().segments().collect::<Vec<_>>();
    let numeric = segments
        .iter()
        .flat_map(|s| [s.start, s.end])
        .flatten()
        .map(|x| (x.next_up() - x).abs())
        .fold(0_f64, f64::max);
    let transported = Profile::from_transported_segments(
        &segments,
        multicad_rust_kernel::head_profile::TransportError {
            numeric_conversion_mm: numeric,
            arithmetic_mm: 0.,
        },
    )
    .unwrap();
    let transported = PhillipsScrew::new(
        *source.thread(),
        transported,
        source.spec(),
        source.cavity().clone(),
        100_000,
        &|| false,
    )
    .unwrap();
    assert_eq!(
        Solid::rust_phillips_source_for_test(
            transported,
            feature(),
            PlaneFrame::default(),
            Default::default(),
            &|| false
        )
        .err()
        .unwrap()
        .kind,
        GeometryErrorKind::Unsupported
    );
    let bare = Solid::rust_helical_thread(
        imported_request(source.thread().spec()),
        None,
        feature(),
        Default::default(),
        Default::default(),
    )
    .unwrap();
    assert!(bare.rust_phillips_screw_geometry().is_none());
    assert!(bare.rust_phillips_mesh_volume_pair().is_none());
    assert!(bare.faces().is_ok());
    assert!(bare.volume_mm3().unwrap() > 0.);
    for (cause, kind) in [
        (source::Error::Invalid, GeometryErrorKind::Input),
        (source::Error::Unresolved, GeometryErrorKind::Unsupported),
        (source::Error::Limit, GeometryErrorKind::Limit),
        (source::Error::Cancelled, GeometryErrorKind::Cancelled),
    ] {
        assert_eq!(failure(cause).kind, kind);
    }
}

#[test]
fn final_cancellation_and_source_precision_refuse_without_publishing_geometry_or_receipts() {
    use std::cell::Cell;
    let source = fixture(1, 1., false);
    let policy = GeometryPolicy::default();
    let calls = Cell::new(0);
    let _ = Body::new(source.clone(), policy, &|| {
        calls.set(calls.get() + 1);
        false
    })
    .unwrap();
    let final_call = calls.get();
    let now = Cell::new(0);
    assert_eq!(
        Body::new(source.clone(), policy, &|| {
            now.set(now.get() + 1);
            now.get() >= final_call
        })
        .err()
        .unwrap()
        .kind,
        GeometryErrorKind::Cancelled
    );
    let body = Body::new(source.clone(), policy, &|| false).unwrap();
    let calls = Cell::new(0);
    assert_eq!(
        body.mesh_with_interrupt(
            Default::default(),
            Solid::rust_mesh_options(Default::default()),
            &|| {
                calls.set(calls.get() + 1);
                calls.get() > 2
            }
        )
        .err()
        .unwrap()
        .kind,
        GeometryErrorKind::Cancelled
    );
    assert!(body.last_mesh_pair.borrow().is_none());
    let head = Head::Cylindrical {
        radius_mm: 1e6,
        height_mm: 4.,
    };
    let authored = Profile::new(head, 1.5, 1.).unwrap();
    let profile =
        Profile::from_retained_segments(&authored.segments().collect::<Vec<_>>()).unwrap();
    let large = PhillipsScrew::new(
        *source.thread(),
        profile,
        SourceSpec {
            head: None,
            ..source.spec()
        },
        source.cavity().clone(),
        100_000,
        &|| false,
    )
    .unwrap();
    assert!(
        Solid::rust_phillips_source_for_test(large, feature(), Default::default(), policy, &|| {
            false
        })
        .is_err()
    );
    let strict = GeometryPolicy {
        minimum_feature_mm: 1e-4,
        ..policy
    };
    assert!(
        Solid::rust_phillips_source_for_test(
            source,
            feature(),
            Default::default(),
            strict,
            &|| false
        )
        .is_err()
    );
}

#[test]
fn actual_large_head_and_admitted_strain_can_refuse_unresolved_world_area_precision() {
    let base = fixture(1, 1., false);
    let head = Head::Crowned {
        radius_mm: 50.,
        height_mm: 25.,
        crown_radius_mm: 75.,
        top_land_radius_mm: 35.,
        neck_radius_mm: 0.1,
    };
    let source = PhillipsScrew::new(
        *base.thread(),
        Profile::new(head, 1.5, 1.).unwrap(),
        SourceSpec {
            head: Some(head),
            ..base.spec()
        },
        base.cavity().clone(),
        100_000,
        &|| false,
    )
    .unwrap();
    let frame = PlaneFrame {
        u: [1. + 4e-11, 0., 0.],
        ..Default::default()
    };
    let solid =
        Solid::rust_phillips_source_for_test(source, feature(), frame, Default::default(), &|| {
            false
        })
        .unwrap();
    assert_eq!(
        solid.faces().err().unwrap().kind,
        GeometryErrorKind::Unsupported
    );
}

#[test]
fn public_source_constructor_and_interrupt_wrapper_keep_the_actual_source() {
    let source = fixture(1, 1., false);
    let spec = source.thread().spec();
    let solid = Solid::rust_phillips_screw_from_source(
        source.clone(),
        feature(),
        Default::default(),
        Default::default(),
    )
    .unwrap();
    assert_eq!(
        solid
            .rust_phillips_screw_geometry()
            .unwrap()
            .0
            .thread()
            .spec(),
        spec
    );
    assert_eq!(
        Solid::rust_phillips_screw_from_source_with_interrupt(
            source.clone(),
            feature(),
            Default::default(),
            Default::default(),
            Some(&|| true)
        )
        .err()
        .unwrap()
        .kind,
        GeometryErrorKind::Cancelled
    );
    assert!(
        Solid::rust_phillips_screw_from_source_with_interrupt(
            source,
            feature(),
            Default::default(),
            Default::default(),
            None
        )
        .is_ok()
    );
}
