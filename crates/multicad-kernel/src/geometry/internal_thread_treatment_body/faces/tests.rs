use super::*;
use multicad_rust_kernel::Rational;
fn solid(frame: PlaneFrame, left: bool, starts: u32) -> Solid {
    let stock = super::super::tests::fixture_stock(Rational::new(4, 1).unwrap(), frame);
    stock
        .rust_cut_authored_internal_thread_treatments_with_interrupt(
            super::super::tests::cut(frame, left, starts),
            EntityId::from_uuid_bits(9),
            GeometryPolicy::default(),
            None,
        )
        .unwrap()
}
#[test]
fn actual_source_faces_match_roles_areas_centroids_and_mfg_ordinals() {
    for (left, starts, frame) in [
        (false, 1, PlaneFrame::default()),
        (
            true,
            2,
            PlaneFrame {
                origin: [10., 20., 30.],
                u: [0.8, 0.6, 0.],
                v: [-0.36, 0.48, 0.8],
            },
        ),
        (
            false,
            1,
            PlaneFrame {
                u: [1. + 2_f64.powi(-44), 0., 0.],
                ..PlaneFrame::default()
            },
        ),
    ] {
        let solid = if frame.u[0] > 1. {
            solid(PlaneFrame::default(), left, starts)
                .transformed(RigidTransform {
                    rotation: [[frame.u[0], 0., 0.], [0., 1., 0.], [0., 0., 1.]],
                    ..Default::default()
                })
                .unwrap()
        } else {
            solid(frame, left, starts)
        };
        let data = solid.rust_internal_thread_treatment_data().unwrap();
        let faces = solid.faces().unwrap();
        assert_eq!(faces.len(), data.material.face_count() as usize);
        let mesh = solid
            .manufacturing_mesh_with_receipt(Default::default())
            .unwrap();
        assert_eq!(mesh.edge_pairs as usize * 2, mesh.mesh.indices.len());
        for face in faces {
            assert!(mesh.mesh.face_indices.contains(&face.index));
            let key = data.material.face_key(face.index).unwrap();
            assert!(face.area_mm2 > 0. && face.area_mm2.is_finite());
            let receipt = solid
                .rust_internal_thread_treatment_world_face_properties_with_interrupt(
                    key,
                    source::FacePropertyOptions {
                        area_tolerance_mm2: 1e-6,
                        first_moment_tolerance_mm3: 1e-3,
                        centroid_tolerance_mm: 1e-6,
                        ..Default::default()
                    },
                    None,
                )
                .unwrap()
                .unwrap();
            assert!(
                receipt.world.area_mm2.lower <= face.area_mm2
                    && face.area_mm2 <= receipt.world.area_mm2.upper
            );
            for (v, b) in face.center_mm.into_iter().zip(receipt.world.centroid_mm) {
                assert!(b.lower <= v && v <= b.upper);
            }
            if let Some((p, n)) = face.plane {
                assert!(matches!(
                    key,
                    FaceKey::StartCap | FaceKey::EndCap | FaceKey::TerminalShoulder
                ));
                assert!((super::super::super::dot(n, n) - 1.).abs() < 1e-12);
                assert!(super::super::super::dot(n, frame.u).abs() < 1e-12);
                assert!(super::super::super::dot(n, frame.v).abs() < 1e-12);
                let delta = std::array::from_fn(|i| face.center_mm[i] - p[i]);
                assert!(
                    super::super::super::dot(n, delta).abs()
                        < GeometryPolicy::default().linear_tolerance_mm
                );
            } else {
                assert!(!matches!(
                    key,
                    FaceKey::StartCap | FaceKey::EndCap | FaceKey::TerminalShoulder
                ));
            }
        }
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("unqualified.step");
        assert!(solid.write_step(&path).is_err());
        assert!(!path.exists());
    }
}
#[test]
fn source_domain_and_precision_failures_publish_no_partial_faces() {
    let solid = solid(PlaneFrame::default(), false, 1);
    let (body, frame) = solid.internal_thread_treatment_body().unwrap();
    let strict = GeometryPolicy {
        linear_tolerance_mm: 1e-15,
        ..Default::default()
    };
    assert!(body.faces(frame, strict).is_err());
    assert!(
        body.faces(
            PlaneFrame {
                origin: [1e7, 0., 0.],
                ..frame
            },
            GeometryPolicy::default()
        )
        .is_err()
    );
    assert!(
        body.faces(
            PlaneFrame {
                u: [0.; 3],
                ..frame
            },
            GeometryPolicy::default()
        )
        .is_err()
    );

    let far = PlaneFrame {
        origin: [1e6; 3],
        u: [0.8, 0.6, 0.],
        v: [-0.36, 0.48, 0.8],
    };
    assert!(
        body.faces(
            far,
            GeometryPolicy {
                linear_tolerance_mm: 2e-9,
                ..GeometryPolicy::default()
            }
        )
        .is_err()
    );
    assert!(midpoint([B::point(f64::NAN); 3]).is_err());
    // Each coordinate could meet a separate width budget while the 3D
    // descriptor's total positional error still exceeds its publication limit.
    assert!(centroid([B::new(-1e-8, 1e-8); 3], GeometryPolicy::default()).is_err());
    assert!(
        plane(
            &body.material,
            PlaneFrame {
                u: [0.; 3],
                ..frame
            },
            FaceKey::StartCap,
            GeometryPolicy::default()
        )
        .is_err()
    );
    assert!(
        plane(
            &body.material,
            PlaneFrame {
                origin: [1e6; 3],
                ..frame
            },
            FaceKey::EndCap,
            GeometryPolicy {
                linear_tolerance_mm: 1e-15,
                ..GeometryPolicy::default()
            }
        )
        .is_err()
    );
    let stock = super::super::tests::fixture_stock(Rational::new(4, 1).unwrap(), frame);
    let mut request = super::super::tests::cut(frame, false, 1);
    request.intent.entry_depth_mm = Some(0.01);
    let small_entry = stock
        .rust_cut_authored_internal_thread_treatments_with_interrupt(
            request,
            EntityId::from_uuid_bits(10),
            GeometryPolicy::default(),
            None,
        )
        .unwrap();
    assert!(small_entry.faces().is_err());
    let before = solid.faces().unwrap();
    assert_eq!(solid.faces().unwrap(), before);
}
#[test]
fn profile_end_cap_keeps_actual_asymmetric_centroid_without_circle_inference() {
    let frame = PlaneFrame::default();
    let stock = super::super::tests::fixture_stock(Rational::new(4, 1).unwrap(), frame);
    let mut request = super::super::tests::cut(frame, false, 1);
    request.intent = source::Intent {
        entry_depth_mm: Some(0.5),
        terminal: None,
    };
    let actual = stock
        .rust_cut_authored_internal_thread_treatments_with_interrupt(
            request,
            EntityId::from_uuid_bits(12),
            GeometryPolicy::default(),
            None,
        )
        .unwrap();
    let faces = actual.faces().unwrap();
    let data = actual.rust_internal_thread_treatment_data().unwrap();
    assert_eq!(faces.len(), data.material.face_count() as usize);
    assert_eq!(faces.iter().filter(|f| f.plane.is_some()).count(), 2);
    let face = &faces[1];
    assert!(face.center_mm[0].hypot(face.center_mm[1]) > 1e-4);
    assert_eq!(face.plane.unwrap().0[2], 4.);
    assert!((face.center_mm[2] - 4.).abs() <= GeometryPolicy::default().linear_tolerance_mm);
}
