use super::*;
use crate::geometry::source_face_measure::{
    Placement, cardinal, local_options, metric, world_point,
};
use multicad_rust_kernel::{
    authored_flat_thread::{AuthoredThread, SourceOptions},
    helical_thread::Spec,
};
use std::cell::Cell;
fn source(left: bool, nb: bool) -> Insert {
    let opts = SourceOptions {
        position_limit_mm: 1e-7,
        normal_angle_limit_rad: 1e-7,
    };
    let outer = Spec {
        major_radius_mm: 5.,
        pitch_mm: 2.,
        length_mm: if nb { 4.3 } else { 4. },
        depth_mm: 1.,
        crest_width_mm: 0.4,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: left,
        root_radius_mm: None,
    };
    let minor = (2.459 + 2.599) / 2.;
    let pd = (2.675 + 2.775) / 2.;
    let inner = Spec {
        major_radius_mm: 1.5,
        pitch_mm: 0.5,
        depth_mm: (3. - minor) / 2.,
        crest_width_mm: 0.25 - (3. - pd) / 3_f64.sqrt(),
        ..outer
    };
    Insert::new(
        AuthoredThread::new(outer, opts).unwrap(),
        AuthoredThread::new(inner, opts).unwrap(),
        0.5,
        None,
    )
    .unwrap()
}
fn options() -> FacePropertyOptions {
    FacePropertyOptions {
        area_tolerance_mm2: 1e-6,
        first_moment_tolerance_mm3: 1e-5,
        centroid_tolerance_mm: 1e-7,
        ..Default::default()
    }
}
fn proper() -> PlaneFrame {
    PlaneFrame {
        origin: [11., -7., 5.],
        u: [0.8, 0.6, 0.],
        v: [-0.36, 0.48, 0.8],
    }
}
fn overlaps(a: Enclosure, b: B) {
    assert!(a.lower <= b.hi && b.lo <= a.upper, "{a:?} vs{b:?}");
}
#[test]
fn actual_source_faces_preserve_asymmetric_centroids_through_proper_affine_pose() {
    let policy = GeometryPolicy::default();
    for left in [false, true] {
        for nb in [false, true] {
            let source = source(left, nb);
            let body = WoodInsertMaterial::new(source, proper(), policy, None).unwrap();
            let before = body;
            let faces = body.faces(policy, None).unwrap();
            assert_eq!(faces.len(), source.face_count() as usize);
            for (index, descriptor) in faces.into_iter().enumerate() {
                let face = source.face_key(index as u32).unwrap();
                let p = body
                    .world_face_properties(face, options(), policy, None)
                    .unwrap();
                assert_eq!(p.local.face, face);
                assert_eq!(p.world.face, face);
                assert_eq!(descriptor.index, index as u32);
                assert!(descriptor.area_mm2 > 0.);
                assert_eq!(descriptor.plane.is_some(), planar(face));
                for k in 0..3 {
                    let expected = B::point(body.frame.origin[k])
                        + (0..3).fold(B::ZERO, |s, j| {
                            s + B::point(Placement::from(body.frame).axes[j][k])
                                * B::from(p.local.centroid_mm[j])
                        });
                    overlaps(p.world.centroid_mm[k], expected);
                    assert!((descriptor.center_mm[k] - p.world.centroid_mm[k].lower).abs() < 1e-6);
                }
                if planar(face) {
                    assert_eq!(p.local_centroid_reweighting_error_mm, [0.; 3]);
                    let (point, n) = descriptor.plane.unwrap();
                    let station = if face == FaceKey::StartCap {
                        0.
                    } else {
                        source.full_length_mm()
                    };
                    let x = [source.outer().maximum_radius_mm(), 0., station];
                    let wp = world_point(Placement::from(body.frame), x.map(B::point));
                    let d = (0..3).fold(B::ZERO, |s, k| {
                        s + B::point(n[k]) * (wp[k] - B::point(point[k]))
                    });
                    assert!(d.lo <= 1e-10 && d.hi >= -1e-10);
                    assert!(
                        p.local.centroid_mm[0].lower.abs() + p.local.centroid_mm[1].lower.abs()
                            > 1e-3
                    );
                }
            }
            assert_eq!(body, before);
        }
    }
}
#[test]
fn actual_small_metric_change_retains_normal_dependent_area_and_explicit_precision_refusal() {
    let src = source(false, false);
    let policy = GeometryPolicy::default();
    let frame = PlaneFrame {
        u: [1. + 2e-12, 0., 0.],
        ..PlaneFrame::default()
    };
    let body = WoodInsertMaterial::new(src, frame, policy, None).unwrap();
    let cap = body
        .world_face_properties(FaceKey::StartCap, options(), policy, None)
        .unwrap();
    assert!(cap.area_jacobian.lower <= 1. + 2e-12 && cap.area_jacobian.upper >= 1. + 2e-12);
    assert_eq!(cap.local_centroid_reweighting_error_mm, [0.; 3]);
    let key = FaceKey::OuterThread {
        start: 0,
        band: multicad_rust_kernel::authored_internal_treatment::ThreadBandRole::Lower,
    };
    let side = body
        .world_face_properties(key, options(), policy, None)
        .unwrap();
    assert!(
        side.local_centroid_reweighting_error_mm
            .into_iter()
            .any(|v| v > 0.)
    );
    let distortion_limited = FacePropertyOptions {
        area_tolerance_mm2: 1e-10,
        ..options()
    };
    assert_eq!(
        body.world_face_properties(key, distortion_limited, policy, None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Unsupported
    );
    let tight = FacePropertyOptions {
        area_tolerance_mm2: 1e-20,
        ..options()
    };
    assert_eq!(
        body.world_face_properties(key, tight, policy, None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert!(body.faces(policy, None).is_ok());
}
#[test]
fn queries_validate_policy_resource_source_and_all_cancel_checkpoints_atomically() {
    let policy = GeometryPolicy::default();
    let body =
        WoodInsertMaterial::new(source(false, false), PlaneFrame::default(), policy, None).unwrap();
    let before = body;
    let key = FaceKey::StartCap;
    assert_eq!(
        body.world_face_properties(
            key,
            FacePropertyOptions {
                maximum_cells: 1,
                ..options()
            },
            policy,
            None
        )
        .unwrap_err()
        .kind,
        GeometryErrorKind::Limit
    );
    assert_eq!(
        body.world_face_properties(
            key,
            FacePropertyOptions {
                area_tolerance_mm2: 0.,
                ..options()
            },
            policy,
            None
        )
        .unwrap_err()
        .kind,
        GeometryErrorKind::Input
    );
    assert_eq!(
        body.world_face_properties(
            FaceKey::OuterThread {
                start: 1,
                band: multicad_rust_kernel::authored_internal_treatment::ThreadBandRole::Root
            },
            options(),
            policy,
            None
        )
        .unwrap_err()
        .kind,
        GeometryErrorKind::Input
    );
    let count = Cell::new(0);
    let poll = || {
        count.set(count.get() + 1);
        false
    };
    body.faces(policy, Some(&poll)).unwrap();
    for at in 1..=count.get() {
        let n = Cell::new(0);
        let cancel = || {
            n.set(n.get() + 1);
            n.get() == at
        };
        assert_eq!(
            body.faces(policy, Some(&cancel)).unwrap_err().kind,
            GeometryErrorKind::Cancelled
        );
    }
    assert_eq!(body, before);
}

#[test]
fn literal_map_primitives_refuse_invalid_singular_reflected_and_exhausted_numeric_bounds() {
    let identity = Placement::from(PlaneFrame::default());
    assert_eq!(metric(identity, false).unwrap().0, B::ONE);
    for axes in [
        [[0.; 3]; 3],
        [[-1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        [[1.01, 0., 0.], [0., 1., 0.], [0., 0., 1.]],
    ] {
        assert!(metric(Placement { axes, ..identity }, false).is_err());
    }
    assert_eq!(
        metric(
            Placement {
                axes: [[f64::NAN; 3]; 3],
                ..identity
            },
            false
        )
        .unwrap_err()
        .kind,
        GeometryErrorKind::Input
    );
    assert!(!cardinal([[1., 1., 0.], [0., 0., 1.], [0., 1., 0.]]));
    assert!(!cardinal([[0.; 3]; 3]));
    assert!(!cardinal([[1., 0., 0.], [1., 0., 0.], [0., 0., 1.]]));
    assert!(positive(B::new(-1., 1.)).is_err());
    assert!(positive(B::new(f64::NAN, 1.)).is_err());
    assert!(
        local_options(
            FacePropertyOptions {
                centroid_tolerance_mm: f64::from_bits(1),
                ..options()
            },
            B::ONE,
            100.,
            2.,
            100.
        )
        .is_err()
    );
    assert!(midpoint([B::new(f64::NAN, 1.); 3]).is_err());
    assert!(center([B::new(-1., 1.); 3], GeometryPolicy::default()).is_err());
    let src = source(false, false);
    assert!(
        plane(
            &src,
            PlaneFrame {
                u: [0.; 3],
                ..PlaneFrame::default()
            },
            FaceKey::StartCap,
            GeometryPolicy::default()
        )
        .is_err()
    );
    assert!(
        plane(
            &src,
            PlaneFrame {
                origin: [1e12; 3],
                ..proper()
            },
            FaceKey::EndCap,
            GeometryPolicy::default()
        )
        .is_err()
    );
    let far = GeometryPolicy {
        maximum_coordinate_mm: f64::MAX,
        ..GeometryPolicy::default()
    };
    assert!(
        transport(
            &src,
            PlaneFrame {
                origin: [f64::MAX; 3],
                ..PlaneFrame::default()
            },
            FaceKey::StartCap,
            options(),
            far,
            None
        )
        .is_err()
    );
    let flip = PlaneFrame {
        u: [-1., 0., 0.],
        ..PlaneFrame::default()
    };
    let actual = WoodInsertMaterial::new(src, flip, GeometryPolicy::default(), None).unwrap();
    assert!(actual.faces(GeometryPolicy::default(), None).is_ok());
}

#[test]
fn world_wall_area_and_moments_agree_with_actual_transformed_tangent_cross_oracle() {
    let source = source(true, true);
    let frame = proper();
    let policy = GeometryPolicy::default();
    let body = WoodInsertMaterial::new(source, frame, policy, None).unwrap();
    let thread = source.outer();
    let spec = thread.spec();
    let rate = -spec.pitch_mm / std::f64::consts::TAU;
    let mut sums = [[0.; 4]; 4];
    for (band, b) in thread.bands().into_iter().enumerate() {
        let role = [1, 2, 3, 0, 1][band];
        let nq = 256;
        let nz = 512;
        let dq = (b.phase_mm[1] - b.phase_mm[0]) / nq as f64;
        let dz = spec.length_mm / nz as f64;
        let m = (b.radius_mm[1] - b.radius_mm[0]) / (b.phase_mm[1] - b.phase_mm[0]);
        let a = Placement::from(frame).axes;
        let map = |p: [f64; 3]| {
            std::array::from_fn::<_, 3, _>(|i| (0..3).map(|j| a[j][i] * p[j]).sum::<f64>())
        };
        for i in 0..nq {
            let q = b.phase_mm[0] + (i as f64 + 0.5) * dq;
            let r = b.radius_mm[0] + m * (q - b.phase_mm[0]);
            for j in 0..nz {
                let z = (j as f64 + 0.5) * dz;
                let (s, c) = ((z - q) / rate).portable_sin_cos();
                let u = map([m * c + r / rate * s, m * s - r / rate * c, 0.]);
                let v = map([-r / rate * s, r / rate * c, 1.]);
                let cross = [
                    u[1] * v[2] - u[2] * v[1],
                    u[2] * v[0] - u[0] * v[2],
                    u[0] * v[1] - u[1] * v[0],
                ];
                let area = cross.into_iter().map(|x| x * x).sum::<f64>().sqrt() * dq * dz;
                let point = map([r * c, r * s, z]);
                for (axis, value) in [
                    area,
                    area * (frame.origin[0] + point[0]),
                    area * (frame.origin[1] + point[1]),
                    area * (frame.origin[2] + point[2]),
                ]
                .into_iter()
                .enumerate()
                {
                    sums[role][axis] += value;
                }
            }
        }
    }
    for (role, band) in [
        multicad_rust_kernel::authored_internal_treatment::ThreadBandRole::Lower,
        multicad_rust_kernel::authored_internal_treatment::ThreadBandRole::Crest,
        multicad_rust_kernel::authored_internal_treatment::ThreadBandRole::Upper,
        multicad_rust_kernel::authored_internal_treatment::ThreadBandRole::Root,
    ]
    .into_iter()
    .enumerate()
    {
        let p = body
            .world_face_properties(
                FaceKey::OuterThread { start: 0, band },
                options(),
                policy,
                None,
            )
            .unwrap();
        assert!((p.world.area_mm2.lower - sums[role][0]).abs() < 5e-4);
        for axis in 0..3 {
            assert!((p.world.first_moment_mm3[axis].lower - sums[role][axis + 1]).abs() < 5e-3);
        }
    }
}
