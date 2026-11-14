use super::source_tests::source;
use super::*;
use multicad_rust_kernel::authored_hanger_bolt::{FaceKey, FacePropertyOptions};
use multicad_rust_kernel::authored_internal_treatment::ThreadBandRole;
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
#[test]
fn every_actual_hanger_face_retains_placed_centroid_area_and_true_planar_station() {
    let policy = GeometryPolicy::default();
    for (left, frame) in [(false, PlaneFrame::default()), (true, proper())] {
        let source = source(left, true);
        let material = HangerMaterial::new(source, frame, policy, None).unwrap();
        let before = material;
        let descriptors = material.faces(policy, None).unwrap();
        assert_eq!(descriptors.len(), 13);
        let axes = [frame.u, frame.v, frame.normal()];
        for (index, descriptor) in descriptors.iter().enumerate() {
            let key = source.face_key(index).unwrap();
            let p = material
                .world_face_properties(key, options(), policy, None)
                .unwrap();
            assert_eq!(p.local.face, key);
            assert_eq!(p.world.face, key);
            assert_eq!(descriptor.index, index as u32);
            assert!(p.world.area_mm2.lower > 0.);
            assert!(
                p.world.area_mm2.lower <= descriptor.area_mm2 + 1e-6
                    && descriptor.area_mm2 <= p.world.area_mm2.upper + 1e-6
            );
            for (i, world_origin) in frame.origin.into_iter().enumerate() {
                let expected = (0..3).fold(I::point(world_origin), |v, j| {
                    v + I::point(axes[j][i]) * I::from(p.local.centroid_mm[j])
                });
                let actual = I::from(p.world.centroid_mm[i]);
                assert!(actual.lo <= expected.hi && expected.lo <= actual.hi);
                assert!((descriptor.center_mm[i] - actual.lo).abs() < 1e-6);
            }
            if let Some((point, normal)) = descriptor.plane {
                assert_eq!(p.local_centroid_reweighting_error_mm, [0.; 3]);
                let z = match key {
                    FaceKey::MetricTip => source.sections()[0].start_mm,
                    FaceKey::MetricShoulder => source.sections()[0].end_mm,
                    FaceKey::WoodShoulder => [0., 0.],
                    FaceKey::WoodTip => [source.wood().full_length_mm(); 2],
                    _ => panic!("A curved wall cannot provide a planar support"),
                };
                let station = z[0].midpoint(z[1]);
                let at = std::array::from_fn::<_, 3, _>(|i| {
                    frame.origin[i] + frame.normal()[i] * station
                });
                let distance = (0..3).map(|i| (at[i] - point[i]) * normal[i]).sum::<f64>();
                assert!(distance.abs() < 1e-10);
                let sign = if key == FaceKey::WoodTip { 1. } else { -1. };
                assert!(
                    (0..3)
                        .map(|i| normal[i] * frame.normal()[i] * sign)
                        .sum::<f64>()
                        > 0.999999999
                );
            }
        }
        assert_eq!(descriptors.iter().filter(|f| f.plane.is_some()).count(), 4);
        let feature = EntityId::new();
        let solid = Solid::rust_hanger(source, feature, frame, policy, None).unwrap();
        assert_eq!(solid.faces().unwrap(), descriptors);
        let mesh = solid
            .manufacturing_mesh_with_receipt(Default::default())
            .unwrap()
            .mesh;
        assert!(
            mesh.face_indices
                .iter()
                .all(|&f| (f as usize) < descriptors.len())
        );
        assert_eq!(
            mesh.face_indices
                .iter()
                .copied()
                .collect::<std::collections::BTreeSet<_>>()
                .len(),
            descriptors.len()
        );
        assert_eq!(material, before);
    }
}
#[test]
fn hanger_world_face_queries_refuse_invalid_and_cancelled_requests_without_changing_source() {
    use std::cell::Cell;
    let policy = GeometryPolicy::default();
    let material = HangerMaterial::new(source(true, true), proper(), policy, None).unwrap();
    let before = material;
    assert!(
        material
            .world_face_properties(
                FaceKey::MetricThread {
                    start: 99,
                    band: ThreadBandRole::Root
                },
                options(),
                policy,
                None
            )
            .is_err()
    );
    for options in [
        FacePropertyOptions {
            maximum_cells: 0,
            ..options()
        },
        FacePropertyOptions {
            maximum_cells: 1,
            ..options()
        },
        FacePropertyOptions {
            area_tolerance_mm2: f64::NAN,
            ..options()
        },
        FacePropertyOptions {
            centroid_tolerance_mm: f64::from_bits(1),
            ..options()
        },
    ] {
        assert!(
            material
                .world_face_properties(FaceKey::MetricTip, options, policy, None)
                .is_err()
        );
    }
    let calls = Cell::new(0usize);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    material.faces(policy, Some(&count)).unwrap();
    let last = calls.get();
    assert!(last > 1);
    for stop in [1, last] {
        calls.set(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() == stop
        };
        assert_eq!(
            material.faces(policy, Some(&cancel)).unwrap_err().kind,
            GeometryErrorKind::Cancelled
        );
    }
    let key = FaceKey::MetricThread {
        start: 0,
        band: ThreadBandRole::Root,
    };
    calls.set(0);
    material
        .world_face_properties(key, options(), policy, Some(&count))
        .unwrap();
    let last = calls.get();
    calls.set(0);
    let cancel = || {
        calls.set(calls.get() + 1);
        calls.get() == last
    };
    assert_eq!(
        material
            .world_face_properties(key, options(), policy, Some(&cancel))
            .unwrap_err()
            .kind,
        GeometryErrorKind::Cancelled
    );
    assert_eq!(calls.get(), last);
    assert_eq!(material, before);
}

const GAUSS: [(f64, f64); 8] = [
    (-0.9602898564975363, 0.1012285362903763),
    (-0.7966664774136267, 0.2223810344533745),
    (-0.525532409916329, 0.3137066458778873),
    (-0.1834346424956498, 0.362683783378362),
    (0.1834346424956498, 0.362683783378362),
    (0.525532409916329, 0.3137066458778873),
    (0.7966664774136267, 0.2223810344533745),
    (0.9602898564975363, 0.1012285362903763),
];
fn independent_root_integral(source: Hanger, frame: PlaneFrame, nq: usize, nz: usize) -> [f64; 4] {
    let thread = source.metric();
    let spec = thread.spec();
    assert_eq!(spec.starts, 1);
    let axes = [frame.u, frame.v, frame.normal()];
    let placed = |v: [f64; 3]| {
        std::array::from_fn::<_, 3, _>(|i| (0..3).map(|j| axes[j][i] * v[j]).sum::<f64>())
    };
    let s = spec.pitch_mm / std::f64::consts::TAU * if spec.left_hand { -1. } else { 1. };
    let endpoints = [0., 1.].map(|t| {
        let p = thread.profile_bounds(2, [t, t]).unwrap().phase_mm;
        p[0].midpoint(p[1])
    });
    let dq = endpoints[1] - endpoints[0];
    let shift = -spec.length_mm - source.spec().plain_length_mm;
    let mut sum = [0.; 4];
    for cell in 0..nq {
        for (u, wu) in GAUSS {
            let t = (cell as f64 + (u + 1.) / 2.) / nq as f64;
            let p = thread.profile_bounds(2, [t, t]).unwrap();
            let r = p.radius_mm[0].midpoint(p.radius_mm[1]);
            let q = p.phase_mm[0].midpoint(p.phase_mm[1]);
            let slope = p.slope[0].midpoint(p.slope[1]);
            for layer in 0..nz {
                for (v, wv) in GAUSS {
                    let z = spec.length_mm * (layer as f64 + (v + 1.) / 2.) / nz as f64;
                    let (sn, cs) = ((z - q) / s).sin_cos();
                    let a = placed([slope * cs + r * sn / s, slope * sn - r * cs / s, 0.]);
                    let b = placed([-r * sn / s, r * cs / s, 1.]);
                    let cross = [
                        a[1] * b[2] - a[2] * b[1],
                        a[2] * b[0] - a[0] * b[2],
                        a[0] * b[1] - a[1] * b[0],
                    ];
                    let jac = cross.into_iter().map(|x| x * x).sum::<f64>().sqrt();
                    let position = placed([r * cs, r * sn, z + shift]);
                    let weight = jac * wu * wv * dq * spec.length_mm / (4. * nq as f64 * nz as f64);
                    sum[0] += weight;
                    for i in 0..3 {
                        sum[1 + i] += weight * (frame.origin[i] + position[i]);
                    }
                }
            }
        }
    }
    sum
}
#[test]
fn actual_placed_metric_root_area_and_moments_match_independent_two_dimensional_tangent_cross_integrals()
 {
    let source = source(true, true);
    let material = HangerMaterial::new(source, proper(), GeometryPolicy::default(), None).unwrap();
    let p = material
        .world_face_properties(
            FaceKey::MetricThread {
                start: 0,
                band: ThreadBandRole::Root,
            },
            options(),
            GeometryPolicy::default(),
            None,
        )
        .unwrap();
    let coarse = independent_root_integral(source, proper(), 16, 64);
    let fine = independent_root_integral(source, proper(), 32, 128);
    let values = std::iter::once(p.world.area_mm2).chain(p.world.first_moment_mm3);
    for ((coarse, fine), bound) in coarse.into_iter().zip(fine).zip(values) {
        assert!(
            (coarse - fine).abs() < 5e-8,
            "numerical oracle did not converge: {coarse} vs{fine}"
        );
        assert!(
            bound.lower - 5e-8 <= fine && fine <= bound.upper + 5e-8,
            "actual world measure {bound:?} vs independent {fine}"
        );
    }
}
