use super::*;
use crate::authored_internal_treatment::{FaceKey as ThreadFace, ThreadBandRole};
use multicad_precision::PortableMath as _;
use std::cell::Cell;
fn width(e: measure::Enclosure) -> f64 {
    e.upper - e.lower
}
fn mid(e: measure::Enclosure) -> f64 {
    (e.upper + e.lower) * 0.5
}
fn contains(e: measure::Enclosure, x: f64, error: f64) {
    assert!(e.lower - error <= x && x <= e.upper + error, "{e:?} vs{x}");
}
#[test]
fn actual_paired_face_measures_retain_phase_asymmetric_caps_and_semantic_roles() {
    for left in [false, true] {
        for starts in [1, 2] {
            for nb in [false, true] {
                let body = crate::wood_insert::tests::fixture(left, starts, nb);
                let before = body;
                let opts = FacePropertyOptions::default();
                let mut caps = Vec::new();
                for index in 0..body.face_count() {
                    let key = body.face_key(index).unwrap();
                    assert_eq!(body.face_index(key), Some(index));
                    let value = body.local_face_properties(key, opts, None).unwrap();
                    assert_eq!(value.face, key);
                    assert!(
                        value.area_mm2.lower > 0.
                            && width(value.area_mm2) <= opts.area_tolerance_mm2
                    );
                    assert!(
                        value
                            .first_moment_mm3
                            .into_iter()
                            .all(|v| width(v) <= opts.first_moment_tolerance_mm3)
                    );
                    assert!(
                        value
                            .centroid_mm
                            .into_iter()
                            .all(|v| width(v) <= opts.centroid_tolerance_mm)
                    );
                    match key {
                        FaceKey::StartCap | FaceKey::EndCap => caps.push(value),
                        FaceKey::OuterThread { start, band }
                        | FaceKey::InnerThread { start, band } => {
                            let t = if matches!(key, FaceKey::OuterThread { .. }) {
                                body.outer()
                            } else {
                                body.inner()
                            };
                            let bare = t
                                .local_face_properties(
                                    ThreadFace::Thread { start, band },
                                    opts,
                                    None,
                                )
                                .unwrap();
                            assert_eq!(value.area_mm2, bare.area_mm2);
                            assert_eq!(value.first_moment_mm3, bare.first_moment_mm3);
                        }
                    }
                }
                assert_eq!(caps.len(), 2);
                let volume = body
                    .volume(
                        QueryOptions {
                            volume_error_mm3: 1e-7,
                            maximum_cells: 10,
                        },
                        None,
                    )
                    .unwrap();
                contains(
                    caps[0].area_mm2,
                    (volume.material_mm3[0] + volume.material_mm3[1]) * 0.5 / body.full_length_mm(),
                    1e-9,
                );
                contains(caps[0].centroid_mm[2], 0., 1e-12);
                contains(caps[1].centroid_mm[2], body.full_length_mm(), 1e-12);
                if starts == 1 {
                    assert!(
                        mid(caps[0].centroid_mm[0]).abs() + mid(caps[0].centroid_mm[1]).abs()
                            > 1e-3
                    );
                }
                assert_eq!(body.face_key(body.face_count()), None);
                assert_eq!(
                    body.face_index(FaceKey::InnerThread {
                        start: starts,
                        band: ThreadBandRole::Root
                    }),
                    None
                );
                assert_eq!(body, before);
            }
        }
    }
}
#[test]
fn cap_subtraction_matches_independent_source_angular_integral() {
    let body = crate::wood_insert::tests::fixture(false, 1, true);
    for (key, z) in [
        (FaceKey::StartCap, 0.),
        (FaceKey::EndCap, body.full_length_mm()),
    ] {
        let result = body
            .local_face_properties(key, FacePropertyOptions::default(), None)
            .unwrap();
        let count = 200000;
        let mut oracle = [0.; 4];
        for j in 0..count {
            let theta = std::f64::consts::TAU * (j as f64 + 0.5) / count as f64;
            let radius = |thread: &AuthoredThread| {
                let spec = thread.spec();
                let lead = spec.pitch_mm * f64::from(spec.starts);
                let hand = if spec.left_hand { -1. } else { 1. };
                let q = (z - hand * lead * theta / std::f64::consts::TAU).rem_euclid(spec.pitch_mm);
                let b = thread
                    .bands()
                    .into_iter()
                    .find(|b| q <= b.phase_mm[1])
                    .unwrap();
                b.radius_mm[0]
                    + (b.radius_mm[1] - b.radius_mm[0]) * (q - b.phase_mm[0])
                        / (b.phase_mm[1] - b.phase_mm[0])
            };
            let a = radius(body.outer());
            let b = radius(body.inner());
            let da = (a * a - b * b) * 0.5;
            let dm = (a * a * a - b * b * b) / 3.;
            for (i, v) in [
                da,
                dm * theta.portable_cos(),
                dm * theta.portable_sin(),
                z * da,
            ]
            .into_iter()
            .enumerate()
            {
                oracle[i] += v * std::f64::consts::TAU / count as f64;
            }
        }
        contains(result.area_mm2, oracle[0], 1e-7);
        for i in 0..3 {
            contains(result.first_moment_mm3[i], oracle[i + 1], 2e-6);
        }
    }
}
#[test]
fn invalid_precision_resource_and_cancellation_never_publish_partial_properties() {
    let body = crate::wood_insert::tests::fixture(false, 1, false);
    let before = body;
    let key = FaceKey::StartCap;
    let opts = FacePropertyOptions::default();
    assert_eq!(
        body.local_face_properties(
            FaceKey::OuterThread {
                start: 1,
                band: ThreadBandRole::Root
            },
            opts,
            None
        ),
        Err(Error::Invalid)
    );
    for bad in [0., f64::NAN, f64::INFINITY] {
        assert_eq!(
            body.local_face_properties(
                key,
                FacePropertyOptions {
                    area_tolerance_mm2: bad,
                    ..opts
                },
                None
            ),
            Err(Error::Invalid)
        );
    }
    assert_eq!(
        body.local_face_properties(
            key,
            FacePropertyOptions {
                maximum_cells: 1,
                ..opts
            },
            None
        ),
        Err(Error::Limit)
    );
    assert_eq!(
        body.local_face_properties(
            key,
            FacePropertyOptions {
                area_tolerance_mm2: 1e-30,
                ..opts
            },
            None
        ),
        Err(Error::Unresolved)
    );
    let count = Cell::new(0);
    let count_calls = || {
        count.set(count.get() + 1);
        false
    };
    body.local_face_properties(key, opts, Some(&count_calls))
        .unwrap();
    for at in 1..=count.get() {
        let calls = Cell::new(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() == at
        };
        assert_eq!(
            body.local_face_properties(key, opts, Some(&cancel)),
            Err(Error::Cancelled)
        );
    }
    assert_eq!(body, before);
}

#[test]
fn bare_thread_wall_measures_match_independent_tangent_cross_quadrature() {
    for left in [false, true] {
        let body = crate::wood_insert::tests::fixture(left, 1, true);
        for thread in [body.outer(), body.inner()] {
            let spec = thread.spec();
            let signed_rate =
                (if spec.left_hand { -1. } else { 1. }) * spec.pitch_mm / std::f64::consts::TAU;
            let mut reference = [[0.; 4]; 4];
            for (band, b) in thread.bands().into_iter().enumerate() {
                let role = [1, 2, 3, 0, 1][band];
                let m = (b.radius_mm[1] - b.radius_mm[0]) / (b.phase_mm[1] - b.phase_mm[0]);
                let nq = 512;
                let nz = 512;
                let dq = (b.phase_mm[1] - b.phase_mm[0]) / nq as f64;
                let dz = spec.length_mm / nz as f64;
                for i in 0..nq {
                    let q = b.phase_mm[0] + (i as f64 + 0.5) * dq;
                    let r = b.radius_mm[0] + m * (q - b.phase_mm[0]);
                    for j in 0..nz {
                        let z = (j as f64 + 0.5) * dz;
                        let (s, c) = ((z - q) / signed_rate).portable_sin_cos();
                        let u = [m * c + r / signed_rate * s, m * s - r / signed_rate * c, 0.];
                        let v = [-r / signed_rate * s, r / signed_rate * c, 1.];
                        let cross = [u[1] * v[2], -u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
                        let area = cross.into_iter().map(|a| a * a).sum::<f64>().sqrt() * dq * dz;
                        for (axis, value) in [area, area * r * c, area * r * s, area * z]
                            .into_iter()
                            .enumerate()
                        {
                            reference[role][axis] += value;
                        }
                    }
                }
            }
            for (role, band) in [
                ThreadBandRole::Lower,
                ThreadBandRole::Crest,
                ThreadBandRole::Upper,
                ThreadBandRole::Root,
            ]
            .into_iter()
            .enumerate()
            {
                let value = thread
                    .local_face_properties(
                        ThreadFace::Thread { start: 0, band },
                        FacePropertyOptions::default(),
                        None,
                    )
                    .unwrap();
                contains(value.area_mm2, reference[role][0], 1e-4);
                for axis in 0..3 {
                    contains(
                        value.first_moment_mm3[axis],
                        reference[role][axis + 1],
                        1e-3,
                    );
                }
            }
        }
    }
}

#[test]
fn bare_property_domain_protocol_and_actual_canonical_default_are_explicit() {
    let body = crate::wood_insert::tests::fixture(false, 1, false);
    let options = FacePropertyOptions::default();
    assert_eq!(
        body.outer()
            .local_face_properties(ThreadFace::EntryCone, options, None),
        Err(measure::Error::Invalid)
    );
    assert_eq!(
        bare_properties(
            body.outer(),
            Some(body.inner()),
            ThreadFace::Thread {
                start: 0,
                band: ThreadBandRole::Root
            },
            options,
            None
        ),
        Err(measure::Error::Invalid)
    );
    let unequal = crate::wood_insert::tests::fixture(false, 1, true);
    assert_eq!(
        bare_properties(
            body.outer(),
            Some(unequal.inner()),
            ThreadFace::StartCap,
            options,
            None
        ),
        Err(measure::Error::Invalid)
    );
    let mut spec = body.outer().spec();
    spec.crest_width_mm = 0.4;
    spec.lower_flank_deg = 30.;
    spec.upper_flank_deg = 30.;
    let thread = AuthoredThread::new(
        spec,
        crate::authored_flat_thread::SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-7,
        },
    )
    .unwrap();
    let canonical = Insert::new(thread, *body.inner(), 0.5, None).unwrap();
    for index in 0..canonical.face_count() {
        let key = canonical.face_key(index).unwrap();
        assert!(
            canonical
                .local_face_properties(key, options, None)
                .unwrap()
                .area_mm2
                .lower
                > 0.
        );
    }
}
