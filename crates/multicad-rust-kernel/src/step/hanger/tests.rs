use super::*;
use crate::authored_metric_thread::Spec;
fn source(pitch: f64, left_hand: bool) -> AuthoredMetricThread {
    AuthoredMetricThread::new(Spec {
        major_radius_mm: 3. * pitch,
        pitch_mm: pitch,
        length_mm: 24. * pitch,
        crest_width_mm: pitch / 8.,
        starts: 1,
        left_hand,
    })
    .unwrap()
}
fn bezier(mut p: Vec<[f64; 3]>, t: f64) -> [f64; 3] {
    for n in (1..p.len()).rev() {
        for i in 0..n {
            p[i] = std::array::from_fn(|k| p[i][k] * (1. - t) + p[i + 1][k] * t);
        }
    }
    p[0]
}
fn actual(surface: &RootSurface, u: f64, v: f64) -> [f64; 3] {
    let n = surface.angular_knots().len() - 1;
    let span = ((v * n as f64) as usize).min(n - 1);
    let t = v * n as f64 - span as f64;
    let count = surface.poles_mm().len() / 3;
    let rows = (0..3)
        .map(|i| {
            bezier(
                surface.poles_mm()[i * count + 11 * span..i * count + 11 * span + 12].to_vec(),
                t,
            )
        })
        .collect::<Vec<_>>();
    let weights = surface.weights();
    let b = [(1. - u).powi(2), 2. * u * (1. - u), u * u];
    let w = (0..3).map(|i| weights[i] * b[i]).sum::<f64>();
    std::array::from_fn(|k| (0..3).map(|i| rows[i][k] * weights[i] * b[i]).sum::<f64>() / w)
}
#[test]
fn actual_rounded_root_tensor_retains_original_source_and_certifies_both_hands_full_length() {
    for pitch in [0.5, 0.3] {
        for left in [false, true] {
            let s = source(pitch, left);
            let surface =
                RootSurface::prepare(s, [-0.625, 23.625], CurveOptions::default(), None).unwrap();
            assert_eq!(*surface.source(), s);
            assert_eq!(surface.physical_turns(), [-0.625, 23.625]);
            let bound = surface.bounds();
            assert!(
                bound.position_error_mm <= 1e-7
                    && bound.derivative_error_mm.into_iter().all(|v| v <= 1e-6)
            );
            assert!(bound.normal_angle_rad < 0.001);
            assert!(
                surface.work_used() <= 16_000_000 && surface.retained_bytes() <= 64 * 1024 * 1024
            );
            assert_eq!(
                surface
                    .certify_controls(surface.poles_mm(), CurveOptions::default(), None)
                    .unwrap(),
                bound
            );
            assert!(surface.angular_knots().len() > 2);
            assert_eq!(
                surface.certify_controls(
                    surface.poles_mm(),
                    CurveOptions {
                        maximum_spans_per_curve: 1,
                        ..CurveOptions::default()
                    },
                    None
                ),
                Err(Error::Limit)
            );
            for j in 0..=8 {
                for k in 0..=16 {
                    let u = f64::from(j) / 8.;
                    let v = f64::from(k) / 16.;
                    let profile = s
                        .root_meridian()
                        .bounds([u, u], None)
                        .unwrap()
                        .point_mm
                        .map(hull);
                    let q = profile[0].midpoint();
                    let r = profile[1].midpoint();
                    let turn = -0.625 + 24.25 * v;
                    let theta = std::f64::consts::TAU * turn * if left { -1. } else { 1. };
                    let expected = [r * theta.cos(), r * theta.sin(), q + pitch * turn];
                    let actual = actual(&surface, u, v);
                    let error = (0..3).map(|i| (actual[i] - expected[i]).abs()).sum::<f64>();
                    assert!(
                        error <= bound.position_error_mm + 1e-12,
                        "{error} {bound:?}"
                    );
                }
            }
        }
    }
}
#[test]
fn actual_coefficient_quality_ownership_budget_and_cancellation_failures_are_atomic() {
    use std::cell::Cell;
    let options = CurveOptions::default();
    let s = source(0.5, false);
    let surface = RootSurface::prepare(s, [0., 0.25], options, None).unwrap();
    let before = surface.poles_mm().to_vec();
    let mut changed = before.clone();
    changed[0][0] += 0.01;
    assert_eq!(
        surface.certify_controls(&changed, options, None),
        Err(Error::Unresolved)
    );
    changed[0][0] = f64::NAN;
    assert_eq!(
        surface.certify_controls(&changed, options, None),
        Err(Error::Invalid)
    );
    assert_eq!(
        surface.certify_controls(&[], options, None),
        Err(Error::Invalid)
    );
    assert_eq!(
        surface.certify_controls(
            surface.poles_mm(),
            CurveOptions { work: 1, ..options },
            None
        ),
        Err(Error::Limit)
    );
    let bytes = surface.retained_bytes() + SCRATCH;
    assert!(
        surface
            .certify_controls(surface.poles_mm(), CurveOptions { bytes, ..options }, None)
            .is_ok()
    );
    assert_eq!(
        surface.certify_controls(
            surface.poles_mm(),
            CurveOptions {
                bytes: bytes - 1,
                ..options
            },
            None
        ),
        Err(Error::Limit)
    );
    assert_eq!(
        surface.certify_controls(&before, CurveOptions { bytes, ..options }, None),
        Err(Error::Limit)
    );
    for domain in [[f64::NAN, 1.], [1., 0.], [0., 257.]] {
        assert!(matches!(
            RootSurface::prepare(s, domain, options, None),
            Err(Error::Invalid)
        ));
    }
    for opts in [
        CurveOptions {
            maximum_spans_per_curve: 0,
            ..options
        },
        CurveOptions {
            position_error_mm: f64::NAN,
            ..options
        },
    ] {
        assert!(matches!(
            RootSurface::prepare(s, [0., 0.25], opts, None),
            Err(Error::Invalid)
        ));
    }
    assert!(matches!(
        RootSurface::prepare(
            s,
            [0., 0.25],
            CurveOptions {
                bytes: 1,
                ..options
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        RootSurface::prepare(
            s,
            [0., 0.25],
            CurveOptions {
                position_error_mm: f64::from_bits(1),
                ..options
            },
            None
        ),
        Err(Error::Unresolved)
    ));
    assert!(matches!(
        RootSurface::prepare(
            s,
            [0., 0.25],
            CurveOptions {
                position_error_mm: 1e-30,
                maximum_spans_per_curve: 1,
                ..options
            },
            None
        ),
        Err(Error::Unresolved)
    ));
    let count = Cell::new(0);
    surface
        .certify_controls(
            surface.poles_mm(),
            options,
            Some(&|| {
                count.set(count.get() + 1);
                false
            }),
        )
        .unwrap();
    for stop in [0, 1, count.get() / 2, count.get() - 1] {
        let n = Cell::new(0);
        assert_eq!(
            surface.certify_controls(
                surface.poles_mm(),
                options,
                Some(&|| {
                    let at = n.get();
                    n.set(at + 1);
                    at == stop
                })
            ),
            Err(Error::Cancelled)
        );
    }
    let count = Cell::new(0);
    RootSurface::prepare(
        s,
        [0., 0.25],
        options,
        Some(&|| {
            count.set(count.get() + 1);
            false
        }),
    )
    .unwrap();
    for stop in [0, 1, count.get() / 2, count.get() - 1] {
        let n = Cell::new(0);
        assert!(matches!(
            RootSurface::prepare(
                s,
                [0., 0.25],
                options,
                Some(&|| {
                    let at = n.get();
                    n.set(at + 1);
                    at == stop
                })
            ),
            Err(Error::Cancelled)
        ));
    }
    assert_eq!(surface.poles_mm(), before);
    assert_eq!(*surface.source(), s);
}
