use super::*;
fn i(x: [f64; 2]) -> I {
    I::new(x[0], x[1])
}
fn overlaps(a: I, b: I) {
    assert!(
        a.valid() && b.valid() && a.lo <= b.hi && b.lo <= a.hi,
        "{a:?} {b:?}"
    );
}
fn source(r: f64, p: f64, left: bool, starts: u32) -> AuthoredMetricThread {
    AuthoredMetricThread::new(Spec {
        major_radius_mm: r,
        pitch_mm: p,
        length_mm: 12. * p,
        crest_width_mm: p / 8.,
        starts,
        left_hand: left,
    })
    .unwrap()
}
#[test]
fn symbolic_controls_cover_the_same_root_circle_and_original_c2_profile() {
    for (r, p) in [(1.4635, 0.5), (2.933, 1.), (1.7, 0.3), (1e7, 0.001)] {
        for left in [false, true] {
            for starts in [1, 2] {
                let s = source(r, p, left, starts);
                let curve = s.root_meridian();
                assert_eq!(*curve.source(), s);
                let c = curve.homogeneous_controls();
                assert_eq!(c.map(|v| v[2]), [[1., 1.], [0.5, 0.5], [1., 1.]]);
                let a = s.knot(Knot::UpperRoot);
                let b = s.knot(Knot::LowerRoot);
                for (t, k) in [(0., a), (1., b)] {
                    let q = curve.bounds([t, t], None).unwrap();
                    overlaps(i(q.point_mm[0]), i(k.phase_mm));
                    overlaps(i(q.point_mm[1]), i(k.radius_mm));
                    let slope = i(q.first_derivative_mm[1]) / i(q.first_derivative_mm[0]);
                    overlaps(slope, sqrt3() * I::point(if t == 0. { -1. } else { 1. }));
                }
                let whole = curve.bounds([0., 1.], None).unwrap();
                assert!(whole.first_derivative_mm[0][0] > 0.);
                for j in 1..32 {
                    let t = f64::from(j) / 32.;
                    let q = curve.bounds([t, t], None).unwrap();
                    let phase = i(q.point_mm[0]);
                    let radius = i(q.point_mm[1]);
                    let circle = (phase - I::point(p / 2.)).square()
                        + (radius - s.center_radius()).square()
                        - s.root_radius().square();
                    assert!(circle.lo <= 0. && circle.hi >= 0.);
                    let u = (phase - I::point(p) * I::point(3.) / I::point(8.))
                        / (I::point(p) / I::point(4.));
                    let original = s.profile_bounds(2, [u.lo, u.hi]).unwrap();
                    overlaps(radius, i(original.radius_mm));
                    let dq = i(q.first_derivative_mm[0]);
                    let dr = i(q.first_derivative_mm[1]);
                    overlaps(dr, i(original.slope) * dq);
                    overlaps(
                        i(q.second_derivative_mm[1]),
                        i(original.second_derivative_per_mm) * dq.square()
                            + i(original.slope) * i(q.second_derivative_mm[0]),
                    );
                    // Independent homogeneous Bernstein evaluation of THESE
                    // source-owned controls, not a fit through sampled endpoints.
                    let t = I::point(t);
                    let v = I::point(1.) - t;
                    let h: [I; 3] = std::array::from_fn(|k| {
                        i(c[0][k]) * v.square()
                            + I::point(2.) * i(c[1][k]) * t * v
                            + i(c[2][k]) * t.square()
                    });
                    overlaps(phase, h[0] / h[2]);
                    overlaps(radius, h[1] / h[2]);
                }
            }
        }
    }
}
#[test]
fn original_owner_cancellation_and_invalid_trim_are_atomic() {
    use std::cell::Cell;
    let s = source(1.4635, 0.5, false, 1);
    let curve = s.root_meridian();
    let before = curve.homogeneous_controls();
    for domain in [
        [f64::NAN, 0.],
        [0.5, 0.25],
        [-f64::MIN_POSITIVE, 0.],
        [0., 1_f64.next_up()],
    ] {
        assert_eq!(curve.bounds(domain, None), Err(Error::Invalid));
    }
    for stop in 0..2 {
        let at = Cell::new(0);
        assert_eq!(
            curve.bounds(
                [0.25, 0.75],
                Some(&|| {
                    let n = at.get();
                    at.set(n + 1);
                    n == stop
                })
            ),
            Err(Error::Cancelled)
        );
    }
    let count = Cell::new(0);
    assert!(
        curve
            .bounds(
                [0.25, 0.75],
                Some(&|| {
                    count.set(count.get() + 1);
                    false
                })
            )
            .is_ok()
    );
    assert_eq!(count.get(), 2);
    assert_eq!(curve.homogeneous_controls(), before);
    assert_eq!(*curve.source(), s);
}
