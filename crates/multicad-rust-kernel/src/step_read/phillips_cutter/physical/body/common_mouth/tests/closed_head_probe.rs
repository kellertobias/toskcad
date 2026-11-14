//! Source-only feasibility for a separately repaired closed head, not raw mass.
use super::*;
use crate::{head_profile::Curve as HeadCurve, tapered_thread::interval::I};
fn ray(p: [f64; 2], c: [f64; 2]) -> [I; 2] {
    if p[0] == c[0] {
        return [I::point(0.), I::point(if p[1] > c[1] { 1. } else { -1. })];
    }
    if p[1] == c[1] {
        return [I::point(if p[0] > c[0] { 1. } else { -1. }), I::point(0.)];
    }
    let d = std::array::from_fn::<_, 2, _>(|i| I::point(p[i]) - I::point(c[i]));
    let n = (d[0].square() + d[1].square()).sqrt();
    [d[0] / n, d[1] / n]
}
#[test]
fn actual_h0_normalized_minor_arcs_admit_cubic_endpoint_correction_and_z_order() {
    for scale in [1., 1.137] {
        for posed in [false, true] {
            let text =
                super::super::super::super::branches::representative::scaled_text(scale, posed);
            let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
            let profile = body.retained_meridian_header();
            let repaired = profile
                .closed_repaired_meridian(Default::default(), None)
                .unwrap();
            assert!(std::ptr::eq(repaired.profile_proposal(), profile));
            assert_eq!(repaired.segment_count(), profile.segments().count());
            for i in 0..repaired.segment_count() {
                let c = repaired.segment_certificate(i).unwrap();
                assert!(c.minimum_radius_mm > 0. && c.minimum_progress_mm > 0.);
                let full = repaired.evaluate(i, [0., 1.]).unwrap();
                assert!(
                    full.position_mm
                        .into_iter()
                        .chain(full.derivative_mm)
                        .flatten()
                        .all(f64::is_finite)
                );
            }
            let mut arcs = 0;
            for s in profile.segments() {
                let HeadCurve::Arc {
                    center,
                    radius_mm,
                    sweep,
                } = s.curve
                else {
                    continue;
                };
                let a = ray(s.start, center);
                let b = ray(s.end, center);
                let det = a[0] * b[1] - a[1] * b[0];
                let radius = I::point(radius_mm);
                let da = std::array::from_fn::<_, 2, _>(|i| {
                    I::point(s.start[i]) - I::point(center[i]) - radius * a[i]
                });
                let db = std::array::from_fn::<_, 2, _>(|i| {
                    I::point(s.end[i]) - I::point(center[i]) - radius * b[i]
                });
                let delta = (db[1] - da[1]).magnitude();
                let sign = if sweep > 0. { 1. } else { -1. };
                assert!((det * I::point(sign)).lo > 0.);
                let coefficient = if s.start[0] == center[0] {
                    (radius * det * b[0] - I::point(6.) * I::point(delta)).lo
                } else if s.end[0] == center[0] {
                    (radius * det * a[0] - I::point(6.) * I::point(delta)).lo
                } else {
                    let det_lower = if det.lo > 0. { det.lo } else { -det.hi };
                    let x_lower = if sign > 0. {
                        a[0].lo.min(b[0].lo)
                    } else {
                        (-a[0].hi).min(-b[0].hi)
                    };
                    (radius * I::point(det_lower) * I::point(x_lower)
                        - I::point(1.5) * I::point(delta))
                    .lo
                };
                let transfer = da
                    .into_iter()
                    .chain(db)
                    .map(I::magnitude)
                    .fold(0., f64::max);
                eprintln!(
                    "head scale{scale} pose{posed} {:?} rays{a:?}/{b:?} det{det:?} dzCoefficient{coefficient} coordinateTransfer{transfer}",
                    s.role
                );
                assert!(coefficient > 0.);
                assert!(transfer < profile.rounding_error_mm());
                arcs += 1;
            }
            assert_eq!(arcs, 2);
        }
    }
}
