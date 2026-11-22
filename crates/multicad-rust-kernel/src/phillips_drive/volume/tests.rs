use super::*;
use multicad_precision::PortableMath as _;
#[test]
fn benchmark_actual_h1_policy_budget() {
    for index in 1..=5 {
        let (spec, domain) = super::super::tests::fixture(index);
        let c = Cutter::new(spec, domain).unwrap();
        let extent = c.envelope.entry.hi.max(c.envelope.tip.hi);
        let width = 1e-7 / 32. * extent * extent;
        let now = std::time::Instant::now();
        let result = c.volume_bounds_with_cancel(width, 100_000, &|| false);
        eprintln!(
            "H1 width={width:e} result={result:?} elapsed={:?}",
            now.elapsed()
        );
        let receipt = result.unwrap();
        let v = I::new(receipt.volume_mm3[0], receipt.volume_mm3[1]);
        let oracle = oracle(spec);
        eprintln!("independent asin/Simpson oracle={oracle:.16e}");
        assert!(v.lo <= oracle && oracle <= v.hi);
    }
}

// Independent pre-existing model regression's asin section formula, evaluated
// without directed jets. This is a cross-check, not the enclosure authority.
fn oracle(s: Spec) -> f64 {
    let Inner::ShallowV {
        angle_deg,
        e_width_mm: e,
    } = s.inner
    else {
        unreachable!()
    };
    let alpha = (angle_deg / 2.).to_radians();
    let beta = s.tilt_angle_deg.to_radians();
    let a = s.outer_flank_angle_deg.to_radians();
    let c = s.floor_angle_deg.to_radians();
    let theta = (s.outer_corner_angle_deg / 2.).to_radians();
    let m = s.mouth_diameter_mm;
    let r = s.mouth_radius_mm;
    let b = s.b_width_mm;
    let tg = (m - s.g_diameter_mm) / (2. * a.portable_tan());
    let total = tg + s.g_diameter_mm / 2. * c.portable_tan();
    let tb = tg + (s.g_diameter_mm - b) / 2. * c.portable_tan();
    let join = r * (1. - a.portable_sin());
    let entry = m / 2. + r * (1. - a.portable_sin()) / a.portable_cos();
    let de = e / (2. * alpha.portable_tan() * beta.portable_cos());
    let first = 1. / (alpha.portable_tan() * beta.portable_cos());
    let second = 1. / (theta.portable_tan() * beta.portable_cos());
    let area = |d: f64| {
        let rr = if d <= join {
            entry - (r * r - (d - r).powi(2)).max(0.).sqrt()
        } else if d <= tg {
            m / 2. - d * a.portable_tan()
        } else {
            (total - d).max(0.) / c.portable_tan()
        };
        if rr == 0. {
            return 0.;
        }
        let rb = b / 2. + (tb - d) * beta.portable_tan();
        let primitive = |y: f64| {
            let y = y.clamp(0., rr);
            (y * (rr * rr - y * y).max(0.).sqrt() + rr * rr * (y / rr).portable_asin()) / 2.
        };
        let part = |aa: f64, bb: f64, lo: f64, upper: f64| {
            let disc = rr * rr * (1. + bb * bb) - aa * aa;
            if disc <= 0. {
                return 0.;
            }
            let hi = upper.min(rr).min((-aa * bb + disc.sqrt()) / (1. + bb * bb));
            if hi <= lo {
                0.
            } else {
                primitive(hi) - primitive(lo) - aa * (hi - lo) - bb * (hi * hi - lo * lo) / 2.
            }
        };
        PI * rr * rr
            - 8. * (part(rb, first, 0., rr.min(e / 2.))
                + part(rb + de - second * e / 2., second, e / 2., rr))
    };
    [0., join, tg, total]
        .windows(2)
        .map(|w| integral(&area, w[0], w[1], 1e-11, 32))
        .sum()
}
fn integral(f: &impl Fn(f64) -> f64, a: f64, b: f64, eps: f64, n: u32) -> f64 {
    let m = (a + b) / 2.;
    let fa = f(a);
    let fb = f(b);
    let fm = f(m);
    let l = (m - a) * (fa + 4. * f((a + m) / 2.) + fm) / 6.;
    let r = (b - m) * (fm + 4. * f((m + b) / 2.) + fb) / 6.;
    let s = (b - a) * (fa + 4. * fm + fb) / 6.;
    if n == 0 || (l + r - s).abs() < 15. * eps {
        l + r + (l + r - s) / 15.
    } else {
        integral(f, a, m, eps / 2., n - 1) + integral(f, m, b, eps / 2., n - 1)
    }
}

#[test]
fn prototype_refuses_unqualified_finite_domains_and_honors_work_cancel() {
    use std::cell::Cell;
    let (s, d) = super::super::tests::fixture(1);
    let c = Cutter::new(s, d).unwrap();
    assert!(matches!(
        compute(&c, 0., 100, &|| false),
        Err(Error::Invalid)
    ));
    assert!(matches!(compute(&c, 1e-8, 1, &|| false), Err(Error::Limit)));
    assert!(matches!(
        compute(&c, 1e-8, 100000, &|| true),
        Err(Error::Cancelled)
    ));
    let count = Cell::new(0);
    assert!(matches!(
        compute(&c, 1e-8, 100000, &|| {
            count.set(count.get() + 1);
            count.get() > 20
        }),
        Err(Error::Cancelled)
    ));
    for domain in [
        Domain {
            corner_reach_mm: 0.01,
            ..d
        },
        Domain {
            corner_start_depth_mm: 0.,
            ..d
        },
        Domain {
            corner_end_depth_mm: 0.1,
            ..d
        },
    ] {
        let finite = Cutter::new(s, domain).unwrap();
        assert!(matches!(
            compute(&finite, 1e-8, 100000, &|| false),
            Err(Error::Unresolved)
        ));
    }
    let overlap = Spec {
        outer_corner_angle_deg: 170.,
        ..s
    };
    assert!(matches!(
        compute(&Cutter::new(overlap, d).unwrap(), 1e-8, 100000, &|| false),
        Err(Error::Unresolved)
    ));
    assert_eq!(c, Cutter::new(s, d).unwrap());
    assert!(c.volume_bounds_mm3(1e-5, 100000).is_ok());
}

#[test]
fn bounded_work_and_reduction_cancel_are_explicit() {
    use std::cell::Cell;
    let (s, d) = super::super::tests::fixture(1);
    let c = Cutter::new(s, d).unwrap();
    assert_eq!(c.volume_bounds_mm3(1e-8, 1_000_001), Err(Error::Limit));
    assert_eq!(c.volume_bounds_mm3(1e-8, 8), Err(Error::Limit));
    // A broad target enters the final outward reduction immediately. Cancel at
    // its first block, rather than allowing a cached success to bypass cancel.
    let calls = Cell::new(0);
    assert_eq!(
        c.volume_bounds_with_cancel(100., 100, &|| {
            calls.set(calls.get() + 1);
            calls.get() == 3
        }),
        Err(Error::Cancelled)
    );
    let result = c.volume_bounds_mm3(100., 100).unwrap();
    assert!(result.volume_mm3[0] >= 0.);
    assert!(result.volume_mm3[1] > result.volume_mm3[0]);
    assert!(result.midpoint_evaluation_width_mm3 >= 0.);
    assert!(result.curvature_and_event_width_mm3 > 0.);
}

#[test]
fn source_section_events_and_numerical_precision_have_real_refusals() {
    let (base, d) = super::super::tests::fixture(1);
    let broad = Spec {
        mouth_radius_mm: 0.001,
        b_width_mm: 2.9,
        g_diameter_mm: 2.95,
        inner: Inner::ShallowV {
            angle_deg: 170.,
            e_width_mm: 2.85,
        },
        ..base
    };
    let source = Cutter::new(broad, d).unwrap();
    assert!(Sections::new(&source).is_ok());
    let e = source.envelope.entry.midpoint() * 2_f64.sqrt();
    let touching = Cutter::new(
        Spec {
            inner: Inner::ShallowV {
                angle_deg: 170.,
                e_width_mm: e,
            },
            ..broad
        },
        d,
    )
    .unwrap();
    assert!(matches!(Sections::new(&touching), Err(Error::Unresolved)));
    // An actual fixed section join contained in the root enclosure cannot select
    // either jet formula. The range path owns this event, even at tiny widths.
    let r = n(2.);
    let a = n(1.);
    let b = n(0.5);
    let crossing = (-a * b + root(sq(r) * (n(1.) + sq(b)) - sq(a)).unwrap()) / (n(1.) + sq(b));
    assert!(part(r, a, b, crossing.v, None).is_none());
    // Independent exact polynomial integral: integral_0^2 x² dx=8/3.
    let (quadratic, _) =
        midpoint_enclosure(I::point(2.), I::point(2.), I::point(2.), I::point(100.)).unwrap();
    assert!(quadratic.lo <= 8. / 3. && quadratic.hi >= 8. / 3.);
    assert!(quadratic.width() < 1e-13);
    assert!(
        midpoint_enclosure(
            I::point(1.),
            I::point(1e308),
            I::point(1e100),
            I::point(1e308)
        )
        .is_none()
    );
    let lo = source.envelope.join.lo;
    assert_eq!(subdivision_point(lo, lo.next_up()), Err(Error::Unresolved));
    assert_eq!(subdivision_point(0., 2.), Ok(1.));
}

#[test]
fn arithmetic_overflow_preserves_the_whole_area_enclosure_and_nonconvex_sections_integrate() {
    // The generic directed quadrature cell must fall back to its geometric
    // range if the curvature product overflows, rather than emitting NaN mass.
    let fallback = evaluated_cell(
        0.,
        1e100,
        I::point(1.),
        I::point(1e308),
        I::point(1e100),
        I::point(1e308),
    );
    assert_eq!(bounds(fallback.v), [0., 1e308]);
    assert_eq!(fallback.numeric.hi, 0.);
    assert_eq!(fallback.quadrature.hi, 1e308);
    let (s, d) = super::super::tests::fixture(1);
    let s = Spec {
        inner: Inner::ShallowV {
            angle_deg: 91.8,
            e_width_mm: 0.435,
        },
        ..s
    };
    assert!(91.8 / 2. < s.outer_corner_angle_deg / 2.);
    let c = Cutter::new(s, d).unwrap();
    let receipt = c.volume_bounds_mm3(2e-8, 100000).unwrap();
    let expected = oracle(s);
    assert!(receipt.volume_mm3[0] <= expected && expected <= receipt.volume_mm3[1]);
}

// Independent H0 section oracle: bracket the actual ellipse/circle crossing
// directly, without the production quadratic reconstruction or interval jets.
fn rounded_oracle_section(s: Spec, depth: f64) -> f64 {
    let Inner::Rounded { radius_mm: radius } = s.inner else {
        unreachable!()
    };
    let beta = s.tilt_angle_deg.to_radians();
    let theta = (s.outer_corner_angle_deg / 2.).to_radians();
    let a = s.outer_flank_angle_deg.to_radians();
    let c = s.floor_angle_deg.to_radians();
    let flank = (s.mouth_diameter_mm - s.g_diameter_mm) / (2. * a.portable_tan());
    let tip = flank + s.g_diameter_mm / 2. * c.portable_tan();
    let join = s.mouth_radius_mm * (1. - a.portable_sin());
    let entry =
        s.mouth_diameter_mm / 2. + s.mouth_radius_mm * (1. - a.portable_sin()) / a.portable_cos();
    let rr = if depth <= join {
        entry
            - (s.mouth_radius_mm.powi(2) - (depth - s.mouth_radius_mm).powi(2))
                .max(0.)
                .sqrt()
    } else if depth <= flank {
        s.mouth_diameter_mm / 2. - depth * a.portable_tan()
    } else {
        (tip - depth).max(0.) / c.portable_tan()
    };
    if rr == 0. {
        return 0.;
    }
    let tb = flank + (s.g_diameter_mm - s.b_width_mm) / 2. * c.portable_tan();
    let rb = s.b_width_mm / 2. + (tb - depth) * beta.portable_tan();
    if rr <= rb {
        return PI * rr * rr;
    }
    let js = radius * (1. - theta.portable_sin());
    let jy = radius * theta.portable_cos();
    let slope = 1. / (theta.portable_tan() * beta.portable_cos());
    let front = |y: f64| {
        rb + if y <= jy {
            (radius - (radius * radius - y * y).max(0.).sqrt()) / beta.portable_cos()
        } else {
            (js + (y - jy) / theta.portable_tan()) / beta.portable_cos()
        }
    };
    let (mut lo, mut hi) = (0., rr);
    for _ in 0..80 {
        let mid = (lo + hi) / 2.;
        if (rr * rr - mid * mid).max(0.).sqrt() > front(mid) {
            lo = mid
        } else {
            hi = mid
        }
    }
    let crossing = (lo + hi) / 2.;
    let upper = crossing.min(jy);
    let primitive = |rad: f64, y: f64| {
        0.5 * (y * (rad * rad - y * y).max(0.).sqrt()
            + rad * rad * (y / rad).clamp(-1., 1.).portable_asin())
    };
    let inner = primitive(rr, upper) - (rb + radius / beta.portable_cos()) * upper
        + primitive(radius, upper) / beta.portable_cos();
    let outer = if crossing > jy {
        primitive(rr, crossing)
            - primitive(rr, jy)
            - (rb + (js - jy / theta.portable_tan()) / beta.portable_cos()) * (crossing - jy)
            - slope * (crossing * crossing - jy * jy) / 2.
    } else {
        0.
    };
    PI * rr * rr - 8. * (inner + outer)
}
fn rounded_oracle_volume(s: Spec) -> f64 {
    let a = s.outer_flank_angle_deg.to_radians();
    let floor = s.floor_angle_deg.to_radians();
    let flank = (s.mouth_diameter_mm - s.g_diameter_mm) / (2. * a.portable_tan());
    let tip = flank + s.g_diameter_mm / 2. * floor.portable_tan();
    [0., s.mouth_radius_mm * (1. - a.portable_sin()), flank, tip]
        .windows(2)
        .map(|w| integral(&|d| rounded_oracle_section(s, d), w[0], w[1], 1e-12, 34))
        .sum()
}
#[test]
fn actual_h0_ellipse_sections_and_volume_meet_policy_budget() {
    let (base, d) = super::super::tests::fixture(0);
    for tilt in [7., 0.000001] {
        let s = Spec {
            tilt_angle_deg: tilt,
            ..base
        };
        let c = Cutter::new(s, d).unwrap();
        let sections = Sections::new(&c).unwrap();
        for fraction in [0.00001, 0.03, 0.3, 0.6, 0.85, 0.95] {
            let z = c.envelope.tip.midpoint() * fraction;
            let branch = if z < c.envelope.join.lo {
                0
            } else if z < c.envelope.flank.lo {
                1
            } else {
                2
            };
            let actual = sections.area(J::d(I::point(z)), branch).unwrap().v;
            let expected = rounded_oracle_section(s, z);
            assert!(
                actual.lo <= expected && expected <= actual.hi,
                "section {z} {actual:?} vs {expected}"
            );
        }
        let extent = c.envelope.entry.hi.max(c.envelope.tip.hi);
        let width = 1e-7 / 32. * extent * extent;
        let now = std::time::Instant::now();
        let receipt = c.volume_bounds_mm3(width, 100000).unwrap();
        let expected = rounded_oracle_volume(s);
        eprintln!(
            "H0 tilt={tilt} width={width:e} {receipt:?} oracle={expected:.16e} elapsed={:?}",
            now.elapsed()
        );
        assert!(receipt.volume_mm3[0] <= expected && expected <= receipt.volume_mm3[1]);
    }
    let broad = Spec {
        mouth_diameter_mm: 3.,
        mouth_radius_mm: 0.001,
        b_width_mm: 2.9,
        g_diameter_mm: 2.95,
        inner: Inner::Rounded { radius_mm: 1.9 },
        ..base
    };
    assert!(Sections::new(&Cutter::new(broad, d).unwrap()).is_ok());
    let capped = Cutter::new(
        base,
        Domain {
            corner_start_depth_mm: 0.,
            ..d
        },
    )
    .unwrap();
    assert_eq!(
        capped.volume_bounds_mm3(1e-8, 100000),
        Err(Error::Unresolved)
    );
}
