//! Private experiment: certified sections only, no manufacturing admission.
use super::*;
use crate::tapered_thread::interval::J;
use std::cmp::Ordering;
use std::collections::BinaryHeap;
fn j(v: I) -> J {
    J {
        v,
        ..J::constant(0.)
    }
}
fn n(v: f64) -> J {
    J::constant(v)
}
fn sq(v: J) -> J {
    v * v
}
fn root(v: J) -> Option<J> {
    if v.v.lo <= 0. {
        None
    } else {
        v.sqrt().checked().ok()
    }
}
struct Sections<'a> {
    c: &'a Cutter,
    b0: I,
    b1: I,
    offset: I,
}
impl<'a> Sections<'a> {
    fn new(c: &'a Cutter) -> Result<Self, Error> {
        let k = &c.corners[0];
        let SourceInner::ShallowV { sine, cosine } = k.inner else {
            return Err(Error::Unresolved);
        };
        // Section fronts have slopes b0,b1. With both below one, front(t)-t
        // decreases continuously, including a nonconvex shallow-V junction.
        let b0 = cosine / sine / k.cos_beta;
        let b1 = k.cos_theta / k.sin_theta / k.cos_beta;
        let offset = (k.join_s - k.join_y * k.cos_theta / k.sin_theta) / k.cos_beta;
        let e = &c.envelope;
        let tb = k.sin_beta / k.cos_beta;
        let rb = k.b_half + (k.b_depth - e.tip) * tb;
        // Every removed section point has s>=0. Bound its maximum over the
        // WHOLE cavity box, then prove the authored finite closure and both
        // extrusion caps are inactive before using an infinite-front integral.
        let smax = (k.cos_beta * (I::new(-e.entry.hi, e.entry.hi) - k.b_half)
            + k.sin_beta * (I::new(0., e.tip.hi) - k.b_depth))
            .hi
            .max(0.);
        // If adjacent rotated corners overlapped, monotonicity would put an
        // overlap on the diagonal with t<=entry/sqrt(2). Test that worst t and
        // minimum rb together, without losing the correlation f(t)-t.
        let t = e.entry / I::point(2.).sqrt();
        let f = if t.hi < k.join_y.lo {
            t * cosine / sine
        } else if t.lo > k.join_y.hi {
            k.join_s + (t - k.join_y) * k.cos_theta / k.sin_theta
        } else {
            return Err(Error::Unresolved);
        };
        if b0.lo <= 0.
            || b1.lo <= 0.
            || b0.hi >= 1.
            || b1.hi >= 1.
            || rb.lo <= 0.
            || (rb + offset).lo <= 0.
            || (rb + f / k.cos_beta - t).lo <= 0.
            || smax >= k.far_s.lo
            || (-I::point(c.domain.corner_start_depth_mm) - I::point(smax) * k.sin_beta).lo <= 0.
            || e.tip.hi >= c.domain.corner_end_depth_mm
        {
            return Err(Error::Unresolved);
        }
        Ok(Self { c, b0, b1, offset })
    }
    fn area(&self, d: J, branch: usize) -> Option<J> {
        let e = &self.c.envelope;
        let s = self.c.spec;
        let k = &self.c.corners[0];
        let r = match branch {
            0 => {
                j(e.entry)
                    - root(j(I::point(s.mouth_radius_mm).square()) - sq(d - n(s.mouth_radius_mm)))?
            }
            1 => j(I::point(s.mouth_diameter_mm) / I::point(2.)) - d * j(e.sin_a / e.cos_a),
            _ => (j(e.tip) - d) * j(e.cos_c / e.sin_c),
        };
        if r.v.lo <= 0. {
            return None;
        }
        let rb = j(k.b_half) + (j(k.b_depth) - d) * j(k.sin_beta / k.cos_beta);
        let a = part(r, rb, j(self.b0), I::point(0.), Some(k.join_y))?;
        let b = part(r, rb + j(self.offset), j(self.b1), k.join_y, None)?;
        let value = j(I::new(PI.next_down(), PI.next_up())) * sq(r) - n(8.) * (a + b);
        value.checked().ok()
    }
}
fn primitive(r: J, y: J) -> Option<J> {
    if y.v.lo == 0. && y.v.hi == 0. {
        return Some(n(0.));
    }
    let q = root(sq(r) - sq(y))?;
    Some((y * q + sq(r) * (y / q).atan()) / n(2.))
}
fn part(r: J, a: J, b: J, lo: I, upper: Option<I>) -> Option<J> {
    if r.v.hi <= a.v.lo {
        return Some(n(0.));
    }
    if r.v.lo <= a.v.hi {
        return None;
    }
    let disc = sq(r) * (n(1.) + sq(b)) - sq(a);
    let mut hi = (-a * b + root(disc)?) / (n(1.) + sq(b));
    // Whole-cell order tests only. At any unresolved intersection, clamp,
    // tangency, or inactive/active event the caller uses an area range instead
    // of differentiating across a branch change.
    if let Some(u) = upper {
        if hi.v.lo >= u.hi {
            hi = j(u)
        } else if hi.v.hi > u.lo {
            return None;
        }
    }
    if hi.v.hi <= lo.lo {
        return Some(n(0.));
    }
    if hi.v.lo <= lo.hi {
        return None;
    }
    Some(
        primitive(r, hi)?
            - primitive(r, j(lo))?
            - a * (hi - j(lo))
            - b * (sq(hi) - j(lo.square())) / n(2.),
    )
}
#[derive(Clone, Copy)]
struct Cell {
    lo: f64,
    hi: f64,
    v: I,
    numeric: f64,
}
#[derive(Clone, Copy)]
struct Priority(f64, usize);
impl PartialEq for Priority {
    fn eq(&self, r: &Self) -> bool {
        self.0 == r.0 && self.1 == r.1
    }
}
impl Eq for Priority {}
impl PartialOrd for Priority {
    fn partial_cmp(&self, r: &Self) -> Option<Ordering> {
        Some(self.cmp(r))
    }
}
impl Ord for Priority {
    fn cmp(&self, r: &Self) -> Ordering {
        self.0.total_cmp(&r.0).then(self.1.cmp(&r.1))
    }
}
fn cell(s: &Sections<'_>, lo: f64, hi: f64) -> Cell {
    let e = &s.c.envelope;
    let h = I::point(hi) - I::point(lo);
    let mid = (I::point(lo) + I::point(hi)) / I::point(2.);
    let branch = if hi <= e.join.lo {
        Some(0)
    } else if lo >= e.join.hi && hi <= e.flank.lo {
        Some(1)
    } else if lo >= e.flank.hi && hi <= e.tip.lo {
        Some(2)
    } else {
        None
    };
    let max = I::new(PI.next_down(), PI.next_up()) * e.entry.square() * h;
    if let Some(b) = branch
        && let (Some(whole), Some(m)) = (s.area(J::d(I::new(lo, hi)), b), s.area(J::d(mid), b))
    {
        let numeric = m.v * h;
        // The midpoint Peano kernel is nonnegative and integrates to
        // h^3/24. Its signed second-derivative enclosure is valid only on
        // the C2 section branch proved above; midpoint and h stay outward.
        let v = numeric + whole.dd * h * h * h / I::point(24.);
        if v.valid() {
            return Cell {
                lo,
                hi,
                v: I::new(v.lo.max(0.), v.hi.min(max.hi)),
                numeric: numeric.width(),
            };
        }
    }
    Cell {
        lo,
        hi,
        v: I::new(0., max.hi),
        numeric: 0.,
    }
}
fn volume(
    c: &Cutter,
    width: f64,
    limit: usize,
    cancel: &dyn Fn() -> bool,
) -> Result<(I, usize, f64), Error> {
    if !width.is_finite() || width <= 0. {
        return Err(Error::Invalid);
    }
    if cancel() {
        return Err(Error::Cancelled);
    }
    let s = Sections::new(c)?;
    let e = &c.envelope;
    let mut knots = vec![
        0., e.join.lo, e.join.hi, e.flank.lo, e.flank.hi, e.tip.lo, e.tip.hi,
    ];
    knots.sort_by(f64::total_cmp);
    knots.dedup();
    if limit < knots.len() - 1 {
        return Err(Error::Limit);
    }
    let mut cells: Vec<_> = knots.windows(2).map(|w| cell(&s, w[0], w[1])).collect();
    let mut heap: BinaryHeap<_> = cells
        .iter()
        .enumerate()
        .map(|(i, c)| Priority(c.v.width(), i))
        .collect();
    let mut work = cells.len();
    let mut estimate: f64 = cells.iter().map(|c| c.v.width()).sum();
    loop {
        if cancel() {
            return Err(Error::Cancelled);
        }
        if estimate <= width || work % 512 == 0 {
            let total = cells.iter().fold(I::point(0.), |v, c| v + c.v);
            if total.width() <= width {
                return Ok((total, work, cells.iter().map(|c| c.numeric).sum()));
            }
            estimate = cells.iter().map(|c| c.v.width()).sum();
        }
        if work + 2 > limit {
            return Err(Error::Limit);
        }
        let Priority(_, index) = heap.pop().ok_or(Error::Unresolved)?;
        let old = cells[index];
        let mid = old.lo + (old.hi - old.lo) / 2.;
        if mid <= old.lo || mid >= old.hi {
            return Err(Error::Unresolved);
        }
        let left = cell(&s, old.lo, mid);
        let right = cell(&s, mid, old.hi);
        estimate += left.v.width() + right.v.width() - old.v.width();
        cells[index] = left;
        heap.push(Priority(left.v.width(), index));
        heap.push(Priority(right.v.width(), cells.len()));
        cells.push(right);
        work += 2;
    }
}
#[test]
fn benchmark_actual_h1_policy_budget() {
    for index in 1..=5 {
        let (spec, domain) = super::tests::fixture(index);
        let c = Cutter::new(spec, domain).unwrap();
        let extent = c.envelope.entry.hi.max(c.envelope.tip.hi);
        let width = 1e-7 / 32. * extent * extent;
        let now = std::time::Instant::now();
        let result = volume(&c, width, 100_000, &|| false);
        eprintln!(
            "H1 width={width:e} result={result:?} elapsed={:?}",
            now.elapsed()
        );
        let (v, _, _) = result.unwrap();
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
    let tg = (m - s.g_diameter_mm) / (2. * a.tan());
    let total = tg + s.g_diameter_mm / 2. * c.tan();
    let tb = tg + (s.g_diameter_mm - b) / 2. * c.tan();
    let join = r * (1. - a.sin());
    let entry = m / 2. + r * (1. - a.sin()) / a.cos();
    let de = e / (2. * alpha.tan() * beta.cos());
    let first = 1. / (alpha.tan() * beta.cos());
    let second = 1. / (theta.tan() * beta.cos());
    let area = |d: f64| {
        let rr = if d <= join {
            entry - (r * r - (d - r).powi(2)).max(0.).sqrt()
        } else if d <= tg {
            m / 2. - d * a.tan()
        } else {
            (total - d).max(0.) / c.tan()
        };
        if rr == 0. {
            return 0.;
        }
        let rb = b / 2. + (tb - d) * beta.tan();
        let primitive = |y: f64| {
            let y = y.clamp(0., rr);
            (y * (rr * rr - y * y).max(0.).sqrt() + rr * rr * (y / rr).asin()) / 2.
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
    let (s, d) = super::tests::fixture(1);
    let c = Cutter::new(s, d).unwrap();
    assert!(matches!(
        volume(&c, 0., 100, &|| false),
        Err(Error::Invalid)
    ));
    assert!(matches!(volume(&c, 1e-8, 1, &|| false), Err(Error::Limit)));
    assert!(matches!(
        volume(&c, 1e-8, 100000, &|| true),
        Err(Error::Cancelled)
    ));
    let count = Cell::new(0);
    assert!(matches!(
        volume(&c, 1e-8, 100000, &|| {
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
            volume(&finite, 1e-8, 100000, &|| false),
            Err(Error::Unresolved)
        ));
    }
    let (h0, d0) = super::tests::fixture(0);
    assert!(matches!(
        volume(&Cutter::new(h0, d0).unwrap(), 1e-8, 100000, &|| false),
        Err(Error::Unresolved)
    ));
    let overlap = Spec {
        outer_corner_angle_deg: 170.,
        ..s
    };
    assert!(matches!(
        volume(&Cutter::new(overlap, d).unwrap(), 1e-8, 100000, &|| false),
        Err(Error::Unresolved)
    ));
    assert_eq!(c, Cutter::new(s, d).unwrap());
    assert!(volume(&c, 1e-5, 100000, &|| false).is_ok());
}
