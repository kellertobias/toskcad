//! Proven chord bounds for curves given in a two-coordinate chart, either
//! as an explicit graph `P(x)` or as a component of the zero set of a
//! relation `g(x, y)` with a point map `P(x, y)`.
//!
//! Every bound is computed with outward-rounded interval arithmetic and
//! second-order jets in both chart coordinates. For one chord between the
//! published points of two parameters:
//!
//! - each published point is compared with a rigorous enclosure of an
//!   exact curve point next to it (`ends`, see [`settle`]), giving `eps`;
//! - the exact piece is proven to stay inside a box: its native coordinate
//!   range is known by construction and the relation keeps one sign on the
//!   two box edges across the other coordinate, which the connected piece
//!   cannot cross (the edges are padded until that is proven);
//! - where `g_y` (or `g_x`) excludes zero on the box, the zero set inside it
//!   is a graph over `x` (or `y`), so the piece is that graph between its
//!   end coordinates and its second derivative follows from implicit
//!   differentiation of the jets. With `h` the coordinate range, the chord
//!   deviates by at most `h^2 / 8 sup |P''|`;
//! - otherwise (or when smaller) the diameter of the box's image bounds it.
//!
//! The deviation from the published chord is then at most that bound plus
//! `eps`, since moving both chord ends by at most `eps` moves every chord
//! point by at most `eps`.
use super::super::P3;
use super::Path;
use std::ops::{Add, Div, Mul, Neg, Sub};

/// A closed interval with outward-rounded arithmetic.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(in crate::drawing) struct Iv {
    pub lo: f64,
    pub hi: f64,
}

const ENTIRE: Iv = Iv {
    lo: f64::NEG_INFINITY,
    hi: f64::INFINITY,
};

/// Bound of the absolute error of the platform's `sin`, `cos`, `atan` and
/// `acos` for results of magnitude at most 4, with a wide margin over their
/// documented sub-ulp accuracy.
const LIBM: f64 = 4e-15;

impl Iv {
    pub fn new(lo: f64, hi: f64) -> Self {
        if lo.is_nan() || hi.is_nan() || lo > hi {
            return ENTIRE;
        }
        Self { lo, hi }
    }

    pub fn point(x: f64) -> Self {
        Self::new(x, x)
    }

    fn out(lo: f64, hi: f64) -> Self {
        Self::new(lo.next_down(), hi.next_up())
    }

    pub fn hull(self, other: Self) -> Self {
        Self::new(self.lo.min(other.lo), self.hi.max(other.hi))
    }

    pub fn mag(self) -> f64 {
        self.lo.abs().max(self.hi.abs())
    }

    pub fn width(self) -> f64 {
        (self.hi - self.lo).next_up()
    }

    pub fn excludes_zero(self) -> bool {
        self.lo > 0. || self.hi < 0.
    }

    fn widened(self, by: f64) -> Self {
        Self::out(self.lo - by, self.hi + by)
    }

    fn sqr(self) -> Self {
        if self.lo >= 0. {
            Self::out(self.lo * self.lo, self.hi * self.hi)
        } else if self.hi <= 0. {
            Self::out(self.hi * self.hi, self.lo * self.lo)
        } else {
            Self::new(0., (self.mag() * self.mag()).next_up())
        }
    }

    fn sqrt(self) -> Self {
        if self.hi < 0. {
            return ENTIRE;
        }
        // IEEE square roots are correctly rounded.
        Self::new(
            self.lo.max(0.).sqrt().next_down().max(0.),
            self.hi.sqrt().next_up(),
        )
    }

    fn recip(self) -> Self {
        if !self.excludes_zero() {
            return ENTIRE;
        }
        Self::out(1. / self.hi, 1. / self.lo)
    }

    /// Whether `[lo, hi]` may contain `center + j period` for an integer
    /// `j`; generous at the ends, which only widens the result.
    fn reaches(self, center: f64, period: f64) -> bool {
        let j = ((self.lo - center) / period).floor();
        [j, j + 1., j + 2.].into_iter().any(|j| {
            let at = center + j * period;
            at >= self.lo - 1e-9 && at <= self.hi + 1e-9
        })
    }

    fn sin(self) -> Self {
        use std::f64::consts::{FRAC_PI_2, TAU};
        if !(self.lo.is_finite() && self.hi.is_finite()) || self.hi - self.lo >= 6. {
            return Self::new(-1., 1.);
        }
        let [a, b] = [self.lo.sin(), self.hi.sin()];
        let hi = if self.reaches(FRAC_PI_2, TAU) {
            1.
        } else {
            a.max(b)
        };
        let lo = if self.reaches(-FRAC_PI_2, TAU) {
            -1.
        } else {
            a.min(b)
        };
        Self::new((lo - LIBM).max(-1.), (hi + LIBM).min(1.))
    }

    fn cos(self) -> Self {
        use std::f64::consts::{PI, TAU};
        if !(self.lo.is_finite() && self.hi.is_finite()) || self.hi - self.lo >= 6. {
            return Self::new(-1., 1.);
        }
        let [a, b] = [self.lo.cos(), self.hi.cos()];
        let hi = if self.reaches(0., TAU) { 1. } else { a.max(b) };
        let lo = if self.reaches(PI, TAU) { -1. } else { a.min(b) };
        Self::new((lo - LIBM).max(-1.), (hi + LIBM).min(1.))
    }

    pub fn atan(self) -> Self {
        Self::new(self.lo.atan() - LIBM, self.hi.atan() + LIBM)
    }

    /// The arc cosine of the part of the interval within `[-1, 1]`: the
    /// caller knows the exact argument lies there.
    pub fn acos(self) -> Self {
        let [lo, hi] = [self.lo.max(-1.), self.hi.min(1.)];
        if lo > hi {
            let end = if self.lo > 1. { 1_f64 } else { -1. };
            return Self::new(end.acos() - LIBM, end.acos() + LIBM);
        }
        Self::new(hi.acos() - LIBM, lo.acos() + LIBM)
    }

    /// The closed-form constant pi.
    pub fn pi() -> Self {
        Self::out(std::f64::consts::PI, std::f64::consts::PI)
    }
}

impl Add for Iv {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self::out(self.lo + o.lo, self.hi + o.hi)
    }
}

impl Sub for Iv {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self::out(self.lo - o.hi, self.hi - o.lo)
    }
}

impl Neg for Iv {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.hi, -self.lo)
    }
}

impl Mul for Iv {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        let p = [
            self.lo * o.lo,
            self.lo * o.hi,
            self.hi * o.lo,
            self.hi * o.hi,
        ];
        if p.iter().any(|x| x.is_nan()) {
            return ENTIRE;
        }
        let lo = p.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = p.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        Self::out(lo, hi)
    }
}

impl Div for Iv {
    type Output = Self;
    fn div(self, o: Self) -> Self {
        Mul::mul(self, o.recip())
    }
}

/// The arithmetic the charts are written in: intervals for enclosures and
/// second-order interval jets for derivative bounds.
pub(in crate::drawing) trait Num:
    Copy
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
{
    fn k(x: f64) -> Self;
    fn sin(self) -> Self;
    fn cos(self) -> Self;
    fn sqrt(self) -> Self;
    fn sqr(self) -> Self;
}

impl Num for Iv {
    fn k(x: f64) -> Self {
        Self::point(x)
    }
    fn sin(self) -> Self {
        Iv::sin(self)
    }
    fn cos(self) -> Self {
        Iv::cos(self)
    }
    fn sqrt(self) -> Self {
        Iv::sqrt(self)
    }
    fn sqr(self) -> Self {
        Iv::sqr(self)
    }
}

/// A function of the chart coordinates with its first and second partial
/// derivatives, each enclosed over a box.
#[derive(Clone, Copy, Debug)]
pub(in crate::drawing) struct J2 {
    pub v: Iv,
    pub d: [Iv; 2],
    /// `[xx, xy, yy]`.
    pub dd: [Iv; 3],
}

const ZERO: Iv = Iv { lo: 0., hi: 0. };

impl J2 {
    fn variable(value: Iv, index: usize) -> Self {
        let mut d = [ZERO; 2];
        d[index] = Iv::point(1.);
        Self {
            v: value,
            d,
            dd: [ZERO; 3],
        }
    }

    /// `f(self)` from `f`, `f'` and `f''` enclosed over the value.
    fn chain(self, f: Iv, f1: Iv, f2: Iv) -> Self {
        let [x, y] = self.d;
        let [xx, xy, yy] = self.dd;
        Self {
            v: f,
            d: [f1 * x, f1 * y],
            dd: [
                f2 * x.sqr() + f1 * xx,
                f2 * (x * y) + f1 * xy,
                f2 * y.sqr() + f1 * yy,
            ],
        }
    }
}

impl Add for J2 {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self {
            v: self.v + o.v,
            d: [self.d[0] + o.d[0], self.d[1] + o.d[1]],
            dd: [
                self.dd[0] + o.dd[0],
                self.dd[1] + o.dd[1],
                self.dd[2] + o.dd[2],
            ],
        }
    }
}

impl Neg for J2 {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            v: -self.v,
            d: self.d.map(|x| -x),
            dd: self.dd.map(|x| -x),
        }
    }
}

impl Sub for J2 {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        self + -o
    }
}

impl Mul for J2 {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        let (a, b) = (self, o);
        let two = Iv::point(2.);
        Self {
            v: a.v * b.v,
            d: [a.d[0] * b.v + a.v * b.d[0], a.d[1] * b.v + a.v * b.d[1]],
            dd: [
                a.dd[0] * b.v + two * (a.d[0] * b.d[0]) + a.v * b.dd[0],
                a.dd[1] * b.v + a.d[0] * b.d[1] + a.d[1] * b.d[0] + a.v * b.dd[1],
                a.dd[2] * b.v + two * (a.d[1] * b.d[1]) + a.v * b.dd[2],
            ],
        }
    }
}

impl Div for J2 {
    type Output = Self;
    fn div(self, o: Self) -> Self {
        let r = o.v.recip();
        let r2 = r.sqr();
        self * o.chain(r, -r2, Iv::point(2.) * r2 * r)
    }
}

impl Num for J2 {
    fn k(x: f64) -> Self {
        Self {
            v: Iv::point(x),
            d: [ZERO; 2],
            dd: [ZERO; 3],
        }
    }
    fn sin(self) -> Self {
        let (s, c) = (self.v.sin(), self.v.cos());
        self.chain(s, c, -s)
    }
    fn cos(self) -> Self {
        let (s, c) = (self.v.sin(), self.v.cos());
        self.chain(c, -s, -c)
    }
    fn sqrt(self) -> Self {
        let f = self.v.sqrt();
        let f1 = (Iv::point(2.) * f).recip();
        // f'' = -1 / (4 v^(3/2)) = -2 f'^3.
        let f2 = -(Iv::point(2.) * f1 * f1.sqr());
        self.chain(f, f1, f2)
    }
    fn sqr(self) -> Self {
        self.chain(self.v.sqr(), Iv::point(2.) * self.v, Iv::point(2.))
    }
}

/// The narrowest tried interval around `x` across which `f` provably
/// changes sign: it holds a root of `f`.
pub(in crate::drawing) fn bracket(x: f64, f: impl Fn(Iv) -> Iv) -> Option<Iv> {
    let mut w = 1e-15 * x.abs().max(1.);
    while w < 1e-7 {
        let [lo, hi] = [f(Iv::point(x - w)), f(Iv::point(x + w))];
        if lo.excludes_zero() && hi.excludes_zero() && (lo.lo > 0.) != (hi.lo > 0.) {
            return Some(Iv::new(x - w, x + w));
        }
        w *= 2.;
    }
    None
}

/// An enclosure of an exact point of an implicit chart's curve next to the
/// computed chart point `at`: one coordinate keeps its computed value and
/// the other is enclosed by a proven sign change of the relation, choosing
/// the coordinate whose enclosure places tighter.
///
/// A chord bound holds for any exact curve points next to the published
/// ones, not only those of the same parameter. Where the parameter is
/// ill-conditioned (near the axis a meridian line meets the profile at a
/// grazing angle, so the root moves by the angle's error times the radius
/// squared over the offset) the point at the computed value of the other
/// coordinate is enclosed tightly. Each stop gets one enclosure, which the
/// chords on both sides share, so their exact pieces still tile the curve.
pub(in crate::drawing) fn settle<C: Chart>(chart: &C, at: [f64; 2]) -> Option<[Iv; 2]> {
    let along = |axis: usize| {
        let fixed = at.map(Iv::point);
        let root = bracket(at[axis], |x| {
            let mut q = fixed;
            q[axis] = x;
            chart.relation(q[0], q[1])
        })?;
        let mut q = fixed;
        q[axis] = root;
        Some(q)
    };
    let size = |q: &[Iv; 2]| diameter(chart.place(q[0], q[1]));
    [along(1), along(0)]
        .into_iter()
        .flatten()
        .min_by(|p, q| size(p).total_cmp(&size(q)))
}

/// The view point of local coordinates in a placement.
pub(in crate::drawing) fn frame<N: Num>(origin: P3, axes: [P3; 3], [x, y, z]: [N; 3]) -> [N; 3] {
    std::array::from_fn(|i| {
        x * N::k(axes[0][i]) + y * N::k(axes[1][i]) + z * N::k(axes[2][i]) + N::k(origin[i])
    })
}

/// A curve in a two-coordinate chart.
pub(in crate::drawing) trait Chart {
    /// Whether the curve is a zero set of `relation`; otherwise it is the
    /// explicit graph `place(x, _)`.
    fn implicit(&self) -> bool;
    fn relation<N: Num>(&self, x: N, y: N) -> N;
    fn place<N: Num>(&self, x: N, y: N) -> [N; 3];
}

/// Supremum of the distance from a point to the points of a box.
fn reach(p: P3, b: [Iv; 3]) -> f64 {
    let squares: f64 = (0..3)
        .map(|i| {
            let d = (p[i] - b[i].lo).abs().max((b[i].hi - p[i]).abs());
            d * d
        })
        .sum();
    squares.sqrt() * (1. + 1e-12)
}

fn norm(v: [Iv; 3]) -> f64 {
    v.iter().map(|x| x.mag() * x.mag()).sum::<f64>().sqrt() * (1. + 1e-12)
}

fn diameter(b: [Iv; 3]) -> f64 {
    norm(b.map(|x| Iv::new(0., x.width())))
}

/// Proven bound of the deviation of the exact chart piece from the chord
/// between `points`, the published points of its ends. `ends` enclose
/// exact chart points next to them, shared with the neighbouring chords;
/// along `native` the piece is known to stay within their range. Infinite
/// when nothing is proven.
pub(in crate::drawing) fn chord_bound<C: Chart>(
    chart: &C,
    native: usize,
    ends: [[Iv; 2]; 2],
    points: [P3; 2],
) -> f64 {
    let place_iv = |q: [Iv; 2]| chart.place(q[0], q[1]);
    let eps = (0..2)
        .map(|i| reach(points[i], place_iv(ends[i])))
        .fold(0., f64::max);
    if !eps.is_finite() {
        return f64::INFINITY;
    }
    let hull = [0, 1].map(|axis| ends[0][axis].hull(ends[1][axis]));
    if !chart.implicit() {
        let x = J2::variable(hull[0], 0);
        let p = chart.place(x, J2::k(0.));
        let h = hull[0].width();
        let second = h * h / 8. * norm(p.map(|c| c.dd[0]));
        return second.min(diameter(p.map(|c| c.v))) * (1. + 1e-12) + eps;
    }
    let other = 1 - native;
    // On an edge the relation keeps one sign. Its terms can cancel (near a
    // branch end its slope along the edge does), so it is enclosed in the
    // mean-value form as well as directly.
    let edge = |at: f64| {
        let mut q = hull;
        q[other] = Iv::point(at);
        let jet = chart.relation(J2::variable(q[0], 0), J2::variable(q[1], 1));
        if jet.v.excludes_zero() {
            return true;
        }
        let mid = q.map(|x| Iv::point(0.5 * (x.lo + x.hi)));
        let centre = chart.relation(mid[0], mid[1]);
        (centre + jet.d[0] * (q[0] - mid[0]) + jet.d[1] * (q[1] - mid[1])).excludes_zero()
    };
    let span = hull[other].width();
    let scale = 1e-12 * (1. + hull[other].mag());
    // Where the relation hardly changes across the native coordinate's
    // hull, the edges need a pad of at least its change along the native
    // coordinate (and its rounding) over its slope across: near the axis the
    // ends' enclosures are tight in both coordinates but the relation is
    // nearly flat in the other one.
    let drift = {
        let jet = chart.relation(J2::variable(hull[0], 0), J2::variable(hull[1], 1));
        let mid = hull.map(|x| Iv::point(0.5 * (x.lo + x.hi)));
        let slope = jet.d[other];
        let least = slope.lo.abs().min(slope.hi.abs());
        let change =
            jet.d[native].mag() * hull[native].width() + chart.relation(mid[0], mid[1]).width();
        (slope.excludes_zero() && least > 0.).then(|| change / least)
    };
    let pads = [2., 0.5, 8.].map(|factor| factor * span + scale);
    let wide = drift.map_or([f64::NAN; 2], |reach| {
        [2., 8.].map(|factor| span + factor * reach + scale)
    });
    let Some(pad) = pads
        .into_iter()
        .chain(wide.into_iter().filter(|pad| pad.is_finite()))
        .find(|pad| edge(hull[other].lo - pad) && edge(hull[other].hi + pad))
    else {
        return f64::INFINITY;
    };
    let mut cell = hull;
    cell[other] = hull[other].widened(pad);
    let [x, y] = [J2::variable(cell[0], 0), J2::variable(cell[1], 1)];
    let g = chart.relation(x, y);
    let p = chart.place(x, y);
    let mut best = diameter(p.map(|c| c.v));
    for by in [native, other] {
        let of = 1 - by;
        if !g.d[of].excludes_zero() {
            continue;
        }
        // Second partials in (by, of) order.
        let pick = |dd: [Iv; 3]| {
            if by == 0 {
                [dd[0], dd[1], dd[2]]
            } else {
                [dd[2], dd[1], dd[0]]
            }
        };
        let [gpp, gpq, gqq] = pick(g.dd);
        let slope = -(g.d[by] / g.d[of]);
        let two = Iv::point(2.);
        let curve = -((gpp + two * gpq * slope + gqq * slope.sqr()) / g.d[of]);
        let second = p.map(|c| {
            let [ppp, ppq, pqq] = pick(c.dd);
            ppp + two * ppq * slope + pqq * slope.sqr() + c.d[of] * curve
        });
        let h = hull[by].width();
        best = best.min(h * h / 8. * norm(second));
    }
    best * (1. + 1e-12) + eps
}

/// A path whose chords have proven deviation bounds.
pub(in crate::drawing) trait Proven: Path {
    /// Parameters, from 0 to 1, between which `bound` applies.
    fn breaks(&self) -> Vec<f64> {
        vec![0., 1.]
    }

    /// Proven bound of the deviation of the path over `range`, inside one
    /// break interval, from the chord of its end points.
    fn bound(&self, range: [f64; 2]) -> f64;
}

/// Published chords with their bounds.
#[cfg(test)]
type Chords = Vec<([f64; 2], f64)>;

#[cfg(test)]
thread_local! {
    /// When set, every published chord with its bound.
    pub(in crate::drawing) static RECORD: std::cell::RefCell<Option<Chords>> =
        const { std::cell::RefCell::new(None) };
}

/// Parameters whose chords provably stay within half the tolerance: split
/// until each chord's bound is met, into as many parts as the bound
/// predicts. A stretch narrower than any feature the tolerance resolves is
/// accepted as it stands.
pub(in crate::drawing) fn certified(path: &dyn Proven, tolerance: f64) -> Vec<f64> {
    let goal = 0.5 * tolerance;
    let mut output = vec![0.];
    let breaks = path.breaks();
    let mut work = 0usize;
    for pair in breaks.windows(2) {
        if pair[1] <= pair[0] {
            continue;
        }
        let mut stack = vec![(pair[0], pair[1])];
        while let Some((a, b)) = stack.pop() {
            work += 1;
            let bound = path.bound([a, b]);
            let narrow = b - a <= 1e-12 || work > 1 << 20;
            if bound <= goal || narrow {
                #[cfg(test)]
                RECORD.with(|record| {
                    if let Some(record) = record.borrow_mut().as_mut() {
                        record.push(([a, b], bound));
                    }
                });
                output.push(b);
                continue;
            }
            let parts = if bound.is_finite() {
                ((bound / goal).sqrt() * 1.1).ceil().clamp(2., 16.) as usize
            } else {
                2
            };
            for i in (0..parts).rev() {
                let [s, e] = [i, i + 1].map(|j| a + (b - a) * j as f64 / parts as f64);
                stack.push((s, if i + 1 == parts { b } else { e }));
            }
        }
    }
    *output.last_mut().unwrap() = 1.;
    output
}

/// Counts of a sampled check of proven flattening.
#[cfg(test)]
#[derive(Debug, Default)]
pub(in crate::drawing) struct Tally {
    pub paths: usize,
    pub chords: usize,
    /// Chords whose bound does not meet half the tolerance.
    pub unproven: usize,
    /// Largest sampled deviation over its bound, and over half the tolerance.
    pub tightness: f64,
    pub usage: f64,
}

/// Flatten a path, and compare every chord's bound with the deviation of
/// densely sampled path points from the chord.
#[cfg(test)]
pub(in crate::drawing) fn tally(path: &dyn Path, tolerance: f64, tally: &mut Tally, label: &str) {
    RECORD.with(|record| *record.borrow_mut() = Some(Vec::new()));
    let stops = path.flatten(tolerance);
    let record = RECORD.with(|record| record.borrow_mut().take()).unwrap();
    if record.is_empty() {
        return;
    }
    assert_eq!(record.len() + 1, stops.len(), "{label}");
    tally.paths += 1;
    let goal = 0.5 * tolerance;
    for ([a, b], bound) in record {
        tally.chords += 1;
        if bound.is_nan() || bound > goal {
            tally.unproven += 1;
        }
        let [p, q] = [path.point(a), path.point(b)];
        let d: P3 = std::array::from_fn(|i| q[i] - p[i]);
        let length = d.iter().map(|x| x * x).sum::<f64>();
        let mut worst: f64 = 0.;
        for i in 1..64 {
            let m = path.point(a + (b - a) * f64::from(i) / 64.);
            let e: P3 = std::array::from_fn(|k| m[k] - p[k]);
            let s = if length > 0. {
                (e.iter().zip(d).map(|(x, y)| x * y).sum::<f64>() / length).clamp(0., 1.)
            } else {
                0.
            };
            let off = (0..3)
                .map(|k| (e[k] - s * d[k]).powi(2))
                .sum::<f64>()
                .sqrt();
            worst = worst.max(off);
        }
        assert!(
            worst <= bound + 1e-9,
            "{label}: chord [{a}, {b}] deviates {worst} beyond its bound {bound}"
        );
        assert!(
            worst <= goal + 1e-9,
            "{label}: chord [{a}, {b}] deviates {worst}"
        );
        if bound > 0. && bound.is_finite() {
            tally.tightness = tally.tightness.max(worst / bound);
        }
        tally.usage = tally.usage.max(worst / goal);
    }
}
