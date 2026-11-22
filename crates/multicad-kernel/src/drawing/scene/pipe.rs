//! Authored axial pipe threads (NPT flat crests, BSPT rounded crests): the
//! material is `0 < z < L`, `r < F(theta, z) = A(z - k theta) + tau z` with
//! the single-start helical rate `k = p / 2pi` and the radial taper `tau`.
//!
//! In each meridian half-plane the boundary is one fixed chain `W` of exact
//! lines and circular arcs in `(u, radius)` coordinates, translated along the
//! taper direction `(1, tau)` by `s = k theta + m p`: the point `(u, W(u))` of
//! the chain lies at `z = u + s`, `r = W(u) + tau s`. Equivalently
//! `A(phi) = W(phi) - tau phi` over one pitch, and `W(u + p) = W(u) + tau p`.
//! The chain is recomputed from the source's `spec()` formulas and refused
//! unless it agrees with the source's own certified `radial_bounds`.
//!
//! Exactness: rails are conical helices, cap curves follow the exact law,
//! sections are closed-form meets of each meridian line with the chain, and
//! contour curves are, for every angle, the exact meet of the contour
//! condition with the chain: a closed-form point on a flank line, the unique
//! root inside a proven monotone window on a crest or root arc, and the two
//! straight generators of a cone band. Every event (line material, sweep and
//! contour range ends, cap contour points) is isolated with a proven bound.
//! Contour and section branch paths flatten with proven per-chord bounds in
//! their charts (see `helix::proven`).
use super::helix::{
    Branch, Chart, Iv, Meridian, Num, Proven, RESOLVED, Range, Span, certified, chord_bound, frame,
    isolate, level_stop, settle,
};
use super::paths::{Lift, Segment, add, dot3, scaled};
use super::profile::merge;
use super::*;
use multicad_precision::PortableMath as _;
use multicad_rust_kernel::axial_pipe_thread::authored::{AuthoredSection, AuthoredThread};
use multicad_rust_kernel::{normal_sweep, spur_gear::Edge};
use std::f64::consts::{PI, TAU};

/// Placement, chain and bounds of one pipe thread member.
struct Law {
    origin: P3,
    axes: [P3; 3],
    pitch: f64,
    rate: f64,
    taper: f64,
    length: f64,
    /// The meridian chain `W` over `[0, pitch]`.
    spans: Vec<Span>,
    /// Spans that are cone bands: lines of slope exactly `taper` in `W`.
    cones: Vec<bool>,
    /// Kind of the rail at the start of each span, `None` where the chain
    /// continues straight through.
    rails: Vec<Option<Kind>>,
    /// Bounds of `F` over the body.
    minimum: f64,
    maximum: f64,
    /// Bounds of `|A'| = |W' - tau|` and of `|W''|` inside the spans.
    slope: f64,
    curvature: f64,
}

impl Law {
    fn view(&self, [x, y, z]: P3) -> P3 {
        std::array::from_fn(|i| {
            self.origin[i] + x * self.axes[0][i] + y * self.axes[1][i] + z * self.axes[2][i]
        })
    }

    fn local(&self, p: P3) -> P3 {
        let d = add(p, scaled(self.origin, -1.));
        self.axes.map(|axis| dot3(axis, d))
    }

    fn vector(&self, v: P3) -> P3 {
        self.axes.map(|axis| dot3(axis, v))
    }

    fn surface(&self, r: f64, theta: f64, z: f64) -> P3 {
        let (sin, cos) = theta.portable_sin_cos();
        self.view([r * cos, r * sin, z])
    }

    /// Span index and period base of an `A` phase.
    fn locate(&self, phi: f64) -> (usize, f64) {
        let local = phi.rem_euclid(self.pitch);
        let index = self
            .spans
            .iter()
            .position(|span| local <= span.range()[1])
            .unwrap_or(self.spans.len() - 1);
        (index, phi - local)
    }

    /// `A` and `A'` at a phase.
    fn a(&self, phi: f64) -> [f64; 2] {
        let (index, base) = self.locate(phi);
        let local = phi - base;
        let [w, slope] = self.spans[index].eval(local);
        [w - self.taper * local, slope - self.taper]
    }

    /// The boundary radius `F` at an angle and level.
    fn radius(&self, theta: f64, z: f64) -> f64 {
        self.a(z - self.rate * theta)[0] + self.taper * z
    }

    /// The span translated along the taper direction by `s`.
    fn translated(&self, span: &Span, s: f64) -> Span {
        let lift = self.taper * s;
        match *span {
            Span::Line { a, b } => Span::Line {
                a: [a[0] + s, a[1] + lift],
                b: [b[0] + s, b[1] + lift],
            },
            Span::Arc {
                center,
                radius,
                sign,
                range,
            } => Span::Arc {
                center: [center[0] + s, center[1] + lift],
                radius,
                sign,
                range: range.map(|x| x + s),
            },
        }
    }

    /// Copies `m` whose translation `k theta + m p` can reach the material
    /// for some angle in `[0, 2pi]`.
    fn copies(&self) -> std::ops::RangeInclusive<i64> {
        let turn = TAU * self.rate;
        let low = ((-turn - self.pitch) / self.pitch).floor() as i64 - 1;
        let high = (self.length / self.pitch).ceil() as i64 + 1;
        low..=high
    }

    /// Bound of the second derivative of a cap curve point in its angle.
    fn cap_curvature(&self) -> f64 {
        let k = self.rate;
        k * k * self.curvature + 2. * k * self.slope + self.maximum
    }
}

/// A rail: the conical helix of the chain point `(u, w)`, over `s` in range.
struct Rail {
    law: Arc<Law>,
    u: f64,
    w: f64,
    range: [f64; 2],
}

impl Rail {
    fn at(&self, s: f64) -> P3 {
        let law = &self.law;
        law.surface(self.w + law.taper * s, s / law.rate, self.u + s)
    }
}

impl Path for Rail {
    fn point(&self, t: f64) -> P3 {
        self.at(self.range[0] + (self.range[1] - self.range[0]) * t)
    }

    fn flatten(&self, tolerance_mm: f64) -> Vec<f64> {
        // In the angle, |P''| <= r + 2 |dr/dtheta| with dr/dtheta = tau k.
        let law = &self.law;
        let second = law.maximum + 2. * (law.taper * law.rate).abs();
        let step = (8. * tolerance_mm / second).sqrt().min(0.5);
        let sweep = (self.range[1] - self.range[0]).abs() / law.rate;
        let count = ((sweep / step).ceil() as usize).clamp(1, 1 << 20);
        (0..=count).map(|i| i as f64 / count as f64).collect()
    }

    fn paper(&self, _: [f64; 2]) -> Option<Curve> {
        None
    }
}

/// The cap curve at `level` over `theta` in `range`, inside one span of the
/// period based at `base`.
struct Cap {
    law: Arc<Law>,
    span: Span,
    base: f64,
    level: f64,
    range: [f64; 2],
}

impl Cap {
    fn at(&self, theta: f64) -> P3 {
        let law = &self.law;
        let local = self.level - law.rate * theta - self.base;
        let r = self.span.eval(local)[0] - law.taper * local + law.taper * self.level;
        law.surface(r, theta, self.level)
    }
}

impl Path for Cap {
    fn point(&self, t: f64) -> P3 {
        self.at(self.range[0] + (self.range[1] - self.range[0]) * t)
    }

    fn flatten(&self, tolerance_mm: f64) -> Vec<f64> {
        let step = (8. * tolerance_mm / self.law.cap_curvature())
            .sqrt()
            .min(0.5);
        let count = ((self.range[1] - self.range[0]).abs() / step).ceil();
        let count = (count as usize).clamp(1, 1 << 20);
        (0..=count).map(|i| i as f64 / count as f64).collect()
    }

    fn paper(&self, _: [f64; 2]) -> Option<Curve> {
        None
    }
}

/// The view direction's terms of the contour condition at one angle:
/// `n . v = r (P - W' v_z) - k (W' - tau) S` with `P = rho cos psi`,
/// `S = rho sin psi`, `psi = theta - alpha`.
#[derive(Clone, Copy)]
struct Sight {
    rho: f64,
    alpha: f64,
    vz: f64,
}

impl Sight {
    fn new(view: P3) -> Self {
        Self {
            rho: view[0].portable_hypot(view[1]),
            alpha: view[1].portable_atan2(view[0]),
            vz: view[2],
        }
    }

    fn terms(&self, theta: f64) -> [f64; 2] {
        let (sin, cos) = (theta - self.alpha).portable_sin_cos();
        [self.rho * cos, self.rho * sin]
    }
}

/// One contour branch on a non-cone span, parametrised by the angle: for each
/// `theta` the exact contour point of the meridian chain at `s = k theta`.
struct Contour {
    law: Arc<Law>,
    /// The span of the base period; the branch covers every translation.
    span: Span,
    /// On an arc, the root window: `+1` or `-1` side of the dominant root.
    window: f64,
    sight: Sight,
    range: [f64; 2],
}

impl Contour {
    /// The chain parameter `u` of the contour point at an angle, if any.
    fn u(law: &Law, span: &Span, window: f64, sight: &Sight, theta: f64) -> Option<f64> {
        let s = law.rate * theta;
        let [p, q] = sight.terms(theta);
        let (k, tau, vz) = (law.rate, law.taper, sight.vz);
        match *span {
            Span::Line { a, b } => {
                let w = (b[1] - a[1]) / (b[0] - a[0]);
                let denominator = p - w * vz;
                if denominator.abs() <= 1e-300 {
                    return None;
                }
                // (W(u) + tau s)(P - w v_z) = k (w - tau) S.
                let target = k * (w - tau) * q / denominator - tau * s;
                Some(a[0] + (target - a[1]) / w)
            }
            Span::Arc {
                center,
                radius,
                sign,
                ..
            } => {
                let beta = arc_root(law, center, radius, sign, window, sight, theta)?;
                Some(center[0] + radius * beta.portable_sin())
            }
        }
    }

    fn at(&self, theta: f64) -> P3 {
        let law = &self.law;
        let s = law.rate * theta;
        let u = Self::u(law, &self.span, self.window, &self.sight, theta).unwrap_or(f64::NAN);
        let w = self
            .span
            .eval(u.clamp(self.span.range()[0], self.span.range()[1]))[0];
        law.surface(w + law.taper * s, theta, u + s)
    }
}

impl Contour {
    fn theta(&self, t: f64) -> f64 {
        self.range[0] + (self.range[1] - self.range[0]) * t
    }

    /// The exact chart point at a parameter: on an arc `(theta, beta)`, a
    /// root next to the window's computed root, enclosed by a proven sign
    /// change.
    fn end(&self, t: f64) -> Option<[Iv; 2]> {
        let theta = self.theta(t);
        let Span::Arc {
            center,
            radius,
            sign,
            ..
        } = self.span
        else {
            return Some([Iv::point(theta), Iv::point(0.)]);
        };
        let law = &self.law;
        let beta = arc_root(law, center, radius, sign, self.window, &self.sight, theta)?;
        settle(self, [theta, beta])
    }
}

/// On a flank line the contour is the closed-form graph of the angle; on an
/// arc it is the window's root of the arc's condition in `(theta, beta)`,
/// with `u = c0 + R sin beta` and `W = c1 + sign R cos beta`.
impl Chart for Contour {
    fn implicit(&self) -> bool {
        matches!(self.span, Span::Arc { .. })
    }

    fn relation<N: Num>(&self, theta: N, beta: N) -> N {
        let Span::Arc {
            center,
            radius,
            sign,
            ..
        } = self.span
        else {
            return N::k(0.);
        };
        let (law, sight) = (&self.law, &self.sight);
        let psi = theta - N::k(sight.alpha);
        let [p, q] = [N::k(sight.rho) * psi.cos(), N::k(sight.rho) * psi.sin()];
        let c = N::k(center[1]) + N::k(law.taper * law.rate) * theta;
        let s = N::k(law.rate) * q;
        let e = c * p + N::k(law.taper) * s;
        let f = N::k(sign) * (c * N::k(sight.vz) + s);
        let (sin, cos) = (beta.sin(), beta.cos());
        e * cos + f * sin + N::k(sign * radius) * cos * (p * cos + N::k(sign * sight.vz) * sin)
    }

    fn place<N: Num>(&self, theta: N, beta: N) -> [N; 3] {
        let law = &self.law;
        let s = N::k(law.rate) * theta;
        let (u, w) = match self.span {
            Span::Line { a, b } => {
                let w = (b[1] - a[1]) / (b[0] - a[0]);
                let sight = &self.sight;
                let psi = theta - N::k(sight.alpha);
                let [p, q] = [N::k(sight.rho) * psi.cos(), N::k(sight.rho) * psi.sin()];
                let target = N::k(law.rate * (w - law.taper)) * q / (p - N::k(w * sight.vz))
                    - N::k(law.taper) * s;
                let u = N::k(a[0]) + (target - N::k(a[1])) / N::k(w);
                (u, self.span.radius_at(u)[0])
            }
            Span::Arc {
                center,
                radius,
                sign,
                ..
            } => (
                N::k(center[0]) + N::k(radius) * beta.sin(),
                N::k(center[1]) + N::k(sign * radius) * beta.cos(),
            ),
        };
        let r = w + N::k(law.taper) * s;
        frame(
            law.origin,
            law.axes,
            [r * theta.cos(), r * theta.sin(), u + s],
        )
    }
}

impl Proven for Contour {
    fn bound(&self, [a, b]: [f64; 2]) -> f64 {
        let (Some(start), Some(end)) = (self.end(a), self.end(b)) else {
            return f64::INFINITY;
        };
        chord_bound(self, 0, [start, end], [self.point(a), self.point(b)])
    }
}

impl Path for Contour {
    fn point(&self, t: f64) -> P3 {
        self.at(self.theta(t))
    }

    fn flatten(&self, tolerance_mm: f64) -> Vec<f64> {
        certified(self, tolerance_mm)
    }

    fn paper(&self, _: [f64; 2]) -> Option<Curve> {
        None
    }
}

/// The arc terms of the contour condition at an angle, in the circle angle
/// `beta` with `u = c0 + R sin beta`, `W = c1 + sign R cos beta`:
/// `f(beta) = E cos beta + F sin beta + sign R cos beta (P cos beta + sign v_z
/// sin beta)` is `cos beta` times `n . v`.
fn arc_terms(law: &Law, c1: f64, sign: f64, sight: &Sight, theta: f64) -> [f64; 4] {
    let [p, q] = sight.terms(theta);
    let c = c1 + law.taper * law.rate * theta;
    let s = law.rate * q;
    [
        c * p + law.taper * s,
        sign * (c * sight.vz + s),
        p,
        sight.vz,
    ]
}

fn arc_value(radius: f64, sign: f64, [e, f, p, vz]: [f64; 4], beta: f64) -> f64 {
    let (sin, cos) = beta.portable_sin_cos();
    e * cos + f * sin + sign * radius * cos * (p * cos + sign * vz * sin)
}

/// The window root of an arc's contour condition at an angle. The dominant
/// part `N cos(beta - beta0)`, `N = |(E, F)|`, vanishes at `beta0 +- pi/2`;
/// the perturbation and its derivative are at most `R |(P, v_z)| = eps N`.
/// With `eps < 0.4`, within `asin(2 eps)` of each zero the dominant part
/// exceeds the perturbation at the window ends and its slope exceeds the
/// perturbation's (`sqrt(1 - 4 eps^2) > eps`), so each window holds exactly
/// one simple root. Returns it when it lies within a quarter turn.
fn arc_root(
    law: &Law,
    center: P2,
    radius: f64,
    sign: f64,
    window: f64,
    sight: &Sight,
    theta: f64,
) -> Option<f64> {
    let terms = arc_terms(law, center[1], sign, sight, theta);
    let [e, f, p, vz] = terms;
    let norm = e.portable_hypot(f);
    let epsilon = radius * p.portable_hypot(vz) / norm;
    if epsilon.is_nan() || epsilon >= WINDOW {
        return None;
    }
    let half = (2. * epsilon).portable_asin();
    // The window center, as the representative nearest the arc's angles.
    let mut middle = f.portable_atan2(e) + window * PI / 2.;
    middle -= TAU * (middle / TAU).round();
    let [mut lo, mut hi] = [middle - half, middle + half];
    let [mut flo, fhi] = [lo, hi].map(|beta| arc_value(radius, sign, terms, beta));
    if (flo > 0.) == (fhi > 0.) {
        // Rounding at a window end that is itself the root.
        return [lo, hi]
            .into_iter()
            .find(|beta| arc_value(radius, sign, terms, *beta) == 0.)
            .filter(|beta| beta.abs() < PI / 2.);
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if mid <= lo || mid >= hi {
            break;
        }
        let fm = arc_value(radius, sign, terms, mid);
        if (fm > 0.) == (flo > 0.) {
            lo = mid;
            flo = fm;
        } else {
            hi = mid;
        }
    }
    let beta = 0.5 * (lo + hi);
    (beta.abs() < PI / 2.).then_some(beta)
}

/// The largest admitted ratio `R |(P, v_z)| / N` of an arc's contour
/// condition (see `arc_root`).
const WINDOW: f64 = 0.4;

/// A branch of an oblique section on a translated span.
struct Sweep {
    law: Arc<Law>,
    /// The chart's angle origin: the branch's angles are measured from it.
    base: f64,
    branch: Branch,
}

impl Sweep {
    /// The plane's meridian line `alpha W + beta u = gamma` in the chain
    /// coordinates of the copy translated by `k theta` more.
    fn line(law: &Law, (n, c): (P3, f64), theta: f64) -> [f64; 3] {
        Self::meridian(law, (n, c)).line(theta)
    }

    fn meridian(law: &Law, (n, c): (P3, f64)) -> Meridian {
        Meridian {
            n,
            c,
            k: law.rate,
            tau: law.taper,
        }
    }

    /// The chain point of a local point at the chart angle `theta`.
    fn profile(&self, p: P3) -> impl Fn(f64) -> P2 + '_ {
        move |theta| {
            let s = self.law.rate * theta;
            [p[2] - s, p[0].portable_hypot(p[1]) - self.law.taper * s]
        }
    }
}

impl Proven for Sweep {
    fn breaks(&self) -> Vec<f64> {
        self.branch.breaks()
    }

    fn bound(&self, [a, b]: [f64; 2]) -> f64 {
        let placement = (self.law.origin, self.law.axes, self.base);
        self.branch
            .bound(placement, [a, b], [self.point(a), self.point(b)])
    }
}

impl Path for Sweep {
    fn point(&self, t: f64) -> P3 {
        let (theta, p) = self.branch.locate(t);
        let law = &self.law;
        let s = law.rate * theta;
        law.surface(p[1] + law.taper * s, self.base + theta, p[0] + s)
    }

    fn flatten(&self, tolerance_mm: f64) -> Vec<f64> {
        certified(self, tolerance_mm)
    }

    fn paper(&self, _: [f64; 2]) -> Option<Curve> {
        None
    }
}

/// The pipe thread member.
pub(in crate::drawing) struct Pipe {
    law: Arc<Law>,
    tolerance: f64,
}

pub(in crate::drawing) fn member(
    solid: &Solid,
    view: &View,
) -> Option<Result<Box<dyn Member>, GeometryError>> {
    let (source, occurrence) = solid.rust_authored_axial_pipe_thread_source()?;
    let tolerance = solid.geometry_policy().linear_tolerance_mm;
    let construction = source.construction_frame;
    let place = |p: [f64; 3], point: bool| -> [f64; 3] {
        let n = occurrence.normal();
        std::array::from_fn(|i| {
            (if point { occurrence.origin[i] } else { 0. })
                + p[0] * occurrence.u[i]
                + p[1] * occurrence.v[i]
                + p[2] * n[i]
        })
    };
    let origin = view.point(place(construction.origin, true));
    let axes = [construction.u, construction.v, construction.normal()]
        .map(|axis| view.vector(place(axis, false)));
    Some(chain(&source.transport).map(|(spans, cones)| {
        Box::new(Pipe::new(
            &source.transport,
            origin,
            axes,
            spans,
            cones,
            tolerance,
        )) as Box<dyn Member>
    }))
}

/// The meridian chain `W` over one pitch from the authored `spec()` formulas,
/// refused unless it lies within the source's own `radial_bounds`.
fn chain(thread: &AuthoredThread) -> Result<(Vec<Span>, Vec<bool>), GeometryError> {
    let spec = thread.spec();
    let (r, p, l, tau) = (
        spec.major_radius_mm,
        spec.pitch_mm,
        spec.length_mm,
        spec.radial_taper,
    );
    let line = |a: P2, b: P2| Span::Line { a, b };
    let (spans, cones) = match spec.section {
        AuthoredSection::NptFlat {
            depth_mm: d,
            flat_width_mm: width,
        } => {
            let h = width / 2.;
            let m = 3_f64.sqrt();
            let down = h + d / (m + tau);
            let up = p - h - d / (m - tau);
            (
                vec![
                    line([0., r], [h, r + tau * h]),
                    line([h, r + tau * h], [down, r - d + tau * down]),
                    line([down, r - d + tau * down], [up, r - d + tau * up]),
                    line([up, r - d + tau * up], [p - h, r + tau * (p - h)]),
                    line([p - h, r + tau * (p - h)], [p, r + tau * p]),
                ],
                vec![true, false, true, false, true],
            )
        }
        AuthoredSection::BsptRounded { depth_mm: d } => {
            let a = 27.5_f64.to_radians();
            let (sin, cos) = a.portable_sin_cos();
            let m = cos / sin;
            let normal = (1. + tau * tau).sqrt();
            let triangle = p * (m * m - tau * tau) / (2. * m);
            let rho = (triangle - d) / (2. * (1. / sin - normal));
            let valley = triangle / (m + tau);
            let h = rho * cos;
            let crest = r - rho * normal;
            let tangent = crest + rho * sin;
            let root = r - d + rho * normal + tau * valley;
            let arc = |center: P2, sign: f64, range: [f64; 2]| Span::Arc {
                center,
                radius: rho,
                sign,
                range,
            };
            (
                vec![
                    arc([0., crest], 1., [0., h]),
                    line([h, tangent], [valley - h, tangent - m * (valley - 2. * h)]),
                    arc([valley, root], -1., [valley - h, valley + h]),
                    line(
                        [valley + h, tangent + tau * p - m * (p - valley - 2. * h)],
                        [p - h, tangent + tau * p],
                    ),
                    arc([p, crest + tau * p], 1., [p - h, p]),
                ],
                vec![false; 5],
            )
        }
    };
    let invalid = |message: &'static str| err(GeometryErrorKind::InvalidResult, message);
    // The chain must be continuous, close over one pitch with W(p) = W(0) +
    // tau p, and agree with the certified radial bounds of the source.
    let scale = r + p + l + 1.;
    let slack = 256. * f64::EPSILON * scale;
    let w = |span: &Span, u: f64| span.eval(u)[0];
    if spans[0].range()[0] != 0. || spans[4].range()[1] != p {
        return Err(invalid("Pipe thread chain does not span one pitch"));
    }
    for i in 0..spans.len() {
        let end = spans[i].range()[1];
        let (next, shift) = if i + 1 < spans.len() {
            (&spans[i + 1], 0.)
        } else {
            (&spans[0], p)
        };
        let start = next.range()[0] + shift;
        if (end - start).abs() > slack
            || (w(&spans[i], end) - (w(next, next.range()[0]) + tau * shift)).abs() > 64. * slack
        {
            return Err(invalid("Pipe thread chain is not continuous"));
        }
    }
    if spec
        .crest_radius_limit_mm
        .is_some_and(|limit| limit <= r + (tau * l).max(0.) + slack)
    {
        return Err(unsupported(
            "Drawing a pipe thread whose crest radius limit trims its crests is not supported",
        ));
    }
    let samples = 96;
    for span in &spans {
        let [a, b] = span.range();
        for i in 0..=samples {
            let u = a + (b - a) * i as f64 / samples as f64;
            for z in [0., 0.5 * l, l] {
                let mine = w(span, u) - tau * u + tau * z;
                let [lo, hi] = thread
                    .radial_bounds([u, u], [z, z])
                    .map_err(|_| invalid("Pipe thread radial bounds cannot be read"))?;
                if mine < lo - slack || mine > hi + slack {
                    return Err(invalid(
                        "Pipe thread chain disagrees with the source's radial bounds",
                    ));
                }
            }
        }
    }
    Ok((spans, cones))
}

impl Pipe {
    fn new(
        thread: &AuthoredThread,
        origin: P3,
        axes: [P3; 3],
        spans: Vec<Span>,
        cones: Vec<bool>,
        tolerance: f64,
    ) -> Self {
        let spec = thread.spec();
        let (p, tau, l) = (spec.pitch_mm, spec.radial_taper, spec.length_mm);
        let count = spans.len();
        let rails = (0..count)
            .map(|i| {
                let before = spans[(i + count - 1) % count];
                let after = spans[i];
                if cones[i] && cones[(i + count - 1) % count] {
                    return None;
                }
                let [_, s0] = before.eval(before.range()[1]);
                let [_, s1] = after.eval(after.range()[0]);
                let scale = 1. + s0.abs().max(s1.abs());
                Some(if (s0 - s1).abs() <= 1e-7 * scale {
                    Kind::Smooth
                } else {
                    Kind::Edge
                })
            })
            .collect();
        let mut low = f64::INFINITY;
        let mut high = f64::NEG_INFINITY;
        let mut slope: f64 = 0.;
        let mut curvature: f64 = 0.;
        for span in &spans {
            let [a, b] = span.range();
            let mut phases = vec![a, b];
            if let Span::Arc { center, radius, .. } = *span {
                let x = (a - center[0]).abs().max((b - center[0]).abs());
                let h = ((radius - x) * (radius + x)).max(0.).sqrt();
                curvature = curvature.max(radius * radius / h.max(1e-300).powi(3));
                // A' = W' - tau vanishes inside an arc where W' = tau.
                let x = -tau * radius / (1. + tau * tau).sqrt();
                let sign = match *span {
                    Span::Arc { sign, .. } => sign,
                    _ => 1.,
                };
                let phi = center[0] + sign * x;
                if phi > a && phi < b {
                    phases.push(phi);
                }
            }
            for phi in phases {
                let [w, s] = span.eval(phi);
                low = low.min(w - tau * phi);
                high = high.max(w - tau * phi);
                if s.is_finite() {
                    slope = slope.max((s - tau).abs());
                }
            }
        }
        let law = Law {
            origin,
            axes,
            pitch: p,
            rate: p / TAU,
            taper: tau,
            length: l,
            spans,
            cones,
            rails,
            minimum: (low + (tau * l).min(0.)) * (1. - 1e-12),
            maximum: (high + (tau * l).max(0.)) * (1. + 1e-12),
            slope: slope * (1. + 1e-12) + 1e-12,
            curvature: curvature * (1. + 1e-12) + 1e-12,
        };
        Self {
            law: Arc::new(law),
            tolerance,
        }
    }
}

impl Pipe {
    /// Material parameters along a unit-direction local line.
    fn line_local(&self, s: P3, d: P3) -> Result<Vec<[f64; 2]>, GeometryError> {
        let law = &self.law;
        let slab = if d[2].abs() <= 1e-15 {
            if s[2] <= 0. || s[2] >= law.length {
                return Ok(Vec::new());
            }
            [-1e9, 1e9]
        } else {
            let ends = [0., law.length].map(|z| (z - s[2]) / d[2]);
            [ends[0].min(ends[1]), ends[0].max(ends[1])]
        };
        let rho = d[0].portable_hypot(d[1]);
        let within = |radius: f64| -> Option<[f64; 2]> {
            if rho <= 1e-15 {
                return (s[0].portable_hypot(s[1]) < radius)
                    .then_some([f64::NEG_INFINITY, f64::INFINITY]);
            }
            let a = rho * rho;
            let b = s[0] * d[0] + s[1] * d[1];
            let c = s[0] * s[0] + s[1] * s[1] - radius * radius;
            let disc = b * b - a * c;
            (disc > 0.).then(|| {
                let root = disc.sqrt();
                [(-b - root) / a, (-b + root) / a]
            })
        };
        let core = within(law.minimum);
        let Some(outer) = within(law.maximum) else {
            return Ok(Vec::new());
        };
        let span = [slab[0].max(outer[0]), slab[1].min(outer[1])];
        if span[1] <= span[0] {
            return Ok(Vec::new());
        }
        let g = |t: f64| {
            let p = add(s, scaled(d, t));
            law.radius(p[1].portable_atan2(p[0]), p[2]) - p[0].portable_hypot(p[1])
        };
        let mut stretches = vec![span];
        let mut output = Vec::new();
        if let Some([a, b]) = core {
            let [a, b] = [a.max(span[0]), b.min(span[1])];
            if b > a {
                // Inside the minimum radius the line is in material.
                stretches = vec![[span[0], a], [b, span[1]]];
                output.push([a, b]);
            }
        }
        for [a, b] in stretches {
            if b <= a {
                continue;
            }
            let roots = self.crossings(s, d, [a, b])?;
            let mut stops = vec![a];
            stops.extend(roots.into_iter().filter(|t| *t > a && *t < b));
            stops.push(b);
            for pair in stops.windows(2) {
                if pair[1] > pair[0] && g(0.5 * (pair[0] + pair[1])) > 0. {
                    output.push([pair[0], pair[1]]);
                }
            }
        }
        Ok(merge(output, 0.))
    }

    /// Parameters in `[a, b]`, where the line's radius is at least the
    /// minimum of `F`, at which the unit-direction line crosses the thread
    /// surface: `g = A(phi) + tau z - r` with `|g'|` and `|g''|` bounded.
    fn crossings(&self, s: P3, d: P3, [a, b]: [f64; 2]) -> Result<Vec<f64>, GeometryError> {
        let law = &self.law;
        let k = law.rate;
        let rho = d[0].portable_hypot(d[1]);
        if rho <= 1e-15 {
            return Ok(self.axial_crossings(s, d, [a, b]));
        }
        let g = |t: f64| {
            let p = add(s, scaled(d, t));
            law.radius(p[1].portable_atan2(p[0]), p[2]) - p[0].portable_hypot(p[1])
        };
        let slope = |t: f64| {
            let p = add(s, scaled(d, t));
            let r2 = p[0] * p[0] + p[1] * p[1];
            let phi = p[2] - k * p[1].portable_atan2(p[0]);
            let (index, base) = law.locate(phi);
            let a = law.a(phi)[1];
            let dphi = d[2] - k * (p[0] * d[1] - p[1] * d[0]) / r2;
            (
                a * dphi + law.taper * d[2] - (p[0] * d[0] + p[1] * d[1]) / r2.sqrt(),
                phi,
                index,
                base,
            )
        };
        let phase_rate = d[2].abs() + k * rho / law.minimum;
        let lipschitz =
            (law.slope * phase_rate + law.taper.abs() * d[2].abs() + rho) * (1. + 1e-6) + 1e-12;
        let second = law.curvature * phase_rate * phase_rate
            + law.slope * k * 2. * rho * rho / (law.minimum * law.minimum)
            + rho * rho / law.minimum;
        let second = second * (1. + 1e-6) + 1e-12;
        let width = 1e-13 * a.abs().max(b.abs()).max(1.);
        let mut output = Vec::new();
        let mut pending = vec![(a, g(a), b, g(b))];
        let mut work = 0usize;
        while let Some((a, fa, b, fb)) = pending.pop() {
            work += 1;
            if work > 4_000_000 {
                return Err(err(
                    GeometryErrorKind::Limit,
                    "Pipe thread drawing root isolation exceeded its work budget",
                ));
            }
            let h = b - a;
            let change = (fa <= 0.) != (fb <= 0.);
            if !change && fa.abs() + fb.abs() > lipschitz * h {
                continue;
            }
            // Second order, from either end, within one span of the chain.
            let excluded = [(a, fa), (b, fb)].into_iter().any(|(t, value)| {
                let (derivative, phi, index, base) = slope(t);
                let [lo, hi] = law.spans[index].range().map(|x| x + base);
                let room = (phi - lo).min(hi - phi);
                phase_rate * h < room && value.abs() > derivative.abs() * h + 0.5 * second * h * h
            });
            if excluded && !change {
                continue;
            }
            if h <= width || lipschitz * h <= RESOLVED {
                if change {
                    output.push(0.5 * (a + b));
                }
                continue;
            }
            let m = 0.5 * (a + b);
            let fm = g(m);
            pending.push((a, fa, m, fm));
            pending.push((m, fm, b, fb));
        }
        output.sort_by(f64::total_cmp);
        output.dedup_by(|x, y| (*x - *y).abs() <= 4. * width);
        Ok(output)
    }

    /// Exact crossings of a line parallel to the axis: on each translated
    /// span, the levels where the chain reaches the line's radius.
    fn axial_crossings(&self, s: P3, d: P3, [a, b]: [f64; 2]) -> Vec<f64> {
        let law = &self.law;
        let r = s[0].portable_hypot(s[1]);
        let theta = s[1].portable_atan2(s[0]);
        let [p, q] = [s[2] + a * d[2], s[2] + b * d[2]];
        let [lo, hi] = [p.min(q), p.max(q)];
        let shift = law.rate * theta;
        let first = ((lo - shift) / law.pitch).floor() as i64 - 2;
        let last = ((hi - shift) / law.pitch).ceil() as i64 + 1;
        let mut output = Vec::new();
        for m in first..=last {
            for span in &law.spans {
                let span = law.translated(span, shift + m as f64 * law.pitch);
                let [ra, rb] = span.range();
                let mut push = |z: f64| {
                    if z >= ra && z <= rb && z > lo && z < hi {
                        output.push((z - s[2]) / d[2]);
                    }
                };
                match span {
                    Span::Line { a, b } => {
                        if (a[1] - r) * (b[1] - r) < 0. {
                            push(a[0] + (b[0] - a[0]) * (r - a[1]) / (b[1] - a[1]));
                        }
                    }
                    Span::Arc {
                        center,
                        radius,
                        sign,
                        ..
                    } => {
                        let dr = r - center[1];
                        if sign * dr > 0. && dr.abs() < radius {
                            let w = ((radius - dr.abs()) * (radius + dr.abs())).sqrt();
                            push(center[0] - w);
                            push(center[0] + w);
                        }
                    }
                }
            }
        }
        output.sort_by(f64::total_cmp);
        output.dedup_by(|x, y| (*x - *y).abs() <= 1e-13);
        output
    }

    fn own_line(&self, start: P3, direction: P3) -> Result<Vec<[f64; 2]>, GeometryError> {
        let s = self.law.local(start);
        let d = self.law.vector(direction);
        let length = dot3(d, d).sqrt();
        if length <= 1e-300 {
            return Ok(Vec::new());
        }
        let unit = scaled(d, 1. / length);
        Ok(self
            .line_local(s, unit)?
            .into_iter()
            .map(|[a, b]| [a / length, b / length])
            .collect())
    }

    fn view_vector(&self, [x, y, z]: P3) -> P3 {
        let axes = self.law.axes;
        std::array::from_fn(|i| x * axes[0][i] + y * axes[1][i] + z * axes[2][i])
    }

    /// The rail or stationary helix of chain point `u`, split at quarter
    /// turns and where its radius turns perpendicular to the sight.
    fn rail(&self, u: f64, kind: Kind, output: &mut Vec<Candidate>) {
        let law = &self.law;
        let (index, base) = law.locate(u);
        let w = law.spans[index].eval(u - base)[0] + law.taper * base;
        let [a, b] = [-u, law.length - u];
        let k = law.rate;
        let view = law.vector([0., 0., 1.]);
        let mut stops = vec![a, b];
        if view[0].portable_hypot(view[1]) > 1e-12 {
            let fold = view[1].portable_atan2(view[0]) + PI / 2.;
            let mut j = ((a / k - fold) / PI).ceil();
            while k * (fold + j * PI) < b {
                stops.push(k * (fold + j * PI));
                j += 1.;
            }
        }
        let quarter = k * PI / 2.;
        let count = ((b - a) / quarter).ceil().max(1.) as usize;
        stops.extend((1..count).map(|i| a + (b - a) * i as f64 / count as f64));
        stops.sort_by(f64::total_cmp);
        stops.dedup_by(|x, y| (*x - *y).abs() <= 1e-12);
        for pair in stops.windows(2) {
            output.push(Candidate {
                kind,
                path: Arc::new(Rail {
                    law: law.clone(),
                    u,
                    w,
                    range: [pair[0], pair[1]],
                }),
            });
        }
    }

    /// Pieces `(span, base, angles)` of the cap curve at `level` over
    /// `[0, 2pi]`, each inside one span.
    fn cap_pieces(&self, level: f64) -> Vec<(Span, f64, [f64; 2])> {
        let law = &self.law;
        let k = law.rate;
        let mut cuts = vec![0., TAU];
        for m in law.copies() {
            for span in &law.spans {
                let theta = (level - span.range()[0] - m as f64 * law.pitch) / k;
                if theta > 0. && theta < TAU {
                    cuts.push(theta);
                }
            }
        }
        cuts.sort_by(f64::total_cmp);
        cuts.dedup_by(|a, b| (*a - *b).abs() <= 1e-14);
        cuts.windows(2)
            .filter(|pair| pair[1] - pair[0] > 1e-14)
            .map(|pair| {
                let (index, base) = law.locate(level - k * 0.5 * (pair[0] + pair[1]));
                (law.spans[index], base, [pair[0], pair[1]])
            })
            .collect()
    }

    fn cap(&self, level: f64) -> Vec<Arc<dyn Path>> {
        let mut output: Vec<Arc<dyn Path>> = Vec::new();
        for (span, base, [a, b]) in self.cap_pieces(level) {
            let count = ((b - a) / (PI / 2.)).ceil().max(1.) as usize;
            for i in 0..count {
                output.push(Arc::new(Cap {
                    law: self.law.clone(),
                    span,
                    base,
                    level,
                    range: [
                        a + (b - a) * i as f64 / count as f64,
                        a + (b - a) * (i + 1) as f64 / count as f64,
                    ],
                }));
            }
        }
        output
    }

    /// Angles on each cap curve where the surface normal is perpendicular to
    /// the sight: there contour curves end on a cap. Isolated with a proven
    /// Lipschitz bound of `n . v` along the cap curve.
    fn contour_on_caps(&self, sight: &Sight) -> Result<Vec<(f64, f64)>, GeometryError> {
        let law = &self.law;
        let (k, tau) = (law.rate, law.taper);
        let (rmax, slope, curvature) = (law.maximum, law.slope, law.curvature);
        let (rho, vz) = (sight.rho, sight.vz.abs());
        let lipschitz = (k * slope * rho + rmax * rho)
            + (k * k * curvature * rho + k * slope * rho)
            + vz * (k * slope * (slope + tau.abs()) + rmax * k * curvature)
            + 1e-12;
        let mut output = Vec::new();
        for level in [0., law.length] {
            for (span, base, range) in self.cap_pieces(level) {
                let g = |theta: f64| {
                    let local = level - k * theta - base;
                    let [w, ws] = span.eval(local);
                    let r = w - tau * local + tau * level;
                    let [p, q] = sight.terms(theta);
                    r * p - k * (ws - tau) * q - r * ws * sight.vz
                };
                let mut roots = Vec::new();
                isolate(&g, range, lipschitz, &mut roots)?;
                output.extend(roots.into_iter().map(|theta| (level, theta)));
            }
        }
        Ok(output)
    }

    fn outlines(&self, output: &mut Vec<Candidate>) -> Result<(), GeometryError> {
        let law = &self.law;
        let view = law.vector([0., 0., 1.]);
        let sight = Sight::new(view);
        if sight.rho <= 1e-9 {
            // Along the axis the outlines are the helices of stationary W.
            for span in &law.spans {
                if let Span::Arc { .. } = span
                    && let Some(u) = span.stationary()
                    && u < law.pitch
                {
                    self.rail(u, Kind::Outline, output);
                }
            }
            return Ok(());
        }
        let ends = self.contour_on_caps(&sight)?;
        let k = law.rate;
        for (index, span) in law.spans.iter().enumerate() {
            let [a, b] = span.range();
            if law.cones[index] {
                self.cone_generators(span, &sight, output);
                continue;
            }
            let range = [-b / k, (law.length - a) / k];
            let windows: &[f64] = match span {
                Span::Line { .. } => &[0.],
                Span::Arc { .. } => &[-1., 1.],
            };
            let mut events = vec![range[0], range[1]];
            match *span {
                Span::Line { a: pa, b: pb } => {
                    let w = (pb[1] - pa[1]) / (pb[0] - pa[0]);
                    // The contour meets a rail where n . v vanishes there.
                    let reach = law.taper.abs() * (range[0].abs().max(range[1].abs())) * k;
                    for end in [pa, pb] {
                        let e = |theta: f64| {
                            let [p, q] = sight.terms(theta);
                            (end[1] + law.taper * k * theta) * (p - w * sight.vz)
                                - k * (w - law.taper) * q
                        };
                        let c = end[1].abs() + reach;
                        let lipschitz = law.taper.abs() * k * (sight.rho + (w * sight.vz).abs())
                            + c * sight.rho
                            + k * (w - law.taper).abs() * sight.rho
                            + 1e-12;
                        isolate(&e, range, lipschitz, &mut events)?;
                    }
                    // Poles of the closed form, where P = w v_z.
                    let ratio = w * sight.vz / sight.rho;
                    if ratio.abs() <= 1. {
                        for sign in [-1., 1.] {
                            let theta = sight.alpha + sign * ratio.portable_acos();
                            let first = ((range[0] - theta) / TAU).floor() as i64;
                            let last = ((range[1] - theta) / TAU).ceil() as i64;
                            for j in first..=last {
                                let at = theta + j as f64 * TAU;
                                if at > range[0] && at < range[1] {
                                    events.push(at);
                                }
                            }
                        }
                    }
                }
                Span::Arc {
                    center,
                    radius,
                    sign,
                    range: arc,
                } => {
                    self.certify_windows(center, radius, sign, &sight, range)?;
                    let lipschitz = self.arc_lipschitz(center, radius, &sight, range);
                    for end in arc {
                        let beta = ((end - center[0]) / radius).clamp(-1., 1.).portable_asin();
                        let e = |theta: f64| {
                            arc_value(
                                radius,
                                sign,
                                arc_terms(law, center[1], sign, &sight, theta),
                                beta,
                            )
                        };
                        isolate(&e, range, lipschitz, &mut events)?;
                    }
                }
            }
            events.sort_by(f64::total_cmp);
            events.dedup_by(|x, y| (*x - *y).abs() <= 1e-13);
            for &window in windows {
                let u = |theta: f64| Contour::u(law, span, window, &sight, theta);
                for pair in events.windows(2) {
                    let mid = 0.5 * (pair[0] + pair[1]);
                    if pair[1] - pair[0] <= 1e-13 || !u(mid).is_some_and(|u| u > a && u < b) {
                        continue;
                    }
                    let level = |theta: f64| u(theta).map(|u| u + k * theta);
                    let mut stops = vec![pair[0], pair[1]];
                    for &(cap, theta) in &ends {
                        // Cap contour points lie on the branch at their angle.
                        for turn in [0., TAU, -TAU] {
                            let theta = theta + turn;
                            let shift = ((mid - theta) / TAU).round() * TAU;
                            let theta = theta + shift;
                            if theta > pair[0]
                                && theta < pair[1]
                                && level(theta).is_some_and(|z| (z - cap).abs() < 1e-7)
                            {
                                stops.push(theta);
                            }
                        }
                    }
                    stops.sort_by(f64::total_cmp);
                    stops.dedup_by(|x, y| (*x - *y).abs() <= 1e-13);
                    let width = pair[1] - pair[0];
                    for piece in stops.windows(2) {
                        let centre = 0.5 * (piece[0] + piece[1]);
                        let inside = level(centre).is_some_and(|z| z > 0. && z < law.length);
                        if !inside || piece[1] - piece[0] <= 1e-13 {
                            continue;
                        }
                        output.push(Candidate {
                            kind: Kind::Outline,
                            path: Arc::new(Range {
                                path: Arc::new(Contour {
                                    law: law.clone(),
                                    span: *span,
                                    window,
                                    sight,
                                    range: [pair[0], pair[1]],
                                }),
                                range: [(piece[0] - pair[0]) / width, (piece[1] - pair[0]) / width],
                            }),
                        });
                    }
                }
            }
        }
        Ok(())
    }

    /// Contour generators of a cone band: where `rho cos psi = tau v_z`.
    fn cone_generators(&self, span: &Span, sight: &Sight, output: &mut Vec<Candidate>) {
        let law = &self.law;
        let ratio = law.taper * sight.vz / sight.rho;
        if ratio.abs() >= 1. {
            return;
        }
        let [a, b] = span.range();
        let k = law.rate;
        for sign in [-1., 1.] {
            let theta = sight.alpha + sign * ratio.portable_acos();
            let (sin, cos) = theta.portable_sin_cos();
            let first = ((-b / k - theta) / TAU).floor() as i64 - 1;
            let last = (((law.length - a) / k - theta) / TAU).ceil() as i64 + 1;
            for j in first..=last {
                let s = k * (theta + j as f64 * TAU);
                let [lo, hi] = [(a + s).max(0.), (b + s).min(law.length)];
                if hi - lo <= 1e-12 {
                    continue;
                }
                let at = |z: f64| {
                    let r = span.eval(z - s)[0] + law.taper * s;
                    law.view([r * cos, r * sin, z])
                };
                output.push(Candidate {
                    kind: Kind::Outline,
                    path: Arc::new(Segment(at(lo), at(hi))),
                });
            }
        }
    }

    /// Bounds of the derivatives in the angle of `E`, `F` and the
    /// perturbation of an arc's contour condition, summed: a Lipschitz bound
    /// of the condition at a fixed circle angle.
    fn arc_lipschitz(&self, center: P2, radius: f64, sight: &Sight, range: [f64; 2]) -> f64 {
        let law = &self.law;
        let (k, tau, rho, vz) = (law.rate, law.taper.abs(), sight.rho, sight.vz.abs());
        let c = center[1].abs() + tau * k * range[0].abs().max(range[1].abs());
        let de = tau * k * rho + c * rho + tau * k * rho;
        let df = tau * k * vz + k * rho;
        de + df + radius * rho + 1e-12
    }

    /// Prove `N^2 - (R / WINDOW)^2 (P^2 + v_z^2) > 0` over the angles, so
    /// each window of an arc's contour condition holds exactly one root.
    fn certify_windows(
        &self,
        center: P2,
        radius: f64,
        sign: f64,
        sight: &Sight,
        range: [f64; 2],
    ) -> Result<(), GeometryError> {
        let law = &self.law;
        let (k, tau, rho, vz) = (law.rate, law.taper.abs(), sight.rho, sight.vz.abs());
        let c = center[1].abs() + tau * k * range[0].abs().max(range[1].abs());
        let e = c * rho + tau * k * rho;
        let f = c * vz + k * rho;
        let de = tau * k * rho + c * rho + tau * k * rho;
        let df = tau * k * vz + k * rho;
        let factor = (radius / WINDOW).powi(2);
        let lipschitz = 2. * e * de + 2. * f * df + 2. * factor * rho * rho + 1e-12;
        let margin = 1e-9 * (e * e + f * f).max(1e-300);
        let q = |theta: f64| {
            let [e, f, p, vz] = arc_terms(law, center[1], sign, sight, theta);
            e * e + f * f - factor * (p * p + vz * vz) - margin
        };
        let mut roots = Vec::new();
        isolate(&q, range, lipschitz, &mut roots)?;
        if !roots.is_empty() || q(range[0]) <= 0. {
            return Err(unsupported(
                "Drawing a pipe thread's rounded contour in a view this close to its lead angle is not supported",
            ));
        }
        Ok(())
    }
}

impl Pipe {
    /// The chords of the cap region at `level` cut by the plane, as local
    /// end points.
    fn cap_chord(&self, n: P3, c: f64, level: f64) -> Result<Vec<[P3; 2]>, GeometryError> {
        let law = &self.law;
        let nxy = n[0].portable_hypot(n[1]);
        if nxy <= 1e-12 {
            return Ok(Vec::new());
        }
        let d = [n[1] / nxy, -n[0] / nxy, 0.];
        let offset = (c - n[2] * level) / (nxy * nxy);
        let s = [n[0] * offset, n[1] * offset, level];
        let g = |t: f64| {
            let p = add(s, scaled(d, t));
            law.radius(p[1].portable_atan2(p[0]), level) - p[0].portable_hypot(p[1])
        };
        let foot = s[0].portable_hypot(s[1]);
        if foot >= law.maximum {
            return Ok(Vec::new());
        }
        let half = (law.maximum * law.maximum - foot * foot).sqrt();
        let parts: Vec<[f64; 2]> = if foot < law.minimum {
            let inner = (law.minimum * law.minimum - foot * foot).sqrt();
            vec![[-half, -inner], [inner, half]]
        } else {
            vec![[-half, half]]
        };
        let mut stops = vec![-half, half];
        for part in parts {
            stops.extend(part);
            stops.extend(self.crossings(s, d, part)?);
        }
        stops.sort_by(f64::total_cmp);
        stops.dedup_by(|x, y| (*x - *y).abs() <= 1e-13);
        let mut spans = Vec::new();
        for pair in stops.windows(2) {
            if pair[1] - pair[0] > 1e-13 && g(0.5 * (pair[0] + pair[1])) > 0. {
                spans.push([pair[0], pair[1]]);
            }
        }
        Ok(merge(spans, 0.)
            .into_iter()
            .map(|[a, b]| [add(s, scaled(d, a)), add(s, scaled(d, b))])
            .collect())
    }

    /// The section boundary by the local plane `n . p = c`.
    fn section_loop(&self, n: P3, c: f64) -> Result<Loop, GeometryError> {
        let law = &self.law;
        let k = law.rate;
        let nxy = n[0].portable_hypot(n[1]);
        let mut output: Loop = Vec::new();
        if nxy <= 1e-12 {
            // Across the axis: the cap curve at that level.
            let level = c / n[2];
            let probe = level - n[2].signum() * 4. * self.tolerance;
            if probe > 0. && probe < law.length {
                output.extend(self.cap(level));
            }
            return Ok(output);
        }
        if n[2].abs() <= 1e-12 && c.abs() <= 4. * self.tolerance {
            // Containing the axis: both meridian chains and the cap lines.
            let h = [-n[1] / nxy, n[0] / nxy, 0.];
            let theta = h[1].portable_atan2(h[0]);
            for (side, angle) in [(1., theta), (-1., theta + PI)] {
                for m in law.copies() {
                    let shift = k * angle + m as f64 * law.pitch;
                    for span in &law.spans {
                        let local = law.translated(span, shift);
                        let [a, b] = local.range();
                        let [lo, hi] = [a.max(0.), b.min(law.length)];
                        if hi - lo <= 1e-12 {
                            continue;
                        }
                        let point = |z: f64| {
                            let r = local.eval(z)[0];
                            law.view([side * r * h[0], side * r * h[1], z])
                        };
                        let path: Arc<dyn Path> = match local {
                            Span::Line { .. } => Arc::new(Segment(point(lo), point(hi))),
                            Span::Arc { center, sign, .. } => {
                                // A circular arc in the plane: an exact conic.
                                let lift = Lift {
                                    origin: law.view([0.; 3]),
                                    x: self.view_vector([0., 0., 1.]),
                                    y: self.view_vector([side * h[0], side * h[1], 0.]),
                                };
                                let end = |z: f64| [z, local.eval(z)[0]];
                                Arc::new(super::paths::ProfilePath {
                                    edge: Edge::Planar(normal_sweep::Edge::Arc {
                                        center,
                                        a: end(lo),
                                        b: end(hi),
                                        clockwise: sign > 0.,
                                    }),
                                    range: [0., 1.],
                                    lift,
                                    straight: false,
                                })
                            }
                        };
                        output.push(path);
                    }
                }
            }
            for level in [0., law.length] {
                let at = |s: f64| law.view([s * h[0], s * h[1], level]);
                let [right, left] = [law.radius(theta, level), -law.radius(theta + PI, level)];
                output.push(Arc::new(Segment(at(left), at(right))));
            }
            return Ok(output);
        }
        // Chords of the caps; their ends are where section branches meet a cap.
        let mut chords = Vec::new();
        for level in [0., law.length] {
            chords.extend(self.cap_chord(n, c, level)?);
        }
        let chord_ends: Vec<P3> = chords.iter().flatten().copied().collect();
        let tau = law.taper.abs();
        // Sweep the meridian half-plane once around the axis, in two charts
        // of half a turn each centred where the meridian half-plane lies in
        // the plane's direction. There the plane meets the meridian lines at
        // the most grazing angle (near the axis a profile point moves by the
        // angle's error times the radius squared over the offset), so the
        // angle is measured from the chart centre at full precision.
        let base = n[0].portable_atan2(-n[1]);
        let [low, high] = [*law.copies().start() - 1, *law.copies().end() + 1];
        for base in [base, base + PI] {
            let (sin, cos) = base.portable_sin_cos();
            let n = [n[0] * cos + n[1] * sin, n[1] * cos - n[0] * sin, n[2]];
            let shift = k * base;
            let chart = [-PI / 2., PI / 2.];
            for m in low..=high {
                for span in &law.spans {
                    let span = law.translated(span, m as f64 * law.pitch + shift);
                    let mut events = chart.to_vec();
                    let mut touches = Vec::new();
                    let [a, b] = span.range();
                    let drift = nxy * tau * k * (TAU + 1.) + n[2].abs() * k;
                    for e in [span.point(a), span.point(b)] {
                        let f = |theta: f64| {
                            let [alpha, beta, gamma] = Sweep::line(law, (n, c), theta);
                            alpha * e[1] + beta * e[0] - gamma
                        };
                        let lip = nxy * e[1].abs() + drift + 1e-12;
                        isolate(&f, chart, lip, &mut events)?;
                    }
                    if let Span::Arc { center, radius, .. } = span {
                        // The meridian line touches the circle where its signed
                        // distance from the center, `u / |(alpha, beta)|`, is
                        // `+-R`: the roots of `u -+ R |(alpha, beta)|`, isolated
                        // one sign at a time so that a pair of close touches
                        // never needs resolving as one near-double root.
                        let lu = nxy * center[1].abs() + drift;
                        for side in [-1., 1.] {
                            let f = |theta: f64| {
                                let [alpha, beta, gamma] = Sweep::line(law, (n, c), theta);
                                alpha * center[1] + beta * center[0] - gamma
                                    + side * radius * alpha.portable_hypot(beta)
                            };
                            let lip = lu + radius * nxy + 1e-12;
                            isolate(&f, chart, lip, &mut touches)?;
                        }
                    }
                    events.extend(touches.iter().copied());
                    events.sort_by(f64::total_cmp);
                    events.dedup_by(|x, y| (*x - *y).abs() <= 1e-13);
                    let touch = |theta: f64| touches.iter().any(|t| (t - theta).abs() <= 1e-13);
                    for pair in events.windows(2) {
                        if pair[1] - pair[0] <= 1e-13 {
                            continue;
                        }
                        let mid = 0.5 * (pair[0] + pair[1]);
                        let [alpha, beta, gamma] = Sweep::line(law, (n, c), mid);
                        for (index, _) in span.meet(alpha, beta, gamma) {
                            let branch = Sweep {
                                law: law.clone(),
                                base,
                                branch: Branch::new(
                                    span,
                                    index,
                                    Sweep::meridian(law, (n, c)),
                                    [pair[0], pair[1]],
                                    [touch(pair[0]), touch(pair[1])],
                                ),
                            };
                            // The branch leaves the material only through a
                            // cap, at an end of that cap's exact chord.
                            let mut stops = vec![0., 1.];
                            for point in &chord_ends {
                                let theta = point[1].portable_atan2(point[0]) - base;
                                let theta = theta - TAU * (theta / TAU).round();
                                let profile = branch.profile(*point);
                                let Some(t) = branch.branch.parameter(theta, &profile) else {
                                    continue;
                                };
                                if t <= 0. || t >= 1. {
                                    continue;
                                }
                                // Where the cap is crossed, exactly: near the
                                // axis the chord end's angle alone places the
                                // branch point poorly along a span that the
                                // meridian lines meet at a grazing angle.
                                let z = |t: f64| law.local(branch.point(t))[2];
                                let t = level_stop(&z, t, point[2]);
                                let q = law.local(branch.point(t));
                                let d = add(q, scaled(*point, -1.));
                                if dot3(d, d).sqrt() < 1e-6 {
                                    stops.push(t);
                                }
                            }
                            let sweep: Arc<dyn Proven> = Arc::new(branch);
                            stops.sort_by(f64::total_cmp);
                            stops.dedup_by(|x, y| (*x - *y).abs() <= 1e-13);
                            for piece in stops.windows(2) {
                                let z = law.local(sweep.point(0.5 * (piece[0] + piece[1])))[2];
                                if z > 0. && z < law.length && piece[1] - piece[0] > 1e-13 {
                                    output.push(Arc::new(Range {
                                        path: sweep.clone(),
                                        range: [piece[0], piece[1]],
                                    }));
                                }
                            }
                        }
                    }
                }
            }
        }
        for [a, b] in chords {
            output.push(Arc::new(Segment(law.view(a), law.view(b))));
        }
        Ok(output)
    }
}

impl Member for Pipe {
    fn bounds(&self) -> [P3; 2] {
        let law = &self.law;
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for z in [0., law.length] {
            let center = law.view([0., 0., z]);
            for i in 0..3 {
                let axis = law.axes[2][i];
                let extent = law.maximum * (1. - axis * axis).max(0.).sqrt();
                lo[i] = lo[i].min(center[i] - extent);
                hi[i] = hi[i].max(center[i] + extent);
            }
        }
        [lo, hi]
    }

    fn candidates(&self, _: f64) -> Result<Vec<Candidate>, GeometryError> {
        let law = &self.law;
        let mut output = Vec::new();
        for (i, kind) in law.rails.iter().enumerate() {
            if let Some(kind) = *kind {
                self.rail(law.spans[i].range()[0], kind, &mut output);
            }
        }
        for level in [0., law.length] {
            output.extend(self.cap(level).into_iter().map(|path| Candidate {
                kind: Kind::Edge,
                path,
            }));
        }
        self.outlines(&mut output)?;
        Ok(output)
    }

    fn sight(&self, [x, y]: P2) -> Result<Vec<[f64; 2]>, GeometryError> {
        self.own_line([x, y, 0.], [0., 0., 1.])
    }

    fn section(&self, depth: f64, _: f64) -> Result<Vec<Loop>, GeometryError> {
        let n = self.law.vector([0., 0., 1.]);
        let c = depth - self.law.origin[2];
        let output = self.section_loop(n, c)?;
        Ok(if output.is_empty() {
            vec![]
        } else {
            vec![output]
        })
    }

    fn line(&self, start: P3, direction: P3) -> Option<Result<Vec<[f64; 2]>, GeometryError>> {
        Some(self.own_line(start, direction))
    }
}
