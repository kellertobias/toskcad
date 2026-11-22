//! Straight single- or multi-start screw threads: material inside the
//! helicoid `r <= R(z - k * theta)` between two flat caps, with an optional
//! coaxial cylindrical neck beyond the far cap (a terminal thread relief),
//! an entry cone, or a screw's shank and head. Threaded bores (a rounded
//! thread bore, an internal thread treatment) are the stock cylinder less
//! the helicoid's void, which a treatment's entry cone widens to
//! `r < max(R, C(z))`.
//!
//! The axial meridian profile `R` is a periodic chain of exact line segments
//! and circular arcs read from the thread source. Every answer below is
//! evaluated on that profile and the helical law; nothing is sampled from a
//! display or manufacturing mesh.
//!
//! Exactness: rails, cap curves, outline curves and section curves are exact
//! evaluable paths (closed-form points, or closed-form roots of the meridian
//! line and the profile). Every event is isolated with a proven bound: line
//! material, section sweep events, outline existence (roots of the
//! discriminant), contour ends on the caps, and branch ends on the caps (the
//! exact cap chord ends). Helices and cap curves flatten with proven
//! curvature bounds. Outline, seam, cone section and section branch paths
//! flatten with proven per-chord bounds (see [`proven`]): outlines and
//! branches are components of the zero set of an exact relation in a
//! two-coordinate chart, so near a branch end, where the native
//! parametrisation has a square-root singularity, the bound is taken in the
//! other chart coordinate, in which the curve is smooth.
use multicad_precision::PortableMath as _;
mod fastener;
mod proven;

use super::paths::{Lift, Segment, add, dot3, scaled};
use super::prism::{Band, Prism};
use super::profile::merge;
use super::*;
use multicad_rust_kernel::{helical_thread, normal_sweep, rounded_thread, spur_gear::Edge};
pub(super) use proven::{Chart, Iv, Num, Proven, bracket, certified, chord_bound, frame, settle};
#[cfg(test)]
pub(super) use proven::{Tally, tally};
use std::f64::consts::{PI, TAU};

/// One piece of the axial meridian profile in `(phase, radius)` coordinates.
#[derive(Clone, Copy, Debug)]
pub(super) enum Span {
    Line {
        a: P2,
        b: P2,
    },
    /// `r = center[1] + sign * sqrt(radius^2 - (phi - center[0])^2)` over the
    /// phase range.
    Arc {
        center: P2,
        radius: f64,
        sign: f64,
        range: [f64; 2],
    },
}

impl Span {
    pub(super) fn range(&self) -> [f64; 2] {
        match *self {
            Self::Line { a, b } => [a[0], b[0]],
            Self::Arc { range, .. } => range,
        }
    }

    fn shifted(self, offset: f64) -> Self {
        match self {
            Self::Line { a, b } => Self::Line {
                a: [a[0] + offset, a[1]],
                b: [b[0] + offset, b[1]],
            },
            Self::Arc {
                center,
                radius,
                sign,
                range,
            } => Self::Arc {
                center: [center[0] + offset, center[1]],
                radius,
                sign,
                range: range.map(|x| x + offset),
            },
        }
    }

    /// Radius and slope at a phase inside the span.
    pub(super) fn eval(&self, phi: f64) -> [f64; 2] {
        match *self {
            Self::Line { a, b } => {
                let slope = (b[1] - a[1]) / (b[0] - a[0]);
                [a[1] + slope * (phi - a[0]), slope]
            }
            Self::Arc {
                center,
                radius,
                sign,
                ..
            } => {
                let x = phi - center[0];
                let h = ((radius - x) * (radius + x)).max(0.).sqrt();
                [center[1] + sign * h, -sign * x / h.max(1e-300)]
            }
        }
    }

    pub(super) fn point(&self, phi: f64) -> P2 {
        [phi, self.eval(phi)[0]]
    }

    /// Radius and slope at a phase, in chart arithmetic, as `eval`.
    pub(super) fn radius_at<N: Num>(&self, phi: N) -> [N; 2] {
        match *self {
            Self::Line { a, b } => {
                let slope = (b[1] - a[1]) / (b[0] - a[0]);
                [N::k(a[1]) + N::k(slope) * (phi - N::k(a[0])), N::k(slope)]
            }
            Self::Arc {
                center,
                radius,
                sign,
                ..
            } => {
                let x = phi - N::k(center[0]);
                let h = ((N::k(radius) - x) * (N::k(radius) + x)).sqrt();
                [N::k(center[1]) + N::k(sign) * h, -(N::k(sign) * x / h)]
            }
        }
    }

    /// Phases inside the span where the radius is stationary.
    pub(super) fn stationary(&self) -> Option<f64> {
        match *self {
            Self::Line { a, b } => (a[1] == b[1]).then_some(0.5 * (a[0] + b[0])),
            Self::Arc { center, range, .. } => {
                (range[0] <= center[0] && center[0] <= range[1]).then_some(center[0])
            }
        }
    }

    /// Phases inside the span where the radius equals `radius`.
    pub(super) fn at_radius(&self, radius: f64) -> Vec<f64> {
        let [lo, hi] = self.range();
        let within = |phi: &f64| *phi > lo && *phi < hi;
        match *self {
            Self::Line { a, b } => {
                if a[1] == b[1] || (a[1] - radius) * (b[1] - radius) >= 0. {
                    return Vec::new();
                }
                vec![a[0] + (b[0] - a[0]) * (radius - a[1]) / (b[1] - a[1])]
                    .into_iter()
                    .filter(within)
                    .collect()
            }
            Self::Arc {
                center,
                radius: rho,
                sign,
                ..
            } => {
                let dr = radius - center[1];
                if sign * dr < 0. || dr.abs() >= rho {
                    return Vec::new();
                }
                let w = ((rho - dr.abs()) * (rho + dr.abs())).sqrt();
                [center[0] - w, center[0] + w]
                    .into_iter()
                    .filter(within)
                    .collect()
            }
        }
    }

    fn flat(&self) -> Option<f64> {
        match *self {
            Self::Line { a, b } if a[1] == b[1] => Some(a[1]),
            _ => None,
        }
    }

    /// Up to two intersections of the meridian line `alpha r + beta phi = gamma`
    /// with the span, as `(phi, r)`, tagged by root index.
    pub(super) fn meet(&self, alpha: f64, beta: f64, gamma: f64) -> Vec<(usize, P2)> {
        match *self {
            Self::Line { a, b } => {
                let [fa, fb] = [a, b].map(|p| alpha * p[1] + beta * p[0] - gamma);
                if fa == fb || (fa > 0.) == (fb > 0.) && fa != 0. && fb != 0. {
                    return Vec::new();
                }
                let t = (fa / (fa - fb)).clamp(0., 1.);
                vec![(0, [a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])])]
            }
            Self::Arc { .. } => [0, 1]
                .into_iter()
                .filter_map(|index| {
                    let p = self.arc_root(alpha, beta, gamma, index, false)?;
                    Some((index, p))
                })
                .collect(),
        }
    }

    /// The `index`-th intersection of the meridian line with the full circle
    /// of an arc; with `clamp`, the nearest point when the line misses it.
    pub(super) fn arc_root(
        &self,
        alpha: f64,
        beta: f64,
        gamma: f64,
        index: usize,
        clamp: bool,
    ) -> Option<P2> {
        let Self::Arc {
            center,
            radius,
            sign,
            range,
        } = *self
        else {
            return None;
        };
        let norm = alpha.portable_hypot(beta);
        if norm <= 1e-300 {
            return None;
        }
        // Normal (beta, alpha) in (phi, r); direction along the line.
        let n = [beta / norm, alpha / norm];
        let u = [-n[1], n[0]];
        let delta = (beta * center[0] + alpha * center[1] - gamma) / norm;
        let reach = (radius - delta.abs()) * (radius + delta.abs());
        if reach < 0. && !clamp {
            return None;
        }
        // Index 2 is the foot point, where both roots meet at a touch.
        let h = match index {
            0 => -reach.max(0.).sqrt(),
            1 => reach.max(0.).sqrt(),
            _ => 0.,
        };
        let p = [
            center[0] - delta * n[0] + h * u[0],
            center[1] - delta * n[1] + h * u[1],
        ];
        if clamp {
            return Some(p);
        }
        let on_side = sign * (p[1] - center[1]) >= -1e-12 * radius;
        let within = p[0] >= range[0] - 1e-12 * radius && p[0] <= range[1] + 1e-12 * radius;
        (on_side && within).then_some(p)
    }
}

/// The periodic profile over one pitch, with the kind of each break.
#[derive(Clone, Debug)]
struct Profile {
    pitch: f64,
    spans: Vec<Span>,
    /// Kind of the rail at the start of each span; `None` where the profile
    /// continues straight through (no edge at all).
    rails: Vec<Option<Kind>>,
    minimum: f64,
    maximum: f64,
    slope: f64,
    /// Bound of the profile's second derivative inside its spans.
    curvature: f64,
}

impl Profile {
    fn new(pitch: f64, spans: Vec<Span>) -> Self {
        let count = spans.len();
        let rails = (0..count)
            .map(|i| {
                let before = spans[(i + count - 1) % count];
                let after = spans[i];
                let phi = after.range()[0];
                // The previous span ends where this one starts, a pitch
                // later for the first span.
                let back = before.range()[1];
                let [r0, s0] = before.eval(back);
                let [r1, s1] = after.eval(phi);
                let scale = 1. + s0.abs().max(s1.abs());
                if (r0 - r1).abs() > 1e-9 * r0.abs().max(1.) {
                    return Some(Kind::Edge);
                }
                match (before.flat(), after.flat()) {
                    (Some(_), Some(_)) => None,
                    _ if (s0 - s1).abs() <= 1e-7 * scale => Some(Kind::Smooth),
                    _ => Some(Kind::Edge),
                }
            })
            .collect();
        let mut minimum = f64::INFINITY;
        let mut maximum = f64::NEG_INFINITY;
        let mut slope: f64 = 0.;
        let mut curvature: f64 = 0.;
        for span in &spans {
            if let Span::Arc {
                center,
                radius,
                range,
                ..
            } = *span
            {
                let x = (range[0] - center[0])
                    .abs()
                    .max((range[1] - center[0]).abs());
                let h = ((radius - x) * (radius + x)).max(0.).sqrt();
                curvature = curvature.max(radius * radius / h.max(1e-300).powi(3));
            }
            let [a, b] = span.range();
            for phi in [a, b] {
                let [r, s] = span.eval(phi);
                minimum = minimum.min(r);
                maximum = maximum.max(r);
                if s.is_finite() {
                    slope = slope.max(s.abs());
                }
            }
            if let Some(phi) = span.stationary() {
                let r = span.eval(phi)[0];
                minimum = minimum.min(r);
                maximum = maximum.max(r);
            }
        }
        Self {
            pitch,
            spans,
            rails,
            minimum,
            maximum,
            // Line slopes are constant and arc slopes are monotone in the
            // phase, so the extreme lies at a span end.
            slope: slope * (1. + 1e-12) + 1e-12,
            curvature: curvature * (1. + 1e-12) + 1e-12,
        }
    }

    fn locate(&self, phi: f64) -> (usize, f64) {
        let local = phi.rem_euclid(self.pitch);
        let index = self
            .spans
            .iter()
            .position(|span| local <= span.range()[1])
            .unwrap_or(self.spans.len() - 1);
        (index, phi - local)
    }

    fn radius(&self, phi: f64) -> f64 {
        let (index, base) = self.locate(phi);
        self.spans[index].eval(phi - base)[0]
    }

    fn from_helical(thread: &helical_thread::Thread) -> Option<Self> {
        let breaks = thread.profile_breaks();
        let spec = thread.spec();
        let pitch = spec.pitch_mm;
        let r = |phi: f64| thread.radius_at_phase(phi);
        let line = |a: f64, b: f64| Span::Line {
            a: [a, r(a)],
            b: [b, r(b)],
        };
        let root = if spec.root_radius_mm.is_some() {
            // The tangent root circle through three exact profile points.
            let [a, b] = [breaks[2], breaks[3]];
            let m = 0.5 * (a + b);
            let center = circle([a, r(a)], [m, r(m)], [b, r(b)])?;
            Span::Arc {
                center,
                radius: (a - center[0]).portable_hypot(r(a) - center[1]),
                sign: -1.,
                range: [a, b],
            }
        } else {
            line(breaks[2], breaks[3])
        };
        let spans = vec![
            Span::Line {
                a: [breaks[0], spec.major_radius_mm],
                b: [breaks[1], spec.major_radius_mm],
            },
            line(breaks[1], breaks[2]),
            root,
            line(breaks[3], breaks[4]),
            Span::Line {
                a: [breaks[4], spec.major_radius_mm],
                b: [breaks[5], spec.major_radius_mm],
            },
        ];
        (breaks[5] == pitch).then(|| Self::new(pitch, spans))
    }

    fn from_rounded(thread: &rounded_thread::Thread) -> Option<Self> {
        let breaks = thread.profile_breaks();
        let pitch = thread.spec().pitch_mm;
        let spans = (0..5)
            .map(|band| {
                let [a, b] = [breaks[band], breaks[band + 1]];
                match thread.profile_circle(band) {
                    Some(c) => Some(Span::Arc {
                        center: [c.center_phase_mm, c.center_radius_mm],
                        radius: c.radius_mm,
                        sign: c.radial_sign,
                        range: [a, b],
                    }),
                    None => Some(Span::Line {
                        a: [a, thread.radius_at_phase(a)],
                        b: [b, thread.radius_at_phase(b)],
                    }),
                }
            })
            .collect::<Option<Vec<_>>>()?;
        Some(Self::new(pitch, spans))
    }
}

/// Phases (z in the meridian chart of a placed span) where a span meets the
/// entry cone's generator `r = entry + slope z`, inside the cone's length.
fn cone_crossings(span: &Span, cone: &Cone) -> Vec<f64> {
    let [lo, hi] = span.range();
    let mut output = Vec::new();
    match *span {
        Span::Line { a, b } => {
            let fa = a[1] - cone.radius(a[0]);
            let fb = b[1] - cone.radius(b[0]);
            if fa != fb && (fa > 0.) != (fb > 0.) {
                output.push(a[0] + (b[0] - a[0]) * fa / (fa - fb));
            }
        }
        Span::Arc {
            center,
            radius,
            sign,
            ..
        } => {
            // (z - zc)^2 + (e + s z - rc)^2 = radius^2.
            let k = cone.slope;
            let w = cone.entry - center[1];
            let qa = 1. + k * k;
            let qb = -2. * center[0] + 2. * k * w;
            let qc = center[0] * center[0] + w * w - radius * radius;
            let disc = qb * qb - 4. * qa * qc;
            if disc >= 0. {
                let root = disc.sqrt();
                for z in [(-qb - root) / (2. * qa), (-qb + root) / (2. * qa)] {
                    if sign * (cone.radius(z) - center[1]) >= 0. {
                        output.push(z);
                    }
                }
            }
        }
    }
    output
        .into_iter()
        .filter(|z| *z >= lo && *z <= hi && *z >= 0. && *z <= cone.length)
        .collect()
}

/// A counterclockwise hexagon ring of the given width across flats.
fn hexagon(across_flats_mm: f64) -> Vec<Edge> {
    let corners = multicad_rust_kernel::headed_screw::hex_vertices(across_flats_mm);
    (0..6)
        .map(|i| {
            Edge::Planar(normal_sweep::Edge::Line {
                a: corners[i],
                b: corners[(i + 1) % 6],
            })
        })
        .collect()
}

/// The center of the circle through three points.
fn circle(a: P2, b: P2, c: P2) -> Option<P2> {
    let d = 2. * (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]));
    if d.abs() <= 1e-300 {
        return None;
    }
    let s = |p: P2| p[0] * p[0] + p[1] * p[1];
    Some([
        (s(a) * (b[1] - c[1]) + s(b) * (c[1] - a[1]) + s(c) * (a[1] - b[1])) / d,
        (s(a) * (c[0] - b[0]) + s(b) * (a[0] - c[0]) + s(c) * (b[0] - a[0])) / d,
    ])
}

/// Distance from the surface, in millimetres, below which a line is not
/// resolved further. Harmless: an unresolved stretch is one where the line
/// stays within this distance of the surface on both sides, so any material
/// it misses is thinner than 1e-10 mm across the line and at most as long as
/// the stretch, where the Lipschitz bound times the length is below the same
/// distance. The engine hides only behind more than 1e-6 of the scene scale,
/// and the linear tolerance of every source is at least 1e-7 mm, so neither
/// visibility nor any published boundary depends on such material.
pub(super) const RESOLVED: f64 = 1e-10;

/// Roots of `f` on `[a, b]` with `|f'| <= lipschitz`: a stretch is excluded
/// when its end values prove it free of zeros; sign changes are bisected.
pub(super) fn isolate(
    f: &dyn Fn(f64) -> f64,
    [a, b]: [f64; 2],
    lipschitz: f64,
    output: &mut Vec<f64>,
) -> Result<(), GeometryError> {
    let width = 1e-13 * a.abs().max(b.abs()).max(1.);
    // A function as steep as its bound reaches equality in the exclusion
    // test at a simple root; rounding must never turn that into exclusion.
    let lipschitz = lipschitz * (1. + 1e-6) + 1e-12;
    let mut pending = vec![(a, f(a), b, f(b))];
    let mut work = 0usize;
    while let Some((a, fa, b, fb)) = pending.pop() {
        work += 1;
        if work > 2_000_000 {
            return Err(err(
                GeometryErrorKind::Limit,
                "Thread drawing root isolation exceeded its work budget",
            ));
        }
        // A sign change proves a root; rounding of the two values must never
        // exclude it.
        let change = (fa <= 0.) != (fb <= 0.);
        if !change && fa.abs() + fb.abs() > lipschitz * (b - a) {
            continue;
        }
        // Where the function cannot leave a band of `RESOLVED` around zero,
        // the line stays that close to the surface: keep only sign changes.
        if b - a <= width || lipschitz * (b - a) <= RESOLVED {
            if change {
                output.push(bisect(f, [a, b], fa));
            }
            continue;
        }
        let m = 0.5 * (a + b);
        let fm = f(m);
        pending.push((a, fa, m, fm));
        pending.push((m, fm, b, fb));
    }
    output.sort_by(f64::total_cmp);
    output.dedup_by(|x, y| (*x - *y).abs() <= 4. * width);
    Ok(())
}

/// A sign change of `f` on `[a, b]`, with `fa = f(a)`, bisected to the
/// last representable bracket. Near an axis a plane's meridian lines meet
/// the profile at a grazing angle, so the profile point moves by the root's
/// error times a large factor: the root is resolved fully, not just to the
/// isolation width.
fn bisect(f: &dyn Fn(f64) -> f64, [mut a, mut b]: [f64; 2], mut fa: f64) -> f64 {
    for _ in 0..80 {
        let m = 0.5 * (a + b);
        if m <= a || m >= b {
            break;
        }
        let fm = f(m);
        if (fm <= 0.) == (fa <= 0.) {
            a = m;
            fa = fm;
        } else {
            b = m;
        }
    }
    0.5 * (a + b)
}

/// The parameter near `t` at which the path's height `z` crosses `level`,
/// bisected to the last representable bracket; `t` itself when no crossing
/// is bracketed within a hundredth of the parameter range.
pub(super) fn level_stop(z: &dyn Fn(f64) -> f64, t: f64, level: f64) -> f64 {
    let f = |t: f64| z(t) - level;
    let ft = f(t);
    if ft == 0. {
        return t;
    }
    let mut step = 1e-12;
    while step < 1e-2 {
        for other in [t - step, t + step] {
            if (0. ..=1.).contains(&other) && (f(other) <= 0.) != (ft <= 0.) {
                let [a, b] = if other < t { [other, t] } else { [t, other] };
                return bisect(&f, [a, b], f(a));
            }
        }
        step *= 4.;
    }
    t
}

/// Shared placement and the helical law of one member.
#[derive(Clone)]
struct Law {
    origin: P3,
    axes: [P3; 3],
    profile: Profile,
    rate: f64,
    starts: u32,
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

    fn surface(&self, r: f64, theta: f64, z: f64) -> P3 {
        let (sin, cos) = theta.portable_sin_cos();
        self.view([r * cos, r * sin, z])
    }
}

/// A helix of constant radius: `theta` from `range`, `z = base + rate * theta`.
struct Helix {
    law: Arc<Law>,
    radius: f64,
    base: f64,
    range: [f64; 2],
}

impl Path for Helix {
    fn point(&self, t: f64) -> P3 {
        let theta = self.range[0] + (self.range[1] - self.range[0]) * t;
        self.law
            .surface(self.radius, theta, self.base + self.law.rate * theta)
    }

    fn flatten(&self, tolerance_mm: f64) -> Vec<f64> {
        // The helix curvature is at most 1 / radius in every coordinate set.
        let step = (8. * tolerance_mm / self.radius.max(1e-9)).sqrt().min(0.5);
        let count = ((self.range[1] - self.range[0]).abs() / step).ceil() as usize;
        let count = count.clamp(1, 1 << 20);
        (0..=count).map(|i| i as f64 / count as f64).collect()
    }

    fn paper(&self, [from, to]: [f64; 2]) -> Option<Curve> {
        // Seen along the axis, a helix projects onto a circle.
        if self.law.axes[2][0].portable_hypot(self.law.axes[2][1]) > 1e-12 {
            return None;
        }
        let c = self.law.view([0.; 3]);
        let [x, y] = [self.law.axes[0], self.law.axes[1]];
        let span = self.range[1] - self.range[0];
        paths::conic(
            [c[0], c[1]],
            [self.radius * x[0], self.radius * x[1]],
            [self.radius * y[0], self.radius * y[1]],
            self.range[0] + from * span,
            (to - from) * span,
        )
    }
}

/// The flat cap curve at level `level` over `theta` in `range`, inside one
/// shifted span.
struct Cap {
    law: Arc<Law>,
    span: Span,
    level: f64,
    range: [f64; 2],
}

impl Cap {
    fn theta(&self, t: f64) -> f64 {
        self.range[0] + (self.range[1] - self.range[0]) * t
    }
}

impl Path for Cap {
    fn point(&self, t: f64) -> P3 {
        let theta = self.theta(t);
        let r = self.span.eval(self.level - self.law.rate * theta)[0];
        self.law.surface(r, theta, self.level)
    }

    fn flatten(&self, tolerance_mm: f64) -> Vec<f64> {
        // |d2P/dtheta2| <= k^2 |R''| + 2 |k R'| + R inside one span, so a
        // chord over an angle h deviates by at most that bound times h^2 / 8.
        let profile = &self.law.profile;
        let k = self.law.rate.abs();
        let second = k * k * profile.curvature + 2. * k * profile.slope + profile.maximum;
        let step = (8. * tolerance_mm / second).sqrt().min(0.5);
        let count = ((self.range[1] - self.range[0]).abs() / step).ceil();
        let count = (count as usize).clamp(1, 1 << 20);
        (0..=count).map(|i| i as f64 / count as f64).collect()
    }

    fn paper(&self, [from, to]: [f64; 2]) -> Option<Curve> {
        let r = self.span.flat()?;
        let c = self.law.view([0., 0., self.level]);
        let [x, y] = [self.law.axes[0], self.law.axes[1]];
        let [a, b] = [self.theta(from), self.theta(to)];
        paths::conic(
            [c[0], c[1]],
            [r * x[0], r * x[1]],
            [r * y[0], r * y[1]],
            a,
            b - a,
        )
    }
}

/// One branch of the view outline on a span, over a phase range, for the
/// axial copy `copy`: the surface normal is perpendicular to the sight.
struct Outline {
    law: Arc<Law>,
    span: Span,
    sign: f64,
    copy: f64,
    range: [f64; 2],
    view: P3,
}

impl Outline {
    /// The outline angle at a phase, if the branch exists there.
    fn angle(law: &Law, span: &Span, view: P3, sign: f64, phi: f64) -> Option<f64> {
        Self::angle_within(law, span, view, sign, phi, 0.)
    }

    /// As `angle`, admitting the branch up to a relative `slack` beyond its
    /// end, where the acos term is evaluated at its limit.
    fn angle_within(
        law: &Law,
        span: &Span,
        view: P3,
        sign: f64,
        phi: f64,
        slack: f64,
    ) -> Option<f64> {
        let [r, slope] = span.eval(phi);
        let rho = view[0].portable_hypot(view[1]);
        let alpha = view[1].portable_atan2(view[0]);
        let a = -r * rho;
        let b = law.rate * slope * rho;
        let c = -r * slope * view[2];
        let norm = a.portable_hypot(b);
        if norm <= 1e-300 || c.abs() > norm * (1. + slack) {
            return None;
        }
        // atan2(b, a) is continuous as pi + atan2(-b, -a) since a < 0.
        Some(
            alpha + PI + (-b).portable_atan2(-a) + sign * (c / norm).clamp(-1., 1.).portable_acos(),
        )
    }

    fn at(&self, phi: f64) -> P3 {
        // Branch ends are isolated roots of the discriminant, where rounding
        // can put the acos argument just beyond one: take its limit there.
        let theta = Self::angle_within(&self.law, &self.span, self.view, self.sign, phi, 1e-9)
            .unwrap_or(f64::NAN);
        let r = self.span.eval(phi)[0];
        self.law
            .surface(r, theta, phi + self.law.rate * theta + self.copy)
    }

    fn phi(&self, t: f64) -> f64 {
        self.range[0] + (self.range[1] - self.range[0]) * t
    }

    /// `A`, `B` and `C` of the contour condition `A cos psi + B sin psi = C`.
    fn terms<N: Num>(&self, phi: N) -> [N; 3] {
        let [r, s] = self.span.radius_at(phi);
        let rho = self.view[0].portable_hypot(self.view[1]);
        [
            -(r * N::k(rho)),
            N::k(self.law.rate * rho) * s,
            -(r * s * N::k(self.view[2])),
        ]
    }

    /// The exact chart point `(phi, theta)` at a parameter. A branch end
    /// inside the span is an isolated root of the discriminant: its phase
    /// is enclosed by a proven sign change, and the outline angle there by
    /// the closed form over that enclosure.
    fn end(&self, t: f64) -> Option<[Iv; 2]> {
        let phi = self.phi(t);
        let [a, b] = self.span.range();
        let root = (t == 0. && self.range[0] != a) || (t == 1. && self.range[1] != b);
        let phi = if root {
            bracket(phi, |x| {
                let [p, q, r] = self.terms(x);
                p.sqr() + q.sqr() - r.sqr()
            })?
        } else {
            Iv::point(phi)
        };
        if !root {
            // A proven sign change next to the computed point, when it is
            // not at a fold of the branch.
            let theta =
                Self::angle_within(&self.law, &self.span, self.view, self.sign, phi.lo, 1e-9);
            if let Some(point) = theta.and_then(|theta| settle(self, [phi.lo, theta])) {
                return Some(point);
            }
        }
        let [p, q, r] = self.terms(phi);
        let norm = (p.sqr() + q.sqr()).sqrt();
        let alpha = self.view[1].portable_atan2(self.view[0]);
        // As `angle_within`: atan2(-B, -A) = atan(B / A) since A < 0.
        let theta =
            Iv::point(alpha) + Iv::pi() + (q / p).atan() + Iv::point(self.sign) * (r / norm).acos();
        Some([phi, theta])
    }
}

/// The outline is a component of `A cos psi + B sin psi - C = 0` in the
/// chart `(phi, theta)`, `psi = theta - alpha`.
impl Chart for Outline {
    fn implicit(&self) -> bool {
        true
    }

    fn relation<N: Num>(&self, phi: N, theta: N) -> N {
        let [a, b, c] = self.terms(phi);
        let psi = theta - N::k(self.view[1].portable_atan2(self.view[0]));
        a * psi.cos() + b * psi.sin() - c
    }

    fn place<N: Num>(&self, phi: N, theta: N) -> [N; 3] {
        let [r, _] = self.span.radius_at(phi);
        let law = &self.law;
        let z = phi + N::k(law.rate) * theta + N::k(self.copy);
        frame(law.origin, law.axes, [r * theta.cos(), r * theta.sin(), z])
    }
}

impl Proven for Outline {
    fn bound(&self, [a, b]: [f64; 2]) -> f64 {
        let (Some(start), Some(end)) = (self.end(a), self.end(b)) else {
            return f64::INFINITY;
        };
        chord_bound(self, 0, [start, end], [self.point(a), self.point(b)])
    }
}

impl Path for Outline {
    fn point(&self, t: f64) -> P3 {
        self.at(self.phi(t))
    }

    fn flatten(&self, tolerance_mm: f64) -> Vec<f64> {
        certified(self, tolerance_mm)
    }

    fn paper(&self, _: [f64; 2]) -> Option<Curve> {
        None
    }
}

/// A range of another path.
pub(super) struct Range {
    pub(super) path: Arc<dyn Proven>,
    pub(super) range: [f64; 2],
}

impl Range {
    fn map(&self, t: f64) -> f64 {
        self.range[0] + (self.range[1] - self.range[0]) * t
    }
}

impl Proven for Range {
    fn breaks(&self) -> Vec<f64> {
        let [from, to] = self.range;
        let mut output = vec![0.];
        output.extend(
            self.path
                .breaks()
                .into_iter()
                .map(|t| (t - from) / (to - from))
                .filter(|t| *t > 0. && *t < 1.),
        );
        output.push(1.);
        output
    }

    fn bound(&self, [a, b]: [f64; 2]) -> f64 {
        self.path.bound([self.map(a), self.map(b)])
    }
}

impl Path for Range {
    fn point(&self, t: f64) -> P3 {
        self.path.point(self.map(t))
    }

    fn flatten(&self, tolerance_mm: f64) -> Vec<f64> {
        certified(self, tolerance_mm)
    }

    fn paper(&self, [a, b]: [f64; 2]) -> Option<Curve> {
        let span = self.range[1] - self.range[0];
        self.path
            .paper([self.range[0] + a * span, self.range[0] + b * span])
    }
}

/// The meridian lines of a plane `n . p = c` cut through a law that
/// advances `k` per radian axially and `tau k` per radian radially: at
/// `theta` the plane meets the meridian half-plane in `alpha W + beta u =
/// gamma`, with `W` the profile radius and `u` the profile phase.
#[derive(Clone, Copy, Debug)]
pub(super) struct Meridian {
    pub n: P3,
    pub c: f64,
    pub k: f64,
    pub tau: f64,
}

impl Meridian {
    pub(super) fn line(&self, theta: f64) -> [f64; 3] {
        let (sin, cos) = theta.portable_sin_cos();
        let n = self.n;
        let alpha = n[0] * cos + n[1] * sin;
        [
            alpha,
            n[2],
            self.c - (alpha * self.tau + n[2]) * self.k * theta,
        ]
    }

    /// The plane's signed value `n . P - c` at the surface point over the
    /// profile point `q = (u, W)` at `theta`.
    fn value(&self, theta: f64, q: P2) -> f64 {
        let [alpha, beta, gamma] = self.line(theta);
        alpha * q[1] + beta * q[0] - gamma
    }

    /// `d value / d theta` at a fixed profile point.
    fn slope(&self, theta: f64, q: P2) -> f64 {
        let (sin, cos) = theta.portable_sin_cos();
        let n = self.n;
        let alpha = n[0] * cos + n[1] * sin;
        let turn = n[1] * cos - n[0] * sin;
        turn * (q[1] + self.tau * self.k * theta) + (alpha * self.tau + n[2]) * self.k
    }
}

/// A touching end of a section branch on an arc: there the meridian line is
/// tangent to the circle and the branch's root in `theta` has a square-root
/// singularity, so near it the branch is parametrised by the circle angle
/// `beta` instead. `theta(beta)` is the simple root of the plane's value at
/// the circle point, which behaves as `theta_t + (beta - beta_t)^2 / Q`: both
/// branches leaving a touch meet exactly at the touch point.
#[derive(Clone, Copy, Debug)]
struct Touch {
    /// Circle angles of the touch point and of the switch to `theta`.
    beta: [f64; 2],
    /// The bracket `[no-root side of the touch, switch]` in `theta`, on which
    /// the value is certified strictly monotone for every `beta` between the
    /// two angles, with opposite signs at its ends.
    bracket: [f64; 2],
}

fn wrap(angle: f64) -> f64 {
    angle - TAU * (angle / TAU).round()
}

impl Touch {
    /// The touch end at `touch` of the `index`-th root branch on an arc,
    /// whose other end is at `other`. `None` when the touch cannot be
    /// certified simple; the branch then keeps the foot point as its end.
    fn new(meridian: &Meridian, span: &Span, index: usize, touch: f64, other: f64) -> Option<Self> {
        let Span::Arc { center, radius, .. } = *span else {
            return None;
        };
        let at = |beta: f64| {
            let (sin, cos) = beta.portable_sin_cos();
            [center[0] + radius * cos, center[1] + radius * sin]
        };
        let angle = |p: P2| (p[1] - center[1]).portable_atan2(p[0] - center[0]);
        let reach = |theta: f64| {
            let [alpha, beta, gamma] = meridian.line(theta);
            let norm = alpha.portable_hypot(beta);
            let delta = (beta * center[0] + alpha * center[1] - gamma) / norm;
            (radius - delta.abs()) * (radius + delta.abs())
        };
        let [alpha, beta, gamma] = meridian.line(touch);
        let foot = span.arc_root(alpha, beta, gamma, 2, true)?;
        let beta_t = angle(foot);
        let s0 = meridian.slope(touch, at(beta_t));
        let direction = (other - touch).signum();
        // The no-root side: the line misses the circle just before the touch.
        let mut mu = 1e-12 * touch.abs().max(1.);
        while reach(touch - direction * mu) >= 0. {
            mu *= 2.;
            if mu > 1e-3 * (other - touch).abs() {
                return None;
            }
        }
        let outside = touch - direction * mu;
        let nxy = meridian.n[0].portable_hypot(meridian.n[1]);
        let lead = meridian.tau.abs() * meridian.k;
        let mut fraction = 0.45;
        for _ in 0..48 {
            let switch = touch + fraction * (other - touch);
            fraction *= 0.5;
            let [alpha, beta, gamma] = meridian.line(switch);
            let Some(p) = span.arc_root(alpha, beta, gamma, index, false) else {
                continue;
            };
            let beta_w = beta_t + wrap(angle(p) - beta_t);
            if (beta_w - beta_t).abs() <= 1e-9 {
                break;
            }
            // `|d slope / d theta| <= nxy (|W| + tau k |theta|) + 2 nxy tau k`
            // and `|d slope / d beta| <= nxy R` over the box.
            let theta_max = outside.abs().max(switch.abs());
            let w_max = center[1].abs() + radius + lead * theta_max;
            let reach_theta = (switch - touch).abs().max(mu);
            let drift = (nxy * w_max + 2. * nxy * lead) * reach_theta
                + nxy * radius * (beta_w - beta_t).abs();
            if drift > 0.5 * s0.abs() - 1e-12 * (1. + s0.abs()) {
                continue;
            }
            // Opposite signs across the bracket for the touch and the switch
            // angles (the latter at its own root, so check just inside).
            let inner = beta_t + 0.5 * (beta_w - beta_t);
            let [lo, hi] = [outside, switch].map(|theta| meridian.value(theta, at(inner)));
            if (lo > 0.) == (hi > 0.) {
                continue;
            }
            return Some(Self {
                beta: [beta_t, beta_w],
                bracket: [outside, switch],
            });
        }
        None
    }

    /// The root `theta` and the circle point at `beta`.
    fn locate(&self, meridian: &Meridian, span: &Span, beta: f64) -> (f64, P2) {
        let Span::Arc { center, radius, .. } = *span else {
            unreachable!("touches are on arcs")
        };
        let (sin, cos) = beta.portable_sin_cos();
        let q = [center[0] + radius * cos, center[1] + radius * sin];
        let f = |theta: f64| meridian.value(theta, q);
        let [mut lo, mut hi] = self.bracket;
        let [mut flo, fhi] = [f(lo), f(hi)];
        if (flo > 0.) == (fhi > 0.) {
            // Only at the switch angle itself, whose root is the bracket end.
            return (if flo.abs() < fhi.abs() { lo } else { hi }, q);
        }
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if mid == lo || mid == hi {
                break;
            }
            let fm = f(mid);
            if (fm > 0.) == (flo > 0.) {
                lo = mid;
                flo = fm;
            } else {
                hi = mid;
            }
        }
        (0.5 * (lo + hi), q)
    }
}

/// The parametrisation of a section branch: for each `theta` of `range` the
/// `index`-th meet of the plane's meridian line with the span, switching to
/// the circle angle near touching ends (see [`Touch`]).
pub(super) struct Branch {
    span: Span,
    index: usize,
    meridian: Meridian,
    range: [f64; 2],
    /// Ends at which the plane's meridian line touches the arc: both roots
    /// meet there, so the end is the touching point itself.
    touching: [bool; 2],
    touches: [Option<Touch>; 2],
}

/// The share of the branch parameter given to a touching end's angle range.
const TOUCH_SHARE: f64 = 0.25;

impl Branch {
    pub(super) fn new(
        span: Span,
        index: usize,
        meridian: Meridian,
        range: [f64; 2],
        touching: [bool; 2],
    ) -> Self {
        let touches = [0, 1].map(|end| {
            touching[end]
                .then(|| Touch::new(&meridian, &span, index, range[end], range[1 - end]))
                .flatten()
        });
        Self {
            span,
            index,
            meridian,
            range,
            touching,
            touches,
        }
    }

    fn shares(&self) -> [f64; 2] {
        self.touches
            .map(|touch| if touch.is_some() { TOUCH_SHARE } else { 0. })
    }

    /// The `theta` range swept by plain roots.
    fn plain(&self) -> [f64; 2] {
        [0, 1].map(|end| self.touches[end].map_or(self.range[end], |touch| touch.bracket[1]))
    }

    fn root(&self, theta: f64, touch: bool) -> P2 {
        let [alpha, beta, gamma] = self.meridian.line(theta);
        match self.span {
            Span::Line { a, b } => {
                let [fa, fb] = [a, b].map(|p| alpha * p[1] + beta * p[0] - gamma);
                let t = if fa == fb {
                    0.5
                } else {
                    (fa / (fa - fb)).clamp(0., 1.)
                };
                [a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]
            }
            Span::Arc { .. } => self
                .span
                .arc_root(alpha, beta, gamma, if touch { 2 } else { self.index }, true)
                .unwrap_or([0.; 2]),
        }
    }

    /// The angle and profile point `(u, W)` at the branch parameter `t`.
    pub(super) fn locate(&self, t: f64) -> (f64, P2) {
        let [a, b] = self.shares();
        let [from, to] = self.plain();
        if let Some(touch) = self.touches[0].filter(|_| t < a) {
            let s = (t / a).max(0.);
            let beta = touch.beta[0] + s * (touch.beta[1] - touch.beta[0]);
            return touch.locate(&self.meridian, &self.span, beta);
        }
        if let Some(touch) = self.touches[1].filter(|_| t > 1. - b) {
            let s = ((1. - t) / b).max(0.);
            let beta = touch.beta[0] + s * (touch.beta[1] - touch.beta[0]);
            return touch.locate(&self.meridian, &self.span, beta);
        }
        if t <= 0. {
            return (self.range[0], self.root(self.range[0], self.touching[0]));
        }
        if t >= 1. {
            return (self.range[1], self.root(self.range[1], self.touching[1]));
        }
        let s = ((t - a) / (1. - a - b)).clamp(0., 1.);
        let theta = from + s * (to - from);
        (theta, self.root(theta, false))
    }

    /// The branch parameter of a point at `theta` (any turn) whose profile
    /// point there is `profile(theta)`, when it lies on the branch's range.
    pub(super) fn parameter(&self, theta: f64, profile: &dyn Fn(f64) -> P2) -> Option<f64> {
        let [a, b] = self.shares();
        let [from, to] = self.plain();
        let [lo, hi] = [from.min(to), from.max(to)];
        for turn in [-1., 0., 1.] {
            let theta = theta + TAU * (turn + ((lo - theta) / TAU).floor());
            if theta > lo && theta < hi {
                return Some(a + (theta - from) / (to - from) * (1. - a - b));
            }
        }
        let Span::Arc { center, .. } = self.span else {
            return None;
        };
        for end in [0, 1] {
            let Some(touch) = self.touches[end] else {
                continue;
            };
            let near = touch.bracket[0];
            let q = profile(theta + TAU * ((near - theta) / TAU).round());
            let beta = (q[1] - center[1]).portable_atan2(q[0] - center[0]);
            let s = wrap(beta - touch.beta[0]) / (touch.beta[1] - touch.beta[0]);
            if s > 0. && s <= 1. {
                return Some(if end == 0 { s * a } else { 1. - s * b });
            }
        }
        None
    }

    /// Parameters between the touch parts and the plain part.
    pub(super) fn breaks(&self) -> Vec<f64> {
        let [a, b] = self.shares();
        let mut output = vec![0.];
        if a > 0. {
            output.push(a);
        }
        if b > 0. {
            output.push(1. - b);
        }
        output.push(1.);
        output
    }

    /// `alpha`, `gamma`, the signed distance of the arc's center from the
    /// meridian line and the line's normal length, at an angle.
    fn line_iv(&self, theta: Iv) -> [Iv; 2] {
        let m = &self.meridian;
        let alpha = Iv::point(m.n[0]) * theta.cos() + Iv::point(m.n[1]) * theta.sin();
        let gamma = Iv::point(m.c)
            - (alpha * Iv::point(m.tau) + Iv::point(m.n[2])) * Iv::point(m.k) * theta;
        [alpha, gamma]
    }

    /// The signed distance of the arc's center from the meridian line, with
    /// the line's normal length and `alpha`.
    fn offset_iv(&self, theta: Iv) -> Option<[Iv; 3]> {
        let Span::Arc { center, .. } = self.span else {
            return None;
        };
        let [alpha, gamma] = self.line_iv(theta);
        let beta = Iv::point(self.meridian.n[2]);
        let norm = (alpha.sqr() + beta.sqr()).sqrt();
        let delta = (beta * Iv::point(center[0]) + alpha * Iv::point(center[1]) - gamma) / norm;
        Some([delta, norm, alpha])
    }

    /// The exact chart point `(theta, u)` of a touching end: the touch
    /// angle is enclosed by a proven sign change of the reach `R^2 -
    /// delta^2`, and the foot point over that enclosure.
    fn touching_end(&self, end: usize) -> Option<[Iv; 2]> {
        let Span::Arc { center, radius, .. } = self.span else {
            return None;
        };
        let r2 = Iv::point(radius).sqr();
        let theta = bracket(self.range[end], |theta| {
            let [delta, ..] = self.offset_iv(theta).expect("arc");
            r2 - delta.sqr()
        })?;
        let [delta, norm, _] = self.offset_iv(theta)?;
        let u = Iv::point(center[0]) - delta * Iv::point(self.meridian.n[2]) / norm;
        Some([theta, u])
    }

    /// The exact profile phase of the branch's root at an angle, in closed
    /// form.
    fn root_iv(&self, theta: Iv) -> Option<Iv> {
        match self.span {
            Span::Line { a, b } => {
                let [alpha, gamma] = self.line_iv(theta);
                let n2 = Iv::point(self.meridian.n[2]);
                let value = |p: P2| alpha * Iv::point(p[1]) + n2 * Iv::point(p[0]) - gamma;
                let [fa, fb] = [value(a), value(b)];
                let t = fa / (fa - fb);
                Some(Iv::point(a[0]) + t * (Iv::point(b[0]) - Iv::point(a[0])))
            }
            Span::Arc { center, radius, .. } => {
                let [delta, norm, alpha] = self.offset_iv(theta)?;
                let h = (Iv::point(radius).sqr() - delta.sqr()).sqrt();
                let side = Iv::point(if self.index == 0 { -1. } else { 1. });
                Some(
                    Iv::point(center[0])
                        - delta * Iv::point(self.meridian.n[2]) / norm
                        - side * h * alpha / norm,
                )
            }
        }
    }

    /// An exact chart point `(theta, u)` next to the published point at `t`,
    /// in the mode of `part`.
    fn enclose(&self, chart: &BranchChart, t: f64, part: Option<usize>) -> Option<[Iv; 2]> {
        let [a, b] = self.shares();
        if let Some(end) = part {
            let touch = self.touches[end]?;
            let s = if end == 0 { t / a } else { (1. - t) / b }.max(0.);
            if s == 0. {
                return self.touching_end(end);
            }
            let Span::Arc { center, radius, .. } = self.span else {
                return None;
            };
            // As `locate`: the circle point is exact up to its cosine, and
            // the certified monotone bracket holds one root for it.
            let beta = touch.beta[0] + s * (touch.beta[1] - touch.beta[0]);
            let (theta, _) = touch.locate(&self.meridian, &self.span, beta);
            let u = Iv::point(center[0]) + Iv::point(radius) * Iv::point(beta).cos();
            let theta = bracket(theta, |theta| self.relation(theta, u))?;
            return Some([theta, u]);
        }
        for end in [0, 1] {
            if self.touching[end] && (if end == 0 { t <= 0. } else { t >= 1. }) {
                return self.touching_end(end);
            }
        }
        // A proven sign change next to the computed point, else the root at
        // the exact angle in closed form.
        let (theta, p) = self.locate(t);
        settle(chart, [theta, p[0]]).or_else(|| {
            let theta = Iv::point(theta);
            Some([theta, self.root_iv(theta)?])
        })
    }

    /// Proven bound of the deviation of the branch, placed in a member's
    /// frame at a chart base angle, over `[ta, tb]` inside one break
    /// interval, from the chord of `points`.
    pub(super) fn bound(
        &self,
        (origin, axes, base): (P3, [P3; 3], f64),
        [ta, tb]: [f64; 2],
        points: [P3; 2],
    ) -> f64 {
        let [a, b] = self.shares();
        // The part of a parameter, as in `locate`.
        let part = |t: f64| {
            if self.touches[0].is_some() && t < a {
                Some(0)
            } else if self.touches[1].is_some() && t > 1. - b {
                Some(1)
            } else {
                None
            }
        };
        let chart = BranchChart {
            branch: self,
            origin,
            axes,
            base,
        };
        // Each stop is enclosed in the part of its own parameter, so the
        // chords on both sides of a switch between parts share its exact
        // point and their pieces tile the branch.
        let enclose = |t| self.enclose(&chart, t, part(t));
        let (Some(start), Some(end)) = (enclose(ta), enclose(tb)) else {
            return f64::INFINITY;
        };
        // Plain parts are graphs over the angle, touch parts over the
        // circle angle and so over the phase.
        let native = if part(0.5 * (ta + tb)).is_some() {
            1
        } else {
            0
        };
        chord_bound(&chart, native, [start, end], points)
    }

    /// The plane's value at the profile point of phase `u` at `theta`.
    fn relation<N: Num>(&self, theta: N, u: N) -> N {
        let m = &self.meridian;
        let alpha = N::k(m.n[0]) * theta.cos() + N::k(m.n[1]) * theta.sin();
        let gamma = N::k(m.c) - (alpha * N::k(m.tau) + N::k(m.n[2])) * N::k(m.k) * theta;
        let [w, _] = self.span.radius_at(u);
        alpha * w + N::k(m.n[2]) * u - gamma
    }
}

/// A section branch in the chart `(theta, u)`: the zero set of the plane's
/// value at the profile point of phase `u`, placed at `r = W + tau k theta`,
/// angle `base + theta`, `z = u + k theta` in the member's frame.
struct BranchChart<'a> {
    branch: &'a Branch,
    origin: P3,
    axes: [P3; 3],
    base: f64,
}

impl Chart for BranchChart<'_> {
    fn implicit(&self) -> bool {
        true
    }

    fn relation<N: Num>(&self, theta: N, u: N) -> N {
        self.branch.relation(theta, u)
    }

    fn place<N: Num>(&self, theta: N, u: N) -> [N; 3] {
        let m = &self.branch.meridian;
        let [w, _] = self.branch.span.radius_at(u);
        let s = N::k(m.k) * theta;
        let r = w + N::k(m.tau) * s;
        let angle = theta + N::k(self.base);
        frame(
            self.origin,
            self.axes,
            [r * angle.cos(), r * angle.sin(), u + s],
        )
    }
}

/// A branch of an oblique section on a shifted span.
struct Sweep {
    law: Arc<Law>,
    /// The chart's angle origin: the branch's angles are measured from it.
    base: f64,
    branch: Branch,
}

impl Sweep {
    fn line(law: &Law, (n, c): (P3, f64), theta: f64) -> [f64; 3] {
        Self::meridian(law, (n, c)).line(theta)
    }

    fn meridian(law: &Law, (n, c): (P3, f64)) -> Meridian {
        Meridian {
            n,
            c,
            k: law.rate,
            tau: 0.,
        }
    }

    /// The profile point of a local point at the chart angle `theta`.
    fn profile(&self, p: P3) -> impl Fn(f64) -> P2 + '_ {
        move |theta| [p[2] - self.law.rate * theta, p[0].portable_hypot(p[1])]
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
        let z = p[0] + self.law.rate * theta;
        self.law.surface(p[1], self.base + theta, z)
    }

    fn flatten(&self, tolerance_mm: f64) -> Vec<f64> {
        certified(self, tolerance_mm)
    }

    fn paper(&self, _: [f64; 2]) -> Option<Curve> {
        None
    }
}

/// A screw thread member, with an optional terminal neck.
pub(in crate::drawing) struct Thread {
    law: Arc<Law>,
    /// Material z range of the threaded part.
    z: [f64; 2],
    neck: Option<(f64, [f64; 2], Prism)>,
    /// For a threaded bore: the outer stock radius and its cylinder. The
    /// material is then outside the helicoid and inside the stock.
    stock: Option<(f64, Prism)>,
    /// For an entry chamfer: the cone trimming the crests near the entry.
    cone: Option<Cone>,
    /// For a screw: the shank and head, a prism ending at the thread's start
    /// level whose cross-section there contains the thread's.
    base: Option<Prism>,
    /// For a treated bore with a terminal groove: the void past the station
    /// is the tail cylinder, and the shoulder annulus closes it.
    groove: Option<Groove>,
    /// For a hex nut: the hexagonal stock replacing the stock cylinder.
    hex: Option<multicad_rust_kernel::hex_nut::Nut>,
    /// For a partial-depth internal thread: the stock top beyond the void,
    /// whose floor at `z[1]` is the full stock disc.
    floor: Option<f64>,
    tolerance: f64,
}

/// A treated bore's terminal groove: the void is `r < radius` for
/// `station <= z <= top`, beyond every radius of the thread's void.
#[derive(Clone, Copy, Debug)]
struct Groove {
    radius: f64,
    station: f64,
    top: f64,
}

/// The entry cone `r = entry + slope z` on `0 <= z <= length`; beyond it the
/// thread is untrimmed. On an external thread the material is within the
/// cone; on an internal one (a treated bore) the void is the union of the
/// helicoid's and the cone's, so the material is outside the cone.
#[derive(Clone, Copy, Debug)]
struct Cone {
    entry: f64,
    length: f64,
    slope: f64,
    internal: bool,
}

impl Cone {
    fn radius(&self, z: f64) -> f64 {
        self.entry + self.slope * z
    }

    /// The seam station where the helicoid radius `r` meets the cone.
    fn station(&self, r: f64) -> f64 {
        ((r - self.entry) / self.slope).clamp(0., self.length)
    }

    /// Whether the helicoid point at radius `r` and level `z` bounds the
    /// material, i.e. is not trimmed away by the cone.
    fn keeps(&self, r: f64, z: f64, tolerance: f64) -> bool {
        z >= self.length - tolerance
            || if self.internal {
                r >= self.radius(z) - tolerance
            } else {
                r <= self.radius(z) + tolerance
            }
    }

    /// Whether a helix of radius `r` is trimmed near the entry: it then
    /// starts at its seam station.
    fn trims(&self, r: f64) -> bool {
        if self.internal {
            r < self.entry
        } else {
            r > self.entry
        }
    }

    /// Whether the cone's own surface at radius `cone` bounds the material
    /// where the helicoid's radius is `thread`.
    fn face(&self, cone: f64, thread: f64) -> bool {
        if self.internal {
            cone > thread
        } else {
            cone < thread
        }
    }

    /// Parameters of the unit-direction local line `s + t d` inside the
    /// cone's solid: past its end, or within its radius.
    fn line(&self, s: P3, d: P3) -> Vec<[f64; 2]> {
        let far = 1e9;
        let mut output = Vec::new();
        // Past the cone's end.
        if d[2].abs() <= 1e-15 {
            if s[2] >= self.length {
                output.push([-far, far]);
            }
        } else {
            let t = (self.length - s[2]) / d[2];
            output.push(if d[2] > 0. { [t, far] } else { [-far, t] });
        }
        output.extend(self.inside(s, d));
        merge(output, 0.)
    }

    /// Parameters of the local line within the cone frustum `0 < z < length`,
    /// `r < C(z)`: the cone's part of a treated bore's void.
    fn frustum(&self, s: P3, d: P3) -> Vec<[f64; 2]> {
        let far = 1e9;
        let slab = if d[2].abs() <= 1e-15 {
            if s[2] > 0. && s[2] < self.length {
                [-far, far]
            } else {
                return Vec::new();
            }
        } else {
            let [a, b] = [-s[2] / d[2], (self.length - s[2]) / d[2]];
            [a.min(b), a.max(b)]
        };
        profile::intersect(&self.inside(s, d), &[slab])
    }

    /// Parameters of the local line where `r < C(z)` on the cone's nappe of
    /// positive radius.
    fn inside(&self, s: P3, d: P3) -> Vec<[f64; 2]> {
        let far = 1e9;
        // r^2 < C(z)^2 with C(z) > 0.
        let k = self.slope;
        let c0 = self.entry + k * s[2];
        let a = d[0] * d[0] + d[1] * d[1] - k * k * d[2] * d[2];
        let b = 2. * (s[0] * d[0] + s[1] * d[1]) - 2. * k * d[2] * c0;
        let c = s[0] * s[0] + s[1] * s[1] - c0 * c0;
        let mut inside: Vec<[f64; 2]> = Vec::new();
        if a.abs() <= 1e-15 {
            if b.abs() <= 1e-15 {
                if c < 0. {
                    inside.push([-far, far]);
                }
            } else {
                let t = -c / b;
                inside.push(if b > 0. { [-far, t] } else { [t, far] });
            }
        } else {
            let disc = b * b - 4. * a * c;
            if disc > 0. {
                let root = disc.sqrt();
                let [t1, t2] = [(-b - root) / (2. * a), (-b + root) / (2. * a)];
                let [t1, t2] = [t1.min(t2), t1.max(t2)];
                if a > 0. {
                    inside.push([t1, t2]);
                } else {
                    inside.push([-far, t1]);
                    inside.push([t2, far]);
                }
            } else if a < 0. {
                inside.push([-far, far]);
            }
        }
        // Keep the nappe with a positive radius.
        let positive = if d[2].abs() <= 1e-15 {
            if c0 > 0. { [-far, far] } else { [0., 0.] }
        } else {
            let t = -c0 / (k * d[2]);
            if k * d[2] > 0. { [t, far] } else { [-far, t] }
        };
        profile::intersect(&merge(inside, 0.), &[positive])
    }
}

/// A seam where the helicoid meets the entry cone: for each phase in range
/// the point at the cone station of its radius.
struct Seam {
    law: Arc<Law>,
    span: Span,
    cone: Cone,
    copy: f64,
    range: [f64; 2],
}

impl Seam {
    fn phi(&self, t: f64) -> f64 {
        self.range[0] + (self.range[1] - self.range[0]) * t
    }

    fn at(&self, phi: f64) -> (P3, f64, f64) {
        let r = self.span.eval(phi)[0];
        let z = self.cone.station(r);
        let theta = (z - phi - self.copy) / self.law.rate;
        ([r, theta, z], r, theta)
    }
}

/// The seam is the explicit graph of its phase.
impl Chart for Seam {
    fn implicit(&self) -> bool {
        false
    }

    fn relation<N: Num>(&self, _: N, _: N) -> N {
        N::k(0.)
    }

    fn place<N: Num>(&self, phi: N, _: N) -> [N; 3] {
        let [r, _] = self.span.radius_at(phi);
        let z = (r - N::k(self.cone.entry)) / N::k(self.cone.slope);
        let theta = (z - phi - N::k(self.copy)) / N::k(self.law.rate);
        frame(
            self.law.origin,
            self.law.axes,
            [r * theta.cos(), r * theta.sin(), z],
        )
    }
}

impl Proven for Seam {
    fn bound(&self, [a, b]: [f64; 2]) -> f64 {
        let end = |t: f64| [Iv::point(self.phi(t)), Iv::point(0.)];
        chord_bound(self, 0, [end(a), end(b)], [self.point(a), self.point(b)])
    }
}

impl Path for Seam {
    fn point(&self, t: f64) -> P3 {
        let ([r, theta, z], _, _) = self.at(self.phi(t));
        self.law.surface(r, theta, z)
    }

    fn flatten(&self, tolerance_mm: f64) -> Vec<f64> {
        certified(self, tolerance_mm)
    }

    fn paper(&self, [from, to]: [f64; 2]) -> Option<Curve> {
        // On a cylindrical band the seam is a circle at the cone's end.
        let r = self.span.flat()?;
        let law = &self.law;
        let station = self.cone.station(r);
        let c = law.view([0., 0., station]);
        let [x, y] = [law.axes[0], law.axes[1]];
        let theta = |t: f64| {
            let phi = self.range[0] + (self.range[1] - self.range[0]) * t;
            (station - phi - self.copy) / law.rate
        };
        let [a, b] = [theta(from), theta(to)];
        paths::conic(
            [c[0], c[1]],
            [r * x[0], r * x[1]],
            [r * y[0], r * y[1]],
            a,
            b - a,
        )
    }
}

/// The section of the entry cone by the local plane `n . p = c`, over the
/// angles in `range`: `z = (c - q entry) / (q slope + n_z)`.
struct ConeCut {
    law: Arc<Law>,
    cone: Cone,
    plane: (P3, f64),
    range: [f64; 2],
}

impl ConeCut {
    fn local(cone: &Cone, (n, c): (P3, f64), theta: f64) -> Option<P3> {
        let q = n[0] * theta.portable_cos() + n[1] * theta.portable_sin();
        let denominator = q * cone.slope + n[2];
        if denominator.abs() <= 1e-14 {
            return None;
        }
        let z = (c - q * cone.entry) / denominator;
        let r = cone.radius(z);
        Some([r * theta.portable_cos(), r * theta.portable_sin(), z])
    }
}

/// The cone section is the explicit graph of its angle.
impl Chart for ConeCut {
    fn implicit(&self) -> bool {
        false
    }

    fn relation<N: Num>(&self, _: N, _: N) -> N {
        N::k(0.)
    }

    fn place<N: Num>(&self, theta: N, _: N) -> [N; 3] {
        let (n, c) = self.plane;
        let cone = &self.cone;
        let (sin, cos) = (theta.sin(), theta.cos());
        let q = N::k(n[0]) * cos + N::k(n[1]) * sin;
        let z = (N::k(c) - q * N::k(cone.entry)) / (q * N::k(cone.slope) + N::k(n[2]));
        let r = N::k(cone.entry) + N::k(cone.slope) * z;
        frame(self.law.origin, self.law.axes, [r * cos, r * sin, z])
    }
}

impl ConeCut {
    fn theta(&self, t: f64) -> f64 {
        self.range[0] + (self.range[1] - self.range[0]) * t
    }
}

impl Proven for ConeCut {
    fn bound(&self, [a, b]: [f64; 2]) -> f64 {
        let end = |t: f64| [Iv::point(self.theta(t)), Iv::point(0.)];
        chord_bound(self, 0, [end(a), end(b)], [self.point(a), self.point(b)])
    }
}

impl Path for ConeCut {
    fn point(&self, t: f64) -> P3 {
        let p = Self::local(&self.cone, self.plane, self.theta(t)).unwrap_or([f64::NAN; 3]);
        self.law.view(p)
    }

    fn flatten(&self, tolerance_mm: f64) -> Vec<f64> {
        certified(self, tolerance_mm)
    }

    fn paper(&self, _: [f64; 2]) -> Option<Curve> {
        None
    }
}

pub(in crate::drawing) fn member(
    solid: &Solid,
    view: &View,
) -> Option<Result<Box<dyn Member>, GeometryError>> {
    build(solid, view).map(|built| built.map(|thread| Box::new(thread) as Box<dyn Member>))
}

fn build(solid: &Solid, view: &View) -> Option<Result<Thread, GeometryError>> {
    let invalid = || {
        err(
            GeometryErrorKind::InvalidResult,
            "Rust thread profile cannot be read for drawing",
        )
    };
    let tolerance = solid.geometry_policy().linear_tolerance_mm;
    let place = |frame: PlaneFrame| -> (P3, [P3; 3]) {
        (
            view.point(frame.origin),
            [frame.u, frame.v, frame.normal()].map(|axis| view.vector(axis)),
        )
    };
    if let Some(built) = fastener::build(solid, view, tolerance) {
        return Some(built);
    }
    if let Some((relief, frame)) = solid.rust_thread_end_relief_data() {
        let thread = relief.material_thread();
        let Some(profile) = Profile::from_helical(thread) else {
            return Some(Err(invalid()));
        };
        let (origin, axes) = place(frame);
        let neck = relief.neck_radius_mm();
        let span = [relief.material_length_mm(), relief.full_length_mm()];
        return Some((|| {
            let prism = Prism::new(
                origin,
                axes,
                vec![Band {
                    rings: vec![vec![Edge::Planar(normal_sweep::Edge::Circle {
                        center: [0.; 2],
                        radius_mm: neck,
                        clockwise: false,
                    })]],
                    axis: span,
                }],
                tolerance,
            )?;
            Ok(Thread {
                law: Arc::new(Law {
                    origin,
                    axes,
                    profile,
                    rate: thread.helical_rate(),
                    starts: thread.spec().starts,
                }),
                z: [0., span[0]],
                neck: Some((neck, span, prism)),
                stock: None,
                cone: None,
                base: None,
                groove: None,
                hex: None,
                floor: None,
                tolerance,
            })
        })());
    }
    if let Some((screw, frame)) = solid.rust_machine_screw_geometry() {
        let thread = screw.thread();
        let Some(profile) = Profile::from_helical(thread) else {
            return Some(Err(invalid()));
        };
        let (origin, axes) = place(frame);
        let spec = screw.spec();
        let top = screw.head_top_z();
        let seat = screw.head_seat_z();
        let outline = |head: multicad_rust_kernel::headed_screw::ScrewHead| -> Vec<Edge> {
            use multicad_rust_kernel::headed_screw::ScrewHead;
            match head {
                ScrewHead::SocketCap { radius_mm, .. } => {
                    vec![Edge::Planar(normal_sweep::Edge::Circle {
                        center: [0.; 2],
                        radius_mm,
                        clockwise: false,
                    })]
                }
                ScrewHead::ExternalHex {
                    across_flats_mm, ..
                } => hexagon(across_flats_mm),
            }
        };
        let head = outline(spec.head);
        let mut bands = Vec::new();
        match spec.head.socket() {
            Some((width, depth)) => {
                bands.push(Band {
                    rings: vec![head.clone(), hexagon(width)],
                    axis: [top, top + depth],
                });
                bands.push(Band {
                    rings: vec![head],
                    axis: [top + depth, seat],
                });
            }
            None => bands.push(Band {
                rings: vec![head],
                axis: [top, seat],
            }),
        }
        if spec.shank_length_mm > 0. {
            bands.push(Band {
                rings: vec![vec![Edge::Planar(normal_sweep::Edge::Circle {
                    center: [0.; 2],
                    radius_mm: spec.shank_radius_mm,
                    clockwise: false,
                })]],
                axis: [seat, 0.],
            });
        }
        return Some((|| {
            let base = Prism::new(origin, axes, bands, tolerance)?;
            Ok(Thread {
                law: Arc::new(Law {
                    origin,
                    axes,
                    profile,
                    rate: thread.helical_rate(),
                    starts: thread.spec().starts,
                }),
                z: [0., thread.spec().length_mm],
                neck: None,
                stock: None,
                cone: None,
                base: Some(base),
                groove: None,
                hex: None,
                floor: None,
                tolerance,
            })
        })());
    }
    if let Some(data) = solid.rust_thread_entry_chamfer_data() {
        let chamfer = data.source;
        let tool = chamfer.source_thread();
        let Some(profile) = Profile::from_helical(tool) else {
            return Some(Err(invalid()));
        };
        let (origin, axes) = place(data.frame);
        let spec = tool.spec();
        let material = chamfer.material_length_mm();
        let cone = Cone {
            entry: chamfer.entry_radius_mm(),
            length: chamfer.entry_length_mm(),
            slope: (spec.major_radius_mm - chamfer.entry_radius_mm()) / chamfer.entry_length_mm(),
            internal: false,
        };
        return Some((|| {
            let neck = match chamfer.terminal_relief() {
                Some(relief) => {
                    let radius = relief.neck_radius_mm();
                    let span = [material, chamfer.full_length_mm()];
                    let prism = Prism::new(
                        origin,
                        axes,
                        vec![Band {
                            rings: vec![vec![Edge::Planar(normal_sweep::Edge::Circle {
                                center: [0.; 2],
                                radius_mm: radius,
                                clockwise: false,
                            })]],
                            axis: span,
                        }],
                        tolerance,
                    )?;
                    Some((radius, span, prism))
                }
                None => None,
            };
            Ok(Thread {
                law: Arc::new(Law {
                    origin,
                    axes,
                    profile,
                    rate: tool.helical_rate(),
                    starts: spec.starts,
                }),
                z: [0., material],
                neck,
                stock: None,
                cone: Some(cone),
                base: None,
                groove: None,
                hex: None,
                floor: None,
                tolerance,
            })
        })());
    }
    if let Some(data) = solid.rust_internal_thread_treatment_data() {
        return treatment(&data, view, tolerance);
    }
    if let Some(data) = solid.rust_rounded_thread_bore_data() {
        let thread = data.bore.thread();
        let Some(profile) = Profile::from_rounded(thread) else {
            return Some(Err(invalid()));
        };
        let (origin, axes) = place(data.frame);
        let spec = thread.spec();
        let radius = data.bore.radius_mm();
        return Some((|| {
            let prism = Prism::new(
                origin,
                axes,
                vec![Band {
                    rings: vec![vec![Edge::Planar(normal_sweep::Edge::Circle {
                        center: [0.; 2],
                        radius_mm: radius,
                        clockwise: false,
                    })]],
                    axis: [0., spec.length_mm],
                }],
                tolerance,
            )?;
            Ok(Thread {
                law: Arc::new(Law {
                    origin,
                    axes,
                    profile,
                    rate: thread.helical_rate(),
                    starts: spec.starts,
                }),
                z: [0., spec.length_mm],
                neck: None,
                stock: Some((radius, prism)),
                cone: None,
                base: None,
                groove: None,
                hex: None,
                floor: None,
                tolerance,
            })
        })());
    }
    let data = solid.rust_rounded_thread_data()?;
    let Some(profile) = Profile::from_rounded(&data.thread) else {
        return Some(Err(invalid()));
    };
    let (origin, axes) = place(data.frame);
    let spec = data.thread.spec();
    Some(Ok(Thread {
        law: Arc::new(Law {
            origin,
            axes,
            profile,
            rate: data.thread.helical_rate(),
            starts: spec.starts,
        }),
        z: [0., spec.length_mm],
        neck: None,
        stock: None,
        cone: None,
        base: None,
        groove: None,
        hex: None,
        floor: None,
        tolerance,
    }))
}

/// Actual cylindrical stock minus an authored internal thread cavity: the
/// void is `r < max(R(z - k theta), C(z))` with the flat-band profile `R` and
/// the entry cone `C`, between the stock's caps.
fn treatment(
    data: &crate::InternalThreadTreatmentData,
    view: &View,
    tolerance: f64,
) -> Option<Result<Thread, GeometryError>> {
    let unsupported = |message: &str| Some(Err(err(GeometryErrorKind::Unsupported, message)));
    let cavity = data.material.cavity();
    let thread = cavity.source_thread();
    let spec = thread.spec();
    let pitch = spec.pitch_mm;
    let bands = thread.bands();
    if bands[0].phase_mm[0] != 0. || bands[4].phase_mm[1] != pitch {
        return Some(Err(err(
            GeometryErrorKind::InvalidResult,
            "Rust internal thread bands do not close over one pitch",
        )));
    }
    let profile = Profile::new(
        pitch,
        bands
            .iter()
            .map(|band| Span::Line {
                a: [band.phase_mm[0], band.radius_mm[0]],
                b: [band.phase_mm[1], band.radius_mm[1]],
            })
            .collect(),
    );
    let origin = view.point(data.frame.origin);
    let axes = [data.frame.u, data.frame.v, data.frame.normal()].map(|axis| view.vector(axis));
    let length = cavity.full_length_mm();
    let radius = data.material.stock().circle().radius().as_f64();
    let cone = cavity
        .entry_radius_mm()
        .zip(cavity.intent().entry_depth_mm)
        .map(|(entry, depth)| Cone {
            entry,
            length: depth,
            slope: (thread.minimum_radius_mm() - entry) / depth,
            internal: true,
        });
    let sign = if spec.left_hand { -1. } else { 1. };
    let circle = |radius_mm: f64, clockwise: bool| {
        vec![Edge::Planar(normal_sweep::Edge::Circle {
            center: [0.; 2],
            radius_mm,
            clockwise,
        })]
    };
    // A terminal groove: past its station the void is the tail cylinder,
    // which must clear the whole thread and stay inside the stock.
    let groove = match cavity
        .terminal_radius_mm()
        .zip(cavity.terminal_station_mm())
    {
        Some((tail, station)) => {
            if !(tail > thread.maximum_radius_mm()
                && tail < radius
                && station > 0.
                && station < length)
            {
                return unsupported(
                    "Drawing a treated bore whose terminal groove does not clear the thread inside the stock is not supported",
                );
            }
            if cone.is_some_and(|cone| cone.length >= station) {
                return unsupported(
                    "Drawing a treated bore whose entry cone reaches its terminal groove is not supported",
                );
            }
            Some(Groove {
                radius: tail,
                station,
                top: length,
            })
        }
        None => None,
    };
    let bands = match groove {
        Some(groove) => vec![
            Band {
                rings: vec![circle(radius, false)],
                axis: [0., groove.station],
            },
            Band {
                rings: vec![circle(radius, false), circle(groove.radius, true)],
                axis: [groove.station, length],
            },
        ],
        None => vec![Band {
            rings: vec![circle(radius, false)],
            axis: [0., length],
        }],
    };
    Some((|| {
        let prism = Prism::new(origin, axes, bands, tolerance)?;
        Ok(Thread {
            law: Arc::new(Law {
                origin,
                axes,
                profile,
                rate: sign * pitch * f64::from(spec.starts) / TAU,
                starts: spec.starts,
            }),
            z: [0., groove.map_or(length, |groove| groove.station)],
            neck: None,
            stock: Some((radius, prism)),
            cone,
            base: None,
            groove,
            hex: None,
            floor: None,
            tolerance,
        })
    })())
}

impl Thread {
    fn rate(&self) -> f64 {
        self.law.rate
    }

    /// The level where a helix or band of radius `r` starts: its seam
    /// station on an entry cone that trims it, else the first cap.
    fn start(&self, r: f64) -> f64 {
        match self.cone {
            Some(cone) if cone.trims(r) => self.z[0].max(cone.station(r)),
            _ => self.z[0],
        }
    }

    fn pitch(&self) -> f64 {
        self.law.profile.pitch
    }

    /// Integer copies `m` of the profile whose phases `phi + m * pitch` can
    /// reach the material for some angle in `[0, 2pi]`.
    fn copies(&self) -> std::ops::RangeInclusive<i64> {
        let turn = TAU * self.rate().abs();
        let low = ((self.z[0] - turn) / self.pitch()).floor() as i64 - 1;
        let high = ((self.z[1] + turn) / self.pitch()).ceil() as i64 + 1;
        low..=high
    }

    /// Material parameters of the threaded part along a local line.
    fn threaded(&self, s: P3, d: P3) -> Result<Vec<[f64; 2]>, GeometryError> {
        let profile = &self.law.profile;
        // The axial slab.
        let slab = if d[2].abs() <= 1e-15 {
            if s[2] <= self.z[0] || s[2] >= self.z[1] {
                return Ok(Vec::new());
            }
            [-1e9, 1e9]
        } else {
            let ends = self.z.map(|z| (z - s[2]) / d[2]);
            [ends[0].min(ends[1]), ends[0].max(ends[1])]
        };
        let rho = d[0].portable_hypot(d[1]);
        // Radial range: r(t)^2 = |s_xy + t d_xy|^2.
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
        let core = within(profile.minimum);
        let Some(outer) = within(profile.maximum) else {
            return Ok(Vec::new());
        };
        let span = [slab[0].max(outer[0]), slab[1].min(outer[1])];
        if span[1] <= span[0] {
            return Ok(Vec::new());
        }
        let g = |t: f64| {
            let p = add(s, scaled(d, t));
            let r = p[0].portable_hypot(p[1]);
            profile.radius(p[2] - self.rate() * p[1].portable_atan2(p[0])) - r
        };
        let mut stretches = vec![span];
        if let Some([a, b]) = core {
            // Inside the minimum radius the line is in material.
            let [a, b] = [a.max(span[0]), b.min(span[1])];
            if b > a {
                stretches = vec![[span[0], a], [b, span[1]]];
            }
        }
        let mut output = Vec::new();
        if let Some([a, b]) = core {
            let [a, b] = [a.max(span[0]), b.min(span[1])];
            if b > a {
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

    /// Parameters in `[a, b]`, where the radius of the unit-direction line
    /// `s + t d` is at least the profile minimum, at which the line crosses
    /// the helicoid.
    fn crossings(&self, s: P3, d: P3, [a, b]: [f64; 2]) -> Result<Vec<f64>, GeometryError> {
        let profile = &self.law.profile;
        let rate = self.rate();
        let rho = d[0].portable_hypot(d[1]);
        if rho <= 1e-15 {
            // Along the axis: a fixed angle and radius, a linear phase.
            return Ok(self.axial_crossings(s, d, [a, b]));
        }
        let g = |t: f64| {
            let p = add(s, scaled(d, t));
            profile.radius(p[2] - rate * p[1].portable_atan2(p[0])) - p[0].portable_hypot(p[1])
        };
        // First derivative, exact on the span containing the point.
        let slope = |t: f64| {
            let p = add(s, scaled(d, t));
            let r2 = p[0] * p[0] + p[1] * p[1];
            let r = r2.sqrt();
            let phi = p[2] - rate * p[1].portable_atan2(p[0]);
            let (index, base) = profile.locate(phi);
            let [_, profile_slope] = profile.spans[index].eval(phi - base);
            let dphi = d[2] - rate * (p[0] * d[1] - p[1] * d[0]) / r2;
            (
                profile_slope * dphi - (p[0] * d[0] + p[1] * d[1]) / r,
                phi,
                index,
                base,
            )
        };
        let phase_rate = d[2].abs() + rate.abs() * rho / profile.minimum;
        let lipschitz = (profile.slope * phase_rate + rho) * (1. + 1e-6) + 1e-12;
        let second = profile.curvature * phase_rate * phase_rate
            + profile.slope * rate.abs() * 2. * rho * rho / (profile.minimum * profile.minimum)
            + rho * rho / profile.minimum;
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
                    "Thread drawing root isolation exceeded its work budget",
                ));
            }
            let h = b - a;
            let change = (fa <= 0.) != (fb <= 0.);
            if !change && fa.abs() + fb.abs() > lipschitz * h {
                continue;
            }
            // Second order: from either end, within one span of the profile.
            let excluded = [(a, fa), (b, fb)].into_iter().any(|(t, value)| {
                let (derivative, phi, index, base) = slope(t);
                let [lo, hi] = profile.spans[index].range().map(|x| x + base);
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

    /// Exact crossings of a line parallel to the axis.
    fn axial_crossings(&self, s: P3, d: P3, [a, b]: [f64; 2]) -> Vec<f64> {
        let profile = &self.law.profile;
        let r = s[0].portable_hypot(s[1]);
        let theta = s[1].portable_atan2(s[0]);
        // phi(t) = s_z - rate * theta + t * d_z with d_z = +-1.
        let phi0 = s[2] - self.rate() * theta;
        let [p, q] = [phi0 + a * d[2], phi0 + b * d[2]];
        let [lo, hi] = [p.min(q), p.max(q)];
        let first = (lo / profile.pitch).floor() as i64 - 1;
        let last = (hi / profile.pitch).ceil() as i64 + 1;
        let mut output = Vec::new();
        for m in first..=last {
            for span in &profile.spans {
                let span = span.shifted(m as f64 * profile.pitch);
                let [ra, rb] = span.range();
                let mut push = |phi: f64| {
                    if phi >= ra && phi <= rb && phi > lo && phi < hi {
                        output.push((phi - phi0) / d[2]);
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

    fn line_local(&self, s: P3, d: P3) -> Result<Vec<[f64; 2]>, GeometryError> {
        if let Some((radius, [junction, _], _)) = &self.neck
            && d[2].abs() <= 1e-15
            && (s[2] - junction).abs() <= self.tolerance
        {
            // In the shoulder plane the union is open only where both the
            // thread below and the neck above have material: inside the neck.
            let rho = d[0].portable_hypot(d[1]);
            let b = s[0] * d[0] + s[1] * d[1];
            let c = s[0] * s[0] + s[1] * s[1] - radius * radius;
            if rho <= 1e-15 {
                return Ok(if c < 0. {
                    vec![[-1e9, 1e9]]
                } else {
                    Vec::new()
                });
            }
            let disc = b * b - rho * rho * c;
            return Ok(if disc > 0. {
                let root = disc.sqrt();
                vec![[(-b - root) / (rho * rho), (-b + root) / (rho * rho)]]
            } else {
                Vec::new()
            });
        }
        if let (Some(groove), Some((stock, _))) = (&self.groove, &self.stock)
            && d[2].abs() <= 1e-15
            && (s[2] - groove.station).abs() <= self.tolerance
        {
            // In the shoulder plane the material on both sides is the
            // annulus between the tail cylinder and the stock.
            let disc = |radius: f64| -> Vec<[f64; 2]> {
                let rho = d[0].portable_hypot(d[1]);
                let b = s[0] * d[0] + s[1] * d[1];
                let c = s[0] * s[0] + s[1] * s[1] - radius * radius;
                if rho <= 1e-15 {
                    return if c < 0. {
                        vec![[-1e9, 1e9]]
                    } else {
                        Vec::new()
                    };
                }
                let disc = b * b - rho * rho * c;
                if disc > 0. {
                    let root = disc.sqrt();
                    vec![[(-b - root) / (rho * rho), (-b + root) / (rho * rho)]]
                } else {
                    Vec::new()
                }
            };
            return Ok(profile::exclusive(&disc(*stock), &disc(groove.radius), 0.));
        }
        if let Some(spans) = self.floor_line(s, d) {
            return spans;
        }
        if self.base.is_some() && d[2].abs() <= 1e-15 && (s[2] - self.z[0]).abs() <= self.tolerance
        {
            // In the plane where the thread starts on a screw's shank or head,
            // whose cross-section contains the thread's, the union is open
            // exactly within the thread's cap region.
            return self.cap_line(s, d, self.z[0]);
        }
        let mut spans = self.threaded(s, d)?;
        if let Some(prism) = &self.base
            && let Some(more) = prism.line(self.law.view(s), self.view_vector(d))
        {
            spans.extend(more?);
            spans = merge(spans, 0.);
        }
        match &self.cone {
            // A treated bore's void: the helicoid's or the cone frustum's.
            Some(cone) if cone.internal => {
                spans.extend(cone.frustum(s, d));
                spans = merge(spans, 0.);
            }
            Some(cone) => spans = profile::intersect(&spans, &cone.line(s, d)),
            None => {}
        }
        if let Some((_, prism)) = &self.stock {
            // Inside the stock cylinder and outside the thread's helicoid.
            let outer = match prism.line(self.law.view(s), self.view_vector(d)) {
                Some(outer) => outer?,
                None => Vec::new(),
            };
            // An axis-parallel line on a cylindrical band of the void lies in
            // its surface over the band: that stretch is boundary, not
            // material, so a stretch whose middle is not beyond the void's
            // radius there is dropped.
            let axial = d[0].portable_hypot(d[1]) <= 1e-15;
            let r = s[0].portable_hypot(s[1]);
            let theta = s[1].portable_atan2(s[0]);
            let profile = &self.law.profile;
            let beyond = |t: f64| {
                let z = s[2] + t * d[2];
                if self.floor.is_some() && z > self.z[1] {
                    // Past a partial-depth void: the floor's full disc.
                    return true;
                }
                let mut void = profile.radius(z - self.rate() * theta);
                if let Some(groove) = self.groove.filter(|groove| z >= groove.station) {
                    void = groove.radius;
                }
                if let Some(cone) = self.cone.filter(|cone| z > 0. && z < cone.length) {
                    void = void.max(cone.radius(z));
                }
                r > void + 1e-12 * void.max(1.)
            };
            return Ok(profile::exclusive(&outer, &spans, 0.)
                .into_iter()
                .filter(|[a, b]| b > a && (!axial || beyond(0.5 * (a + b))))
                .collect());
        }
        if let Some((_, _, prism)) = &self.neck
            && let Some(neck) = prism.line(self.law.view(s), self.view_vector(d))
        {
            spans.extend(neck?);
        }
        Ok(merge(spans, 0.))
    }

    fn view_vector(&self, [x, y, z]: P3) -> P3 {
        std::array::from_fn(|i| {
            x * self.law.axes[0][i] + y * self.law.axes[1][i] + z * self.law.axes[2][i]
        })
    }

    fn local_vector(&self, v: P3) -> P3 {
        self.law.axes.map(|axis| dot3(axis, v))
    }

    fn rails(&self, output: &mut Vec<Candidate>) {
        let profile = &self.law.profile;
        let rate = self.rate();
        for (i, kind) in profile.rails.iter().enumerate() {
            let Some(kind) = *kind else { continue };
            let phi = profile.spans[i].range()[0];
            let radius = profile.spans[i].eval(phi)[0];
            self.helix(phi, radius, kind, output);
        }
        let _ = rate;
    }

    /// Every copy of the helix of constant phase `phi` within the material.
    fn helix(&self, phi: f64, radius: f64, kind: Kind, output: &mut Vec<Candidate>) {
        let rate = self.rate();
        for m in 0..self.law.starts {
            let base = phi + f64::from(m) * self.pitch();
            // On an entry chamfer the helix starts at its cone station.
            let low = self.start(radius);
            if self.z[1] - low <= 1e-12 {
                continue;
            }
            let ends = [low, self.z[1]].map(|z| (z - base) / rate);
            let [a, b] = [ends[0].min(ends[1]), ends[0].max(ends[1])];
            // Split where the radius turns perpendicular to the sight: there
            // the helix touches the contour of the cylinder through it, and
            // its visibility can change without a crossing in the paper.
            let view = self.local_vector([0., 0., 1.]);
            let mut stops = vec![a, b];
            if view[0].portable_hypot(view[1]) > 1e-12 {
                let fold = view[1].portable_atan2(view[0]) + PI / 2.;
                let mut k = ((a - fold) / PI).ceil();
                while fold + k * PI < b {
                    stops.push(fold + k * PI);
                    k += 1.;
                }
            }
            let count = ((b - a) / (PI / 2.)).ceil().max(1.) as usize;
            stops.extend((1..count).map(|k| a + (b - a) * k as f64 / count as f64));
            stops.sort_by(f64::total_cmp);
            stops.dedup_by(|x, y| (*x - *y).abs() <= 1e-12);
            for pair in stops.windows(2) {
                let range = [pair[0], pair[1]];
                output.push(Candidate {
                    kind,
                    path: Arc::new(Helix {
                        law: self.law.clone(),
                        radius,
                        base,
                        range,
                    }),
                });
            }
        }
    }

    /// The cap curve at `level` as pieces of at most a quarter turn.
    fn cap(&self, level: f64) -> Vec<Arc<dyn Path>> {
        let mut output: Vec<Arc<dyn Path>> = Vec::new();
        for (span, [a, b], _) in self.cap_pieces(level) {
            let count = ((b - a) / (PI / 2.)).ceil().max(1.) as usize;
            for k in 0..count {
                let range = [
                    a + (b - a) * k as f64 / count as f64,
                    a + (b - a) * (k + 1) as f64 / count as f64,
                ];
                output.push(Arc::new(Cap {
                    law: self.law.clone(),
                    span,
                    level,
                    range,
                }));
            }
        }
        output
    }

    fn outlines(&self, output: &mut Vec<Candidate>) -> Result<(), GeometryError> {
        let view = self.local_vector([0., 0., 1.]);
        let rho = view[0].portable_hypot(view[1]);
        let profile = &self.law.profile;
        if rho <= 1e-9 {
            // Along the axis the outlines are the helices of stationary radius.
            for span in &profile.spans {
                if span.flat().is_none()
                    && let Some(phi) = span.stationary()
                {
                    self.helix(phi, span.eval(phi)[0], Kind::Outline, output);
                }
            }
            return Ok(());
        }
        let mut ends = self.contour_on_caps(view)?;
        if self.cone.is_some() {
            ends.extend(self.contour_on_seams(view)?);
        }
        for span in &profile.spans {
            if span.flat().is_some() {
                // A cylindrical crest or root band: its outlines are the
                // generators where the radius is perpendicular to the sight.
                self.flat_outlines(span, view, output);
                continue;
            }
            // Phase ranges where the branches exist: D = A^2 + B^2 - C^2 >= 0
            // with A = -R rho, B = k R' rho, C = -R R' v_z, isolated with a
            // proven Lipschitz bound of D.
            let [a, b] = span.range();
            let k = self.rate();
            let (rmax, slope, curvature) = (profile.maximum, profile.slope, profile.curvature);
            let vz2 = view[2] * view[2];
            let discriminant = |phi: f64| {
                let [r, d] = span.eval(phi);
                r * r * rho * rho + k * k * d * d * rho * rho - r * r * d * d * vz2
            };
            let lipschitz = 2. * rmax * slope * rho * rho
                + 2. * k * k * slope * curvature * rho * rho
                + 2. * rmax * slope * (slope * slope + rmax * curvature) * vz2
                + 1e-12;
            let mut roots = Vec::new();
            isolate(&discriminant, [a, b], lipschitz, &mut roots)?;
            let mut stops = vec![a];
            stops.extend(roots.into_iter().filter(|x| *x > a && *x < b));
            stops.push(b);
            let ranges: Vec<[f64; 2]> = stops
                .windows(2)
                .filter(|pair| pair[1] > pair[0] && discriminant(0.5 * (pair[0] + pair[1])) >= 0.)
                .map(|pair| [pair[0], pair[1]])
                .collect();
            for sign in [-1., 1.] {
                for &range in &ranges {
                    self.outline_copies(span, sign, range, view, &ends, output);
                }
            }
        }
        Ok(())
    }

    /// Spans of the cap curve at `level` over `theta` in `[0, 2pi]`, each
    /// within one span of the profile.
    fn cap_pieces(&self, level: f64) -> Vec<(Span, [f64; 2], bool)> {
        let profile = &self.law.profile;
        let rate = self.rate();
        let trim = self
            .cone
            .filter(|cone| level < cone.length)
            .map(|cone| cone.radius(level));
        let mut cuts = vec![0., TAU];
        for m in self.copies() {
            for span in &profile.spans {
                let shifted = span.shifted(m as f64 * self.pitch());
                let mut phases = vec![shifted.range()[0]];
                if let Some(radius) = trim {
                    phases.extend(shifted.at_radius(radius));
                }
                for phi in phases {
                    let theta = (level - phi) / rate;
                    if theta > 0. && theta < TAU {
                        cuts.push(theta);
                    }
                }
            }
        }
        cuts.sort_by(f64::total_cmp);
        cuts.dedup_by(|a, b| (*a - *b).abs() <= 1e-14);
        cuts.windows(2)
            .map(|pair| {
                let mid = 0.5 * (pair[0] + pair[1]);
                let (index, base) = profile.locate(level - rate * mid);
                let span = profile.spans[index].shifted(base);
                match trim {
                    // Where the thread reaches beyond the cone, the cap ends
                    // on the cone's circle.
                    Some(radius)
                        if self.cone.is_some_and(|cone| {
                            cone.face(radius, span.eval(level - rate * mid)[0])
                        }) =>
                    {
                        let [a, b] = span.range();
                        (
                            Span::Line {
                                a: [a, radius],
                                b: [b, radius],
                            },
                            [pair[0], pair[1]],
                            true,
                        )
                    }
                    _ => (span, [pair[0], pair[1]], false),
                }
            })
            .collect()
    }

    /// Angles on each cap curve where the surface normal is perpendicular to
    /// the sight: there contour curves end on the cap. Isolated with a proven
    /// Lipschitz bound of `n . v` along the cap curve.
    fn contour_on_caps(&self, view: P3) -> Result<Vec<(f64, f64)>, GeometryError> {
        let profile = &self.law.profile;
        let k = self.rate();
        let rho = view[0].portable_hypot(view[1]);
        let alpha = view[1].portable_atan2(view[0]);
        let (rmax, slope, curvature) = (profile.maximum, profile.slope, profile.curvature);
        let lipschitz = rho * (k.abs() * slope + rmax)
            + rho * k.abs() * (k.abs() * curvature + slope)
            + view[2].abs() * k.abs() * (slope * slope + rmax * curvature)
            + 1e-12;
        let mut output = Vec::new();
        for level in self.z {
            for (span, range, cone) in self.cap_pieces(level) {
                if cone {
                    continue;
                }
                let g = |theta: f64| {
                    let [r, d] = span.eval(level - k * theta);
                    let psi = theta - alpha;
                    -r * rho * psi.portable_cos()
                        + k * d * rho * psi.portable_sin()
                        + r * d * view[2]
                };
                let mut roots = Vec::new();
                isolate(&g, range, lipschitz, &mut roots)?;
                output.extend(roots.into_iter().map(|theta| (level, theta)));
            }
        }
        Ok(output)
    }

    /// Phase sub-ranges of each span, with the span, where the helicoid
    /// reaches beyond the cone's entry radius: there it meets the cone.
    fn seam_ranges(&self) -> Vec<(Span, [f64; 2])> {
        let Some(cone) = self.cone else {
            return Vec::new();
        };
        let mut output = Vec::new();
        for span in &self.law.profile.spans {
            let [a, b] = span.range();
            let mut stops = vec![a];
            stops.extend(span.at_radius(cone.entry));
            stops.push(b);
            stops.sort_by(f64::total_cmp);
            for pair in stops.windows(2) {
                if pair[1] - pair[0] > 1e-13 && cone.trims(span.eval(0.5 * (pair[0] + pair[1]))[0])
                {
                    output.push((*span, [pair[0], pair[1]]));
                }
            }
        }
        output
    }

    /// Seam curves between the helicoid and the entry cone.
    fn seams(&self) -> Vec<Seam> {
        let Some(cone) = self.cone else {
            return Vec::new();
        };
        let mut output = Vec::new();
        for (span, range) in self.seam_ranges() {
            for m in 0..self.law.starts {
                output.push(Seam {
                    law: self.law.clone(),
                    span,
                    cone,
                    copy: f64::from(m) * self.pitch(),
                    range,
                });
            }
        }
        output
    }

    /// Bounds of the seam point's derivative and of the derivative of its
    /// angle with respect to the phase.
    fn seam_rates(&self, cone: &Cone) -> (f64, f64) {
        let profile = &self.law.profile;
        let turn = (profile.slope / cone.slope.abs() + 1.) / self.rate().abs();
        (
            profile.slope + profile.maximum * turn + profile.slope / cone.slope.abs(),
            turn,
        )
    }

    /// Points on the seams where the helicoid's normal is perpendicular to
    /// the sight: contour curves end there. Proven Lipschitz isolation.
    fn contour_on_seams(&self, view: P3) -> Result<Vec<(f64, f64)>, GeometryError> {
        let Some(cone) = self.cone else {
            return Ok(Vec::new());
        };
        let profile = &self.law.profile;
        let k = self.rate();
        let rho = view[0].portable_hypot(view[1]);
        let alpha = view[1].portable_atan2(view[0]);
        let (_, turn) = self.seam_rates(&cone);
        let (rmax, slope, curvature) = (profile.maximum, profile.slope, profile.curvature);
        let lipschitz = rho * (slope + rmax * turn + k.abs() * curvature + k.abs() * slope * turn)
            + view[2].abs() * (slope * slope + rmax * curvature)
            + 1e-12;
        let mut output = Vec::new();
        for seam in self.seams() {
            let g = |phi: f64| {
                let (_, _, theta) = seam.at(phi);
                let [r, d] = seam.span.eval(phi);
                let psi = theta - alpha;
                -r * rho * psi.portable_cos() + k * d * rho * psi.portable_sin() + r * d * view[2]
            };
            let mut roots = Vec::new();
            isolate(&g, seam.range, lipschitz, &mut roots)?;
            for phi in roots {
                let ([_, theta, z], _, _) = seam.at(phi);
                output.push((z, theta));
            }
        }
        Ok(output)
    }

    /// Contour generators of the entry cone, where it bounds the material.
    fn cone_contours(&self, output: &mut Vec<Candidate>) -> Result<(), GeometryError> {
        let Some(cone) = self.cone else {
            return Ok(());
        };
        let view = self.local_vector([0., 0., 1.]);
        let rho = view[0].portable_hypot(view[1]);
        if rho <= 1e-9 {
            return Ok(());
        }
        let ratio = cone.slope * view[2] / rho;
        if ratio.abs() >= 1. {
            return Ok(());
        }
        let alpha = view[1].portable_atan2(view[0]);
        let profile = &self.law.profile;
        let k = self.rate();
        for sign in [-1., 1.] {
            let theta = alpha + sign * ratio.portable_acos();
            let (sin, cos) = theta.portable_sin_cos();
            // The cone bounds the material where it trims the thread: roots
            // of R(z - k theta) - C(z) along the generator.
            let f = |z: f64| profile.radius(z - k * theta) - cone.radius(z);
            let span = [self.z[0], cone.length.min(self.z[1])];
            let mut stops = vec![span[0], span[1]];
            isolate(&f, span, profile.slope + cone.slope.abs(), &mut stops)?;
            stops.sort_by(f64::total_cmp);
            stops.dedup_by(|x, y| (*x - *y).abs() <= 1e-13);
            for pair in stops.windows(2) {
                let mid = 0.5 * (pair[0] + pair[1]);
                if pair[1] - pair[0] > 1e-13
                    && cone.face(cone.radius(mid), profile.radius(mid - k * theta))
                {
                    let at = |z: f64| {
                        let r = cone.radius(z);
                        self.law.view([r * cos, r * sin, z])
                    };
                    output.push(Candidate {
                        kind: Kind::Outline,
                        path: Arc::new(Segment(at(pair[0]), at(pair[1]))),
                    });
                }
            }
        }
        Ok(())
    }

    fn outline_copies(
        &self,
        span: &Span,
        sign: f64,
        range: [f64; 2],
        view: P3,
        ends: &[(f64, f64)],
        output: &mut Vec<Candidate>,
    ) {
        let law = &self.law;
        let k = law.rate;
        let angle = |phi: f64| Outline::angle(law, span, view, sign, phi);
        // theta lies within 3pi/2 of alpha + pi, which bounds z = phi + k theta.
        let alpha = view[1].portable_atan2(view[0]);
        let reach = k.abs() * 1.5 * PI;
        let lo = range[0] + k * (alpha + PI) - reach;
        let hi = range[1] + k * (alpha + PI) + reach;
        let first = ((self.z[0] - hi) / self.pitch()).floor() as i64 - 1;
        let last = ((self.z[1] - lo) / self.pitch()).ceil() as i64 + 1;
        let middle = 0.5 * (range[0] + range[1]);
        for m in first..=last {
            let copy = m as f64 * self.pitch();
            let z = |phi: f64| angle(phi).map(|theta| phi + k * theta + copy);
            // Ends on a cap: contour points of the cap curves on this branch.
            let mut stops = vec![range[0], range[1]];
            for &(level, theta) in ends {
                let phi = level - k * theta;
                let phi = phi - self.pitch() * ((phi - middle) / self.pitch()).round();
                if phi <= range[0] || phi >= range[1] {
                    continue;
                }
                let Some(own) = angle(phi) else { continue };
                let turn = (own - theta).rem_euclid(TAU);
                let same = turn.min(TAU - turn) < 1e-7;
                if same && z(phi).is_some_and(|value| (value - level).abs() < 1e-7) {
                    stops.push(phi);
                }
            }
            stops.sort_by(f64::total_cmp);
            stops.dedup_by(|x, y| (*x - *y).abs() <= 1e-13);
            for pair in stops.windows(2) {
                let mid = 0.5 * (pair[0] + pair[1]);
                let radius = span.eval(mid)[0];
                let inside = z(mid).is_some_and(|value| {
                    value > self.z[0]
                        && value < self.z[1]
                        && self
                            .cone
                            .is_none_or(|cone| cone.keeps(radius, value, 1e-12))
                });
                if !inside || pair[1] - pair[0] <= 1e-13 {
                    continue;
                }
                let width = range[1] - range[0];
                output.push(Candidate {
                    kind: Kind::Outline,
                    path: Arc::new(Range {
                        path: Arc::new(Outline {
                            law: law.clone(),
                            span: *span,
                            sign,
                            copy,
                            range,
                            view,
                        }),
                        range: [(pair[0] - range[0]) / width, (pair[1] - range[0]) / width],
                    }),
                });
            }
        }
    }

    /// Generators of a cylindrical band where the radius is perpendicular to
    /// the sight: a constant angle, and phases from the band crossing it.
    fn flat_outlines(&self, span: &Span, view: P3, output: &mut Vec<Candidate>) {
        let r = span.flat().expect("flat band");
        let alpha = view[1].portable_atan2(view[0]);
        let rate = self.rate();
        let [a, b] = span.range();
        for theta in [alpha + PI / 2., alpha - PI / 2.] {
            // Points at this angle: z = phi + rate * theta + m * pitch.
            for m in self.copies() {
                let shift = rate * theta + m as f64 * self.pitch();
                let [z0, z1] = [a + shift, b + shift];
                let low = self.start(r);
                let [lo, hi] = [z0.max(low), z1.min(self.z[1])];
                if hi - lo > 1e-12 {
                    let (sin, cos) = theta.portable_sin_cos();
                    output.push(Candidate {
                        kind: Kind::Outline,
                        path: Arc::new(Segment(
                            self.law.view([r * cos, r * sin, lo]),
                            self.law.view([r * cos, r * sin, hi]),
                        )),
                    });
                }
            }
        }
    }

    /// Section of the threaded part by the local plane `n . p = c`.
    fn thread_section(&self, n: P3, c: f64, behind: f64) -> Result<Loop, GeometryError> {
        let law = &self.law;
        let profile = &law.profile;
        let rate = self.rate();
        let nxy = n[0].portable_hypot(n[1]);
        let mut output: Loop = Vec::new();
        if nxy <= 1e-12 {
            // Across the axis: the cap curve at that level.
            let level = c / n[2];
            let probe = level - behind * n[2].signum() * 4. * self.tolerance;
            if probe > self.z[0] && probe < self.z[1] {
                output.extend(self.cap(level));
            }
            return Ok(output);
        }
        if n[2].abs() <= 1e-12 && c.abs() <= 4. * self.tolerance {
            // Containing the axis: both meridian profiles and the cap lines.
            let h = [-n[1] / nxy, n[0] / nxy, 0.];
            let theta = h[1].portable_atan2(h[0]);
            for (side, angle) in [(1., theta), (-1., theta + PI)] {
                for m in self.copies() {
                    let shift = rate * angle + m as f64 * self.pitch();
                    for span in &profile.spans {
                        let [a, b] = span.range().map(|phi| phi + shift);
                        let [lo, hi] = [a.max(self.z[0]), b.min(self.z[1])];
                        if hi - lo <= 1e-12 {
                            continue;
                        }
                        let local = span.shifted(m as f64 * self.pitch() + rate * angle);
                        // On an entry chamfer only the part inside the cone
                        // bounds the material.
                        let mut stops = vec![lo, hi];
                        if let Some(cone) = &self.cone {
                            stops.extend(
                                cone_crossings(&local, cone)
                                    .into_iter()
                                    .filter(|z| *z > lo && *z < hi),
                            );
                            // The cone ends there even where rounding puts
                            // its crossing with the band just past its end.
                            if cone.length > lo && cone.length < hi {
                                stops.push(cone.length);
                            }
                        }
                        stops.sort_by(f64::total_cmp);
                        for pair in stops.windows(2) {
                            let [lo, hi] = [pair[0], pair[1]];
                            let mid = 0.5 * (lo + hi);
                            if hi - lo <= 1e-12
                                || self
                                    .cone
                                    .is_some_and(|cone| !cone.keeps(local.eval(mid)[0], mid, 1e-12))
                            {
                                continue;
                            }
                            let point = |z: f64| {
                                let r = local.eval(z)[0];
                                law.view([side * r * h[0], side * r * h[1], z])
                            };
                            let path: Arc<dyn Path> = match local {
                                Span::Line { .. } => Arc::new(Segment(point(lo), point(hi))),
                                Span::Arc {
                                    center,
                                    radius,
                                    sign,
                                    ..
                                } => {
                                    // A circular arc in the plane: exact conic.
                                    let lift = Lift {
                                        origin: law.view([0.; 3]),
                                        x: self.view_vector([0., 0., 1.]),
                                        y: self.view_vector([side * h[0], side * h[1], 0.]),
                                    };
                                    let start = (lo - center[0]) / radius;
                                    let end = (hi - center[0]) / radius;
                                    let angle_of = |x: f64| {
                                        let x = x.clamp(-1., 1.);
                                        sign * (1. - x * x).max(0.).sqrt()
                                    };
                                    let a = [
                                        center[0] + radius * start,
                                        center[1] + radius * angle_of(start),
                                    ];
                                    let b = [
                                        center[0] + radius * end,
                                        center[1] + radius * angle_of(end),
                                    ];
                                    Arc::new(paths::ProfilePath {
                                        edge: Edge::Planar(normal_sweep::Edge::Arc {
                                            center,
                                            a,
                                            b,
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
                // The entry cone's generator where it bounds the material.
                if let Some(cone) = &self.cone {
                    let top = cone.length.min(self.z[1]);
                    let f = |z: f64| profile.radius(z - rate * angle) - cone.radius(z);
                    let mut stops = vec![self.z[0], top];
                    let lipschitz = profile.slope + cone.slope.abs();
                    isolate(&f, [self.z[0], top], lipschitz, &mut stops)?;
                    stops.sort_by(f64::total_cmp);
                    for pair in stops.windows(2) {
                        let mid = 0.5 * (pair[0] + pair[1]);
                        let face = cone.face(cone.radius(mid), profile.radius(mid - rate * angle));
                        if pair[1] - pair[0] > 1e-13 && face {
                            let at = |z: f64| {
                                let r = side * cone.radius(z);
                                law.view([r * h[0], r * h[1], z])
                            };
                            output.push(Arc::new(Segment(at(pair[0]), at(pair[1]))));
                        }
                    }
                }
            }
            for level in self.z {
                let trim = self.cone.filter(|cone| level < cone.length);
                let reach = |angle: f64| {
                    let r = profile.radius(level - rate * angle);
                    match trim {
                        Some(cone) if cone.internal => r.max(cone.radius(level)),
                        Some(cone) => r.min(cone.radius(level)),
                        None => r,
                    }
                };
                let at = |s: f64| law.view([s * h[0], s * h[1], level]);
                let [right, left] = [reach(theta), -reach(theta + PI)];
                match &self.stock {
                    // A bore's cap is the annulus between thread and stock;
                    // at a terminal groove's shoulder, the tail cylinder.
                    // A partial-depth bore's floor closes the void's chord.
                    Some(_) if self.is_floor(level) => {
                        output.push(Arc::new(Segment(at(left), at(right))))
                    }
                    Some(_) => {
                        let [outer, back] = [theta, theta + PI].map(|a| self.stock_reach(level, a));
                        output.push(Arc::new(Segment(at(right), at(outer))));
                        output.push(Arc::new(Segment(at(-back), at(left))));
                    }
                    None => output.push(Arc::new(Segment(at(left), at(right)))),
                }
            }
            return Ok(output);
        }
        // A plane parallel to the axis meets a cylindrical band of radius
        // R along generators at the angles where `n_xy . (cos, sin) R = c`;
        // the sweep below has no extent there, so publish them directly.
        if n[2].abs() <= 1e-12 {
            for span in &profile.spans {
                let Some(radius) = span.flat() else { continue };
                let ratio = c / (radius * nxy);
                if ratio.abs() >= 1. {
                    continue;
                }
                let [a, b] = span.range();
                for sign in [-1., 1.] {
                    let theta = n[1].portable_atan2(n[0]) + sign * ratio.portable_acos();
                    let (sin, cos) = theta.portable_sin_cos();
                    for m in self.copies() {
                        let shift = rate * theta + m as f64 * self.pitch();
                        // On an entry chamfer the band starts at its station.
                        let low = self.start(radius);
                        let [lo, hi] = [(a + shift).max(low), (b + shift).min(self.z[1])];
                        if hi - lo > 1e-12 {
                            output.push(Arc::new(Segment(
                                law.view([radius * cos, radius * sin, lo]),
                                law.view([radius * cos, radius * sin, hi]),
                            )));
                        }
                    }
                }
            }
        }
        // Chords of the caps; their ends are where section branches meet a cap.
        let mut chords = Vec::new();
        for level in self.z {
            let mut level_chords = self.cap_chord(n, c, level)?;
            if let Some(cone) = self.cone.filter(|cone| level < cone.length)
                && let Some(disc) = self.disc_chord(n, c, level, cone.radius(level))
            {
                // On an entry chamfer the cap is also within the cone; on a
                // treated bore the void's cap also includes the cone's disc.
                let d = [n[1] / nxy, -n[0] / nxy, 0.];
                let offset = (c - n[2] * level) / (nxy * nxy);
                let foot = [n[0] * offset, n[1] * offset, level];
                let base = add(foot, scaled(d, -dot3(foot, d)));
                let mut stations: Vec<[f64; 2]> = level_chords
                    .iter()
                    .map(|[a, b]| [dot3(*a, d), dot3(*b, d)])
                    .collect();
                stations = if cone.internal {
                    stations.push(disc);
                    merge(stations, 0.)
                } else {
                    profile::intersect(&stations, &[disc])
                };
                level_chords = stations
                    .into_iter()
                    .map(|[a, b]| [add(base, scaled(d, a)), add(base, scaled(d, b))])
                    .collect();
            } else if self
                .cone
                .is_some_and(|cone| !cone.internal && level < cone.length)
            {
                level_chords.clear();
            }
            chords.extend(level_chords);
        }
        let mut chord_ends: Vec<P3> = chords.iter().flatten().copied().collect();
        // Where section branches meet the entry cone: the plane's crossings
        // of the seams, isolated with a proven Lipschitz bound.
        let mut seam_points = Vec::new();
        if let Some(cone) = &self.cone {
            let (lipschitz, _) = self.seam_rates(cone);
            for seam in self.seams() {
                let h = |phi: f64| {
                    let ([r, theta, z], _, _) = seam.at(phi);
                    dot3(n, [r * theta.portable_cos(), r * theta.portable_sin(), z]) - c
                };
                let mut roots = Vec::new();
                isolate(&h, seam.range, lipschitz + 1e-12, &mut roots)?;
                for phi in roots {
                    let ([r, theta, z], _, _) = seam.at(phi);
                    seam_points.push([r * theta.portable_cos(), r * theta.portable_sin(), z]);
                }
            }
            chord_ends.extend(seam_points.iter().copied());
            // The cone's own section between its ends and the seams.
            let mut events = vec![0., TAU];
            let base = n[1].portable_atan2(n[0]);
            for q in [
                c / cone.entry,
                (c - n[2] * cone.length) / cone.radius(cone.length),
                -n[2] / cone.slope,
            ] {
                let ratio = q / nxy;
                if ratio.abs() < 1. {
                    for sign in [-1., 1.] {
                        events.push((base + sign * ratio.portable_acos()).rem_euclid(TAU));
                    }
                }
            }
            events.extend(
                seam_points
                    .iter()
                    .map(|p| p[1].portable_atan2(p[0]).rem_euclid(TAU)),
            );
            events.sort_by(f64::total_cmp);
            events.dedup_by(|x, y| (*x - *y).abs() <= 1e-13);
            for pair in events.windows(2) {
                let mid = 0.5 * (pair[0] + pair[1]);
                let Some(p) = ConeCut::local(cone, (n, c), mid) else {
                    continue;
                };
                let inside = p[2] > self.z[0]
                    && p[2] < cone.length.min(self.z[1])
                    && cone.face(cone.radius(p[2]), profile.radius(p[2] - rate * mid));
                if inside && pair[1] - pair[0] > 1e-13 {
                    output.push(Arc::new(ConeCut {
                        law: law.clone(),
                        cone: *cone,
                        plane: (n, c),
                        range: [pair[0], pair[1]],
                    }));
                }
            }
        }
        if self.stock.is_some() {
            // A bore's caps are annuli: the stock chord less the thread's,
            // at a terminal groove's shoulder the tail cylinder's.
            let d = [n[1] / nxy, -n[0] / nxy, 0.];
            let mut annular = Vec::new();
            for level in self.z {
                let offset = (c - n[2] * level) / (nxy * nxy);
                let foot = [n[0] * offset, n[1] * offset, level];
                let base = add(foot, scaled(d, -dot3(foot, d)));
                let inner: Vec<[f64; 2]> = chords
                    .iter()
                    .filter(|[a, _]| (a[2] - level).abs() <= 1e-12)
                    .map(|[a, b]| [dot3(*a, d), dot3(*b, d)])
                    .collect();
                if self.is_floor(level) {
                    // The floor bounds only the void's cap.
                    for [a, b] in merge(inner, 0.) {
                        annular.push([add(base, scaled(d, a)), add(base, scaled(d, b))]);
                    }
                    continue;
                }
                let outer: Vec<[f64; 2]> = self.stock_chord(n, c, level).into_iter().collect();
                for [a, b] in profile::exclusive(&outer, &inner, 1e-12) {
                    annular.push([add(base, scaled(d, a)), add(base, scaled(d, b))]);
                }
            }
            chords = annular;
        }
        // General plane: sweep the meridian half-plane once around the axis,
        // in two charts of half a turn each centred where the meridian
        // half-plane lies in the plane's direction. There the plane meets the
        // meridian lines at the most grazing angle (near the axis a profile
        // point moves by the angle's error times the radius squared over the
        // offset), so the angle is measured from the chart centre at full
        // precision.
        let base = n[0].portable_atan2(-n[1]);
        let [low, high] = [*self.copies().start() - 1, *self.copies().end() + 1];
        for base in [base, base + PI] {
            let (sin, cos) = base.portable_sin_cos();
            let n = [n[0] * cos + n[1] * sin, n[1] * cos - n[0] * sin, n[2]];
            let chart = [-PI / 2., PI / 2.];
            for m in low..=high {
                for span in &profile.spans {
                    if n[2].abs() <= 1e-12 && span.flat().is_some() {
                        continue;
                    }
                    let span = span.shifted(m as f64 * self.pitch() + rate * base);
                    let mut events = chart.to_vec();
                    let mut touches = Vec::new();
                    let [a, b] = span.range();
                    for e in [span.point(a), span.point(b)] {
                        let f = |theta: f64| {
                            let [alpha, beta, gamma] = Sweep::line(law, (n, c), theta);
                            alpha * e[1] + beta * e[0] - gamma
                        };
                        let lip = nxy * e[1].abs() + (n[2] * rate).abs() + 1e-12;
                        isolate(&f, chart, lip, &mut events)?;
                    }
                    if let Span::Arc { center, radius, .. } = span {
                        // The meridian line touches the circle where its signed
                        // distance from the center, `u / |(alpha, beta)|`, is
                        // `+-R`: the roots of `u -+ R |(alpha, beta)|`, isolated
                        // one sign at a time so that a pair of close touches (near
                        // the axis, where `|(alpha, beta)|` is small) never needs
                        // resolving as one near-double root.
                        let lu = nxy * center[1].abs() + (n[2] * rate).abs();
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
                            // The branch leaves the material only through a cap,
                            // at an end of that cap's exact chord.
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
                                // Where a cap is crossed, exactly: near the
                                // axis the chord end's angle alone places the
                                // branch point poorly along a span that the
                                // meridian lines meet at a grazing angle.
                                let t = if self.z.contains(&point[2]) {
                                    let z = |t: f64| law.local(branch.point(t))[2];
                                    level_stop(&z, t, point[2])
                                } else {
                                    t
                                };
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
                                let q = law.local(sweep.point(0.5 * (piece[0] + piece[1])));
                                let z = q[2];
                                let kept = self.cone.is_none_or(|cone| {
                                    cone.keeps(q[0].portable_hypot(q[1]), z, 1e-12)
                                });
                                if kept
                                    && z > self.z[0]
                                    && z < self.z[1]
                                    && piece[1] - piece[0] > 1e-13
                                {
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

    /// The chords of a cap disc cut by the plane, as local end points.
    /// Parameters of a unit-direction local line lying in the plane of
    /// `level` inside the thread's cross-section there.
    fn cap_line(&self, s: P3, d: P3, level: f64) -> Result<Vec<[f64; 2]>, GeometryError> {
        let profile = &self.law.profile;
        let rho = d[0].portable_hypot(d[1]);
        if rho <= 1e-15 {
            let r = s[0].portable_hypot(s[1]);
            let inside = r < profile.radius(level - self.rate() * s[1].portable_atan2(s[0]));
            return Ok(if inside {
                vec![[-1e9, 1e9]]
            } else {
                Vec::new()
            });
        }
        let g = |t: f64| {
            let p = add(s, scaled(d, t));
            profile.radius(level - self.rate() * p[1].portable_atan2(p[0]))
                - p[0].portable_hypot(p[1])
        };
        let b = s[0] * d[0] + s[1] * d[1];
        let foot2 = s[0] * s[0] + s[1] * s[1] - b * b;
        let reach = profile.maximum;
        if foot2 >= reach * reach {
            return Ok(Vec::new());
        }
        let half = (reach * reach - foot2).sqrt();
        let mut stops = vec![-b - half, -b + half];
        let parts: Vec<[f64; 2]> = if foot2 < profile.minimum * profile.minimum {
            let inner = (profile.minimum * profile.minimum - foot2).sqrt();
            vec![[-b - half, -b - inner], [-b + inner, -b + half]]
        } else {
            vec![[-b - half, -b + half]]
        };
        for part in parts {
            stops.extend(part);
            stops.extend(self.crossings(s, d, part)?);
        }
        stops.sort_by(f64::total_cmp);
        stops.dedup_by(|x, y| (*x - *y).abs() <= 1e-13);
        Ok(merge(
            stops
                .windows(2)
                .filter(|pair| pair[1] - pair[0] > 1e-13 && g(0.5 * (pair[0] + pair[1])) > 0.)
                .map(|pair| [pair[0], pair[1]])
                .collect(),
            0.,
        ))
    }

    fn cap_chord(&self, n: P3, c: f64, level: f64) -> Result<Vec<[P3; 2]>, GeometryError> {
        let nxy = n[0].portable_hypot(n[1]);
        if nxy <= 1e-12 {
            return Ok(Vec::new());
        }
        let d = [n[1] / nxy, -n[0] / nxy, 0.];
        let offset = (c - n[2] * level) / (nxy * nxy);
        let s = [n[0] * offset, n[1] * offset, level];
        // In the cap plane the line is in the disc where r < R(level - k theta).
        let profile = &self.law.profile;
        let g = |t: f64| {
            let p = add(s, scaled(d, t));
            profile.radius(level - self.rate() * p[1].portable_atan2(p[0]))
                - p[0].portable_hypot(p[1])
        };
        let reach = profile.maximum;
        let foot = s[0].portable_hypot(s[1]);
        if foot >= reach {
            return Ok(Vec::new());
        }
        let half = (reach * reach - foot * foot).sqrt();
        let core = (foot < profile.minimum).then(|| (profile.minimum.powi(2) - foot * foot).sqrt());
        let mut stops = vec![-half, half];
        let parts: Vec<[f64; 2]> = match core {
            Some(inner) => vec![[-half, -inner], [inner, half]],
            None => vec![[-half, half]],
        };
        for part in parts {
            // The cap lies at one level, so the line's phase is the angle's.
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
}

impl Member for Thread {
    fn bounds(&self) -> [P3; 2] {
        let own = self.own_bounds();
        match &self.base {
            Some(prism) => {
                let other = prism.bounds();
                [
                    std::array::from_fn(|k| own[0][k].min(other[0][k])),
                    std::array::from_fn(|k| own[1][k].max(other[1][k])),
                ]
            }
            None => own,
        }
    }

    fn line(&self, start: P3, direction: P3) -> Option<Result<Vec<[f64; 2]>, GeometryError>> {
        self.own_line(start, direction)
    }

    fn candidates(&self, tolerance_mm: f64) -> Result<Vec<Candidate>, GeometryError> {
        let mut output = self.own_candidates(tolerance_mm)?;
        if let Some(prism) = &self.base {
            output.extend(prism.candidates(tolerance_mm)?);
        }
        Ok(output)
    }

    fn sight(&self, [x, y]: P2) -> Result<Vec<[f64; 2]>, GeometryError> {
        match self.own_line([x, y, 0.], [0., 0., 1.]) {
            Some(result) => result,
            None => unreachable!("threads answer lines"),
        }
    }

    fn section(&self, depth: f64, tolerance_mm: f64) -> Result<Vec<Loop>, GeometryError> {
        let mut loops = self.own_section(depth, tolerance_mm)?;
        let Some(prism) = &self.base else {
            return Ok(loops);
        };
        let n = self.local_vector([0., 0., 1.]);
        let c = depth - self.law.origin[2];
        let nxy = n[0].portable_hypot(n[1]);
        let mut output: Loop = loops.pop().unwrap_or_default();
        let base: Loop = prism
            .section(depth, tolerance_mm)?
            .into_iter()
            .flatten()
            .collect();
        if nxy <= 1e-12 {
            // Across the axis only one part is cut: the one behind the plane.
            let level = c / n[2];
            let probe = level - n[2].signum() * 4. * self.tolerance;
            if probe < self.z[0] {
                output = base;
            }
        } else {
            // At the junction both parts publish chords; the boundary there
            // is the base chord less the thread's.
            let junction = self.z[0];
            let at_junction = |path: &Arc<dyn Path>| {
                let [a, b] = [0., 1.].map(|t| self.law.local(path.point(t))[2]);
                (a - junction).abs() <= 4. * self.tolerance
                    && (b - junction).abs() <= 4. * self.tolerance
                    && matches!(path.paper([0., 1.]), Some(Curve::Line(_)))
            };
            let d = [n[1] / nxy, -n[0] / nxy, 0.];
            let station = |path: &Arc<dyn Path>| -> [f64; 2] {
                let [a, b] = [0., 1.].map(|t| dot3(self.law.local(path.point(t)), d));
                [a.min(b), a.max(b)]
            };
            let thread: Vec<[f64; 2]> = merge(
                output
                    .iter()
                    .filter(|p| at_junction(p))
                    .map(station)
                    .collect(),
                0.,
            );
            let lower: Vec<[f64; 2]> = merge(
                base.iter()
                    .filter(|p| at_junction(p))
                    .map(station)
                    .collect(),
                0.,
            );
            output.retain(|path| !at_junction(path));
            output.extend(base.into_iter().filter(|path| !at_junction(path)));
            let offset = (c - n[2] * junction) / (nxy * nxy);
            let foot = [n[0] * offset, n[1] * offset, junction];
            let origin = add(foot, scaled(d, -dot3(foot, d)));
            for [a, b] in profile::exclusive(&lower, &thread, 1e-12) {
                output.push(Arc::new(Segment(
                    self.law.view(add(origin, scaled(d, a))),
                    self.law.view(add(origin, scaled(d, b))),
                )));
            }
        }
        Ok(if output.is_empty() {
            vec![]
        } else {
            vec![output]
        })
    }
}

impl Thread {
    fn own_bounds(&self) -> [P3; 2] {
        if let Some(bounds) = self.stock_bounds() {
            return bounds;
        }
        let reach = self
            .stock
            .as_ref()
            .map_or(self.law.profile.maximum, |(radius, _)| *radius);
        let top = match (&self.neck, &self.groove) {
            (Some((_, span, _)), _) => span[1],
            (None, Some(groove)) => groove.top,
            (None, None) => self.z[1],
        };
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for z in [self.z[0], top] {
            let center = self.law.view([0., 0., z]);
            for k in 0..3 {
                let axis = self.law.axes[2][k];
                let extent = reach * (1. - axis * axis).max(0.).sqrt();
                lo[k] = lo[k].min(center[k] - extent);
                hi[k] = hi[k].max(center[k] + extent);
            }
        }
        [lo, hi]
    }

    fn own_candidates(&self, tolerance_mm: f64) -> Result<Vec<Candidate>, GeometryError> {
        let mut output = Vec::new();
        self.rails(&mut output);
        for level in self.z {
            output.extend(self.cap(level).into_iter().map(|path| Candidate {
                kind: Kind::Edge,
                path,
            }));
        }
        self.outlines(&mut output)?;
        self.cone_contours(&mut output)?;
        for seam in self.seams() {
            output.push(Candidate {
                kind: Kind::Edge,
                path: Arc::new(seam),
            });
        }
        if let Some((_, _, prism)) = &self.neck {
            output.extend(prism.candidates(tolerance_mm)?);
        }
        if let Some((_, prism)) = &self.stock {
            output.extend(prism.candidates(tolerance_mm)?);
        }
        Ok(output)
    }

    fn own_section(&self, depth: f64, tolerance_mm: f64) -> Result<Vec<Loop>, GeometryError> {
        let n = self.local_vector([0., 0., 1.]);
        let c = depth - self.law.origin[2];
        let mut output = self.thread_section(n, c, 1.)?;
        if let Some((_, prism)) = &self.stock {
            // The stock wall; its cap chords are replaced by the annuli.
            let cap = |path: &Arc<dyn Path>| {
                let [a, b] = [0., 1.].map(|t| self.law.local(path.point(t))[2]);
                self.z.iter().any(|level| {
                    (a - level).abs() <= 4. * self.tolerance
                        && (b - level).abs() <= 4. * self.tolerance
                }) && matches!(path.paper([0., 1.]), Some(Curve::Line(_)))
            };
            let axial = n[0].portable_hypot(n[1]) <= 1e-12;
            output.extend(
                prism
                    .section(depth, tolerance_mm)?
                    .into_iter()
                    .flatten()
                    .filter(|path| axial || !cap(path)),
            );
        }
        if let Some((_, span, prism)) = &self.neck {
            // The shoulder: the thread's chord at the junction less the neck's.
            let junction = span[0];
            let is_junction = |path: &Arc<dyn Path>| {
                let [a, b] = [0., 1.].map(|t| self.law.local(path.point(t))[2]);
                (a - junction).abs() <= 4. * self.tolerance
                    && (b - junction).abs() <= 4. * self.tolerance
            };
            output.retain(|path| {
                !(is_junction(path) && matches!(path.paper([0., 1.]), Some(Curve::Line(_))))
            });
            let neck: Loop = prism
                .section(depth, tolerance_mm)?
                .into_iter()
                .flatten()
                .collect();
            let nxy = n[0].portable_hypot(n[1]);
            let axial = nxy <= 1e-12;
            if axial {
                // Across the axis, only one part is cut at a level.
                let level = c / n[2];
                let probe = level - n[2].signum() * 4. * self.tolerance;
                if probe > junction {
                    output = neck;
                }
            } else {
                output.extend(neck.into_iter().filter(|path| !is_junction(path)));
                let thread = self.cap_chord(n, c, junction)?;
                let neck_chord = self.disc_chord(n, c, junction, self.neck.as_ref().unwrap().0);
                let d = [n[1] / nxy, -n[0] / nxy, 0.];
                let station = |p: P3| dot3(p, d);
                let foot = {
                    let offset = (c - n[2] * junction) / (nxy * nxy);
                    [n[0] * offset, n[1] * offset, junction]
                };
                let outer: Vec<[f64; 2]> = thread
                    .iter()
                    .map(|[a, b]| [station(*a), station(*b)])
                    .collect();
                let inner: Vec<[f64; 2]> = neck_chord.into_iter().collect();
                let base = add(foot, scaled(d, -station(foot)));
                for [a, b] in profile::exclusive(&outer, &inner, 1e-12) {
                    output.push(Arc::new(Segment(
                        self.law.view(add(base, scaled(d, a))),
                        self.law.view(add(base, scaled(d, b))),
                    )));
                }
            }
        }
        Ok(if output.is_empty() {
            vec![]
        } else {
            vec![output]
        })
    }

    fn own_line(&self, start: P3, direction: P3) -> Option<Result<Vec<[f64; 2]>, GeometryError>> {
        let s = self.law.local(start);
        let d = self.local_vector(direction);
        let length = dot3(d, d).sqrt();
        if length <= 1e-300 {
            return Some(Ok(Vec::new()));
        }
        let unit = scaled(d, 1. / length);
        Some(self.line_local(s, unit).map(|spans| {
            spans
                .into_iter()
                .map(|[a, b]| [a / length, b / length])
                .collect()
        }))
    }
}

impl Thread {
    /// The outer radius of a bore's cap annulus at `level`: the stock's, or
    /// at a terminal groove's station the tail cylinder's.
    fn cap_reach(&self, level: f64) -> f64 {
        match (&self.groove, &self.stock) {
            (Some(groove), _) if (level - groove.station).abs() <= self.tolerance => groove.radius,
            (_, Some((radius, _))) => *radius,
            _ => self.law.profile.maximum,
        }
    }

    /// The station interval, along the cap chord direction, of a disc of
    /// radius `radius` at `level` cut by the plane.
    fn disc_chord(&self, n: P3, c: f64, level: f64, radius: f64) -> Option<[f64; 2]> {
        let nxy = n[0].portable_hypot(n[1]);
        let d = [n[1] / nxy, -n[0] / nxy, 0.];
        let offset = (c - n[2] * level) / (nxy * nxy);
        let s = [n[0] * offset, n[1] * offset, level];
        let foot = s[0].portable_hypot(s[1]);
        (foot < radius).then(|| {
            let half = (radius * radius - foot * foot).sqrt();
            let center = dot3(s, d);
            [center - half, center + half]
        })
    }
}
