//! Straight single- or multi-start screw threads: material inside the
//! helicoid `r <= R(z - k * theta)` between two flat caps, with an optional
//! coaxial cylindrical neck beyond the far cap (a terminal thread relief).
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
//! curvature bounds. Contour and section branch paths still flatten
//! adaptively: near a branch end their parametrisation has a square-root
//! singularity, so no uniform derivative bound exists.
use super::paths::{Lift, Segment, add, dot3, scaled};
use super::prism::{Band, Prism};
use super::profile::merge;
use super::*;
use multicad_rust_kernel::{helical_thread, normal_sweep, rounded_thread, spur_gear::Edge};
use std::f64::consts::{PI, TAU};

/// One piece of the axial meridian profile in `(phase, radius)` coordinates.
#[derive(Clone, Copy, Debug)]
enum Span {
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
    fn range(&self) -> [f64; 2] {
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
    fn eval(&self, phi: f64) -> [f64; 2] {
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

    fn point(&self, phi: f64) -> P2 {
        [phi, self.eval(phi)[0]]
    }

    /// Phases inside the span where the radius is stationary.
    fn stationary(&self) -> Option<f64> {
        match *self {
            Self::Line { a, b } => (a[1] == b[1]).then_some(0.5 * (a[0] + b[0])),
            Self::Arc { center, range, .. } => {
                (range[0] <= center[0] && center[0] <= range[1]).then_some(center[0])
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
    fn meet(&self, alpha: f64, beta: f64, gamma: f64) -> Vec<(usize, P2)> {
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
    fn arc_root(&self, alpha: f64, beta: f64, gamma: f64, index: usize, clamp: bool) -> Option<P2> {
        let Self::Arc {
            center,
            radius,
            sign,
            range,
        } = *self
        else {
            return None;
        };
        let norm = alpha.hypot(beta);
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
                radius: (a - center[0]).hypot(r(a) - center[1]),
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
const RESOLVED: f64 = 1e-10;

/// Roots of `f` on `[a, b]` with `|f'| <= lipschitz`: a stretch is excluded
/// when its end values prove it free of zeros; sign changes are bisected.
fn isolate(
    f: &dyn Fn(f64) -> f64,
    [a, b]: [f64; 2],
    lipschitz: f64,
    output: &mut Vec<f64>,
) -> Result<(), GeometryError> {
    let width = 1e-13 * a.abs().max(b.abs()).max(1.);
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
        if fa.abs() + fb.abs() > lipschitz * (b - a) {
            continue;
        }
        // Where the function cannot leave a band of `RESOLVED` around zero,
        // the line stays that close to the surface: keep only sign changes.
        if b - a <= width || lipschitz * (b - a) <= RESOLVED {
            if (fa <= 0.) != (fb <= 0.) {
                output.push(0.5 * (a + b));
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

/// Adaptive certified-in-practice flattening of an exact path.
fn adaptive(point: &dyn Fn(f64) -> P3, tolerance: f64, start: usize) -> Vec<f64> {
    let deviation = |a: f64, b: f64| {
        let [p, q] = [point(a), point(b)];
        [0.25, 0.5, 0.75]
            .into_iter()
            .map(|s| {
                let m = point(a + (b - a) * s);
                let c = add(p, scaled(add(q, scaled(p, -1.)), s));
                let d = add(m, scaled(c, -1.));
                dot3(d, d).sqrt()
            })
            .fold(0., f64::max)
    };
    let mut output = vec![0.];
    let count = start.max(2);
    for i in 0..count {
        let mut stack = vec![(i as f64 / count as f64, (i + 1) as f64 / count as f64, 0)];
        while let Some((a, b, depth)) = stack.pop() {
            if depth < 24 && deviation(a, b) > tolerance {
                let m = 0.5 * (a + b);
                stack.push((m, b, depth + 1));
                stack.push((a, m, depth + 1));
            } else {
                output.push(b);
            }
        }
    }
    *output.last_mut().unwrap() = 1.;
    output
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
        let (sin, cos) = theta.sin_cos();
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
        if self.law.axes[2][0].hypot(self.law.axes[2][1]) > 1e-12 {
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
        let rho = view[0].hypot(view[1]);
        let alpha = view[1].atan2(view[0]);
        let a = -r * rho;
        let b = law.rate * slope * rho;
        let c = -r * slope * view[2];
        let norm = a.hypot(b);
        if norm <= 1e-300 || c.abs() > norm * (1. + slack) {
            return None;
        }
        // atan2(b, a) is continuous as pi + atan2(-b, -a) since a < 0.
        Some(alpha + PI + (-b).atan2(-a) + sign * (c / norm).clamp(-1., 1.).acos())
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
}

impl Path for Outline {
    fn point(&self, t: f64) -> P3 {
        self.at(self.phi(t))
    }

    fn flatten(&self, tolerance_mm: f64) -> Vec<f64> {
        adaptive(&|t| self.point(t), tolerance_mm, 32)
    }

    fn paper(&self, _: [f64; 2]) -> Option<Curve> {
        None
    }
}

/// A range of another path.
struct Range {
    path: Arc<dyn Path>,
    range: [f64; 2],
}

impl Path for Range {
    fn point(&self, t: f64) -> P3 {
        self.path
            .point(self.range[0] + (self.range[1] - self.range[0]) * t)
    }

    fn flatten(&self, tolerance_mm: f64) -> Vec<f64> {
        adaptive(&|t| self.point(t), tolerance_mm, 8)
    }

    fn paper(&self, [a, b]: [f64; 2]) -> Option<Curve> {
        let span = self.range[1] - self.range[0];
        self.path
            .paper([self.range[0] + a * span, self.range[0] + b * span])
    }
}

/// A branch of an oblique section on a shifted span: for each `theta` the
/// `index`-th meet of the plane's meridian line with the span.
struct Sweep {
    law: Arc<Law>,
    span: Span,
    index: usize,
    plane: (P3, f64),
    range: [f64; 2],
    /// Ends at which the plane's meridian line touches the arc: both roots
    /// meet there, so the end is the touching point itself.
    touching: [bool; 2],
}

impl Sweep {
    fn line(law: &Law, (n, c): (P3, f64), theta: f64) -> [f64; 3] {
        let (sin, cos) = theta.sin_cos();
        let alpha = n[0] * cos + n[1] * sin;
        [alpha, n[2], c - n[2] * law.rate * theta]
    }

    fn at(&self, theta: f64, touch: bool) -> P3 {
        let [alpha, beta, gamma] = Self::line(&self.law, self.plane, theta);
        let p = match self.span {
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
        };
        self.law.surface(p[1], theta, p[0] + self.law.rate * theta)
    }
}

impl Path for Sweep {
    fn point(&self, t: f64) -> P3 {
        if t <= 0. {
            return self.at(self.range[0], self.touching[0]);
        }
        if t >= 1. {
            return self.at(self.range[1], self.touching[1]);
        }
        self.at(self.range[0] + (self.range[1] - self.range[0]) * t, false)
    }

    fn flatten(&self, tolerance_mm: f64) -> Vec<f64> {
        let start = ((self.range[1] - self.range[0]).abs() * 32.).ceil() as usize;
        adaptive(&|t| self.point(t), tolerance_mm, start)
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
    tolerance: f64,
}

pub(in crate::drawing) fn member(
    solid: &Solid,
    view: &View,
) -> Option<Result<Box<dyn Member>, GeometryError>> {
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
            Ok(Box::new(Thread {
                law: Arc::new(Law {
                    origin,
                    axes,
                    profile,
                    rate: thread.helical_rate(),
                    starts: thread.spec().starts,
                }),
                z: [0., span[0]],
                neck: Some((neck, span, prism)),
                tolerance,
            }) as Box<dyn Member>)
        })());
    }
    let data = solid.rust_rounded_thread_data()?;
    let Some(profile) = Profile::from_rounded(&data.thread) else {
        return Some(Err(invalid()));
    };
    let (origin, axes) = place(data.frame);
    let spec = data.thread.spec();
    Some(Ok(Box::new(Thread {
        law: Arc::new(Law {
            origin,
            axes,
            profile,
            rate: data.thread.helical_rate(),
            starts: spec.starts,
        }),
        z: [0., spec.length_mm],
        neck: None,
        tolerance,
    })))
}

impl Thread {
    fn rate(&self) -> f64 {
        self.law.rate
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
        let rho = d[0].hypot(d[1]);
        // Radial range: r(t)^2 = |s_xy + t d_xy|^2.
        let within = |radius: f64| -> Option<[f64; 2]> {
            if rho <= 1e-15 {
                return (s[0].hypot(s[1]) < radius).then_some([f64::NEG_INFINITY, f64::INFINITY]);
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
            let r = p[0].hypot(p[1]);
            profile.radius(p[2] - self.rate() * p[1].atan2(p[0])) - r
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
        let rho = d[0].hypot(d[1]);
        if rho <= 1e-15 {
            // Along the axis: a fixed angle and radius, a linear phase.
            return Ok(self.axial_crossings(s, d, [a, b]));
        }
        let g = |t: f64| {
            let p = add(s, scaled(d, t));
            profile.radius(p[2] - rate * p[1].atan2(p[0])) - p[0].hypot(p[1])
        };
        // First derivative, exact on the span containing the point.
        let slope = |t: f64| {
            let p = add(s, scaled(d, t));
            let r2 = p[0] * p[0] + p[1] * p[1];
            let r = r2.sqrt();
            let phi = p[2] - rate * p[1].atan2(p[0]);
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
        let lipschitz = profile.slope * phase_rate + rho + 1e-12;
        let second = profile.curvature * phase_rate * phase_rate
            + profile.slope * rate.abs() * 2. * rho * rho / (profile.minimum * profile.minimum)
            + rho * rho / profile.minimum
            + 1e-12;
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
            if fa.abs() + fb.abs() > lipschitz * h {
                continue;
            }
            // Second order: from either end, within one span of the profile.
            let excluded = [(a, fa), (b, fb)].into_iter().any(|(t, value)| {
                let (derivative, phi, index, base) = slope(t);
                let [lo, hi] = profile.spans[index].range().map(|x| x + base);
                let room = (phi - lo).min(hi - phi);
                phase_rate * h < room && value.abs() > derivative.abs() * h + 0.5 * second * h * h
            });
            if excluded {
                continue;
            }
            if h <= width || lipschitz * h <= RESOLVED {
                if (fa <= 0.) != (fb <= 0.) {
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
        let r = s[0].hypot(s[1]);
        let theta = s[1].atan2(s[0]);
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
            let rho = d[0].hypot(d[1]);
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
        let mut spans = self.threaded(s, d)?;
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
            let ends = self.z.map(|z| (z - base) / rate);
            let [a, b] = [ends[0].min(ends[1]), ends[0].max(ends[1])];
            // Split where the radius turns perpendicular to the sight: there
            // the helix touches the contour of the cylinder through it, and
            // its visibility can change without a crossing in the paper.
            let view = self.local_vector([0., 0., 1.]);
            let mut stops = vec![a, b];
            if view[0].hypot(view[1]) > 1e-12 {
                let fold = view[1].atan2(view[0]) + PI / 2.;
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
        for (span, [a, b]) in self.cap_pieces(level) {
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
        let rho = view[0].hypot(view[1]);
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
        let ends = self.contour_on_caps(view)?;
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
    fn cap_pieces(&self, level: f64) -> Vec<(Span, [f64; 2])> {
        let profile = &self.law.profile;
        let rate = self.rate();
        let mut cuts = vec![0., TAU];
        for m in self.copies() {
            for span in &profile.spans {
                let phi = span.range()[0] + m as f64 * self.pitch();
                let theta = (level - phi) / rate;
                if theta > 0. && theta < TAU {
                    cuts.push(theta);
                }
            }
        }
        cuts.sort_by(f64::total_cmp);
        cuts.dedup_by(|a, b| (*a - *b).abs() <= 1e-14);
        cuts.windows(2)
            .map(|pair| {
                let mid = 0.5 * (pair[0] + pair[1]);
                let (index, base) = profile.locate(level - rate * mid);
                (profile.spans[index].shifted(base), [pair[0], pair[1]])
            })
            .collect()
    }

    /// Angles on each cap curve where the surface normal is perpendicular to
    /// the sight: there contour curves end on the cap. Isolated with a proven
    /// Lipschitz bound of `n . v` along the cap curve.
    fn contour_on_caps(&self, view: P3) -> Result<Vec<(f64, f64)>, GeometryError> {
        let profile = &self.law.profile;
        let k = self.rate();
        let rho = view[0].hypot(view[1]);
        let alpha = view[1].atan2(view[0]);
        let (rmax, slope, curvature) = (profile.maximum, profile.slope, profile.curvature);
        let lipschitz = rho * (k.abs() * slope + rmax)
            + rho * k.abs() * (k.abs() * curvature + slope)
            + view[2].abs() * k.abs() * (slope * slope + rmax * curvature)
            + 1e-12;
        let mut output = Vec::new();
        for level in self.z {
            for (span, range) in self.cap_pieces(level) {
                let g = |theta: f64| {
                    let [r, d] = span.eval(level - k * theta);
                    let psi = theta - alpha;
                    -r * rho * psi.cos() + k * d * rho * psi.sin() + r * d * view[2]
                };
                let mut roots = Vec::new();
                isolate(&g, range, lipschitz, &mut roots)?;
                output.extend(roots.into_iter().map(|theta| (level, theta)));
            }
        }
        Ok(output)
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
        let alpha = view[1].atan2(view[0]);
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
                let inside = z(0.5 * (pair[0] + pair[1]))
                    .is_some_and(|value| value > self.z[0] && value < self.z[1]);
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
        let alpha = view[1].atan2(view[0]);
        let rate = self.rate();
        let [a, b] = span.range();
        for theta in [alpha + PI / 2., alpha - PI / 2.] {
            // Points at this angle: z = phi + rate * theta + m * pitch.
            for m in self.copies() {
                let shift = rate * theta + m as f64 * self.pitch();
                let [z0, z1] = [a + shift, b + shift];
                let [lo, hi] = [z0.max(self.z[0]), z1.min(self.z[1])];
                if hi - lo > 1e-12 {
                    let (sin, cos) = theta.sin_cos();
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
        let nxy = n[0].hypot(n[1]);
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
            let theta = h[1].atan2(h[0]);
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
                                let b =
                                    [center[0] + radius * end, center[1] + radius * angle_of(end)];
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
            for level in self.z {
                let reach = |angle: f64| profile.radius(level - rate * angle);
                let a = law.view([reach(theta) * h[0], reach(theta) * h[1], level]);
                let b = law.view([-reach(theta + PI) * h[0], -reach(theta + PI) * h[1], level]);
                output.push(Arc::new(Segment(b, a)));
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
                    let theta = n[1].atan2(n[0]) + sign * ratio.acos();
                    let (sin, cos) = theta.sin_cos();
                    for m in self.copies() {
                        let shift = rate * theta + m as f64 * self.pitch();
                        let [lo, hi] = [(a + shift).max(self.z[0]), (b + shift).min(self.z[1])];
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
            chords.extend(self.cap_chord(n, c, level)?);
        }
        let chord_ends: Vec<P3> = chords.iter().flatten().copied().collect();
        // General plane: sweep the meridian half-plane once around the axis.
        for m in self.copies() {
            for span in &profile.spans {
                if n[2].abs() <= 1e-12 && span.flat().is_some() {
                    continue;
                }
                let span = span.shifted(m as f64 * self.pitch());
                let mut events = vec![0., TAU];
                let mut touches = Vec::new();
                let [a, b] = span.range();
                for e in [span.point(a), span.point(b)] {
                    let f = |theta: f64| {
                        let [alpha, beta, gamma] = Sweep::line(law, (n, c), theta);
                        alpha * e[1] + beta * e[0] - gamma
                    };
                    let lip = nxy * e[1].abs() + (n[2] * rate).abs() + 1e-12;
                    isolate(&f, [0., TAU], lip, &mut events)?;
                }
                if let Span::Arc { center, radius, .. } = span {
                    let f = |theta: f64| {
                        let [alpha, beta, gamma] = Sweep::line(law, (n, c), theta);
                        let u = alpha * center[1] + beta * center[0] - gamma;
                        u * u - radius * radius * (alpha * alpha + beta * beta)
                    };
                    let lu = nxy * center[1].abs() + (n[2] * rate).abs();
                    let umax = nxy * center[1].abs()
                        + n[2].abs() * center[0].abs()
                        + c.abs()
                        + (n[2] * rate).abs() * TAU;
                    let lip = 2. * umax * lu + 2. * radius * radius * nxy * nxy + 1e-12;
                    isolate(&f, [0., TAU], lip, &mut touches)?;
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
                        let sweep: Arc<dyn Path> = Arc::new(Sweep {
                            law: law.clone(),
                            span,
                            index,
                            plane: (n, c),
                            range: [pair[0], pair[1]],
                            touching: [touch(pair[0]), touch(pair[1])],
                        });
                        // The branch leaves the material only through a cap,
                        // at an end of that cap's exact chord.
                        let mut stops = vec![0., 1.];
                        for point in &chord_ends {
                            let theta = point[1].atan2(point[0]).rem_euclid(TAU);
                            for theta in [theta, theta + TAU, theta - TAU] {
                                if theta <= pair[0] || theta >= pair[1] {
                                    continue;
                                }
                                let t = (theta - pair[0]) / (pair[1] - pair[0]);
                                let q = law.local(sweep.point(t));
                                let d = add(q, scaled(*point, -1.));
                                if dot3(d, d).sqrt() < 1e-6 {
                                    stops.push(t);
                                }
                            }
                        }
                        stops.sort_by(f64::total_cmp);
                        stops.dedup_by(|x, y| (*x - *y).abs() <= 1e-13);
                        for piece in stops.windows(2) {
                            let z = law.local(sweep.point(0.5 * (piece[0] + piece[1])))[2];
                            if z > self.z[0] && z < self.z[1] && piece[1] - piece[0] > 1e-13 {
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
        for [a, b] in chords {
            output.push(Arc::new(Segment(law.view(a), law.view(b))));
        }
        Ok(output)
    }

    /// The chords of a cap disc cut by the plane, as local end points.
    fn cap_chord(&self, n: P3, c: f64, level: f64) -> Result<Vec<[P3; 2]>, GeometryError> {
        let nxy = n[0].hypot(n[1]);
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
            profile.radius(level - self.rate() * p[1].atan2(p[0])) - p[0].hypot(p[1])
        };
        let reach = profile.maximum;
        let foot = s[0].hypot(s[1]);
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
        let reach = self.law.profile.maximum;
        let top = self.neck.as_ref().map_or(self.z[1], |(_, span, _)| span[1]);
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

    fn candidates(&self, tolerance_mm: f64) -> Result<Vec<Candidate>, GeometryError> {
        let mut output = Vec::new();
        self.rails(&mut output);
        for level in self.z {
            output.extend(self.cap(level).into_iter().map(|path| Candidate {
                kind: Kind::Edge,
                path,
            }));
        }
        self.outlines(&mut output)?;
        if let Some((_, _, prism)) = &self.neck {
            output.extend(prism.candidates(tolerance_mm)?);
        }
        Ok(output)
    }

    fn sight(&self, [x, y]: P2) -> Result<Vec<[f64; 2]>, GeometryError> {
        match self.line([x, y, 0.], [0., 0., 1.]) {
            Some(result) => result,
            None => unreachable!("threads answer lines"),
        }
    }

    fn section(&self, depth: f64, tolerance_mm: f64) -> Result<Vec<Loop>, GeometryError> {
        let n = self.local_vector([0., 0., 1.]);
        let c = depth - self.law.origin[2];
        let mut output = self.thread_section(n, c, 1.)?;
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
            let nxy = n[0].hypot(n[1]);
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

    fn line(&self, start: P3, direction: P3) -> Option<Result<Vec<[f64; 2]>, GeometryError>> {
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
    /// The station interval, along the cap chord direction, of a disc of
    /// radius `radius` at `level` cut by the plane.
    fn disc_chord(&self, n: P3, c: f64, level: f64, radius: f64) -> Option<[f64; 2]> {
        let nxy = n[0].hypot(n[1]);
        let d = [n[1] / nxy, -n[0] / nxy, 0.];
        let offset = (c - n[2] * level) / (nxy * nxy);
        let s = [n[0] * offset, n[1] * offset, level];
        let foot = s[0].hypot(s[1]);
        (foot < radius).then(|| {
            let half = (radius * radius - foot * foot).sqrt();
            let center = dot3(s, d);
            [center - half, center + half]
        })
    }
}
