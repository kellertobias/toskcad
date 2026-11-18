//! A perpendicular cross-bore removed from one round band.
//!
//! The band stays an exact cylinder; the void is the finite tool cylinder
//! whose axis crosses the band axis at a right angle. A line of sight is
//! blocked when the band material along it, less the open tool interval,
//! keeps positive length in front of the point. Visibility along a candidate
//! curve can only change where its line of sight is tangent to the tool,
//! where the curve crosses the tool wall or a bottom plane, where the line
//! of sight crosses a bottom circle, or where its paper image crosses the
//! paper image of a bored-wall rim. Every such parameter is an event; each
//! interval between events is classified once at its midpoint.
//!
//! Rim crossings are the roots of `g(θ) = σ A w(θ) + H(θ)` with
//! `w = sqrt(r² - (e + ρ cos θ)²)` along the rim angle θ. Squaring removes
//! the branch sign: `A² w² - H²` is a trigonometric polynomial, solved in
//! two bounded half-angle charts by the shared polynomial root finder.
//! Squared roots of the other branch and repeated-root candidates only add
//! harmless extra events.
use super::*;
use crate::curved_projection::{derivative, roots};
use crate::drawing_events::Harmonic;
use std::f64::consts::{FRAC_PI_2, PI, TAU};

type Poly = Vec<f64>;

fn sub(a: Point, b: Point) -> Point {
    std::array::from_fn(|i| a[i] - b[i])
}

fn dot(a: Point, b: Point) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn add(a: &[f64], b: &[f64]) -> Poly {
    (0..a.len().max(b.len()))
        .map(|i| a.get(i).copied().unwrap_or(0.) + b.get(i).copied().unwrap_or(0.))
        .collect()
}

fn mul(a: &[f64], b: &[f64]) -> Poly {
    let mut product = vec![0.; a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            product[i + j] = x.mul_add(*y, product[i + j]);
        }
    }
    product
}

fn scale(a: &[f64], factor: f64) -> Poly {
    a.iter().map(|x| x * factor).collect()
}

/// Numerator over `1 + τ²` of `h0 + h1 cos θ + h2 sin θ` in the chart
/// `θ = phase + 2 atan τ`, where `sign` is `cos phase` for phase 0 or π.
fn harmonic([h0, h1, h2]: [f64; 3], sign: f64) -> Poly {
    vec![h0 + sign * h1, 2. * sign * h2, h0 - sign * h1]
}

/// Every root angle in `[0, 2π)` of the numerator built for each chart,
/// plus repeated-root candidates and the chart boundaries.
fn trig_roots(build: impl Fn(f64) -> Poly, out: &mut Vec<f64>) {
    for (phase, sign) in [(0., 1.), (PI, -1.)] {
        let polynomial = build(sign);
        let largest = polynomial.iter().map(|c| c.abs()).fold(0., f64::max);
        if !largest.is_finite() || largest == 0. {
            continue;
        }
        let polynomial = scale(&polynomial, 1. / largest);
        let chart = [-1., 1.];
        out.extend(
            roots(&polynomial, chart)
                .into_iter()
                .chain(roots(&derivative(&polynomial), chart))
                .map(|t| (phase + 2. * t.atan()).rem_euclid(TAU)),
        );
        out.push(phase + FRAC_PI_2);
    }
}

/// A finite tool cylinder removed from a band, perpendicular to its axis.
#[derive(Clone, Copy, Debug)]
pub struct CrossVoid {
    origin: Point,
    /// Unit tool travel `x'`, offset direction `y'` and band axis `z`.
    axes: [Point; 3],
    wall: f64,
    offset: f64,
    level: f64,
    radius: f64,
    stations: [f64; 2],
    through: [bool; 2],
}

impl CrossVoid {
    /// `origin` is the band-axis point of the bore chart. All values are in
    /// drawing coordinates; the bored wall radius `wall` must hold the tool
    /// section strictly, and at least one side must pass the wall.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        origin: Point,
        axes: [Point; 3],
        wall: f64,
        offset: f64,
        level: f64,
        radius: f64,
        stations: [f64; 2],
        through: [bool; 2],
    ) -> Result<Self, ProjectionError> {
        let orthonormal = (0..3).all(|i| {
            (0..3).all(|j| {
                let expected = if i == j { 1. } else { 0. };
                (dot(axes[i], axes[j]) - expected).abs() <= 1e-9
            })
        });
        if !bounded(origin)
            || !orthonormal
            || ![wall, offset, level, radius, stations[0], stations[1]]
                .into_iter()
                .all(|v| v.is_finite() && v.abs() <= 1e8)
            || radius <= 0.
            || offset.abs() + radius >= wall
            || stations[0] >= stations[1]
            || through == [false; 2]
        {
            return Err(ProjectionError::Invalid);
        }
        Ok(Self {
            origin,
            axes,
            wall,
            offset,
            level,
            radius,
            stations,
            through,
        })
    }

    fn local(&self, point: Point) -> Point {
        let delta = sub(point, self.origin);
        self.axes.map(|axis| dot(axis, delta))
    }

    /// Components of the direction toward the viewer in the bore chart.
    fn view(&self) -> Point {
        self.axes.map(|axis| axis[2])
    }

    /// The bored-wall rim on a through `side` at the tool angle.
    fn rim(&self, side: usize, angle: f64) -> Point {
        let (sin, cos) = angle.sin_cos();
        let y = self.offset + self.radius * cos;
        let w = ((self.wall - y) * (self.wall + y)).max(0.).sqrt();
        let local = [
            if side == 0 { -w } else { w },
            y,
            self.level + self.radius * sin,
        ];
        std::array::from_fn(|i| {
            self.origin[i] + (0..3).map(|k| self.axes[k][i] * local[k]).sum::<f64>()
        })
    }

    /// Open depth interval along the ray from `start` toward the viewer in
    /// which the ray is strictly inside the tool.
    pub(super) fn interval(&self, start: Point) -> Option<Span> {
        let [x, y, z] = self.local(start);
        let [dx, dy, dz] = self.view();
        let [y, z] = [y - self.offset, z - self.level];
        let a = dy.mul_add(dy, dz * dz);
        let b = 2. * y.mul_add(dy, z * dz);
        let c = y.mul_add(y, z * z) - self.radius * self.radius;
        let mut span = if a <= 1e-24 {
            if c >= 0. {
                return None;
            }
            [f64::NEG_INFINITY, f64::INFINITY]
        } else {
            let disc = b.mul_add(b, -4. * a * c);
            if disc <= 0. {
                return None;
            }
            let root = disc.sqrt();
            [(-b - root) / (2. * a), (-b + root) / (2. * a)]
        };
        if dx == 0. {
            if x <= self.stations[0] || x >= self.stations[1] {
                return None;
            }
        } else {
            let [p, q] = self.stations.map(|s| (s - x) / dx);
            span = [span[0].max(p.min(q)), span[1].min(p.max(q))];
        }
        (span[0] < span[1]).then_some(span)
    }

    /// The parts of a line outside the open tool.
    pub(super) fn trim_line(&self, [a, b]: [Point; 2]) -> Vec<[Point; 2]> {
        let [la, lb] = [self.local(a), self.local(b)];
        let d = sub(lb, la);
        let [y, z] = [la[1] - self.offset, la[2] - self.level];
        let qa = d[1].mul_add(d[1], d[2] * d[2]);
        let qb = 2. * y.mul_add(d[1], z * d[2]);
        let qc = y.mul_add(y, z * z) - self.radius * self.radius;
        let mut inside = if qa <= 1e-24 {
            if qc >= 0. {
                return vec![[a, b]];
            }
            [0., 1.]
        } else {
            let disc = qb.mul_add(qb, -4. * qa * qc);
            if disc <= 0. {
                return vec![[a, b]];
            }
            let root = disc.sqrt();
            [(-qb - root) / (2. * qa), (-qb + root) / (2. * qa)]
        };
        if d[0] == 0. {
            if la[0] <= self.stations[0] || la[0] >= self.stations[1] {
                return vec![[a, b]];
            }
        } else {
            let [p, q] = self.stations.map(|s| (s - la[0]) / d[0]);
            inside = [inside[0].max(p.min(q)), inside[1].min(p.max(q))];
        }
        let [lo, hi] = [inside[0].max(0.), inside[1].min(1.)];
        if lo >= hi {
            return vec![[a, b]];
        }
        let at = |t: f64| -> Point { std::array::from_fn(|i| (b[i] - a[i]).mul_add(t, a[i])) };
        [[0., lo], [hi, 1.]]
            .into_iter()
            .filter(|[from, to]| to > from)
            .map(|[from, to]| [at(from), at(to)])
            .collect()
    }

    /// Rim angles at which the paper image of a bored-wall rim meets the
    /// paper line through `base` normal to `normal`.
    fn line_crossings(&self, base: [f64; 2], normal: [f64; 2], out: &mut Vec<f64>) {
        let paper = |v: Point| normal[0].mul_add(v[0], normal[1] * v[1]);
        let a = paper(self.axes[0]);
        let b = self.radius * paper(self.axes[1]);
        let c = self.radius * paper(self.axes[2]);
        let center: Point = std::array::from_fn(|i| {
            self.origin[i] + self.offset * self.axes[1][i] + self.level * self.axes[2][i]
        });
        let k = paper(center) - normal[0].mul_add(base[0], normal[1] * base[1]);
        let (e, rho, r) = (self.offset, self.radius, self.wall);
        trig_roots(
            |sign| {
                let one = vec![1., 0., 1.];
                let wall = add(
                    &scale(&mul(&one, &one), r * r),
                    &scale(
                        &mul(&harmonic([e, rho, 0.], sign), &harmonic([e, rho, 0.], sign)),
                        -1.,
                    ),
                );
                let h = harmonic([k, b, c], sign);
                add(&scale(&wall, a * a), &scale(&mul(&h, &h), -1.))
            },
            out,
        );
    }

    fn sides(&self) -> impl Iterator<Item = usize> + '_ {
        (0..2).filter(|&side| self.through[side])
    }

    /// Events along the straight edge `[a, b]`, already lifted by the depth
    /// guard, as parameters in `(0, 1)`.
    pub(super) fn line_events(&self, [a, b]: [Point; 2], cuts: &mut Vec<f64>) {
        let [la, lb] = [self.local(a), self.local(b)];
        let [x, y, z] = std::array::from_fn(|i| Affine::new(la[i], lb[i]));
        let (y, z) = (y.offset(-self.offset), z.offset(-self.level));
        let [v0, v1, v2] = self.view();
        let mut push = |t: f64| {
            if 0. < t && t < 1. {
                cuts.push(t);
            }
        };
        let mut level = |f: Affine, value: f64| {
            if f.delta != 0. {
                push((value - f.start) / f.delta);
            }
        };
        let tangent = v1.hypot(v2);
        if tangent > 1e-12 {
            let f = y.combine(-v2 / tangent, z, v1 / tangent);
            level(f, self.radius);
            level(f, -self.radius);
        }
        for s in self.stations {
            level(x, s);
        }
        let mut circle = |p: Affine, q: Affine, radius: f64| {
            let mut events = Events::new();
            events.circle(p, q, radius);
            let sorted = events.sorted();
            for &t in &sorted[1..sorted.len() - 1] {
                push(t);
            }
        };
        circle(y, z, self.radius);
        if v0.abs() > 1e-12 {
            for s in self.stations {
                circle(
                    y.combine(v0, x, -v1).offset(s * v1),
                    z.combine(v0, x, -v2).offset(s * v2),
                    self.radius * v0.abs(),
                );
            }
        }
        let start = [a[0], a[1]];
        let delta = [b[0] - a[0], b[1] - a[1]];
        let length = delta[0].mul_add(delta[0], delta[1] * delta[1]);
        if length == 0. {
            return;
        }
        let mut angles = Vec::new();
        self.line_crossings(start, [-delta[1], delta[0]], &mut angles);
        for angle in angles {
            for side in self.sides() {
                let p = self.rim(side, angle);
                push(((p[0] - start[0]) * delta[0] + (p[1] - start[1]) * delta[1]) / length);
            }
        }
    }

    /// Events along a harmonic rim `origin + cos t * cosine + sin t * sine`,
    /// already lifted by the depth guard, within `trim`.
    pub(super) fn rim_events(
        &self,
        [origin, cosine, sine]: [Point; 3],
        trim: Span,
        cuts: &mut Vec<f64>,
    ) {
        let start = sub(origin, self.origin);
        let [x, y, z] = std::array::from_fn(|k| {
            Harmonic([
                dot(self.axes[k], start),
                dot(self.axes[k], cosine),
                dot(self.axes[k], sine),
            ])
        });
        let (y, z) = (y.offset(-self.offset), z.offset(-self.level));
        let [v0, v1, v2] = self.view();
        let tangent = v1.hypot(v2);
        if tangent > 1e-12 {
            let f = y.combine(-v2 / tangent, z, v1 / tangent);
            f.level(self.radius, trim, cuts);
            f.level(-self.radius, trim, cuts);
        }
        for s in self.stations {
            x.level(s, trim, cuts);
        }
        crate::drawing_events::circle(y, z, self.radius, trim, cuts);
        if v0.abs() > 1e-12 {
            for s in self.stations {
                crate::drawing_events::circle(
                    y.combine(v0, x, -v1).offset(s * v1),
                    z.combine(v0, x, -v2).offset(s * v2),
                    self.radius * v0.abs(),
                    trim,
                    cuts,
                );
            }
        }
        // Crossings of the rim's paper image with each bored-wall rim.
        let o = [origin[0], origin[1]];
        let [a, b] = [[cosine[0], cosine[1]], [sine[0], sine[1]]];
        let det = a[0] * b[1] - a[1] * b[0];
        let size = a[0].hypot(a[1]).max(b[0].hypot(b[1]));
        let mut angles = Vec::new();
        let mut push = |angle: f64| {
            if angle > trim[0] && angle < trim[1] {
                cuts.push(angle);
            }
        };
        if size == 0. {
            return;
        }
        if det.abs() <= 1e-9 * size * size {
            // An edge-on rim is a paper segment along its larger axis.
            let dir = if a[0].hypot(a[1]) >= b[0].hypot(b[1]) {
                a
            } else {
                b
            };
            self.line_crossings(o, [-dir[1], dir[0]], &mut angles);
            for angle in angles {
                for side in self.sides() {
                    let p = self.rim(side, angle);
                    let along = (p[0] - o[0]) * dir[0] + (p[1] - o[1]) * dir[1];
                    let mut found = Vec::new();
                    Harmonic([
                        0.,
                        a[0] * dir[0] + a[1] * dir[1],
                        b[0] * dir[0] + b[1] * dir[1],
                    ])
                    .level(along, trim, &mut found);
                    found.into_iter().for_each(&mut push);
                }
            }
            return;
        }
        // q = M⁻¹ (p - o) with M = [a b]; the rim is |q| = 1.
        let inverse = |v: [f64; 2]| {
            [
                (b[1] * v[0] - b[0] * v[1]) / det,
                (a[0] * v[1] - a[1] * v[0]) / det,
            ]
        };
        let paper = |v: Point| [v[0], v[1]];
        let center: Point = std::array::from_fn(|i| {
            self.origin[i] + self.offset * self.axes[1][i] + self.level * self.axes[2][i]
        });
        let q0 = inverse([center[0] - o[0], center[1] - o[1]]);
        let q1 = inverse(paper(self.axes[0]));
        let q2 = inverse(paper(self.axes[1]).map(|c| c * self.radius));
        let q3 = inverse(paper(self.axes[2]).map(|c| c * self.radius));
        let (e, rho, r) = (self.offset, self.radius, self.wall);
        let width = q1[0].mul_add(q1[0], q1[1] * q1[1]);
        trig_roots(
            |sign| {
                let one = vec![1., 0., 1.];
                let one2 = mul(&one, &one);
                let hw = harmonic([e, rho, 0.], sign);
                let wall = add(&scale(&one2, r * r), &scale(&mul(&hw, &hw), -1.));
                let qx = harmonic([q0[0], q2[0], q3[0]], sign);
                let qy = harmonic([q0[1], q2[1], q3[1]], sign);
                let u = add(
                    &add(&mul(&qx, &qx), &mul(&qy, &qy)),
                    &add(&scale(&wall, width), &scale(&one2, -1.)),
                );
                let v = harmonic(
                    [
                        2. * (q1[0] * q0[0] + q1[1] * q0[1]),
                        2. * (q1[0] * q2[0] + q1[1] * q2[1]),
                        2. * (q1[0] * q3[0] + q1[1] * q3[1]),
                    ],
                    sign,
                );
                add(&mul(&u, &u), &scale(&mul(&wall, &mul(&v, &v)), -1.))
            },
            &mut angles,
        );
        for angle in angles {
            for side in self.sides() {
                let p = self.rim(side, angle);
                let [c, s] = inverse([p[0] - o[0], p[1] - o[1]]);
                push(s.atan2(c).rem_euclid(TAU));
            }
        }
    }
}

impl RoundOccluder {
    /// Remove a perpendicular cross-bore from this band.
    pub fn with_cross_void(mut self, void: CrossVoid) -> Self {
        self.cross = Some(void);
        self
    }

    /// Whether band material less the tool remains along the ray from
    /// `start` toward the viewer, up to `maximum_depth`.
    pub(super) fn cross_ray_hits(
        &self,
        start: Point,
        maximum_depth: f64,
        void: &CrossVoid,
    ) -> bool {
        let w0 = sub(start, self.base);
        let axis = self.axis;
        let a0 = dot(w0, axis);
        let ad = axis[2];
        let mut span = [0., maximum_depth];
        if ad == 0. {
            if !(a0 > 0. && a0 < self.height) {
                return false;
            }
        } else {
            let [p, q] = [-a0 / ad, (self.height - a0) / ad];
            span = [span[0].max(p.min(q)), span[1].min(p.max(q))];
        }
        let wp: Point = std::array::from_fn(|i| w0[i] - a0 * axis[i]);
        let dp: Point = std::array::from_fn(|i| f64::from(u8::from(i == 2)) - ad * axis[i]);
        let qa = dot(dp, dp);
        let qb = 2. * dot(wp, dp);
        let qc = dot(wp, wp) - self.radius * self.radius;
        if qa <= 1e-24 {
            if qc >= 0. {
                return false;
            }
        } else {
            let disc = qb.mul_add(qb, -4. * qa * qc);
            if disc <= 0. {
                return false;
            }
            let root = disc.sqrt();
            span = [
                span[0].max((-qb - root) / (2. * qa)),
                span[1].min((-qb + root) / (2. * qa)),
            ];
        }
        if span[1] <= span[0] {
            return false;
        }
        match void.interval(start) {
            None => true,
            Some([lo, hi]) => lo > span[0] || hi < span[1],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A band along +x through the origin, radius 9, x ∈ [-15, 15], bored
    /// radially along the view axis by a Ø5.5 tool.
    fn bored(axes: [Point; 3]) -> RoundOccluder {
        let band = RoundOccluder::new([-15., 0., 0.], axes[2], 9., 30., None).unwrap();
        band.with_cross_void(
            CrossVoid::new([0.; 3], axes, 9., 0., 0., 2.75, [-50., 50.], [true; 2]).unwrap(),
        )
    }

    #[test]
    fn a_line_behind_the_bore_is_seen_through_it() {
        // Tool travel toward the viewer; band axis along x.
        let round = bored([[0., 0., 1.], [0., 1., 0.], [1., 0., 0.]]);
        let hidden = round
            .hidden_intervals([[-20., 0., -30.], [20., 0., -30.]], 1e-6)
            .unwrap();
        // Hidden behind the band except through the hole |x| < 2.75.
        assert_eq!(hidden.len(), 2);
        let [first, second] = [hidden[0], hidden[1]];
        assert!((first[0] - 5. / 40.).abs() < 1e-9);
        assert!((first[1] - (20. - 2.75) / 40.).abs() < 1e-9);
        assert!((second[0] - (20. + 2.75) / 40.).abs() < 1e-9);
        assert!((second[1] - 35. / 40.).abs() < 1e-9);
        assert!(!round.blocks_point([0., 0., -30.], f64::INFINITY, 1e-6));
        assert!(round.blocks_point([5., 0., -30.], f64::INFINITY, 1e-6));
    }

    #[test]
    fn a_blind_bore_and_an_oblique_view_still_hide_correctly() {
        let axes = [[0., 0., 1.], [0., 1., 0.], [1., 0., 0.]];
        let band = RoundOccluder::new([-15., 0., 0.], axes[2], 9., 30., None).unwrap();
        // Blind from the viewer side to 1 mm past the axis: still opaque.
        let blind = band.with_cross_void(
            CrossVoid::new([0.; 3], axes, 9., 0., 0., 2.75, [-1., 50.], [false, true]).unwrap(),
        );
        assert!(blind.blocks_point([0., 0., -30.], f64::INFINITY, 1e-6));
        assert!(!blind.blocks_point([0., 0., 0.], f64::INFINITY, 1e-6));
        // Tilted view: the straight-through window shrinks and moves, and
        // the events still split exactly where visibility changes.
        let (s, c) = 0.3_f64.sin_cos();
        let tilted = [[0., -s, c], [0., c, s], [1., 0., 0.]];
        let round = bored(tilted);
        let line = [[-20., 0., -30.], [20., 0., -30.]];
        let hidden = round.hidden_intervals(line, 1e-6).unwrap();
        let at = |t: f64| std::array::from_fn(|i| (line[1][i] - line[0][i]).mul_add(t, line[0][i]));
        for t in (1..400).map(|k| f64::from(k) / 400.) {
            let inside = hidden.iter().any(|[a, b]| *a < t && t < *b);
            let near = hidden
                .iter()
                .any(|[a, b]| (t - a).abs() < 1e-6 || (t - b).abs() < 1e-6);
            if !near {
                assert_eq!(
                    inside,
                    round.blocks_point(at(t), f64::INFINITY, 1e-6),
                    "{t}"
                );
            }
        }
    }

    #[test]
    fn rim_events_split_where_the_bore_reveals_a_rim() {
        let (s, c) = 0.25_f64.sin_cos();
        let round = bored([[0., -s, c], [0., c, s], [1., 0., 0.]]);
        // A tilted circle behind the band, crossing the bore window.
        let rim = Rim {
            ellipse: curved_projection::Ellipse {
                center: [0.1, 0.05],
                u: [1., 0.],
                v: [0., 1.],
                major: 4.,
                minor: 0.4,
                trim: [0., TAU],
            },
            depth: -30.,
            depth_cosine: 2.,
            depth_sine: -1.,
        };
        let mut cuts = Vec::new();
        let mut hidden = Vec::new();
        round.hide_rim(rim, 1e-6, &mut cuts, false, |span| hidden.push(span));
        assert!(!hidden.is_empty());
        let mut revealed = 0;
        for k in 1..720 {
            let t = TAU * f64::from(k) / 720.;
            let near = hidden
                .iter()
                .any(|[a, b]| (t - a).abs() < 1e-6 || (t - b).abs() < 1e-6);
            if near {
                continue;
            }
            let inside = hidden.iter().any(|[a, b]| *a < t && t < *b);
            let blocked = round.blocks_point(rim.point(t), f64::INFINITY, 1e-6);
            assert_eq!(inside, blocked, "{t}");
            revealed += usize::from(!blocked && rim.point(t)[0].abs() < 15.);
        }
        assert!(revealed > 0, "part of the rim is seen through the bore");
    }

    #[test]
    fn invalid_voids_are_refused_and_silhouettes_are_trimmed() {
        let axes = [[0., 0., 1.], [0., 1., 0.], [1., 0., 0.]];
        assert!(CrossVoid::new([0.; 3], axes, 9., 7., 0., 2.75, [-1., 1.], [true; 2]).is_err());
        assert!(CrossVoid::new([0.; 3], axes, 9., 0., 0., 2.75, [1., -1.], [true; 2]).is_err());
        assert!(CrossVoid::new([0.; 3], axes, 9., 0., 0., 2.75, [-1., 1.], [false; 2]).is_err());
        let void = CrossVoid::new([0.; 3], axes, 9., 0., 0., 2.75, [-50., 50.], [true; 2]).unwrap();
        // The generator on the near wall loses the hole's chord.
        let pieces = void.trim_line([[-15., 0., 9.], [15., 0., 9.]]);
        assert_eq!(pieces.len(), 2);
        assert!((pieces[0][1][0] + 2.75).abs() < 1e-12 && (pieces[1][0][0] - 2.75).abs() < 1e-12);
        assert_eq!(void.trim_line([[-15., 5., 9.], [15., 5., 9.]]).len(), 1);
    }
}
