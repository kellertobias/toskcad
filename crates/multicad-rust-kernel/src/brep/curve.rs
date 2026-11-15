//! Edge support curves in their exchange parameterisation.
use super::surface::Frame;
use super::{
    Error, NurbsCurve, Point, Result, add, distance, dot, ensure, nearest, scale, sub, turn,
};
use std::f64::consts::{PI, TAU};

/// The widest angular panel used for sampling and quadrature.
pub(crate) const PANEL: f64 = PI / 8.;

#[derive(Clone, Debug, PartialEq)]
pub enum Curve {
    /// `origin + t * direction`; the direction keeps the stated magnitude.
    Line {
        origin: Point,
        direction: Point,
    },
    /// `origin + radius * (cos t * x + sin t * y)`.
    Circle {
        frame: Frame,
        radius: f64,
    },
    /// `origin + a * cos t * x + b * sin t * y`.
    Ellipse {
        frame: Frame,
        semiaxes: [f64; 2],
    },
    Nurbs {
        curve: NurbsCurve,
        closed: bool,
    },
}

impl Curve {
    pub fn nurbs(curve: NurbsCurve, tolerance: f64) -> Self {
        let [a, b] = curve.range();
        let closed = distance(curve.evaluate(a).0, curve.evaluate(b).0) <= tolerance;
        Self::Nurbs { curve, closed }
    }

    pub fn period(&self) -> Option<f64> {
        match self {
            Self::Line { .. } => None,
            Self::Circle { .. } | Self::Ellipse { .. } => Some(TAU),
            Self::Nurbs { curve, closed } => {
                let [a, b] = curve.range();
                closed.then_some(b - a)
            }
        }
    }

    fn wrap(&self, t: f64) -> f64 {
        match self {
            Self::Nurbs {
                curve,
                closed: true,
            } => {
                let [a, b] = curve.range();
                a + (t - a).rem_euclid(b - a)
            }
            _ => t,
        }
    }

    /// Point and first derivative.
    pub fn evaluate(&self, t: f64) -> (Point, Point) {
        match self {
            Self::Line { origin, direction } => (add(*origin, scale(*direction, t)), *direction),
            Self::Circle { frame, radius } => ellipse(frame, [*radius; 2], t),
            Self::Ellipse { frame, semiaxes } => ellipse(frame, *semiaxes, t),
            Self::Nurbs { curve, .. } => curve.evaluate(self.wrap(t)),
        }
    }

    /// Parameter of the curve point closest to `point`. With a seed the
    /// result is the periodic representative nearest to it.
    pub fn invert(&self, point: Point, seed: Option<f64>) -> f64 {
        let t = match self {
            Self::Line { origin, direction } => {
                dot(sub(point, *origin), *direction) / dot(*direction, *direction)
            }
            Self::Circle { frame, .. } => {
                let local = frame.local(point);
                local[1].atan2(local[0]).rem_euclid(TAU)
            }
            Self::Ellipse { frame, semiaxes } => {
                let local = frame.local(point);
                (local[1] / semiaxes[1])
                    .atan2(local[0] / semiaxes[0])
                    .rem_euclid(TAU)
            }
            Self::Nurbs { curve, .. } => {
                let [a, b] = curve.range();
                let start = seed.map_or_else(
                    || Self::nearest_sample(curve, point),
                    |seed| self.wrap(seed),
                );
                let mut t = start;
                for _ in 0..40 {
                    let (p, d) = curve.evaluate(t);
                    let step = dot(sub(point, p), d) / dot(d, d).max(f64::MIN_POSITIVE);
                    let next = if self.period().is_some() {
                        self.wrap(t + step)
                    } else {
                        (t + step).clamp(a, b)
                    };
                    let done = (next - t).abs() <= 1e-13 * (b - a);
                    t = next;
                    if done {
                        break;
                    }
                }
                t
            }
        };
        seed.map_or(t, |seed| nearest(t, seed, self.period()))
    }

    /// The knot-panel sample of `curve` nearest to `point`.
    fn nearest_sample(curve: &NurbsCurve, point: Point) -> f64 {
        let [a, b] = curve.range();
        let cuts = curve.breaks(None, a, b);
        let mut best = (f64::INFINITY, a);
        for pair in cuts.windows(2) {
            for k in 0..=4 {
                let t = pair[0] + (pair[1] - pair[0]) * f64::from(k) / 4.;
                let d = distance(curve.evaluate(t).0, point);
                if d < best.0 {
                    best = (d, t);
                }
            }
        }
        best.1
    }

    /// Panel boundaries from `a` to `b`: knot spans for splines and bounded
    /// angular steps for conics.
    pub(crate) fn breaks(&self, a: f64, b: f64) -> Vec<f64> {
        match self {
            Self::Line { .. } => vec![a, b],
            Self::Circle { .. } | Self::Ellipse { .. } => uniform(a, b, PANEL),
            Self::Nurbs { curve, .. } => curve.breaks(self.period(), a, b),
        }
    }

    /// Parameters from `a` to `b` whose chords stay within `deflection` and
    /// whose tangent turns by at most `angle` along each chord.
    pub(crate) fn sample(&self, a: f64, b: f64, deflection: f64, angle: f64) -> Result<Vec<f64>> {
        let mut parameters = vec![a];
        for pair in self.breaks(a, b).windows(2) {
            let mut pieces = 1usize;
            while !self.flat(pair[0], pair[1], pieces, deflection, angle) {
                pieces *= 2;
                ensure(pieces <= 4096, Error::Limit)?;
            }
            for k in 1..=pieces {
                parameters.push(pair[0] + (pair[1] - pair[0]) * k as f64 / pieces as f64);
            }
        }
        Ok(parameters)
    }

    fn flat(&self, a: f64, b: f64, pieces: usize, deflection: f64, angle: f64) -> bool {
        (0..pieces).all(|k| {
            let at = |s: f64| self.evaluate(a + (b - a) * (k as f64 + s) / pieces as f64);
            let (start, middle, end) = (at(0.), at(0.5), at(1.));
            let chord = scale(add(start.0, end.0), 0.5);
            distance(chord, middle.0) <= deflection && turn(start.1, end.1) <= angle
        })
    }

    /// Largest and smallest coordinate along `axis` between `a` and `b`.
    pub(crate) fn extent(&self, a: f64, b: f64, axis: usize) -> [f64; 2] {
        let mut range = [f64::INFINITY, f64::NEG_INFINITY];
        let mut include = |t: f64| {
            let x = self.evaluate(t).0[axis];
            range = [range[0].min(x), range[1].max(x)];
        };
        include(a);
        include(b);
        for pair in self.breaks(a, b).windows(2) {
            // Each panel is narrow enough that the coordinate has at most
            // one interior extremum, bracketed by its derivative's sign.
            let slope = |t: f64| self.evaluate(t).1[axis];
            let (mut low, mut high) = (pair[0], pair[1]);
            include(high);
            if slope(low) * slope(high) < 0. {
                for _ in 0..80 {
                    let middle = 0.5 * (low + high);
                    if slope(low) * slope(middle) <= 0. {
                        high = middle;
                    } else {
                        low = middle;
                    }
                }
                include(0.5 * (low + high));
            }
        }
        range
    }

    pub(crate) fn mapped(
        &self,
        map: &dyn Fn(Point) -> Point,
        linear: &dyn Fn(Point) -> Point,
    ) -> Self {
        match self {
            Self::Line { origin, direction } => Self::Line {
                origin: map(*origin),
                direction: linear(*direction),
            },
            Self::Circle { frame, radius } => Self::Circle {
                frame: frame.mapped(map, linear),
                radius: *radius,
            },
            Self::Ellipse { frame, semiaxes } => Self::Ellipse {
                frame: frame.mapped(map, linear),
                semiaxes: *semiaxes,
            },
            Self::Nurbs { curve, closed } => Self::Nurbs {
                curve: curve.mapped(map),
                closed: *closed,
            },
        }
    }
}

fn ellipse(frame: &Frame, [a, b]: [f64; 2], t: f64) -> (Point, Point) {
    let (sin, cos) = t.sin_cos();
    (
        add(
            frame.origin,
            add(scale(frame.x, a * cos), scale(frame.y, b * sin)),
        ),
        add(scale(frame.x, -a * sin), scale(frame.y, b * cos)),
    )
}

/// `a..=b` cut into equal steps no wider than `step`.
pub(crate) fn uniform(a: f64, b: f64, step: f64) -> Vec<f64> {
    let pieces = ((b - a).abs() / step).ceil().max(1.) as usize;
    (0..=pieces)
        .map(|k| a + (b - a) * k as f64 / pieces as f64)
        .collect()
}
