//! The common planar arrangement of two line-and-arc profiles.
//!
//! Every profile curve is split at every point where it meets a curve of
//! either profile, including the ends of collinear or co-circular overlaps,
//! so coincident boundary pieces of both profiles become one shared piece.
//! Each piece then knows, for each profile, whether the material lies on its
//! left, its right, both sides or neither. Tangent contacts, near contacts
//! and junctions that do not separate within the tolerance are refused.
use multicad_precision::PortableMath as _;
use std::collections::BTreeMap;
use std::f64::consts::{FRAC_PI_2, PI, TAU};

pub(super) type P2 = [f64; 2];

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) enum Shape {
    Line {
        a: P2,
        b: P2,
    },
    /// `center + radius * (cos, sin)(start + t * sweep)` for `t` in [0, 1];
    /// a positive sweep runs counter-clockwise.
    Arc {
        center: P2,
        radius: f64,
        start: f64,
        sweep: f64,
    },
}

pub(super) fn sub(a: P2, b: P2) -> P2 {
    [a[0] - b[0], a[1] - b[1]]
}

fn dot(a: P2, b: P2) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

fn cross(a: P2, b: P2) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

pub(super) fn distance(a: P2, b: P2) -> f64 {
    let d = sub(a, b);
    d[0].portable_hypot(d[1])
}

fn polar(center: P2, radius: f64, angle: f64) -> P2 {
    let (sin, cos) = angle.portable_sin_cos();
    [center[0] + radius * cos, center[1] + radius * sin]
}

fn angle_of(center: P2, point: P2) -> f64 {
    let d = sub(point, center);
    d[1].portable_atan2(d[0])
}

impl Shape {
    pub(super) fn point(self, t: f64) -> P2 {
        match self {
            Self::Line { a, b } => [a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])],
            Self::Arc {
                center,
                radius,
                start,
                sweep,
            } => polar(center, radius, start + t * sweep),
        }
    }

    pub(super) fn ends(self) -> [P2; 2] {
        match self {
            Self::Line { a, b } => [a, b],
            _ => [self.point(0.), self.point(1.)],
        }
    }

    pub(super) fn length(self) -> f64 {
        match self {
            Self::Line { a, b } => distance(a, b),
            Self::Arc { radius, sweep, .. } => radius * sweep.abs(),
        }
    }

    /// The piece of this curve between two parameters.
    fn between(self, t0: f64, t1: f64, ends: [P2; 2]) -> Self {
        match self {
            Self::Line { .. } => Self::Line {
                a: ends[0],
                b: ends[1],
            },
            Self::Arc {
                center,
                radius,
                start,
                sweep,
            } => Self::Arc {
                center,
                radius,
                start: start + t0 * sweep,
                sweep: (t1 - t0) * sweep,
            },
        }
    }

    /// The parameters of `point` on this curve, if it lies within
    /// `tolerance` of the finite curve.
    fn parameter(self, point: P2, tolerance: f64) -> Option<f64> {
        let [a, b] = self.ends();
        if distance(point, a) <= tolerance {
            return Some(0.);
        }
        if distance(point, b) <= tolerance {
            return Some(1.);
        }
        match self {
            Self::Line { a, b } => {
                let d = sub(b, a);
                let t = dot(sub(point, a), d) / dot(d, d);
                ((0. ..=1.).contains(&t) && distance(self.point(t), point) <= tolerance)
                    .then_some(t)
            }
            Self::Arc {
                center,
                radius,
                start,
                sweep,
            } => {
                if (distance(point, center) - radius).abs() > tolerance {
                    return None;
                }
                let offset = (angle_of(center, point) - start) * sweep.signum();
                let t = offset.rem_euclid(TAU) / sweep.abs();
                (t <= 1.).then_some(t)
            }
        }
    }

    /// The nearest distance from `point` to this finite curve.
    fn distance_to(self, point: P2) -> f64 {
        let [a, b] = self.ends();
        let ends = distance(point, a).min(distance(point, b));
        match self {
            Self::Line { a, b } => {
                let d = sub(b, a);
                let t = dot(sub(point, a), d) / dot(d, d);
                if (0. ..=1.).contains(&t) {
                    cross(d, sub(point, a)).abs() / d[0].portable_hypot(d[1])
                } else {
                    ends
                }
            }
            Self::Arc {
                center,
                radius,
                start,
                sweep,
            } => {
                let offset = (angle_of(center, point) - start) * sweep.signum();
                if offset.rem_euclid(TAU) <= sweep.abs() {
                    (distance(point, center) - radius).abs()
                } else {
                    ends
                }
            }
        }
    }

    /// Twice the signed area this curve contributes to a closed boundary.
    pub(super) fn twice_area(self) -> f64 {
        match self {
            Self::Line { a, b } => cross(a, b),
            Self::Arc {
                center,
                radius,
                sweep,
                ..
            } => {
                let [p, q] = self.ends();
                center[0] * (q[1] - p[1]) - center[1] * (q[0] - p[0]) + radius * radius * sweep
            }
        }
    }

    /// Pieces monotone in `y`, for crossing parity.
    fn monotone(self) -> Vec<Self> {
        let Self::Arc { start, sweep, .. } = self else {
            return vec![self];
        };
        let mut cuts = vec![0.];
        // Extremes of y at odd multiples of a right angle.
        let [low, high] = if sweep > 0. {
            [start, start + sweep]
        } else {
            [start + sweep, start]
        };
        let mut k = ((low - FRAC_PI_2) / PI).floor() + 1.;
        while FRAC_PI_2 + k * PI < high {
            let angle = FRAC_PI_2 + k * PI;
            cuts.push((angle - start) / sweep);
            k += 1.;
        }
        cuts.push(1.);
        cuts.sort_by(f64::total_cmp);
        cuts.windows(2)
            .filter(|w| w[1] > w[0])
            .map(|w| {
                let ends = [self.point(w[0]), self.point(w[1])];
                self.between(w[0], w[1], ends)
            })
            .collect()
    }

    /// Whether the ray from `point` towards +x crosses this y-monotone
    /// curve, half-open in y so shared ends count once.
    fn crosses(self, point: P2) -> bool {
        let [a, b] = self.ends();
        let py = point[1];
        if (a[1] <= py) == (b[1] <= py) {
            return false;
        }
        let x = match self {
            Self::Line { a, b } => a[0] + (py - a[1]) * (b[0] - a[0]) / (b[1] - a[1]),
            Self::Arc {
                center,
                radius,
                start,
                sweep,
            } => {
                let half = (radius * radius - (py - center[1]).powi(2)).max(0.).sqrt();
                let mid = (start + 0.5 * sweep).portable_cos();
                center[0] + half.copysign(mid)
            }
        };
        x > point[0]
    }
}

/// Whether `point` lies inside the closed boundary formed by `curves`.
pub(super) fn contains(curves: impl IntoIterator<Item = Shape>, point: P2) -> bool {
    curves
        .into_iter()
        .flat_map(Shape::monotone)
        .filter(|piece| piece.crosses(point))
        .count()
        % 2
        == 1
}

/// The supporting line or circle shared by coincident curves.
#[derive(Clone, Copy, Debug)]
pub(super) enum Support {
    Line { point: P2, direction: P2 },
    Circle { center: P2, radius: f64 },
}

impl Support {
    fn of(shape: Shape) -> Self {
        match shape {
            Shape::Line { a, b } => {
                let d = sub(b, a);
                let length = d[0].portable_hypot(d[1]);
                Self::Line {
                    point: a,
                    direction: [d[0] / length, d[1] / length],
                }
            }
            Shape::Arc { center, radius, .. } => Self::Circle { center, radius },
        }
    }

    /// Whether the whole of `shape` lies on this support.
    fn carries(self, shape: Shape, tolerance: f64) -> bool {
        match (self, shape) {
            (Self::Line { point, direction }, Shape::Line { a, b }) => [a, b]
                .iter()
                .all(|&p| cross(direction, sub(p, point)).abs() <= tolerance),
            (
                Self::Circle { center, radius },
                Shape::Arc {
                    center: c,
                    radius: r,
                    ..
                },
            ) => distance(center, c) <= tolerance && (radius - r).abs() <= tolerance,
            _ => false,
        }
    }
}

/// One oriented profile curve with its material on the left.
#[derive(Clone, Debug)]
pub(super) struct Segment {
    pub(super) shape: Shape,
    pub(super) operand: usize,
}

/// An elementary piece of the arrangement. `marks` lists the profile
/// segments lying on it and whether each runs in the piece's direction.
#[derive(Clone, Debug)]
pub(super) struct Piece {
    pub(super) shape: Shape,
    pub(super) nodes: [usize; 2],
    pub(super) support: usize,
    pub(super) marks: Vec<(usize, bool)>,
    /// For each operand, whether its material lies [left, right].
    pub(super) sides: [[bool; 2]; 2],
}

pub(super) struct Arrangement {
    pub(super) nodes: Vec<P2>,
    pub(super) supports: Vec<Support>,
    pub(super) pieces: Vec<Piece>,
}

/// Why the arrangement cannot be formed exactly.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Refusal {
    /// Curves touch tangentially or pass within the tolerance of each other
    /// without a separable crossing.
    Contact,
    /// More curves than this provider admits.
    Limit,
}

pub(super) const MAX_SEGMENTS: usize = 512;

/// The node at `point`, added if no node lies within the tolerance; a point
/// near two distinct nodes is an unresolved junction.
fn node(nodes: &mut Vec<P2>, point: P2, tolerance: f64) -> Result<usize, Refusal> {
    let mut near = nodes
        .iter()
        .enumerate()
        .filter(|(_, n)| distance(**n, point) <= tolerance)
        .map(|(i, _)| i);
    match (near.next(), near.next()) {
        (Some(_), Some(_)) => Err(Refusal::Contact),
        (Some(i), None) => Ok(i),
        _ => {
            nodes.push(point);
            Ok(nodes.len() - 1)
        }
    }
}

/// Points where two curves on different supports meet.
fn crossings(a: Shape, b: Shape, tolerance: f64) -> Result<Vec<P2>, Refusal> {
    let points = match (a, b) {
        (Shape::Line { a: p, b: q }, Shape::Line { a: r, b: s }) => {
            let [d, e] = [sub(q, p), sub(s, r)];
            let denominator = cross(d, e);
            let scale = d[0].portable_hypot(d[1]) * e[0].portable_hypot(e[1]);
            if denominator.abs() <= 1e-12 * scale {
                // Parallel lines on different supports never meet.
                return Ok(vec![]);
            }
            let t = cross(sub(r, p), e) / denominator;
            vec![[p[0] + t * d[0], p[1] + t * d[1]]]
        }
        (Shape::Line { a: p, b: q }, Shape::Arc { center, radius, .. })
        | (Shape::Arc { center, radius, .. }, Shape::Line { a: p, b: q }) => {
            let d = sub(q, p);
            let length = d[0].portable_hypot(d[1]);
            let unit = [d[0] / length, d[1] / length];
            let foot_t = dot(sub(center, p), unit);
            let foot = [p[0] + foot_t * unit[0], p[1] + foot_t * unit[1]];
            let offset = distance(center, foot);
            if (offset - radius).abs() <= tolerance {
                return touching(a, b, foot, tolerance);
            }
            if offset > radius {
                return Ok(vec![]);
            }
            let half = (radius * radius - offset * offset).sqrt();
            [-half, half]
                .map(|h| [foot[0] + h * unit[0], foot[1] + h * unit[1]])
                .to_vec()
        }
        (
            Shape::Arc {
                center: c0,
                radius: r0,
                ..
            },
            Shape::Arc {
                center: c1,
                radius: r1,
                ..
            },
        ) => {
            let d = distance(c0, c1);
            if (d - (r0 + r1)).abs() <= tolerance || (d - (r0 - r1).abs()).abs() <= tolerance {
                let towards = sub(c1, c0);
                let unit = [towards[0] / d, towards[1] / d];
                // The tangency point lies on the line of centres.
                let along = if (d - (r0 + r1)).abs() <= tolerance || r0 > r1 {
                    r0
                } else {
                    -r0
                };
                return touching(
                    a,
                    b,
                    [c0[0] + along * unit[0], c0[1] + along * unit[1]],
                    tolerance,
                );
            }
            if d > r0 + r1 || d < (r0 - r1).abs() {
                return Ok(vec![]);
            }
            let along = (d * d + r0 * r0 - r1 * r1) / (2. * d);
            let half = (r0 * r0 - along * along).max(0.).sqrt();
            let towards = sub(c1, c0);
            let unit = [towards[0] / d, towards[1] / d];
            let base = [c0[0] + along * unit[0], c0[1] + along * unit[1]];
            [-half, half]
                .map(|h| [base[0] - h * unit[1], base[1] + h * unit[0]])
                .to_vec()
        }
    };
    Ok(points)
}

/// A tangency is admitted only where it is a shared end of both curves.
fn touching(a: Shape, b: Shape, point: P2, tolerance: f64) -> Result<Vec<P2>, Refusal> {
    let at_end = |shape: Shape| {
        shape
            .ends()
            .iter()
            .any(|&e| distance(e, point) <= tolerance)
    };
    match (a.parameter(point, tolerance), b.parameter(point, tolerance)) {
        (Some(_), Some(_)) if at_end(a) && at_end(b) => Ok(vec![point]),
        (Some(_), Some(_)) => Err(Refusal::Contact),
        _ => Ok(vec![]),
    }
}

impl Arrangement {
    /// Split `segments` of two operands against each other. Each operand's
    /// segments form closed boundaries with the material on their left.
    pub(super) fn new(segments: &[Segment], tolerance: f64) -> Result<Self, Refusal> {
        if segments.len() > MAX_SEGMENTS {
            return Err(Refusal::Limit);
        }
        // Coincident supports are shared.
        let mut supports: Vec<Support> = vec![];
        let support_of: Vec<usize> = segments
            .iter()
            .map(|segment| {
                supports
                    .iter()
                    .position(|s| s.carries(segment.shape, tolerance))
                    .unwrap_or_else(|| {
                        supports.push(Support::of(segment.shape));
                        supports.len() - 1
                    })
            })
            .collect();
        let mut cuts: Vec<Vec<(f64, P2)>> = segments
            .iter()
            .map(|s| {
                let [a, b] = s.shape.ends();
                vec![(0., a), (1., b)]
            })
            .collect();
        for i in 0..segments.len() {
            for j in i + 1..segments.len() {
                let [a, b] = [segments[i].shape, segments[j].shape];
                let points = if support_of[i] == support_of[j] {
                    // Overlap ends: each curve's ends that lie on the other.
                    a.ends().into_iter().chain(b.ends()).collect()
                } else {
                    crossings(a, b, tolerance)?
                };
                for point in points {
                    if let (Some(s), Some(t)) =
                        (a.parameter(point, tolerance), b.parameter(point, tolerance))
                    {
                        cuts[i].push((s, point));
                        cuts[j].push((t, point));
                    }
                }
            }
        }
        let mut nodes: Vec<P2> = vec![];
        let mut pieces: Vec<Piece> = vec![];
        let mut by_ends: BTreeMap<(usize, usize, usize), Vec<usize>> = BTreeMap::new();
        for (index, (segment, points)) in segments.iter().zip(&mut cuts).enumerate() {
            points.sort_by(|p, q| p.0.total_cmp(&q.0));
            let mut stops: Vec<(f64, usize)> = vec![];
            for &(t, p) in points.iter() {
                let id = node(&mut nodes, p, tolerance)?;
                if stops.last().is_none_or(|&(_, last)| last != id) {
                    stops.push((t, id));
                }
            }
            if stops.len() < 2 {
                return Err(Refusal::Contact);
            }
            for w in stops.windows(2) {
                let [(t0, n0), (t1, n1)] = [w[0], w[1]];
                let shape = segment.shape.between(t0, t1, [nodes[n0], nodes[n1]]);
                if shape.length() <= tolerance {
                    return Err(Refusal::Contact);
                }
                let middle = shape.point(0.5);
                let key = (support_of[index], n0.min(n1), n0.max(n1));
                let same = by_ends.get(&key).and_then(|list| {
                    list.iter()
                        .copied()
                        .find(|&p| distance(pieces[p].shape.point(0.5), middle) <= tolerance)
                });
                match same {
                    Some(p) => {
                        let forward = pieces[p].nodes == [n0, n1];
                        pieces[p].marks.push((index, forward));
                    }
                    None => {
                        by_ends.entry(key).or_default().push(pieces.len());
                        pieces.push(Piece {
                            shape,
                            nodes: [n0, n1],
                            support: support_of[index],
                            marks: vec![(index, true)],
                            sides: [[false; 2]; 2],
                        });
                    }
                }
            }
        }
        let mut arrangement = Self {
            nodes,
            supports,
            pieces,
        };
        arrangement.classify(segments, tolerance)?;
        Ok(arrangement)
    }

    /// Where each operand's material lies around every piece.
    fn classify(&mut self, segments: &[Segment], tolerance: f64) -> Result<(), Refusal> {
        // A sample must clear every unrelated boundary by this margin.
        let margin = 16. * tolerance;
        for piece in &mut self.pieces {
            for operand in 0..2 {
                let marks: Vec<bool> = piece
                    .marks
                    .iter()
                    .filter(|&&(s, _)| segments[s].operand == operand)
                    .map(|&(_, forward)| forward)
                    .collect();
                piece.sides[operand] = match marks.as_slice() {
                    [forward] => [*forward, !*forward],
                    [] => {
                        let middle = piece.shape.point(0.5);
                        let own = segments.iter().filter(|s| s.operand == operand);
                        if own.clone().any(|s| s.shape.distance_to(middle) <= margin) {
                            return Err(Refusal::Contact);
                        }
                        let inside = contains(own.map(|s| s.shape), middle);
                        [inside, inside]
                    }
                    // One operand's own boundary never overlaps itself.
                    _ => return Err(Refusal::Contact),
                };
            }
        }
        Ok(())
    }
}

/// The curves of one oriented ring, a full circle split into two halves.
pub(super) fn ring_shapes(
    ring: &[multicad_rust_kernel::normal_sweep::Edge],
    tolerance: f64,
) -> Option<Vec<Shape>> {
    use multicad_rust_kernel::normal_sweep::Edge;
    let mut shapes = vec![];
    for &edge in ring {
        match edge {
            Edge::Line { a, b } => shapes.push(Shape::Line { a, b }),
            Edge::Arc {
                center,
                a,
                b,
                clockwise,
            } => {
                let radius = distance(a, center);
                if (distance(b, center) - radius).abs() > tolerance || distance(a, b) <= tolerance {
                    return None;
                }
                let start = angle_of(center, a);
                let end = angle_of(center, b);
                let sweep = if clockwise {
                    -(start - end).rem_euclid(TAU)
                } else {
                    (end - start).rem_euclid(TAU)
                };
                shapes.push(Shape::Arc {
                    center,
                    radius,
                    start,
                    sweep,
                });
            }
            Edge::Circle {
                center,
                radius_mm,
                clockwise,
            } => {
                let sweep = if clockwise { -PI } else { PI };
                for start in [0., sweep] {
                    shapes.push(Shape::Arc {
                        center,
                        radius: radius_mm,
                        start,
                        sweep,
                    });
                }
            }
        }
    }
    Some(shapes)
}
