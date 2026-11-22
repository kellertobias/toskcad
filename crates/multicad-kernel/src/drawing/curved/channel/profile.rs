//! Decompose a closed straight-sweep section into exact blocker pieces.
//!
//! Paired concentric quarter arcs (bends) become annular sectors and convex
//! corner quarters become solid sectors. Their radial chords cancel against
//! the remaining straight boundary, leaving rectangles (channel walls) and
//! other simple polygons (the body of a rounded plate). Circular and
//! rounded-slot hole rings must each lie clear inside one such piece.
use multicad_precision::PortableMath as _;
use multicad_rust_kernel::normal_sweep::Edge;
use std::f64::consts::{FRAC_PI_2, PI, TAU};

type Point = [f64; 2];
type Segment = [Point; 2];

pub(super) struct Bend {
    pub center: Point,
    /// The centre as authored, which a snapped corner quarter replaces.
    pub authored: Point,
    pub spokes: [Point; 2],
    /// `None` for a solid corner quarter.
    pub inner: Option<f64>,
    pub outer: f64,
}

/// An axial through-hole and the piece that contains it: a circle when
/// both centres coincide, otherwise a rounded slot (stadium).
pub(super) struct Hole {
    pub centers: [Point; 2],
    pub radius: f64,
    /// Index into the rectangles followed by the polygons.
    pub host: usize,
}

pub(super) struct Profile {
    pub bends: Vec<Bend>,
    pub rectangles: Vec<Segment>,
    /// Simple counterclockwise straight-edged pieces that are not rectangles.
    pub polygons: Vec<Vec<Point>>,
    pub holes: Vec<Hole>,
}

#[derive(Clone, Copy)]
struct Arc {
    center: Point,
    ends: Segment,
    spokes: Segment,
    radius: f64,
}

impl Arc {
    fn new(edge: Edge) -> Option<Self> {
        let Edge::Arc { center, a, b, .. } = edge else {
            return None;
        };
        if (edge.signed_sweep_rad()?.abs() - FRAC_PI_2).abs() > 32. * f64::EPSILON {
            return None;
        }
        let radius = edge.radius_mm()?;
        let spokes = [a, b].map(|p| std::array::from_fn(|i| (p[i] - center[i]) / radius));
        Some(Self {
            center,
            ends: [a, b],
            spokes,
            radius,
        })
    }
}

impl Profile {
    /// `linear_tolerance` bounds how far a convex corner arc may deviate
    /// radially from the exact cardinal quarter that replaces it.
    pub(super) fn new(rings: &[Vec<Edge>], linear_tolerance: f64) -> Option<Self> {
        let (ring, hole_rings) = rings.split_first()?;
        if ring.len() > 256 || hole_rings.len() > 256 {
            return None;
        }
        let scale = rings
            .iter()
            .flatten()
            .flat_map(|edge| edge.start())
            .map(f64::abs)
            .fold(1., f64::max);
        let tolerance = 64. * f64::EPSILON * scale;
        let voids = hole_rings
            .iter()
            .map(|ring| match ring[..] {
                [
                    Edge::Circle {
                        center, radius_mm, ..
                    },
                ] => Some(([center; 2], radius_mm)),
                _ => slot(ring, tolerance),
            })
            .collect::<Option<Vec<_>>>()?;
        let mut lines = Vec::new();
        let mut inner = Vec::new();
        let mut outer = Vec::new();
        for &edge in ring {
            match edge {
                Edge::Line { a, b } => lines.push([a, b]),
                Edge::Arc {
                    clockwise: true, ..
                } => inner.push(Arc::new(edge)?),
                Edge::Arc {
                    clockwise: false, ..
                } => outer.push(edge),
                Edge::Circle { .. } => return None,
            }
        }
        if outer.is_empty() {
            return None;
        }
        let mut bends = Vec::new();
        for edge in outer {
            let paired = Arc::new(edge).and_then(|arc| {
                let i = inner.iter().position(|candidate| {
                    same(candidate.center, arc.center, tolerance)
                        && arc.radius - candidate.radius > tolerance
                        && same(candidate.spokes[0], arc.spokes[1], 32. * f64::EPSILON)
                        && same(candidate.spokes[1], arc.spokes[0], 32. * f64::EPSILON)
                })?;
                Some((arc, inner.swap_remove(i)))
            });
            if let Some((arc, nested)) = paired {
                // Opposite orientation to the radial walls of the annular sector.
                lines.extend([[nested.ends[0], arc.ends[1]], [arc.ends[0], nested.ends[1]]]);
                bends.push(Bend {
                    center: arc.center,
                    authored: arc.center,
                    spokes: arc.spokes,
                    inner: Some(nested.radius),
                    outer: arc.radius,
                });
            } else {
                let (bend, ends) = corner(edge, linear_tolerance)?;
                lines.extend([[ends[0], bend.center], [bend.center, ends[1]]]);
                bends.push(bend);
            }
        }
        if !inner.is_empty() {
            return None;
        }

        let mut cycles = Vec::new();
        while let Some(first) = lines.pop() {
            let mut cycle = vec![first[0]];
            let mut endpoint = first[1];
            while !same(endpoint, first[0], tolerance) {
                cycle.push(endpoint);
                let i = lines
                    .iter()
                    .position(|edge| same(edge[0], endpoint, tolerance))?;
                endpoint = lines.swap_remove(i)[1];
            }
            cycles.push(cycle);
        }
        // The internal chords cancel, leaving the original simple CCW
        // boundary. Every piece has positive winding, so their total winding
        // proves exact coverage with no overlapping interiors.
        let mut rectangles = Vec::new();
        let mut polygons = Vec::new();
        for cycle in cycles {
            match rectangle(&cycle, tolerance) {
                Some(region) => rectangles.push(region),
                None if simple_counterclockwise(&cycle, tolerance) => polygons.push(cycle),
                None => return None,
            }
        }
        let pieces: Vec<Vec<Point>> = rectangles
            .iter()
            .map(|&[low, high]| vec![low, [high[0], low[1]], high, [low[0], high[1]]])
            .chain(polygons.iter().cloned())
            .collect();
        let holes = voids
            .iter()
            .map(|&(centers, radius)| {
                // The void must clear every edge of exactly one piece, so the
                // opening lies inside one cap face of one blocker.
                let mut hosts = pieces.iter().enumerate().filter(|(_, piece)| {
                    inside(piece, centers[0])
                        && (0..piece.len()).all(|i| {
                            let edge = [piece[i], piece[(i + 1) % piece.len()]];
                            segments_distance(centers, edge) > radius + tolerance
                        })
                });
                let (host, _) = hosts.next()?;
                hosts.next().is_none().then_some(Hole {
                    centers,
                    radius,
                    host,
                })
            })
            .collect::<Option<Vec<_>>>()?;
        Some(Self {
            bends,
            rectangles,
            polygons,
            holes,
        })
    }
}

/// Replace a counterclockwise convex corner arc by the exact cardinal
/// quarter through its endpoints. Imported outlines can carry an arc
/// reconstructed from a rounded midpoint, whose centre is off the corner by
/// a fraction of the file resolution. The replacement is accepted only when
/// every point of the authored arc lies within `tolerance` of the quarter's
/// circle, and the quarter's radius is resolved above that tolerance.
fn corner(edge: Edge, tolerance: f64) -> Option<(Bend, Segment)> {
    let Edge::Arc {
        center,
        a,
        b,
        clockwise: false,
    } = edge
    else {
        return None;
    };
    let sweep = edge.signed_sweep_rad()?;
    let radius = edge.radius_mm()?;
    if !(tolerance.is_finite() && tolerance >= 0. && sweep > 0. && sweep < FRAC_PI_2 * 1.5) {
        return None;
    }
    let cardinal = |p: Point| {
        let [x, y] = sub(p, center);
        if x.abs() >= y.abs() {
            [x.signum(), 0.]
        } else {
            [0., y.signum()]
        }
    };
    let spokes = [cardinal(a), cardinal(b)];
    // A counterclockwise quarter turns its first spoke by +90 degrees.
    if spokes[1] != [-spokes[0][1], spokes[0][0]] {
        return None;
    }
    // The exact centre lies on the radial line through each endpoint.
    let snapped = if spokes[0][0] == 0. {
        [a[0], b[1]]
    } else {
        [b[0], a[1]]
    };
    let radii = [a, b]
        .map(|p| sub(p, snapped))
        .map(|[x, y]| x.portable_hypot(y));
    let quarter = radii[0].midpoint(radii[1]);
    if radii.iter().any(|&r| r <= tolerance) {
        return None;
    }
    // Distances from the snapped centre to the authored circle are extreme
    // at the endpoints or along the line through both centres.
    let offset = sub(center, snapped);
    let distance = offset[0].portable_hypot(offset[1]);
    let start = (a[1] - center[1]).portable_atan2(a[0] - center[0]);
    let mut deviation = radii.iter().map(|r| (r - quarter).abs()).fold(0., f64::max);
    if distance > 0. {
        let toward = offset[1].portable_atan2(offset[0]);
        for (angle, reach) in [
            (toward, radius + distance),
            (toward + PI, radius - distance),
        ] {
            if (angle - start).rem_euclid(TAU) <= sweep {
                deviation = deviation.max((reach - quarter).abs());
            }
        }
    }
    (deviation <= tolerance).then_some((
        Bend {
            center: snapped,
            authored: center,
            spokes,
            inner: None,
            outer: quarter,
        },
        [a, b],
    ))
}

/// A rounded-slot ring: two semicircles of one radius bulging away from
/// each other, joined by straight walls tangent to both.
fn slot(ring: &[Edge], tolerance: f64) -> Option<([Point; 2], f64)> {
    let arcs: Vec<_> = ring
        .iter()
        .filter(|edge| matches!(edge, Edge::Arc { .. }))
        .collect();
    let [first, second] = arcs[..] else {
        return None;
    };
    if ring.len() != 4 {
        return None;
    }
    let round = |edge: &Edge| match *edge {
        Edge::Arc { center, .. } => Some((center, edge.radius_mm()?, edge.signed_sweep_rad()?)),
        _ => None,
    };
    let [(c0, r0, s0), (c1, r1, s1)] = [round(first)?, round(second)?];
    let radius = r0.midpoint(r1);
    let spine = sub(c1, c0);
    let length = spine[0].portable_hypot(spine[1]);
    if (r0 - r1).abs() > tolerance
        || radius <= tolerance
        || length <= tolerance
        || [s0, s1]
            .iter()
            .any(|sweep| (sweep.abs() - PI).abs() > 32. * f64::EPSILON)
    {
        return None;
    }
    let along = |p: Point, c: Point| (p[0] - c[0]) * spine[0] + (p[1] - c[1]) * spine[1];
    for (edge, center, other) in [(first, c0, c1), (second, c1, c0)] {
        // Each semicircle bulges away from the other centre.
        let middle = edge.point(0.5);
        let outward = (middle[0] - center[0]) * (center[0] - other[0])
            + (middle[1] - center[1]) * (center[1] - other[1]);
        if outward <= 0. {
            return None;
        }
    }
    for edge in ring {
        let Edge::Line { a, b } = *edge else { continue };
        // Wall endpoints sit on the radial lines perpendicular to the spine.
        let ends = [[(a, c0), (b, c1)], [(a, c1), (b, c0)]]
            .into_iter()
            .find(|pair| {
                pair.iter()
                    .all(|&(p, c)| along(p, c).abs() <= tolerance * length)
            })?;
        for (p, center) in ends {
            let [x, y] = sub(p, center);
            if (x.portable_hypot(y) - radius).abs() > tolerance {
                return None;
            }
        }
    }
    Some(([c0, c1], radius))
}

/// A cycle with positive area whose non-adjacent edges never meet.
fn simple_counterclockwise(cycle: &[Point], tolerance: f64) -> bool {
    let n = cycle.len();
    if !(3..=512).contains(&n) {
        return false;
    }
    let edge = |i: usize| [cycle[i], cycle[(i + 1) % n]];
    let area: f64 = (0..n)
        .map(|i| {
            let [p, q] = edge(i);
            p[0] * q[1] - q[0] * p[1]
        })
        .sum();
    area > 0.
        && (0..n).all(|i| {
            (i + 2..n)
                .filter(|&j| (j + 1) % n != i)
                .all(|j| segments_distance(edge(i), edge(j)) > tolerance)
        })
}

/// Even-odd containment of a point strictly inside a simple polygon.
fn inside(polygon: &[Point], [x, y]: Point) -> bool {
    let n = polygon.len();
    let mut crossings = false;
    for i in 0..n {
        let [p, q] = [polygon[i], polygon[(i + 1) % n]];
        if (p[1] > y) != (q[1] > y) && x < p[0] + (y - p[1]) * (q[0] - p[0]) / (q[1] - p[1]) {
            crossings = !crossings;
        }
    }
    crossings
}

fn point_segment_distance(p: Point, [a, b]: Segment) -> f64 {
    let d = sub(b, a);
    let length = d[0] * d[0] + d[1] * d[1];
    let t = if length == 0. {
        0.
    } else {
        (((p[0] - a[0]) * d[0] + (p[1] - a[1]) * d[1]) / length).clamp(0., 1.)
    };
    let [x, y] = sub(p, [a[0] + t * d[0], a[1] + t * d[1]]);
    x.portable_hypot(y)
}

/// Distance between two closed segments, either possibly degenerate.
fn segments_distance(first: Segment, second: Segment) -> f64 {
    let cross = |o: Point, a: Point, b: Point| {
        (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    };
    let [a, b] = first;
    let [c, d] = second;
    let [d1, d2, d3, d4] = [
        cross(c, d, a),
        cross(c, d, b),
        cross(a, b, c),
        cross(a, b, d),
    ];
    if d1 * d2 < 0. && d3 * d4 < 0. {
        return 0.;
    }
    [
        point_segment_distance(a, second),
        point_segment_distance(b, second),
        point_segment_distance(c, first),
        point_segment_distance(d, first),
    ]
    .into_iter()
    .fold(f64::INFINITY, f64::min)
}

fn sub(a: Point, b: Point) -> Point {
    [a[0] - b[0], a[1] - b[1]]
}

fn same(a: Point, b: Point, tolerance: f64) -> bool {
    (0..2).all(|i| (a[i] - b[i]).abs() <= tolerance)
}

fn rectangle(ring: &[Point], tolerance: f64) -> Option<Segment> {
    let [a, b, c, d] = <[Point; 4]>::try_from(ring).ok()?;
    let bounds = [false, true].map(|high| {
        std::array::from_fn(|i| {
            ring.iter()
                .map(|p| p[i])
                .reduce(|a, b| if high { a.max(b) } else { a.min(b) })
                .expect("four corners")
        })
    });
    let [low, high] = bounds;
    let corners = [low, [high[0], low[1]], high, [low[0], high[1]]];
    let ordered = [a, b, c, d];
    if (0..2).any(|i| high[i] - low[i] <= tolerance)
        || !(0..4)
            .any(|shift| (0..4).all(|i| same(ordered[i], corners[(i + shift) % 4], tolerance)))
    {
        return None;
    }
    Some(bounds)
}

#[cfg(test)]
mod tests;
