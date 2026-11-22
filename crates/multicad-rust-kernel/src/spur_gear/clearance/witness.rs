//! Nearest actual boundary points of two placed exact outlines. Both points
//! are evaluations of the analytic edges, so their distance is an upper
//! bound on the true boundary gap. Chords only steer the search and never
//! stand in for either manufacturing boundary.
use super::{Error, Placement, checkpoint, segment_count};
use crate::{
    revolution::{dot, length},
    spur_gear::{Edge, Point, add, scale, sub},
};

const MAX_WITNESS_CHORDS: usize = 65_536;
const MAX_WITNESS_TESTS: usize = 4_000_000;
const REFINEMENT_ROUNDS: usize = 8;
const GOLDEN_STEPS: usize = 48;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Witness {
    /// Profile-local points, evaluated on the exact edges.
    pub local_mm: [Point; 2],
    /// The caller's index of the edge holding each point.
    pub edges: [usize; 2],
    /// Distance between the two placed points.
    pub distance_mm: f64,
}

#[derive(Clone, Copy)]
struct Span {
    edge: usize,
    t: [f64; 2],
    ends: [Point; 2],
    min_x: f64,
}

fn spans(
    edges: &[(usize, Edge)],
    placement: Placement,
    deflection: f64,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Vec<Span>, Error> {
    let mut result = Vec::new();
    for (slot, &(_, edge)) in edges.iter().enumerate() {
        checkpoint(interrupt)?;
        let count = segment_count(edge, deflection)?;
        if result.len() + count > MAX_WITNESS_CHORDS {
            return Err(Error::Limit);
        }
        for i in 0..count {
            let t = [i as f64 / count as f64, (i + 1) as f64 / count as f64];
            let ends = t.map(|t| placement.point(edge.point(t)));
            result.push(Span {
                edge: slot,
                t,
                ends,
                min_x: ends[0][0].min(ends[1][0]),
            });
        }
    }
    Ok(result)
}

fn point_segment(point: Point, [a, b]: [Point; 2]) -> f64 {
    let direction = sub(b, a);
    let squared = dot(direction, direction);
    let t = if squared > 0. {
        (dot(sub(point, a), direction) / squared).clamp(0., 1.)
    } else {
        0.
    };
    length(sub(point, add(a, scale(direction, t))))
}

pub(super) fn segment_distance(a: [Point; 2], b: [Point; 2]) -> f64 {
    // Crossing chords only steer the search; their exact edge points are
    // refined below, so zero is a valid lower estimate here.
    let side = |p: Point, [q, r]: [Point; 2]| {
        let [x, y] = sub(r, q);
        let [u, v] = sub(p, q);
        x * v - y * u
    };
    if side(a[0], b) * side(a[1], b) < 0. && side(b[0], a) * side(b[1], a) < 0. {
        return 0.;
    }
    point_segment(a[0], b)
        .min(point_segment(a[1], b))
        .min(point_segment(b[0], a))
        .min(point_segment(b[1], a))
}

/// Minimize a distance over one parameter interval by golden sections.
fn golden(interval: [f64; 2], distance: impl Fn(f64) -> f64) -> f64 {
    const RATIO: f64 = 0.618_033_988_749_894_9;
    let [mut low, mut high] = interval;
    let mut a = high - RATIO * (high - low);
    let mut b = low + RATIO * (high - low);
    let (mut fa, mut fb) = (distance(a), distance(b));
    for _ in 0..GOLDEN_STEPS {
        if fa <= fb {
            high = b;
            b = a;
            fb = fa;
            a = high - RATIO * (high - low);
            fa = distance(a);
        } else {
            low = a;
            a = b;
            fa = fb;
            b = low + RATIO * (high - low);
            fb = distance(b);
        }
    }
    [low, a, b, high]
        .into_iter()
        .min_by(|x, y| distance(*x).total_cmp(&distance(*y)))
        .expect("four candidates")
}

/// Find a pair of actual boundary points within `reach_mm`, searching the
/// supplied (index, edge) subsets of two placed outlines. `Ok(None)` means
/// no chord pair came within reach of each other, or no exact points did.
pub fn nearest_boundary(
    edges: [&[(usize, Edge)]; 2],
    placements: [Placement; 2],
    reach_mm: f64,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Option<Witness>, Error> {
    if !(reach_mm.is_finite() && reach_mm > 0.)
        || placements.into_iter().any(|p| !p.valid(1e-3))
        || edges.iter().any(|set| set.is_empty())
    {
        return Err(Error::Invalid);
    }
    // segment_count keeps each chord within half its deflection of the
    // exact edge, so two chords near enough include every reachable pair.
    let deflection = (reach_mm * 0.25).clamp(1e-6, 0.01);
    let margin = reach_mm + deflection;
    let first = spans(edges[0], placements[0], deflection, interrupt)?;
    let mut second = spans(edges[1], placements[1], deflection, interrupt)?;
    second.sort_by(|a, b| a.min_x.total_cmp(&b.min_x));
    let span = chord_span(&second);
    let mut best: Option<(f64, Span, Span)> = None;
    let mut tests = 0usize;
    for a in &first {
        checkpoint(interrupt)?;
        let max_x = a.ends[0][0].max(a.ends[1][0]) + margin;
        let min_x = a.min_x - margin - span;
        let start = second.partition_point(|b| b.min_x < min_x);
        for b in second[start..].iter().take_while(|b| b.min_x <= max_x) {
            tests += 1;
            if tests > MAX_WITNESS_TESTS {
                return Err(Error::Limit);
            }
            let distance = segment_distance(a.ends, b.ends);
            if distance <= margin && best.is_none_or(|(d, ..)| distance < d) {
                best = Some((distance, *a, *b));
            }
        }
    }
    let Some((_, a, b)) = best else {
        return Ok(None);
    };
    let curves = [edges[0][a.edge].1, edges[1][b.edge].1];
    let around = |t: [f64; 2]| {
        let width = t[1] - t[0];
        [(t[0] - width).max(0.), (t[1] + width).min(1.)]
    };
    let intervals = [around(a.t), around(b.t)];
    let place = |side: usize, t: f64| placements[side].point(curves[side].point(t));
    let mut t = [a.t, b.t].map(|t| (t[0] + t[1]) / 2.);
    for _ in 0..REFINEMENT_ROUNDS {
        checkpoint(interrupt)?;
        let other = place(1, t[1]);
        t[0] = golden(intervals[0], |s| length(sub(place(0, s), other)));
        let other = place(0, t[0]);
        t[1] = golden(intervals[1], |s| length(sub(place(1, s), other)));
    }
    // Also keep the exact chord endpoints, which are actual points too.
    let candidates = [
        t,
        [a.t[0], b.t[0]],
        [a.t[0], b.t[1]],
        [a.t[1], b.t[0]],
        [a.t[1], b.t[1]],
    ];
    let (distance, t) = candidates
        .into_iter()
        .map(|t| (length(sub(place(0, t[0]), place(1, t[1]))), t))
        .min_by(|x, y| x.0.total_cmp(&y.0))
        .expect("five candidates");
    Ok((distance <= reach_mm).then(|| Witness {
        local_mm: [curves[0].point(t[0]), curves[1].point(t[1])],
        edges: [edges[0][a.edge].0, edges[1][b.edge].0],
        distance_mm: distance,
    }))
}

/// Longest X extent of any chord, so the sorted sweep cannot skip a chord
/// that starts left of the window but still reaches into it.
fn chord_span(spans: &[Span]) -> f64 {
    spans
        .iter()
        .map(|s| (s.ends[0][0] - s.ends[1][0]).abs())
        .fold(0., f64::max)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{normal_sweep, spur_gear::Profile};

    fn square(offset: f64) -> Vec<(usize, Edge)> {
        let p = [[0., 0.], [1., 0.], [1., 1.], [0., 1.]].map(|[x, y]| [x + offset, y]);
        (0..4)
            .map(|i| {
                (
                    10 + i,
                    Edge::Planar(normal_sweep::Edge::Line {
                        a: p[i],
                        b: p[(i + 1) % 4],
                    }),
                )
            })
            .collect()
    }

    #[test]
    fn exact_edge_points_report_the_gap_indices_and_respect_the_reach() {
        let placements = [Placement::default(); 2];
        let [a, b] = [square(0.), square(1.25)];
        let witness = nearest_boundary([&a, &b], placements, 0.3, None)
            .unwrap()
            .unwrap();
        assert!((witness.distance_mm - 0.25).abs() < 1e-12);
        // Bottom, right and top sides tie at the corners; any is actual.
        assert!(witness.edges.iter().all(|e| (10..14).contains(e)));
        assert!((witness.local_mm[0][0] - 1.).abs() < 1e-12);
        assert!((witness.local_mm[1][0] - 1.25).abs() < 1e-12);
        assert_eq!(nearest_boundary([&a, &b], placements, 0.2, None), Ok(None));
    }

    #[test]
    fn curved_edges_refine_to_their_exact_gap() {
        // Two unit circles 0.1 mm apart: their chords lie inside, so the
        // refined points must still reach the exact 0.1 mm gap.
        let circle = |x: f64| {
            vec![(
                0,
                Edge::Planar(normal_sweep::Edge::Circle {
                    center: [x, 0.],
                    radius_mm: 1.,
                    clockwise: false,
                }),
            )]
        };
        let witness = nearest_boundary(
            [&circle(0.), &circle(2.1)],
            [Placement::default(); 2],
            0.2,
            None,
        )
        .unwrap()
        .unwrap();
        assert!((witness.distance_mm - 0.1).abs() < 1e-9, "{witness:?}");
    }

    #[test]
    fn invalid_requests_budgets_and_cancellation_fail_closed() {
        let placements = [Placement::default(); 2];
        let a = square(0.);
        assert_eq!(
            nearest_boundary([&a, &[]], placements, 0.1, None),
            Err(Error::Invalid)
        );
        for reach in [0., f64::NAN] {
            assert_eq!(
                nearest_boundary([&a, &a], placements, reach, None),
                Err(Error::Invalid)
            );
        }
        let mut bad = placements;
        bad[0].uncertainty_mm = 1.;
        assert_eq!(
            nearest_boundary([&a, &a], bad, 0.1, None),
            Err(Error::Invalid)
        );
        assert_eq!(
            nearest_boundary([&a, &a], placements, 0.1, Some(&|| true)),
            Err(Error::Cancelled)
        );
        // A large curved outline at a microscopic reach exceeds the chord budget.
        let pressure = 20_f64.to_radians();
        let big = Profile::new(
            crate::spur_gear::Spec {
                teeth: 64,
                pitch_radius_mm: 64.,
                base_radius_mm: 64. * pressure.cos(),
                tip_radius_mm: 66.,
                root_radius_mm: 61.5,
                root_fillet_radius_mm: 0.4,
                half_base_angle_rad: (std::f64::consts::PI - 0.05) / 128. + pressure.tan()
                    - pressure,
                face_width_mm: 5.,
                bore_radius_mm: 4.,
            },
            1e-6,
            1e7,
        )
        .unwrap();
        let edges: Vec<_> = big
            .outline()
            .iter()
            .enumerate()
            .map(|(i, e)| (i, e.edge))
            .collect();
        assert_eq!(
            nearest_boundary([&edges, &edges], placements, 1e-6, None),
            Err(Error::Limit)
        );
    }
}
