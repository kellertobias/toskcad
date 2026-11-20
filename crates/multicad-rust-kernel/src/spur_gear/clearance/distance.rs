//! Certified minimum distance between two disjoint placed exact outlines.
//!
//! The upper bound is a pair of actual edge evaluations (`nearest_boundary`).
//! The lower bound is the least chord-pair distance minus both chord-error
//! tubes. Coarse chords first exclude every edge pair whose lower bound
//! exceeds the current upper bound; only the remaining chord intervals are
//! refined, level by level, down to the tolerance.
//! A witness is returned only when both bounds agree within the tolerance.
use super::{Error, Placement, checkpoint, segment_count, witness};
use crate::{
    revolution::length,
    spur_gear::{Edge, Point, sub},
};
use std::collections::BTreeSet;

const MAX_CHORDS: usize = 262_144;
const MAX_TESTS: usize = 8_000_000;

struct Chord {
    slot: usize,
    /// Edge parameters of the chord ends.
    t: [f64; 2],
    ends: [Point; 2],
    min_x: f64,
    max_x: f64,
}

/// Chords over parameter intervals `(slot, [t0, t1])` of the given edges.
/// Each interval gets its share of the whole-edge chord count; chord error
/// shrinks with the square of the span, so the deflection bound still holds.
fn chords(
    edges: &[(usize, Edge)],
    intervals: &[(usize, [f64; 2])],
    placement: Placement,
    deflection: f64,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Vec<Chord>, Error> {
    let mut result = Vec::new();
    for &(slot, [t0, t1]) in intervals {
        checkpoint(interrupt)?;
        let edge = edges[slot].1;
        let whole = segment_count(edge, deflection)? as f64;
        let count = (whole * (t1 - t0)).ceil().max(1.) as usize;
        if result.len() + count > MAX_CHORDS {
            return Err(Error::Limit);
        }
        let points: Vec<_> = (0..=count)
            .map(|i| {
                let t = if i == count {
                    t1
                } else {
                    t0 + (t1 - t0) * i as f64 / count as f64
                };
                (t, placement.point(edge.point(t)))
            })
            .collect();
        result.extend(points.windows(2).map(|pair| Chord {
            slot,
            t: [pair[0].0, pair[1].0],
            ends: [pair[0].1, pair[1].1],
            min_x: pair[0].1[0].min(pair[1].1[0]),
            max_x: pair[0].1[0].max(pair[1].1[0]),
        }));
    }
    result.sort_by(|a, b| a.min_x.total_cmp(&b.min_x));
    Ok(result)
}

/// Every chord pair within `reach` of each other, by an X sweep.
fn near_pairs(
    first: &[Chord],
    second: &[Chord],
    reach: f64,
    mut visit: impl FnMut(&Chord, &Chord, f64),
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<(), Error> {
    let span = second.iter().map(|c| c.max_x - c.min_x).fold(0., f64::max);
    let mut tests = 0usize;
    for a in first {
        checkpoint(interrupt)?;
        let start = second.partition_point(|b| b.min_x < a.min_x - reach - span);
        for b in second[start..]
            .iter()
            .take_while(|b| b.min_x <= a.max_x + reach)
        {
            tests += 1;
            if tests > MAX_TESTS {
                return Err(Error::Limit);
            }
            let d = witness::segment_distance(a.ends, b.ends);
            if d <= reach {
                visit(a, b, d);
            }
        }
    }
    Ok(())
}

/// Parameters in [0, 1] of a closest point pair of two segments.
fn closest_parameters(a: [Point; 2], b: [Point; 2]) -> [f64; 2] {
    let dot = |p: Point, q: Point| p[0] * q[0] + p[1] * q[1];
    let (d1, d2, r) = (sub(a[1], a[0]), sub(b[1], b[0]), sub(a[0], b[0]));
    let (aa, ee, f) = (dot(d1, d1), dot(d2, d2), dot(d2, r));
    let clamp = |v: f64| if v.is_finite() { v.clamp(0., 1.) } else { 0. };
    let candidates = [0., 1.]
        .into_iter()
        .map(|s| [s, clamp((f + s * dot(d1, d2)) / ee)])
        .chain(
            [0., 1.]
                .into_iter()
                .map(|u| [clamp((u * dot(d1, d2) - dot(d1, r)) / aa), u]),
        )
        .chain({
            let (c, bb) = (dot(d1, r), dot(d1, d2));
            let denominator = aa * ee - bb * bb;
            let s = clamp((bb * f - c * ee) / denominator);
            [[s, clamp((bb * s + f) / ee)]]
        });
    let at = |p: [Point; 2], t: f64| {
        [
            p[0][0] + (p[1][0] - p[0][0]) * t,
            p[0][1] + (p[1][1] - p[0][1]) * t,
        ]
    };
    candidates
        .min_by(|x, y| {
            length(sub(at(a, x[0]), at(b, x[1]))).total_cmp(&length(sub(at(a, y[0]), at(b, y[1]))))
        })
        .expect("five candidates")
}

/// Chord-to-edge error at a deflection: half the deflection by
/// `segment_count`, the placement uncertainty and a rounding guard.
fn tube(deflection: f64, placement: Placement, scale: f64) -> f64 {
    deflection * 0.5 + placement.uncertainty_mm + 1024. * f64::EPSILON * scale
}

/// Certify the minimum distance between the two boundaries. The caller must
/// separately prove the enclosed regions disjoint, so that the region gap is
/// this boundary gap. `Ok(None)` means the bounds did not meet in budget.
pub fn certify_boundary_distance(
    edges: [&[(usize, Edge)]; 2],
    placements: [Placement; 2],
    tolerance_mm: f64,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Option<witness::Witness>, Error> {
    if !(1e-7..=1e-3).contains(&tolerance_mm)
        || placements.into_iter().any(|p| !p.valid(1e-3))
        || edges.iter().any(|set| set.is_empty())
    {
        return Err(Error::Invalid);
    }
    let scale = placements
        .iter()
        .flat_map(|p| p.center_mm)
        .map(f64::abs)
        .fold(1., f64::max)
        * 4.;
    let mut candidates: [Vec<(usize, [f64; 2])>; 2] =
        edges.map(|set| (0..set.len()).map(|slot| (slot, [0., 1.])).collect());
    // A first actual-point upper bound from the least coarse chord gap.
    let coarse = 0.01;
    let tubes = placements.map(|p| tube(coarse, p, scale));
    let [a, b] = [0, 1].map(|i| chords(edges[i], &candidates[i], placements[i], coarse, interrupt));
    let (a, b) = (a?, b?);
    let mut least = f64::INFINITY;
    for window in [0.05, 0.5, 5., 50., 500.] {
        near_pairs(&a, &b, window, |_, _, d| least = least.min(d), interrupt)?;
        if least.is_finite() {
            break;
        }
    }
    if !least.is_finite() {
        return Ok(None);
    }
    let reach = least + tubes[0] + tubes[1] + coarse;
    let Some(mut witness) = witness::nearest_boundary(edges, placements, reach, interrupt)? else {
        return Ok(None);
    };
    // Each level keeps only edges whose chords come within the current
    // upper bound plus both tubes: all others are provably farther. The
    // final chords lie within a quarter tolerance of their edges, so exact
    // points at the closest chord pair's parameters exceed the final lower
    // bound by at most the tolerance.
    let mut lower = 0.;
    for deflection in [coarse, 1e-3, 1e-4, 1e-5, 1e-6, tolerance_mm * 0.45] {
        let deflection = deflection.max(tolerance_mm * 0.45);
        let tubes = placements.map(|p| tube(deflection, p, scale));
        let [a, b] = [0, 1].map(|i| {
            chords(
                edges[i],
                &candidates[i],
                placements[i],
                deflection,
                interrupt,
            )
        });
        let (a, b) = (a?, b?);
        let reach = witness.distance_mm + tubes[0] + tubes[1];
        let mut kept = [BTreeSet::new(), BTreeSet::new()];
        let mut least = f64::INFINITY;
        let mut best = None;
        near_pairs(
            &a,
            &b,
            reach,
            |x, y, d| {
                kept[0].insert((x.slot, x.t.map(f64::to_bits)));
                kept[1].insert((y.slot, y.t.map(f64::to_bits)));
                if d < least {
                    least = d;
                    best = Some((x.slot, x.t, x.ends, y.slot, y.t, y.ends));
                }
            },
            interrupt,
        )?;
        let Some((x, xt, xe, y, yt, ye)) = best else {
            return Ok(None);
        };
        lower = (least - tubes[0] - tubes[1]).max(0.);
        let [s, u] = closest_parameters(xe, ye);
        let parameters = [xt[0] + (xt[1] - xt[0]) * s, yt[0] + (yt[1] - yt[0]) * u];
        let local = [
            edges[0][x].1.point(parameters[0]),
            edges[1][y].1.point(parameters[1]),
        ];
        let gap = length(sub(
            placements[0].point(local[0]),
            placements[1].point(local[1]),
        ));
        if gap < witness.distance_mm {
            witness = witness::Witness {
                local_mm: local,
                edges: [edges[0][x].0, edges[1][y].0],
                distance_mm: gap,
            };
        }
        candidates = kept.map(|set| {
            set.into_iter()
                .map(|(slot, t)| (slot, t.map(f64::from_bits)))
                .collect()
        });
    }
    Ok((witness.distance_mm - lower <= tolerance_mm).then_some(witness))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::normal_sweep;

    fn square(offset: f64) -> Vec<(usize, Edge)> {
        let p = [[0., 0.], [1., 0.], [1., 1.], [0., 1.]].map(|[x, y]| [x + offset, y]);
        (0..4)
            .map(|i| {
                (
                    i,
                    Edge::Planar(normal_sweep::Edge::Line {
                        a: p[i],
                        b: p[(i + 1) % 4],
                    }),
                )
            })
            .collect()
    }

    fn disc(center: [f64; 2], radius: f64) -> Vec<(usize, Edge)> {
        let at = |a: f64| [center[0] + radius * a.cos(), center[1] + radius * a.sin()];
        (0..4)
            .map(|i| {
                let a = std::f64::consts::FRAC_PI_2 * i as f64;
                (
                    i,
                    Edge::Planar(normal_sweep::Edge::Arc {
                        center,
                        a: at(a),
                        b: at(a + std::f64::consts::FRAC_PI_2),
                        clockwise: false,
                    }),
                )
            })
            .collect()
    }

    #[test]
    fn straight_and_curved_gaps_are_certified_within_tolerance() {
        let p = [Placement::default(); 2];
        let found = certify_boundary_distance([&square(0.), &square(1.25)], p, 1e-7, None)
            .unwrap()
            .unwrap();
        assert!((found.distance_mm - 0.25).abs() <= 1e-7);
        let circle = disc([4., 0.5], 1.);
        let found = certify_boundary_distance([&square(0.), &circle], p, 1e-7, None)
            .unwrap()
            .unwrap();
        // The arc's nearest point lies straight across from the right side.
        assert!(
            (found.distance_mm - 2.).abs() <= 1e-7,
            "{}",
            found.distance_mm
        );
        assert!((found.local_mm[1][0] - 3.).abs() < 1e-6);
    }

    #[test]
    fn invalid_requests_fail_closed() {
        let p = [Placement::default(); 2];
        assert_eq!(
            certify_boundary_distance([&square(0.), &[]], p, 1e-7, None),
            Err(Error::Invalid)
        );
        assert_eq!(
            certify_boundary_distance([&square(0.), &square(2.)], p, 1., None),
            Err(Error::Invalid)
        );
    }
}
