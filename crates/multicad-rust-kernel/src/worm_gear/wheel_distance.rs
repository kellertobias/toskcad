//! Certified distance from a point to one generated wheel patch.
//!
//! A patch point is `S(u, v) = (P(u, v), v)` over `[0, 1] × [-h, h]` in
//! sector-0 module units. Branch and bound covers the domain with parameter
//! cells. On a cell `B` with center `c` and half-widths `(ρu, ρv)` the mean
//! value theorem bounds the squared distance `g = |Q - S|²` from below by
//!
//! `g(c) - ρu max|∂g/∂u (B)| - ρv max|∂g/∂v (B)|`,
//!
//! with `S(c)` and the derivatives enclosed by outward-rounded interval
//! Taylor arithmetic (`interval.rs`), including the flank exit radius
//! (`wheel_bounds.rs`). The bound is second order where the gradient
//! vanishes, so cells near a nondegenerate closest point stay few. The
//! enclosure of `S(B)` itself gives a zero-order bound used when larger.
//! Every cell center evaluates an actual patch point, whose distance bounds
//! the minimum from above. The search ends when the smallest open lower
//! bound meets that upper bound within the requested gap.
use super::interval::Interval;
use super::jet::Jet;
use super::wheel::{Piece, Wheel, turn};
use super::wheel_bounds::{Pair, cells};
use std::cmp::Ordering;
use std::collections::BinaryHeap;

/// Cells evaluated before the search reports a limit.
const MAX_CELLS: usize = 400_000;
const START: [usize; 2] = [8, 8];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    Invalid,
    /// The patch bounds could not be proved, or the gap did not close.
    Unresolved,
    Limit,
}

/// Closest patch point within a certified gap, in module units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PatchDistance {
    /// Distance from the query to `point`, an actual patch point.
    pub distance: f64,
    /// Proven lower bound of the distance from the query to the patch.
    pub lower: f64,
    /// Patch point in sector-0 coordinates, `z = v`.
    pub point: [f64; 3],
}

struct Cell {
    /// Lower bound of the squared distance over the cell.
    lower: f64,
    u: Interval,
    v: Interval,
    /// Mean-value contributions along u and v, which choose the split.
    spread: [f64; 2],
}

impl PartialEq for Cell {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}
impl Eq for Cell {}
impl PartialOrd for Cell {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Cell {
    // The heap pops the smallest lower bound first.
    fn cmp(&self, other: &Self) -> Ordering {
        other.lower.total_cmp(&self.lower)
    }
}

fn middle(i: Interval) -> f64 {
    (i.lo + i.hi) / 2.
}

/// Interval distance between two enclosures, zero when they meet.
fn gap(a: Interval, b: Interval) -> f64 {
    (b.lo - a.hi).max(a.lo - b.hi).max(0.)
}

impl Wheel {
    /// Patch point and its cutter-angle at one parameter pair, in f64.
    fn patch_point(&self, piece: Piece, u: f64, v: f64) -> ([f64; 3], f64) {
        let crossing = piece.flank_side().map_or(0., |side| self.crossing(v, side));
        let (q, angle) = self.polar_patch(piece, u, v, crossing);
        let [x, y] = turn(q, angle);
        ([x, y, v], angle)
    }

    /// Gauss-Newton steps toward the closest patch point from `(u, v)`,
    /// clamped to the domain. Every accepted step is an actual patch point
    /// no farther from the query.
    fn polish(
        &self,
        piece: Piece,
        query: [f64; 3],
        mut at: [f64; 2],
        mut best: f64,
    ) -> ([f64; 2], f64) {
        let h = self.half_width();
        for _ in 0..16 {
            let [u, v] = at;
            let crossing = piece.flank_side().map_or(
                Jet {
                    value: 0.,
                    d: [0.; 2],
                },
                |side| {
                    let slope = self.crossing_slope(v, side).0;
                    Jet {
                        value: slope[0],
                        d: [0., slope[1]],
                    }
                },
            );
            let [x, y] = self.patch(piece, Jet::variable(u, 0), Jet::variable(v, 1), crossing);
            let r = [query[0] - x.value, query[1] - y.value, query[2] - v];
            let su = [x.d[0], y.d[0], 0.];
            let sv = [x.d[1], y.d[1], 1.];
            let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
            let [a, b, c] = [dot(su, su), dot(su, sv), dot(sv, sv)];
            let [p, q] = [dot(su, r), dot(sv, r)];
            let det = a * c - b * b;
            if !(det.is_finite() && det > 0.) {
                break;
            }
            let next = [
                (u + (c * p - b * q) / det).clamp(0., 1.),
                (v + (a * q - b * p) / det).clamp(-h, h),
            ];
            let point = self.patch_point(piece, next[0], next[1]).0;
            let distance = (0..3)
                .map(|k| (query[k] - point[k]).powi(2))
                .sum::<f64>()
                .sqrt();
            if distance.is_nan() || distance >= best {
                break;
            }
            (at, best) = (next, distance);
        }
        (at, best)
    }

    /// Lower bound of `g` over one cell, its split weights and the center
    /// point, or `None` when the enclosures are not proved.
    fn bound_cell(
        &self,
        piece: Piece,
        query: [f64; 3],
        u: Interval,
        v: Interval,
    ) -> Option<(Cell, [f64; 3])> {
        let mid = [middle(u), middle(v)];
        let (point, base) = self.patch_point(piece, mid[0], mid[1]);
        if !point.iter().all(|x| x.is_finite()) {
            return None;
        }
        // Rotate the query, not the patch, so the enclosures stay tight.
        let turned = turn(
            [Interval::point(query[0]), Interval::point(query[1])],
            Interval::point(-base),
        );
        let q = [turned[0], turned[1], Interval::point(query[2])];
        let center = self.coefficients::<2, 2>(
            piece,
            mid.map(Interval::point),
            self.exit(piece, Interval::point(mid[1]))?,
            base,
        );
        let cell: [Pair<2, 2>; 2] = self.coefficients(piece, [u, v], self.exit(piece, v)?, base);
        let value = |p: &Pair<2, 2>| p.0[0].0[0];
        let c = [
            value(&center[0]),
            value(&center[1]),
            Interval::point(mid[1]),
        ];
        let g = (0..3).fold(Interval::point(0.), |sum, k| {
            let d = q[k] - c[k];
            sum + d * d
        });
        let s = [value(&cell[0]), value(&cell[1]), v];
        let su = [cell[0].0[1].0[0], cell[1].0[1].0[0], Interval::point(0.)];
        let sv = [cell[0].0[0].0[1], cell[1].0[0].0[1], Interval::point(1.)];
        let slope = |ds: [Interval; 3]| {
            (0..3)
                .fold(Interval::point(0.), |sum, k| sum + (q[k] - s[k]) * ds[k])
                .magnitude()
                * 2.
        };
        let radius = [u, v].map(|x| (x.hi - x.lo) / 2.);
        let spread = [slope(su) * radius[0], slope(sv) * radius[1]];
        let mean = g.lo - (spread[0] + spread[1]).next_up();
        let zero = (0..3)
            .map(|k| {
                let d = gap(q[k], s[k]);
                d * d
            })
            .sum::<f64>()
            .next_down();
        let lower = mean.max(zero).max(0.);
        if !lower.is_finite() || !spread.iter().all(|x| x.is_finite()) {
            return None;
        }
        let point = [point[0], point[1], point[2]];
        Some((
            Cell {
                lower,
                u,
                v,
                spread,
            },
            point,
        ))
    }

    /// Closest point of one sector-0 patch to `query` (module units), with
    /// a proven lower bound no more than `gap` below the returned distance.
    pub fn patch_distance(
        &self,
        piece: Piece,
        query: [f64; 3],
        target_gap: f64,
    ) -> Result<PatchDistance, Error> {
        if !query.iter().all(|x| x.is_finite()) || target_gap.is_nan() || target_gap <= 0. {
            return Err(Error::Invalid);
        }
        let h = self.half_width();
        let mut best: Option<(f64, [f64; 2])> = None;
        let mut heap = BinaryHeap::new();
        let mut evaluated = 0;
        let mut add = |heap: &mut BinaryHeap<Cell>,
                       best: &mut Option<(f64, [f64; 2])>,
                       u: Interval,
                       v: Interval|
         -> Result<(), Error> {
            evaluated += 1;
            if evaluated > MAX_CELLS {
                return Err(Error::Limit);
            }
            let (cell, point) = self
                .bound_cell(piece, query, u, v)
                .ok_or(Error::Unresolved)?;
            let distance = (0..3)
                .map(|k| (query[k] - point[k]).powi(2))
                .sum::<f64>()
                .sqrt();
            if best.is_none_or(|(d, _)| distance < d) {
                *best = Some((distance, [middle(u), middle(v)]));
            }
            heap.push(cell);
            Ok(())
        };
        for v in cells(-h, h, START[1]) {
            for u in cells(0., 1., START[0]) {
                add(&mut heap, &mut best, u, v)?;
            }
        }
        while let Some(cell) = heap.pop() {
            let (distance, at) = best.ok_or(Error::Unresolved)?;
            let lower = cell.lower.sqrt();
            if lower >= distance - target_gap {
                let (at, distance) = self.polish(piece, query, at, distance);
                let point = self.patch_point(piece, at[0], at[1]).0;
                return Ok(PatchDistance {
                    distance,
                    lower: lower.min(distance),
                    point,
                });
            }
            let halves = |i: Interval| {
                let m = middle(i);
                [Interval::new(i.lo, m), Interval::new(m, i.hi)]
            };
            if cell.spread[0] >= cell.spread[1] {
                for u in halves(cell.u) {
                    add(&mut heap, &mut best, u, cell.v)?;
                }
            } else {
                for v in halves(cell.v) {
                    add(&mut heap, &mut best, cell.u, v)?;
                }
            }
        }
        Err(Error::Unresolved)
    }
}

#[cfg(test)]
mod tests {
    use super::super::wheel::{
        PIECES,
        tests::{spec, wheel},
    };
    use super::*;

    #[test]
    fn offset_patch_points_measure_their_offset_and_far_points_stay_bounded() {
        let wheel = wheel(spec(1., 32, 2., 2., 0.05));
        for piece in PIECES {
            for (u, v) in [(0.3, -0.4), (0.5, 0.), (0.8, 0.7)] {
                let crossing = piece.flank_side().map_or(0., |s| wheel.crossing(v, s));
                let p = |u: f64, v: f64| {
                    let c = piece.flank_side().map_or(0., |s| wheel.crossing(v, s));
                    let [x, y] = wheel.patch(piece, u, v, c);
                    [x, y, v]
                };
                let [x, y] = wheel.patch(piece, u, v, crossing);
                let at = [x, y, v];
                // Unit normal from central differences of the exact patch.
                let e = 1e-6;
                let du: [f64; 3] = std::array::from_fn(|k| (p(u + e, v)[k] - p(u - e, v)[k]) / 2.);
                let dv: [f64; 3] = std::array::from_fn(|k| (p(u, v + e)[k] - p(u, v - e)[k]) / 2.);
                let n = [
                    du[1] * dv[2] - du[2] * dv[1],
                    du[2] * dv[0] - du[0] * dv[2],
                    du[0] * dv[1] - du[1] * dv[0],
                ];
                let length = n.iter().map(|x| x * x).sum::<f64>().sqrt();
                for side in [1., -1.] {
                    let query = std::array::from_fn(|k| at[k] + side * 0.01 * n[k] / length);
                    let found = wheel.patch_distance(piece, query, 1e-7).unwrap();
                    assert!(found.lower <= found.distance);
                    assert!(found.distance - found.lower <= 1e-7);
                    // An interior offset is the distance when no other part
                    // of the patch comes closer, so the certified bracket
                    // holds it (the difference normal is good to ~1e-10).
                    assert!(found.lower <= 0.01 + 1e-9, "{piece:?} {u} {v} {found:?}");
                    assert!(found.distance >= 0.01 - 1e-9, "{piece:?} {u} {v} {found:?}");
                    let moved = (0..3)
                        .map(|k| (found.point[k] - at[k]).powi(2))
                        .sum::<f64>()
                        .sqrt();
                    assert!(moved < 1e-7, "{piece:?} {u} {v}: {moved}");
                }
            }
        }
        // A point beyond the patch reaches its boundary within the gap.
        let far = wheel
            .patch_distance(Piece::Torus, [0., 0., 0.], 1e-7)
            .unwrap();
        assert!(far.distance - far.lower <= 1e-7);
        assert!(far.distance > 10.);
        assert_eq!(
            wheel.patch_distance(Piece::Torus, [f64::NAN, 0., 0.], 1e-7),
            Err(Error::Invalid)
        );
        assert_eq!(
            wheel.patch_distance(Piece::Torus, [0.; 3], 0.),
            Err(Error::Invalid)
        );
    }
}
