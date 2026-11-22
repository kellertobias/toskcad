//! Exact line intersections with planar profiles of lines, circular arcs and
//! involutes. A family splits every edge where its distance along one fixed
//! normal turns, so that each piece meets a parallel line at most once.
use super::paths::{dot2, edge_point, round};
use super::*;
use multicad_rust_kernel::spur_gear::Edge;
use std::f64::consts::PI;

pub(in crate::drawing) type Rings = Vec<Vec<Edge>>;

/// Fractions strictly inside the edge where its tangent is perpendicular to
/// `normal`, ascending.
pub(in crate::drawing) fn folds(edge: Edge, normal: P2) -> Vec<f64> {
    let (first, start, span) = match edge {
        Edge::Planar(planar) => {
            let Some((_, _, start, sweep)) = round(planar) else {
                return Vec::new();
            };
            (normal[1].atan2(normal[0]), start, sweep)
        }
        Edge::Involute(curve) => {
            let [x, y] = curve.basis().map(|axis| dot2(axis, normal));
            let [a, b] = curve.interval();
            ((-x).atan2(y), a, b - a)
        }
    };
    let [lo, hi] = [start.min(start + span), start.max(start + span)];
    let mut fractions = Vec::new();
    let mut k = ((lo - first) / PI).ceil();
    while first + k * PI <= hi {
        let fraction = (first + k * PI - start) / span;
        if fraction > 1e-12 && fraction < 1. - 1e-12 {
            fractions.push(fraction);
        }
        k += 1.;
    }
    fractions.sort_by(f64::total_cmp);
    fractions
}

#[derive(Clone, Copy, Debug)]
pub(in crate::drawing) struct Piece {
    pub ring: usize,
    pub edge: usize,
    pub range: [f64; 2],
    /// Distance along the family normal at both ends.
    pub level: [f64; 2],
}

#[derive(Clone, Copy, Debug)]
pub(in crate::drawing) struct Hit {
    pub fraction: f64,
    pub point: P2,
}

pub(in crate::drawing) struct Family {
    pub normal: P2,
    pub pieces: Vec<Piece>,
}

impl Family {
    pub(in crate::drawing) fn new(rings: &Rings, normal: P2) -> Self {
        let mut pieces = Vec::new();
        for (ring, edges) in rings.iter().enumerate() {
            for (index, &edge) in edges.iter().enumerate() {
                let mut stops = vec![0.];
                stops.extend(folds(edge, normal));
                stops.push(1.);
                for pair in stops.windows(2) {
                    pieces.push(Piece {
                        ring,
                        edge: index,
                        range: [pair[0], pair[1]],
                        level: [pair[0], pair[1]].map(|f| dot2(edge_point(edge, f), normal)),
                    });
                }
            }
        }
        Self { normal, pieces }
    }

    /// The fraction on a monotone piece at the given level.
    pub(in crate::drawing) fn solve(&self, rings: &Rings, piece: &Piece, level: f64) -> Hit {
        let edge = rings[piece.ring][piece.edge];
        let value = |f: f64| dot2(edge_point(edge, f), self.normal) - level;
        let [mut a, mut b] = piece.range;
        let [mut fa, mut fb] = piece.level.map(|l| l - level);
        let mut fraction = a + (b - a) * fa / (fa - fb);
        if !matches!(edge, Edge::Planar(planar) if round(planar).is_none()) {
            // Illinois false position on the exact curve.
            for _ in 0..200 {
                let f = value(fraction);
                if f.abs() <= 1e-15 * (1. + level.abs())
                    || (b - a).abs() <= 4. * f64::EPSILON * (a.abs() + b.abs())
                {
                    break;
                }
                if (f < 0.) == (fa < 0.) {
                    a = fraction;
                    fa = f;
                    fb *= 0.5;
                } else {
                    b = fraction;
                    fb = f;
                    fa *= 0.5;
                }
                let next = a + (b - a) * fa / (fa - fb);
                fraction = if next.is_finite() && (next - a) * (next - b) <= 0. {
                    next
                } else {
                    0.5 * (a + b)
                };
            }
        }
        Hit {
            fraction,
            point: edge_point(edge, fraction),
        }
    }

    /// Every crossing of the profile with the line at `level` along the
    /// normal. A piece counts when its ends are on opposite sides, the lower
    /// side exclusive, so closed rings always cross an even number of times.
    pub(in crate::drawing) fn crossings(&self, rings: &Rings, level: f64) -> Vec<Hit> {
        self.pieces
            .iter()
            .filter(|piece| (piece.level[0] < level) != (piece.level[1] < level))
            .map(|piece| self.solve(rings, piece, level))
            .collect()
    }

    /// Material intervals along the line at `level`, measured by `along`.
    pub(in crate::drawing) fn spans(
        &self,
        rings: &Rings,
        level: f64,
        along: impl Fn(P2) -> f64,
    ) -> Result<Vec<[f64; 2]>, GeometryError> {
        let mut stops: Vec<f64> = self
            .crossings(rings, level)
            .into_iter()
            .map(|hit| along(hit.point))
            .collect();
        if !stops.len().is_multiple_of(2) {
            return Err(err(
                GeometryErrorKind::InvalidResult,
                "A drawing profile is not a closed boundary",
            ));
        }
        stops.sort_by(f64::total_cmp);
        Ok(stops.as_chunks::<2>().0.to_vec())
    }
}

/// Even-odd containment by a ray along +x.
pub(in crate::drawing) fn contains(family: &Family, rings: &Rings, point: P2) -> bool {
    debug_assert!(family.normal == [0., 1.]);
    family
        .crossings(rings, point[1])
        .into_iter()
        .filter(|hit| hit.point[0] > point[0])
        .count()
        % 2
        == 1
}

/// Ascending disjoint union of intervals, joining those closer than `gap`.
pub(in crate::drawing) fn merge(mut spans: Vec<[f64; 2]>, gap: f64) -> Vec<[f64; 2]> {
    spans.retain(|[a, b]| b > a);
    spans.sort_by(|a, b| a[0].total_cmp(&b[0]));
    let mut merged: Vec<[f64; 2]> = Vec::with_capacity(spans.len());
    for [a, b] in spans {
        match merged.last_mut() {
            Some(last) if a <= last[1] + gap => last[1] = last[1].max(b),
            _ => merged.push([a, b]),
        }
    }
    merged
}

/// Intersection of two ascending disjoint interval lists.
pub(in crate::drawing) fn intersect(a: &[[f64; 2]], b: &[[f64; 2]]) -> Vec<[f64; 2]> {
    let mut output = Vec::new();
    for x in a {
        for y in b {
            let [lo, hi] = [x[0].max(y[0]), x[1].min(y[1])];
            if hi > lo {
                output.push([lo, hi]);
            }
        }
    }
    output
}

/// Stretches covered by exactly one of two ascending disjoint lists.
pub(in crate::drawing) fn exclusive(a: &[[f64; 2]], b: &[[f64; 2]], gap: f64) -> Vec<[f64; 2]> {
    let mut stops: Vec<(f64, usize)> = Vec::new();
    for (list, spans) in [a, b].into_iter().enumerate() {
        stops.extend(spans.iter().flatten().map(|&stop| (stop, list)));
    }
    stops.sort_by(|x, y| x.0.total_cmp(&y.0));
    let mut inside = [false; 2];
    let mut output = Vec::new();
    for pair in 0..stops.len() {
        inside[stops[pair].1] ^= true;
        if inside[0] != inside[1]
            && let Some(next) = stops.get(pair + 1)
            && next.0 - stops[pair].0 > gap
        {
            output.push([stops[pair].0, next.0]);
        }
    }
    output
}
