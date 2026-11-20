//! Certified minimum distance between a placed worm and generated wheel
//! that `clearance::certify` has proved disjoint.
//!
//! The wheel boundary consists of the generated space patches, the bore,
//! and the land and end-face material `M` on the convex blank `K`. Worm
//! points are never wheel material, because the solids are disjoint.
//!
//! * A worm point inside `K` is nearest to a patch or the bore: the segment
//!   from a land, end-face or rim point toward it stays in `K` and starts in
//!   material, so such a point is never the nearest.
//! * For a worm point `p` outside `K` with projection `p*` onto `K` not in
//!   `M`, a nearest point of `M` with `p - q` in its normal cone would
//!   project to itself, so it lies on an edge of `M`, which is a patch or
//!   bore edge, unless `p - q` points into `K` and the segment crosses the
//!   whole blank, at least its width `2h`. If `p*` is in `M`, the distance
//!   is at least `|p - p*|`.
//!
//! A box of worm points has a projection free of material when the
//! projection avoids every edge of `M` (patch exits and end edges), the
//! box avoids the patches and the bore, and the projected center lies in a
//! tooth space: the projection is connected and can only enter material
//! across those edges.
//!
//! Branch and bound covers the worm boundary with parameter boxes. Each box
//! keeps the wheel patch cells, over all sectors, that may still come
//! within the best distance. On a pair of a worm box `X` and a wheel cell
//! `Y` three bounds hold, each enclosed by outward-rounded intervals: the
//! gap between centered enclosures of `P(X)` and `S(Y)`; the mean value
//! form of `g = |P - S|²`, `g(c) - Σ ρᵢ max|∂g/∂xᵢ|`; and the offset along
//! the unit cell normal `n` at its center, `(P - S)·n ≥ f(c) - Σ ρᵢ
//! max|∂f/∂xᵢ|`. The last is second order with a curvature-sized constant
//! where the flanks face each other. Cell centers and box centers, and the
//! material projections of box centers, are actual boundary points, whose
//! distances bound the minimum from above.
use super::clearance::Placement;
use super::interval::Interval;
use super::jet::Real;
use super::wheel::{PIECES, Piece, Wheel};
use super::wheel_bounds::{Pair, cells};
use super::worm::Worm;
use super::*;
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::f64::consts::TAU;

const MAX_BOXES: usize = 2_000_000;
const MAX_PAIRS: usize = 1_000_000_000;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Error {
    Invalid,
    /// An enclosure could not be proved.
    Unresolved,
    Limit,
    Cancelled,
}

/// Closest worm and wheel points within a certified gap.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PairDistance {
    pub distance_mm: f64,
    pub lower_mm: f64,
    /// World points on the worm and on the wheel.
    pub points_mm: [Point; 2],
}

pub(super) type V3 = [Interval; 3];

pub(super) fn k(value: f64) -> Interval {
    Interval::point(value)
}

pub(super) fn mid(i: Interval) -> f64 {
    (i.lo + i.hi) / 2.
}

fn halves(i: Interval) -> [Interval; 2] {
    let m = mid(i);
    [Interval::new(i.lo, m), Interval::new(m, i.hi)]
}

fn gap(a: Interval, b: Interval) -> f64 {
    (b.lo - a.hi).max(a.lo - b.hi).max(0.)
}

pub(super) fn square(i: Interval) -> Interval {
    let m = i.lo.abs().min(i.hi.abs());
    let lo = if i.lo <= 0. && i.hi >= 0. { 0. } else { m * m };
    let top = i.lo.abs().max(i.hi.abs());
    Interval::new(lo.next_down().max(0.), (top * top).next_up())
}

/// Intersection of two enclosures of the same value.
fn meet(a: Interval, b: Interval) -> Interval {
    Interval::new(a.lo.max(b.lo), a.hi.min(b.hi))
}

/// Centered form `f(c) + f_p (p - c_p) + f_q (q - c_q)` met with a direct
/// enclosure of the same function over the box.
pub(super) fn centered(direct: V3, center: V3, d: [V3; 2], x: [Interval; 2]) -> V3 {
    let offset = x.map(|i| i - k(mid(i)));
    std::array::from_fn(|i| {
        meet(
            direct[i],
            center[i] + d[0][i] * offset[0] + d[1][i] * offset[1],
        )
    })
}

fn dot3(a: V3, b: V3) -> Interval {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Worm boundary patches in the thread frame (mm).
#[derive(Clone, Copy, Debug)]
pub(super) enum Patch {
    /// Profile band with radius `r0 + slope (s - s0)` over `(s, z)`.
    Band { r0: f64, s0: f64, slope: f64 },
    /// End cap at height `z` over `(p, q)`, outer radius `a + b q`.
    Cap { z: f64, a: f64, b: f64 },
    /// Bore wall over `(angle, z)`.
    Bore,
}

pub(super) struct WormMap<'a> {
    pub(super) worm: &'a Worm,
    rate: f64,
    bore: f64,
    /// Thread frame to wheel local module units: `M t + c`.
    matrix: [[f64; 3]; 3],
    offset: Point,
    pub(super) to_world: (Placement, Placement),
}

impl<'a> WormMap<'a> {
    pub(super) fn new(
        worm: &'a Worm,
        worm_at: Placement,
        wheel: &Wheel,
        wheel_at: Placement,
    ) -> Self {
        let m = wheel.module_mm();
        // Compose thread → worm local → world → wheel local, in module units.
        let image = |t: Point| {
            wheel_at
                .local(worm_at.world(worm.to_local(t)))
                .map(|v| v / m)
        };
        let offset = image([0.; 3]);
        let columns = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]].map(|e| sub(image(e), offset));
        Self {
            worm,
            rate: worm.thread().helical_rate(),
            bore: worm.bore_radius_mm().unwrap_or(0.),
            matrix: std::array::from_fn(|i| [columns[0][i], columns[1][i], columns[2][i]]),
            offset,
            to_world: (worm_at, wheel_at),
        }
    }

    /// Thread point and its two parameter derivatives over a box.
    fn thread(&self, patch: Patch, [p, q]: [Interval; 2]) -> [V3; 3] {
        match patch {
            Patch::Band { r0, s0, slope } => {
                let r = k(r0) + k(slope) * (p - k(s0));
                let theta = (q - p) / k(self.rate);
                let (s, c) = theta.sin_cos();
                let w = k(1. / self.rate);
                [
                    [r * c, r * s, q],
                    [k(slope) * c + r * s * w, k(slope) * s - r * c * w, k(0.)],
                    [-(r * s * w), r * c * w, k(1.)],
                ]
            }
            Patch::Cap { z, a, b } => {
                let outer = k(a) + k(b) * q;
                let reach = outer - k(self.bore);
                let rho = k(self.bore) + reach * p;
                let (s, c) = q.sin_cos();
                let drho = k(b) * p;
                [
                    [rho * c, rho * s, k(z)],
                    [reach * c, reach * s, k(0.)],
                    [drho * c - rho * s, drho * s + rho * c, k(0.)],
                ]
            }
            Patch::Bore => {
                let (s, c) = p.sin_cos();
                let r = k(self.bore);
                [
                    [r * c, r * s, q],
                    [-(r * s), r * c, k(0.)],
                    [k(0.), k(0.), k(1.)],
                ]
            }
        }
    }

    fn apply(&self, t: V3, affine: bool) -> V3 {
        std::array::from_fn(|i| {
            let row = self.matrix[i];
            let linear = k(row[0]) * t[0] + k(row[1]) * t[1] + k(row[2]) * t[2];
            if affine {
                linear + k(self.offset[i])
            } else {
                linear
            }
        })
    }

    /// Wheel-local value and derivatives over a box, in module units.
    pub(super) fn local(&self, patch: Patch, x: [Interval; 2]) -> [V3; 3] {
        let [t, tp, tq] = self.thread(patch, x);
        [
            self.apply(t, true),
            self.apply(tp, false),
            self.apply(tq, false),
        ]
    }

    pub(super) fn thread_point(&self, patch: Patch, x: [f64; 2]) -> Point {
        self.thread(patch, x.map(k))[0].map(mid)
    }

    pub(super) fn point(&self, patch: Patch, x: [f64; 2]) -> (Point, Point) {
        let [t, ..] = self.thread(patch, x.map(k));
        let t = t.map(mid);
        let local = self.apply(t.map(k), true).map(mid);
        let (worm_at, _) = self.to_world;
        (local, worm_at.world(self.worm.to_local(t)))
    }
}

/// Cached enclosures of one sector-0 patch cell, shared by every sector.
struct Cell {
    piece: Piece,
    u: Interval,
    v: Interval,
    /// Center patch point (sector 0, module units) and its rotation angle.
    point: Point,
    base: f64,
    /// Enclosures in the frame turned back by `base`: center value, value
    /// over the cell, and the two parameter derivatives.
    center: V3,
    value: V3,
    du: V3,
    dv: V3,
    /// Unit normal at the center, shrunk below unit length for rounding.
    normal: Point,
    children: Option<[usize; 2]>,
    /// Image extents along u and v, which choose the split.
    extent: [f64; 2],
}

struct Cells<'a> {
    wheel: &'a Wheel,
    all: Vec<Cell>,
}

impl Cells<'_> {
    fn make(&mut self, piece: Piece, u: Interval, v: Interval) -> Result<usize, Error> {
        let wheel = self.wheel;
        let mid_uv = [mid(u), mid(v)];
        let crossing = piece
            .flank_side()
            .map_or(0., |side| wheel.crossing(mid_uv[1], side));
        let (q, base) = wheel.polar_patch(piece, mid_uv[0], mid_uv[1], crossing);
        let [x, y] = super::wheel::turn(q, base);
        let enclose = |u: Interval, v: Interval| -> Option<[Pair<2, 2>; 2]> {
            Some(wheel.coefficients(piece, [u, v], wheel.exit(piece, v)?, base))
        };
        let center = enclose(k(mid_uv[0]), k(mid_uv[1])).ok_or(Error::Unresolved)?;
        let over = enclose(u, v).ok_or(Error::Unresolved)?;
        let value = |p: &Pair<2, 2>| p.0[0].0[0];
        let du = [over[0].0[1].0[0], over[1].0[1].0[0], k(0.)];
        let dv = [over[0].0[0].0[1], over[1].0[0].0[1], k(1.)];
        let length = |d: V3| d.iter().map(|c| c.magnitude().powi(2)).sum::<f64>().sqrt();
        let extent = [length(du) * (u.hi - u.lo), length(dv) * (v.hi - v.lo)];
        let tangent = |p: &Pair<2, 2>, along_u: bool| {
            if along_u {
                mid(p.0[1].0[0])
            } else {
                mid(p.0[0].0[1])
            }
        };
        let su = [tangent(&center[0], true), tangent(&center[1], true), 0.];
        let sv = [tangent(&center[0], false), tangent(&center[1], false), 1.];
        let normal = cross(su, sv);
        let size = norm(normal) * (1. + 1e-12);
        let normal = normal.map(|x| x / size);
        self.all.push(Cell {
            normal,
            piece,
            u,
            v,
            point: [x, y, mid_uv[1]],
            base,
            center: [value(&center[0]), value(&center[1]), k(mid_uv[1])],
            value: centered(
                [value(&over[0]), value(&over[1]), v],
                [value(&center[0]), value(&center[1]), k(mid_uv[1])],
                [du, dv],
                [u, v],
            ),
            du,
            dv,
            children: None,
            extent,
        });
        Ok(self.all.len() - 1)
    }

    fn children(&mut self, index: usize) -> Result<[usize; 2], Error> {
        if let Some(children) = self.all[index].children {
            return Ok(children);
        }
        let Cell {
            piece,
            u,
            v,
            extent,
            ..
        } = self.all[index];
        let pair = if extent[0] >= extent[1] {
            let [a, b] = halves(u);
            [self.make(piece, a, v)?, self.make(piece, b, v)?]
        } else {
            let [a, b] = halves(v);
            [self.make(piece, u, a)?, self.make(piece, u, b)?]
        };
        self.all[index].children = Some(pair);
        Ok(pair)
    }
}

/// What a target bounds: a patch cell measured directly, or a cell on an
/// edge of the land and end-face material that the blank projection of a
/// box must avoid.
#[derive(Clone, Copy, PartialEq)]
enum Kind {
    Patch,
    Edge,
}

/// One wheel cell of one sector still open for a worm box.
#[derive(Clone, Copy)]
struct Target {
    cell: usize,
    sector: u32,
    kind: Kind,
}

struct WormBox {
    patch: Patch,
    x: [Interval; 2],
    lower: f64,
    /// Proved free of material after projection.
    open: bool,
    /// Every projected point is at least the best distance from the blank.
    settled: bool,
    targets: Vec<Target>,
    /// Worm-side spreads of the closest target, or image extents.
    spread: [f64; 2],
}

impl PartialEq for WormBox {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}
impl Eq for WormBox {}
impl PartialOrd for WormBox {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for WormBox {
    fn cmp(&self, other: &Self) -> Ordering {
        other.lower.total_cmp(&self.lower)
    }
}

/// One smooth map of the worm box toward the blank: the identity, the end
/// clamp `z ↦ ±h`, the radial scaling `(x, y) ↦ B (x, y) / r`, or both.
/// `value` encloses the map over the box points it projects; the center
/// and derivatives hold for the map formula over the whole box.
struct Map {
    value: V3,
    center: V3,
    dp: V3,
    dq: V3,
}

/// Evaluated worm box.
struct Image {
    /// The box lies inside the closed blank.
    inside: bool,
    identity: Map,
    /// Blank projections of the box points outside the blank.
    maps: Vec<Map>,
    /// Distance from the box to the blank, module units.
    outside: f64,
    radius: [f64; 2],
    x: [Interval; 2],
    /// Image extents along each worm parameter.
    extent: [f64; 2],
    /// Actual center point, wheel local module units and world mm.
    local: Point,
    world: Point,
}

struct Search<'a> {
    map: WormMap<'a>,
    cells: Cells<'a>,
    wheel: &'a Wheel,
    blank: f64,
    h: f64,
    /// Rounding slack of the composed placement map, module units.
    slack: f64,
    best: f64,
    witness: Option<[Point; 2]>,
    pairs: usize,
}

impl Search<'_> {
    fn image(&self, patch: Patch, x: [Interval; 2]) -> Image {
        let [value, dp, dq] = self.map.local(patch, x);
        let [center, ..] = self.map.local(patch, x.map(|i| k(mid(i))));
        let value = centered(value, center, [dp, dq], x);
        let (local, world) = self.map.point(patch, x.map(mid));
        let radius = (square(value[0]) + square(value[1])).sqrt();
        let z = value[2];
        let (blank, h) = (self.blank, self.h);
        let beyond = (z.lo - h).max(-h - z.hi).max(0.);
        let outside = (radius.lo - blank).max(0.).hypot(beyond);
        let within = radius.lo <= blank;
        let beyond_radius = radius.hi > blank;
        let slab = z.lo <= h && z.hi >= -h;
        let identity = Map {
            value,
            center,
            dp,
            dq,
        };
        let flat = |m: &Map, plane: f64| Map {
            value: [m.value[0], m.value[1], k(plane)],
            center: [m.center[0], m.center[1], k(plane)],
            dp: [m.dp[0], m.dp[1], k(0.)],
            dq: [m.dq[0], m.dq[1], k(0.)],
        };
        let radial = |m: &Map| {
            // Over the projected points r ≥ B; derivatives hold for all r.
            let clipped = Interval::new(radius.lo.max(blank), radius.hi.max(blank));
            let unit = |r: Interval| [value[0] / r, value[1] / r];
            let [ex, ey] = unit(radius);
            let [cx, cy] = unit(clipped);
            let rc = (square(m.center[0]) + square(m.center[1])).sqrt();
            let scale = k(blank) / radius;
            let turn = |d: V3| {
                let along = ex * d[0] + ey * d[1];
                [
                    scale * (d[0] - ex * along),
                    scale * (d[1] - ey * along),
                    d[2],
                ]
            };
            Map {
                value: [k(blank) * cx, k(blank) * cy, m.value[2]],
                center: [
                    k(blank) * m.center[0] / rc,
                    k(blank) * m.center[1] / rc,
                    m.center[2],
                ],
                dp: turn(m.dp),
                dq: turn(m.dq),
            }
        };
        let inside = radius.hi <= blank && z.lo >= -h && z.hi <= h;
        let mut maps = vec![];
        if !inside {
            let radial_map = beyond_radius.then(|| radial(&identity));
            for (plane, reach) in [(h, z.hi > h), (-h, z.lo < -h)] {
                if reach && within {
                    maps.push(flat(&identity, plane));
                }
                if let (true, Some(r)) = (reach, &radial_map) {
                    maps.push(flat(r, plane));
                }
            }
            if let (true, Some(r)) = (slab, radial_map) {
                maps.push(r);
            }
        }
        let length = |d: V3| d.iter().map(|c| c.magnitude().powi(2)).sum::<f64>().sqrt();
        let radius = x.map(|i| (i.hi - i.lo) / 2.);
        Image {
            inside,
            identity,
            outside,
            extent: [length(dp) * radius[0], length(dq) * radius[1]].map(|e| 2. * e),
            maps,
            radius,
            x,
            local,
            world,
        }
    }

    /// Lower bound of the gap between a worm box image and one sector
    /// cell, with the worm-side and wheel-side spreads of that bound.
    fn bound(
        &mut self,
        image: &Image,
        cell: usize,
        sector: u32,
        projected: bool,
    ) -> (f64, [f64; 2], [f64; 2]) {
        self.pairs += 1;
        let c = &self.cells.all[cell];
        let angle = self.wheel.sector_angle(sector) + c.base;
        let (s, co) = k(-angle).sin_cos();
        let turn = |p: V3| [co * p[0] - s * p[1], s * p[0] + co * p[1], p[2]];
        let offsets = [image.x[0], image.x[1], c.u, c.v].map(|i| i - k(mid(i)));
        let mut result = (f64::INFINITY, [0.; 2], [0.; 2]);
        let maps = if projected {
            &image.maps[..]
        } else {
            std::slice::from_ref(&image.identity)
        };
        for map in maps {
            let value = turn(map.value);
            let center = turn(map.center);
            let (dp, dq) = (turn(map.dp), turn(map.dq));
            let diff: V3 = std::array::from_fn(|i| center[i] - c.center[i]);
            // Joint centered form of P - S, met with the enclosure difference.
            let d: V3 = std::array::from_fn(|i| {
                meet(
                    value[i] - c.value[i],
                    diff[i] + dp[i] * offsets[0] + dq[i] * offsets[1]
                        - c.du[i] * offsets[2]
                        - c.dv[i] * offsets[3],
                )
            });
            let zero = d
                .iter()
                .map(|i| gap(*i, k(0.)).powi(2))
                .sum::<f64>()
                .next_down()
                .max(0.)
                .sqrt();
            let g = dot3(diff, diff);
            let slope = |e: V3| (dot3(d, e).magnitude() * 2.).next_up();
            let worm = [slope(dp) * image.radius[0], slope(dq) * image.radius[1]];
            let wheel = [
                slope(c.du) * (c.u.hi - c.u.lo) / 2.,
                slope(c.dv) * (c.v.hi - c.v.lo) / 2.,
            ];
            let mean = (g.lo - (worm[0] + worm[1] + wheel[0] + wheel[1]).next_up())
                .max(0.)
                .sqrt();
            // Separating plane through the cell center normal: the signed
            // offset along a unit vector bounds the distance from below.
            let n = c.normal;
            let sign = if dot(n, diff.map(mid)) < 0. { -1. } else { 1. };
            let n = n.map(|x| k(sign * x));
            let along = |e: V3| dot3(e, n).magnitude().next_up();
            let tilt_worm = [along(dp) * image.radius[0], along(dq) * image.radius[1]];
            let tilt_wheel = [
                along(c.du) * (c.u.hi - c.u.lo) / 2.,
                along(c.dv) * (c.v.hi - c.v.lo) / 2.,
            ];
            let plane = dot3(diff, n).lo
                - (tilt_worm[0] + tilt_worm[1] + tilt_wheel[0] + tilt_wheel[1]).next_up();
            let lower = zero.max(mean).max(plane).max(0.);
            let lower = if lower.is_finite() { lower } else { 0. };
            if lower < result.0 {
                result = if plane >= mean {
                    (lower, tilt_worm, tilt_wheel)
                } else {
                    (lower, worm, wheel)
                };
            }
        }
        if result.0.is_infinite() {
            result.0 = 0.;
        }
        result
    }

    fn observe(&mut self, image: &Image, cell: usize, sector: u32) {
        let c = &self.cells.all[cell];
        let point = rotate_z(c.point, self.wheel.sector_angle(sector));
        let distance = norm(sub(point, image.local));
        if distance < self.best {
            self.best = distance;
            let m = self.wheel.module_mm();
            let (_, wheel_at) = self.map.to_world;
            self.witness = Some([image.world, wheel_at.world(point.map(|v| v * m))]);
        }
    }

    /// Lower bound of the box distance in module units, refining its
    /// targets, and the worm-side spread of its closest target.
    fn evaluate(
        &mut self,
        patch: Patch,
        x: [Interval; 2],
        inherited: &[Target],
        (parent_open, parent_settled): (bool, bool),
        cut: f64,
    ) -> Result<WormBox, Error> {
        let image = self.image(patch, x);
        let settled = parent_settled || image.outside - self.slack >= self.best - cut;
        let projected = !image.inside && !settled && !parent_open;
        let mut open = Vec::with_capacity(inherited.len());
        let mut pending: Vec<Target> = inherited
            .iter()
            .filter(|t| t.kind == Kind::Patch || projected)
            .copied()
            .collect();
        let mut spread = [0.; 2];
        let mut direct = f64::INFINITY;
        let mut edge = f64::INFINITY;
        while let Some(target) = pending.pop() {
            if self.pairs > MAX_PAIRS {
                return Err(Error::Limit);
            }
            let through = target.kind == Kind::Edge;
            let (lower, worm, wheel) = self.bound(&image, target.cell, target.sector, through);
            if !through {
                self.observe(&image, target.cell, target.sector);
            }
            // An edge only has to be avoided by the projection.
            if through && lower > 0. || !through && lower - self.slack >= self.best - cut {
                continue;
            }
            let c = &self.cells.all[target.cell];
            let small = c.u.hi - c.u.lo < 1e-9 && c.v.hi - c.v.lo < 1e-9;
            if !small && wheel[0] + wheel[1] > worm[0] + worm[1] {
                for child in self.cells.children(target.cell)? {
                    pending.push(Target {
                        cell: child,
                        ..target
                    });
                }
                continue;
            }
            let slot = if through { &mut edge } else { &mut direct };
            if lower < *slot {
                *slot = lower;
            }
            if !through && lower <= direct {
                spread = worm;
            }
            open.push(target);
        }
        if !image.inside {
            self.observe_material(&image);
        }
        let (bore, bore_edge) = self.bore_gaps(&image);
        let direct = direct.min(bore);
        let edge = edge.min(bore_edge);
        let free = image.inside
            || settled
            || parent_open
            || (direct > 0. && edge > 0. && self.projected_space(&image));
        // With the projection free of material, a nearest land or end-face
        // point is an edge point, a patch point, unless the segment to it
        // crosses the whole blank, at least its width 2h.
        let material = if image.inside || settled {
            f64::INFINITY
        } else if free {
            2. * self.h
        } else {
            image.outside
        };
        let lower = direct.min(material).min(self.best) - self.slack;
        let spread = if spread[0] + spread[1] > 0. {
            spread
        } else {
            image.extent
        };
        Ok(WormBox {
            patch,
            x,
            lower,
            open: free && !image.inside,
            settled,
            targets: if free {
                open.into_iter().filter(|t| t.kind == Kind::Patch).collect()
            } else {
                open
            },
            spread,
        })
    }

    /// Lower bounds to the wheel bore wall for the box, and to the bore end
    /// circles for its projected points.
    fn bore_gaps(&self, image: &Image) -> (f64, f64) {
        let Some(bore) = self.wheel.bore_radius_mm() else {
            return (f64::INFINITY, f64::INFINITY);
        };
        let bore = k(bore / self.wheel.module_mm());
        let radius = |m: &Map| (square(m.value[0]) + square(m.value[1])).sqrt();
        let wall = gap(radius(&image.identity), bore);
        let circle = image
            .maps
            .iter()
            .map(|m| {
                let ends = gap(m.value[2], k(self.h)).min(gap(m.value[2], k(-self.h)));
                gap(radius(m), bore).hypot(ends)
            })
            .fold(f64::INFINITY, f64::min);
        (wall, circle)
    }

    /// The blank projection of the box center, when it is land or end-face
    /// material, is an actual wheel point.
    fn observe_material(&mut self, image: &Image) {
        let p = image.local;
        let r = p[0].hypot(p[1]);
        let scale = if r > self.blank { self.blank / r } else { 1. };
        let q = [p[0] * scale, p[1] * scale, p[2].clamp(-self.h, self.h)];
        let wheel = self.wheel;
        let m = wheel.module_mm();
        let rq = q[0].hypot(q[1]);
        if wheel.bore_radius_mm().is_some_and(|b| rq * m <= b)
            || wheel.space_margin(q[0], q[1], q[2]) <= 0.
        {
            return;
        }
        let distance = norm(sub(p, q));
        if distance < self.best {
            self.best = distance;
            let (_, wheel_at) = self.map.to_world;
            self.witness = Some([image.world, wheel_at.world(q.map(|v| v * m))]);
        }
    }

    /// Whether the projected box center lies strictly in a tooth space.
    fn projected_space(&self, image: &Image) -> bool {
        let [x, y, z] = image.local;
        let r = x.hypot(y);
        let scale = if r > self.blank { self.blank / r } else { 1. };
        let (x, y, z) = (x * scale, y * scale, z.clamp(-self.h, self.h));
        let wheel = self.wheel;
        let m = wheel.module_mm();
        if wheel.bore_radius_mm().is_some_and(|b| (x.hypot(y)) * m < b) {
            return true;
        }
        wheel.space_margin(x, y, z) < -1e-12
    }
}

/// Worm boundary patches and their parameter boxes.
pub(super) fn worm_patches(worm: &Worm) -> Vec<(Patch, [Interval; 2])> {
    let thread = worm.thread();
    let length = worm.length_mm();
    let phases = worm.rail_phases();
    let mut patches = vec![];
    for i in 0..4 {
        let (lo, hi) = (phases[i + 1], phases[i]);
        let s0 = (lo + hi) / 2.;
        let [r0, slope] = thread.radius_and_slope_at_phase(s0);
        patches.push((
            Patch::Band { r0, s0, slope },
            [Interval::new(lo, hi), Interval::new(0., length)],
        ));
    }
    let rate = thread.helical_rate();
    let pitch = thread.spec().pitch_mm;
    let breaks = thread.profile_breaks();
    for cap in 0..2 {
        let z = cap as f64 * length;
        let start = z.rem_euclid(pitch);
        // phase(q) = start - rate q falls by one pitch over a turn.
        let mut cuts: Vec<f64> = vec![0., TAU];
        for b in breaks {
            for n in -2..=2 {
                let q = (start - b - f64::from(n) * pitch) / rate;
                if q > 0. && q < TAU {
                    cuts.push(q);
                }
            }
        }
        cuts.sort_by(f64::total_cmp);
        cuts.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
        for pair in cuts.windows(2) {
            let q = (pair[0] + pair[1]) / 2.;
            let [r, slope] = thread.radius_and_slope_at_phase(start - rate * q);
            // r(q') = r + slope (phase(q') - phase(q)) = r - slope rate (q' - q).
            let b = -slope * rate;
            patches.push((
                Patch::Cap { z, a: r - b * q, b },
                [Interval::new(0., 1.), Interval::new(pair[0], pair[1])],
            ));
        }
    }
    if worm.bore_radius_mm().is_some() {
        patches.push((
            Patch::Bore,
            [Interval::new(0., TAU), Interval::new(0., length)],
        ));
    }
    patches
}

/// Certified minimum distance of a placed worm and wheel whose clearance
/// certificate is `Clear`. The returned distance is attained by actual
/// boundary points and exceeds the proved lower bound by at most `gap_mm`.
pub fn minimum_distance(
    worm: &Worm,
    worm_at: Placement,
    wheel: &Wheel,
    wheel_at: Placement,
    gap_mm: f64,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<PairDistance, Error> {
    if !worm_at.valid() || !wheel_at.valid() || !(gap_mm.is_finite() && gap_mm > 0.) {
        return Err(Error::Invalid);
    }
    let m = wheel.module_mm();
    let map = WormMap::new(worm, worm_at, wheel, wheel_at);
    let scale = norm(map.offset) + wheel.blank() + worm.length_mm() / m;
    let mut search = Search {
        map,
        cells: Cells { wheel, all: vec![] },
        wheel,
        blank: wheel.blank(),
        h: wheel.half_width(),
        slack: 1e-12 * scale,
        best: f64::INFINITY,
        witness: None,
        pairs: 0,
    };
    let h = wheel.half_width();
    let mut roots = vec![];
    for piece in PIECES {
        for v in cells(-h, h, 4) {
            for u in cells(0., 1., 2) {
                roots.push((search.cells.make(piece, u, v)?, Kind::Patch));
            }
        }
        // Edges on the blank: the flank exits and every piece at both ends.
        if let Some(side) = piece.flank_side() {
            let u = k(if side > 0. { 0. } else { 1. });
            for v in cells(-h, h, 4) {
                roots.push((search.cells.make(piece, u, v)?, Kind::Edge));
            }
        }
        for v in [-h, h] {
            for u in cells(0., 1., 2) {
                roots.push((search.cells.make(piece, u, k(v))?, Kind::Edge));
            }
        }
    }
    let mut targets = vec![];
    for sector in 0..wheel.teeth() {
        for &(cell, kind) in &roots {
            targets.push(Target { cell, sector, kind });
        }
    }
    let cut = gap_mm / m;
    let mut heap = BinaryHeap::new();
    for (patch, x) in worm_patches(worm) {
        heap.push(search.evaluate(patch, x, &targets, (false, false), cut)?);
    }
    let mut boxes = 0;
    while let Some(b) = heap.pop() {
        if b.lower >= search.best - cut {
            let lower = b.lower.min(search.best);
            let points = search.witness.ok_or(Error::Unresolved)?;
            return Ok(PairDistance {
                distance_mm: norm(sub(points[0], points[1])),
                lower_mm: lower * m,
                points_mm: points,
            });
        }
        boxes += 1;
        if boxes > MAX_BOXES {
            return Err(Error::Limit);
        }
        if boxes % 256 == 0 && interrupt.is_some_and(|stop| stop()) {
            return Err(Error::Cancelled);
        }
        let axis = usize::from(b.spread[1] > b.spread[0]);
        for half in halves(b.x[axis]) {
            let mut x = b.x;
            x[axis] = half;
            let child = search.evaluate(b.patch, x, &b.targets, (b.open, b.settled), cut)?;
            if child.lower < search.best - cut {
                heap.push(child);
            }
        }
    }
    let points = search.witness.ok_or(Error::Unresolved)?;
    Ok(PairDistance {
        distance_mm: norm(sub(points[0], points[1])),
        lower_mm: search.best * m - gap_mm,
        points_mm: points,
    })
}

#[cfg(test)]
mod tests {
    use super::super::clearance::{ResultKind, certify};
    use super::super::{wheel::tests::spec as wheel_spec, worm::tests::spec as worm_spec};
    use super::*;

    /// Canonical mating pose, as in the clearance tests.
    fn pose(teeth: u32, drive: f64, extra: f64) -> [Placement; 2] {
        let center = f64::from(teeth) / 2. + 5.;
        let (sw, cw) = drive.sin_cos();
        let turn = -drive / f64::from(teeth) + extra;
        let (s, c) = turn.sin_cos();
        [
            Placement {
                origin_mm: [center, 0., 0.],
                axes: [[cw, 0., -sw], [-sw, 0., -cw], [0., 1., 0.]],
            },
            Placement {
                origin_mm: [0.; 3],
                axes: [[c, s, 0.], [-s, c, 0.], [0., 0., 1.]],
            },
        ]
    }

    #[test]
    fn clear_mating_pairs_measure_a_certified_flank_gap() {
        for (teeth, drive) in [(32, 0.), (40, 1.3)] {
            let worm = Worm::new(worm_spec(1., 16., 1., 0.05), 1e-3).unwrap();
            let wheel = Wheel::new(wheel_spec(1., teeth, 2., 2., 0.05), 1e-3, 1e7).unwrap();
            let [a, b] = pose(teeth, drive, 0.);
            let clear = certify(&worm, a, &wheel, b, 1e-7, 0.2, None).unwrap();
            assert_eq!(clear.result, ResultKind::Clear);
            let found = minimum_distance(&worm, a, &wheel, b, 1e-7, None).unwrap();
            assert!(found.lower_mm <= found.distance_mm);
            assert!(found.distance_mm - found.lower_mm <= 1e-7 + 1e-12);
            let measured = norm(sub(found.points_mm[0], found.points_mm[1]));
            assert!((measured - found.distance_mm).abs() < 1e-12);
            // The clearance witness is an actual pair, so never closer.
            assert!(found.distance_mm <= clear.nearest.unwrap().distance_mm + 1e-12);
            // Both witnesses lie on the boundaries.
            assert!(worm.margin(a.local(found.points_mm[0])).abs() < 1e-9);
            assert!(wheel.margin(b.local(found.points_mm[1])).abs() < 1e-9);
        }
    }
}
