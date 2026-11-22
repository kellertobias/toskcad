//! An independent exact evaluator for axis-aligned box assemblies with
//! straight prism and cylinder cuts, used as a test oracle. Material is the
//! union of boxes minus removed boxes, optionally intersected with a crop
//! box, evaluated exactly on the grid of all box coordinates. Convex prisms
//! and cylinders are then subtracted in closed form against that grid:
//! polygon/rectangle clipping for prisms, circle/band segments for
//! cylinders. No geometry kernel is involved.
#![allow(dead_code)]

use multicad_core::{
    document::{Document, EntityId},
    model::{self, PrimitiveParameters},
    parameters::{self, Unit},
};
use multicad_precision::PortableMath as _;

pub type Aabb = [[f64; 3]; 2];

pub struct Grid {
    pub lines: [Vec<f64>; 3],
    solid: Vec<bool>,
}

fn inside(b: &Aabb, p: [f64; 3]) -> bool {
    (0..3).all(|i| b[0][i] < p[i] && p[i] < b[1][i])
}

impl Grid {
    pub fn new(material: &[Aabb], removed: &[Aabb], crop: Option<Aabb>) -> Self {
        let lines: [Vec<f64>; 3] = std::array::from_fn(|axis| {
            let mut v: Vec<f64> = material
                .iter()
                .chain(removed)
                .chain(crop.iter())
                .flat_map(|b| [b[0][axis], b[1][axis]])
                .collect();
            v.sort_by(f64::total_cmp);
            v.dedup();
            v
        });
        let n = lines.each_ref().map(|l| l.len() - 1);
        let mut solid = vec![false; n[0] * n[1] * n[2]];
        for i in 0..n[0] {
            for j in 0..n[1] {
                for k in 0..n[2] {
                    let c = [
                        (lines[0][i] + lines[0][i + 1]) / 2.,
                        (lines[1][j] + lines[1][j + 1]) / 2.,
                        (lines[2][k] + lines[2][k + 1]) / 2.,
                    ];
                    solid[(i * n[1] + j) * n[2] + k] = material.iter().any(|b| inside(b, c))
                        && !removed.iter().any(|b| inside(b, c))
                        && crop.is_none_or(|b| inside(&b, c));
                }
            }
        }
        Self { lines, solid }
    }
    fn n(&self) -> [usize; 3] {
        self.lines.each_ref().map(|l| l.len() - 1)
    }
    fn at(&self, c: [isize; 3]) -> bool {
        let n = self.n();
        if (0..3).any(|i| c[i] < 0 || c[i] as usize >= n[i]) {
            return false;
        }
        self.solid[(c[0] as usize * n[1] + c[1] as usize) * n[2] + c[2] as usize]
    }
    /// Every solid cell as its box.
    pub fn cells(&self) -> Vec<Aabb> {
        let n = self.n();
        let mut out = vec![];
        for i in 0..n[0] {
            for j in 0..n[1] {
                for k in 0..n[2] {
                    if self.at([i as isize, j as isize, k as isize]) {
                        let idx = [i, j, k];
                        out.push([
                            std::array::from_fn(|a| self.lines[a][idx[a]]),
                            std::array::from_fn(|a| self.lines[a][idx[a] + 1]),
                        ]);
                    }
                }
            }
        }
        out
    }
    /// Every boundary face between solid and empty cells: its axis, plane
    /// coordinate and rectangle.
    pub fn boundary(&self) -> Vec<(usize, f64, Aabb)> {
        let n = self.n();
        let mut out = vec![];
        for axis in 0..3 {
            for i in 0..n[0] as isize + 1 {
                for j in 0..n[1] as isize + 1 {
                    for k in 0..n[2] as isize + 1 {
                        let c = [i, j, k];
                        if (0..3).any(|a| a != axis && c[a] as usize >= n[a]) {
                            continue;
                        }
                        let mut below = c;
                        below[axis] -= 1;
                        if self.at(c) == self.at(below) {
                            continue;
                        }
                        let idx = c.map(|x| x as usize);
                        let face: Aabb = [
                            std::array::from_fn(|a| self.lines[a][idx[a]]),
                            std::array::from_fn(|a| {
                                if a == axis {
                                    self.lines[a][idx[a]]
                                } else {
                                    self.lines[a][idx[a] + 1]
                                }
                            }),
                        ];
                        out.push((axis, self.lines[axis][idx[axis]], face));
                    }
                }
            }
        }
        out
    }
    pub fn volume_mm3(&self) -> f64 {
        self.cells()
            .iter()
            .map(|b| (0..3).map(|i| b[1][i] - b[0][i]).product::<f64>())
            .sum()
    }
    pub fn area_mm2(&self) -> f64 {
        self.boundary()
            .iter()
            .map(|(axis, _, f)| {
                (0..3)
                    .filter(|a| a != axis)
                    .map(|a| f[1][a] - f[0][a])
                    .product::<f64>()
            })
            .sum()
    }
    pub fn bounds(&self) -> Aabb {
        let mut b = [[f64::INFINITY; 3], [f64::NEG_INFINITY; 3]];
        for c in self.cells() {
            for i in 0..3 {
                b[0][i] = b[0][i].min(c[0][i]);
                b[1][i] = b[1][i].max(c[1][i]);
            }
        }
        b
    }
}

/// A tool end may lie on a grid plane only away from material there, so the
/// strict cap and mouth tests below stay exact.
fn ends_clear(grid: &Grid, axis: usize, span: [f64; 2], bounds: &Aabb) {
    for end in span {
        if !grid.lines[axis].contains(&end) {
            continue;
        }
        for cell in grid.cells() {
            if cell[0][axis] == end || cell[1][axis] == end {
                assert!(
                    (0..3)
                        .any(|i| i != axis
                            && (cell[1][i] <= bounds[0][i] || bounds[1][i] <= cell[0][i])),
                    "a tool ends exactly on a material face"
                );
            }
        }
    }
}

fn others(axis: usize) -> [usize; 2] {
    match axis {
        0 => [1, 2],
        1 => [0, 2],
        _ => [0, 1],
    }
}
fn overlap(a: [f64; 2], b: [f64; 2]) -> f64 {
    (a[1].min(b[1]) - a[0].max(b[0])).max(0.)
}
fn shoelace(p: &[[f64; 2]]) -> f64 {
    (0..p.len())
        .map(|i| {
            let [a, b] = [p[i], p[(i + 1) % p.len()]];
            a[0] * b[1] - b[0] * a[1]
        })
        .sum::<f64>()
        .abs()
        / 2.
}
/// Sutherland–Hodgman clip of a polygon to a rectangle.
fn clip(polygon: &[[f64; 2]], rect: [[f64; 2]; 2]) -> Vec<[f64; 2]> {
    let mut out = polygon.to_vec();
    for (axis, bound, keep_below) in [
        (0, rect[0][0], false),
        (0, rect[1][0], true),
        (1, rect[0][1], false),
        (1, rect[1][1], true),
    ] {
        let input = std::mem::take(&mut out);
        let kept = |p: [f64; 2]| {
            if keep_below {
                p[axis] <= bound
            } else {
                p[axis] >= bound
            }
        };
        for i in 0..input.len() {
            let (a, b) = (input[i], input[(i + 1) % input.len()]);
            let crossing = || {
                let t = (bound - a[axis]) / (b[axis] - a[axis]);
                [a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]
            };
            match (kept(a), kept(b)) {
                (true, true) => out.push(b),
                (true, false) => out.push(crossing()),
                (false, true) => {
                    out.push(crossing());
                    out.push(b);
                }
                (false, false) => {}
            }
        }
        if out.is_empty() {
            break;
        }
    }
    out
}
/// Length of a segment inside a rectangle (Liang–Barsky).
fn clipped_length(a: [f64; 2], b: [f64; 2], rect: [[f64; 2]; 2]) -> f64 {
    let d = [b[0] - a[0], b[1] - a[1]];
    let (mut t0, mut t1) = (0_f64, 1_f64);
    for axis in 0..2 {
        if d[axis] == 0. {
            if a[axis] < rect[0][axis] || a[axis] > rect[1][axis] {
                return 0.;
            }
            continue;
        }
        let ta = (rect[0][axis] - a[axis]) / d[axis];
        let tb = (rect[1][axis] - a[axis]) / d[axis];
        t0 = t0.max(ta.min(tb));
        t1 = t1.min(ta.max(tb));
    }
    (t1 - t0).max(0.) * d[0].portable_hypot(d[1])
}

/// A convex polygon prism swept along one world axis.
pub struct Prism {
    pub axis: usize,
    pub span: [f64; 2],
    /// Profile in the two other axes, in increasing axis order.
    pub profile: Vec<[f64; 2]>,
}
/// A cylinder along one world axis.
pub struct Cylinder {
    pub axis: usize,
    pub span: [f64; 2],
    pub center: [f64; 2],
    pub radius: f64,
}

impl Prism {
    pub fn bounds(&self) -> Aabb {
        let [b, c] = others(self.axis);
        let mut out = [[0.; 3]; 2];
        out[0][self.axis] = self.span[0];
        out[1][self.axis] = self.span[1];
        for (k, i) in [b, c].into_iter().enumerate() {
            out[0][i] = self
                .profile
                .iter()
                .map(|p| p[k])
                .fold(f64::INFINITY, f64::min);
            out[1][i] = self
                .profile
                .iter()
                .map(|p| p[k])
                .fold(f64::NEG_INFINITY, f64::max);
        }
        out
    }
    fn section(&self, k: usize, value: f64) -> Option<[f64; 2]> {
        let mut hits = vec![];
        for i in 0..self.profile.len() {
            let (a, b) = (self.profile[i], self.profile[(i + 1) % self.profile.len()]);
            if (a[k] - value) * (b[k] - value) < 0. {
                let t = (value - a[k]) / (b[k] - a[k]);
                hits.push(a[1 - k] + t * (b[1 - k] - a[1 - k]));
            }
        }
        (hits.len() == 2).then(|| [hits[0].min(hits[1]), hits[0].max(hits[1])])
    }
    /// Volume and surface area change when this prism is cut from the grid.
    pub fn cut(&self, grid: &Grid) -> (f64, f64) {
        let [b, c] = others(self.axis);
        // The profile must be convex and in general position on the grid.
        for i in 0..self.profile.len() {
            let [p, q, r] = [0, 1, 2].map(|d| self.profile[(i + d) % self.profile.len()]);
            let turn = (q[0] - p[0]) * (r[1] - q[1]) - (q[1] - p[1]) * (r[0] - q[0]);
            assert!(
                turn.abs() > 1e-12
                    && turn.signum() == {
                        let [p, q, r] = [0, 1, 2].map(|d| self.profile[d]);
                        ((q[0] - p[0]) * (r[1] - q[1]) - (q[1] - p[1]) * (r[0] - q[0])).signum()
                    }
            );
            for (k, axis) in [b, c].into_iter().enumerate() {
                assert!(
                    !grid.lines[axis].contains(&p[k]),
                    "profile vertex on a grid line"
                );
            }
        }
        ends_clear(grid, self.axis, self.span, &self.bounds());
        let rect = |cell: &Aabb| [[cell[0][b], cell[0][c]], [cell[1][b], cell[1][c]]];
        let [mut volume, mut area] = [0., 0.];
        for cell in grid.cells() {
            let along = overlap([cell[0][self.axis], cell[1][self.axis]], self.span);
            let clipped = shoelace(&clip(&self.profile, rect(&cell)));
            volume -= clipped * along;
            // New groove walls.
            for i in 0..self.profile.len() {
                let (p, q) = (self.profile[i], self.profile[(i + 1) % self.profile.len()]);
                area += clipped_length(p, q, rect(&cell)) * along;
            }
            // End caps inside material.
            for end in self.span {
                if cell[0][self.axis] < end && end < cell[1][self.axis] {
                    area += clipped;
                }
            }
        }
        for (axis, value, face) in grid.boundary() {
            if axis == self.axis {
                if self.span[0] < value && value < self.span[1] {
                    area -= shoelace(&clip(&self.profile, rect(&face)));
                }
            } else {
                let k = usize::from(axis == c);
                let other = if k == 0 { c } else { b };
                if let Some(range) = self.section(k, value) {
                    area -= overlap(range, [face[0][other], face[1][other]])
                        * overlap([face[0][self.axis], face[1][self.axis]], self.span);
                }
            }
        }
        (volume, area)
    }
}

/// Area of a disc between two chords τ0 < τ1 (offsets from the center
/// along one axis), and the arc length inside that band.
fn band(radius: f64, t0: f64, t1: f64) -> (f64, f64) {
    let [t0, t1] = [t0, t1].map(|t| t.clamp(-radius, radius));
    let f = |t: f64| {
        radius * radius * (t / radius).portable_asin() + t * (radius * radius - t * t).sqrt()
    };
    (
        f(t1) - f(t0),
        2. * radius * ((t1 / radius).portable_asin() - (t0 / radius).portable_asin()),
    )
}

impl Cylinder {
    pub fn bounds(&self) -> Aabb {
        let [b, c] = others(self.axis);
        let mut out = [[0.; 3]; 2];
        out[0][self.axis] = self.span[0];
        out[1][self.axis] = self.span[1];
        for (k, i) in [b, c].into_iter().enumerate() {
            out[0][i] = self.center[k] - self.radius;
            out[1][i] = self.center[k] + self.radius;
        }
        out
    }
    /// The disc inside a rectangle that cuts it along at most one axis.
    fn in_rect(&self, rect: [[f64; 2]; 2]) -> (f64, f64) {
        let r = self.radius;
        let cuts: Vec<usize> = (0..2)
            .filter(|&k| rect[0][k] > self.center[k] - r || rect[1][k] < self.center[k] + r)
            .collect();
        match cuts.as_slice() {
            [] => (std::f64::consts::PI * r * r, 2. * std::f64::consts::PI * r),
            [k] => band(
                r,
                rect[0][*k] - self.center[*k],
                rect[1][*k] - self.center[*k],
            ),
            _ => panic!("a rectangle corner cuts the bore disc"),
        }
    }
    pub fn cut(&self, grid: &Grid) -> (f64, f64) {
        let (volume, removed, inside) = self.parts(grid);
        (volume, inside - removed)
    }
    /// Volume and surface area change when this cylinder is joined to the grid.
    pub fn join(&self, grid: &Grid) -> (f64, f64) {
        let (volume, removed, inside) = self.parts(grid);
        let r = self.radius;
        let length = self.span[1] - self.span[0];
        let pi = std::f64::consts::PI;
        (
            pi * r * r * length + volume,
            2. * pi * r * length + 2. * pi * r * r - removed - inside,
        )
    }
    /// The material volume inside the cylinder (negative), the material
    /// boundary area inside it, and the cylinder surface inside material.
    fn parts(&self, grid: &Grid) -> (f64, f64, f64) {
        let [b, c] = others(self.axis);
        ends_clear(grid, self.axis, self.span, &self.bounds());
        let rect = |cell: &Aabb| [[cell[0][b], cell[0][c]], [cell[1][b], cell[1][c]]];
        let [mut volume, mut inside, mut removed] = [0., 0., 0.];
        for cell in grid.cells() {
            let along = overlap([cell[0][self.axis], cell[1][self.axis]], self.span);
            let rect = rect(&cell);
            let disjoint = (0..2).any(|k| {
                rect[1][k] <= self.center[k] - self.radius
                    || rect[0][k] >= self.center[k] + self.radius
            });
            if disjoint {
                continue;
            }
            let (disc, arc) = self.in_rect(rect);
            volume -= disc * along;
            inside += arc * along;
            for end in self.span {
                if cell[0][self.axis] < end && end < cell[1][self.axis] {
                    inside += disc;
                }
            }
        }
        for (axis, value, face) in grid.boundary() {
            let rect = rect(&face);
            if axis == self.axis {
                let disjoint = (0..2).any(|k| {
                    rect[1][k] <= self.center[k] - self.radius
                        || rect[0][k] >= self.center[k] + self.radius
                });
                if self.span[0] < value && value < self.span[1] && !disjoint {
                    removed += self.in_rect(rect).0;
                }
            } else {
                let k = usize::from(axis == c);
                let d = value - self.center[k];
                if d.abs() < self.radius {
                    let half = (self.radius * self.radius - d * d).sqrt();
                    let chord = [self.center[1 - k] - half, self.center[1 - k] + half];
                    let extent = [rect[0][1 - k], rect[1][1 - k]];
                    removed += overlap(chord, extent)
                        * overlap([face[0][self.axis], face[1][self.axis]], self.span);
                }
            }
        }
        (volume, removed, inside)
    }
}

pub enum Tool {
    Box(Aabb),
    Prism(Prism),
    Cylinder(Cylinder),
}

/// The saved primitive or straight extrusion behind `feature`, in world
/// millimeters from its authored dimensions and frame.
pub fn operand(document: &Document, feature: EntityId) -> Tool {
    let values = parameters::evaluate(document).unwrap();
    let mm = |d: &multicad_core::constraint_sketch::Dimension| {
        d.value_in(Unit::Millimeter, Some(&values)).unwrap()
    };
    let axis_of = |d: [f64; 3]| {
        let axis = (0..3).find(|&i| d[i].abs() == 1.).unwrap();
        assert!(
            (0..3).all(|i| i == axis || d[i] == 0.),
            "axis-aligned frames only"
        );
        (axis, d[axis])
    };
    if let Ok(primitive) = model::read_primitive(document, feature) {
        let (origin, u, v, extents) = match &primitive.parameters {
            PrimitiveParameters::Box { size, frame } => (
                frame.origin.each_ref().map(mm),
                frame.u,
                frame.v,
                size.each_ref().map(mm),
            ),
            PrimitiveParameters::Cylinder {
                radius,
                height,
                frame,
            } => (
                frame.origin.each_ref().map(mm),
                frame.u,
                frame.v,
                [mm(radius), mm(radius), mm(height)],
            ),
            other => panic!("Unexpected primitive {other:?}"),
        };
        let n = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        let directions = [u, v, n].map(axis_of);
        let range = |(axis, sign): (usize, f64), length: f64| {
            let [a, b] = [origin[axis], origin[axis] + sign * length];
            (axis, [a.min(b), a.max(b)])
        };
        return match &primitive.parameters {
            PrimitiveParameters::Box { .. } => {
                let mut b = [[0.; 3]; 2];
                for (direction, length) in directions.into_iter().zip(extents) {
                    let (axis, [lo, hi]) = range(direction, length);
                    [b[0][axis], b[1][axis]] = [lo, hi];
                }
                Tool::Box(b)
            }
            _ => {
                let (axis, span) = range(directions[2], extents[2]);
                let [p, q] = [0, 1, 2]
                    .map(|i| i)
                    .into_iter()
                    .filter(|&i| i != axis)
                    .collect::<Vec<_>>()[..]
                else {
                    unreachable!()
                };
                Tool::Cylinder(Cylinder {
                    axis,
                    span,
                    center: [origin[p], origin[q]],
                    radius: extents[0],
                })
            }
        };
    }
    let prism = super::saved_prism::extruded_prism(document, feature);
    assert!(
        prism.arcs.iter().all(Option::is_none),
        "polygon prisms only"
    );
    let rectangle = prism.loop_mm.len() == 4
        && (0..4).all(|i| {
            let [a, b] = [prism.loop_mm[i], prism.loop_mm[(i + 1) % 4]];
            (0..3).filter(|&k| a[k] != b[k]).count() == 1
        });
    if rectangle {
        let mut b = [[f64::INFINITY; 3], [f64::NEG_INFINITY; 3]];
        for p in &prism.loop_mm {
            for q in [*p, std::array::from_fn(|i| p[i] + prism.sweep[i])] {
                for i in 0..3 {
                    b[0][i] = b[0][i].min(q[i]);
                    b[1][i] = b[1][i].max(q[i]);
                }
            }
        }
        return Tool::Box(b);
    }
    let (axis, _) = axis_of(prism.n);
    let others: Vec<usize> = (0..3).filter(|&i| i != axis).collect();
    let start = prism.loop_mm[0][axis];
    assert!(prism.loop_mm.iter().all(|p| p[axis] == start));
    let end = start + prism.sweep[axis];
    Tool::Prism(Prism {
        axis,
        span: [start.min(end), start.max(end)],
        profile: prism
            .loop_mm
            .iter()
            .map(|p| [p[others[0]], p[others[1]]])
            .collect(),
    })
}
