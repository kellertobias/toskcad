//! Regularized subtraction of one finite circular cylinder from a solid
//! stepped shaft, with the tool axis perpendicular to the shaft axis: a
//! radial or offset cross-bore, through or blind.
//!
//! The shaft is a full revolution whose meridian has only axis-parallel and
//! axis-perpendicular edges. The tool's circular section must lie within one
//! cylindrical wall of the shaft, both across the axis and along it. The
//! result keeps its actual supports: every face of the shaft, the tool
//! cylinder and the planar tool caps. The bored wall meets the tool wall in
//! the exact quartic curve `x = ±sqrt(r² - y²)` over the tool angle, and a
//! blind end meets it in a circle. No mesh, polygonal cylinder or fitted
//! curve defines the boundary.
//!
//! Every admitted classification clears `RESOLUTION` times the largest
//! source magnitude. Operands that are the literal source values decide
//! true contact exactly in rational arithmetic, so a tool that only touches
//! the shaft leaves it unchanged. Tangency of the tool wall to the shaft
//! wall, a cap that crosses the shaft wall, a tool that crosses a shoulder
//! and an enclosed void are refused as unsupported; anything closer to a
//! decision boundary than the margin is unresolved.
//!
//! Volume, face areas and centroids are periodic integrals over the tool
//! angle of functions analytic in a strip of known width. The trapezoid rule
//! is used with the strip error bound of Trefethen and Weideman (SIAM Review
//! 56, 2014, Theorem 3.2), so each value carries a certified error bound.
use crate::display::{
    DisplayMesh, MeshError, Tessellation, triangulate_region, triangulate_region_with_interior,
};
use crate::revolution::Revolution;
use crate::revolution_mesh;
use num_rational::BigRational;
use spade::Point2;
use std::{
    cmp::Ordering,
    f64::consts::{PI, TAU},
};

/// Relative margin that every admitted classification must clear.
pub const RESOLUTION: f64 = 1. / (1_u64 << 30) as f64;
/// Largest trapezoid node count; tighter strips are unresolved.
const MAXIMUM_NODES: usize = 1 << 21;
const MAXIMUM_MM: f64 = 1e7;

type V = [f64; 3];
/// Appends one triangle of a face, oriented by its outward normal.
type Emit = dyn Fn(&mut DisplayMesh, [usize; 3], u32, V);

fn dot(a: V, b: V) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn mix(a: V, scale: f64, b: V) -> V {
    std::array::from_fn(|i| a[i] + scale * b[i])
}

fn cross(a: V, b: V) -> V {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn exact(value: f64) -> BigRational {
    BigRational::from_float(value).expect("validated finite operand")
}

/// Finite tool cylinder in the shaft chart, whose `z` axis is the shaft axis.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Tool {
    pub base_mm: V,
    /// Unit travel direction, perpendicular to the shaft axis.
    pub axis: V,
    pub radius_mm: f64,
    pub height_mm: f64,
    /// The operands are the exact source values, so true contact can be
    /// decided exactly; otherwise every decision must clear the margin.
    pub literal: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    Invalid,
    Unsupported,
    Unresolved,
    Cancelled,
}

/// How the tool leaves the shaft on one side of the axis.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum End {
    /// The tool passes the shaft wall: the bored wall carries a hole.
    Through,
    /// The planar tool cap lies inside the shaft: a blind bottom.
    Cap,
}

#[derive(Debug)]
pub enum Cut {
    /// The tool misses the material or only touches it.
    Unchanged,
    Body(Box<CrossBore>),
}

/// Area and first moment of one face in the shaft chart.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Moments {
    pub area_mm2: f64,
    pub centroid_mm: V,
}

/// A stepped shaft with one cross-bore.
#[derive(Clone, Debug)]
pub struct CrossBore {
    shaft: Revolution,
    /// Meridian edge of the bored cylindrical wall.
    wall: usize,
    radius: f64,
    span: [f64; 2],
    /// Tool travel `x'` and `y' = z × x'`, both perpendicular to the axis.
    across: [V; 2],
    /// Signed distance of the tool axis from the shaft axis along `y'`.
    offset: f64,
    /// Axial level of the tool axis.
    level: f64,
    rho: f64,
    /// Start and end of the tool along `x'`.
    stations: [f64; 2],
    ends: [End; 2],
    resolution: f64,
    nodes: usize,
    removed: f64,
    removed_error: f64,
    holes: [Option<Moments>; 2],
    tool_wall: Moments,
    /// Largest second derivative of the bored-wall curves over the angle.
    bend: f64,
}

/// Exact sign of `value` decided by margin or, for literal operands, by the
/// rational evaluation `exact`.
fn decide(
    value: f64,
    margin: f64,
    literal: bool,
    exact: impl FnOnce() -> BigRational,
) -> Result<Ordering, Error> {
    if value.abs() > margin {
        return Ok(value.total_cmp(&0.));
    }
    if literal {
        return Ok(exact().cmp(&BigRational::from_integer(0.into())));
    }
    Err(Error::Unresolved)
}

/// The cylindrical walls `(edge, radius, [low, high])` of a full, solid,
/// stepped meridian: one axis edge and a staircase that rises in `z`.
pub fn stepped_walls(shaft: &Revolution) -> Option<Vec<(usize, f64, [f64; 2])>> {
    let [ring] = shaft.rings() else {
        return None;
    };
    if !shaft.is_full() {
        return None;
    }
    let count = ring.len();
    let edge = |i: usize| [ring[i], ring[(i + 1) % count]];
    let axis: Vec<usize> = (0..count)
        .filter(|&i| {
            let [a, b] = edge(i);
            a[0] == 0. && b[0] == 0.
        })
        .collect();
    let [axis] = axis.as_slice() else {
        return None;
    };
    let mut walls = Vec::new();
    for step in 1..count {
        let i = (axis + step) % count;
        let [a, b] = edge(i);
        if a[0] == b[0] {
            if b[1] <= a[1] || a[0] <= 0. {
                return None;
            }
            walls.push((i, a[0], [a[1], b[1]]));
        } else if a[1] != b[1] {
            return None;
        }
    }
    let [low, high] = [edge(*axis)[1][1], edge(*axis)[0][1]];
    let spans: Vec<[f64; 2]> = walls.iter().map(|wall| wall.2).collect();
    (spans.first()?[0] == low
        && spans.last()?[1] == high
        && spans.windows(2).all(|pair| pair[0][1] == pair[1][0]))
    .then_some(walls)
}

/// Bounds the trapezoid sums of one periodic integrand.
struct Strip {
    /// Half width of the analytic strip used by the bound.
    width: f64,
    /// Bound of `|r² - y²|` from above and below over the strip.
    upper: f64,
    lower: f64,
    /// Bound of `|y|` over the strip.
    reach: f64,
}

impl Strip {
    /// Error of an `n`-node trapezoid sum of a function bounded by `bound`.
    fn error(&self, n: usize, bound: f64) -> f64 {
        4. * PI * bound / (self.width * n as f64).exp_m1()
    }

    fn nodes(&self, bound: f64, target: f64) -> Option<usize> {
        let needed = ((4. * PI * bound / target).ln_1p() / self.width).ceil();
        (needed.is_finite() && needed <= MAXIMUM_NODES as f64)
            .then(|| (needed as usize).max(16).next_multiple_of(4))
    }
}

impl CrossBore {
    /// Point `x' * across[0] + y' * across[1] + z * ẑ` of the bore chart.
    fn chart(&self, [x, y, z]: V) -> V {
        let [a, b] = self.across;
        let mut point = mix(mix([0.; 3], x, a), y, b);
        point[2] += z;
        point
    }

    fn local(&self, point: V) -> V {
        let [a, b] = self.across;
        [dot(point, a), dot(point, b), point[2]]
    }

    fn half_chord(&self, y: f64) -> f64 {
        ((self.radius - y) * (self.radius + y)).max(0.).sqrt()
    }

    /// `x'` range of removed material over the tool angle, as `[low, high]`.
    fn reach(&self, angle: f64) -> [f64; 2] {
        let w = self.half_chord(self.offset + self.rho * angle.cos());
        [
            match self.ends[0] {
                End::Through => -w,
                End::Cap => self.stations[0],
            },
            match self.ends[1] {
                End::Through => w,
                End::Cap => self.stations[1],
            },
        ]
    }

    pub fn shaft(&self) -> &Revolution {
        &self.shaft
    }

    /// Meridian edge, radius and axial span of the bored wall.
    pub fn wall(&self) -> (usize, f64, [f64; 2]) {
        (self.wall, self.radius, self.span)
    }

    pub fn ends(&self) -> [End; 2] {
        self.ends
    }

    /// The bore chart: unit travel `x'` and offset direction `y'` in the
    /// shaft chart, the axis offset along `y'`, the axial level and the
    /// tool start and end along `x'`.
    pub fn bore_chart(&self) -> ([V; 2], f64, f64, [f64; 2]) {
        (self.across, self.offset, self.level, self.stations)
    }

    /// Tool base, unit travel, radius and height in the shaft chart.
    pub fn tool(&self) -> (V, V, f64, f64) {
        (
            self.chart([self.stations[0], self.offset, self.level]),
            self.across[0],
            self.rho,
            self.stations[1] - self.stations[0],
        )
    }

    /// Absolute clearance that every admitted classification cleared.
    pub fn resolution_mm(&self) -> f64 {
        self.resolution
    }

    /// Trapezoid nodes over the tool angle used by every integral.
    pub fn nodes(&self) -> usize {
        self.nodes
    }

    /// Exact shaft volume less the certified removed volume.
    pub fn volume_mm3(&self) -> f64 {
        self.shaft.volume_mm3() - self.removed
    }

    /// Certified bound of `|volume_mm3 - exact volume|`, including floating
    /// evaluation of the sum.
    pub fn volume_error_mm3(&self) -> f64 {
        self.removed_error + 8. * f64::EPSILON * self.shaft.volume_mm3()
    }

    /// Volume the tool removes from the shaft.
    pub fn removed_mm3(&self) -> f64 {
        self.removed
    }

    /// The bored wall less its holes, in the shaft chart.
    pub fn bored_wall(&self) -> Moments {
        let height = self.span[1] - self.span[0];
        let mut area = TAU * self.radius * height;
        let mut moment = [0., 0., area * (self.span[0] + self.span[1]) / 2.];
        for hole in self.holes.iter().flatten() {
            area -= hole.area_mm2;
            moment = mix(moment, -hole.area_mm2, hole.centroid_mm);
        }
        Moments {
            area_mm2: area,
            centroid_mm: moment.map(|m| m / area),
        }
    }

    /// The retained part of the tool cylinder, in the shaft chart.
    pub fn tool_wall(&self) -> Moments {
        self.tool_wall
    }

    /// Center, outward normal and radius of a planar blind bottom.
    pub fn cap(&self, side: usize) -> Option<(V, V, f64)> {
        (self.ends[side] == End::Cap).then(|| {
            let center = self.chart([self.stations[side], self.offset, self.level]);
            let sign = if side == 0 { 1. } else { -1. };
            (center, self.across[0].map(|c| sign * c), self.rho)
        })
    }

    /// Point of the curve where the tool wall leaves material on `side`, at
    /// the tool angle measured from `y'` toward `z`.
    pub fn rim_point(&self, side: usize, angle: f64) -> V {
        let x = self.reach(angle)[side];
        let (sin, cos) = angle.sin_cos();
        self.chart([x, self.offset + self.rho * cos, self.level + self.rho * sin])
    }

    /// Angle stations at which every rim of the bore is sampled together.
    pub fn rim_stations(&self, options: Tessellation) -> Result<usize, MeshError> {
        options.check()?;
        let step = options
            .angle
            .min((8. * options.deflection / self.bend).sqrt());
        let count = (TAU / step).ceil().max(8.);
        if count > options.triangle_limit as f64 {
            return Err(MeshError::BudgetExceeded);
        }
        Ok(count as usize)
    }

    /// Whether a shaft-chart point lies in the closed result, within
    /// `tolerance` of its boundary.
    pub fn contains(&self, point: V, tolerance: f64) -> bool {
        if !self.shaft.contains(point, tolerance) {
            return false;
        }
        let [x, y, z] = self.local(point);
        let inner = self.rho - tolerance;
        let deep = inner > 0.
            && (y - self.offset).hypot(z - self.level) < inner
            && x > self.stations[0] + tolerance
            && x < self.stations[1] - tolerance;
        !deep
    }

    /// Shared tessellation: the shaft mesh with its bored wall replaced, the
    /// retained tool wall as face `faces` and blind bottoms after it.
    pub fn mesh(&self, options: Tessellation) -> Result<DisplayMesh, MeshError> {
        let base = revolution_mesh::tessellate(&self.shaft, options)?;
        let bored = self
            .shaft
            .wall_face(0, self.wall)
            .expect("validated bored wall");
        let faces = self.shaft_faces();
        let n = self.rim_stations(options)?;
        let mut mesh = DisplayMesh {
            vertices: base.vertices.clone(),
            indices: Vec::with_capacity(base.indices.len()),
            face_indices: Vec::with_capacity(base.face_indices.len()),
        };
        // Stations of the bored wall's two rims, by angle.
        let segments = {
            let count = base
                .face_indices
                .iter()
                .filter(|face| **face == bored)
                .count();
            count / 2
        };
        let mut rims = [vec![usize::MAX; segments], vec![usize::MAX; segments]];
        let step = TAU / segments as f64;
        for (triangle, &face) in base
            .indices
            .as_chunks::<3>()
            .0
            .iter()
            .zip(&base.face_indices)
        {
            if face != bored {
                mesh.indices.extend_from_slice(triangle);
                mesh.face_indices.push(face);
                continue;
            }
            for &id in triangle {
                let id = id as usize;
                let [x, y, z] = [0, 1, 2].map(|i| base.vertices[3 * id + i]);
                let station = ((y.atan2(x) / step).round() as i64).rem_euclid(segments as i64);
                let rim = usize::from(z == self.span[1]);
                if z != self.span[rim] {
                    return Err(MeshError::UnresolvedBoundary);
                }
                rims[rim][station as usize] = id;
            }
        }
        if segments < 3 || rims.iter().flatten().any(|id| *id == usize::MAX) {
            return Err(MeshError::UnresolvedBoundary);
        }
        let angles: Vec<f64> = (0..n).map(|k| TAU * k as f64 / n as f64).collect();
        let push = |mesh: &mut DisplayMesh, point: V| {
            mesh.vertices.extend(point);
            mesh.vertices.len() / 3 - 1
        };
        let lines: [Vec<usize>; 2] = std::array::from_fn(|side| {
            angles
                .iter()
                .map(|&angle| push(&mut mesh, self.rim_point(side, angle)))
                .collect()
        });
        let point = |mesh: &DisplayMesh, id: usize| -> V {
            std::array::from_fn(|i| mesh.vertices[3 * id + i])
        };
        // Orient each triangle by the outward normal of its support.
        let emit = move |mesh: &mut DisplayMesh, [a, b, c]: [usize; 3], face: u32, normal: V| {
            let [pa, pb, pc] = [a, b, c].map(|id| point(mesh, id));
            let turn = cross(mix(pb, -1., pa), mix(pc, -1., pa));
            if dot(turn, normal) >= 0. {
                mesh.triangle([a, b, c], face);
            } else {
                mesh.triangle([a, c, b], face);
            }
        };
        let tool_face = faces;
        let [y_axis, z_axis] = [self.across[1], [0., 0., 1.]];
        for k in 0..n {
            let next = (k + 1) % n;
            let middle = TAU * (k as f64 + 0.5) / n as f64;
            let (sin, cos) = middle.sin_cos();
            let normal = mix(mix([0.; 3], -cos, y_axis), -sin, z_axis);
            let [l0, l1, h0, h1] = [lines[0][k], lines[0][next], lines[1][k], lines[1][next]];
            emit(&mut mesh, [l0, l1, h1], tool_face, normal);
            emit(&mut mesh, [l0, h1, h0], tool_face, normal);
        }
        let mut next_face = tool_face + 1;
        for (side, line) in lines.iter().enumerate() {
            let Some((_, normal, _)) = self.cap(side) else {
                continue;
            };
            let ring: Vec<[usize; 2]> = (0..n).map(|k| [k, (k + 1) % n]).collect();
            let points = angles
                .iter()
                .map(|angle| Point2::new(angle.cos(), angle.sin()))
                .collect();
            for triangle in triangulate_region(points, &ring, 0)?.triangles {
                emit(&mut mesh, triangle.map(|k| line[k]), next_face, normal);
            }
            next_face += 1;
        }
        self.bored_wall_mesh(&mut mesh, &rims, &lines, bored, &emit)?;
        if mesh.face_indices.len() > options.triangle_limit {
            return Err(MeshError::BudgetExceeded);
        }
        Ok(mesh)
    }

    /// Triangulate the bored wall in its unrolled chart from the shaft's
    /// rim stations and the sampled hole rims.
    fn bored_wall_mesh(
        &self,
        mesh: &mut DisplayMesh,
        rims: &[Vec<usize>; 2],
        lines: &[Vec<usize>; 2],
        bored: u32,
        emit: &Emit,
    ) -> Result<(), MeshError> {
        let segments = rims[0].len();
        let step = TAU / segments as f64;
        let angle_of = |id: usize| {
            let [x, y] = [mesh.vertices[3 * id], mesh.vertices[3 * id + 1]];
            y.atan2(x).rem_euclid(TAU)
        };
        // Unroll from the rim station farthest from every hole sample.
        let samples: Vec<f64> = (0..2)
            .filter(|&side| self.ends[side] == End::Through)
            .flat_map(|side| lines[side].iter().map(|&id| angle_of(id)))
            .collect();
        let gap = |angle: f64| {
            samples
                .iter()
                .map(|sample| {
                    let d = (sample - angle).rem_euclid(TAU);
                    d.min(TAU - d)
                })
                .fold(f64::INFINITY, f64::min)
        };
        let seam = (0..segments)
            .max_by(|a, b| gap(*a as f64 * step).total_cmp(&gap(*b as f64 * step)))
            .expect("rim stations");
        let start = seam as f64 * step;
        if gap(start) <= step / 4. {
            return Err(MeshError::UnresolvedBoundary);
        }
        let unrolled = |angle: f64| start + (angle - start).rem_euclid(TAU);
        let mut points: Vec<Point2<f64>> = Vec::new();
        let mut ids: Vec<usize> = Vec::new();
        let mut edges: Vec<[usize; 2]> = Vec::new();
        let r = self.radius;
        // Outer boundary, counterclockwise in (r * angle, z).
        let mut outer: Vec<(f64, f64, usize)> = Vec::new();
        for k in 0..=segments {
            let station = (seam + k) % segments;
            outer.push((
                r * (start + k as f64 * step),
                self.span[0],
                rims[0][station],
            ));
        }
        for k in (0..=segments).rev() {
            let station = (seam + k) % segments;
            outer.push((
                r * (start + k as f64 * step),
                self.span[1],
                rims[1][station],
            ));
        }
        let first = points.len();
        for (u, v, id) in &outer {
            points.push(Point2::new(*u, *v));
            ids.push(*id);
        }
        let count = outer.len();
        edges.extend((0..count).map(|i| [first + i, first + (i + 1) % count]));
        let mut holes = 0;
        for (line, end) in lines.iter().zip(self.ends) {
            if end != End::Through {
                continue;
            }
            let ring: Vec<(Point2<f64>, usize)> = line
                .iter()
                .map(|&id| {
                    let z = mesh.vertices[3 * id + 2];
                    (Point2::new(r * unrolled(angle_of(id)), z), id)
                })
                .collect();
            let area: f64 = (0..ring.len())
                .map(|i| {
                    let [a, b] = [ring[i].0, ring[(i + 1) % ring.len()].0];
                    a.x * b.y - a.y * b.x
                })
                .sum();
            let first = points.len();
            let ordered: Vec<_> = if area > 0. {
                ring.into_iter().rev().collect()
            } else {
                ring
            };
            let length = ordered.len();
            for (point, id) in ordered {
                points.push(point);
                ids.push(id);
            }
            edges.extend((0..length).map(|i| [first + i, first + (i + 1) % length]));
            holes += 1;
        }
        // Steiner columns at every rim station across the hole band keep
        // each triangle within about one station of angle. Points within a
        // station of a hole stay out so no sliver meets its rim.
        let band = [
            (self.level - self.rho - 2. * r * step).max(self.span[0]),
            (self.level + self.rho + 2. * r * step).min(self.span[1]),
        ];
        let rows = ((band[1] - band[0]) / (r * step)).ceil().max(1.) as usize;
        let clearance = self.rho + 1.5 * r * step;
        let mut interior = 0;
        for k in 0..segments {
            let angle = start + k as f64 * step;
            for row in 0..=rows {
                let z = band[0] + (band[1] - band[0]) * row as f64 / rows as f64;
                if z - self.span[0] < 0.5 * r * step || self.span[1] - z < 0.5 * r * step {
                    continue;
                }
                let (sin, cos) = angle.sin_cos();
                let [x, y, _] = self.local([r * cos, r * sin, z]);
                let side = usize::from(x > 0.);
                if self.ends[side] == End::Through
                    && (y - self.offset).hypot(z - self.level) < clearance
                {
                    continue;
                }
                if k == 0 {
                    // The seam column is the outer boundary itself.
                    continue;
                }
                points.push(Point2::new(r * angle, z));
                ids.push(mesh.vertices.len() / 3);
                mesh.vertices.extend([r * cos, r * sin, z]);
                interior += 1;
            }
        }
        if mesh.vertices.len() / 3 > u32::MAX as usize {
            return Err(MeshError::BudgetExceeded);
        }
        let region = triangulate_region_with_interior(points, &edges, holes, interior)?;
        for triangle in region.triangles {
            let [a, b, c] = triangle.map(|i| ids[i]);
            let centroid: V = std::array::from_fn(|i| {
                (mesh.vertices[3 * a + i] + mesh.vertices[3 * b + i] + mesh.vertices[3 * c + i])
                    / 3.
            });
            let normal = [centroid[0], centroid[1], 0.];
            emit(mesh, [a, b, c], bored, normal);
        }
        Ok(())
    }

    /// Number of faces of the unbored shaft; the tool wall follows them.
    pub fn shaft_faces(&self) -> u32 {
        (0..self.shaft.rings()[0].len())
            .filter_map(|edge| self.shaft.wall_face(0, edge))
            .max()
            .map_or(0, |last| last + 1)
    }

    /// Samples of the bored-wall rim on a through `side`, closed at both ends.
    pub fn sample_rim(&self, side: usize, options: Tessellation) -> Result<Vec<V>, MeshError> {
        let n = self.rim_stations(options)?;
        Ok((0..=n)
            .map(|k| self.rim_point(side, TAU * k as f64 / n as f64))
            .collect())
    }

    /// Range of the solid along a direction in the shaft chart. A cross-bore
    /// inside one wall's span never removes an extreme point of the shaft.
    pub fn extent(&self, direction: V) -> [f64; 2] {
        [
            -self.shaft.support(direction.map(|c| -c)),
            self.shaft.support(direction),
        ]
    }
}

/// Classify the tool against the shaft and build the exact cut.
pub fn subtract(
    shaft: &Revolution,
    tool: Tool,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Cut, Error> {
    let check = || {
        if interrupt.is_some_and(|stop| stop()) {
            Err(Error::Cancelled)
        } else {
            Ok(())
        }
    };
    check()?;
    let walls = stepped_walls(shaft).ok_or(Error::Unsupported)?;
    let finite = |v: f64| v.is_finite() && v.abs() <= MAXIMUM_MM;
    let [ax, ay, az] = tool.axis;
    let norm = ax.hypot(ay);
    if !tool.base_mm.into_iter().all(finite)
        || az != 0.
        || !(1. - 1e-12..=1. + 1e-12).contains(&norm)
        || !(finite(tool.radius_mm) && tool.radius_mm > 0.)
        || !(finite(tool.height_mm) && tool.height_mm > 0.)
    {
        return Err(Error::Invalid);
    }
    let travel = [ax / norm, ay / norm, 0.];
    let side = cross([0., 0., 1.], travel);
    let [s0, e, c] = [
        dot(tool.base_mm, travel),
        dot(tool.base_mm, side),
        tool.base_mm[2],
    ];
    let rho = tool.radius_mm;
    let s1 = s0 + tool.height_mm;
    let low = walls[0].2[0];
    let high = walls[walls.len() - 1].2[1];
    let scale = walls
        .iter()
        .flat_map(|wall| [wall.1, wall.2[0].abs(), wall.2[1].abs()])
        .chain([s0.abs(), s1.abs(), e.abs(), c.abs(), rho])
        .fold(1., f64::max);
    let margin = RESOLUTION * scale;
    let area_margin = 4. * margin * scale;
    let literal = tool.literal && norm == 1. && (travel[0] == 0. || travel[1] == 0.);
    // Exact operands in the bore chart.
    let base = |i: usize| exact(tool.base_mm[i]);
    let signed = |v: f64| exact(v);
    let xe = || {
        if travel[0] != 0. {
            base(0) * signed(travel[0])
        } else {
            base(1) * signed(travel[1])
        }
    };
    let ye = || {
        if travel[0] != 0. {
            base(1) * signed(side[1])
        } else {
            base(0) * signed(side[0])
        }
    };
    let [re, he, ce] = [exact(rho), exact(tool.height_mm), base(2)];

    // Along the axis: miss, touch, or one wall that holds the whole band.
    let below = decide(c - rho - high, margin, literal, || {
        ce.clone() - re.clone() - exact(high)
    })?;
    let above = decide(low - c - rho, margin, literal, || {
        exact(low) - ce.clone() - re.clone()
    })?;
    if below != Ordering::Less || above != Ordering::Less {
        return Ok(Cut::Unchanged);
    }
    let abs_y = || {
        let y = ye();
        if y < BigRational::from_integer(0.into()) {
            -y
        } else {
            y
        }
    };
    let reach = e.abs();
    // The widest wall the tool section can meet bounds the shaft in that
    // band: beyond it, radially or along the tool, nothing is touched.
    let widest_wall = walls
        .iter()
        .filter(|wall| wall.2[0] < c + rho && wall.2[1] > c - rho)
        .map(|wall| wall.1)
        .fold(0., f64::max);
    let apart = decide(reach - rho - widest_wall, margin, literal, || {
        abs_y() - re.clone() - exact(widest_wall)
    })?;
    let before = decide(s0 + tool.height_mm + widest_wall, margin, literal, || {
        xe() + he.clone() + exact(widest_wall)
    })?;
    let after = decide(widest_wall - s0, margin, literal, || {
        exact(widest_wall) - xe()
    })?;
    if apart != Ordering::Less || before != Ordering::Greater || after != Ordering::Greater {
        return Ok(Cut::Unchanged);
    }
    let (wall, radius, span) = walls
        .iter()
        .copied()
        .find(|wall| wall.2[0] <= c && c <= wall.2[1])
        .ok_or(Error::Unsupported)?;
    for gap in [c - rho - span[0], span[1] - c - rho] {
        if gap <= margin {
            return Err(if gap < -margin {
                Error::Unsupported
            } else {
                Error::Unresolved
            });
        }
    }
    check()?;
    // Across the axis: miss or touch, tangency, or a section inside the wall.
    let apart = decide(reach - rho - radius, margin, literal, || {
        abs_y() - re.clone() - exact(radius)
    })?;
    if apart != Ordering::Less {
        return Ok(Cut::Unchanged);
    }
    let inside = decide(radius - reach - rho, margin, literal, || {
        exact(radius) - abs_y() - re.clone()
    })?;
    if inside != Ordering::Greater {
        return Err(Error::Unsupported);
    }
    // Along the tool: compare each station with the narrowest and widest
    // half chord of the shaft over the tool section, by their squares.
    let near = (reach - rho).max(0.);
    let far = reach + rho;
    let station_exact = |end: usize| {
        if end == 0 { xe() } else { xe() + he.clone() }
    };
    let widest = |s: f64, end: usize| {
        decide(
            s * s + near * near - radius * radius,
            area_margin,
            literal,
            || {
                let s = station_exact(end);
                let m = abs_y() - re.clone();
                let m = if m < BigRational::from_integer(0.into()) {
                    BigRational::from_integer(0.into())
                } else {
                    m
                };
                &s * &s + &m * &m - exact(radius) * exact(radius)
            },
        )
    };
    let narrowest = |s: f64, end: usize| {
        decide(
            s * s + far * far - radius * radius,
            area_margin,
            literal,
            || {
                let s = station_exact(end);
                let m = abs_y() + re.clone();
                &s * &s + &m * &m - exact(radius) * exact(radius)
            },
        )
    };
    let mut ends = [End::Through; 2];
    for (end, s) in [s0, s1].into_iter().enumerate() {
        // Outward is negative travel at the start and positive at the end.
        let outward = if end == 0 { s < 0. } else { s > 0. };
        let wide = widest(s, end)?;
        if wide != Ordering::Less {
            if outward && wide == Ordering::Greater {
                continue;
            }
            if !outward {
                // The tool lies entirely beyond the shaft or touches it.
                return Ok(Cut::Unchanged);
            }
            return Err(Error::Unsupported);
        }
        if narrowest(s, end)? == Ordering::Less {
            ends[end] = End::Cap;
            continue;
        }
        return Err(Error::Unsupported);
    }
    if ends == [End::Cap; 2] {
        // An enclosed void needs a second shell.
        return Err(Error::Unsupported);
    }
    check()?;
    let mut bore = CrossBore {
        shaft: shaft.clone(),
        wall,
        radius,
        span,
        across: [travel, side],
        offset: e,
        level: c,
        rho,
        stations: [s0, s1],
        ends,
        resolution: margin,
        nodes: 0,
        removed: 0.,
        removed_error: 0.,
        holes: [None; 2],
        tool_wall: Moments {
            area_mm2: 0.,
            centroid_mm: [0.; 3],
        },
        bend: 0.,
    };
    integrate(&mut bore, check)?;
    Ok(Cut::Body(Box::new(bore)))
}

/// Volume, face areas and centroids by certified periodic trapezoid sums.
fn integrate(bore: &mut CrossBore, check: impl Fn() -> Result<(), Error>) -> Result<(), Error> {
    let (r, e, rho, c) = (bore.radius, bore.offset, bore.rho, bore.level);
    let q = (r - e.abs()) / rho;
    if q.is_nan() || q <= 1. {
        return Err(Error::Unresolved);
    }
    let width = q.acosh() / 2.;
    let stretch = width.cosh();
    let reach = e.abs() + rho * stretch;
    let strip = Strip {
        width,
        upper: r * r + reach * reach,
        lower: (r - reach) * (r + reach),
        reach,
    };
    if strip.lower.is_nan() || strip.lower <= 0. {
        return Err(Error::Unresolved);
    }
    // Strip bounds of each integrand: sin² by cosh², |w| and 1/|w| by the
    // bounds of |r² - y²|, |y| by its reach.
    let trig = stretch * stretch;
    let root_upper = strip.upper.sqrt();
    let root_lower = strip.lower.sqrt();
    let bounds = [
        trig * root_upper,                   // sin² w
        root_upper,                          // w
        trig * r / root_lower,               // r sin² / w
        trig * r * strip.reach / root_lower, // r sin² y / w
        strip.reach * root_upper,            // y w
        root_upper * root_upper,             // w²
    ];
    let target = 4. * f64::EPSILON * TAU * r.max(rho);
    let mut n = 16;
    for bound in bounds {
        n = n.max(strip.nodes(bound, target).ok_or(Error::Unresolved)?);
    }
    check()?;
    let h = TAU / n as f64;
    let mut sums = [0.; 6];
    let mut magnitude = [0.; 6];
    for k in 0..n {
        if k % 4096 == 0 {
            check()?;
        }
        let (sin, cos) = (h * k as f64).sin_cos();
        let y = e + rho * cos;
        let w = ((r - y) * (r + y)).sqrt();
        let terms = [
            sin * sin * w,
            w,
            r * sin * sin / w,
            r * sin * sin * y / w,
            y * w,
            w * w,
        ];
        for i in 0..6 {
            sums[i] += terms[i];
            magnitude[i] += terms[i].abs();
        }
    }
    // Truncation by the strip bound, rounding by the summed magnitudes.
    let value = |i: usize| h * sums[i];
    let error =
        |i: usize| strip.error(n, bounds[i]) + (n as f64 + 16.) * f64::EPSILON * h * magnitude[i];
    let through = bore.ends.iter().filter(|end| **end == End::Through).count() as f64;
    let caps = [
        if bore.ends[0] == End::Cap {
            -bore.stations[0]
        } else {
            0.
        },
        if bore.ends[1] == End::Cap {
            bore.stations[1]
        } else {
            0.
        },
    ];
    let disc = PI * rho * rho;
    let cap_sum = caps[0] + caps[1];
    bore.removed = through * rho * rho * value(0) + disc * cap_sum;
    bore.removed_error = through * rho * rho * error(0)
        + 8. * f64::EPSILON * (bore.removed.abs() + disc * cap_sum.abs());
    // Holes in the bored wall: area element r / w over the tool section.
    for side in 0..2 {
        if bore.ends[side] != End::Through {
            continue;
        }
        let sign = if side == 0 { -1. } else { 1. };
        let area = rho * rho * value(2);
        let first = [sign * r * PI * rho * rho, rho * rho * value(3), c * area];
        bore.holes[side] = Some(Moments {
            area_mm2: area,
            centroid_mm: bore.chart(first.map(|m| m / area)),
        });
    }
    // The retained tool wall: x' from the low to the high reach.
    let width_sum = through * value(1) + TAU * cap_sum;
    let area = rho * width_sum;
    let y_moment = rho * (through * value(4) + cap_sum * TAU * e);
    let x_moment = rho
        * [0, 1]
            .into_iter()
            .map(|side| {
                let sign = if side == 0 { -1. } else { 1. };
                match bore.ends[side] {
                    End::Through => sign * value(5) / 2.,
                    End::Cap => sign * TAU * bore.stations[side].powi(2) / 2.,
                }
            })
            .sum::<f64>();
    bore.tool_wall = Moments {
        area_mm2: area,
        centroid_mm: bore.chart([x_moment / area, y_moment / area, c]),
    };
    let w_min = root(r * r - (e.abs() + rho).powi(2));
    bore.bend = rho + (r * rho).powi(2) / w_min.powi(3) + (e.abs() + rho) * rho / w_min;
    bore.nodes = n;
    Ok(())
}

fn root(value: f64) -> f64 {
    value.max(0.).sqrt()
}

#[cfg(test)]
mod tests;
