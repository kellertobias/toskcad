//! Watertight tessellation of the retained boundary.
//!
//! Every edge is sampled once and each face triangulates its parameter
//! region from exactly those samples, so neighbouring faces share their
//! boundary vertices. Interior edges are then split until every interior
//! chord and every triangle centre lies within the requested deflection of
//! the exact surface.
use super::domain::{self, Domain, Key, Segment};
use super::polygon::{self, Triangle};
use super::{
    Error, Interrupt, Point, Result, Solid, Surface, add, checkpoint, cross, distance, ensure,
    norm, scale, turn,
};
use crate::display::DisplayMesh;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Options {
    /// Largest admitted distance between a chord or triangle centre and
    /// the exact boundary.
    pub deflection_mm: f64,
    /// Largest admitted turn of the surface along one chord.
    pub angle_rad: f64,
    pub triangle_limit: usize,
}

/// One triangle of a face, with the surface parameters of its corners.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Cell {
    pub vertices: [u32; 3],
    pub uv: [[f64; 2]; 3],
}

impl Cell {
    /// Parameters a fraction `s` of the way along side `i`.
    pub(crate) fn along(&self, i: usize, s: f64) -> [f64; 2] {
        let (a, b) = (self.uv[i], self.uv[(i + 1) % 3]);
        [a[0] + (b[0] - a[0]) * s, a[1] + (b[1] - a[1]) * s]
    }

    fn centre(&self) -> [f64; 2] {
        let c = self.uv;
        [
            (c[0][0] + c[1][0] + c[2][0]) / 3.,
            (c[0][1] + c[1][1] + c[2][1]) / 3.,
        ]
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Mesh {
    pub positions: Vec<Point>,
    /// Outward counter-clockwise triangles.
    pub triangles: Vec<[u32; 3]>,
    /// The boundary face of each triangle.
    pub faces: Vec<u32>,
    /// Largest distance measured between an interior chord midpoint or
    /// triangle centre and the exact surface.
    pub maximum_deviation_mm: f64,
    pub(crate) cells: Vec<Cell>,
    /// The edge piece under each triangle side that lies on the boundary.
    pub(crate) sides: Vec<[Option<Segment>; 3]>,
}

impl Mesh {
    /// The mesh in the kernel's flat display layout.
    pub fn display(&self) -> DisplayMesh {
        DisplayMesh {
            vertices: self.positions.iter().flatten().copied().collect(),
            indices: self.triangles.iter().flatten().copied().collect(),
            face_indices: self.faces.clone(),
        }
    }
}

fn pair(a: u32, b: u32) -> (u32, u32) {
    (a.min(b), a.max(b))
}

/// A vertex of one face while it is refined. Copies of a pole along its
/// row are separate here and share one mesh vertex.
#[derive(Clone, Copy)]
struct Local {
    position: Point,
    /// Unit surface normal, or zero where the surface collapses.
    normal: Point,
    vertex: u32,
    pole: bool,
}

fn normal(surface: &Surface, uv: [f64; 2]) -> Point {
    let (_, [su, sv]) = surface.evaluate(uv);
    let n = cross(su, sv);
    let (short, long) = (norm(su).min(norm(sv)), norm(su).max(norm(sv)));
    // Next to a pole the direction of the normal is rounding noise.
    if short > 1e-6 * long && norm(n) > 0. {
        scale(n, 1. / norm(n))
    } else {
        [0.; 3]
    }
}

struct Refiner<'a> {
    surface: &'a Surface,
    positions: &'a mut Vec<Point>,
    locals: Vec<Local>,
    cells: Vec<Cell>,
    fixed: BTreeSet<(u32, u32)>,
    deflection: f64,
    angle: f64,
    deviation: f64,
}

impl Refiner<'_> {
    fn insert(&mut self, uv: [f64; 2]) -> u32 {
        let position = self.surface.evaluate(uv).0;
        let vertex = self.positions.len() as u32;
        self.positions.push(position);
        self.locals.push(Local {
            position,
            normal: normal(self.surface, uv),
            vertex,
            pole: false,
        });
        self.locals.len() as u32 - 1
    }

    /// One pass: split the interior sides that stray from the surface or
    /// turn too far, at most one per triangle. Sides on boundary edges and
    /// along a pole row are never split.
    fn pass(&mut self) -> Result<bool> {
        let mut sides = BTreeMap::<(u32, u32), Vec<(usize, usize)>>::new();
        for (index, cell) in self.cells.iter().enumerate() {
            for i in 0..3 {
                sides
                    .entry(pair(cell.vertices[i], cell.vertices[(i + 1) % 3]))
                    .or_default()
                    .push((index, i));
            }
        }
        // Chord length of every free side, and whether it must go.
        let mut lengths = BTreeMap::new();
        let mut wanted = BTreeSet::new();
        self.deviation = 0.;
        for (&(a, b), uses) in &sides {
            let [from, to] = [a, b].map(|v| self.locals[v as usize]);
            let row = from.pole && to.pole;
            if self.fixed.contains(&(a, b)) || row {
                continue;
            }
            ensure(uses.len() == 2, Error::Unsupported)?;
            lengths.insert((a, b), distance(from.position, to.position));
            let (index, i) = uses[0];
            let cell = &self.cells[index];
            let [from, to] = [i, (i + 1) % 3].map(|k| self.locals[cell.vertices[k] as usize]);
            // The chord against the surface along it; three stations see
            // a side that crosses the surface at its middle.
            let mut deviation: f64 = 0.;
            for s in [0.25, 0.5, 0.75] {
                let chord = add(scale(from.position, 1. - s), scale(to.position, s));
                deviation =
                    deviation.max(distance(self.surface.evaluate(cell.along(i, s)).0, chord));
            }
            // A chord far below the deflection cannot stray by it.
            let resolved = lengths[&(a, b)] <= 1e-3 * self.deflection;
            if !resolved
                && (deviation > self.deflection || turn(from.normal, to.normal) > self.angle)
            {
                wanted.insert((a, b));
            } else {
                self.deviation = self.deviation.max(deviation);
            }
        }
        // A triangle is only ever cut across its longest free side, so no
        // triangle thins out while a longer side waits.
        let mut chosen = BTreeSet::new();
        for cell in &self.cells {
            let keys = [0, 1, 2].map(|i| pair(cell.vertices[i], cell.vertices[(i + 1) % 3]));
            if keys.iter().any(|key| wanted.contains(key)) {
                let longest = keys
                    .iter()
                    .filter_map(|key| lengths.get(key).map(|&length| (length, *key)))
                    .fold(None, |best: Option<(f64, (u32, u32))>, side| match best {
                        Some(best) if best.0 >= side.0 => Some(best),
                        _ => Some(side),
                    });
                chosen.extend(longest.map(|(_, key)| key));
            }
        }
        let mut wanted: Vec<_> = chosen
            .into_iter()
            .map(|(a, b)| (lengths[&(a, b)], a, b))
            .collect();
        wanted.sort_by(|x, y| y.0.total_cmp(&x.0).then(x.1.cmp(&y.1)).then(x.2.cmp(&y.2)));
        let mut taken = vec![false; self.cells.len()];
        let mut changed = !wanted.is_empty();
        for (_, a, b) in wanted {
            let uses = &sides[&(a, b)];
            if uses.iter().any(|&(index, _)| taken[index]) {
                continue;
            }
            let (index, i) = uses[0];
            let vertex = self.insert(self.cells[index].along(i, 0.5));
            for &(index, i) in uses {
                taken[index] = true;
                let cell = self.cells[index];
                let uv = cell.along(i, 0.5);
                let j = (i + 1) % 3;
                let mut first = cell;
                (first.vertices[j], first.uv[j]) = (vertex, uv);
                let mut second = cell;
                (second.vertices[i], second.uv[i]) = (vertex, uv);
                self.cells[index] = first;
                self.cells.push(second);
            }
        }
        // A triangle whose sides all conform can still bulge at its centre.
        for index in 0..taken.len() {
            let cell = self.cells[index];
            let centre = cell.centre();
            let exact = self.surface.evaluate(centre).0;
            let [a, b, c] = cell.vertices.map(|v| self.locals[v as usize].position);
            let deviation = distance(exact, scale(add(add(a, b), c), 1. / 3.));
            // A triangle on a pole row has no extent on the surface.
            let on_row = cell
                .vertices
                .iter()
                .filter(|&&v| self.locals[v as usize].pole)
                .count()
                > 1;
            if taken[index] || on_row || deviation <= self.deflection {
                if !taken[index] && !on_row {
                    self.deviation = self.deviation.max(deviation);
                }
                continue;
            }
            changed = true;
            let vertex = self.insert(centre);
            for i in 0..3 {
                let mut part = cell;
                (part.vertices[i], part.uv[i]) = (vertex, centre);
                if i == 0 {
                    self.cells[index] = part;
                } else {
                    self.cells.push(part);
                }
            }
        }
        Ok(changed)
    }
}

struct Builder<'a> {
    solid: &'a Solid,
    samples: Vec<Vec<f64>>,
    vertices: BTreeMap<Key, u32>,
    mesh: Mesh,
}

impl Builder<'_> {
    fn vertex(&mut self, corner: &polygon::Corner, surface: &Surface) -> u32 {
        let boundary = self.solid.boundary();
        let next = self.mesh.positions.len() as u32;
        let id = *self.vertices.entry(corner.key).or_insert(next);
        if id == next {
            self.mesh.positions.push(match corner.key {
                Key::Vertex(vertex) => boundary.vertices[vertex],
                Key::Sample(edge, k) => {
                    boundary.edges[edge].curve.evaluate(self.samples[edge][k]).0
                }
                Key::Extra(..) => surface.evaluate(corner.uv).0,
            });
        }
        id
    }

    fn face(&mut self, index: usize, options: Options, outward: bool) -> Result<()> {
        let solid = self.solid;
        let surface = &solid.boundary().faces[index].surface;
        let domain = domain::build(solid.boundary(), index, &self.samples, solid.tolerance())?;
        let triangles: Vec<Triangle> = polygon::triangulate(&domain, surface, index)?;
        // Local vertices: one per key, and one per parameter along a pole row.
        let mut locals = Vec::new();
        let mut local_of = BTreeMap::new();
        let mut cells = Vec::with_capacity(triangles.len());
        for triangle in &triangles {
            let vertices = [0, 1, 2].map(|i| {
                let corner = &triangle[i];
                let vertex = self.vertex(corner, surface);
                *local_of.entry((vertex, corner.copy)).or_insert_with(|| {
                    locals.push(Local {
                        position: self.mesh.positions[vertex as usize],
                        normal: normal(surface, corner.uv),
                        vertex,
                        pole: corner.pole,
                    });
                    locals.len() as u32 - 1
                })
            });
            cells.push(Cell {
                vertices,
                uv: [0, 1, 2].map(|i| triangle[i].uv),
            });
        }
        // Boundary edge pieces by mesh vertex pair, with their parameters:
        // a pole or a seam point meets the same vertices elsewhere too.
        let periods = surface.periods();
        let mut pieces = BTreeMap::<(u32, u32), Vec<(Segment, [[f64; 2]; 2])>>::new();
        if let Domain::Trimmed { loops, .. } = &domain {
            for trim in loops {
                for (i, segment) in trim.segments.iter().enumerate() {
                    let Some(segment) = segment else { continue };
                    let next = (i + 1) % trim.points.len();
                    let ends = [i, next].map(|k| self.vertices[&trim.points[k].key]);
                    let uv = [trim.points[i].uv, trim.at(i + 1, periods)];
                    let back = Segment {
                        edge: segment.edge,
                        t: [segment.t[1], segment.t[0]],
                    };
                    pieces
                        .entry((ends[0], ends[1]))
                        .or_default()
                        .push((*segment, uv));
                    pieces
                        .entry((ends[1], ends[0]))
                        .or_default()
                        .push((back, [uv[1], uv[0]]));
                }
            }
        }
        let same = |a: [f64; 2], b: [f64; 2]| {
            (0..2).all(|axis| {
                let turns =
                    periods[axis].map_or(0., |period| ((a[axis] - b[axis]) / period).round());
                (a[axis] - b[axis] - turns * periods[axis].unwrap_or(0.)).abs() <= 1e-9
            })
        };
        let piece = |vertices: (u32, u32), from: [f64; 2], to: [f64; 2]| {
            let found = pieces
                .get(&vertices)?
                .iter()
                .find(|(_, uv)| same(uv[0], from) && same(uv[1], to));
            found.map(|&(segment, _)| segment)
        };
        let mut fixed = BTreeSet::new();
        for cell in &cells {
            for i in 0..3 {
                let (a, b) = (cell.vertices[i], cell.vertices[(i + 1) % 3]);
                let vertices = (locals[a as usize].vertex, locals[b as usize].vertex);
                if piece(vertices, cell.uv[i], cell.uv[(i + 1) % 3]).is_some() {
                    fixed.insert(pair(a, b));
                }
            }
        }
        let mut refiner = Refiner {
            surface,
            positions: &mut self.mesh.positions,
            locals,
            cells,
            fixed,
            deflection: options.deflection_mm,
            angle: options.angle_rad,
            deviation: 0.,
        };
        let mut passes = 0;
        while !surface.is_plane() && refiner.pass()? {
            passes += 1;
            ensure(
                passes < 4096 && refiner.cells.len() <= options.triangle_limit,
                Error::Limit,
            )?;
        }
        let Refiner {
            locals,
            cells,
            deviation,
            ..
        } = refiner;
        self.mesh.maximum_deviation_mm = self.mesh.maximum_deviation_mm.max(deviation);
        for mut cell in cells {
            cell.vertices = cell.vertices.map(|v| locals[v as usize].vertex);
            // A triangle left on a pole row has no extent on the surface.
            let [a, b, c] = cell.vertices;
            if a == b || b == c || c == a {
                continue;
            }
            if !outward {
                cell.vertices.swap(1, 2);
                cell.uv.swap(1, 2);
            }
            let side = |i: usize| {
                let j = (i + 1) % 3;
                piece((cell.vertices[i], cell.vertices[j]), cell.uv[i], cell.uv[j])
            };
            self.mesh.sides.push([side(0), side(1), side(2)]);
            self.mesh.triangles.push(cell.vertices);
            self.mesh.faces.push(index as u32);
            self.mesh.cells.push(cell);
        }
        ensure(
            self.mesh.triangles.len() <= options.triangle_limit,
            Error::Limit,
        )
    }
}

/// Tessellate the solid within `options.deflection_mm`.
pub fn tessellate(solid: &Solid, options: Options, interrupt: Interrupt<'_>) -> Result<Mesh> {
    ensure(
        options.deflection_mm.is_finite()
            && options.deflection_mm > 0.
            && options.angle_rad > 0.
            && options.triangle_limit > 0,
        Error::Invalid,
    )?;
    let boundary = solid.boundary();
    let mut builder = Builder {
        solid,
        samples: boundary.sample(options.deflection_mm, options.angle_rad)?,
        vertices: BTreeMap::new(),
        mesh: Mesh::default(),
    };
    let signs = boundary.signs();
    for index in 0..boundary.faces.len() {
        checkpoint(interrupt)?;
        builder.face(index, options, signs[index] > 0.)?;
    }
    Ok(builder.mesh)
}
