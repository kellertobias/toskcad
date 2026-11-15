//! The retained boundary of one solid and its integrated properties.
use super::domain::{self, Domain};
use super::{
    Curve, Error, Interrupt, Point, Result, Surface, checkpoint, distance, ensure, extent, mass,
};
use std::collections::BTreeMap;

/// The largest admitted number of faces, edges or vertices in one solid.
pub const MAX_ELEMENTS: usize = 100_000;

/// An edge runs along its curve from `range[0]` at `vertices[0]` to
/// `range[1]` at `vertices[1]`.
#[derive(Clone, Debug, PartialEq)]
pub struct Edge {
    pub curve: Curve,
    pub vertices: [usize; 2],
    pub range: [f64; 2],
}

impl Edge {
    /// Bound `curve` between two vertices, in the curve's own direction.
    pub fn new(
        curve: Curve,
        vertices: [usize; 2],
        points: [Point; 2],
        tolerance: f64,
    ) -> Result<Self> {
        let start = curve.invert(points[0], None);
        let mut end = curve.invert(points[1], None);
        if let Some(period) = curve.period() {
            let closed = vertices[0] == vertices[1];
            end = start + (end - start).rem_euclid(period);
            if closed || end - start <= 1e-12 * period {
                ensure(closed, Error::Invalid)?;
                end = start + period;
            }
        }
        ensure(end > start, Error::Invalid)?;
        let edge = Self {
            curve,
            vertices,
            range: [start, end],
        };
        for side in 0..2 {
            let on_curve = edge.curve.evaluate(edge.range[side]).0;
            ensure(
                distance(on_curve, points[side]) <= tolerance,
                Error::Invalid,
            )?;
        }
        Ok(edge)
    }
}

/// One trim loop: oriented edge uses, or a lone vertex.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct Loop {
    /// Edge index and whether the loop follows the edge direction.
    pub uses: Vec<(usize, bool)>,
    pub vertex: Option<usize>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Face {
    pub surface: Surface,
    /// Whether the face normal agrees with the surface normal.
    pub same_sense: bool,
    pub loops: Vec<Loop>,
}

/// A closed shell. The first shell of a solid is its exterior and the
/// others bound voids. `reversed` flips every face normal of the shell.
#[derive(Clone, Debug, PartialEq)]
pub struct Shell {
    pub faces: Vec<usize>,
    pub reversed: bool,
}

#[derive(Clone, Debug, PartialEq, Default)]
pub struct Boundary {
    pub vertices: Vec<Point>,
    pub edges: Vec<Edge>,
    pub faces: Vec<Face>,
    pub shells: Vec<Shell>,
}

impl Boundary {
    /// Every edge of a closed oriented shell is used once in each direction
    /// by that shell alone, loops are vertex-continuous, and each shell is
    /// one connected set of faces.
    pub(crate) fn validate(&self) -> Result<()> {
        let counts = [self.vertices.len(), self.edges.len(), self.faces.len()];
        ensure(counts.iter().all(|&n| n <= MAX_ELEMENTS), Error::Limit)?;
        ensure(!self.shells.is_empty(), Error::Invalid)?;
        let mut owner = BTreeMap::new();
        let mut claimed = vec![false; self.faces.len()];
        for (shell_index, shell) in self.shells.iter().enumerate() {
            ensure(!shell.faces.is_empty(), Error::Invalid)?;
            let mut uses = BTreeMap::<usize, [u32; 2]>::new();
            let mut group: Vec<usize> = (0..shell.faces.len()).collect();
            let mut first_face = BTreeMap::new();
            for (local, &face_index) in shell.faces.iter().enumerate() {
                let face = self.faces.get(face_index).ok_or(Error::Invalid)?;
                ensure(
                    !std::mem::replace(&mut claimed[face_index], true),
                    Error::Invalid,
                )?;
                for trim in &face.loops {
                    ensure(
                        trim.vertex.is_none_or(|v| v < self.vertices.len()),
                        Error::Invalid,
                    )?;
                    for (k, &(edge_index, forward)) in trim.uses.iter().enumerate() {
                        let edge = self.edges.get(edge_index).ok_or(Error::Invalid)?;
                        let (next_edge, next_forward) = trim.uses[(k + 1) % trim.uses.len()];
                        let next = self.edges.get(next_edge).ok_or(Error::Invalid)?;
                        ensure(
                            edge.vertices[usize::from(forward)]
                                == next.vertices[usize::from(!next_forward)],
                            Error::Invalid,
                        )?;
                        uses.entry(edge_index).or_default()[usize::from(forward)] += 1;
                        let other = *first_face.entry(edge_index).or_insert(local);
                        let (a, b) = (find(&mut group, local), find(&mut group, other));
                        group[a] = b;
                    }
                }
            }
            for (&edge, &count) in &uses {
                ensure(count == [1, 1], Error::Invalid)?;
                ensure(owner.insert(edge, shell_index).is_none(), Error::Invalid)?;
            }
            let root = find(&mut group, 0);
            ensure(
                (0..group.len()).all(|i| find(&mut group, i) == root),
                Error::Invalid,
            )?;
        }
        ensure(
            claimed.iter().all(|&c| c) && owner.len() == self.edges.len(),
            Error::Invalid,
        )
    }

    pub(crate) fn mapped(
        &self,
        map: &dyn Fn(Point) -> Point,
        linear: &dyn Fn(Point) -> Point,
    ) -> Self {
        Self {
            vertices: self.vertices.iter().map(|&p| map(p)).collect(),
            edges: self
                .edges
                .iter()
                .map(|edge| Edge {
                    curve: edge.curve.mapped(map, linear),
                    ..edge.clone()
                })
                .collect(),
            faces: self
                .faces
                .iter()
                .map(|face| Face {
                    surface: face.surface.mapped(map, linear),
                    ..face.clone()
                })
                .collect(),
            shells: self.shells.clone(),
        }
    }

    /// Edge parameters whose chords stay within `deflection` and turn by
    /// at most `angle`.
    pub(crate) fn sample(&self, deflection: f64, angle: f64) -> Result<Vec<Vec<f64>>> {
        self.edges
            .iter()
            .map(|edge| {
                edge.curve
                    .sample(edge.range[0], edge.range[1], deflection, angle)
            })
            .collect()
    }

    /// Outward sign of each face: surface sense combined with a void shell.
    pub(crate) fn signs(&self) -> Vec<f64> {
        let mut signs = vec![1.; self.faces.len()];
        for shell in &self.shells {
            for &face in &shell.faces {
                let outward = self.faces[face].same_sense != shell.reversed;
                signs[face] = if outward { 1. } else { -1. };
            }
        }
        signs
    }
}

fn find(group: &mut [usize], mut i: usize) -> usize {
    while group[i] != i {
        group[i] = group[group[i]];
        i = group[i];
    }
    i
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FaceProperties {
    pub area_mm2: f64,
    pub centroid_mm: Point,
    /// The face's signed share of the enclosed volume.
    pub flux_mm3: f64,
}

/// A validated solid: one outer shell and any number of void shells.
#[derive(Clone, Debug)]
pub struct Solid {
    boundary: Boundary,
    tolerance: f64,
    domains: Vec<Domain>,
    faces: Vec<FaceProperties>,
    volume_mm3: f64,
    volume_error_mm3: f64,
    bounds_mm: [Point; 2],
}

impl Solid {
    /// Validate the boundary and integrate its properties. `tolerance` is
    /// the largest admitted gap between a vertex, its edges and their faces.
    pub fn new(boundary: Boundary, tolerance: f64, interrupt: Interrupt<'_>) -> Result<Self> {
        ensure(tolerance.is_finite() && tolerance > 0., Error::Invalid)?;
        let mut boundary = boundary;
        let samples = boundary.sample(f64::INFINITY, f64::INFINITY)?;
        for index in 0..boundary.faces.len() {
            checkpoint(interrupt)?;
            let counter_clockwise = domain::orientation(&boundary, index, &samples, tolerance)?;
            let face = &mut boundary.faces[index];
            // A reversed face whose loops alone do not show their direction
            // could be either of its surface's two sides.
            ensure(
                counter_clockwise.is_some() || face.same_sense,
                Error::Unsupported,
            )?;
            if counter_clockwise.is_some_and(|ccw| ccw != face.same_sense) {
                // Keep one convention: loops run about the face normal.
                for trim in &mut face.loops {
                    trim.uses.reverse();
                    for used in &mut trim.uses {
                        used.1 = !used.1;
                    }
                }
            }
        }
        boundary.validate()?;
        let anchor = *boundary.vertices.first().ok_or(Error::Invalid)?;
        let signs = boundary.signs();
        let mut domains = Vec::with_capacity(boundary.faces.len());
        let mut faces = Vec::with_capacity(boundary.faces.len());
        let mut error = 0.;
        for (index, face) in boundary.faces.iter().enumerate() {
            checkpoint(interrupt)?;
            let domain = domain::build(&boundary, index, &samples, tolerance)?;
            let coarse = mass::integrate(&boundary, face, &domain, anchor, 1, tolerance)?;
            let fine = mass::integrate(&boundary, face, &domain, anchor, 2, tolerance)?;
            ensure(fine[1].is_finite() && fine[1] > 0., Error::Invalid)?;
            error += (fine[0] - coarse[0]).abs();
            faces.push(FaceProperties {
                area_mm2: fine[1],
                centroid_mm: [fine[2] / fine[1], fine[3] / fine[1], fine[4] / fine[1]],
                flux_mm3: signs[index] * fine[0],
            });
            domains.push(domain);
        }
        let mut volume_mm3 = 0.;
        for (index, shell) in boundary.shells.iter().enumerate() {
            let enclosed: f64 = shell.faces.iter().map(|&face| faces[face].flux_mm3).sum();
            // The outer shell encloses material; every further shell is a
            // void whose outward normals point into it.
            ensure(
                enclosed.is_finite() && (enclosed > 0.) == (index == 0),
                Error::Invalid,
            )?;
            volume_mm3 += enclosed;
        }
        ensure(volume_mm3 > 0., Error::Invalid)?;
        let area: f64 = faces.iter().map(|face| face.area_mm2).sum();
        let bounds_mm = extent::bounds(&boundary, &domains)?;
        Ok(Self {
            boundary,
            tolerance,
            domains,
            faces,
            volume_mm3,
            volume_error_mm3: error + tolerance * area,
            bounds_mm,
        })
    }

    pub fn boundary(&self) -> &Boundary {
        &self.boundary
    }

    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }

    pub fn volume_mm3(&self) -> f64 {
        self.volume_mm3
    }

    /// Bound on the volume error: the measured quadrature convergence plus
    /// the boundary's positional tolerance swept over its area.
    pub fn volume_error_mm3(&self) -> f64 {
        self.volume_error_mm3
    }

    pub fn bounds_mm(&self) -> [Point; 2] {
        self.bounds_mm
    }

    pub fn faces(&self) -> &[FaceProperties] {
        &self.faces
    }

    /// The same solid under a rigid motion. `map` moves points and `linear`
    /// turns directions. Volume and areas are invariant and are kept.
    pub fn moved(
        &self,
        map: &dyn Fn(Point) -> Point,
        linear: &dyn Fn(Point) -> Point,
    ) -> Result<Self> {
        let boundary = self.boundary.mapped(map, linear);
        let bounds_mm = extent::bounds(&boundary, &self.domains)?;
        Ok(Self {
            faces: self
                .faces
                .iter()
                .map(|face| FaceProperties {
                    centroid_mm: map(face.centroid_mm),
                    ..*face
                })
                .collect(),
            boundary,
            bounds_mm,
            domains: self.domains.clone(),
            ..*self
        })
    }
}
