//! The retained boundary of one solid and its integrated properties.
use super::domain::{self, Domain};
use super::{
    Curve, Error, Interrupt, Point, Result, Surface, checkpoint, distance, ensure, extent, mass,
};
use super::{crosses_ray, mesh, point_triangle};
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

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
        for (side, &point) in points.iter().enumerate() {
            let on_curve = edge.curve.evaluate(edge.range[side]).0;
            ensure(distance(on_curve, point) <= tolerance, Error::Invalid)?;
        }
        Ok(edge)
    }
}

impl Edge {
    /// Points along the edge whose chords stay within `deflection` and
    /// turn by at most `angle`.
    pub fn polyline(&self, deflection: f64, angle: f64) -> Result<Vec<Point>> {
        let parameters = self
            .curve
            .sample(self.range[0], self.range[1], deflection, angle)?;
        Ok(parameters
            .into_iter()
            .map(|t| self.curve.evaluate(t).0)
            .collect())
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

/// Vertices on the loops of one face.
fn face_vertices(boundary: &Boundary, face: usize) -> Vec<usize> {
    let mut vertices = Vec::new();
    for trim in &boundary.faces[face].loops {
        vertices.extend(trim.vertex);
        for &(edge, _) in &trim.uses {
            vertices.extend(boundary.edges[edge].vertices);
        }
    }
    vertices
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

/// Check earned boundary incidence without returning geometry. Non-analytic
/// support orientations remain unresolved, so topology that depends on them is
/// not classified. No area, volume, extent or nesting permission is returned.
pub(crate) fn validate_incidence(
    boundary: Boundary,
    tolerance: f64,
    interrupt: Interrupt<'_>,
) -> Result<()> {
    prepare_boundary(boundary, tolerance, interrupt, Preparation::Refusal).map(|_| ())
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Preparation {
    Full,
    Refusal,
}

fn prepare_boundary(
    boundary: Boundary,
    tolerance: f64,
    interrupt: Interrupt<'_>,
    mode: Preparation,
) -> Result<(Boundary, Vec<Vec<f64>>)> {
    ensure(tolerance.is_finite() && tolerance > 0., Error::Invalid)?;
    let mut boundary = boundary;
    let samples = boundary.sample(f64::INFINITY, f64::INFINITY)?;
    let mut normalized = true;
    for index in 0..boundary.faces.len() {
        checkpoint(interrupt)?;
        if mode == Preparation::Refusal
            && matches!(
                &boundary.faces[index].surface,
                Surface::Nurbs { .. } | Surface::Extrusion { .. } | Surface::Revolution { .. }
            )
        {
            domain::incidence_only(&boundary, index, &samples, tolerance)?;
            normalized = false;
            continue;
        }
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
    // Edge-use senses may only be checked after ALL actual face loops
    // were normalized. A skipped orientation grants no topology authority.
    if normalized {
        boundary.validate()?;
    }
    // A lone vertex bounding a face, as at a sphere's pole, lies on it.
    for face in &boundary.faces {
        for vertex in face.loops.iter().filter_map(|trim| trim.vertex) {
            let point = *boundary.vertices.get(vertex).ok_or(Error::Invalid)?;
            let uv = face.surface.invert(point, None);
            ensure(
                distance(face.surface.evaluate(uv).0, point) <= tolerance,
                Error::Invalid,
            )?;
        }
    }
    Ok((boundary, samples))
}

/// Faces below which one thread integrates them all.
pub(super) const PARALLEL_FACES: usize = 8;

/// Threads that integrate `count` faces: one below `2 * PARALLEL_FACES`,
/// otherwise up to eight.
fn face_workers(count: usize) -> usize {
    std::thread::available_parallelism()
        .map_or(1, usize::from)
        .min(8)
        .min(count / PARALLEL_FACES.max(1))
}

/// `f` of every face index, computed on `workers` threads. Each result is
/// stored at its index, so the caller combines them in the same order as a
/// sequential loop and every face's own arithmetic is unchanged: the
/// results are bit-identical whatever the thread count. The interrupt is
/// polled on the calling thread before each face it takes, as the
/// sequential loop did. Once a face fails, faces not yet started stay
/// `None`. Once the interrupt fires, faces not yet started stay `None` on
/// one thread; on several, helpers may already have finished every face,
/// so all results are dropped and the caller reports cancellation.
pub(super) fn each_face<T: Send>(
    count: usize,
    workers: usize,
    interrupt: Interrupt<'_>,
    f: &(dyn Fn(usize) -> Result<T> + Sync),
) -> Vec<Option<Result<T>>> {
    let mut output: Vec<Option<Result<T>>> = (0..count).map(|_| None).collect();
    if workers < 2 {
        for (index, slot) in output.iter_mut().enumerate() {
            if checkpoint(interrupt).is_err() {
                break;
            }
            let result = f(index);
            let failed = result.is_err();
            *slot = Some(result);
            if failed {
                break;
            }
        }
        return output;
    }
    let next = AtomicUsize::new(0);
    let stop = AtomicBool::new(false);
    let take = || {
        (!stop.load(Ordering::Relaxed))
            .then(|| next.fetch_add(1, Ordering::Relaxed))
            .filter(|&index| index < count)
    };
    // The caller stops at the first failed face, and every lower index
    // was handed out before it, so no further face is needed.
    let run = |index: usize| {
        let result = f(index);
        if result.is_err() {
            stop.store(true, Ordering::Relaxed);
        }
        (index, result)
    };
    let work = || {
        let mut done = Vec::new();
        while let Some(index) = take() {
            done.push(run(index));
        }
        done
    };
    if checkpoint(interrupt).is_err() {
        return output;
    }
    let cancelled = AtomicBool::new(false);
    let finished = std::thread::scope(|scope| {
        let helpers: Vec<_> = (1..workers).map(|_| scope.spawn(work)).collect();
        let mut done = Vec::new();
        loop {
            if checkpoint(interrupt).is_err() {
                cancelled.store(true, Ordering::Relaxed);
                stop.store(true, Ordering::Relaxed);
                break;
            }
            let Some(index) = take() else { break };
            done.push(run(index));
        }
        for helper in helpers {
            match helper.join() {
                Ok(more) => done.extend(more),
                Err(panic) => std::panic::resume_unwind(panic),
            }
        }
        done
    });
    if cancelled.into_inner() {
        return output;
    }
    for (index, result) in finished {
        output[index] = Some(result);
    }
    output
}

impl Solid {
    /// Validate the boundary and integrate its properties. `tolerance` is
    /// the largest admitted gap between a vertex, its edges and their faces.
    pub fn new(boundary: Boundary, tolerance: f64, interrupt: Interrupt<'_>) -> Result<Self> {
        let workers = face_workers(boundary.faces.len());
        Self::build(boundary, tolerance, interrupt, workers)
    }

    /// [`Solid::new`] integrating its faces on `workers` threads.
    pub(super) fn build(
        boundary: Boundary,
        tolerance: f64,
        interrupt: Interrupt<'_>,
        workers: usize,
    ) -> Result<Self> {
        let (boundary, samples) =
            prepare_boundary(boundary, tolerance, interrupt, Preparation::Full)?;
        let anchor = *boundary.vertices.first().ok_or(Error::Invalid)?;
        let signs = boundary.signs();
        let mut domains = Vec::with_capacity(boundary.faces.len());
        let mut faces = Vec::with_capacity(boundary.faces.len());
        let mut error = 0.;
        let integrate = |index: usize| {
            let face = &boundary.faces[index];
            domain::build(&boundary, index, &samples, tolerance).and_then(|domain| {
                let coarse = mass::integrate(&boundary, face, &domain, anchor, 1, tolerance)?;
                let fine = mass::integrate(&boundary, face, &domain, anchor, 2, tolerance)?;
                ensure(fine[1].is_finite() && fine[1] > 0., Error::Invalid)?;
                Ok((domain, coarse, fine))
            })
        };
        let mut integrated = each_face(boundary.faces.len(), workers, interrupt, &integrate);
        for index in 0..boundary.faces.len() {
            // A face left unintegrated was stopped by the interrupt.
            let Some(integrated) = integrated[index].take() else {
                return Err(Error::Cancelled);
            };
            let (domain, coarse, fine) = match integrated {
                // A face whose edges coincide within the tolerance bounds no
                // region. Writers emit such slivers between faces that meet
                // within their uncertainty; they are not invalid, but the
                // boundary would need them merged away.
                Err(Error::Invalid)
                    if domain::degenerate(&boundary, index, &samples, tolerance)? =>
                {
                    return Err(Error::Unsupported);
                }
                other => other?,
            };
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
        let solid = Self {
            boundary,
            tolerance,
            domains,
            faces,
            volume_mm3,
            volume_error_mm3: error + tolerance * area,
            bounds_mm,
        };
        solid.nest(interrupt)?;
        Ok(solid)
    }

    /// Voids lie strictly inside the exterior shell and outside each other.
    /// Shells share no vertex, and every vertex of a void is classified
    /// against the other shells on a tessellation fine against the part:
    /// a vertex outside where it belongs is invalid, and one closer to
    /// another shell than that tessellation resolves is refused as
    /// unsupported rather than guessed.
    fn nest(&self, interrupt: Interrupt<'_>) -> Result<()> {
        let boundary = &self.boundary;
        if boundary.shells.len() < 2 {
            return Ok(());
        }
        let mut owner = vec![usize::MAX; boundary.vertices.len()];
        let mut face_shell = vec![0; boundary.faces.len()];
        for (shell_index, shell) in boundary.shells.iter().enumerate() {
            for &face in &shell.faces {
                face_shell[face] = shell_index;
                for vertex in face_vertices(boundary, face) {
                    let slot = &mut owner[vertex];
                    ensure(*slot == usize::MAX || *slot == shell_index, Error::Invalid)?;
                    *slot = shell_index;
                }
            }
        }
        let diagonal = distance(self.bounds_mm[0], self.bounds_mm[1]);
        let deflection = (1e-4 * diagonal).max(10. * self.tolerance);
        let options = mesh::Options {
            deflection_mm: deflection,
            angle_rad: 0.5,
            triangle_limit: 2_000_000,
        };
        let tessellation = mesh::tessellate(self, options, interrupt)?;
        let triangles: Vec<(usize, [Point; 3])> = tessellation
            .triangles
            .iter()
            .zip(&tessellation.faces)
            .map(|(triangle, &face)| {
                (
                    face_shell[face as usize],
                    triangle.map(|i| tessellation.positions[i as usize]),
                )
            })
            .collect();
        for (vertex, &shell) in owner.iter().enumerate() {
            if shell == 0 || shell == usize::MAX {
                continue;
            }
            checkpoint(interrupt)?;
            let point = boundary.vertices[vertex];
            for other in (0..boundary.shells.len()).filter(|&other| other != shell) {
                let near = triangles
                    .iter()
                    .filter(|(owner, _)| *owner == other)
                    .map(|(_, triangle)| point_triangle(point, *triangle))
                    .fold(f64::INFINITY, f64::min);
                ensure(near > deflection, Error::Unsupported)?;
                let crossings = triangles
                    .iter()
                    .filter(|(owner, triangle)| *owner == other && crosses_ray(point, triangle))
                    .count();
                // Inside the exterior, outside every other void.
                ensure((crossings % 2 == 1) == (other == 0), Error::Invalid)?;
            }
        }
        Ok(())
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

    /// Axis-aligned bounds of each trimmed face, in face order.
    pub fn face_bounds_mm(&self) -> Result<Vec<[Point; 2]>> {
        extent::face_bounds(&self.boundary, &self.domains)
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

#[cfg(test)]
mod refusal_tests;
