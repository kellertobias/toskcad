//! Hex stock joined to original finite thread rims. This stage earns source
//! distance, normals and indexed closure; a separate literal intersection
//! proof is required before manufacturing admission.
use super::*;
use crate::{
    display::{DisplayMesh, Tessellation},
    manufacturing::{self, MeshProof},
    radial_tensor::TracedMesh,
};
mod sleeve;
pub use sleeve::{join_sleeve, sleeve_segments};

#[derive(Clone, Copy)]
pub enum ThreadBoundaryMesh<'a, 's> {
    Clipped(&'a StockClipTracedMesh<'s>),
    Retained(&'a TracedMesh<'s>),
}
#[derive(Clone, Copy)]
pub struct MeshOptions {
    pub tessellation: Tessellation,
    /// New joining work, separate from immutable input qualification.
    pub maximum_work: usize,
    pub maximum_bytes: usize,
    /// Other live caller/source state besides the borrowed mesh and trace.
    /// In particular include retained RadialBoundary source allocations here.
    pub additional_resident_bytes: usize,
}
pub struct JoinedMeshCandidate<'a, 's> {
    source: ThreadBoundaryMesh<'a, 's>,
    mesh: DisplayMesh,
    connectivity: MeshProof,
    distance: f64,
    work: usize,
    peak: usize,
}
impl JoinedMeshCandidate<'_, '_> {
    pub fn source(&self) -> ThreadBoundaryMesh<'_, '_> {
        self.source
    }
    pub fn mesh(&self) -> &DisplayMesh {
        &self.mesh
    }
    pub fn into_mesh(self) -> DisplayMesh {
        self.mesh
    }
    pub fn connectivity(&self) -> MeshProof {
        self.connectivity
    }
    pub fn source_distance_mm(&self) -> f64 {
        self.distance
    }
    pub fn joining_work_used(&self) -> usize {
        self.work
    }
    pub fn inherited_source_work_used(&self) -> usize {
        self.source.work()
    }
    /// Borrowed mesh/trace, caller-declared other residents and joining peak.
    /// This is not process RSS or an undeclared whole-scene heap bound.
    pub fn accounted_peak_bytes(&self) -> usize {
        self.peak
    }
}
fn mesh_bytes(mesh: &DisplayMesh) -> usize {
    mesh.vertices.capacity() * 8 + (mesh.indices.capacity() + mesh.face_indices.capacity()) * 4
}
impl<'a, 's> ThreadBoundaryMesh<'a, 's> {
    fn mesh(self) -> &'a DisplayMesh {
        match self {
            Self::Clipped(t) => t.candidate().mesh(),
            Self::Retained(t) => t.mesh(),
        }
    }
    fn work(self) -> usize {
        match self {
            Self::Clipped(t) => t.work_used(),
            Self::Retained(t) => t.work_used(),
        }
    }
    fn bytes(self) -> usize {
        match self {
            Self::Clipped(t) => t.accounted_peak_bytes(),
            Self::Retained(t) => mesh_bytes(t.mesh()) + t.retained_metadata_bytes(),
        }
    }
}
struct Work<'a> {
    used: usize,
    held: usize,
    peak: usize,
    options: MeshOptions,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Work<'_> {
    fn charge(&mut self, n: usize) -> Result<(), Error> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(Error::Cancelled);
        }
        self.used = self
            .used
            .checked_add(n)
            .filter(|n| *n <= self.options.maximum_work)
            .ok_or(Error::Limit)?;
        Ok(())
    }
    fn storage(&mut self, n: usize) -> Result<(), Error> {
        let total = self
            .held
            .checked_add(n)
            .filter(|n| *n <= self.options.maximum_bytes)
            .ok_or(Error::Limit)?;
        self.peak = self.peak.max(total);
        Ok(())
    }
    fn triangle(
        &mut self,
        mesh: &mut DisplayMesh,
        ids: [u32; 3],
        face: u32,
        normal: [I; 3],
    ) -> Result<(), Error> {
        self.charge(256)?;
        if mesh.face_indices.len() >= self.options.tessellation.triangle_limit {
            return Err(Error::Limit);
        }
        let p = ids.map(|id| {
            std::array::from_fn::<_, 3, _>(|k| I::point(mesh.vertices[3 * id as usize + k]))
        });
        let a = std::array::from_fn::<_, 3, _>(|k| p[1][k] - p[0][k]);
        let b = std::array::from_fn::<_, 3, _>(|k| p[2][k] - p[0][k]);
        let cross = std::array::from_fn::<_, 3, _>(|k| {
            a[(k + 1) % 3] * b[(k + 2) % 3] - a[(k + 2) % 3] * b[(k + 1) % 3]
        });
        let length = cross
            .into_iter()
            .map(I::square)
            .fold(I::point(0.), |a, b| a + b)
            .sqrt();
        if !length.lo.is_finite() || !length.hi.is_finite() || length.lo <= 0. {
            return Err(Error::Unresolved);
        }
        let error = (0..3)
            .map(|k| (cross[k] / length - normal[k]).square())
            .fold(I::point(0.), |a, b| a + b)
            .sqrt();
        let allowed =
            (I::point(2.) * (I::point(self.options.tessellation.angle) / I::point(2.)).sin()).lo;
        if !error.hi.is_finite() || error.hi > allowed {
            return Err(Error::Unresolved);
        }
        mesh.indices.extend(ids);
        mesh.face_indices.push(face);
        Ok(())
    }
}
struct Rims {
    vertices: [Vec<u32>; 2],
    caps: [u32; 2],
    height: [f64; 2],
    // Every unchanged source face, including a retained head-tip cap fan,
    // retains the source's proven distance bound.
    // This may be a guaranteed requested quality, not a measured residual.
    unchanged_source_bound: f64,
    // Actual finite original cap-curve to emitted rim chord correspondence.
    distance: f64,
}
fn rims(source: ThreadBoundaryMesh<'_, '_>, work: &mut Work<'_>) -> Result<Rims, Error> {
    match source {
        ThreadBoundaryMesh::Clipped(trace) => {
            let count = trace.candidate().mesh().vertices.len() / 3;
            let mut vertices = [Vec::new(), Vec::new()];
            let unchanged_source_bound = trace.candidate().source_distance_mm();
            let mut distance = 0f64;
            let allowed = (I::point(2.)
                * (I::point(work.options.tessellation.angle) / I::point(2.)).sin())
            .lo;
            if unchanged_source_bound > work.options.tessellation.deflection
                || trace.candidate().normal_chord() > allowed
            {
                return Err(Error::Unsupported);
            }
            work.storage(count * 16 + trace.cap_edges().len() * 8)?;
            let mut next = vec![u32::MAX; count];
            for (end, station) in [StockClipStation::Start, StockClipStation::End]
                .into_iter()
                .enumerate()
            {
                next.fill(u32::MAX);
                let mut first = None;
                let mut n = 0;
                for edge in trace.cap_edges().iter().filter(|e| e.station == station) {
                    work.charge(16)?;
                    first = Some(edge.vertices[0]);
                    next[edge.vertices[0] as usize] = edge.vertices[1];
                    distance = distance.max(edge.chord_error_mm);
                    n += 1;
                }
                // SAME sealed trace gives one original cap cycle; all edges
                // and all endpoint indices were earned by its source builder.
                let first = first.expect("sealed source cap has edges");
                vertices[end].reserve_exact(n);
                let mut current = first;
                for _ in 0..n {
                    work.charge(1)?;
                    vertices[end].push(current);
                    current = next[current as usize];
                }
                assert_eq!(current, first, "sealed original source cap cycle");
                if end == 0 {
                    vertices[end].reverse();
                }
            }
            Ok(Rims {
                vertices,
                caps: [5, 6],
                height: [0., trace.source().source().source().spec().length_mm],
                unchanged_source_bound,
                distance,
            })
        }
        ThreadBoundaryMesh::Retained(trace) => {
            let q = trace.tessellation();
            if q.deflection > work.options.tessellation.deflection
                || q.angle > work.options.tessellation.angle
            {
                return Err(Error::Unsupported);
            }
            let axes = trace.source().cap_axes();
            if axes.iter().any(|a| {
                a.origin[0] != 0.
                    || a.origin[1] != 0.
                    || a.reference != [1., 0., 0.]
                    || a.normal[0] != 0.
                    || a.normal[1] != 0.
                    || a.normal[2].abs() != 1.
            }) {
                return Err(Error::Unsupported);
            }
            let rings = trace.cap_rims();
            work.storage((rings[0].len() + rings[1].len()) * 4)?;
            let mut distance = 0f64;
            let vertices = std::array::from_fn(|end| {
                rings[end]
                    .iter()
                    .map(|p| {
                        distance = distance.max(p.outgoing_chord_error_mm);
                        p.vertex
                    })
                    .collect()
            });
            work.charge(rings[0].len() + rings[1].len())?;
            Ok(Rims {
                vertices,
                caps: rings.map(|r| r[0].cap_face as u32),
                height: axes.map(|a| a.origin[2]),
                unchanged_source_bound: q.deflection,
                distance,
            })
        }
    }
}
fn plane(
    mesh: &mut DisplayMesh,
    outer: &[u32],
    inner: &[u32],
    positive: bool,
    face: u32,
    work: &mut Work<'_>,
) -> Result<(), Error> {
    let n = outer.len() + inner.len();
    work.charge(n * 32)?;
    // Pinned Spade2.15.1 planar workspace: original points only, <=3n
    // edges/2n faces;8192n+fixed1MiB covers its DCEL, constraint/flood paths
    // and temporary vectors conservatively. Input/output held separately.
    work.storage(mesh_bytes(mesh) + 1024 * 1024 + 8192 * (n + 8))?;
    let vertices = outer.iter().chain(inner).copied().collect::<Vec<_>>();
    let points = vertices
        .iter()
        .map(|id| {
            spade::Point2::new(
                mesh.vertices[3 * *id as usize],
                mesh.vertices[3 * *id as usize + 1],
            )
        })
        .collect();
    let edges = (0..outer.len())
        .map(|i| [i, (i + 1) % outer.len()])
        .chain((0..inner.len()).map(|i| [outer.len() + (i + 1) % inner.len(), outer.len() + i]))
        .collect::<Vec<_>>();
    let region = crate::display::triangulate_region(points, &edges, usize::from(!inner.is_empty()))
        .map_err(|_| Error::Unresolved)?;
    for triangle in region.triangles {
        let mut ids = triangle.map(|i| vertices[i]);
        if !positive {
            ids.swap(1, 2);
        }
        work.triangle(
            mesh,
            ids,
            face,
            [
                I::point(0.),
                I::point(0.),
                I::point(if positive { 1. } else { -1. }),
            ],
        )?;
    }
    Ok(())
}
/// Join to the actual qualified thread source, without assigning an authored
/// recipe to independently retained geometry. Cap two-sided C0 follows from
/// its actual paired boundary homotopy; straight hex edges carry their own
/// directed literal-rounding bound. Remaining source/facet normals are checked
/// against the exact planar supports. No manufacturing flag is returned here.
pub fn join<'a, 's>(
    source: ThreadBoundaryMesh<'a, 's>,
    kind: Kind,
    options: MeshOptions,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<JoinedMeshCandidate<'a, 's>, Error> {
    options.tessellation.check().map_err(|_| Error::Invalid)?;
    if options.maximum_work == 0
        || options.maximum_work > 512_000_000
        || options.maximum_bytes == 0
        || options.maximum_bytes > 64 * 1024 * 1024
    {
        return Err(Error::Invalid);
    }
    let (af, head) = match kind {
        Kind::HexHead {
            across_flats_mm,
            height_mm,
        } => (across_flats_mm, Some(height_mm)),
        Kind::HexNut { across_flats_mm } => (across_flats_mm, None),
    };
    if !af.is_finite() || af <= 0. || head.is_some_and(|h| !h.is_finite() || h <= 0.) {
        return Err(Error::Invalid);
    }
    let held = source
        .bytes()
        .checked_add(options.additional_resident_bytes)
        .ok_or(Error::Limit)?;
    let mut work = Work {
        used: 0,
        held,
        peak: 0,
        options,
        interrupt,
    };
    work.storage(0)?;
    work.charge(1)?;
    let mut rings = rims(source, &mut work)?;
    work.held = work
        .held
        .checked_add(
            rings
                .vertices
                .iter()
                .map(|r| r.capacity() * 4)
                .sum::<usize>(),
        )
        .ok_or(Error::Limit)?;
    work.storage(0)?;
    if !rings.height[0].is_finite()
        || !rings.height[1].is_finite()
        || rings.height[0] >= rings.height[1]
    {
        return Err(Error::Unsupported);
    }
    let source_mesh = source.mesh();
    // Both sealed trace providers emit one triangle per cap coedge. Replacing
    // one fan by a hex annulus, six walls and the four-triangle top adds22;
    // replacing both fans by two annuli and six walls adds24. Refuse a caller
    // triangle frontier before copying the already-qualified input arrays.
    let expected_triangles = source_mesh
        .face_indices
        .len()
        .checked_add(if head.is_some() { 22 } else { 24 })
        .ok_or(Error::Limit)?;
    if expected_triangles > options.tessellation.triangle_limit {
        return Err(Error::Limit);
    }
    let mut cap_height_error = 0f64;
    for (end, ring) in rings.vertices.iter().enumerate() {
        for &id in ring {
            work.charge(1)?;
            cap_height_error = cap_height_error.max(
                (I::point(source_mesh.vertices[3 * id as usize + 2]) - I::point(rings.height[end]))
                    .magnitude(),
            );
        }
    }
    // Declared sewn source vertices need not lie exactly on the cap support.
    // Preserve those SAME literal coordinates and charge their full plane lift
    // in the annulus transfer; every new actual facet normal is still checked.
    rings.distance = (I::point(rings.distance) + I::point(cap_height_error)).hi;
    work.charge(4)?;
    let containment_distance = rings.unchanged_source_bound.max(rings.distance);
    let root3 = I::point(3.).sqrt();
    let half = I::point(af) / I::point(2.);
    let radius = I::point(af) / root3;
    let zero = I::point(0.);
    let hex = [
        [radius, zero],
        [radius / I::point(2.), half],
        [-radius / I::point(2.), half],
        [-radius, zero],
        [-radius / I::point(2.), -half],
        [radius / I::point(2.), -half],
    ];
    for p in source_mesh.vertices.as_chunks::<3>().0 {
        work.charge(16)?;
        // Convexity covers complete literal triangles. The actual complete
        // source is within the inherited Euclidean C0 tube; these halfspaces
        // have gradient norms2 and1 respectively. Thus this strict margin
        // contains the complete actual thread source, not merely its vertices.
        if (root3 * I::point(p[0].abs())
            + I::point(p[1].abs())
            + I::point(2.) * I::point(containment_distance))
        .hi >= I::point(af).lo
            || (I::point(p[1].abs()) + I::point(containment_distance)).hi >= half.lo
        {
            return Err(Error::Unsupported);
        }
    }
    let lower = if head.is_some() {
        I::point(rings.height[1])
    } else {
        I::point(rings.height[0])
    };
    let upper = if let Some(h) = head {
        I::point(rings.height[1]) + I::point(h)
    } else {
        I::point(rings.height[1])
    };
    if !lower.lo.is_finite() || !upper.hi.is_finite() || lower.hi >= upper.lo {
        return Err(Error::Unsupported);
    }
    let extra_triangles =
        source_mesh.face_indices.len() + rings.vertices[0].len() + rings.vertices[1].len() + 40;
    let count = source_mesh.vertices.len() / 3 + 12;
    work.storage(
        count * 24 + extra_triangles * 16 + count * 4 + 52 * extra_triangles + 1024 * 1024,
    )?;
    let mut mesh = DisplayMesh {
        vertices: Vec::with_capacity(3 * count),
        indices: Vec::with_capacity(3 * extra_triangles),
        face_indices: Vec::with_capacity(extra_triangles),
    };
    mesh.vertices.extend_from_slice(&source_mesh.vertices);
    for (triangle, &face) in source_mesh
        .indices
        .as_chunks::<3>()
        .0
        .iter()
        .zip(&source_mesh.face_indices)
    {
        work.charge(16)?;
        if face == rings.caps[1] || (head.is_none() && face == rings.caps[0]) {
            continue;
        }
        let mut ids = *triangle;
        if head.is_none() {
            ids.swap(1, 2);
        }
        mesh.indices.extend(ids);
        mesh.face_indices.push(face);
    }
    let base = mesh.vertices.len() / 3;
    let mut rounding = 0f64;
    for z in [lower, upper] {
        for xy in hex {
            work.charge(32)?;
            let p = [xy[0], xy[1], z];
            let literal = p.map(|v| v.lo + (v.hi - v.lo) * 0.5);
            if literal.iter().any(|v| !v.is_finite()) {
                return Err(Error::Unsupported);
            }
            let error = (0..3)
                .map(|k| (I::point(literal[k]) - p[k]).square())
                .fold(zero, |a, b| a + b)
                .sqrt()
                .hi;
            if !error.is_finite() {
                return Err(Error::Unsupported);
            }
            rounding = rounding.max(error);
            mesh.vertices.extend(literal);
        }
    }
    // The retained lateral triangles are unchanged. New planar hex/cap
    // patches alone incur the actual rim/plane/literal-hex transfer. A bound
    // for the complete disjoint face union is their MAXIMUM, not adding new
    // cap error to every unchanged lateral face's already-earned quality.
    let cap_distance = (I::point(rings.distance) + I::point(rounding)).hi;
    let distance = rings.unchanged_source_bound.max(cap_distance).max(rounding);
    if distance > options.tessellation.deflection {
        return Err(Error::Unsupported);
    }
    let face_base = source_mesh
        .face_indices
        .iter()
        .copied()
        .max()
        .expect("qualified mesh has faces")
        + 1;
    for i in 0..6 {
        let j = (i + 1) % 6;
        let dx = hex[j][0] - hex[i][0];
        let dy = hex[j][1] - hex[i][1];
        let length = (dx.square() + dy.square()).sqrt();
        let normal = [dy / length, -dx / length, zero];
        let a = (base + i) as u32;
        let b = (base + j) as u32;
        let c = b + 6;
        let d = a + 6;
        work.triangle(&mut mesh, [a, b, c], face_base + i as u32, normal)?;
        work.triangle(&mut mesh, [a, c, d], face_base + i as u32, normal)?;
    }
    let outer0 = (base as u32..base as u32 + 6).collect::<Vec<_>>();
    let outer1 = (base as u32 + 6..base as u32 + 12).collect::<Vec<_>>();
    if head.is_some() {
        plane(
            &mut mesh,
            &outer0,
            &rings.vertices[1],
            false,
            face_base + 6,
            &mut work,
        )?;
        plane(&mut mesh, &outer1, &[], true, face_base + 7, &mut work)?;
    } else {
        plane(
            &mut mesh,
            &outer0,
            &rings.vertices[0],
            false,
            face_base + 6,
            &mut work,
        )?;
        plane(
            &mut mesh,
            &outer1,
            &rings.vertices[1],
            true,
            face_base + 7,
            &mut work,
        )?;
    }
    // Remove former cap apices only by actual index use; no coordinate weld.
    work.storage(mesh_bytes(&mesh) * 2 + count * 4 + 52 * mesh.face_indices.len())?;
    let mut remap = vec![u32::MAX; count];
    let mut compact = DisplayMesh {
        vertices: Vec::with_capacity(mesh.vertices.len()),
        indices: Vec::with_capacity(mesh.indices.len()),
        face_indices: mesh.face_indices.clone(),
    };
    let old_mesh_bytes = mesh_bytes(&mesh);
    for id in mesh.indices {
        work.charge(1)?;
        let new = &mut remap[id as usize];
        if *new == u32::MAX {
            *new = (compact.vertices.len() / 3) as u32;
            compact
                .vertices
                .extend_from_slice(&mesh.vertices[3 * id as usize..3 * id as usize + 3]);
        }
        compact.indices.push(*new);
    }

    assert_eq!(
        compact.face_indices.len(),
        expected_triangles,
        "sealed original cap fan and planar region triangle counts"
    );
    let resident = old_mesh_bytes + mesh_bytes(&compact) + remap.capacity() * 4;
    let connectivity = manufacturing::certify_counted(&compact, &mut |n, scratch| {
        if let Some(scratch) = scratch {
            work.storage(resident.checked_add(scratch).ok_or(Error::Limit)?)?;
        }
        work.charge(n)
    })
    .map_err(|error| match error {
        manufacturing::CountedError::Mesh(_) => Error::Unresolved,
        manufacturing::CountedError::Interrupted(error) => error,
    })?;
    Ok(JoinedMeshCandidate {
        source,
        mesh: compact,
        connectivity,
        distance,
        work: work.used,
        peak: work.peak,
    })
}

#[cfg(test)]
mod tests;
