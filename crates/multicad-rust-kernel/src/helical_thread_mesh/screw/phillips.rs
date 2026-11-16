//! Actual indexed cavity sewing. Global literal embedding remains separate.
use super::*;
use crate::{
    manufacturing::{self, MeshProof},
    phillips_drive::{Error, FaceKey, PhillipsScrew, ScrewFace},
    tapered_thread::interval::I,
};
use std::collections::{BTreeMap, BTreeSet};
#[derive(Clone, Copy, Debug)]
pub enum PhillipsThreadVertexKind {
    Grid {
        fraction: f64,
        phase: f64,
        period: i64,
    },
    CapColumn {
        fraction: f64,
        cap: usize,
        rail: Option<(i64, f64)>,
    },
    CapRail {
        phase: f64,
        period: i64,
        cap: usize,
    },
}
#[derive(Clone, Copy, Debug)]
pub struct PhillipsThreadVertex {
    pub vertex: u32,
    pub kind: PhillipsThreadVertexKind,
}
#[derive(Clone, Copy, Debug)]
pub struct PhillipsHeadVertex {
    pub vertex: u32,
    pub segment: usize,
    pub parameter: f64,
    pub angle: f64,
}
#[derive(Clone, Copy, Debug)]
pub struct PhillipsCavityVertex {
    pub vertex: u32,
    pub original: [f64; 3],
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PhillipsMeshError {
    Geometry(MeshError),
    Cancelled,
}
impl From<MeshError> for PhillipsMeshError {
    fn from(e: MeshError) -> Self {
        Self::Geometry(e)
    }
}
fn source_error(e: Error) -> PhillipsMeshError {
    match e {
        Error::Cancelled => PhillipsMeshError::Cancelled,
        Error::Limit => MeshError::BudgetExceeded.into(),
        Error::Invalid | Error::Unresolved => MeshError::UnresolvedBoundary.into(),
    }
}
fn cancel(f: &dyn Fn() -> bool) -> Result<(), PhillipsMeshError> {
    if f() {
        Err(PhillipsMeshError::Cancelled)
    } else {
        Ok(())
    }
}
/// Closed source candidate; connectivity and face ownership are certified.
/// This is not a native manufacturing receipt or independent file mass proof.
#[derive(Clone, Debug)]
pub struct PhillipsScrewMesh {
    mesh: DisplayMesh,
    connectivity: MeshProof,
    faces: Vec<ScrewFace>,
    source_volume: crate::phillips_drive::SourceVolumeBounds,
    literal_volume_mm3: [f64; 2],
    volume_difference_bound_mm3: f64,
    work: usize,
    handoffs: [usize; 3],
    thread_vertices: Vec<PhillipsThreadVertex>,
    head_vertices: Vec<PhillipsHeadVertex>,
    cavity_vertices: Vec<PhillipsCavityVertex>,
    cavity_distance: f64,
    proposal: PhillipsScrew,
}
impl PhillipsScrewMesh {
    pub fn proposal(&self) -> &PhillipsScrew {
        &self.proposal
    }
    pub fn retained_budget_bytes(&self) -> usize {
        std::mem::size_of::<Self>()
            + self.mesh.vertices.capacity() * 8
            + (self.mesh.indices.capacity() + self.mesh.face_indices.capacity()) * 4
            + self.faces.capacity() * std::mem::size_of::<ScrewFace>()
            + self.thread_vertices.capacity() * std::mem::size_of::<PhillipsThreadVertex>()
            + self.head_vertices.capacity() * std::mem::size_of::<PhillipsHeadVertex>()
            + self.cavity_vertices.capacity() * std::mem::size_of::<PhillipsCavityVertex>()
            + 65536
    }
    pub fn thread_vertices(&self) -> &[PhillipsThreadVertex] {
        &self.thread_vertices
    }
    pub fn head_vertices(&self) -> &[PhillipsHeadVertex] {
        &self.head_vertices
    }
    pub fn cavity_vertices(&self) -> &[PhillipsCavityVertex] {
        &self.cavity_vertices
    }
    pub fn cavity_source_distance_mm(&self) -> f64 {
        self.cavity_distance
    }
    /// New bounded-path logical cell, trace, copy and proof units. This is
    /// neither a CPU instruction count nor an allocator/RSS measurement.
    pub fn work_used(&self) -> usize {
        self.work
    }
    /// Cumulative logical debits after Cutter refinement, before planar sewing,
    /// and before connectivity/integration. They share the same total ledger.
    pub fn work_handoffs(&self) -> [usize; 3] {
        self.handoffs
    }
    pub fn mesh(&self) -> &DisplayMesh {
        &self.mesh
    }
    pub fn connectivity(&self) -> MeshProof {
        self.connectivity
    }
    pub fn faces(&self) -> &[ScrewFace] {
        &self.faces
    }
    /// Source-only analytic interval paired during this candidate's creation.
    pub fn source_volume(&self) -> crate::phillips_drive::SourceVolumeBounds {
        self.source_volume
    }
    /// Directed signed-volume evaluation of the actual literal triangles.
    pub fn literal_volume_mm3(&self) -> [f64; 2] {
        self.literal_volume_mm3
    }
    pub fn volume_difference_bound_mm3(&self) -> f64 {
        self.volume_difference_bound_mm3
    }
}
// The caller has just certified this owned mesh's finite coordinates, indices
// and closure. Source dimensions and the five-million-facet resource ceiling
// bound these local cubic terms far below floating overflow. A local literal
// vertex is the origin, reducing cancellation without moving any geometry.
fn literal_volume(mesh: &DisplayMesh) -> I {
    let origin = mesh.vertices.as_chunks::<3>().0[0];
    let mut volume = I::point(0.);
    for tri in mesh.indices.as_chunks::<3>().0 {
        let [a, b, c] = tri.map(|id| {
            std::array::from_fn::<_, 3, _>(|i| I::point(point(mesh, id)[i]) - I::point(origin[i]))
        });
        let cross = [
            b[1] * c[2] - b[2] * c[1],
            b[2] * c[0] - b[0] * c[2],
            b[0] * c[1] - b[1] * c[0],
        ];
        volume = volume + (a[0] * cross[0] + a[1] * cross[1] + a[2] * cross[2]) / I::point(6.);
    }
    volume
}
fn point(mesh: &DisplayMesh, id: u32) -> [f64; 3] {
    mesh.vertices.as_chunks::<3>().0[id as usize]
}
fn translated_normal(
    old: [[f64; 3]; 3],
    new: [[f64; 3]; 3],
    reserve: f64,
) -> Result<(), MeshError> {
    let normal = |p: [[f64; 3]; 3]| {
        let edges = [1, 2]
            .map(|j| std::array::from_fn::<_, 3, _>(|i| I::point(p[j][i]) - I::point(p[0][i])));
        std::array::from_fn::<_, 3, _>(|i| {
            edges[0][(i + 1) % 3] * edges[1][(i + 2) % 3]
                - edges[0][(i + 2) % 3] * edges[1][(i + 1) % 3]
        })
    };
    let a = normal(old);
    let b = normal(new);
    let dot = (0..3).fold(I::point(0.), |s, i| s + a[i] * b[i]);
    let cross = (0..3).fold(I::point(0.), |s, i| {
        s + I::point(
            (a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3]).magnitude(),
        )
    });
    let ratio = cross / dot;
    if dot.lo <= 0. || !ratio.valid() || ratio.hi > reserve {
        return Err(MeshError::UnresolvedBoundary);
    }
    Ok(())
}
// Preserve the source cap's exact edge uses, reversing only its boundary
// orientation to obtain the cavity's CCW inner ring. Never sort/weld points.
fn mouth_ring(
    mesh: &DisplayMesh,
    cap: u32,
    source_order: crate::phillips_drive::MouthSweepProof,
) -> Result<Vec<u32>, MeshError> {
    let mut uses = BTreeMap::<(u32, u32), Vec<(u32, u32)>>::new();
    for (tri, &face) in mesh
        .indices
        .as_chunks::<3>()
        .0
        .iter()
        .zip(&mesh.face_indices)
    {
        if face != cap {
            continue;
        }
        for (a, b) in [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
            uses.entry((a.min(b), a.max(b))).or_default().push((a, b));
        }
    }
    let mut outgoing = BTreeMap::new();
    for edges in uses.values() {
        match *edges.as_slice() {
            [(a, b)] => {
                if outgoing.insert(a, b).is_some() {
                    return Err(MeshError::UnresolvedBoundary);
                }
            }
            [(a, b), (c, d)] if a == d && b == c => {}
            _ => return Err(MeshError::UnresolvedBoundary),
        }
    }
    let start = *outgoing
        .keys()
        .next()
        .ok_or(MeshError::UnresolvedBoundary)?;
    let count = outgoing.len();
    let mut ring = Vec::with_capacity(count);
    let mut at = start;
    loop {
        ring.push(at);
        at = outgoing.remove(&at).ok_or(MeshError::UnresolvedBoundary)?;
        if at == start {
            break;
        }
    }
    if !outgoing.is_empty() || ring.len() < source_order.source_edge_count() {
        return Err(MeshError::UnresolvedBoundary);
    }
    ring.reverse();
    let angle = |id| {
        let [x, y, _] = point(mesh, id);
        y.atan2(x).rem_euclid(TAU)
    };
    let first = (0..ring.len())
        .min_by(|&a, &b| angle(ring[a]).total_cmp(&angle(ring[b])))
        .ok_or(MeshError::UnresolvedBoundary)?;
    ring.rotate_left(first);
    if angle(ring[0]) != 0. {
        return Err(MeshError::UnresolvedBoundary);
    }
    for i in 0..ring.len() {
        let [x, y, z] = point(mesh, ring[i]);
        let [u, v, w] = point(mesh, ring[(i + 1) % ring.len()]);
        let cross = I::point(x) * I::point(v) - I::point(y) * I::point(u);
        let dot = I::point(x) * I::point(u) + I::point(y) * I::point(v);
        if z != 0.
            || w != 0.
            || cross.lo <= 0.
            || dot.lo <= 0.
            || (i + 1 < ring.len() && angle(ring[i]) >= angle(ring[i + 1]))
        {
            return Err(MeshError::UnresolvedBoundary);
        }
    }
    Ok(ring)
}
fn top_triangle(shell: &mut HeadMesh<'_, '_>, tri: [u32; 3], face: u32) -> Result<(), MeshError> {
    let [a, b, c] = tri.map(|id| point(&shell.builder.mesh, id));
    let dx = I::point(b[0]) - I::point(a[0]);
    let dy = I::point(b[1]) - I::point(a[1]);
    let ux = I::point(c[0]) - I::point(a[0]);
    let uy = I::point(c[1]) - I::point(a[1]);
    // Positive projected triangles plus the disjoint ordered source rings
    // prevent a planar fold; never repair a failed zipper by winding flips.
    if a[2] != b[2] || a[2] != c[2] || (dx * uy - dy * ux).lo <= 0. {
        return Err(MeshError::UnresolvedBoundary);
    }
    shell.triangle(tri, face, true)
}
fn top_annulus(
    shell: &mut HeadMesh<'_, '_>,
    inner: &[u32],
    outer: &[u32],
    face: u32,
    cancelled: &dyn Fn() -> bool,
) -> Result<(), PhillipsMeshError> {
    use spade::Point2;
    cancel(cancelled)?;
    let vertices: Vec<_> = outer.iter().chain(inner).copied().collect();
    let points: Vec<_> = vertices
        .iter()
        .map(|&id| {
            let [x, y, _] = point(&shell.builder.mesh, id);
            Point2::new(x, y)
        })
        .collect();
    let mut edges = Vec::with_capacity(vertices.len());
    for i in 0..outer.len() {
        edges.push([i, (i + 1) % outer.len()]);
    }
    // The constrained region has material on its left: outer CCW, hole CW.
    for i in 0..inner.len() {
        edges.push([outer.len() + (i + 1) % inner.len(), outer.len() + i]);
    }
    // This existing checked region API preserves input vertex indices and
    // literal coordinates; duplicates/intersecting constraints are refused.
    // Prospective planar point pairs bound constrained insertion/flood units;
    // reserve before delegating to the existing literal-index CDT contract.
    shell.builder.charge(
        vertices
            .len()
            .checked_mul(vertices.len())
            .and_then(|n| n.checked_mul(16))
            .ok_or(MeshError::BudgetExceeded)?,
    )?;
    let region = crate::display::triangulate_region(points, &edges, 1)?;
    for tri in region.triangles {
        cancel(cancelled)?;
        top_triangle(shell, tri.map(|id| vertices[id]), face)?;
    }
    Ok(())
}
pub fn tessellate_phillips_screw(
    screw: &PhillipsScrew,
    options: Tessellation,
    cancelled: &dyn Fn() -> bool,
) -> Result<PhillipsScrewMesh, PhillipsMeshError> {
    tessellate_phillips_screw_with_work(screw, options, 512_000_000, cancelled)
}
/// Same actual indexed producer with one caller-lowerable logical work ledger.
/// Counts proposed phase cells (also skipped finite cells), traced Cutter
/// refinement, emitted facets, prospective planar pairs and connectivity input.
/// Existing triangle and delegated allocator contracts remain unchanged.
pub fn tessellate_phillips_screw_with_work(
    screw: &PhillipsScrew,
    options: Tessellation,
    maximum_work: usize,
    cancelled: &dyn Fn() -> bool,
) -> Result<PhillipsScrewMesh, PhillipsMeshError> {
    let work = super::super::Counted {
        remaining: std::cell::Cell::new(maximum_work),
        cancelled: std::cell::Cell::new(false),
        interrupt: cancelled,
    };
    if maximum_work > 512_000_000 {
        return Err(MeshError::InvalidOptions.into());
    }
    let result = phillips_counted(screw, options, cancelled, &work, maximum_work);
    if work.cancelled.get() {
        Err(PhillipsMeshError::Cancelled)
    } else {
        result
    }
}
fn phillips_counted(
    screw: &PhillipsScrew,
    options: Tessellation,
    cancelled: &dyn Fn() -> bool,
    work: &super::super::Counted<'_>,
    maximum_work: usize,
) -> Result<PhillipsScrewMesh, PhillipsMeshError> {
    check_options(options)?;
    work.step(32768)?;
    cancel(cancelled)?;
    let source_volume = screw.source_volume_bounds().map_err(source_error)?;
    let source = screw.cavity().source();
    let witness = screw.cavity().mouth_sweep();
    // Cavity volume admission proves the shared positive front intercepts;
    // this sealed witness accompanies the subsequent literal ring checks.
    let extent = source
        .bounding_box_mm()
        .into_iter()
        .flatten()
        .map(f64::abs)
        .fold(1., f64::max)
        .max(screw.thread().spec().length_mm);
    let reserve = (I::point(4.) * I::point(screw.profile().rounding_error_mm())
        + I::point(64. * f64::EPSILON) * I::point(extent))
    .hi;
    if options.deflection <= reserve {
        return Err(MeshError::UnresolvedBoundary.into());
    }
    let working = Tessellation {
        deflection: (options.deflection - reserve) / 2.,
        angle: options.angle * 0.99,
        ..options
    };
    let radius = screw.profile().support([1., 0., 0.]);
    let step = (working.angle * 0.5).min((working.deflection / radius).sqrt());
    let count = (TAU / (12. * step)).ceil() * 12.;
    if count * 2. > working.triangle_limit as f64 {
        return Err(MeshError::BudgetExceeded.into());
    }
    work.step(32768)?;
    let atlas = source
        .surface_patches_with_cancel(cancelled)
        .map_err(source_error)?;
    let traced = atlas
        .boundary_mesh_with_mouth_trace(
            working,
            crate::phillips_drive::MouthTraceOptions {
                maximum_work: work.remaining.get(),
                maximum_bytes: 32 * 1024 * 1024,
            },
            cancelled,
        )
        .map_err(source_error)?;
    work.step(traced.work_used())?;
    let after_cavity = maximum_work - work.remaining.get();
    let cavity = traced.boundary_mesh();
    work.step(
        cavity
            .mesh()
            .indices
            .len()
            .checked_mul(32)
            .ok_or(MeshError::BudgetExceeded)?,
    )?;
    let topology = source.boundary_topology().map_err(source_error)?;
    let cap = topology
        .faces()
        .iter()
        .position(|f| f.source_identity() == FaceKey::Mouth)
        .expect("qualified source mouth") as u32;
    let ring = mouth_ring(cavity.mesh(), cap, witness)?;
    let strips =
        super::super::phase_strips_counted(screw.thread(), TAU / count, working, Some(work))?;
    let mut builder = Builder::new(screw.thread(), count as usize, strips, working);
    builder.work = Some(work);
    builder.lower = screw.visible_thread_interval_mm()[0];
    builder.insert_cap_stations()?;
    builder.build_wall()?;
    cancel(cancelled)?;
    let [mut shoulder, tip] = std::mem::take(&mut builder.caps);
    builder.add_cap(1, tip)?;
    Builder::sort_ring(&mut shoulder)?;
    work.step(shoulder.len().saturating_mul(32))?;
    let faces = screw.faces();
    let face = |key| {
        faces
            .iter()
            .position(|&f| f == key)
            .expect("sealed Phillips face") as u32
    };
    let mut shell = HeadMesh {
        builder: &mut builder,
        angles: shoulder.iter().map(|&(angle, _)| angle).collect(),
    };
    let shaft: Vec<_> = shoulder.into_iter().map(|(_, id)| id).collect();
    let mut head_vertices = Vec::new();
    let outer = shell.profile_outer_with_witnesses(
        screw.profile(),
        screw.spec().shank_length_mm,
        screw.visible_thread_interval_mm()[0],
        &shaft,
        (face(ScrewFace::ThreadShoulder), |i| {
            face(ScrewFace::ProfileWall(i))
        }),
        working,
        Some(&mut head_vertices),
    )?;
    let mut used = BTreeSet::new();
    for (tri, &f) in cavity
        .mesh()
        .indices
        .as_chunks::<3>()
        .0
        .iter()
        .zip(&cavity.mesh().face_indices)
    {
        if f != cap {
            used.extend(tri.iter().copied());
        }
    }
    let mut ids = vec![None; cavity.mesh().vertices.len() / 3];
    let mut cavity_vertices = Vec::new();
    for old in used {
        cancel(cancelled)?;
        let [x, y, d] = point(cavity.mesh(), old);
        let vertex = shell.push([x, y, screw.head_top_z() + d]);
        cavity_vertices.push(PhillipsCavityVertex {
            vertex,
            original: [x, y, d],
        });
        ids[old as usize] = Some(vertex);
    }
    for (tri, &f) in cavity
        .mesh()
        .indices
        .as_chunks::<3>()
        .0
        .iter()
        .zip(&cavity.mesh().face_indices)
    {
        if f == cap {
            continue;
        }
        cancel(cancelled)?;
        let key = topology.faces()[f as usize].source_identity();
        work.step(256)?;
        let original = tri.map(|id| point(cavity.mesh(), id));
        let tri = [tri[0], tri[2], tri[1]].map(|id| ids[id as usize].expect("used source vertex"));
        let translated = tri.map(|id| point(&shell.builder.mesh, id));
        // Reverse both copies when comparing: the material uses cavity sides
        // with the opposite original winding. Translation gets a fresh literal
        // facet-normal proof, including large-coordinate cancellation.
        translated_normal(
            [original[0], original[2], original[1]],
            translated,
            options.angle / 200.,
        )?;
        shell.builder.triangle(tri, face(ScrewFace::Cavity(key)))?;
    }
    let inner: Vec<_> = ring
        .into_iter()
        .map(|id| ids[id as usize].expect("mouth vertex belongs to cavity side"))
        .collect();
    let before_planar = maximum_work - work.remaining.get();
    top_annulus(
        &mut shell,
        &inner,
        &outer,
        face(ScrewFace::HeadTop),
        cancelled,
    )?;
    cancel(cancelled)?;
    let mut thread_vertices = Vec::with_capacity(builder.vertex_ids.len());
    work.step(
        builder.vertex_ids.len().saturating_mul(32)
            + head_vertices.len().saturating_mul(32)
            + cavity_vertices.len().saturating_mul(32),
    )?;
    for (&key, &vertex) in &builder.vertex_ids {
        let kind = match key {
            super::super::VertexKey::Grid(column, rail) => PhillipsThreadVertexKind::Grid {
                fraction: builder.fraction(column),
                phase: builder.rail_phase(rail),
                period: rail.div_euclid(builder.strips.len() as i64),
            },
            super::super::VertexKey::CapColumn(column, cap) => {
                PhillipsThreadVertexKind::CapColumn {
                    fraction: builder.fraction(column),
                    cap,
                    rail: builder.cap_rails.get(&(column, cap)).map(|&rail| {
                        (
                            rail.div_euclid(builder.strips.len() as i64),
                            builder.rail_phase(rail),
                        )
                    }),
                }
            }
            super::super::VertexKey::CapRail(rail, cap) => PhillipsThreadVertexKind::CapRail {
                phase: builder.rail_phase(rail),
                period: rail.div_euclid(builder.strips.len() as i64),
                cap,
            },
        };
        thread_vertices.push(PhillipsThreadVertex { vertex, kind });
    }
    thread_vertices.sort_unstable_by_key(|p| p.vertex);
    head_vertices.sort_unstable_by_key(|p| (p.vertex, p.segment));
    cavity_vertices.sort_unstable_by_key(|p| p.vertex);
    let mesh = builder.mesh;
    let before_proof = maximum_work - work.remaining.get();
    work.step(
        mesh.indices
            .len()
            .saturating_mul(64)
            .saturating_add(mesh.vertices.len()),
    )?;
    let connectivity = manufacturing::certify(&mesh)?;
    let literal = literal_volume(&mesh);
    let source = source_volume.volume_mm3();
    let difference = literal - I::new(source[0], source[1]);
    Ok(PhillipsScrewMesh {
        mesh,
        connectivity,
        faces,
        source_volume,
        literal_volume_mm3: [literal.lo, literal.hi],
        volume_difference_bound_mm3: difference.magnitude(),
        work: maximum_work - work.remaining.get(),
        handoffs: [after_cavity, before_planar, before_proof],
        thread_vertices,
        head_vertices,
        cavity_vertices,
        cavity_distance: cavity.source_distance_mm(),
        proposal: screw.clone(),
    })
}

#[cfg(test)]
mod tests;
