//! Actual indexed cavity sewing. Global literal embedding remains separate.
use super::*;
use crate::{
    manufacturing::{self, MeshProof},
    phillips_drive::{Error, FaceKey, PhillipsScrew, ScrewFace},
    tapered_thread::interval::I,
};
use std::collections::{BTreeMap, BTreeSet};
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
}
impl PhillipsScrewMesh {
    pub fn mesh(&self) -> &DisplayMesh {
        &self.mesh
    }
    pub fn connectivity(&self) -> MeshProof {
        self.connectivity
    }
    pub fn faces(&self) -> &[ScrewFace] {
        &self.faces
    }
}
fn point(mesh: &DisplayMesh, id: u32) -> [f64; 3] {
    mesh.vertices.as_chunks::<3>().0[id as usize]
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
    check_options(options)?;
    cancel(cancelled)?;
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
        ..options
    };
    let radius = screw.profile().support([1., 0., 0.]);
    let step = (working.angle * 0.5).min((working.deflection / radius).sqrt());
    let count = (TAU / (12. * step)).ceil() * 12.;
    if count * 2. > working.triangle_limit as f64 {
        return Err(MeshError::BudgetExceeded.into());
    }
    let cavity = source
        .surface_patches_with_cancel(cancelled)
        .map_err(source_error)?
        .boundary_mesh_with_cancel(working, cancelled)
        .map_err(source_error)?;
    let topology = source.boundary_topology().map_err(source_error)?;
    let cap = topology
        .faces()
        .iter()
        .position(|f| f.source_identity() == FaceKey::Mouth)
        .expect("qualified source mouth") as u32;
    let ring = mouth_ring(cavity.mesh(), cap, witness)?;
    let strips = phase_strips(screw.thread(), TAU / count, working)?;
    let mut builder = Builder::new(screw.thread(), count as usize, strips, working);
    builder.lower = screw.visible_thread_interval_mm()[0];
    builder.build_wall()?;
    cancel(cancelled)?;
    let [mut shoulder, tip] = std::mem::take(&mut builder.caps);
    builder.add_cap(1, tip)?;
    Builder::sort_ring(&mut shoulder)?;
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
    let outer = shell.profile_outer(
        screw.profile(),
        screw.spec().shank_length_mm,
        screw.visible_thread_interval_mm()[0],
        &shaft,
        (face(ScrewFace::ThreadShoulder), |i| {
            face(ScrewFace::ProfileWall(i))
        }),
        working,
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
    for old in used {
        cancel(cancelled)?;
        let [x, y, d] = point(cavity.mesh(), old);
        ids[old as usize] = Some(shell.push([x, y, screw.head_top_z() + d]));
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
        let tri = [tri[0], tri[2], tri[1]].map(|id| ids[id as usize].expect("used source vertex"));
        shell.builder.triangle(tri, face(ScrewFace::Cavity(key)))?;
    }
    let inner: Vec<_> = ring
        .into_iter()
        .map(|id| ids[id as usize].expect("mouth vertex belongs to cavity side"))
        .collect();
    top_annulus(
        &mut shell,
        &inner,
        &outer,
        face(ScrewFace::HeadTop),
        cancelled,
    )?;
    cancel(cancelled)?;
    let mesh = builder.mesh;
    let connectivity = manufacturing::certify(&mesh)?;
    Ok(PhillipsScrewMesh {
        mesh,
        connectivity,
        faces,
    })
}

#[cfg(test)]
mod tests;
