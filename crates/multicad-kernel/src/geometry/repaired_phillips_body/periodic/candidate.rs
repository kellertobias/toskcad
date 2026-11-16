//! Independently closed indexed approximation, with explicit transfer to the
//! authoritative owned periodic source. Global literal embedding follows in
//! the shared Native cold pipeline; this candidate does not skip that gate.
use super::*;
mod normals;
use multicad_rust_kernel::{
    display::{DisplayMesh, Tessellation},
    manufacturing::MeshProof,
    phillips_drive::ScrewFace,
    step_read::{
        PhillipsJoinedDegree11MaterialCandidate, PhillipsJoinedMaterialCandidate,
        PhillipsJoinedMeshOptions, PhillipsMaterialFaceSource,
    },
};
pub(in crate::geometry::repaired_phillips_body) fn literal_volume(
    mesh: &Mesh,
) -> Result<[f64; 2], GeometryError> {
    let origin = &mesh.vertices[..3];
    let mut sum = I::ZERO;
    for tri in mesh.indices.as_chunks::<3>().0 {
        let [a, b, c] = tri.map(|id| {
            std::array::from_fn::<_, 3, _>(|i| {
                I::point(mesh.vertices[3 * id as usize + i]) - I::point(origin[i])
            })
        });
        let cross = std::array::from_fn::<_, 3, _>(|i| {
            b[(i + 1) % 3] * c[(i + 2) % 3] - b[(i + 2) % 3] * c[(i + 1) % 3]
        });
        sum = sum + (0..3).fold(I::ZERO, |s, i| s + a[i] * cross[i]) / I::point(6.);
    }
    if !sum.valid() {
        return Err(unsupported(
            "Periodic literal signed-volume enclosure is unresolved",
        ));
    }
    Ok([sum.lo, sum.hi])
}
pub(in crate::geometry::repaired_phillips_body) enum Candidate<'s> {
    Retained(PhillipsJoinedMaterialCandidate<'s>),
    Periodic(PeriodicCandidate<'s>),
    Degree11(PhillipsJoinedDegree11MaterialCandidate<'s>),
}
pub(in crate::geometry::repaired_phillips_body) struct PeriodicCandidate<'s> {
    _source: &'s OwnedPeriodicMaterial,
    mesh: DisplayMesh,
    proof: MeshProof,
    used: usize,
    peak: usize,
    distance: f64,
}
impl Candidate<'_> {
    pub(in crate::geometry::repaired_phillips_body) fn periodic_surface_error(
        &self,
    ) -> Option<f64> {
        match self {
            Self::Retained(_) => None,
            Self::Periodic(c) => Some(c.distance),
            Self::Degree11(c) => Some(c.new_surface_error_mm()),
        }
    }
    pub(in crate::geometry::repaired_phillips_body) fn mesh(&self) -> &DisplayMesh {
        match self {
            Self::Retained(c) => c.mesh(),
            Self::Periodic(c) => &c.mesh,
            Self::Degree11(c) => c.mesh(),
        }
    }
    pub(in crate::geometry::repaired_phillips_body) fn connectivity(&self) -> MeshProof {
        match self {
            Self::Retained(c) => c.connectivity(),
            Self::Periodic(c) => c.proof,
            Self::Degree11(c) => c.connectivity(),
        }
    }
    pub(in crate::geometry::repaired_phillips_body) fn work_used(&self) -> usize {
        match self {
            Self::Retained(c) => c.work_used(),
            Self::Periodic(c) => c.used,
            Self::Degree11(c) => c.work_used(),
        }
    }
    pub(in crate::geometry::repaired_phillips_body) fn combined_live_upper_bound_bytes(
        &self,
    ) -> usize {
        match self {
            Self::Retained(c) => c.combined_live_upper_bound_bytes(),
            Self::Periodic(c) => c.peak,
            Self::Degree11(c) => c.combined_live_upper_bound_bytes(),
        }
    }
}
fn connected(mesh: &DisplayMesh) -> Result<(), GeometryError> {
    let mut parents: Vec<usize> = (0..mesh.vertices.len() / 3).collect();
    fn root(p: &mut [usize], mut i: usize) -> usize {
        while p[i] != i {
            p[i] = p[p[i]];
            i = p[i];
        }
        i
    }
    for t in mesh.indices.as_chunks::<3>().0 {
        let a = root(&mut parents, t[0] as usize);
        for &b in &t[1..] {
            let b = root(&mut parents, b as usize);
            parents[b] = a;
        }
    }
    let anchor = root(&mut parents, mesh.indices[0] as usize);
    if mesh
        .indices
        .iter()
        .any(|&i| root(&mut parents, i as usize) != anchor)
    {
        return Err(unsupported(
            "Periodic approximation has more than one boundary component",
        ));
    }
    Ok(())
}
pub(in crate::geometry::repaired_phillips_body) fn build<'s>(
    source: &'s OwnedPeriodicMaterial,
    o: PhillipsJoinedMeshOptions,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<PeriodicCandidate<'s>, GeometryError> {
    checkpoint(interrupt)?;
    let approximation = source
        .mesh_approximation(o.maximum_work.min(16_000_000), interrupt)
        .map_err(source_failure)?;
    let mut remaining = o
        .maximum_work
        .checked_sub(approximation.work_used())
        .ok_or_else(exhausted)?;
    let normal = approximation.maximum_normal_angle_radians();
    // The delegated producer receives only the unused local source error.
    let q = Tessellation {
        deflection: (I::point(o.tessellation.deflection) / I::point(2.)).lo,
        angle: (I::point(o.tessellation.angle) - I::point(normal)).lo,
        ..o.tessellation
    };
    if q.deflection <= 0. || q.angle < 0.01 {
        return Err(unsupported(
            "Periodic approximation cannot fit source quality reserve",
        ));
    }
    let cancel = || interrupt.is_some_and(|f| f());
    let literal = multicad_rust_kernel::helical_thread_mesh::tessellate_phillips_screw_with_work(
        approximation.proposal(),
        q,
        remaining,
        &cancel,
    )
    .map_err(|e| match e {
        multicad_rust_kernel::helical_thread_mesh::PhillipsMeshError::Cancelled => {
            failure(GeometryErrorKind::Cancelled, "Periodic mesh cancelled")
        }
        multicad_rust_kernel::helical_thread_mesh::PhillipsMeshError::Geometry(
            multicad_rust_kernel::display::MeshError::BudgetExceeded,
        ) => exhausted(),
        _ => unsupported("Periodic approximation producer refused source sewing"),
    })?;
    remaining = remaining
        .checked_sub(literal.work_used())
        .ok_or_else(exhausted)?;
    let measured = approximation
        .certify_literal(
            &literal,
            o.tessellation.deflection,
            remaining,
            o.maximum_bytes,
            interrupt,
        )
        .map_err(source_failure)?;
    remaining = remaining
        .checked_sub(measured.work_used())
        .ok_or_else(exhausted)?;
    let position = measured.maximum_outer_vertex_error_mm();
    let measured_peak = measured.accounted_peak_bytes();
    let surface_distance = measured.maximum_surface_error_mm();
    let original = literal.mesh();
    let triangles = original.indices.len() / 3;
    let vertices = original.vertices.len() / 3;
    let copies = original.vertices.len() + original.indices.len() + original.face_indices.len();
    remaining = remaining
        .checked_sub(copies + triangles * 16 + vertices * 4)
        .ok_or_else(exhausted)?;
    let actual_bytes = original.vertices.capacity() * 8
        + (original.indices.capacity() + original.face_indices.capacity()) * 4;
    let owned = source
        .retained_budget_bytes()
        .saturating_add(literal.retained_budget_bytes())
        .saturating_add(actual_bytes)
        .saturating_add(vertices * 8 + 65536);
    if owned > o.maximum_bytes {
        return Err(failure(
            GeometryErrorKind::Limit,
            "Periodic mesh consumer buffer bound exceeded",
        ));
    }
    let faces = approximation.proposal().faces();
    let mut remap = Vec::with_capacity(faces.len());
    for face in faces {
        let target = match face {
            ScrewFace::ThreadShoulder => PhillipsMaterialFaceSource::JoinAnnulus,
            ScrewFace::ThreadTip => PhillipsMaterialFaceSource::ThreadTip,
            ScrewFace::ThreadWall(i) => PhillipsMaterialFaceSource::ThreadSide {
                face: [2, 3, 0, 1][i as usize],
            },
            ScrewFace::HeadTop => PhillipsMaterialFaceSource::TopAnnulus,
            ScrewFace::ProfileWall(i) => PhillipsMaterialFaceSource::HeadMeridian {
                segment: i as usize,
            },
            ScrewFace::Cavity(face) => PhillipsMaterialFaceSource::ReversedCavity { face },
        };
        let id = (0..source.face_count())
            .find(|&i| source.face_source(i) == Some(target))
            .ok_or_else(|| {
                unsupported("Periodic approximation lost an actual source face owner")
            })?;
        remap.push(id as u32);
    }
    let mut mesh = original.clone();
    for face in &mut mesh.face_indices {
        *face = remap[*face as usize];
    }
    connected(&mesh)?;
    let used = normals::certify(
        &mesh,
        source,
        approximation.periodic_outward_normal_coefficients(),
        position,
        o.tessellation.angle,
        remaining,
        interrupt,
    )?;
    remaining -= used;
    checkpoint(interrupt)?;
    // Existing producer delegated workspaces remain separate from caller's
    // consumer bytes: Cutter refinement/trace, BTree/HashMap growth, planar
    // DCEL and connectivity. Final facet/vertex counts dominate the conforming
    // intermediate grids;32MiB also covers the counted trace metadata cap.
    let delegated = triangles
        .checked_mul(16384)
        .and_then(|n| n.checked_add(36 * 1024 * 1024))
        .ok_or_else(exhausted)?;
    let peak = owned
        .max(measured_peak)
        .checked_add(delegated)
        .ok_or_else(exhausted)?;
    Ok(PeriodicCandidate {
        _source: source,
        mesh,
        proof: literal.connectivity(),
        used: o.maximum_work - remaining,
        peak,
        distance: surface_distance,
    })
}
