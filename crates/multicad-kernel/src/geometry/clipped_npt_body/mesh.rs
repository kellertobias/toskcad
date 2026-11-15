use super::*;
use multicad_rust_kernel::manufacturing::MeshProof;
const WORK: usize = 512_000_000;
const BYTES: usize = 64 * 1024 * 1024;
fn grid(body: &ClippedNptBody, q: Tessellation) -> Result<usize, GeometryError> {
    quality(q).map_err(|_| failure(Cause::Invalid))?;
    let spec = body.source.source().spec();
    let r = body.source.stock_radius_limit_mm();
    // Steering only. Full source/facet/normal certification below decides
    // acceptance; no grid estimate serves as geometry authority.
    for n in [128usize, 256, 512] {
        let omega = std::f64::consts::TAU / spec.pitch_mm;
        let ds = spec.pitch_mm / n as f64;
        let dz = spec.length_mm / 128.;
        let bound = 0.5 * ((omega * omega * r + 4. * omega) * ds * ds + 4. * omega * ds * dz);
        if bound < q.deflection * 0.9 && 24. / (n as f64) < q.angle {
            return Ok(n);
        }
    }
    Err(failure(Cause::Limit))
}
impl ClippedNptBody {
    fn proved(
        &self,
        q: Tessellation,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<(Mesh, MeshProof), GeometryError> {
        checkpoint(interrupt)?;
        let n = grid(self, q)?;
        // Construction includes this fixed steering allowance. The returned
        // candidate work already includes its original cell-plan prerequisite.
        let available = WORK - 1024;
        let plan = self
            .source
            .cell_plan(n, 128, available, BYTES, interrupt)
            .map_err(failure)?;
        let c = plan
            .mesh_candidate(q, available, BYTES, interrupt)
            .map_err(failure)?;
        let proof = c.connectivity();
        let remaining = available - c.work_used();
        let literal = c.into_mesh();
        drop(plan);
        let mesh = Mesh {
            vertices: literal.vertices,
            indices: literal.indices,
            face_indices: literal.face_indices,
        };
        // Source64MiB scope is complete before this unchanged bounded literal
        // embedding phase. Do not present its triangle-limited workspace as
        // part of the source64MiB accounting or as measured process RSS.
        mesh.prove_literal_embedding(remaining.min(64_000_000), interrupt)?;
        checkpoint(interrupt)?;
        Ok((mesh, proof))
    }
    pub(in crate::geometry) fn mesh(&self, q: Tessellation) -> Result<Mesh, GeometryError> {
        self.proved(q, None).map(|x| x.0)
    }
    pub(in crate::geometry) fn manufacturing(
        &self,
        q: Tessellation,
    ) -> Result<crate::ManufacturingMesh, GeometryError> {
        let (mesh, p) = self.proved(q, None)?;
        Ok(crate::ManufacturingMesh {
            mesh,
            maximum_displacement_mm: 0.,
            edge_pairs: p.edge_pairs,
            vertex_links: p.vertex_links,
            removed_native_pole_triangles: 0,
            removed_exact_zero_triangles: 0,
            exact_vertex_pairs: 0,
            internal_constraint_edges: 0,
            removed_meshing_internal_edges: 0,
            removed_meshing_constraint_wire_containers: 0,
            removed_meshing_internal_vertices: 0,
            remeshed_without_internal_constraints: false,
        })
    }
}
