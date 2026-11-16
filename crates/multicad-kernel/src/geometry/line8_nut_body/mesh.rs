//! Fresh actual-world source/facet proof, then independent literal embedding.
use super::*;
use multicad_rust_kernel::{display::Tessellation, manufacturing::MeshProof};
fn facet(v: &[f64], ids: [u32; 3]) -> Result<[I; 3], GeometryError> {
    let points = ids.map(|id| std::array::from_fn(|i| v[id as usize * 3 + i]));
    multicad_rust_kernel::literal_triangle::normal_bounds(points)
        .map(|n| n.map(iv))
        .map_err(|_| unresolved())
}
impl Line8NutBody {
    fn proved(&self, q: Tessellation) -> Result<(Mesh, MeshProof), GeometryError> {
        // The source already reserves steering headroom; actual world C0
        // and normal checks below independently decide the requested quality.
        let inner = q;
        let c = self
            .source
            .mesh_candidate(inner, 64_000_000, None)
            .map_err(fail)?;
        let proof = c.connectivity();
        let mut remaining = 128_000_000usize
            .checked_sub(c.work_used())
            .ok_or_else(unresolved)?;
        let bytes = c.retained_bytes();
        let local_chord = c.source_distance_mm();
        let local_normal = c.normal_chord();
        let mut literal = c.into_mesh();
        let cost = literal.vertices.len() / 3 * 128
            + literal.indices.len() / 3
                * (2 * multicad_rust_kernel::literal_triangle::NORMAL_WORK + 512);
        remaining = remaining.checked_sub(cost).ok_or_else(unresolved)?;
        let extra = literal.vertices.len() * 8 + 16384;
        if bytes + extra > 64 * 1024 * 1024 {
            return Err(fail(multicad_rust_kernel::line8_nut::Error::Limit));
        }
        let mut world = Vec::new();
        world
            .try_reserve_exact(literal.vertices.len())
            .map_err(|_| fail(multicad_rust_kernel::line8_nut::Error::Limit))?;
        if bytes + world.capacity() * 8 + 16384 > 64 * 1024 * 1024 {
            return Err(fail(multicad_rust_kernel::line8_nut::Error::Limit));
        }
        let mut rounding = 0_f64;
        for &p in literal.vertices.as_chunks::<3>().0 {
            let exact = self.point_bounds(p.map(|v| [v, v]))?;
            let point = exact.map(|v| v[0].midpoint(v[1]));
            rounding =
                rounding.max(norm(std::array::from_fn(|i| iv(exact[i]) - I::point(point[i]))).hi);
            world.extend(point);
        }
        let outer = self.outer();
        let scale = matrix_norm(self.stock) * matrix_norm(outer);
        if (scale * I::point(local_chord) + I::point(rounding)).hi > q.deflection {
            return Err(unresolved());
        }
        let condition = condition(self.stock) * condition(outer);
        let ca = cofactor(self.stock.matrix().map(|r| r.map(I::point)));
        let cb = cofactor(outer.matrix().map(|r| r.map(I::point)));
        let angle = I::point(q.angle) / I::point(2.);
        let allowed = (I::point(2.) * (angle - angle * angle * angle / I::point(6.))).lo;
        for &ids in literal.indices.as_chunks::<3>().0 {
            let normal = facet(&literal.vertices, ids)?;
            let a = unit(ca.map(|r| {
                r.into_iter()
                    .zip(normal)
                    .fold(I::ZERO, |s, (a, b)| s + a * b)
            }))?;
            let exact =
                unit(cb.map(|r| r.into_iter().zip(a).fold(I::ZERO, |s, (a, b)| s + a * b)))?;
            let actual = facet(&world, ids)?;
            if (condition * I::point(local_normal)
                + norm(std::array::from_fn(|i| exact[i] - actual[i])))
            .hi > allowed
            {
                return Err(unresolved());
            }
        }
        literal.vertices = world;
        let mesh = Mesh {
            vertices: literal.vertices,
            indices: literal.indices,
            face_indices: literal.face_indices,
        };
        // The finite64MiB source/placement scope ends before the unchanged embedding
        // workspace. Keep its independent work inside the same128M operation limit.
        mesh.prove_literal_embedding(remaining.min(64_000_000), None)?;
        Ok((mesh, proof))
    }
    pub(in crate::geometry) fn mesh(&self, q: Tessellation) -> Result<Mesh, GeometryError> {
        self.proved(q).map(|v| v.0)
    }
    pub(in crate::geometry) fn manufacturing(
        &self,
        q: Tessellation,
    ) -> Result<crate::ManufacturingMesh, GeometryError> {
        let (mesh, p) = self.proved(q)?;
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
