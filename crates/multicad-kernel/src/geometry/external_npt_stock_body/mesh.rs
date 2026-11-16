use super::*;
use multicad_rust_kernel::display::Tessellation;
impl ExternalNptStockBody {
    fn proved(
        &self,
        q: Tessellation,
    ) -> Result<(Mesh, multicad_rust_kernel::manufacturing::MeshProof), GeometryError> {
        let c = self
            .source
            .mesh_candidate(q, 512_000_000, 64 * 1024 * 1024, None)
            .map_err(failure)?;
        let proof = c.connectivity();
        let mut remaining = 512_000_000 - c.work_used();
        let chord = c.source_distance_mm();
        let normal = c.normal_chord();
        let retained = c.retained_budget_bytes();
        let mut literal = c.into_mesh();
        if self.pose != RigidTransform::default() {
            let calibration = self.calibration();
            let extra = literal
                .vertices
                .len()
                .checked_mul(8)
                .and_then(|v| v.checked_add(16384))
                .ok_or_else(|| failure(Cause::Limit))?;
            if retained
                .checked_add(extra)
                .is_none_or(|v| v > 64 * 1024 * 1024)
            {
                return Err(failure(Cause::Limit));
            }
            let cost = literal.vertices.len() / 3 * 128 + literal.indices.len() / 3 * 640;
            remaining = remaining
                .checked_sub(cost)
                .ok_or_else(|| failure(Cause::Limit))?;
            let mut world = Vec::new();
            world
                .try_reserve_exact(literal.vertices.len())
                .map_err(|_| failure(Cause::Limit))?;
            if retained + world.capacity() * 8 + 16384 > 64 * 1024 * 1024 {
                return Err(failure(Cause::Limit));
            }
            let mut rounding = 0_f64;
            for &p in literal.vertices.as_chunks::<3>().0 {
                let exact = calibration
                    .point_bounds(p.map(|v| [v, v]))
                    .map_err(failure)?;
                let v = exact.map(|x| x[0].midpoint(x[1]));
                rounding =
                    rounding.max(norm(std::array::from_fn(|i| iv(exact[i]) - I::point(v[i]))).hi);
                world.extend(v);
            }
            let delta = calibration.metric_error();
            let scale = (I::point(1.) + I::point(delta)).sqrt();
            if (scale * I::point(chord) + I::point(rounding)).hi > q.deflection {
                return Err(failure(Cause::Unresolved));
            }
            let m = self.pose.rotation.map(|r| r.map(I::point));
            let pairs = [(1, 2), (2, 0), (0, 1)];
            let cofactor: [[I; 3]; 3] = std::array::from_fn(|i| {
                std::array::from_fn(|j| {
                    let (a, b) = pairs[i];
                    let (c, d) = pairs[j];
                    m[a][c] * m[b][d] - m[a][d] * m[b][c]
                })
            });
            let condition = (I::point(1.) + I::point(delta)) / (I::point(1.) - I::point(delta));
            // sin(x)>=x-x^3/6 for 0<=x<=0.5: an outward lower
            // bound on the allowed unit-normal chord without libm authority.
            let angle = I::point(q.angle) / I::point(2.);
            let allowed = (I::point(2.) * (angle - angle * angle * angle / I::point(6.))).lo;
            for &ids in literal.indices.as_chunks::<3>().0 {
                let n = facet(&literal.vertices, ids)?;
                let transformed = unit(
                    cofactor.map(|r| r.into_iter().zip(n).fold(I::ZERO, |s, (a, b)| s + a * b)),
                )?;
                let actual = facet(&world, ids)?;
                let err = (condition * I::point(normal)
                    + norm(std::array::from_fn(|i| actual[i] - transformed[i])))
                .hi;
                if err > allowed {
                    return Err(failure(Cause::Unresolved));
                }
            }
            // Same indices preserve the already-earned edge and cyclic link
            // identities; all actual world facets above have nonzero area.
            literal.vertices = world;
        }
        let mesh = Mesh {
            vertices: literal.vertices,
            indices: literal.indices,
            face_indices: literal.face_indices,
        };
        remaining = remaining.min(64_000_000);
        mesh.prove_literal_embedding(remaining, None)?;
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

fn cross(a: [I; 3], b: [I; 3]) -> [I; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn unit(a: [I; 3]) -> Result<[I; 3], GeometryError> {
    let n = norm(a);
    if n.lo <= 0. || !n.hi.is_finite() {
        return Err(failure(Cause::Unresolved));
    }
    Ok(a.map(|v| v / n))
}
fn facet(v: &[f64], ids: [u32; 3]) -> Result<[I; 3], GeometryError> {
    let [a, b, c]: [[I; 3]; 3] =
        ids.map(|id| std::array::from_fn(|i| I::point(v[id as usize * 3 + i])));
    unit(cross(
        std::array::from_fn(|i| b[i] - a[i]),
        std::array::from_fn(|i| c[i] - a[i]),
    ))
}
