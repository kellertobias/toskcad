//! Actual source/world candidate plus independent literal separation, sharing
//! the original1.6B balance. No geometry fallback or mesh repair is performed.
use super::*;
const WORK: usize = 1_600_000_000;
struct Built {
    mesh: Mesh,
    proof: multicad_rust_kernel::manufacturing::MeshProof,
    record: MeshRecord,
}
impl ActualClippedHeadBody {
    fn build_mesh(
        &self,
        q: Tessellation,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Built, GeometryError> {
        let options = read::full_head::mesh::HeadMeshOptions {
            tessellation: q,
            maximum_work: WORK
                .checked_sub(self.preparation)
                .ok_or_else(|| failure(read::Error::Limit))?,
            maximum_bytes: BYTES,
            additional_resident_bytes: self.external_bytes(),
        };
        let candidate = if let Some(a) = &self.occurrence {
            self.source.placed_mesh_candidate(a, options, interrupt)
        } else {
            self.source.mesh_candidate(options, interrupt)
        }
        .map_err(failure)?;
        let receipt = candidate.receipt();
        let connectivity = candidate.connectivity();
        let raw = candidate.into_mesh();
        let mesh = Mesh {
            vertices: raw.vertices,
            indices: raw.indices,
            face_indices: raw.face_indices,
        };
        if interrupt.is_some_and(|f| f()) {
            return Err(failure(read::Error::Cancelled));
        }
        let native_work = mesh
            .face_indices
            .len()
            .checked_add(128)
            .ok_or_else(|| failure(read::Error::Limit))?;
        let embedding_balance = receipt
            .remaining_work
            .checked_sub(native_work)
            .ok_or_else(|| failure(read::Error::Limit))?;

        // Fixed17actual face groups. BVH boxes48n, nodes capacity2n*96B,
        // IDs8n, fan labels16n, grouping/rim/tree bookkeeping and traversal
        // are conservatively bounded by512n+2MiB. This preflight is a NEW
        // post-child phase, not a reinterpretation of the inherited peak.
        let workspace = mesh
            .face_indices
            .len()
            .checked_mul(512)
            .and_then(|n| n.checked_add(2 * 1024 * 1024))
            .ok_or_else(|| failure(read::Error::Limit))?;
        let peak = self
            .resident_bytes()
            .checked_add(mesh_bytes(&mesh))
            .and_then(|n| n.checked_add(workspace))
            .ok_or_else(|| failure(read::Error::Limit))?;
        if peak > BYTES || mesh.face_indices.iter().any(|f| *f >= 17) {
            return Err(failure(read::Error::Limit));
        }
        let embedding = mesh.prove_literal_embedding(embedding_balance, interrupt)?;
        let work = receipt
            .aggregate_work
            .checked_add(self.preparation)
            .and_then(|n| n.checked_add(native_work))
            .and_then(|n| n.checked_add(embedding.work_used()))
            .ok_or_else(|| failure(read::Error::Limit))?;
        if work > WORK {
            return Err(failure(read::Error::Limit));
        }
        let inherited = self.source.receipt().aggregate_work
            + self.occurrence.as_ref().map_or(0, |a| a.work_used())
            + self.preparation;
        let fresh = work
            .checked_sub(inherited)
            .ok_or_else(|| failure(read::Error::Limit))?;
        let record = MeshRecord {
            work,
            fresh,
            peak: peak.max(receipt.accounted_peak_bytes),
            source_error: receipt.source_position_error_mm,
            normal: receipt.normal_tangent_cone,
            pairs: embedding.pairs(),
        };
        Ok(Built {
            mesh,
            proof: connectivity,
            record,
        })
    }
    fn copy_mesh(
        &self,
        input: &Mesh,
        record: MeshRecord,
        already_cached: bool,
    ) -> Result<(Mesh, MeshRecord), GeometryError> {
        let size = mesh_bytes(input);
        let peak = self
            .resident_bytes()
            .checked_add(if already_cached { 0 } else { size })
            .and_then(|n| n.checked_add(size))
            .ok_or_else(|| failure(read::Error::Limit))?;
        if peak > BYTES {
            return Err(failure(read::Error::Limit));
        }
        let copy_work = input
            .vertices
            .len()
            .checked_add(input.indices.len())
            .and_then(|n| n.checked_add(input.face_indices.len()))
            .and_then(|n| n.checked_add(32))
            .ok_or_else(|| failure(read::Error::Limit))?;
        let work = record
            .work
            .checked_add(copy_work)
            .filter(|n| *n <= WORK)
            .ok_or_else(|| failure(read::Error::Limit))?;
        let mut mesh = Mesh {
            vertices: Vec::new(),
            indices: Vec::new(),
            face_indices: Vec::new(),
        };
        mesh.vertices
            .try_reserve_exact(input.vertices.len())
            .map_err(|_| failure(read::Error::Limit))?;
        mesh.indices
            .try_reserve_exact(input.indices.len())
            .map_err(|_| failure(read::Error::Limit))?;
        mesh.face_indices
            .try_reserve_exact(input.face_indices.len())
            .map_err(|_| failure(read::Error::Limit))?;
        if mesh_bytes(&mesh) > size {
            return Err(failure(read::Error::Limit));
        }
        mesh.vertices.extend_from_slice(&input.vertices);
        mesh.indices.extend_from_slice(&input.indices);
        mesh.face_indices.extend_from_slice(&input.face_indices);
        Ok((
            mesh,
            MeshRecord {
                work,
                fresh: if already_cached {
                    copy_work
                } else {
                    record.fresh + copy_work
                },
                peak: peak.max(record.peak),
                ..record
            },
        ))
    }
    pub(in crate::geometry) fn mesh(&self, q: Tessellation) -> Result<Mesh, GeometryError> {
        let key = [
            q.deflection.to_bits(),
            q.angle.to_bits(),
            q.triangle_limit as u64,
        ];
        if let Some(cache) = self.display.borrow().as_ref().filter(|c| c.key == key) {
            let record = cache.record;
            let (mesh, record) = self.copy_mesh(&cache.mesh, record, true)?;
            self.last_mesh.set(Some(record));
            return Ok(mesh);
        }
        let built = self.build_mesh(q, None)?;
        let (copy, record) = self.copy_mesh(&built.mesh, built.record, false)?;
        *self.display.borrow_mut() = Some(CachedDisplay {
            key,
            mesh: built.mesh,
            record: built.record,
        });
        self.last_mesh.set(Some(record));
        self.mesh_builds.set(self.mesh_builds.get() + 1);
        Ok(copy)
    }
    pub(in crate::geometry) fn manufacturing(
        &self,
        q: Tessellation,
    ) -> Result<crate::ManufacturingMesh, GeometryError> {
        let built = self.build_mesh(q, None)?;
        self.last_mesh.set(Some(built.record));
        self.mesh_builds.set(self.mesh_builds.get() + 1);
        Ok(crate::ManufacturingMesh {
            mesh: built.mesh,
            maximum_displacement_mm: 0.,
            edge_pairs: built.proof.edge_pairs,
            vertex_links: built.proof.vertex_links,
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
