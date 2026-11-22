//! One source/hex-join/actual-world/literal-embedding work ledger. Display
//! reuse belongs to this immutable owner; manufacturing always proves afresh.
use super::*;
use multicad_rust_kernel::{
    axial_pipe_thread::authored::hardware::mesh::{self as joined, ThreadBoundaryMesh},
    display::DisplayMesh,
    manufacturing::MeshProof,
};
const WORK: usize = 512_000_000;
const BYTES: usize = 64 * 1024 * 1024;
pub(super) struct Cached {
    q: Tessellation,
    mesh: Mesh,
}
fn same(a: Tessellation, b: Tessellation) -> bool {
    a.deflection.to_bits() == b.deflection.to_bits()
        && a.angle.to_bits() == b.angle.to_bits()
        && a.triangle_limit == b.triangle_limit
}
fn bytes(mesh: &Mesh) -> usize {
    mesh.vertices.capacity() * 8 + (mesh.indices.capacity() + mesh.face_indices.capacity()) * 4
}
fn copy(mesh: &Mesh) -> Mesh {
    Mesh {
        vertices: mesh.vertices.clone(),
        indices: mesh.indices.clone(),
        face_indices: mesh.face_indices.clone(),
    }
}
struct Work<'a> {
    remaining: usize,
    held: usize,
    peak: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Work<'_> {
    fn budget(&self) -> Result<usize, GeometryError> {
        if self.remaining == 0 {
            Err(failure(Cause::Limit))
        } else {
            Ok(self.remaining)
        }
    }
    fn charge(&mut self, n: usize) -> Result<(), GeometryError> {
        checkpoint(self.interrupt)?;
        self.remaining = self
            .remaining
            .checked_sub(n)
            .ok_or_else(|| failure(Cause::Limit))?;
        Ok(())
    }
    fn storage(&mut self, n: usize) -> Result<(), GeometryError> {
        let used = self
            .held
            .checked_add(n)
            .filter(|v| *v <= BYTES)
            .ok_or_else(|| failure(Cause::Limit))?;
        self.peak = self.peak.max(used);
        Ok(())
    }
}
impl PipeHardwareBody {
    fn resident_bytes(&self) -> usize {
        // The fresh authored facade has only its identity placement/cache and
        // <=34faces/112edges/80vertices.1MiB conservatively covers both Native
        // headers, scalar descriptor caches, topology/role strings and maps,
        // plus added12vertices/18edges/8faces. Original pole/normal owners are
        // charged by the actual immutable RadialBoundary storage receipt.
        let boundary = self
            .source
            .thread
            .rust_radial_boundary()
            .map_or(0, |(b, _)| b.retained_budget_bytes());
        1024 * 1024 + boundary + self.placement.retained_bytes()
    }
    fn local_candidate(
        &self,
        q: Tessellation,
        work: &mut Work<'_>,
    ) -> Result<(DisplayMesh, MeshProof, f64), GeometryError> {
        // The join checks its actual cap homotopy/rounding separately from
        // unchanged lateral source quality. The final maximum and world C0
        // checks retain q; do not overspecify source refinement by a fixed1%.
        let source_q = q;
        quality(source_q)?;
        let finish = |trace: ThreadBoundaryMesh<'_, '_>,
                      work: &mut Work<'_>|
         -> Result<(DisplayMesh, MeshProof, f64), GeometryError> {
            let trace_work = match trace {
                ThreadBoundaryMesh::Clipped(t) => t.work_used(),
                ThreadBoundaryMesh::Retained(t) => t.work_used(),
            };
            work.charge(trace_work)?;
            #[cfg(test)]
            eprintln!(
                "PIPE_MESH before_stock_join source_work={trace_work} remaining={}",
                work.remaining
            );
            let candidate = joined::join(
                trace,
                self.source.intent.kind(),
                joined::MeshOptions {
                    tessellation: q,
                    maximum_work: work.budget()?,
                    maximum_bytes: BYTES,
                    additional_resident_bytes: work.held,
                },
                work.interrupt,
            )
            .map_err(|e| {
                error(
                    "pipe hardware stock join",
                    failure(e).kind,
                    format!(
                        "Actual stock join refused with {} work remaining: {e:?}",
                        work.remaining
                    ),
                )
            })?;
            #[cfg(test)]
            eprintln!(
                "PIPE_MESH source_work={} join_work={} remaining_before_join_debit={}",
                trace_work,
                candidate.joining_work_used(),
                work.remaining
            );
            work.charge(candidate.joining_work_used())?;
            work.peak = work.peak.max(candidate.accounted_peak_bytes());
            let proof = candidate.connectivity();
            let distance = candidate.source_distance_mm();
            Ok((candidate.into_mesh(), proof, distance))
        };
        if let Some(source) = self.source.thread.rust_clipped_npt_source() {
            let spec = source.source().spec();
            let mut proposed = None;
            for n in [128usize, 256, 512] {
                work.charge(256)?;
                let omega = std::f64::consts::TAU / spec.pitch_mm;
                let ds = spec.pitch_mm / n as f64;
                let dz = spec.length_mm / 128.;
                let bound = 0.5
                    * ((omega * omega * source.stock_radius_limit_mm() + 4. * omega) * ds * ds
                        + 4. * omega * ds * dz);
                if bound < source_q.deflection * 0.9 && 24. / (n as f64) < source_q.angle {
                    proposed = Some(n);
                    break;
                }
            }
            let n = proposed.ok_or_else(|| failure(Cause::Limit))?;
            let available = BYTES
                .checked_sub(work.held)
                .filter(|n| *n > 0)
                .ok_or_else(|| failure(Cause::Limit))?;
            let plan = source
                .cell_plan(n, 128, work.budget()?, available, work.interrupt)
                .map_err(failure)?;
            let trace = plan
                .mesh_candidate_with_cap_trace(source_q, work.budget()?, available, work.interrupt)
                .map_err(failure)?;
            finish(ThreadBoundaryMesh::Clipped(&trace), work)
        } else {
            let (source, frame) = self
                .source
                .thread
                .rust_radial_boundary()
                .ok_or_else(|| failure(Cause::Unsupported))?;
            if frame != PlaneFrame::default() {
                return Err(failure(Cause::Unsupported));
            }
            let trace = source
                .tessellate_with_cap_trace_with_work_budget(
                    source_q,
                    Default::default(),
                    work.budget()?,
                    work.interrupt,
                )
                .map_err(|e| {
                    error(
                        "pipe hardware retained thread mesh",
                        radial_failure(e).kind,
                        format!("Actual retained thread mesh refused: {e:?}"),
                    )
                })?;
            finish(ThreadBoundaryMesh::Retained(&trace), work)
        }
    }
    fn fresh(
        &self,
        q: Tessellation,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<(Mesh, MeshProof), GeometryError> {
        quality(q)?;
        checkpoint(interrupt)?;
        self.mesh_builds
            .set(self.mesh_builds.get().saturating_add(1));
        let cached = self.display.borrow().as_ref().map_or(0, |c| bytes(&c.mesh));
        let mut work = Work {
            remaining: WORK,
            held: self
                .resident_bytes()
                .checked_add(cached)
                .ok_or_else(|| failure(Cause::Limit))?,
            peak: 0,
            interrupt,
        };
        work.storage(0)?;
        let local = self.placement.local_quality(q)?;
        #[cfg(test)]
        eprintln!(
            "PIPE_MESH requested {:?} local {:?} residents {}",
            q, local, work.held
        );
        let (mut original, mut proof, distance) = self.local_candidate(local, &mut work)?;
        for face in &mut original.face_indices {
            work.charge(1)?;
            *face = *self
                .source
                .face_map
                .get(face)
                .expect("same finite source/hex face owner");
        }
        let mesh = if self.placement.is_identity() {
            Mesh {
                vertices: original.vertices,
                indices: original.indices,
                face_indices: original.face_indices,
            }
        } else {
            let original_bytes = original.vertices.capacity() * 8
                + (original.indices.capacity() + original.face_indices.capacity()) * 4;
            work.storage(original_bytes * 2)?;
            let mut target = DisplayMesh {
                vertices: Vec::with_capacity(original.vertices.len()),
                indices: original.indices.clone(),
                face_indices: original.face_indices.clone(),
            };
            let mut rounding = 0f64;
            for &p in original.vertices.as_chunks::<3>().0 {
                work.charge(self.placement.point_work())?;
                let (point, error) = self.placement.point(p)?;
                rounding = rounding.max(error);
                target.vertices.extend(point);
            }
            let bound = I::point(distance) * I::point(self.placement.sigma[1]) + I::point(rounding);
            if !bound.valid() || bound.hi > q.deflection {
                return Err(failure(Cause::Unresolved));
            }
            let reserve = (I::point(q.angle) / I::point(4096.)).lo;
            let anchor = std::array::from_fn::<_, 3, _>(|k| I::point(target.vertices[k]));
            let mut flux = I::ZERO;
            for ids in original.indices.as_chunks::<3>().0 {
                work.charge(self.placement.normal_work())?;
                let a =
                    ids.map(|id| std::array::from_fn(|k| original.vertices[3 * id as usize + k]));
                let b = ids.map(|id| std::array::from_fn(|k| target.vertices[3 * id as usize + k]));
                let actual_normal = self
                    .placement
                    .require_normal(a, b, reserve, |n| work.charge(n))?;
                // The already earned normal encloses cross(p1-p0,p2-p0)
                // of these SAME literal world points. Exact determinant
                // multilinearity gives dot(p0-anchor,normal)/6 for their
                // anchored signed tetrahedral volume; no local mass reuse.
                work.charge(24)?;
                let at = std::array::from_fn::<_, 3, _>(|k| I::point(b[0][k]) - anchor[k]);
                flux =
                    flux + (0..3).fold(I::ZERO, |v, k| v + at[k] * actual_normal[k]) / I::point(6.);
            }
            work.storage(
                original_bytes
                    + target.vertices.capacity() * 8
                    + (target.indices.capacity() + target.face_indices.capacity()) * 4,
            )?;
            // These are the SAME indices, vertex count and face-array length.
            // Finite mapped points and the strictly positive actual facet
            // projection above earn nonzero consistently oriented world faces.
            // Therefore indexed edge pairing and cyclic vertex links are
            // unchanged: do not sort/union the identical graph a second time.
            // The old local flux is NOT transported as a world scalar. Earn
            // positive actual signed flux around a shared actual world anchor;
            // the unchanged independent literal embedding proof follows below.
            if !flux.valid() || flux.lo <= 0. {
                return Err(failure(Cause::Unresolved));
            }
            proof = MeshProof {
                edge_pairs: proof.edge_pairs,
                vertex_links: proof.vertex_links,
                signed_volume_mm3: flux.lo.midpoint(flux.hi),
            };
            Mesh {
                vertices: target.vertices,
                indices: target.indices,
                face_indices: target.face_indices,
            }
        };
        work.storage(bytes(&mesh))?;
        // Original retained preparation and the generic intersection validator
        // have their existing separately bounded delegated workspaces. The64MiB
        // receipt above covers the Native owner/cache/transfer buffers and the
        // actual source/join buffer phases, not whole-process allocator RSS.
        work.charge(1)?;
        #[cfg(test)]
        eprintln!(
            "PIPE_MESH before_embedding remaining={} facets={}",
            work.remaining,
            mesh.indices.len() / 3
        );
        let embedding = mesh.prove_literal_embedding(work.remaining, interrupt)?;
        let used = embedding.work_used();
        #[cfg(test)]
        eprintln!(
            "PIPE_MESH embedding_work={used} pairs={}",
            embedding.pairs()
        );
        work.charge(used)?;
        checkpoint(interrupt)?;
        Ok((mesh, proof))
    }
    pub(in crate::geometry) fn mesh(&self, q: Tessellation) -> Result<Mesh, GeometryError> {
        quality(q)?;
        if let Some(cached) = self.display.borrow().as_ref().filter(|c| same(c.q, q)) {
            if self
                .resident_bytes()
                .saturating_add(2 * bytes(&cached.mesh))
                > BYTES
            {
                return Err(failure(Cause::Limit));
            }
            return Ok(copy(&cached.mesh));
        }
        let (mesh, _) = self.fresh(q, None)?;
        // Atomically publish only a complete source/world/embedding result.
        // A failed request or fresh manufacturing never clears the old cache.
        if self.resident_bytes().saturating_add(2 * bytes(&mesh)) <= BYTES {
            *self.display.borrow_mut() = Some(Cached {
                q,
                mesh: copy(&mesh),
            });
        }
        Ok(mesh)
    }
    pub(in crate::geometry) fn manufacturing(
        &self,
        q: Tessellation,
    ) -> Result<crate::ManufacturingMesh, GeometryError> {
        let (mesh, proof) = self.fresh(q, None)?;
        Ok(crate::ManufacturingMesh {
            mesh,
            maximum_displacement_mm: 0.,
            edge_pairs: proof.edge_pairs,
            vertex_links: proof.vertex_links,
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
