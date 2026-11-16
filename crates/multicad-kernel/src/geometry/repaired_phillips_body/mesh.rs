//! Fresh manufacturing proof; one immutable-owner display receipt may be reused.
use super::*;
use multicad_rust_kernel::{
    display::{DisplayMesh, Tessellation},
    manufacturing::{self, MeshProof},
    step_read::PhillipsJoinedMeshOptions as JoinedMeshOptions,
};
pub(super) const SCRATCH: usize = 65536;
const CACHE_BYTES: usize = 32 * 1024 * 1024;
// Per triangle: 26 interval primitives, 27 scalar indexing operations,
// 19 point constructions, <=24 loads/copies and <=20 control actions.
// Round 116 up to 128, plus 16 fixed setup/finalization actions. As elsewhere
// in this ledger, an interval primitive is one action, not one CPU instruction.
const LITERAL_VOLUME_TRIANGLE_WORK: usize = 128;
const LITERAL_VOLUME_FIXED_WORK: usize = 16;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RepairedPhillipsMeshBudget {
    pub maximum_work: usize,
    /// Native owner/result/transfer/cache buffers. Delegated source/provider
    /// temporaries retain their separately reported conservative bounds.
    pub maximum_bytes: usize,
}
impl Default for RepairedPhillipsMeshBudget {
    fn default() -> Self {
        Self {
            maximum_work: WORK,
            maximum_bytes: BYTES,
        }
    }
}
#[derive(Clone, Copy, PartialEq, Eq)]
struct Key {
    owner: usize,
    effective: [u64; 2],
    quality: [u64; 2],
    triangles: usize,
    budget: RepairedPhillipsMeshBudget,
    pose: [[u64; 3]; 4],
}
impl Key {
    fn new(
        source: &Source,
        frame: PlaneFrame,
        quality: Tessellation,
        local: Tessellation,
        budget: RepairedPhillipsMeshBudget,
    ) -> Self {
        Self {
            owner: source.identity(),
            effective: [local.deflection.to_bits(), local.angle.to_bits()],
            quality: [quality.deflection.to_bits(), quality.angle.to_bits()],
            triangles: quality.triangle_limit,
            budget,
            pose: [frame.origin, frame.u, frame.v, frame.normal()].map(|p| p.map(f64::to_bits)),
        }
    }
}
pub(super) fn local_quality(
    options: Tessellation,
    placement: Placement,
) -> Result<Tessellation, GeometryError> {
    if placement.is_identity() {
        return Ok(options);
    }
    let local = Tessellation {
        deflection: (I::point(options.deflection) / (I::point(2.) * I::point(placement.sigma_max)))
            .lo,
        angle: (I::point(0.99) * I::point(options.angle) / I::point(placement.condition)).lo,
        ..options
    };
    if local.angle < 0.01 {
        return Err(unsupported(
            "World angular request leaves no qualified source reserve",
        ));
    }
    Ok(local)
}
#[derive(Clone, Copy)]
struct Metadata {
    quality: Tessellation,
    occurrence: PlaneFrame,
    proof: MeshProof,
    volume: [f64; 2],
    transfer: f64,
    embedding_pairs: usize,
    source_work: usize,
    proof_work: usize,
    proof_remaining: usize,
    query_work: usize,
    remaining: usize,
    producer_peak: usize,
    consumer_peak: usize,
    validator_peak: usize,
    combined: usize,
    cached: bool,
    periodic_literal_volume: Option<[f64; 2]>,
    periodic_mass_difference: Option<f64>,
    periodic_surface_error: Option<f64>,
    #[cfg(test)]
    timings: [std::time::Duration; 4],
}
/// Diagnostic/quality receipt bound to this immutable repaired source. The
/// source interval and literal signed volume have no asserted equality or
/// source-to-mesh mass error. This is not a raw STEP material admission token.
pub struct RepairedPhillipsMeshReceipt<'s> {
    source: &'s OwnedDerivedMaterial,
    metadata: Metadata,
}
/// Fresh or exact-key cached proof attached to its immutable typed source.
/// The default source is periodic; the degree-eleven alias remains distinct.
pub struct PeriodicRepairedPhillipsMeshReceipt<'s, S = OwnedPeriodicMaterial> {
    source: &'s S,
    metadata: Metadata,
}
pub type Degree11RepairedPhillipsMeshReceipt<'s> =
    PeriodicRepairedPhillipsMeshReceipt<'s, OwnedDegree11Material>;
impl<S> PeriodicRepairedPhillipsMeshReceipt<'_, S> {
    pub fn source(&self) -> &S {
        self.source
    }
    pub fn occurrence(&self) -> PlaneFrame {
        self.metadata.occurrence
    }
    pub fn quality(&self) -> MeshOptions {
        let q = self.metadata.quality;
        MeshOptions {
            linear_deflection_mm: q.deflection,
            angular_deflection_rad: q.angle,
            maximum_triangles: q.triangle_limit as u32,
        }
    }
    pub fn new_source_volume_mm3(&self) -> [f64; 2] {
        self.metadata.volume
    }
    pub fn literal_signed_volume_mm3(&self) -> f64 {
        self.metadata.proof.signed_volume_mm3
    }
    pub fn embedding_pairs(&self) -> usize {
        self.metadata.embedding_pairs
    }
    pub fn query_work_used(&self) -> usize {
        self.metadata.query_work
    }
    pub fn remaining_work(&self) -> usize {
        self.metadata.remaining
    }
    pub fn source_work_used(&self) -> usize {
        self.metadata.source_work
    }
    pub fn consumer_accounted_peak_bytes(&self) -> usize {
        self.metadata.consumer_peak
    }
    pub fn producer_combined_upper_bound_bytes(&self) -> usize {
        self.metadata.producer_peak
    }
    pub fn validator_workspace_upper_bound_bytes(&self) -> usize {
        self.metadata.validator_peak
    }
    pub fn combined_live_upper_bound_bytes(&self) -> usize {
        self.metadata.combined
    }
    pub fn reused_display_receipt(&self) -> bool {
        self.metadata.cached
    }
    pub fn literal_volume_enclosure_mm3(&self) -> [f64; 2] {
        self.metadata
            .periodic_literal_volume
            .expect("periodic cold path and exact-owner cache retain its directed literal integral")
    }
    /// Direct interval difference of this actual literal mesh and the distinct
    /// authoritative periodic source; no scalar equality is asserted.
    pub fn source_volume_error_bound_mm3(&self) -> f64 {
        self.metadata.periodic_mass_difference.expect(
            "periodic cold path measures the same owned source against its literal triangles",
        )
    }
    pub fn local_source_surface_error_mm(&self) -> f64 {
        self.metadata
            .periodic_surface_error
            .expect("periodic candidate retains its directed whole-source surface certificate")
    }
}
impl RepairedPhillipsMeshReceipt<'_> {
    pub fn source(&self) -> &OwnedDerivedMaterial {
        self.source
    }
    pub fn quality(&self) -> MeshOptions {
        let q = self.metadata.quality;
        MeshOptions {
            linear_deflection_mm: q.deflection,
            angular_deflection_rad: q.angle,
            maximum_triangles: q.triangle_limit as u32,
        }
    }
    pub fn occurrence(&self) -> PlaneFrame {
        self.metadata.occurrence
    }
    pub fn new_source_volume_mm3(&self) -> [f64; 2] {
        self.metadata.volume
    }
    pub fn literal_signed_volume_mm3(&self) -> f64 {
        self.metadata.proof.signed_volume_mm3
    }
    pub fn point_transfer_bound_mm(&self) -> f64 {
        self.metadata.transfer
    }
    pub fn embedding_pairs(&self) -> usize {
        self.metadata.embedding_pairs
    }
    pub fn source_work_used(&self) -> usize {
        self.metadata.source_work
    }
    pub fn proof_work_used(&self) -> usize {
        self.metadata.proof_work
    }
    pub fn proof_remaining_work(&self) -> usize {
        self.metadata.proof_remaining
    }
    #[cfg(test)]
    pub(super) fn timings(&self) -> [std::time::Duration; 4] {
        self.metadata.timings
    }
    pub fn query_work_used(&self) -> usize {
        self.metadata.query_work
    }
    pub fn remaining_work(&self) -> usize {
        self.metadata.remaining
    }
    pub fn consumer_accounted_peak_bytes(&self) -> usize {
        self.metadata.consumer_peak
    }
    pub fn producer_combined_upper_bound_bytes(&self) -> usize {
        self.metadata.producer_peak
    }
    /// Sequential Native manifold/embedding workspace allocation requests,
    /// outside the caller-lowerable consumer buffer budget; not RSS.
    pub fn validator_workspace_upper_bound_bytes(&self) -> usize {
        self.metadata.validator_peak
    }
    /// Conservative sum, including sequential delegated workspaces; not RSS.
    pub fn combined_live_upper_bound_bytes(&self) -> usize {
        self.metadata.combined
    }
    pub fn reused_display(&self) -> bool {
        self.metadata.cached
    }
}
pub(super) struct Cached {
    key: Key,
    mesh: Mesh,
    metadata: Metadata,
}
fn bytes(mesh: &Mesh) -> usize {
    buffers(
        mesh.vertices.capacity(),
        mesh.indices.capacity(),
        mesh.face_indices.capacity(),
    )
}
fn display_bytes(mesh: &DisplayMesh) -> usize {
    buffers(
        mesh.vertices.capacity(),
        mesh.indices.capacity(),
        mesh.face_indices.capacity(),
    )
}
fn buffers(v: usize, i: usize, f: usize) -> usize {
    v.saturating_mul(8)
        .saturating_add(i.saturating_add(f).saturating_mul(4))
}
// This validation layer accepts literal mesh data; it grants no source
// correspondence. Production reaches it only after full affine transport.
pub(super) fn placed_connectivity(target: &DisplayMesh) -> Result<MeshProof, GeometryError> {
    manufacturing::certify(target).map_err(|e| {
        failure(
            GeometryErrorKind::InvalidResult,
            format!("Placed repaired manifold refused: {e:?}"),
        )
    })
}
struct Work<'a> {
    remaining: usize,
    held: usize,
    peak: usize,
    maximum: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Work<'_> {
    fn charge(&mut self, n: usize) -> Result<(), GeometryError> {
        checkpoint(self.interrupt)?;
        self.remaining = self.remaining.checked_sub(n).ok_or_else(|| {
            failure(
                GeometryErrorKind::Limit,
                "Native repaired mesh work exhausted",
            )
        })?;
        Ok(())
    }
    fn storage(&mut self, extra: usize) -> Result<(), GeometryError> {
        self.charge(1)?;
        let total = self.held.saturating_add(extra).saturating_add(SCRATCH);
        if total > self.maximum {
            return Err(failure(
                GeometryErrorKind::Limit,
                "Native repaired mesh storage exhausted",
            ));
        }
        self.peak = self.peak.max(total);
        Ok(())
    }
}
// A cancellation predicate may issue a nested query. Nested reads/builds are
// allowed, but cannot replace the outer query's live cache allocation. This
// keeps its resident-byte charge stable and never holds RefCell borrows across
// caller code. Normal top-level publication remains atomic after the last poll.
struct QueryScope<'a> {
    depth: &'a std::cell::Cell<usize>,
    previous: usize,
}
impl Drop for QueryScope<'_> {
    fn drop(&mut self) {
        self.depth.set(self.previous);
    }
}
impl RepairedPhillipsBody {
    pub(super) fn resident_bytes(&self) -> usize {
        // Source and body each have one Rc allocation. Source retained bytes
        // already include its payload/header; charge both Rc counters/slack.
        // Body includes cache descriptors even when the entry is empty.
        self.source
            .retained_budget_bytes()
            .saturating_add(std::mem::size_of::<Self>() + 8 * std::mem::size_of::<usize>())
            .saturating_add(self.display.borrow().as_ref().map_or(0, |c| bytes(&c.mesh)))
            .saturating_add(
                self.degree11_properties
                    .borrow()
                    .capacity()
                    .saturating_mul(std::mem::size_of::<Option<Degree11LocalFace>>()),
            )
            .saturating_add(self.periodic_properties.borrow().as_ref().map_or(0, |p| {
                p.capacity()
                    * std::mem::size_of::<
                        multicad_rust_kernel::step_read::PhillipsPeriodicFaceProperties,
                    >()
            }))
    }
    pub(in crate::geometry) fn mesh(
        &self,
        frame: PlaneFrame,
        options: Tessellation,
        usage: MeshUse,
    ) -> Result<Mesh, GeometryError> {
        self.controlled(frame, options, usage, Default::default(), None)
            .map(|(mesh, _)| mesh)
    }
    fn controlled(
        &self,
        frame: PlaneFrame,
        options: Tessellation,
        usage: MeshUse,
        budget: RepairedPhillipsMeshBudget,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<(Mesh, Metadata), GeometryError> {
        let previous = self.query_depth.get();
        self.query_depth.set(previous.saturating_add(1));
        let _scope = QueryScope {
            depth: &self.query_depth,
            previous,
        };
        checkpoint(interrupt)?;
        if !(self.policy.linear_tolerance_mm.max(0.001)..=self.policy.maximum_coordinate_mm)
            .contains(&options.deflection)
            || !(0.01..=1.).contains(&options.angle)
            || !(1..=5_000_000).contains(&options.triangle_limit)
        {
            return Err(input(format!(
                "Invalid repaired mesh quality: deflection {} mm, angle {} rad, triangle limit {}",
                options.deflection, options.angle, options.triangle_limit
            )));
        }
        // The caller supplies an upper bound, not a required facet count.
        // Keep this provider's hard ceiling while accepting normal application
        // budgets; chord and angular accuracy are unchanged.
        let options = Tessellation {
            triangle_limit: options.triangle_limit.min(1_000_000),
            ..options
        };
        if budget.maximum_work > WORK
            || budget.maximum_bytes > BYTES
            || options.triangle_limit > 1_000_000
        {
            return Err(input("Invalid repaired mesh hard budget"));
        }
        let mut work = Work {
            remaining: budget.maximum_work,
            held: self.resident_bytes(),
            peak: 0,
            maximum: budget.maximum_bytes,
            interrupt,
        };
        work.storage(0)?;
        work.charge(1024)?; // Fixed two-matrix metric/quality/metadata work.
        let placement = self.placement(frame)?;
        let local = local_quality(options, placement)?;
        let key = Key::new(&self.source, frame, options, local, budget);
        if matches!(usage, MeshUse::Display) {
            // Never hold a RefCell borrow while invoking caller cancellation.
            let descriptor = self
                .display
                .borrow()
                .as_ref()
                .filter(|c| c.key == key)
                .map(|c| {
                    (
                        bytes(&c.mesh),
                        c.mesh
                            .vertices
                            .len()
                            .saturating_add(c.mesh.indices.len())
                            .saturating_add(c.mesh.face_indices.len()),
                    )
                });
            if let Some((held, copy_work)) = descriptor {
                work.storage(held)?;
                work.charge(copy_work)?;
                let saved = self
                    .display
                    .borrow()
                    .as_ref()
                    .filter(|c| c.key == key)
                    .map(|c| (c.mesh.clone(), c.metadata));
                // The matching descriptor was earned before polling. No
                // nested query can publish while this QueryScope is live.
                let (mesh, mut metadata) =
                    saved.expect("QueryScope preserves the matched cache entry");
                work.storage(bytes(&mesh))?;
                checkpoint(interrupt)?;
                metadata.cached = true;
                metadata.query_work = budget.maximum_work - work.remaining;
                metadata.remaining = work.remaining;
                metadata.consumer_peak = work.peak;
                // The cold receipt checked this sum with maximum_bytes;
                // this exact cache key includes the same byte budget.
                metadata.combined = metadata.producer_peak + work.peak + metadata.validator_peak;
                return Ok((mesh, metadata));
            }
        }
        let (mesh, mut metadata) =
            self.cold(frame, options, local, budget, placement, &mut work)?;
        // Optional bounded retention cannot reject an otherwise proved mesh.
        // Clone only when both its peak and copy work fit the existing balance.
        let copy = mesh
            .vertices
            .len()
            .saturating_add(mesh.indices.len())
            .saturating_add(mesh.face_indices.len());
        let mut retain = None;
        if matches!(usage, MeshUse::Display)
            && self.query_depth.get() == 1
            && bytes(&mesh) <= CACHE_BYTES
            && work
                .held
                .saturating_add(2usize.saturating_mul(bytes(&mesh)))
                .saturating_add(SCRATCH)
                <= work.maximum
            && copy.saturating_add(3) <= work.remaining
        {
            work.storage(2usize.saturating_mul(bytes(&mesh)))?;
            work.charge(copy)?;
            let retained = mesh.clone();
            work.storage(bytes(&mesh).saturating_add(bytes(&retained)))?;
            checkpoint(interrupt)?;
            retain = Some(retained);
        }
        checkpoint(interrupt)?;
        metadata.query_work = budget.maximum_work - work.remaining;
        metadata.remaining = work.remaining;
        metadata.consumer_peak = work.peak;
        // cold() checked this sum with maximum_bytes >= every actual peak.
        metadata.combined = metadata.producer_peak + work.peak + metadata.validator_peak;
        // No fallible operation follows publication; failed/cancelled queries
        // therefore retain the old successful display entry unchanged.
        if let Some(retained) = retain {
            *self.display.borrow_mut() = Some(Cached {
                key,
                mesh: retained,
                metadata,
            });
        }
        Ok((mesh, metadata))
    }
    fn cold(
        &self,
        frame: PlaneFrame,
        options: Tessellation,
        local: Tessellation,
        budget: RepairedPhillipsMeshBudget,
        placement: Placement,
        work: &mut Work<'_>,
    ) -> Result<(Mesh, Metadata), GeometryError> {
        #[cfg(test)]
        self.builds.set(self.builds.get() + 1);
        #[cfg(test)]
        let start = std::time::Instant::now();
        // storage(0) earned maximum >= held+SCRATCH. held contains this
        // exact immutable shared source once; QueryScope keeps the old cache
        // resident throughout a nested callback. Reserve Native-only live
        // storage BEFORE the provider allocates its own source buffers.
        let native_only = work.held - self.source.retained_budget_bytes();
        let provider_bytes = work.maximum - native_only - SCRATCH;
        let candidate = self.source.joined_mesh_candidate(
            JoinedMeshOptions {
                tessellation: local,
                maximum_work: work.remaining,
                maximum_bytes: provider_bytes,
            },
            work.interrupt,
        )?;
        work.charge(candidate.work_used())?;
        #[cfg(test)]
        let source_time = start.elapsed();
        #[cfg(test)]
        let start = std::time::Instant::now();
        let producer_peak = candidate.combined_live_upper_bound_bytes();
        let source_work = candidate.work_used();
        let original = candidate.mesh();
        let validator_peak =
            workspace::validators(original.vertices.len() / 3, original.indices.len() / 3)?;
        // Check representability before Native validator allocations. Later
        // actual consumer peaks cannot exceed this request's byte cap.
        workspace::combined(producer_peak, work.maximum, validator_peak)?;
        #[cfg(test)]
        let placement_time;
        #[cfg(test)]
        let mut manifold_time = std::time::Duration::ZERO;
        let (target, proof, transfer) = if placement.is_identity() {
            // The FULL material construction and occurrence maps are exactly
            // identity. This clone preserves every source-proved byte; it is
            // the same manifold, not a detached scalar-based admission.
            work.storage(display_bytes(original).saturating_mul(2))?;
            work.charge(
                original
                    .vertices
                    .len()
                    .saturating_add(original.indices.len())
                    .saturating_add(original.face_indices.len()),
            )?;
            let target = original.clone();
            work.storage(display_bytes(original).saturating_add(display_bytes(&target)))?;
            #[cfg(test)]
            {
                placement_time = start.elapsed();
            }
            (target, candidate.connectivity(), 0.)
        } else {
            let vertices = original.vertices.len() / 3;
            let world_bytes = vertices.saturating_mul(std::mem::size_of::<(Point, f64)>());
            let target_bound = buffers(
                original.vertices.len(),
                original.indices.len(),
                original.face_indices.len(),
            );
            work.storage(
                display_bytes(original)
                    .saturating_add(world_bytes)
                    .saturating_add(target_bound),
            )?;
            let mut world = Vec::with_capacity(vertices);
            let mut target = DisplayMesh {
                vertices: Vec::with_capacity(original.vertices.len()),
                indices: original.indices.clone(),
                face_indices: original.face_indices.clone(),
            };
            work.storage(
                display_bytes(original)
                    .saturating_add(
                        world
                            .capacity()
                            .saturating_mul(std::mem::size_of::<(Point, f64)>()),
                    )
                    .saturating_add(display_bytes(&target)),
            )?;
            let mut transfer = 0f64;
            let point_reserve = (I::point(options.deflection) / I::point(4.)).lo;
            for &point in original.vertices.as_chunks::<3>().0 {
                work.charge(128)?;
                let mapped = placement.point_within(point, point_reserve)?;
                transfer = transfer.max(mapped.1);
                world.push(mapped);
            }
            let normal_reserve = (I::point(options.angle) / I::point(200.)).lo;
            work.charge(affine::NORMAL_PREPARATION_WORK)?;
            let normals = placement.prepare_normals();
            for ids in original.indices.as_chunks::<3>().0 {
                work.charge(affine::NORMAL_TRIANGLE_WORK)?;
                let local =
                    ids.map(|i| std::array::from_fn(|k| original.vertices[3 * i as usize + k]));
                let actual = ids.map(|i| world[i as usize].0);
                normals.require_normal(local, actual, normal_reserve)?;
            }
            for (p, _) in world {
                target.vertices.extend(p);
            }
            work.charge(
                target
                    .face_indices
                    .len()
                    .saturating_mul(52)
                    .saturating_add(vertices),
            )?;
            #[cfg(test)]
            {
                placement_time = start.elapsed();
            }
            #[cfg(test)]
            let start = std::time::Instant::now();
            let proof = placed_connectivity(&target)?;
            #[cfg(test)]
            {
                manifold_time = start.elapsed();
            }
            (target, proof, transfer)
        };
        #[cfg(test)]
        let start = std::time::Instant::now();
        let mesh = Mesh {
            vertices: target.vertices,
            indices: target.indices,
            face_indices: target.face_indices,
        };
        // Source and every placement unit have already consumed this same
        // balance. The literal validator receives ONLY what remains.
        work.charge(1)?; // final-publication unit, polled again below
        let (embedding_pairs, used) = {
            let embedding = mesh.prove_literal_embedding(work.remaining, work.interrupt)?;
            (embedding.pairs(), embedding.work_used())
        };
        work.charge(used)?;
        let volume = self.volume(frame)?;
        let (periodic_literal_volume, periodic_mass_difference) =
            if self.source.periodic().is_some() || self.source.degree11().is_some() {
                work.charge(
                    (mesh.indices.len() / 3)
                        .saturating_mul(LITERAL_VOLUME_TRIANGLE_WORK)
                        .saturating_add(LITERAL_VOLUME_FIXED_WORK),
                )?;
                let literal = periodic::candidate::literal_volume(&mesh)?;
                let difference = I::new(literal[0], literal[1]) - I::new(volume[0], volume[1]);
                (Some(literal), Some(difference.magnitude()))
            } else {
                (None, None)
            };
        checkpoint(work.interrupt)?;
        let used = budget.maximum_work - work.remaining;
        Ok((
            mesh,
            Metadata {
                quality: options,
                occurrence: frame,
                proof,
                volume,
                transfer,
                embedding_pairs,
                source_work,
                proof_work: used,
                proof_remaining: work.remaining,
                query_work: used,
                remaining: work.remaining,
                producer_peak,
                consumer_peak: work.peak,
                validator_peak,
                combined: producer_peak + work.peak + validator_peak,
                cached: false,
                periodic_literal_volume,
                periodic_mass_difference,
                periodic_surface_error: candidate.periodic_surface_error(),
                #[cfg(test)]
                timings: [source_time, placement_time, manifold_time, start.elapsed()],
            },
        ))
    }
    pub(in crate::geometry) fn manufacturing(
        &self,
        frame: PlaneFrame,
        options: MeshOptions,
    ) -> Result<crate::ManufacturingMesh, GeometryError> {
        let quality = Tessellation {
            deflection: options.linear_deflection_mm,
            angle: options.angular_deflection_rad,
            triangle_limit: options.maximum_triangles as usize,
        };
        let (mesh, metadata) = self.controlled(
            frame,
            quality,
            MeshUse::Manufacturing,
            Default::default(),
            None,
        )?;
        Ok(crate::ManufacturingMesh {
            mesh,
            maximum_displacement_mm: 0.,
            edge_pairs: metadata.proof.edge_pairs,
            vertex_links: metadata.proof.vertex_links,
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
impl Solid {
    fn periodic_repaired_mesh_query(
        &self,
        options: MeshOptions,
        budget: RepairedPhillipsMeshBudget,
        usage: MeshUse,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<(Mesh, PeriodicRepairedPhillipsMeshReceipt<'_>), GeometryError> {
        let (body, frame) = self
            .repaired_phillips_body()
            .ok_or_else(|| unsupported("Solid has no periodic repaired source"))?;
        let source = body
            .source
            .periodic()
            .ok_or_else(|| unsupported("Use the retained repaired source mesh API"))?;
        let q = Tessellation {
            deflection: options.linear_deflection_mm,
            angle: options.angular_deflection_rad,
            triangle_limit: options.maximum_triangles as usize,
        };
        let (mesh, metadata) = body.controlled(frame, q, usage, budget, interrupt)?;
        Ok((
            mesh,
            PeriodicRepairedPhillipsMeshReceipt { source, metadata },
        ))
    }
    pub fn rust_periodic_repaired_phillips_mesh_with_receipt(
        &self,
        options: MeshOptions,
        budget: RepairedPhillipsMeshBudget,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<(Mesh, PeriodicRepairedPhillipsMeshReceipt<'_>), GeometryError> {
        self.periodic_repaired_mesh_query(options, budget, MeshUse::Manufacturing, interrupt)
    }
    pub fn rust_periodic_repaired_phillips_display_mesh_with_receipt(
        &self,
        options: MeshOptions,
        budget: RepairedPhillipsMeshBudget,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<(Mesh, PeriodicRepairedPhillipsMeshReceipt<'_>), GeometryError> {
        self.periodic_repaired_mesh_query(options, budget, MeshUse::Display, interrupt)
    }
    fn degree11_repaired_mesh_query(
        &self,
        options: MeshOptions,
        budget: RepairedPhillipsMeshBudget,
        usage: MeshUse,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<(Mesh, Degree11RepairedPhillipsMeshReceipt<'_>), GeometryError> {
        let (body, frame) = self
            .repaired_phillips_body()
            .ok_or_else(|| unsupported("Solid has no degree-eleven repaired source"))?;
        let source = body
            .source
            .degree11()
            .ok_or_else(|| unsupported("Use the matching source mesh API"))?;
        let q = Tessellation {
            deflection: options.linear_deflection_mm,
            angle: options.angular_deflection_rad,
            triangle_limit: options.maximum_triangles as usize,
        };
        let (mesh, metadata) = body.controlled(frame, q, usage, budget, interrupt)?;
        Ok((
            mesh,
            PeriodicRepairedPhillipsMeshReceipt { source, metadata },
        ))
    }
    pub fn rust_degree11_repaired_phillips_mesh_with_receipt(
        &self,
        options: MeshOptions,
        budget: RepairedPhillipsMeshBudget,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<(Mesh, Degree11RepairedPhillipsMeshReceipt<'_>), GeometryError> {
        self.degree11_repaired_mesh_query(options, budget, MeshUse::Manufacturing, interrupt)
    }
    pub fn rust_degree11_repaired_phillips_display_mesh_with_receipt(
        &self,
        options: MeshOptions,
        budget: RepairedPhillipsMeshBudget,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<(Mesh, Degree11RepairedPhillipsMeshReceipt<'_>), GeometryError> {
        self.degree11_repaired_mesh_query(options, budget, MeshUse::Display, interrupt)
    }
    fn repaired_mesh_query(
        &self,
        options: MeshOptions,
        budget: RepairedPhillipsMeshBudget,
        usage: MeshUse,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<(Mesh, RepairedPhillipsMeshReceipt<'_>), GeometryError> {
        let (body, frame) = self
            .repaired_phillips_body()
            .ok_or_else(|| unsupported("Solid has no repaired Phillips source"))?;
        let source = body
            .source
            .retained()
            .ok_or_else(|| unsupported("Use the explicit periodic source mesh API"))?;
        let quality = Tessellation {
            deflection: options.linear_deflection_mm,
            angle: options.angular_deflection_rad,
            triangle_limit: options.maximum_triangles as usize,
        };
        let (mesh, metadata) = body.controlled(frame, quality, usage, budget, interrupt)?;
        Ok((mesh, RepairedPhillipsMeshReceipt { source, metadata }))
    }
    /// Fresh complete repaired-source/placement/embedding proof. Display reuse
    /// is deliberately unavailable on this manufacturing query.
    pub fn rust_repaired_phillips_mesh_with_receipt(
        &self,
        options: MeshOptions,
        budget: RepairedPhillipsMeshBudget,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<(Mesh, RepairedPhillipsMeshReceipt<'_>), GeometryError> {
        self.repaired_mesh_query(options, budget, MeshUse::Manufacturing, interrupt)
    }
    /// Exact-key display reuse retains the original proof ledger separately
    /// from this query's bounded lookup/copy work. Never a manufacturing query.
    pub fn rust_repaired_phillips_display_mesh_with_receipt(
        &self,
        options: MeshOptions,
        budget: RepairedPhillipsMeshBudget,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<(Mesh, RepairedPhillipsMeshReceipt<'_>), GeometryError> {
        self.repaired_mesh_query(options, budget, MeshUse::Display, interrupt)
    }
}
