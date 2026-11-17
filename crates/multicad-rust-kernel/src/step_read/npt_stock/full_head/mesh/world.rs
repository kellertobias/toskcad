//! Actual file B, then ONLY additional Native occurrence A. In particular B
//! already contains any encoded intrinsic screw F(L); never append it twice.
use super::*;
use crate::step_read::npt_stock::{
    LiteralOccurrencePlacement, full_head::transport::MappedActualHead,
};
#[cfg(test)]
mod tests;
#[derive(Clone, Copy, Debug)]
pub struct MappedHeadMeshReceipt {
    pub local: HeadMeshReceipt,
    pub inherited_mapped_material_work: usize,
    pub inherited_occurrence_work: usize,
    pub world_transfer_work: usize,
    pub aggregate_work: usize,
    pub accounted_peak_bytes: usize,
    pub source_position_error_mm: f64,
    pub literal_world_rounding_mm: f64,
    pub normal_tangent_cone: f64,
    pub remaining_work: usize,
}
pub struct MappedHeadMeshCandidate<'source> {
    source: &'source MappedActualHead,
    occurrence: Option<&'source LiteralOccurrencePlacement>,
    mesh: DisplayMesh,
    connectivity: MeshProof,
    receipt: MappedHeadMeshReceipt,
}
impl<'source> MappedHeadMeshCandidate<'source> {
    pub fn source(&self) -> &'source MappedActualHead {
        self.source
    }
    pub fn occurrence(&self) -> Option<&'source LiteralOccurrencePlacement> {
        self.occurrence
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
    pub fn receipt(&self) -> MappedHeadMeshReceipt {
        self.receipt
    }
}
impl MappedActualHead {
    pub fn mesh_candidate(
        &self,
        options: HeadMeshOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<MappedHeadMeshCandidate<'_>> {
        self.mesh_with_occurrence(None, options, interrupt)
    }
    pub fn placed_mesh_candidate<'source>(
        &'source self,
        occurrence: &'source LiteralOccurrencePlacement,
        options: HeadMeshOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<MappedHeadMeshCandidate<'source>> {
        self.mesh_with_occurrence(Some(occurrence), options, interrupt)
    }
    fn mesh_with_occurrence<'source>(
        &'source self,
        occurrence: Option<&'source LiteralOccurrencePlacement>,
        options: HeadMeshOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<MappedHeadMeshCandidate<'source>> {
        options.tessellation.check().map_err(|_| Error::Invalid)?;
        if options.maximum_work == 0
            || options.maximum_work > MAX_MESH_WORK
            || options.maximum_bytes == 0
            || options.maximum_bytes > MAX_BYTES
        {
            return Err(Error::Invalid);
        }
        let source = self.local_geometry();
        let shared = source.receipt().aggregate_work;
        let old = self.receipt();
        let occurrence_work = occurrence.map_or(0, |a| a.work_used());
        let occurrence_bytes = occurrence.map_or(0, |a| a.retained_budget_bytes());
        let inherited = old
            .aggregate_work
            .checked_add(occurrence_work)
            .ok_or(Error::Limit)?;
        let header = size_of::<MappedHeadMeshCandidate<'_>>() + HEAD_SCRATCH;
        let base = old
            .retained_bytes
            .checked_add(occurrence_bytes)
            .and_then(|n| n.checked_add(options.additional_resident_bytes))
            .and_then(|n| n.checked_add(header))
            .ok_or(Error::Limit)?;
        let inherited_peak = old
            .accounted_peak_bytes
            .checked_sub(old.additional_resident_bytes)
            .and_then(|n| n.checked_add(occurrence_bytes))
            .and_then(|n| n.checked_add(options.additional_resident_bytes))
            .and_then(|n| n.checked_add(header))
            .ok_or(Error::Limit)?;
        let occurrence_peak = occurrence.map_or(Ok(0), |a| {
            a.accounted_peak_bytes()
                .checked_add(old.retained_bytes)
                .and_then(|n| n.checked_add(options.additional_resident_bytes))
                .and_then(|n| n.checked_add(header))
                .ok_or(Error::Limit)
        })?;
        if base > options.maximum_bytes
            || inherited_peak.max(occurrence_peak) > options.maximum_bytes
        {
            return Err(Error::Limit);
        }
        let mut w = Work {
            used: inherited,
            remaining: options
                .maximum_work
                .checked_sub(inherited)
                .ok_or(Error::Limit)?,
            base,
            held: 0,
            peak: base.max(inherited_peak).max(occurrence_peak),
            maximum_bytes: options.maximum_bytes,
            interrupt,
        };
        w.step(128)?;
        let delta = occurrence.map_or(0., |a| a.metric_error());
        if !delta.is_finite() || !(0.0..1.0).contains(&delta) {
            return Err(Error::Unsupported);
        }
        let frame = self.actual_frame();
        let scale = frame.scale_bounds();
        let sigma_a = occurrence.map_or(1., |a| a.linear_l2_norm_upper());
        let sigma = (I::point(sigma_a) * I::point(scale[1])).hi;
        let local_deflection = (I::point(options.tessellation.deflection) / I::point(sigma)).lo;
        let low = (I::point(1.) - I::point(delta)).square();
        let high = (I::point(1.) + I::point(delta)).square();
        let wanted = I::point(options.tessellation.angle)
            - I::point(options.tessellation.angle) / I::point(4096.);
        let local_angle = (low * wanted / (high + (high - low) * wanted / I::point(2.)))
            .lo
            .max(0.01)
            .min(options.tessellation.angle);
        let extra_work = inherited
            .checked_sub(shared)
            .and_then(|n| n.checked_add(w.used - inherited))
            .ok_or(Error::Limit)?;
        let outside = base
            .checked_sub(source.receipt().retained_bytes)
            .ok_or(Error::Limit)?;
        let local = source.mesh_candidate(
            HeadMeshOptions {
                tessellation: Tessellation {
                    deflection: local_deflection,
                    angle: local_angle,
                    ..options.tessellation
                },
                maximum_work: options
                    .maximum_work
                    .checked_sub(extra_work)
                    .ok_or(Error::Limit)?,
                maximum_bytes: options.maximum_bytes,
                additional_resident_bytes: outside,
            },
            interrupt,
        )?;
        let lr = local.receipt();
        w.step(lr.aggregate_work.checked_sub(shared).ok_or(Error::Limit)?)?;
        w.peak = w.peak.max(lr.accounted_peak_bytes);
        let original = local.mesh();
        let input_bytes = size_of::<ActualHeadMeshCandidate<'_>>()
            + original.vertices.capacity() * 8
            + (original.indices.capacity() + original.face_indices.capacity()) * 4;
        w.held = input_bytes;
        w.storage(0)?;
        let mut world = DisplayMesh {
            vertices: w.allocate(original.vertices.len())?,
            indices: Vec::new(),
            face_indices: Vec::new(),
        };
        let mut rounding = 0f64;
        for &p in original.vertices.as_chunks::<3>().0 {
            w.step(128)?;
            let mut exact = frame.map_box_mm(p.map(|x| [x, x]))?;
            if let Some(a) = occurrence {
                exact = a.map_box_mm(exact)?;
            }
            let interval = exact.map(|x| I::new(x[0], x[1]));
            let point = interval.map(I::midpoint);
            if point.iter().any(|x| !x.is_finite()) {
                return Err(Error::Unsupported);
            }
            let error = (0..3).fold(I::point(0.), |s, k| {
                s + I::point((interval[k] - I::point(point[k])).magnitude())
            });
            if !error.valid() {
                return Err(Error::Unsupported);
            }
            rounding = rounding.max(error.hi);
            world.vertices.extend(point);
        }
        let distance =
            (I::point(sigma) * I::point(lr.source_position_error_mm) + I::point(rounding)).hi;
        if !distance.is_finite() || distance > options.tessellation.deflection {
            return Err(Error::Unsupported);
        }
        let transported = normal::transported_metric(delta, lr.normal_tangent_cone, &mut w)?;
        let mut angle = 0f64;
        let mut flux = I::point(0.);
        let anchor = std::array::from_fn::<_, 3, _>(|i| I::point(world.vertices[i]));
        for &ids in original.indices.as_chunks::<3>().0 {
            let raw = [
                vertex(original, ids[0])?,
                vertex(original, ids[1])?,
                vertex(original, ids[2])?,
            ];
            let actual = [
                vertex(&world, ids[0])?,
                vertex(&world, ids[1])?,
                vertex(&world, ids[2])?,
            ];
            let (round, normal) =
                normal::world_rounding(raw, actual, frame.axes(), occurrence, &mut w)?;
            let bound = normal::combine(transported, round, &mut w)?;
            if bound > options.tessellation.angle {
                return Err(Error::Unsupported);
            }
            angle = angle.max(bound);
            w.step(32)?;
            let at = std::array::from_fn::<_, 3, _>(|i| I::point(actual[0][i]) - anchor[i]);
            flux = flux + (0..3).fold(I::point(0.), |s, i| s + at[i] * normal[i]) / I::point(6.);
        }
        if !flux.valid() || flux.lo <= 0. {
            return Err(Error::Unsupported);
        }
        let proof = MeshProof {
            edge_pairs: local.connectivity().edge_pairs,
            vertex_links: local.connectivity().vertex_links,
            signed_volume_mm3: flux.midpoint(),
        };
        // SAME indices and finite nonzero positively mapped actual facets
        // retain indexed edge pairing/links; world flux above was earned from
        // SAME world points. Whole literal embedding is STILL required.
        let raw = local.into_mesh();
        world.indices = raw.indices;
        world.face_indices = raw.face_indices;
        w.held -= size_of::<ActualHeadMeshCandidate<'_>>() + raw.vertices.capacity() * 8;
        drop(raw.vertices);
        w.step(0)?;
        let world_work = w
            .used
            .checked_sub(inherited)
            .and_then(|n| n.checked_sub(lr.aggregate_work - shared))
            .ok_or(Error::Limit)?;
        Ok(MappedHeadMeshCandidate {
            source: self,
            occurrence,
            mesh: world,
            connectivity: proof,
            receipt: MappedHeadMeshReceipt {
                local: lr,
                inherited_mapped_material_work: old.aggregate_work,
                inherited_occurrence_work: occurrence_work,
                world_transfer_work: world_work,
                aggregate_work: w.used,
                accounted_peak_bytes: w.peak,
                source_position_error_mm: distance,
                literal_world_rounding_mm: rounding,
                normal_tangent_cone: angle,
                remaining_work: w.remaining,
            },
        })
    }
}
