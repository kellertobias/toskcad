//! Owned actual repaired source; no raw-exact or manufacturing admission.
use super::*;
use crate::head_profile::OwnedRepairedHead;
use std::sync::Arc;
const SCRATCH: usize = 16384;
mod mesh;
pub use mesh::{DelegatedWorkspaceBounds, JoinedMaterialCandidate, JoinedMeshOptions};
#[derive(Clone, Copy)]
pub struct OwnedMaterialOptions {
    pub maximum_work: usize,
    pub maximum_bytes: usize,
}
impl Default for OwnedMaterialOptions {
    fn default() -> Self {
        Self {
            maximum_work: 512_000_000,
            maximum_bytes: BYTES,
        }
    }
}
/// Independent immutable geometry and provenance. No borrowed parser arena,
/// repaired source view or prerequisite receipt survives in this value.
/// Manufacturing, finite facet embedding and Native admission remain separate.
pub struct OwnedDerivedMaterial {
    raw: Arc<str>,
    head: OwnedRepairedHead,
    thread: RadialBoundary,
    cavity: Cutter,
    chart: CommonMouthChart,
    determinant: [f64; 2],
    join_margin: f64,
    containment: [f64; 2],
    cavity_faces: [CavityFace; 88],
    cavity_count: usize,
    local_volume: [f64; 2],
    local_enclosure: [[f64; 2]; 3],
    components: [[f64; 2]; 3],
    work: usize,
}
impl OwnedDerivedMaterial {
    pub fn original_input(&self) -> &str {
        &self.raw
    }
    pub fn head(&self) -> &OwnedRepairedHead {
        &self.head
    }
    pub fn thread(&self) -> &RadialBoundary {
        &self.thread
    }
    pub fn cavity(&self) -> &Cutter {
        &self.cavity
    }
    pub fn chart(&self) -> CommonMouthChart {
        self.chart
    }
    pub fn determinant_bounds(&self) -> [f64; 2] {
        self.determinant
    }
    pub fn join_radial_margin_mm(&self) -> f64 {
        self.join_margin
    }
    /// Certified local radial and axial-floor margins, respectively.
    pub fn containment_clearances_mm(&self) -> [f64; 2] {
        self.containment
    }
    pub fn local_volume_mm3(&self) -> [f64; 2] {
        self.local_volume
    }
    /// Conservative finite envelope of the actual repaired head/thread union.
    /// The contained cavity only removes material. Construction earned this
    /// from actual head Part fields and typed thread Bernstein/cap hulls.
    pub fn local_enclosure_bounds(&self) -> [[f64; 2]; 3] {
        self.local_enclosure
    }
    pub fn component_volumes_mm3(&self) -> [[f64; 2]; 3] {
        self.components
    }
    /// Original common-chart volume, recomputed from LOCAL authority and the
    /// actual four-column determinant. This API has no occurrence-pose setter;
    /// a future placed owner must certify its full composed affine determinant.
    pub fn world_volume_mm3(&self) -> [f64; 2] {
        let v = I::new(self.local_volume[0], self.local_volume[1])
            * I::new(self.determinant[0], self.determinant[1]);
        [v.lo, v.hi]
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        std::mem::size_of::<Self>()
            .saturating_add(
                self.thread
                    .retained_budget_bytes()
                    .saturating_sub(std::mem::size_of::<RadialBoundary>()),
            )
            .saturating_add(self.raw.len())
            .saturating_add(2 * std::mem::size_of::<usize>() + 2 * std::mem::align_of::<usize>())
    }
    pub fn face_count(&self) -> usize {
        self.head.view().segment_count() + self.thread.side_faces().len() + 3 + self.cavity_count
    }
    pub fn face_source(&self, mut index: usize) -> Option<MaterialFaceSource> {
        let h = self.head.view().segment_count();
        if index < h {
            return Some(MaterialFaceSource::HeadMeridian { segment: index });
        }
        index -= h;
        let sides = self.thread.side_faces();
        if index < sides.len() {
            return Some(MaterialFaceSource::ThreadSide { face: sides[index] });
        }
        index -= sides.len();
        match index {
            0 => Some(MaterialFaceSource::ThreadTip),
            1 => Some(MaterialFaceSource::JoinAnnulus),
            2 => Some(MaterialFaceSource::TopAnnulus),
            _ => self.cavity_faces[..self.cavity_count]
                .get(index - 3)
                .copied()
                .map(|face| MaterialFaceSource::ReversedCavity { face }),
        }
    }
    pub fn join_inner_coedges(&self) -> &[(usize, bool)] {
        &self.thread.input().faces[0].outer
    }
}
fn check_peak(
    prefix: usize,
    payload: usize,
    maximum: usize,
) -> std::result::Result<(), PhillipsReadError> {
    if prefix.saturating_add(payload) > maximum {
        Err(PhillipsReadError::Limit)
    } else {
        Ok(())
    }
}
// Final publication validates the ACTUAL retained capacities independently of
// the construction preflight, then polls cancellation before returning ownership.
fn publish(
    mut result: OwnedDerivedMaterial,
    prefix: usize,
    options: OwnedMaterialOptions,
    work: &mut Work<'_>,
) -> std::result::Result<OwnedDerivedMaterial, PhillipsReadError> {
    check_peak(
        prefix,
        result.retained_budget_bytes(),
        options.maximum_bytes,
    )?;
    work.step(1).map_err(|e| read_error(e, work))?;
    result.work = options.maximum_work - work.remaining;
    Ok(result)
}
impl MaterialVolume<'_, '_, '_> {
    /// Conservative full live peak while the original source and receipt still
    /// exist. Includes copied source payload, original input Arc and fixed
    /// receipt/clone scratch; logical accounted bytes, not allocator RSS.
    pub fn owned_peak_budget_bytes(&self) -> usize {
        self.source()
            .retained_budget_bytes()
            .saturating_add(SCRATCH)
            .saturating_add(std::mem::size_of::<OwnedDerivedMaterial>())
            .saturating_add(2usize.saturating_mul(self.source().thread.retained_budget_bytes()))
            .saturating_add(self.source().body.literal_input().len())
            .saturating_add(64)
    }

    /// Consume an earned same-source volume receipt and copy the actual source
    /// payload. The original source remains unchanged on any refusal/cancel.
    pub fn into_owned(
        self,
        options: OwnedMaterialOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> std::result::Result<OwnedDerivedMaterial, PhillipsReadError> {
        let source = self.source();
        let mut work = Work {
            remaining: options.maximum_work,
            interrupt,
            cancelled: Cell::new(false),
        };
        work.step(1).map_err(|e| read_error(e, &work))?;
        if options.maximum_work > 512_000_000 || options.maximum_bytes > BYTES {
            return Err(PhillipsReadError::Invalid);
        }
        let raw = source.body.literal_input();
        let thread_bytes = source.thread.retained_budget_bytes();
        // Clone retains Arc pole owners but copies nested Vec payloads. Charging
        // twice the entire old source covers every cloned Vec's capacity and
        // transient allocation conservatively, including headers/shared poles.
        // Raw Arc<str> has exact slice storage plus counters/alignment slack.
        let candidate_bound = std::mem::size_of::<OwnedDerivedMaterial>()
            .saturating_add(2usize.saturating_mul(thread_bytes))
            .saturating_add(raw.len())
            .saturating_add(64);
        let prefix = source.retained_budget_bytes().saturating_add(SCRATCH);
        check_peak(prefix, candidate_bound, options.maximum_bytes)?;
        // Source construction bounds controls/counts; a byte-linear charge also
        // covers all nested descriptor copies and immutable raw input transfer.
        work.step(self.work_used())
            .and_then(|_| work.step(candidate_bound))
            .map_err(|e| read_error(e, &work))?;
        let mut local_enclosure = source.thread_enclosure;
        for part in 0..source.head.segment_count() {
            // Sealed valid index/whole-unit domain handoff; construction
            // already qualified finite Part fields. This is no nominal Profile.
            work.step(257).map_err(|e| read_error(e, &work))?;
            let field = source.head.certified_interval_bounds(part, [0., 1.]);
            let radius = field.position_mm[0][1];
            local_enclosure[0][0] = local_enclosure[0][0].min(-radius);
            local_enclosure[0][1] = local_enclosure[0][1].max(radius);
            local_enclosure[1][0] = local_enclosure[1][0].min(-radius);
            local_enclosure[1][1] = local_enclosure[1][1].max(radius);
            local_enclosure[2][0] = local_enclosure[2][0].min(field.position_mm[1][0]);
            local_enclosure[2][1] = local_enclosure[2][1].max(field.position_mm[1][1]);
        }
        let head = source.head.owned_snapshot();
        let thread = source.thread.clone();
        work.step(1).map_err(|e| read_error(e, &work))?;
        let raw = Arc::<str>::from(raw);
        let result = OwnedDerivedMaterial {
            raw,
            head,
            thread,
            cavity: source.cavity.clone(),
            chart: source.chart,
            determinant: source.determinant,
            join_margin: source.join_margin,
            containment: [
                source.containment.minimum_radial_clearance_mm(),
                source.containment.axial_floor_clearance_mm(),
            ],
            cavity_faces: source.cavity_faces,
            cavity_count: source.cavity_count,
            local_volume: self.local_volume_mm3(),
            local_enclosure,
            components: self.component_volumes_mm3(),
            work: 0,
        };
        publish(result, prefix, options, &mut work)
    }
}
#[cfg(test)]
mod tests;

#[cfg(test)]
mod mesh_probe;
