//! DISTINCT analytic periodic repair material, not the occupied raw STEP set.
//! The head and thread have disjoint open axial interiors. The shared-crest
//! shoulder removes their virtual joining disks and supplies exactly the
//! noncollapsed disk difference, with its original owning complementary arc.
//! The already embedded actual derived cavity is strictly contained in the
//! repaired head except its entire shaped mouth; replace the filled top disk
//! with that shaped-mouth annulus and reverse every non-mouth cavity face.
//! Positive full-rank common affine transport preserves this closure and
//! orientation. Raw material, mass and manufacturing-file admission are separate.
use super::*;
mod volume;
use crate::phillips_drive::{
    FaceKey as CavityFace, RepairedHeadContainment, RepairedHeadContainmentOptions,
};
pub use volume::{
    AuthoredPeriodicOptions, OwnedPeriodicMaterial, PeriodicAuthoredApproximation,
    PeriodicCurveOwner, PeriodicCurveQuery, PeriodicFaceProperties, PeriodicGeometryApproximation,
    PeriodicLiteralCertificate, PeriodicMaterialProvenance, PeriodicMaterialVolume,
    PeriodicOwnedMaterialQuery, PeriodicOwnedOptions, PeriodicProposalOrigin,
    PeriodicVolumeOptions,
};
#[cfg(test)]
mod tests;
const CELLS: usize = 65536;
const METADATA: usize = 32768;
#[derive(Clone, Copy)]
pub struct PeriodicMaterialOptions {
    pub minimum_clearance_mm: f64,
    pub maximum_cells: usize,
    pub maximum_work: usize,
    pub maximum_total_work: usize,
    pub maximum_bytes: usize,
}
impl Default for PeriodicMaterialOptions {
    fn default() -> Self {
        Self {
            minimum_clearance_mm: 1e-7,
            maximum_cells: CELLS,
            maximum_work: WORK,
            maximum_total_work: TOTAL,
            maximum_bytes: BYTES,
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PeriodicMaterialFace {
    HeadSegment(usize),
    ThreadBand(usize),
    ShoulderBand(usize),
    TipCap,
    TopAnnulus,
    ReversedCavity(CavityFace),
}
pub struct PeriodicClosedMaterial<'s, 'input> {
    join: &'s PeriodicShoulderRepair<'s, 's, 's, 'input>,
    containment: RepairedHeadContainment<'s, 's, 's>,
    chart: CommonMouthChart,
    determinant: [f64; 2],
    cavity_faces: [CavityFace; 88],
    count: usize,
    new_work: usize,
    retained: usize,
}
impl<'s, 'input> PeriodicClosedMaterial<'s, 'input> {
    pub fn shoulder(&self) -> &'s PeriodicShoulderRepair<'s, 's, 's, 'input> {
        self.join
    }
    pub fn containment(&self) -> &RepairedHeadContainment<'s, 's, 's> {
        &self.containment
    }
    pub fn chart(&self) -> CommonMouthChart {
        self.chart
    }
    pub fn determinant_bounds(&self) -> [f64; 2] {
        self.determinant
    }
    pub fn new_work_used(&self) -> usize {
        self.new_work
    }
    pub fn inherited_work_used(&self) -> usize {
        self.join.total_work_used()
    }
    pub fn total_work_used(&self) -> usize {
        self.new_work + self.inherited_work_used()
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.retained
    }
    pub fn boundary_faces(&self) -> impl Iterator<Item = PeriodicMaterialFace> + '_ {
        (0..self.join.head_source().segment_count())
            .map(PeriodicMaterialFace::HeadSegment)
            .chain((0..4).map(PeriodicMaterialFace::ThreadBand))
            .chain((0..3).map(PeriodicMaterialFace::ShoulderBand))
            .chain([
                PeriodicMaterialFace::TipCap,
                PeriodicMaterialFace::TopAnnulus,
            ])
            .chain(
                self.cavity_faces[..self.count]
                    .iter()
                    .copied()
                    .map(PeriodicMaterialFace::ReversedCavity),
            )
    }
}
fn cavity_error(e: crate::phillips_drive::Error) -> PhillipsReadError {
    match e {
        crate::phillips_drive::Error::Limit => PhillipsReadError::Limit,
        crate::phillips_drive::Error::Cancelled => PhillipsReadError::Cancelled,
        _ => PhillipsReadError::Unsupported,
    }
}
impl<'h, 'q, 'r, 'a, 'input> ClosedRepairedHead<'h, 'q, 'r, 'a, 'input> {
    pub fn periodic_closed_material<'s>(
        &'s self,
        join: &'s PeriodicShoulderRepair<'s, 's, 'a, 'input>,
        o: PeriodicMaterialOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> FResult<PeriodicClosedMaterial<'s, 'input>>
    where
        'a: 's,
        'r: 's,
    {
        let mut work = Work {
            remaining: o.maximum_work,
            interrupt,
            cancelled: Cell::new(false),
        };
        let result = (|| {
            work.step(1)?;
            if !o.minimum_clearance_mm.is_finite()
                || o.minimum_clearance_mm <= 0.
                || o.maximum_cells == 0
                || o.maximum_cells > CELLS
                || o.maximum_work > WORK
                || o.maximum_total_work > TOTAL
                || o.maximum_bytes > BYTES
            {
                return Err(PhillipsReadError::Invalid);
            }
            if !std::ptr::eq(self.body(), join.finite_curves().thread().body())
                || !std::ptr::eq(self.source(), join.head_source())
                || self.chart() != join.finite_curves().thread().chart()
            {
                return Err(PhillipsReadError::Invalid);
            }
            let allowance = o.maximum_work.min(
                o.maximum_total_work
                    .checked_sub(join.total_work_used())
                    .ok_or(PhillipsReadError::Limit)?,
            );
            work.remaining = allowance.checked_sub(1).ok_or(PhillipsReadError::Limit)?;
            let retained = join
                .retained_budget_bytes()
                .saturating_add(std::mem::size_of::<PeriodicClosedMaterial<'_, '_>>());
            if retained.saturating_add(RepairedHeadContainment::scratch_budget_bytes(
                o.maximum_cells,
            )) > o.maximum_bytes
            {
                return Err(PhillipsReadError::Limit);
            }
            work.step(METADATA)?;
            let closed = self.common_mouth().closed_cavity();
            let local = closed.local_source_certificate();
            let head: &'s RepairedHeadMeridian<'s> = self.source();
            // The join's closed-head proof already includes BOTH immutable
            // child prerequisites: common_mouth -> closed_cavity -> local,
            // and closed_repaired_head -> this exact head source constructor.
            // Preserve the standalone child's16M cap while giving it only the
            // caller's remaining NEW work after those already-paid prefixes.
            let child_prefix = local.work_used() + head.work_used();
            let child_limit = work.remaining.saturating_add(child_prefix).min(WORK);
            let containment = local
                .certify_repaired_head(
                    head,
                    RepairedHeadContainmentOptions {
                        minimum_clearance_mm: o.minimum_clearance_mm,
                        maximum_cells: o.maximum_cells,
                        maximum_work: child_limit,
                        maximum_bytes: o.maximum_bytes,
                    },
                    interrupt,
                )
                .map_err(cavity_error)?;
            let incremental = containment
                .work_used()
                .checked_sub(child_prefix)
                .expect("containment receipt includes its validated immutable prerequisites");
            work.step(incremental)?;
            let mut cavity_faces = [CavityFace::Mouth; 88];
            let mut count = 0;
            for domain in local.patch_domains() {
                work.step(1)?;
                let face = domain.source_face;
                if face != CavityFace::Mouth && !cavity_faces[..count].contains(&face) {
                    cavity_faces[count] = face;
                    count += 1;
                }
            }
            work.step(1)?;
            Ok(PeriodicClosedMaterial {
                join,
                containment,
                chart: self.chart(),
                determinant: closed.determinant_bounds(),
                cavity_faces,
                count,
                new_work: allowance - work.remaining,
                retained,
            })
        })();
        if work.cancelled.get() {
            Err(PhillipsReadError::Cancelled)
        } else {
            result
        }
    }
}
