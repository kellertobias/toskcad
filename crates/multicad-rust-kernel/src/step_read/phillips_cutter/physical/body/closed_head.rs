//! Same-Body repaired head source and independently measured raw support transfer.
//! Infinite analytic supports are not finite face fillings. Original finite trim
//! matching remains a distinct inherited prerequisite; both closing disks are
//! virtual source geometry. No complete outer/thread/material or mass admission.
use super::*;
use crate::{
    head_profile::{RepairedHeadError, RepairedHeadMeridian, RepairedHeadOptions, Role},
    tapered_thread::interval::I,
};
mod raw_correspondence;
mod supports;
pub use raw_correspondence::{RawHeadCorrespondence, RawHeadFaceCorrespondence};
pub use supports::RawHeadSupportKind;
#[cfg(test)]
pub(super) mod tests;
const WORK: usize = 16_000_000;
const BYTES: usize = 64 * 1024 * 1024;
// All data reads are fixed simple-record/three-coordinate axis reads, not spline
// collectors or cloned ASTs. At most one raw support and its <=4 phase/sign
// candidates are live; the owned source/descriptors are counted separately.
const SCRATCH: usize = 128 * 1024;
#[derive(Clone, Copy)]
pub struct ClosedHeadOptions {
    pub maximum_distance_mm: f64,
    pub maximum_work: usize,
    pub maximum_bytes: usize,
}
impl Default for ClosedHeadOptions {
    fn default() -> Self {
        Self {
            maximum_distance_mm: 1e-7,
            maximum_work: WORK,
            maximum_bytes: BYTES,
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct HeadSupportTransfer {
    /// None denotes the virtual filled top disk on the actual annular plane.
    pub meridian_role: Option<Role>,
    pub material_face_entity: u32,
    pub support_entity: u32,
    /// Raw support class only. The repaired curved face is generally NOT an
    /// exact sphere/torus, especially under the actual near-rigid affine map.
    pub raw_kind: RawHeadSupportKind,
    pub distance_mm: f64,
}
pub struct ClosedRepairedHead<'h, 'q, 'r, 'a, 'input> {
    common: &'h CommonMouthHeadEnvelope<'q, 'r, 'a, 'input>,
    source: RepairedHeadMeridian<'a>,
    bounds: [Option<HeadSupportTransfer>; 6],
    count: usize,
    work: usize,
    retained: usize,
}
impl<'h, 'q, 'r, 'a, 'input> ClosedRepairedHead<'h, 'q, 'r, 'a, 'input> {
    pub fn common_mouth(&self) -> &'h CommonMouthHeadEnvelope<'q, 'r, 'a, 'input> {
        self.common
    }
    pub fn body(&self) -> &'a Body<'input> {
        self.common.closed_cavity().representative().body()
    }
    pub fn source(&self) -> &RepairedHeadMeridian<'a> {
        &self.source
    }
    pub fn chart(&self) -> CommonMouthChart {
        self.common.chart()
    }
    pub fn support_transfers(&self) -> impl Iterator<Item = HeadSupportTransfer> + '_ {
        self.bounds[..self.count].iter().flatten().copied()
    }
    pub fn inherited_finite_matching_admission_mm(&self) -> f64 {
        self.body().profile.rounding_error_mm()
    }
    /// The bottom disk closes ONLY this derived head. It has no associated raw
    /// material face or assertion of a coherent thread/head junction.
    pub fn virtual_bottom_station(&self) -> [f64; 2] {
        self.source.station(self.source.segment_count()).unwrap()
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.retained
    }
}
pub(super) fn repaired_error(e: RepairedHeadError) -> PhillipsReadError {
    match e {
        // Public options were validated before invoking the source constructor.
        RepairedHeadError::Invalid | RepairedHeadError::Unresolved => {
            PhillipsReadError::Unsupported
        }
        RepairedHeadError::Limit => PhillipsReadError::Limit,
        RepairedHeadError::Cancelled => PhillipsReadError::Cancelled,
    }
}
fn read_error(e: ReadError, work: &Work<'_>) -> PhillipsReadError {
    if work.cancelled.get() {
        PhillipsReadError::Cancelled
    } else {
        e.into()
    }
}
impl<'q, 'r, 'a, 'input> CommonMouthHeadEnvelope<'q, 'r, 'a, 'input> {
    pub fn closed_repaired_head<'h>(
        &'h self,
        options: ClosedHeadOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> std::result::Result<ClosedRepairedHead<'h, 'q, 'r, 'a, 'input>, PhillipsReadError> {
        let body = self.closed_cavity().representative().body();
        let mut work = Work {
            remaining: options.maximum_work,
            interrupt,
            cancelled: Cell::new(false),
        };
        let preflight = (|| {
            work.step(1)?;
            if !options.maximum_distance_mm.is_finite()
                || options.maximum_distance_mm <= 0.
                || options.maximum_distance_mm > 1e-7
                || options.maximum_work > WORK
                || options.maximum_bytes > BYTES
            {
                return Err(ReadError::Invalid);
            }
            work.step(self.work_used())?;
            // Common already counts the immutable original Body and all required
            // borrowed receipt payloads once. Self includes the fixed source
            // cache; its Profile borrow refers to that already-counted Body.
            let retained = self
                .retained_budget_bytes()
                .saturating_add(std::mem::size_of::<ClosedRepairedHead<'_, '_, '_, '_, '_>>());
            if retained.saturating_add(SCRATCH) > options.maximum_bytes {
                return Err(ReadError::Limit);
            }
            Ok(retained)
        })()
        .map_err(|e| read_error(e, &work))?;
        let source = body
            .profile
            .closed_repaired_meridian(
                RepairedHeadOptions {
                    maximum_transfer_mm: options.maximum_distance_mm,
                    maximum_work: work.remaining,
                    maximum_bytes: options.maximum_bytes,
                },
                interrupt,
            )
            .map_err(repaired_error)?;
        let result = (|| {
            work.step(source.work_used())?;
            require(body.profile_faces.len() == source.segment_count())?;
            let shell = Shell::new(&body.entities, body.tolerance, true);
            let mut bounds = [None; 6];
            let cap = &body.cap;
            work.step(4096)?;
            let raw = supports::read(&shell, cap.surface)?;
            let top = body.profile.top();
            let distance = supports::top_cap_distance(&raw, self.chart(), top[0], top[1])?;
            require(distance <= options.maximum_distance_mm)?;
            bounds[0] = Some(HeadSupportTransfer {
                meridian_role: None,
                material_face_entity: cap.material_face,
                support_entity: cap.surface,
                raw_kind: RawHeadSupportKind::Plane,
                distance_mm: distance,
            });
            for (index, (s, &face_id)) in
                body.profile.segments().zip(&body.profile_faces).enumerate()
            {
                // Includes fixed record/AXIS2 lookups, complete source interval
                // evaluation, and all <=4 actual whole-support comparisons.
                work.step(4096)?;
                let (surface, _, _) = helical::face(&shell, face_id)?;
                let raw = supports::read(&shell, surface)?;
                let correction = source
                    .segment_certificate(index)
                    .unwrap()
                    .support_transfer_mm;
                let (kind, distance) = supports::distance(&raw, self.chart(), s, correction)?;
                require(distance <= options.maximum_distance_mm)?;
                bounds[index + 1] = Some(HeadSupportTransfer {
                    meridian_role: Some(s.role),
                    material_face_entity: face_id,
                    support_entity: surface,
                    raw_kind: kind,
                    distance_mm: distance,
                });
            }
            work.step(1)?;
            Ok(bounds)
        })();
        let bounds = result.map_err(|e| read_error(e, &work))?;
        Ok(ClosedRepairedHead {
            common: self,
            count: source.segment_count() + 1,
            source,
            bounds,
            work: options.maximum_work - work.remaining,
            retained: preflight,
        })
    }
}
