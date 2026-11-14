//! A DISTINCT closed repaired material source, not exact closure of raw STEP.
//!
//! The repaired head's bottom disk contains the occupied thread's entire star
//! entry disk at their exact common station. Their interiors occupy opposite
//! axial halfspaces. Removing both virtual join disks leaves the plane annulus
//! Gamma(q)*((1-lambda)+lambda*R_head/|Gamma(q)|), with the ORIGINAL polynomial
//! thread rim Gamma, not a fitted circle. Its strict radial gap and the occupied
//! cap's one-wrap polar proof make this an embedded annulus. Its inner uses
//! inherit the removed thread cap's winding, outer uses the removed head cap's.
//! The actual thread tip cap is retained unchanged.
//!
//! The same head's actual derived cavity is strictly contained except its mouth
//! disk. Replace the head top disk by its concentric source annulus, and reverse
//! the cavity side atlas. Its exact mouth radius remains the Cutter expression.
//! The occupied set is the regular-closed difference of (head union thread)
//! and cavity. Their only head/thread intersection is the entry disk; the
//! cavity contacts the outer source only at its mouth disk. Both disks have
//! zero volume and the remaining relative boundary is exactly this inventory.
//! Shared identities/consistent inherited winding thus sew an embedded sphere;
//! outwardness follows the actual head/thread/cavity anchors, +Z on the join
//! annulus and -Z on the top annulus. Contacts/cusps already qualified by those
//! source atlases are retained, without a smooth-boundary claim. This is NOT a
//! raw-exact body, finite manufacturing approximation, Native Solid or file mass.
use super::*;
use crate::{
    head_profile::RepairedHeadMeridian,
    phillips_drive::{
        Cutter, FaceKey as CavityFace, RepairedHeadContainment, RepairedHeadContainmentOptions,
    },
    radial_tensor::RadialBoundary,
    tapered_thread::interval::I,
};
const WORK: usize = 64_000_000;
const BYTES: usize = 64 * 1024 * 1024;
const CELLS: usize = 65536;
// Fixed <=102 pole-owner metadata comparisons plus the complete88-cell source
// face-key inventory. No coordinate scan, retained AST clone or new heap data.
const METADATA_WORK: usize = 32768;
#[derive(Clone, Copy)]
/// Local radial/axial margins only, not Euclidean wall thickness after affine
/// transport. This aggregate proof has a distinct64M hard work cap; the existing
/// carrier16M, containment16M and occupied providers keep their own hard caps.
pub struct ClosedMaterialOptions {
    pub minimum_clearance_mm: f64,
    pub maximum_cells: usize,
    pub maximum_work: usize,
    pub maximum_bytes: usize,
}
impl Default for ClosedMaterialOptions {
    fn default() -> Self {
        Self {
            minimum_clearance_mm: 1e-7,
            maximum_cells: CELLS,
            maximum_work: WORK,
            maximum_bytes: BYTES,
        }
    }
}
/// Actual component role; these local query indices are not persistent names.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MaterialFaceSource {
    HeadMeridian {
        segment: usize,
    },
    ThreadSide {
        face: usize,
    },
    ThreadTip,
    JoinAnnulus,
    TopAnnulus,
    /// Material uses the reverse of the original cavity outward orientation.
    ReversedCavity {
        face: CavityFace,
    },
}
pub struct DerivedClosedMaterial<'s, 'input> {
    body: &'s Body<'input>,
    head: &'s RepairedHeadMeridian<'s>,
    thread: &'s RadialBoundary,
    cavity: &'s Cutter,
    chart: CommonMouthChart,
    determinant: [f64; 2],
    containment: RepairedHeadContainment<'s, 's, 's>,
    join_margin: f64,
    cavity_faces: [CavityFace; 88],
    cavity_count: usize,
    work: usize,
    retained: usize,
}
impl<'s, 'input> DerivedClosedMaterial<'s, 'input> {
    pub fn body(&self) -> &'s Body<'input> {
        self.body
    }
    pub fn head(&self) -> &'s RepairedHeadMeridian<'s> {
        self.head
    }
    pub fn thread(&self) -> &'s RadialBoundary {
        self.thread
    }
    pub fn cavity(&self) -> &'s Cutter {
        self.cavity
    }
    pub fn chart(&self) -> CommonMouthChart {
        self.chart
    }
    pub fn determinant_bounds(&self) -> [f64; 2] {
        self.determinant
    }
    pub fn containment(&self) -> &RepairedHeadContainment<'s, 's, 's> {
        &self.containment
    }
    pub fn join_radial_margin_mm(&self) -> f64 {
        self.join_margin
    }
    /// Known carrier/common-head prerequisite charges plus this query's work.
    /// The already completed occupied provider is a separate bounded operation,
    /// not hidden in this figure or granted a renewed budget by this query.
    pub fn work_used(&self) -> usize {
        self.work
    }
    /// Includes shared Body/common-head once, occupied source/metadata once and
    /// this fixed receipt. Query scratch/RSS/allocator fragmentation are separate.
    pub fn retained_budget_bytes(&self) -> usize {
        self.retained
    }
    pub fn face_count(&self) -> usize {
        self.head.segment_count() + self.thread.side_faces().len() + 3 + self.cavity_count
    }
    pub fn face_source(&self, mut index: usize) -> Option<MaterialFaceSource> {
        let h = self.head.segment_count();
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
    /// The actual input's complete entry coedges, in the removed cap's winding.
    /// The private occupied factory always builds that source cap as face0;
    /// no rounded circle or independently reconstructed wire replaces it.
    pub fn join_inner_coedges(&self) -> &[(usize, bool)] {
        &self.thread.input().faces[0].outer
    }
}
fn cavity_error(e: crate::phillips_drive::Error) -> PhillipsReadError {
    match e {
        crate::phillips_drive::Error::Cancelled => PhillipsReadError::Cancelled,
        crate::phillips_drive::Error::Limit => PhillipsReadError::Limit,
        crate::phillips_drive::Error::Invalid | crate::phillips_drive::Error::Unresolved => {
            PhillipsReadError::Unsupported
        }
    }
}
fn read_error(e: ReadError, work: &Work<'_>) -> PhillipsReadError {
    if work.cancelled.get() {
        PhillipsReadError::Cancelled
    } else {
        e.into()
    }
}
impl<'h, 'q, 'r, 'a, 'input> ClosedRepairedHead<'h, 'q, 'r, 'a, 'input> {
    pub fn closed_repaired_material<'s>(
        &'s self,
        thread: &'s OccupiedDerivedThread<'a, 'input>,
        options: ClosedMaterialOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> std::result::Result<DerivedClosedMaterial<'s, 'input>, PhillipsReadError>
    where
        'a: 's,
        'r: 's,
    {
        let mut work = Work {
            remaining: options.maximum_work,
            interrupt,
            cancelled: Cell::new(false),
        };
        work.step(1).map_err(|e| read_error(e, &work))?;
        if !options.minimum_clearance_mm.is_finite()
            || options.minimum_clearance_mm <= 0.
            || options.maximum_cells == 0
            || options.maximum_cells > CELLS
            || options.maximum_work > WORK
            || options.maximum_bytes > BYTES
        {
            return Err(PhillipsReadError::Invalid);
        }
        let body = self.body();
        let chart = self.chart();
        if !std::ptr::eq(body, thread.body()) || chart != thread.chart() {
            return Err(PhillipsReadError::Unsupported);
        }
        let bottom = self.virtual_bottom_station();
        let join_margin = (I::point(bottom[0]) - I::point(thread.entry_radius_upper_mm())).lo;
        if thread.cap_stations_mm()[0] != bottom[1] || join_margin <= options.minimum_clearance_mm {
            return Err(PhillipsReadError::Unsupported);
        }
        work.step(self.work_used())
            .and_then(|_| work.step(thread.construction_work_used()))
            .and_then(|_| work.step(METADATA_WORK))
            .map_err(|e| read_error(e, &work))?;
        let header = std::mem::size_of::<DerivedClosedMaterial<'_, '_>>();
        let scratch = RepairedHeadContainment::scratch_budget_bytes(options.maximum_cells);
        // Refuse the known prefix before the metadata-only occupied query.
        let prefix = self.retained_budget_bytes().saturating_add(header);
        if prefix.saturating_add(scratch) > options.maximum_bytes {
            return Err(PhillipsReadError::Limit);
        }
        let retained = prefix.saturating_add(thread.additional_budget_bytes());
        if retained.saturating_add(scratch) > options.maximum_bytes {
            return Err(PhillipsReadError::Limit);
        }
        let closed = self.common_mouth().closed_cavity();
        let local = closed.local_source_certificate();
        // This is the head's OWN qualified cavity, never an arbitrary supplied
        // Cutter or independently fitted header. All refs shrink covariantly to
        // the one query lifetime while retaining exact immutable Body identity.
        let head: &'s RepairedHeadMeridian<'s> = self.source();
        let containment: RepairedHeadContainment<'s, 's, 's> = local
            .certify_repaired_head(
                head,
                RepairedHeadContainmentOptions {
                    minimum_clearance_mm: options.minimum_clearance_mm,
                    maximum_cells: options.maximum_cells,
                    maximum_work: work.remaining.min(16_000_000),
                    maximum_bytes: options.maximum_bytes,
                },
                interrupt,
            )
            .map_err(cavity_error)?;
        work.step(containment.work_used())
            .map_err(|e| read_error(e, &work))?;
        let mut cavity_faces = [CavityFace::Mouth; 88];
        let mut count = 0;
        // The earned complete88-domain atlas covers every finite source face.
        // Coalesce only identical source FaceKeys; retain all geometric domains
        // and quotient endpoints through the original closed-source authority.
        for domain in local.patch_domains() {
            work.step(1).map_err(|e| read_error(e, &work))?;
            let face = domain.source_face;
            if face != CavityFace::Mouth && !cavity_faces[..count].contains(&face) {
                cavity_faces[count] = face;
                count += 1;
            }
        }
        work.step(1).map_err(|e| read_error(e, &work))?;
        Ok(DerivedClosedMaterial {
            body,
            head,
            thread: thread.boundary(),
            cavity: local.source(),
            chart,
            determinant: closed.determinant_bounds(),
            containment,
            join_margin,
            cavity_faces,
            cavity_count: count,
            work: options.maximum_work - work.remaining,
            retained,
        })
    }
}

#[cfg(test)]
mod tests;

mod volume;
pub use volume::{MaterialVolume, MaterialVolumeOptions};
