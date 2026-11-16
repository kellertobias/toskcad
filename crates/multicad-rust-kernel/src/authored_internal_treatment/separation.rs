//! Whole-material disjointness in an explicitly COMMON coaxial source chart.
//! No pose, nominal-fit or sampled-point assertion supplies this authority.
use super::*;
use crate::{
    helical_thread::{ProfileContainmentError, ProfileContainmentOptions},
    thread_entry_chamfer::Chamfer,
};

/// Opaque proof bound to both unchanged actual material operands. Placements
/// are deliberately absent: a consumer must earn their exact common chart.
#[derive(Debug)]
pub struct ChamferSeparation<'a> {
    external: &'a Chamfer,
    receiving: &'a Material,
    reference_gap_mm: f64,
    radial_transfer_mm: f64,
    terminal_gap_mm: Option<f64>,
    gap_mm: f64,
    overlap_end_mm: Enclosure,
    cells: usize,
}
impl ChamferSeparation<'_> {
    pub fn external(&self) -> &Chamfer {
        self.external
    }
    pub fn receiving(&self) -> &Material {
        self.receiving
    }
    pub fn reference_profile_gap_mm(&self) -> f64 {
        self.reference_gap_mm
    }
    pub fn radial_transfer_error_mm(&self) -> f64 {
        self.radial_transfer_mm
    }
    pub fn terminal_radial_gap_mm(&self) -> Option<f64> {
        self.terminal_gap_mm
    }
    /// Source-local radial separation, not a closest3D distance or strength.
    pub fn minimum_radial_gap_mm(&self) -> f64 {
        self.gap_mm
    }
    /// Upper end of the actual common finite interval, whose lower end is0.
    /// The receiving end is its retained exact rational, not a display float.
    pub fn overlap_end_mm(&self) -> Enclosure {
        self.overlap_end_mm
    }
    pub fn cells_visited(&self) -> usize {
        self.cells
    }
}
fn checkpoint(interrupt: Option<&dyn Fn() -> bool>) -> Result<(), Error> {
    if interrupt.is_some_and(|f| f()) {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn profile_error(error: ProfileContainmentError) -> Error {
    match error {
        ProfileContainmentError::Invalid => Error::Invalid,
        ProfileContainmentError::Unresolved => Error::Unresolved,
        ProfileContainmentError::Limit => Error::Limit,
        ProfileContainmentError::Cancelled => Error::Cancelled,
    }
}
impl Material {
    /// Prove empty occupied intersection with an actual external Chamfer in
    /// this SAME coaxial chart. The male source is its authoritative Thread,
    /// so it needs no approximate-profile transfer. The female shared source
    /// is compared through its already-earned, retained legacy correspondence.
    /// A strictly positive gap is required; unresolved contact never admits.
    ///
    /// Thread::prove_profile_containment supplies its bounded65-cell stack and
    /// caller cell cap. Additional work is fixed scalar arithmetic, with no
    /// heap allocation or material reconstruction. Cancellation is checked
    /// around this arithmetic as well as inside every profile cell query.
    pub fn prove_coaxial_chamfer_separation<'a>(
        &'a self,
        external: &'a Chamfer,
        options: ProfileContainmentOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<ChamferSeparation<'a>, Error> {
        checkpoint(interrupt)?;
        let male = external.source_thread();
        let female = self.cavity().source_thread();
        let a = male.spec();
        let b = female.spec();
        // These are exact source-law bindings, never approximate fit labels.
        // The visible prefix uses only the earned finite transfer domain.
        // Terminal material below is checked separately against its ACTUAL
        // cylinder, including any exact receiving end above that domain.
        if a.pitch_mm != b.pitch_mm
            || a.starts != b.starts
            || a.left_hand != b.left_hand
            || a.length_mm != b.length_mm
        {
            return Err(Error::Invalid);
        }
        let (reference, transfer) = female.legacy_approximation();
        let proof = male
            .prove_profile_containment(reference, options, interrupt)
            .map_err(profile_error)?;
        checkpoint(interrupt)?;
        // The ACTUAL authored radial profile linearly interpolates shared
        // endpoints, so it is continuous through every crease and the period
        // seam. Compute its exact-source slope enclosures, not a rounded
        // legacy derivative (legacy one-sided joins can differ by an ulp).
        let mut slope = 0_f64;
        for band in female.bands() {
            checkpoint(interrupt)?;
            let width = I::point(band.phase_mm[1]) - I::point(band.phase_mm[0]);
            let change = I::point(band.radius_mm[1]) - I::point(band.radius_mm[0]);
            let derivative = change / width;
            slope = slope.max(derivative.magnitude());
        }
        // Earned position correspondence bounds radial difference at paired
        // parameters by the reverse triangle inequality. Move the CONTINUOUS
        // actual authored point through its earned angular mismatch to the
        // SAME physical angle as the reference point. Its global angular
        // Lipschitz bound is max|dr/dq| times its TRUE helical inverse gain.
        let rate = I::point(1.) / female.angular_gain();
        let radial = I::point(transfer.maximum_position_error_mm)
            + I::point(slope) * I::point(rate.magnitude()) * I::point(transfer.angular_error_rad);
        let gap = I::point(proof.minimum_gap_mm()) - radial;
        if !radial.valid() || !gap.valid() || gap.lo <= 0. {
            return Err(Error::Unresolved);
        }
        // Never extrapolate the finite authored transfer past its originalL.
        // The constructor's actual terminal cylinder replaces that source
        // after its retained station. Compare it directly to the whole male
        // source's maximum radius; the male MIN/neck cannot enlarge that.
        let terminal_gap = self
            .cavity()
            .tail
            .map(|tail| I::point(tail.radius) - I::point(a.major_radius_mm));
        if terminal_gap.is_some_and(|g| !g.valid() || g.lo <= 0.) {
            return Err(Error::Unresolved);
        }
        let end = self.cavity().end_interval();
        let overlap = Enclosure {
            lower: a.length_mm.min(end.lo),
            upper: a.length_mm.min(end.hi),
        };
        // Constructor-earned set relations complete the material argument:
        // Chamfer's entry MIN and optional neck only REMOVE the original male
        // thread. Cavity's entry MAX and terminal cylinder only ENLARGE its
        // authored void. Material is precisely its finite stock minus that
        // cavity; both source domains begin at localz=0. Outside the actual
        // stock there is no receiving material, even if the male continues.
        checkpoint(interrupt)?;
        Ok(ChamferSeparation {
            external,
            receiving: self,
            reference_gap_mm: proof.minimum_gap_mm(),
            radial_transfer_mm: radial.hi,
            terminal_gap_mm: terminal_gap.map(|g| g.lo),
            gap_mm: gap.lo.min(terminal_gap.map_or(f64::INFINITY, |g| g.lo)),
            overlap_end_mm: overlap,
            cells: proof.cells_visited(),
        })
    }
}

#[cfg(test)]
mod tests;
