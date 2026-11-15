//! Explicit approximation proposal for the existing independent indexed mesh
//! producer. The authoritative material and mass remain the owned periodic
//! source. No old source, reader or import receipt is substituted for it.
use super::*;
use crate::phillips_drive::{Cavity, PhillipsScrew};
mod literal;
mod properties;
pub use literal::PeriodicLiteralCertificate;
pub use properties::PeriodicFaceProperties;
#[cfg(test)]
mod tests;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PeriodicProposalOrigin {
    AuthoredSeed,
    DecodedSeed,
}
/// Numerical producer proposal tied to the exact owning periodic source.
/// Authored/decoded seed identity is explicit; neither seed defines geometry.
pub struct PeriodicGeometryApproximation<'s> {
    owner: &'s OwnedPeriodicMaterial,
    origin: PeriodicProposalOrigin,
    proposal: PhillipsScrew,
    position: f64,
    derivative: [f64; 2],
    angle: f64,
    work: usize,
    normal_coefficients: [[[f64; 2]; 2]; 4],
    head_position: f64,
}
/// Compatibility name for the authored-only factory. A decoded proposal is
/// obtained explicitly through mesh_approximation and reports DecodedSeed.
pub type PeriodicAuthoredApproximation<'s> = PeriodicGeometryApproximation<'s>;
impl PeriodicGeometryApproximation<'_> {
    pub fn proposal_origin(&self) -> PeriodicProposalOrigin {
        self.origin
    }
    pub fn source(&self) -> &OwnedPeriodicMaterial {
        self.owner
    }
    pub fn proposal(&self) -> &PhillipsScrew {
        &self.proposal
    }
    /// Continuous proposal fields only. Actual emitted arithmetic is measured
    /// by certify_literal and is never inferred from this bound.
    pub fn maximum_position_error_mm(&self) -> f64 {
        self.position
    }
    pub fn maximum_thread_derivative_error_mm(&self) -> [f64; 2] {
        self.derivative
    }
    pub fn maximum_normal_angle_radians(&self) -> f64 {
        self.angle
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    /// Actual periodic source N=e_r+a*e_theta+b*e_z on each whole band.
    /// These enclosures are conditional on the caller's actual angular cover;
    /// they do not infer an angle from a display mesh or a source header.
    pub fn periodic_outward_normal_coefficients(&self) -> [[[f64; 2]; 2]; 4] {
        self.normal_coefficients
    }
}
fn upper_norm(v: &[I]) -> I {
    v.iter().fold(I::point(0.), |s, x| s + x.square()).sqrt()
}
impl OwnedPeriodicMaterial {
    pub fn authored_mesh_approximation(
        &self,
        maximum_work: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> FResult<PeriodicAuthoredApproximation<'_>> {
        self.geometry_approximation_impl(maximum_work, interrupt, true)
    }
    /// A neutral producer input. The complete source/field correspondence and
    /// actual emitted-literal qualification still run against this owner.
    pub fn mesh_approximation(
        &self,
        maximum_work: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> FResult<PeriodicGeometryApproximation<'_>> {
        self.geometry_approximation_impl(maximum_work, interrupt, false)
    }
    fn geometry_approximation_impl(
        &self,
        maximum_work: usize,
        interrupt: Option<&dyn Fn() -> bool>,
        authored_only: bool,
    ) -> FResult<PeriodicGeometryApproximation<'_>> {
        let mut work = Work {
            remaining: maximum_work,
            interrupt,
            cancelled: Cell::new(false),
        };
        let result = (|| {
            work.step(1)?;
            if maximum_work > WORK {
                return Err(PhillipsReadError::Invalid);
            }
            if authored_only
                && !matches!(self.provenance, PeriodicMaterialProvenance::Authored { .. })
            {
                return Err(PhillipsReadError::Unsupported);
            }
            let (thread, profile, mut source_spec, origin) = match &self.provenance {
                PeriodicMaterialProvenance::Authored {
                    thread,
                    profile,
                    source_spec,
                } => (
                    *thread,
                    *profile,
                    *source_spec,
                    PeriodicProposalOrigin::AuthoredSeed,
                ),
                PeriodicMaterialProvenance::Decoded {
                    thread_proposal, ..
                } => (
                    *thread_proposal,
                    *self.head.view().profile_proposal(),
                    crate::profiled_screw::SourceSpec {
                        head: None,
                        shank_length_mm: 0.,
                        // This lower bound was earned against the actual
                        // decoded head/cavity. It changes no geometric law.
                        minimum_wall_mm: self.containment[0].min(self.containment[1]),
                    },
                    PeriodicProposalOrigin::DecodedSeed,
                ),
            };
            let head = self
                .head
                .view()
                .proposal_transfer(work.remaining.min(WORK), interrupt)
                .map_err(
                    crate::step_read::phillips_cutter::physical::body::closed_head::repaired_error,
                )?;
            work.step(head.work_used)?;
            if origin == PeriodicProposalOrigin::DecodedSeed {
                // This is only a proposal-constructor request, not a material
                // thickness or same-height clearance claim. Leave room for
                // the measured head deviation, then ask for half the remainder.
                // The old constructor independently proves its COMPLETE raw
                // profile/thread/cavity host; this hint cannot bypass that gate.
                let available = I::point(source_spec.minimum_wall_mm) - I::point(head.position_mm);
                require(available.valid() && available.lo > 0.)?;
                source_spec.minimum_wall_mm = (available / I::point(2.)).lo;
            }
            let rotation = rotation::Rotation::new(&mut work)?;
            let real_rate =
                I::point(2.) * rotation.pi(&mut work)? / I::point(self.thread.pitch_mm());
            let old_rate = I::point(std::f64::consts::TAU) / I::point(self.thread.pitch_mm());
            let rate_error = (real_rate - old_rate).magnitude();
            let p = self.thread.phases();
            let r = self.thread.radius_nodes_mm();
            let z = expression(self.caps[0]);
            let z = I::new(z.lo, expression(self.caps[1]).hi);
            let mut position = head.position_mm;
            let mut derivative = [0_f64; 2];
            let mut angle = head.normal_angle_radians;
            let mut normal_coefficients = [[[0.; 2]; 2]; 4];
            for band in 0..4 {
                work.step(2048)?;
                let dp = p[band + 1] - p[band];
                let dr = I::point(r[band + 1]) - I::point(r[band]);
                let exact_slope = dr / dp;
                let radius = I::new(r[band].min(r[band + 1]), r[band].max(r[band + 1]));
                let tangent = I::point(if self.thread.left_hand() { -1. } else { 1. })
                    * exact_slope
                    / (real_rate * radius);
                normal_coefficients[band] = [tangent, -exact_slope].map(|x| [x.lo, x.hi]);
                let slope = I::point(thread.slope_at_phase(p[band].midpoint(), band + 1));
                let anchor = match band {
                    0 => p[0],
                    2 => p[3],
                    _ => p[band],
                };
                let anchor_r = if band == 1 { r[1] } else { r[0] };
                let old = [p[band], p[band + 1]].map(|q| I::point(anchor_r) + slope * (q - anchor));
                let radial_error = (old[0] - I::point(r[band]))
                    .magnitude()
                    .max((old[1] - I::point(r[band + 1])).magnitude());
                let old_r = I::new(old[0].lo.min(old[1].lo), old[0].hi.max(old[1].hi));
                require(old_r.valid() && old_r.lo > 0.)?;
                let old_dr = slope * dp;
                let phi = I::new(p[band].lo, p[band + 1].hi);
                let angle_error = I::point(rate_error) * I::point((z - phi).magnitude());
                let du_angle = real_rate * dp;
                let old_du_angle = old_rate * dp;
                let du_error = I::point((dr - old_dr).magnitude())
                    + I::point(old_dr.magnitude()) * angle_error
                    + I::point(radial_error) * I::point(du_angle.magnitude())
                    + I::point(old_r.hi)
                        * (I::point(rate_error) * dp
                            + I::point(old_du_angle.magnitude()) * angle_error);
                let dz_error = I::point(radial_error) * real_rate
                    + I::point(old_r.hi) * (I::point(rate_error) + old_rate * angle_error);
                let pz = upper_norm(&[I::point(r[0]) * real_rate, I::point(1.)]);
                let old_pu = upper_norm(&[old_dr, I::point(old_r.hi) * old_du_angle]);
                let cross_error = du_error * pz + old_pu * dz_error;
                let min_cross = I::point(old_r.lo) * old_du_angle;
                let q = cross_error / min_cross;
                require(q.valid() && q.hi < 0.25)?;
                let normal = I::point(2.) * q / (I::point(1.) - q);
                let c0 = I::point(radial_error) + I::point(old_r.hi) * angle_error;
                // The immediately dominating valid q<1/4 already proves all
                // component errors finite: cross_error contains each Du/Dz
                // error with strictly positive radial/rate/phase factors.
                // Infinite angle_error also enters Du with a positive radius
                // and positive angular speed. It therefore cannot leave c0
                // finite checking as a new source-facing refusal. Finally
                // normal=2q/(1-q) has a denominator bounded below3/4.
                assert!(
                    [c0, du_error, dz_error, normal]
                        .into_iter()
                        .all(|x| x.valid()),
                    "bounded periodic normal ratio retains finite component errors"
                );
                position = position.max(c0.hi);
                derivative[0] = derivative[0].max(du_error.hi);
                derivative[1] = derivative[1].max(dz_error.hi);
                angle = angle.max(normal.hi);
            }
            // This compares continuous proposal fields only. Actual cached
            // remainder/FMA/trigonometric/cap evaluations are NOT assigned an
            // epsilon allowance; certify_literal measures their emitted data.
            require(position.is_finite() && angle.is_finite())?;
            // These existing producers do not report exact internal operations;
            // debit their entire explicit caller limits before each invocation.
            let cancel = || interrupt.is_some_and(|f| f());
            work.step(1_000_000)?;
            let cavity = Cavity::with_cancel(self.cavity.clone(), 1e-7, 1_000_000, &cancel)
                .map_err(cavity_error)?;
            work.step(1_000_000)?;
            let proposal =
                PhillipsScrew::new(thread, profile, source_spec, cavity, 1_000_000, &cancel)
                    .map_err(cavity_error)?;
            work.step(1)?;
            Ok(PeriodicGeometryApproximation {
                origin,
                owner: self,
                proposal,
                position,
                derivative,
                angle,
                work: maximum_work - work.remaining,
                normal_coefficients,
                head_position: head.position_mm,
            })
        })();
        result.map_err(|e| {
            if work.cancelled.get() {
                PhillipsReadError::Cancelled
            } else {
                e
            }
        })
    }
}
