//! Distinct endpoint-coupled head geometry, not exact closure of raw circles.
//!
//! Each retained station is an exact binary64 source coordinate. A circular
//! piece uses the same literal support and its normalized endpoint-ray minor
//! arc, plus c0+(c1-c0)(3t²-2t³). Its endpoints are therefore the shared stations
//! as an expression identity; correction derivatives vanish at both ends.
//! Lines are unchanged. Constructor proofs establish positive radius and strict
//! axial order (with explicit cardinal endpoint factors); horizontal inward
//! bearing segments remain separate. Ordered axial slabs and their shared
//! stations give a simple meridian. Revolution and the two positive-radius
//! cap disks form a closed embedded sphere with the axis quotient at each disk
//! center. For P(t,theta)=(r cos(theta),r sin(theta),z), P_theta x P_t
//! is r*(z' e_r-r' e_z): outward radial sign follows the certified axial
//! order; an inward bearing has outward +Z. Top/bottom disks point -Z/+Z.
//! This is a NEW source, not a raw-file, outer-thread or mass receipt.
use super::{Curve, Profile, Role, Segment};
use crate::tapered_thread::interval::I;
mod arc;
mod query;
use arc::ArcSource;
pub use query::{RepairedHeadLocation, RepairedHeadQuery, RepairedHeadQueryOptions};
#[cfg(test)]
mod tests;
const WORK: usize = 16_000_000;
const BYTES: usize = 64 * 1024 * 1024;
const SCRATCH: usize = 16 * 1024;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RepairedHeadError {
    Invalid,
    Unresolved,
    Limit,
    Cancelled,
}
#[derive(Clone, Copy)]
pub struct RepairedHeadOptions {
    pub maximum_transfer_mm: f64,
    pub maximum_work: usize,
    pub maximum_bytes: usize,
}
impl Default for RepairedHeadOptions {
    fn default() -> Self {
        Self {
            maximum_transfer_mm: 1e-7,
            maximum_work: WORK,
            maximum_bytes: BYTES,
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AxialProgress {
    Strict,
    ZeroAtStart,
    RadialBearing,
}
#[derive(Clone, Copy, Debug)]
pub struct RepairedSegmentCertificate {
    pub role: Role,
    pub stations: [usize; 2],
    /// Lower bound for dz/dt, or dz/dt divided by t for ZeroAtStart.
    /// RadialBearing instead certifies -dr/dt on its horizontal interval.
    pub progress: AxialProgress,
    pub minimum_progress_mm: f64,
    pub minimum_radius_mm: f64,
    /// Measured bound to the same literal circle support; zero for a line.
    pub support_transfer_mm: f64,
    /// Additional finite-arc endpoint-angle transfer, not a C1/unit-t claim.
    pub finite_trim_transfer_mm: f64,
}
#[derive(Clone, Copy)]
enum Part {
    Line(Segment),
    Arc(ArcSource),
}
impl Part {
    fn evaluate(self, t: I) -> [[I; 2]; 2] {
        match self {
            Self::Line(s) => {
                let d =
                    std::array::from_fn::<_, 2, _>(|i| I::point(s.end[i]) - I::point(s.start[i]));
                [std::array::from_fn(|i| I::point(s.start[i]) + t * d[i]), d]
            }
            Self::Arc(a) => a.evaluate(t),
        }
    }
}
/// Immutable source/proof tied to the actual retained Profile proposal. All
/// interval fields are evaluation caches for the explicit source expression;
/// they are never chosen replacement coefficients or caller error allowances.
pub struct RepairedHeadMeridian<'p> {
    profile: &'p Profile,
    parts: [Option<Part>; 5],
    certificates: [Option<RepairedSegmentCertificate>; 5],
    count: usize,
    work: usize,
}
#[derive(Clone, Copy, Debug)]
pub struct RepairedMeridianBounds {
    pub parameter: [f64; 2],
    pub position_mm: [[f64; 2]; 2],
    pub derivative_mm: [[f64; 2]; 2],
}
fn finite(i: I) -> bool {
    i.valid() && i.lo.is_finite() && i.hi.is_finite()
}
fn ensure(b: bool) -> Result<(), RepairedHeadError> {
    if b {
        Ok(())
    } else {
        Err(RepairedHeadError::Unresolved)
    }
}
fn poll(interrupt: Option<&dyn Fn() -> bool>) -> Result<(), RepairedHeadError> {
    if interrupt.is_some_and(|f| f()) {
        Err(RepairedHeadError::Cancelled)
    } else {
        Ok(())
    }
}
impl<'p> RepairedHeadMeridian<'p> {
    pub fn profile_proposal(&self) -> &'p Profile {
        self.profile
    }
    pub fn segment_count(&self) -> usize {
        self.count
    }
    pub fn closed_face_count(&self) -> usize {
        self.count + 2
    }
    pub fn segment_certificate(&self, index: usize) -> Option<RepairedSegmentCertificate> {
        self.certificates.get(index).copied().flatten()
    }
    pub fn station(&self, index: usize) -> Option<[f64; 2]> {
        if index == 0 {
            Some(self.profile.top())
        } else {
            self.profile.segments().nth(index - 1).map(|s| s.end)
        }
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        std::mem::size_of::<Self>() + std::mem::size_of::<Profile>()
    }
    pub fn evaluate(
        &self,
        index: usize,
        parameter: [f64; 2],
    ) -> Result<RepairedMeridianBounds, RepairedHeadError> {
        self.evaluate_with_interrupt(index, parameter, WORK, BYTES, None)
    }
    /// Bounded whole-interval source evaluation. Work conservatively includes
    /// the construction proof; caller limits can only lower the hard caps.
    pub fn evaluate_with_interrupt(
        &self,
        index: usize,
        parameter: [f64; 2],
        maximum_work: usize,
        maximum_bytes: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<RepairedMeridianBounds, RepairedHeadError> {
        poll(interrupt)?;
        if index >= self.count
            || !parameter.into_iter().all(f64::is_finite)
            || parameter[0] < 0.
            || parameter[0] > parameter[1]
            || parameter[1] > 1.
            || maximum_work > WORK
            || maximum_bytes > BYTES
        {
            return Err(RepairedHeadError::Invalid);
        }
        if maximum_work < self.work.saturating_add(257)
            || maximum_bytes < self.retained_budget_bytes() + SCRATCH
        {
            return Err(RepairedHeadError::Limit);
        }
        let value = self.certified_interval_bounds(index, parameter);
        // Whole-unit constructor qualification plus the validated restriction
        // proves finiteness here; repeating that same immutable proof cannot
        // produce a new numerical failure.
        poll(interrupt)?;
        Ok(value)
    }
    /// Internal metered-query handoff. Caller must own a valid segment index
    /// and closed subinterval of[0,1], and charge/poll its fixed evaluation.
    /// Construction qualified Part::evaluate over the WHOLE unit interval;
    /// its outward interval extension, including unit-domain intersections,
    /// is inclusion-monotone. Every such restriction therefore has finite
    /// position/first derivatives. Shared endpoint overrides are exact finite
    /// retained stations, not approximate repairs of evaluated coefficients.
    pub(crate) fn certified_interval_bounds(
        &self,
        index: usize,
        parameter: [f64; 2],
    ) -> RepairedMeridianBounds {
        let [mut p, d] = self.parts[index]
            .unwrap()
            .evaluate(I::new(parameter[0], parameter[1]));
        if parameter == [0., 0.] {
            p = self.station(index).unwrap().map(I::point)
        }
        if parameter == [1., 1.] {
            p = self.station(index + 1).unwrap().map(I::point)
        }
        RepairedMeridianBounds {
            parameter,
            position_mm: p.map(|v| [v.lo, v.hi]),
            derivative_mm: d.map(|v| [v.lo, v.hi]),
        }
    }
}
impl Profile {
    /// Certify a separately named positive-rim closed head source. Exact shared
    /// station identities and whole-domain order earn closure; successful
    /// construction of Profile alone is not used as an embedding assertion.
    pub fn closed_repaired_meridian(
        &self,
        options: RepairedHeadOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<RepairedHeadMeridian<'_>, RepairedHeadError> {
        poll(interrupt)?;
        if !options.maximum_transfer_mm.is_finite()
            || options.maximum_transfer_mm <= 0.
            || options.maximum_transfer_mm > 1e-7
            || options.maximum_work > WORK
            || options.maximum_bytes > BYTES
        {
            return Err(RepairedHeadError::Invalid);
        }
        let retained =
            std::mem::size_of::<RepairedHeadMeridian<'_>>() + std::mem::size_of::<Profile>();
        if options.maximum_bytes < retained + SCRATCH {
            return Err(RepairedHeadError::Limit);
        }
        let mut result = RepairedHeadMeridian {
            profile: self,
            parts: [None; 5],
            certificates: [None; 5],
            count: 0,
            work: 1,
        };
        let separated_caps =
            self.top()[0] > 0. && self.end()[0] > 0. && self.top()[1] < self.end()[1];
        ensure(separated_caps)?;
        // Earn the exact source chain for both retained and authored Profile
        // paths. The canonical query station is never a repair of a mismatch.
        let mut previous = self.top();
        for (index, s) in self.segments().enumerate() {
            poll(interrupt)?;
            // Each fixed source evaluation includes two unit-ray normalizations,
            // two finite-angle enclosures, correction/progress bounds and a
            // whole-interval jet. Profile owns <=5 segments; total cost <=2562.
            result.work += 512;
            if result.work >= options.maximum_work {
                return Err(RepairedHeadError::Limit);
            }
            ensure(s.start == previous)?;
            previous = s.end;
            let (part, mut certificate) = match s.curve {
                Curve::Line => {
                    let dz = I::point(s.end[1]) - I::point(s.start[1]);
                    let dr = I::point(s.start[0]) - I::point(s.end[0]);
                    let (progress, minimum) = if s.start[1] == s.end[1] {
                        ensure(s.role == Role::Bearing && dr.lo > 0.)?;
                        (AxialProgress::RadialBearing, dr.lo)
                    } else {
                        ensure(dz.lo > 0.)?;
                        (AxialProgress::Strict, dz.lo)
                    };
                    let radius = s.start[0].min(s.end[0]);
                    ensure(radius > 0.)?;
                    (
                        Part::Line(s),
                        RepairedSegmentCertificate {
                            role: s.role,
                            stations: [index, index + 1],
                            progress,
                            minimum_progress_mm: minimum,
                            minimum_radius_mm: radius,
                            support_transfer_mm: 0.,
                            finite_trim_transfer_mm: 0.,
                        },
                    )
                }
                Curve::Arc {
                    center,
                    radius_mm,
                    sweep,
                } => {
                    let (arc, certificate) = ArcSource::new(s, center, radius_mm, sweep)?;
                    (Part::Arc(arc), certificate)
                }
            };
            // An initial inward horizontal bearing would overlap the top disk.
            // Internal/end bearings instead meet disjoint axial slabs/cap annuli.
            ensure(index != 0 || certificate.progress != AxialProgress::RadialBearing)?;
            certificate.stations = [index, index + 1];
            ensure(
                (I::point(certificate.support_transfer_mm)
                    + I::point(certificate.finite_trim_transfer_mm))
                .hi <= options.maximum_transfer_mm,
            )?;
            // Qualify every intermediate denominator/result over the complete
            // unit interval once. Subsequent subinterval arithmetic is inclusive.
            let finite_jet = part
                .evaluate(I::new(0., 1.))
                .into_iter()
                .flatten()
                .all(finite);
            ensure(finite_jet)?;
            result.parts[index] = Some(part);
            result.certificates[index] = Some(certificate);
            result.count += 1;
        }
        ensure(previous == self.end())?;
        result.work += 1;
        poll(interrupt)?;
        Ok(result)
    }
}

mod volume;
pub use volume::{RepairedHeadVolume, RepairedHeadVolumeOptions};

/// Owned snapshot of the actual qualified repaired source expression. It does
/// not reconstruct a meridian from a recipe or retain a self-referential view.
pub struct OwnedRepairedHead {
    profile: Profile,
    parts: [Option<Part>; 5],
    certificates: [Option<RepairedSegmentCertificate>; 5],
    count: usize,
    work: usize,
}
impl RepairedHeadMeridian<'_> {
    pub(crate) fn owned_snapshot(&self) -> OwnedRepairedHead {
        OwnedRepairedHead {
            profile: *self.profile,
            parts: self.parts,
            certificates: self.certificates,
            count: self.count,
            work: self.work,
        }
    }
}
impl OwnedRepairedHead {
    /// Temporary immutable view through the original source evaluator. All
    /// stations, corrected arc expressions and certificates are copied exactly.
    pub fn view(&self) -> RepairedHeadMeridian<'_> {
        RepairedHeadMeridian {
            profile: &self.profile,
            parts: self.parts,
            certificates: self.certificates,
            count: self.count,
            work: self.work,
        }
    }
}
