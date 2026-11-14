//! Actual canonical metric-root STEP SUPPORT carrier. This is not a closed
//! body, trimmed-face, independently decoded material or export certificate.
use super::insert::{
    AngularCarrier, fit_angular,
    resources::{Memory, Work},
};
pub use super::insert::{CurveOptions, Error};
use crate::{
    authored_metric_thread::{AuthoredMetricThread, RationalMeridian},
    tapered_thread::interval::I,
};
mod domain;
pub use domain::{DomainOptions, RootFaceDomain, TrimBounds};
mod periodic;
pub use periodic::PeriodicAngular;
mod encoded;
pub use encoded::{
    CapEndpoint, CapTrace, EncodedCell, EncodedControl, EncodedField, EncodedMeridian,
    EncodedScalar,
};
mod body_field;
pub use body_field::{BodyField, CellSurface, FieldBounds, FieldControl};
mod field;
pub use field::{
    BandPiece, CellSeam, HomogeneousControl, MeridianRow, MeridianRowKey, PhaseNode, SeamPoint,
    SourceField, StationKey,
};

pub struct RootSurface {
    source: RationalMeridian,
    turns: [f64; 2],
    angular: AngularCarrier,
    poles: Vec<[f64; 3]>,
    bound: SupportBounds,
    work: usize,
    bytes: usize,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SupportBounds {
    pub position_error_mm: f64,
    /// Whole normalized root-u and angular-v parameter derivatives.
    pub derivative_error_mm: [f64; 2],
    /// Angle between actual literal/source tangent crosses in the SAME order.
    pub normal_angle_rad: f64,
}
const SCRATCH: usize = 16 * 1024;
const WEIGHTS: [f64; 3] = [1., 0.5, 1.];
fn p(x: f64) -> I {
    I::point(x)
}
fn hull(x: [f64; 2]) -> I {
    I::new(x[0], x[1])
}
fn l1(v: [I; 3]) -> I {
    v.into_iter().fold(p(0.), |sum, x| sum + p(x.magnitude()))
}
fn ideal(h: [I; 3], angular: [f64; 3], lead: I) -> [I; 3] {
    [
        h[1] * p(angular[0]),
        h[1] * p(angular[1]),
        h[0] + h[2] * lead * p(angular[2]),
    ]
}
impl RootSurface {
    /// Root coordinate uses the same owned rational meridian; the second
    /// coordinate t lifts physical turns=v0+(v1-v0)t. Support outside the
    /// source's end planes must later be trimmed; it is not extra material.
    pub fn prepare(
        source: AuthoredMetricThread,
        turns: [f64; 2],
        options: CurveOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        let mut w = Work::new(options.work, interrupt)?;
        w.charge(1)?;
        options.validate()?;
        let mut memory = Memory::new(size_of::<Self>(), options.bytes)?;
        memory.add(SCRATCH)?;
        let meridian = source.root_meridian();
        w.charge(256)?;
        let profile = meridian.whole_bounds();
        let lead = p(source.spec().pitch_mm) * p(f64::from(source.spec().starts));
        let radius = p(source.maximum_radius_mm());
        let dr = hull(profile.first_derivative_mm[1]).magnitude();
        let scale = (radius + lead + p(dr) + p(1.)).hi;
        let angular_options = CurveOptions {
            position_error_mm: options.position_error_mm / (16. * scale),
            derivative_error_mm: options.derivative_error_mm / (16. * scale),
            ..options
        };
        if angular_options.position_error_mm == 0. || angular_options.derivative_error_mm == 0. {
            return Err(Error::Unresolved);
        }
        let angular = fit_angular(
            turns,
            source.spec().left_hand,
            angular_options,
            &mut memory,
            &mut w,
        )?;
        // Sealed fitter emits <=11*1024+1 controls. Even a 16-bit usize
        // contains all three rows; allocation itself remains checked below.
        let count = 3 * angular.poles.len();
        let mut poles = memory.allocate(count, &mut w)?;
        w.charge(256)?;
        let controls = meridian.homogeneous_controls().map(|row| row.map(hull));
        for h in controls {
            for &a in &angular.poles {
                w.charge(32)?;
                poles.push(ideal(h, a, lead).map(|x| (x / h[2]).midpoint()));
            }
        }
        let mut result = Self {
            source: meridian,
            turns,
            angular,
            poles,
            bound: SupportBounds {
                position_error_mm: 0.,
                derivative_error_mm: [0.; 2],
                normal_angle_rad: 0.,
            },
            work: 0,
            bytes: 0,
        };
        result.bound = result.certify(&result.poles, options, &mut w)?;
        w.charge(1)?;
        memory.used -= SCRATCH;
        result.work = w.used;
        result.bytes = memory.used;
        Ok(result)
    }
    fn certify(
        &self,
        poles: &[[f64; 3]],
        options: CurveOptions,
        w: &mut Work<'_>,
    ) -> Result<SupportBounds, Error> {
        let count = self.angular.poles.len();
        if poles.len() != 3 * count {
            return Err(Error::Invalid);
        }
        w.charge(3 * count)?;
        if poles.iter().flatten().any(|x| !x.is_finite()) {
            return Err(Error::Invalid);
        }
        let source = self.source.source();
        let lead = p(source.spec().pitch_mm) * p(f64::from(source.spec().starts));
        w.charge(256)?;
        let controls = self.source.homogeneous_controls().map(|row| row.map(hull));
        let delta = |row: usize, col: usize| -> [I; 3] {
            let expected = ideal(controls[row], self.angular.poles[col], lead);
            std::array::from_fn(|k| {
                p(WEIGHTS[row])
                    * I::new(
                        poles[row * count + col][k].next_down(),
                        poles[row * count + col][k].next_up(),
                    )
                    - expected[k]
            })
        };
        let mut e = 0_f64;
        let mut eu = 0_f64;
        let mut ev = 0_f64;
        let spans = self.angular.knots.len() - 1;
        for row in 0..3 {
            for col in 0..count {
                w.charge(256)?;
                let d = delta(row, col);
                e = e.max(l1(d).hi);
                if row < 2 {
                    let next = delta(row + 1, col);
                    eu = eu.max((p(2.) * l1(std::array::from_fn(|k| next[k] - d[k]))).hi);
                }
                if col + 1 < count {
                    let next = delta(row, col + 1);
                    ev = ev.max(
                        (p((11 * spans) as f64) * l1(std::array::from_fn(|k| next[k] - d[k]))).hi,
                    );
                }
            }
        }
        w.charge(256)?;
        let profile = self.source.whole_bounds();
        let radial = p(source.maximum_radius_mm());
        let dr = p(hull(profile.first_derivative_mm[1]).magnitude());
        let dq = hull(profile.first_derivative_mm[0]);
        // Positive Bernstein weights give W>=3/4 and |Wu|<=1. These
        // quotient bounds include original algebraic coefficient uncertainty,
        // actual multiplication/serialization rounding, and BOTH derivatives.
        let ep = p(e) / p(0.75) + (radial + lead) * p(self.angular.position_error);
        let du = p(eu) / p(0.75) + p(e) / p(0.75).square() + dr * p(self.angular.position_error);
        let dv = p(ev) / p(0.75) + (radial + lead) * p(self.angular.derivative_error);
        let tau = I::new(
            std::f64::consts::TAU.next_down(),
            std::f64::consts::TAU.next_up(),
        );
        let span = p(self.turns[1]) - p(self.turns[0]);
        let source_u = dr + p(dq.magnitude());
        let source_v = (tau * radial + lead) * span;
        let cross_error = du * source_v + dv * source_u + du * dv;
        let normal_lower = tau * hull(source.minimum_radius_mm()) * p(dq.lo) * span;
        if ep.hi > options.position_error_mm
            || du.hi > options.derivative_error_mm
            || dv.hi > options.derivative_error_mm
            || normal_lower.lo <= 0.
            || (p(2.) * cross_error).hi >= normal_lower.lo
        {
            return Err(Error::Unresolved);
        }
        // |delta N|/|N|<1/2: positive dot is retained and angle<=2 ratio.
        let angle = p(2.) * cross_error / p(normal_lower.lo);
        Ok(SupportBounds {
            position_error_mm: ep.hi,
            derivative_error_mm: [du.hi, dv.hi],
            normal_angle_rad: angle.hi,
        })
    }
    pub fn source(&self) -> &AuthoredMetricThread {
        self.source.source()
    }
    pub fn physical_turns(&self) -> [f64; 2] {
        self.turns
    }
    pub fn poles_mm(&self) -> &[[f64; 3]] {
        &self.poles
    }
    pub fn angular_knots(&self) -> &[f64] {
        &self.angular.knots
    }
    pub fn weights(&self) -> [f64; 3] {
        WEIGHTS
    }
    pub fn bounds(&self) -> SupportBounds {
        self.bound
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_bytes(&self) -> usize {
        self.bytes
    }
    /// Same-source check of independently supplied literal CONTROL DATA.
    /// Does not decode a file or grant material/topology authority.
    pub fn certify_controls(
        &self,
        poles: &[[f64; 3]],
        options: CurveOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<SupportBounds, Error> {
        let mut w = Work::new(options.work, interrupt)?;
        w.charge(1)?;
        options.validate()?;
        if self.angular.knots.len() - 1 > options.maximum_spans_per_curve {
            return Err(Error::Limit);
        }
        let mut memory = Memory::new(self.bytes, options.bytes)?;
        memory.add(SCRATCH)?;
        if poles.as_ptr() != self.poles.as_ptr() {
            memory.add(std::mem::size_of_val(poles))?;
        }
        let result = self.certify(poles, options, &mut w)?;
        w.charge(1)?;
        Ok(result)
    }
}
#[cfg(test)]
mod tests;
