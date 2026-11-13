//! Authored shared-boundary flat threads, distinct from legacy literal pieces.
//!
//! A canonical Spec creates common phase/radius endpoints. Every flank is
//! normalized interpolation of those same endpoints, so a source join owns one
//! expression. An independently evaluated legacy Thread is only an approximation
//! with a derived per-band surface correspondence and never source authority.
use crate::helical_thread::{Spec, Thread, ThreadError};
use crate::tapered_thread::interval::I;
#[cfg(test)]
mod tests;
mod transfer;
pub use crate::authored_internal_treatment::{
    ShellEnd, ThreadBoundary, ThreadRimEdge, ThreadShell,
};
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    Invalid,
    Unsupported,
    Limit,
    Unresolved,
    Cancelled,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SourceOptions {
    pub position_limit_mm: f64,
    pub normal_angle_limit_rad: f64,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Band {
    pub phase_mm: [f64; 2],
    pub radius_mm: [f64; 2],
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BandTransfer {
    pub position_error_mm: f64,
    pub slope_error: f64,
    pub normal_angle_error_rad: f64,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransferBounds {
    pub bands: [BandTransfer; 5],
    pub maximum_position_error_mm: f64,
    pub maximum_normal_angle_error_rad: f64,
    /// Error of the angular parameterization over the entire finite q/z chart.
    pub angular_error_rad: f64,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PointSample {
    pub point_mm: [f64; 3],
    pub error_mm: f64,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AuthoredThread {
    spec: Spec,
    phases: [f64; 6],
    radii: [f64; 6],
    approximation: Thread,
    transfer: TransferBounds,
}
impl AuthoredThread {
    pub fn new(spec: Spec, options: SourceOptions) -> Result<Self, Error> {
        Self::new_with_interrupt(spec, options, None)
    }
    pub fn new_with_interrupt(
        spec: Spec,
        options: SourceOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        checkpoint(interrupt)?;
        if !options.position_limit_mm.is_finite()
            || options.position_limit_mm <= 0.
            || !options.normal_angle_limit_rad.is_finite()
            || options.normal_angle_limit_rad <= 0.
        {
            return Err(Error::Invalid);
        }
        if spec.root_radius_mm.is_some() {
            return Err(Error::Unsupported);
        }
        let approximation = Thread::new(spec).map_err(|e| match e {
            ThreadError::Invalid => Error::Invalid,
            ThreadError::Limit => Error::Limit,
        })?;
        // These are authored finite shared expressions from the canonical
        // dimensions, not a recipe inferred from a retained Thread object.
        let half = spec.crest_width_mm * 0.5;
        let upper = half + spec.depth_mm * spec.upper_flank_deg.to_radians().tan();
        let crest = spec.pitch_mm - half;
        let lower = crest - spec.depth_mm * spec.lower_flank_deg.to_radians().tan();
        let phases = [0., half, upper, lower, crest, spec.pitch_mm];
        let minor = spec.major_radius_mm - spec.depth_mm;
        let radii = [
            spec.major_radius_mm,
            spec.major_radius_mm,
            minor,
            minor,
            spec.major_radius_mm,
            spec.major_radius_mm,
        ];
        // The same finite canonical phase expressions have just passed
        // Thread admission; transfer::certify independently requires every
        // normalized source width and literal radius enclosure to be positive.
        let transfer = transfer::certify(spec, phases, radii, &approximation, options, interrupt)?;
        checkpoint(interrupt)?;
        Ok(Self {
            spec,
            phases,
            radii,
            approximation,
            transfer,
        })
    }
    pub fn spec(&self) -> Spec {
        self.spec
    }
    pub fn bands(&self) -> [Band; 5] {
        std::array::from_fn(|i| Band {
            phase_mm: [self.phases[i], self.phases[i + 1]],
            radius_mm: [self.radii[i], self.radii[i + 1]],
        })
    }
    pub fn transfer_bounds(&self) -> TransferBounds {
        self.transfer
    }
    /// This evaluated object is not the shared-boundary manufacturing source.
    pub fn legacy_approximation(&self) -> (&Thread, TransferBounds) {
        (&self.approximation, self.transfer)
    }
    pub fn minimum_radius_mm(&self) -> f64 {
        self.radii[2]
    }
    pub fn maximum_radius_mm(&self) -> f64 {
        self.spec.major_radius_mm
    }
    pub fn full_length_mm(&self) -> f64 {
        self.spec.length_mm
    }
    /// Nominal point plus an outward Euclidean evaluation enclosure.
    /// The authoritative chart closes by true 2pi after exactly one full lead.
    pub fn point_in_band(
        &self,
        band: usize,
        start: u32,
        u: f64,
        z: f64,
    ) -> Result<PointSample, Error> {
        if band >= 5
            || start >= self.spec.starts
            || !u.is_finite()
            || !(0. ..=1.).contains(&u)
            || !z.is_finite()
            || !(0. ..=self.spec.length_mm).contains(&z)
        {
            return Err(Error::Invalid);
        }
        let b = self.bands()[band];
        let q = I::point(b.phase_mm[0])
            + (I::point(b.phase_mm[1]) - I::point(b.phase_mm[0])) * I::point(u)
            + I::point(f64::from(start)) * I::point(self.spec.pitch_mm);
        let r = I::point(b.radius_mm[0])
            + (I::point(b.radius_mm[1]) - I::point(b.radius_mm[0])) * I::point(u);
        let angle = (I::point(z) - q) * self.angular_gain();
        let exact = [r * angle.cos(), r * angle.sin(), I::point(z)];
        let [qbase, rn] = if u == 0. {
            [b.phase_mm[0], b.radius_mm[0]]
        } else if u == 1. {
            [b.phase_mm[1], b.radius_mm[1]]
        } else {
            [
                (b.phase_mm[1] - b.phase_mm[0]).mul_add(u, b.phase_mm[0]),
                (b.radius_mm[1] - b.radius_mm[0]).mul_add(u, b.radius_mm[0]),
            ]
        };
        let qn = qbase + f64::from(start) * self.spec.pitch_mm;
        let theta = (if self.spec.left_hand { -1. } else { 1. }) * std::f64::consts::TAU * (z - qn)
            / (self.spec.pitch_mm * f64::from(self.spec.starts));
        let (s, c) = theta.sin_cos();
        let point_mm = [rn * c, rn * s, z];
        let error = point_mm
            .into_iter()
            .zip(exact)
            .fold(I::point(0.), |sum, (p, b)| {
                sum + I::point((I::point(p) - b).magnitude())
            });
        // Validated source dimensions, u in [0,1] and the finite turn bound
        // keep all point/derivative arithmetic finite (true angle < 1k rad).
        // The explicit enclosure remains separate from constructor transfer.
        Ok(PointSample {
            point_mm,
            error_mm: error.hi,
        })
    }
    pub(crate) fn angular_gain(&self) -> I {
        angular_gain(self.spec)
    }
}
fn angular_gain(spec: Spec) -> I {
    let pi = I::new(
        std::f64::consts::PI.next_down(),
        std::f64::consts::PI.next_up(),
    );
    let gain = I::point(2.) * pi / (I::point(spec.pitch_mm) * I::point(f64::from(spec.starts)));
    if spec.left_hand { -gain } else { gain }
}
fn checkpoint(interrupt: Option<&dyn Fn() -> bool>) -> Result<(), Error> {
    if interrupt.is_some_and(|check| check()) {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
