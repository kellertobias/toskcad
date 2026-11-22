//! Internal thread cavities retain the original phase law and add material removal.
//!
//! The entry is the MAX of that law and a decreasing cone through the literal
//! `(minor + c, 0)` and `(minor, c)` endpoints. The terminal cylinder enlarges
//! the cavity over its last authored interval; neither is an external mask.
use crate::helical_thread::Thread;
use multicad_precision::PortableMath as _;
#[cfg(test)]
mod tests;
mod volume;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    Invalid,
    Limit,
    Unresolved,
    Cancelled,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Terminal {
    pub width_mm: f64,
    pub radial_depth_mm: f64,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Intent {
    pub entry_depth_mm: Option<f64>,
    pub terminal: Option<Terminal>,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Enclosure {
    pub lower: f64,
    pub upper: f64,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QueryOptions {
    pub tolerance_mm3: f64,
    pub maximum_cells: usize,
}
impl Default for QueryOptions {
    fn default() -> Self {
        Self {
            tolerance_mm3: 1e-7,
            maximum_cells: 1_000_000,
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Cavity {
    thread: Thread,
    intent: Intent,
    entry_radius: f64,
    entry_slope: Option<[f64; 2]>,
    terminal_station: f64,
    terminal_radius: f64,
    maximum_radius_bound: f64,
}
impl Cavity {
    pub fn new(thread: Thread, intent: Intent) -> Result<Self, Error> {
        let spec = thread.spec();
        let c = intent.entry_depth_mm.unwrap_or(0.);
        let w = intent.terminal.map_or(0., |tail| tail.width_mm);
        let d = intent.terminal.map_or(0., |tail| tail.radial_depth_mm);
        if spec.root_radius_mm.is_some()
            || (intent.entry_depth_mm.is_none() && intent.terminal.is_none())
            || intent
                .entry_depth_mm
                .is_some_and(|v| !v.is_finite() || v < 0.001)
            || intent.terminal.is_some_and(|tail| {
                [tail.width_mm, tail.radial_depth_mm]
                    .into_iter()
                    .any(|v| !v.is_finite() || v < 0.001)
            })
            || c + w + spec.pitch_mm > spec.length_mm
        {
            return Err(Error::Invalid);
        }
        let entry_radius = thread.root_radius_mm() + c;
        let terminal_station = spec.length_mm - w;
        let terminal_radius = spec.major_radius_mm + d;
        if ![entry_radius, terminal_radius]
            .into_iter()
            .all(|v| v.is_finite() && v <= 1e7)
        {
            return Err(Error::Limit);
        }
        // The original endpoint expressions define a finite cone slope.
        // Qualify this immutable source coefficient once for every query.
        let entry_slope = intent.entry_depth_mm.map(|depth| {
            use crate::tapered_thread::interval::I;
            let slope =
                (I::point(entry_radius) - I::point(thread.root_radius_mm())) / I::point(depth);
            [slope.lo, slope.hi]
        });
        // Explicit finite-query invariants: factor pi/p<4000, all source
        // radii/stations<=1e7, affine slope magnitude<58 and finite band hulls.
        // The endpoint-derived cone slope stays inside [.5,2]. These bounds
        // keep all volume products (including 1M-cell accumulation) finite.
        let mut valid = spec.pitch_mm >= 0.001
            && spec.pitch_mm <= 1e7
            && spec.length_mm <= 1e7
            && spec.major_radius_mm <= 1e7
            && terminal_station > c
            && (intent.terminal.is_none() || terminal_station < spec.length_mm)
            && (c == 0. || entry_radius > thread.root_radius_mm())
            && (d == 0. || terminal_radius > spec.major_radius_mm)
            && entry_slope.is_none_or(|[lo, hi]| {
                lo.is_finite() && hi.is_finite() && lo >= 0.5 && hi <= 2. && lo <= hi
            });
        let mut source = Self {
            thread,
            intent,
            entry_radius,
            entry_slope,
            terminal_station,
            terminal_radius,
            maximum_radius_bound: spec
                .major_radius_mm
                .max(if intent.entry_depth_mm.is_some() {
                    entry_radius
                } else {
                    0.
                })
                .max(if intent.terminal.is_some() {
                    terminal_radius
                } else {
                    0.
                }),
        };
        // Validate the actual retained affine bands, rather than assume that
        // a positive nominal root survives all literal source arithmetic.
        let knots = thread.profile_breaks();
        for band in 0..5 {
            let (radius, slope) = source.profile(
                band,
                crate::tapered_thread::interval::I::new(knots[band], knots[band + 1]),
            );
            valid &= radius.valid()
                && radius.lo > 0.
                && radius.hi <= 1e9
                && radius.width() <= 1e9
                && slope.valid()
                && slope.magnitude() <= 58.;
            source.maximum_radius_bound = source.maximum_radius_bound.max(radius.hi);
        }
        if !valid {
            return Err(Error::Unresolved);
        }
        Ok(source)
    }
    pub fn source_thread(&self) -> &Thread {
        &self.thread
    }
    pub fn intent(&self) -> Intent {
        self.intent
    }
    pub fn full_length_mm(&self) -> f64 {
        self.thread.spec().length_mm
    }
    pub fn entry_radius_mm(&self) -> Option<f64> {
        self.intent.entry_depth_mm.map(|_| self.entry_radius)
    }
    pub fn terminal_station_mm(&self) -> Option<f64> {
        self.intent.terminal.map(|_| self.terminal_station)
    }
    pub fn terminal_radius_mm(&self) -> Option<f64> {
        self.intent.terminal.map(|_| self.terminal_radius)
    }
    /// Outward source bound, including the actual affine flank coefficients.
    /// Receiving wall admission must use this bound, not nominal major radius.
    pub fn maximum_radius_bound_mm(&self) -> f64 {
        self.maximum_radius_bound
    }
    /// Nominal evaluation of the retained cone through its literal endpoints.
    pub fn entry_cone_radius_mm(&self, z: f64) -> Option<f64> {
        let c = self.intent.entry_depth_mm?;
        if !z.is_finite() || !(0. ..=c).contains(&z) {
            return None;
        }
        Some(if z == 0. {
            self.entry_radius
        } else if z == c {
            self.thread.root_radius_mm()
        } else {
            ((self.thread.root_radius_mm() - self.entry_radius) / c).mul_add(z, self.entry_radius)
        })
    }
    /// Station of the cone/profile crossing for an unchanged source phase.
    pub fn entry_seam_station_mm(&self, phase: f64) -> Option<f64> {
        let c = self.intent.entry_depth_mm?;
        if !phase.is_finite() {
            return None;
        }
        let r = self.thread.radius_at_phase(phase);
        Some(
            (c * ((self.entry_radius - r) / (self.entry_radius - self.thread.root_radius_mm())))
                .clamp(0., c),
        )
    }
    pub fn radius_mm(&self, theta: f64, z: f64) -> f64 {
        if !theta.is_finite() || !z.is_finite() || !(0. ..=self.full_length_mm()).contains(&z) {
            return f64::NAN;
        }
        if self.intent.terminal.is_some() && z >= self.terminal_station {
            return self.terminal_radius;
        }
        self.entry_cone_radius_mm(z).map_or_else(
            || self.thread.radius(theta, z),
            |cone| self.thread.radius(theta, z).max(cone),
        )
    }
    pub fn contains_void(&self, point: [f64; 3], tolerance: f64) -> bool {
        if !tolerance.is_finite() || tolerance < 0. || !point.into_iter().all(f64::is_finite) {
            return false;
        }
        let [x, y, z] = point;
        if z < -tolerance || z > self.full_length_mm() + tolerance {
            return false;
        }
        x.portable_hypot(y)
            <= self.radius_mm(y.portable_atan2(x), z.clamp(0., self.full_length_mm())) + tolerance
    }
}
