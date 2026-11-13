//! Actual MAX cavities over an explicitly authored shared-boundary thread.
//! The legacy literal Cavity is a different source and is never admitted here.
use crate::authored_flat_thread::AuthoredThread;
use crate::tapered_thread::interval::I;
pub use crate::thread_internal_treatment::{Enclosure, Error, Intent, QueryOptions, Terminal};
mod material;
mod mesh;
pub(crate) mod properties;
pub use properties::{FaceKey, FaceProperties, FacePropertyOptions, ThreadBandRole};
#[cfg(test)]
mod tests;
mod volume;
pub use material::Material;
pub(crate) use mesh::metric_shell;
pub use mesh::{Candidate, ShellEnd, ThreadRimEdge, ThreadShell, tessellate_candidate};
pub use mesh::{MetricRimEdge, MetricShell};
#[derive(Clone, Copy, Debug, PartialEq)]
struct Entry {
    radius: f64,
    depth: f64,
    slope: [f64; 2],
}
#[derive(Clone, Copy, Debug, PartialEq)]
struct Tail {
    radius: f64,
    station: f64,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Cavity {
    thread: AuthoredThread,
    intent: Intent,
    entry: Option<Entry>,
    tail: Option<Tail>,
}
impl Cavity {
    pub fn new(thread: AuthoredThread, intent: Intent) -> Result<Self, Error> {
        let spec = thread.spec();
        let c = intent.entry_depth_mm.unwrap_or(0.);
        let w = intent.terminal.map_or(0., |t| t.width_mm);
        if (intent.entry_depth_mm.is_none() && intent.terminal.is_none())
            || intent
                .entry_depth_mm
                .is_some_and(|v| !v.is_finite() || v < 0.001)
            || intent.terminal.is_some_and(|t| {
                [t.width_mm, t.radial_depth_mm]
                    .into_iter()
                    .any(|v| !v.is_finite() || v < 0.001)
            })
            || c + w + spec.pitch_mm > spec.length_mm
        {
            return Err(Error::Invalid);
        }
        let entry = intent.entry_depth_mm.map(|depth| {
            let radius = thread.minimum_radius_mm() + depth;
            let slope = (I::point(radius) - I::point(thread.minimum_radius_mm())) / I::point(depth);
            Entry {
                radius,
                depth,
                slope: [slope.lo, slope.hi],
            }
        });
        let tail = intent.terminal.map(|t| Tail {
            radius: spec.major_radius_mm + t.radial_depth_mm,
            station: spec.length_mm - t.width_mm,
        });
        if entry.is_some_and(|e| !e.radius.is_finite() || e.radius > 1e7)
            || tail.is_some_and(|t| !t.radius.is_finite() || t.radius > 1e7)
        {
            return Err(Error::Limit);
        }
        // Enabled dimensions are >=.001 while every admitted radius/station
        // is <=1e7 (ulp<2e-9). Endpoint addition/subtraction therefore cannot
        // collapse an enabled radius or interval, and the endpoint-derived
        // slope stays inside [.5,2]. The retained central pitch is >=.001.
        Ok(Self {
            thread,
            intent,
            entry,
            tail,
        })
    }
    pub fn source_thread(&self) -> &AuthoredThread {
        &self.thread
    }
    pub fn intent(&self) -> Intent {
        self.intent
    }
    pub fn full_length_mm(&self) -> f64 {
        self.thread.full_length_mm()
    }
    pub fn entry_radius_mm(&self) -> Option<f64> {
        self.entry.map(|e| e.radius)
    }
    pub fn terminal_radius_mm(&self) -> Option<f64> {
        self.tail.map(|t| t.radius)
    }
    pub fn terminal_station_mm(&self) -> Option<f64> {
        self.tail.map(|t| t.station)
    }
    pub fn maximum_radius_mm(&self) -> f64 {
        self.thread
            .maximum_radius_mm()
            .max(self.entry_radius_mm().unwrap_or(0.))
            .max(self.terminal_radius_mm().unwrap_or(0.))
    }
    /// Nominal source evaluation; shared band endpoints remain exact data.
    pub fn radius_mm(&self, theta: f64, z: f64) -> f64 {
        if !theta.is_finite() || !z.is_finite() || !(0. ..=self.full_length_mm()).contains(&z) {
            return f64::NAN;
        }
        if let Some(tail) = self.tail
            && z >= tail.station
        {
            return tail.radius;
        }
        let spec = self.thread.spec();
        let rate = (if spec.left_hand { -1. } else { 1. }) * spec.pitch_mm * f64::from(spec.starts)
            / std::f64::consts::TAU;
        let phase = (z - rate * theta).rem_euclid(spec.pitch_mm);
        let bands = self.thread.bands();
        let band = bands
            .iter()
            .find(|b| phase <= b.phase_mm[1])
            .unwrap_or(&bands[4]);
        let u = ((phase - band.phase_mm[0]) / (band.phase_mm[1] - band.phase_mm[0])).clamp(0., 1.);
        let r = (band.radius_mm[1] - band.radius_mm[0]).mul_add(u, band.radius_mm[0]);
        self.entry_cone_radius_mm(z).map_or(r, |cone| r.max(cone))
    }
    pub fn entry_cone_radius_mm(&self, z: f64) -> Option<f64> {
        let e = self.entry?;
        if !z.is_finite() || !(0. ..=e.depth).contains(&z) {
            return None;
        }
        Some(if z == 0. {
            e.radius
        } else if z == e.depth {
            self.thread.minimum_radius_mm()
        } else {
            ((self.thread.minimum_radius_mm() - e.radius) / e.depth).mul_add(z, e.radius)
        })
    }
    /// Semantic seam in the source band's normalized coordinate, not a guess
    /// from a retained Thread or a coordinate-nearest edge correspondence.
    pub fn entry_seam_station_in_band(&self, band: usize, u: f64) -> Result<Option<f64>, Error> {
        if band >= 5 || !u.is_finite() || !(0. ..=1.).contains(&u) {
            return Err(Error::Invalid);
        }
        let Some(e) = self.entry else {
            return Ok(None);
        };
        let b = self.thread.bands()[band];
        let r = if u == 0. {
            b.radius_mm[0]
        } else if u == 1. {
            b.radius_mm[1]
        } else {
            (b.radius_mm[1] - b.radius_mm[0]).mul_add(u, b.radius_mm[0])
        };
        Ok(Some(if r == self.thread.minimum_radius_mm() {
            e.depth
        } else {
            (e.depth * ((e.radius - r) / (e.radius - self.thread.minimum_radius_mm())))
                .clamp(0., e.depth)
        }))
    }
    pub fn contains_void(&self, point: [f64; 3], tolerance: f64) -> bool {
        if !tolerance.is_finite() || tolerance < 0. || !point.into_iter().all(f64::is_finite) {
            return false;
        }
        let [x, y, z] = point;
        if z < -tolerance || z > self.full_length_mm() + tolerance {
            return false;
        }
        x.hypot(y) <= self.radius_mm(y.atan2(x), z.clamp(0., self.full_length_mm())) + tolerance
    }
}

pub use mesh::boundary::{
    Boundary, BoundaryCell, BoundaryCellSource, BoundaryCoedge, BoundaryConstraint, BoundaryCurve,
    BoundaryOptions, BoundaryPoint, BoundaryPointKey, BoundaryPointRole, BoundaryPolyline,
    ExactBoundaryData, ExactBoundaryDomain, ExactBoundaryRational, ThreadBoundary,
};
