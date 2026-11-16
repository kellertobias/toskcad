//! Retained prismatic Line 8 stock minus an original authored metric tap and
//! the finite cylindrical crown relief. The relief may open through the crown.
use crate::{
    authored_flat_thread::AuthoredThread, linear_sweep::Sweep, tapered_thread::interval::I,
};
mod curves;
pub use curves::PointBounds;
mod mesh;
mod topology;
pub use topology::{Boundary, Coedge, Edge, Face, Vertex};
mod properties;
pub use mesh::Candidate;
pub use properties::FaceMeasure;
#[cfg(test)]
mod tests;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    Invalid,
    Unsupported,
    Limit,
    Unresolved,
    Cancelled,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Relief {
    pub start_mm: f64,
    pub radius_mm: f64,
    pub tool_length_mm: f64,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Classification {
    Inside,
    Outside,
    Undetermined,
}
#[derive(Clone, Debug)]
pub struct Nut {
    stock: Sweep,
    bore_from_start_mm: f64,
    thread: AuthoredThread,
    relief: Relief,
    phase_rotation: [[f64; 2]; 2],
    width: [f64; 4],
    height: [f64; 4],
    determinant: I,
    minimum_thread: I,
    maximum_thread: I,
}
#[derive(Clone, Copy, Debug)]
pub struct Volume {
    bounds: [f64; 2],
    work: usize,
}
impl Volume {
    pub fn bounds_mm3(&self) -> [f64; 2] {
        self.bounds
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
}
fn checkpoint(interrupt: Option<&dyn Fn() -> bool>) -> Result<(), Error> {
    if interrupt.is_some_and(|f| f()) {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn valid(v: I) -> Result<I, Error> {
    if v.lo.is_finite() && v.hi.is_finite() && v.lo <= v.hi {
        Ok(v)
    } else {
        Err(Error::Unresolved)
    }
}
fn point(v: f64) -> I {
    I::point(v)
}
impl Nut {
    /// Stock section coordinates are [Y,Z], extruded along X. In this chart
    /// original sweep X=0 is at x=-bore_from_start_mm. The supplied phase map
    /// is literal source geometry, neither normalized nor fitted to a circle.
    pub fn new(
        stock: Sweep,
        bore_from_start_mm: f64,
        thread: AuthoredThread,
        phase_rotation: [[f64; 2]; 2],
        relief: Relief,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        checkpoint(interrupt)?;
        if !bore_from_start_mm.is_finite()
            || [relief.start_mm, relief.radius_mm, relief.tool_length_mm]
                .iter()
                .any(|v| !v.is_finite() || *v <= 0.)
            || phase_rotation.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(Error::Invalid);
        }
        if thread.spec().starts != 1 {
            return Err(Error::Unsupported);
        }
        if stock.base() != stock.top() || stock.base().len() != 1 || stock.base()[0].len() != 8 {
            return Err(Error::Unsupported);
        }
        if stock.retained_bytes() > 64 * 1024 * 1024 - std::mem::size_of::<Self>() {
            return Err(Error::Limit);
        }
        let ring = &stock.base()[0];
        let width = [ring[1][0], ring[2][0], ring[3][0], ring[4][0]];
        let height = [ring[1][1], ring[2][1], ring[3][1], ring[4][1]];
        let expected = [
            [-width[0], height[0]],
            [width[0], height[0]],
            [width[1], height[1]],
            [width[2], height[2]],
            [width[3], height[3]],
            [-width[3], height[3]],
            [-width[2], height[2]],
            [-width[1], height[1]],
        ];
        if ring.as_slice() != expected
            || height[0] != 0.
            || width.iter().any(|v| *v <= 0.)
            || !width.windows(2).all(|v| v[0] > v[1])
            || !height.windows(2).all(|v| v[0] < v[1])
        {
            return Err(Error::Unsupported);
        }
        // The whole-box half-space classifier below uses convex stock.
        // Prove it from the original retained literal vertices, rather than
        // inferring convexity from a family name or descending widths.
        for i in 0..8 {
            let a = ring[i];
            let b = ring[(i + 1) % 8];
            let c = ring[(i + 2) % 8];
            let turn = (point(b[0]) - point(a[0])) * (point(c[1]) - point(b[1]))
                - (point(b[1]) - point(a[1])) * (point(c[0]) - point(b[0]));
            if turn.lo <= 0. {
                return Err(Error::Unsupported);
            }
        }
        if [
            bore_from_start_mm,
            relief.start_mm,
            relief.radius_mm,
            relief.tool_length_mm,
        ]
        .iter()
        .chain(width.iter())
        .chain(height.iter())
        .any(|v| v.abs() > 1e7)
        {
            return Err(Error::Limit);
        }
        let m = phase_rotation.map(|r| r.map(point));
        let determinant = m[0][0] * m[1][1] - m[0][1] * m[1][0];
        let a = m[0][0] * m[0][0] + m[1][0] * m[1][0];
        let b = m[0][1] * m[0][1] + m[1][1] * m[1][1];
        let c = m[0][0] * m[0][1] + m[1][0] * m[1][1];
        let delta = point((a - point(1.)).magnitude())
            + point((b - point(1.)).magnitude())
            + point(2.) * point(c.magnitude());
        if determinant.lo <= 0. || delta.hi > 1e-6 {
            return Err(Error::Unsupported);
        }
        let minimum_thread =
            point(thread.minimum_radius_mm()) * (point(1.) - point(delta.hi)).sqrt();
        let maximum_thread =
            point(thread.maximum_radius_mm()) * (point(1.) + point(delta.hi)).sqrt();
        // The full original tap/relief union equals a prefix ending at the relief
        // start, then this cylinder: the cylinder strictly contains the entire
        // mapped original tap in the actual overlap. No nominal length is changed.
        if relief.start_mm < thread.spec().pitch_mm
            || relief.start_mm >= thread.full_length_mm()
            || relief.start_mm < height[1]
            || relief.start_mm >= height[2]
            || relief.radius_mm <= maximum_thread.hi
            || relief.radius_mm <= width[3]
            || relief.radius_mm >= width[2]
            || (point(relief.start_mm) + point(relief.tool_length_mm)).lo <= height[3]
        {
            return Err(Error::Unsupported);
        }
        let at_join = point(width[1])
            + (point(width[2]) - point(width[1])) * (point(relief.start_mm) - point(height[1]))
                / (point(height[2]) - point(height[1]));
        if at_join.lo <= maximum_thread.hi
            || bore_from_start_mm <= relief.radius_mm
            || (point(stock.height_mm()) - point(bore_from_start_mm)).lo <= relief.radius_mm
        {
            return Err(Error::Unsupported);
        }
        checkpoint(interrupt)?;
        Ok(Self {
            stock,
            bore_from_start_mm,
            thread,
            relief,
            phase_rotation,
            width,
            height,
            determinant: valid(determinant)?,
            minimum_thread: valid(minimum_thread)?,
            maximum_thread: valid(maximum_thread)?,
        })
    }
    pub fn retained_bytes(&self) -> usize {
        std::mem::size_of::<Self>() + self.stock.retained_bytes() - std::mem::size_of::<Sweep>()
    }
    pub fn stock(&self) -> &Sweep {
        &self.stock
    }
    pub fn thread(&self) -> &AuthoredThread {
        &self.thread
    }
    pub fn relief(&self) -> Relief {
        self.relief
    }
    pub fn phase_rotation(&self) -> [[f64; 2]; 2] {
        self.phase_rotation
    }
    pub fn bore_from_start_mm(&self) -> f64 {
        self.bore_from_start_mm
    }
    /// Euclidean lower separation of the entire finite mapped tap envelope
    /// from the original six sloping side planes and both extrusion caps.
    /// Axial entry/relief openings are intentionally not wall constraints.
    pub fn minimum_thread_wall_lower_bound_mm(&self) -> f64 {
        let r = self.maximum_thread;
        let mut lower = (point(self.bore_from_start_mm) - r)
            .lo
            .min((point(self.stock.height_mm()) - point(self.bore_from_start_mm) - r).lo);
        let y = I::new(-r.hi, r.hi);
        let z = I::new(0., self.relief.start_mm);
        let ring = &self.stock.base()[0];
        for side in [1usize, 2, 3, 5, 6, 7] {
            let a = ring[side];
            let b = ring[(side + 1) % 8];
            let dy = point(b[0]) - point(a[0]);
            let dz = point(b[1]) - point(a[1]);
            let signed = dy * (z - point(a[1])) - dz * (y - point(a[0]));
            let normal = (dy * dy + dz * dz).sqrt();
            lower = lower.min((signed / normal).lo);
        }
        lower.max(0.)
    }
    pub fn bounds_mm(&self) -> [[f64; 2]; 3] {
        let x = point(self.stock.height_mm()) - point(self.bore_from_start_mm);
        [
            [-self.bore_from_start_mm, x.hi],
            [-self.width[0], self.width[0]],
            [0., self.height[3]],
        ]
    }
    pub fn volume_bounds(
        &self,
        maximum_width_mm3: f64,
        maximum_work: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Volume, Error> {
        checkpoint(interrupt)?;
        if !maximum_width_mm3.is_finite() || maximum_width_mm3 <= 0. {
            return Err(Error::Invalid);
        }
        if maximum_work < 2048 {
            return Err(Error::Limit);
        }
        let pi = I::new(
            std::f64::consts::PI.next_down(),
            std::f64::consts::PI.next_up(),
        );
        let mut blank = point(0.);
        for i in 0..3 {
            blank = blank
                + (point(self.width[i]) + point(self.width[i + 1]))
                    * (point(self.height[i + 1]) - point(self.height[i]));
        }
        blank = blank * point(self.stock.height_mm());
        let mut square_radius = point(0.);
        for b in self.thread.bands() {
            let a = point(b.radius_mm[0]);
            let c = point(b.radius_mm[1]);
            square_radius = square_radius
                + (point(b.phase_mm[1]) - point(b.phase_mm[0])) * (a * a + a * c + c * c)
                    / point(3.);
        }
        let tap = pi * square_radius / point(self.thread.spec().pitch_mm)
            * point(self.relief.start_mm)
            * self.determinant;
        let r = point(self.relief.radius_mm);
        let a = point(self.width[3]);
        let d = valid(r * r - a * a)?.sqrt();
        let asin = (a / d).atan();
        // Integral over the actual crown opening: the sloping roof removes two
        // circular wedges. This is a disk/section intersection, never full-disk
        // subtraction above the point where the crown width falls below radius.
        let integral = d * d * d / point(3.)
            - a * (pi * r * r / point(4.) - (a * d + r * r * asin) / point(2.));
        let slope = (point(self.height[3]) - point(self.height[2]))
            / (point(self.width[2]) - point(self.width[3]));
        let relief = pi * r * r * (point(self.height[3]) - point(self.relief.start_mm))
            - point(4.) * slope * integral;
        let volume = valid(blank - tap - relief)?;
        if volume.lo <= 0. || volume.width() > maximum_width_mm3 {
            return Err(Error::Unresolved);
        }
        checkpoint(interrupt)?;
        Ok(Volume {
            bounds: [volume.lo, volume.hi],
            work: 2048,
        })
    }
    pub fn classify_box(
        &self,
        bounds: [[f64; 2]; 3],
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Classification, Error> {
        checkpoint(interrupt)?;
        if bounds
            .iter()
            .any(|v| !v[0].is_finite() || !v[1].is_finite() || v[0] > v[1])
        {
            return Err(Error::Invalid);
        }
        if bounds.iter().flatten().any(|v| v.abs() > 1e7) {
            return Err(Error::Limit);
        }
        let [x, y, z] = bounds.map(|v| I::new(v[0], v[1]));
        let relative = x + point(self.bore_from_start_mm);
        if relative.hi < 0. || relative.lo > self.stock.height_mm() {
            return Ok(Classification::Outside);
        }
        let mut all = relative.lo > 0. && relative.hi < self.stock.height_mm();
        let ring = &self.stock.base()[0];
        for i in 0..8 {
            let a = ring[i];
            let b = ring[(i + 1) % 8];
            let side = (point(b[0]) - point(a[0])) * (z - point(a[1]))
                - (point(b[1]) - point(a[1])) * (y - point(a[0]));
            if side.hi < 0. {
                return Ok(Classification::Outside);
            }
            all &= side.lo > 0.;
        }
        // Interval square handles boxes crossing the axis.
        let square = |v: I| {
            let lower = if v.lo <= 0. && v.hi >= 0. {
                0.
            } else {
                v.lo.abs().min(v.hi.abs()).powi(2).next_down().max(0.)
            };
            I::new(lower, v.magnitude().powi(2).next_up())
        };
        let radius = valid(square(x) + square(y))?.sqrt();
        let (inside_cavity, outside_cavity) = if z.lo >= self.relief.start_mm {
            (
                radius.hi < self.relief.radius_mm,
                radius.lo > self.relief.radius_mm,
            )
        } else if z.hi < self.relief.start_mm {
            (
                radius.hi < self.minimum_thread.lo,
                radius.lo > self.maximum_thread.hi,
            )
        } else {
            (
                radius.hi < self.minimum_thread.lo.min(self.relief.radius_mm),
                radius.lo > self.maximum_thread.hi.max(self.relief.radius_mm),
            )
        };
        if inside_cavity {
            return Ok(Classification::Outside);
        }
        if all && outside_cavity {
            Ok(Classification::Inside)
        } else {
            Ok(Classification::Undetermined)
        }
    }
}
