//! A finite circular entry seat over the unchanged original straight thread.
use crate::helical_thread::{Thread, ThreadError};
#[cfg(test)]
mod tests;
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Counterbore {
    thread: Thread,
    radius: f64,
    depth: f64,
}
impl Counterbore {
    pub fn new(thread: Thread, radius_mm: f64, depth_mm: f64) -> Result<Self, ThreadError> {
        if !radius_mm.is_finite()
            || !depth_mm.is_finite()
            || radius_mm <= thread.spec().major_radius_mm
            || depth_mm <= 0.
            || depth_mm >= thread.spec().length_mm
        {
            return Err(ThreadError::Invalid);
        }
        if radius_mm > 1e7 {
            return Err(ThreadError::Limit);
        }
        Ok(Self {
            thread,
            radius: radius_mm,
            depth: depth_mm,
        })
    }
    pub fn thread(&self) -> &Thread {
        &self.thread
    }
    pub fn radius_mm(&self) -> f64 {
        self.radius
    }
    pub fn depth_mm(&self) -> f64 {
        self.depth
    }
    pub fn visible_thread_interval_mm(&self) -> [f64; 2] {
        [self.depth, self.thread.spec().length_mm]
    }
    pub fn cylinder_face(&self) -> u32 {
        2 + 4 * self.thread.spec().starts
    }
    pub fn shoulder_face(&self) -> u32 {
        self.cylinder_face() + 1
    }
    /// Thread floor/wall ordinals are unchanged; the stock mouth is face zero.
    pub fn cavity_face_count(&self) -> u32 {
        self.shoulder_face() + 1
    }
    pub fn removed_volume_mm3(&self) -> f64 {
        let length = self.thread.spec().length_mm;
        std::f64::consts::PI * self.radius * self.radius * self.depth
            + self.thread.volume_mm3() / length * (length - self.depth)
    }
    pub fn contains_void(&self, point: [f64; 3], tolerance: f64) -> bool {
        if !tolerance.is_finite() || tolerance < 0. || !point.into_iter().all(f64::is_finite) {
            return false;
        }
        let [x, y, z] = point;
        if z < -tolerance || z > self.thread.spec().length_mm + tolerance {
            return false;
        }
        let radius = if z <= self.depth {
            self.radius
        } else {
            self.thread.radius(
                y.atan2(x),
                z.clamp(self.depth, self.thread.spec().length_mm),
            )
        };
        x.hypot(y) <= radius + tolerance
    }
}
