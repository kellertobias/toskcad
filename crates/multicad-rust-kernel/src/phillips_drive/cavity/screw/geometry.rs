//! Queries retain the original thread/head union and actual finite cutter.
use super::*;
impl PhillipsScrew {
    /// Support of the occupied outer union. Host qualification proves a
    /// positive radial wall and floor, so the finite recess removes none of
    /// its directional extrema, including the remaining top annular land.
    pub fn support(&self, direction: [f64; 3]) -> f64 {
        if !direction.into_iter().all(f64::is_finite) {
            return f64::NAN;
        }
        let scale = direction.into_iter().map(f64::abs).fold(0., f64::max);
        if scale == 0. {
            return 0.;
        }
        let direction = direction.map(|x| x / scale);
        scale
            * self
                .thread
                .support(direction)
                .max(self.profile.support(direction) - direction[2] * self.spec.shank_length_mm)
    }
    /// Tolerant outer-source membership with conservative finite-void refusal.
    /// A point whose actual translated cutter classification is unresolved is
    /// never admitted as occupied material. This is not a distance certificate.
    pub fn contains(&self, point: [f64; 3], tolerance: f64) -> bool {
        if !point.into_iter().all(f64::is_finite) || !tolerance.is_finite() || tolerance < 0. {
            return false;
        }
        let [x, y, z] = point;
        let local_z = z + self.spec.shank_length_mm;
        let [top, end] = [self.profile.top()[1], self.profile.end()[1]];
        let in_profile = local_z >= top - tolerance
            && local_z <= end + tolerance
            && x.hypot(y)
                <= self
                    .profile
                    .radius_at(local_z.clamp(top, end))
                    .expect("qualified original profile domain")
                    + tolerance;
        if !(self.thread.contains(point, tolerance) || in_profile) {
            return false;
        }
        let depth = I::point(z) - I::point(self.head_top_z());
        // Keep the subtraction enclosure; rounding z-head_top to one point can
        // incorrectly resolve an actual cavity boundary into solid material.
        classify::enclosure(self.cavity.source(), [I::point(x), I::point(y), depth])
            == Ok(Classification::Outside)
    }
}
#[cfg(test)]
mod tests;
