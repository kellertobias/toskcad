//! Directed analytic source mass, distinct from literal-file mass decoding.
use super::*;
use crate::helical_thread::RootCircle;
use std::f64::consts::PI;
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SourceVolumeBounds {
    volume_mm3: [f64; 2],
    profile_precision_error_mm3: f64,
    cavity_source_width_mm3: f64,
}
impl SourceVolumeBounds {
    pub fn volume_mm3(&self) -> [f64; 2] {
        self.volume_mm3
    }
    pub fn profile_precision_error_mm3(&self) -> f64 {
        self.profile_precision_error_mm3
    }
    pub fn cavity_source_width_mm3(&self) -> f64 {
        self.cavity_source_width_mm3
    }
}
fn cube(x: I) -> I {
    x * x * x
}
fn root_primitive(root: RootCircle, s: f64) -> Option<I> {
    let x = I::point(s) - I::point(root.phase);
    let r = I::point(root.radius);
    let c = I::point(root.center_radius);
    let height = (r.square() - x.square()).sqrt();
    if !height.valid() || height.lo <= 0. {
        return None;
    }
    // The admitted root trim stays strictly inside its circle's vertical
    // tangencies. asin(x/r)=atan(x/sqrt(r²-x²)) on this actual finite domain.
    let angle = (x / height).atan();
    Some(
        (c.square() + r.square()) * x
            - cube(x) / I::point(3.)
            - c * (x * height + r.square() * angle),
    )
}
fn thread_volume(thread: &Thread, start: f64) -> Option<I> {
    let breaks = thread.profile_breaks();
    let spec = thread.spec();
    let major = I::point(spec.major_radius_mm);
    let mut integral = I::point(0.);
    for (band, w) in breaks.windows(2).enumerate() {
        let [a, b] = [w[0], w[1]];
        let value = if band == 2
            && let Some(root) = thread.root_circle()
        {
            root_primitive(root, b)? - root_primitive(root, a)?
        } else {
            let radius = |phase| match band {
                1 => {
                    major
                        + (I::point(phase) - I::point(breaks[1]))
                            * I::point(thread.slope_at_phase(phase, 1))
                }
                3 => {
                    major
                        + (I::point(phase) - I::point(breaks[4]))
                            * I::point(thread.slope_at_phase(phase, 3))
                }
                2 => major - I::point(spec.depth_mm),
                _ => major,
            };
            let [ra, rb] = [radius(a), radius(b)];
            (I::point(b) - I::point(a)) * (ra.square() + ra * rb + rb.square()) / I::point(3.)
        };
        integral = integral + value;
    }
    let pi = I::new(PI.next_down(), PI.next_up());
    let result =
        pi * (I::point(spec.length_mm) - I::point(start)) / I::point(spec.pitch_mm) * integral;
    (result.valid() && result.lo > 0.).then_some(result)
}
impl PhillipsScrew {
    /// Enclosure of this retained authored SOURCE's material volume. Head
    /// precision, numerical integration and actual cavity-source bounds remain
    /// explicit. It does not measure an independently decoded STEP file.
    pub fn source_volume_bounds(&self) -> Result<SourceVolumeBounds, Error> {
        let head = crate::head_profile::source_volume::volume_bounds(&self.profile)
            .ok_or(Error::Unresolved)?;
        let precision = crate::head_profile::source_volume::precision_envelope(&self.profile)
            .ok_or(Error::Unresolved)?;
        let thread = thread_volume(&self.thread, self.visible_start_mm).ok_or(Error::Unresolved)?;
        let cavity = self.cavity.volume.volume_mm3;
        let cavity = I::new(cavity[0], cavity[1]);
        let volume = head + thread - cavity;
        if !volume.valid() || volume.lo <= 0. {
            return Err(Error::Unresolved);
        }
        Ok(SourceVolumeBounds {
            volume_mm3: [volume.lo, volume.hi],
            profile_precision_error_mm3: precision.volume.hi,
            cavity_source_width_mm3: cavity.width(),
        })
    }
}

#[cfg(test)]
mod tests;
