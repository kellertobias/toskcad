//! Actual helical screws with a retained, variable-section head and hex drive.
use crate::{
    head_profile::{Curve, Head, Profile},
    headed_screw::{hex_area, hex_radius},
    helical_thread::{Thread, ThreadError},
};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Spec {
    pub head: Head,
    pub shank_length_mm: f64,
    pub socket_across_flats_mm: f64,
    pub socket_depth_mm: f64,
    pub minimum_wall_mm: f64,
}

/// Actual drive/junction data; authored head metadata is absent on a physical
/// import. Shape decisions always use the retained profile, never this hint.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SourceSpec {
    pub head: Option<Head>,
    pub shank_length_mm: f64,
    pub socket_across_flats_mm: f64,
    pub socket_depth_mm: f64,
    pub minimum_wall_mm: f64,
}

impl From<Spec> for SourceSpec {
    fn from(spec: Spec) -> Self {
        Self {
            head: Some(spec.head),
            shank_length_mm: spec.shank_length_mm,
            socket_across_flats_mm: spec.socket_across_flats_mm,
            socket_depth_mm: spec.socket_depth_mm,
            minimum_wall_mm: spec.minimum_wall_mm,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Face {
    ThreadShoulder,
    ThreadTip,
    ThreadWall(u32),
    HeadTop,
    ProfileWall(u32),
    SocketFloor,
    SocketWall(u32),
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Screw {
    thread: Thread,
    spec: SourceSpec,
    profile: Profile,
    visible_start_mm: f64,
}

impl Screw {
    pub fn new(thread: Thread, spec: Spec) -> Result<Self, ThreadError> {
        let profile = Profile::new(
            spec.head,
            thread.spec().major_radius_mm,
            spec.shank_length_mm,
        )?;
        Self::admit_source(thread, profile, spec.into(), false)
    }

    /// Physical meridian/drive admission without guessing an authoring recipe.
    pub fn from_source(
        thread: Thread,
        profile: Profile,
        spec: SourceSpec,
    ) -> Result<Self, ThreadError> {
        if spec.head.is_some() {
            return Err(ThreadError::Invalid);
        }
        Self::admit_source(thread, profile, spec, false)
    }

    /// Opt in to the actual head/thread union; the original thread law and
    /// authored length remain unchanged, including the covered prefix.
    pub fn new_with_head_thread_union(thread: Thread, spec: Spec) -> Result<Self, ThreadError> {
        let profile = Profile::new_for_thread_union(
            spec.head,
            thread.spec().major_radius_mm,
            spec.shank_length_mm,
        )?;
        Self::admit_source(thread, profile, spec.into(), true)
    }

    pub fn from_source_with_union(
        thread: Thread,
        profile: Profile,
        spec: SourceSpec,
    ) -> Result<Self, ThreadError> {
        if spec.head.is_some() {
            return Err(ThreadError::Invalid);
        }
        Self::admit_source(thread, profile, spec, true)
    }

    fn admit_source(
        thread: Thread,
        profile: Profile,
        spec: SourceSpec,
        union: bool,
    ) -> Result<Self, ThreadError> {
        let major = thread.spec().major_radius_mm;
        let visible_start_mm = profile.end()[1] - spec.shank_length_mm;
        if !spec.shank_length_mm.is_finite()
            || spec.shank_length_mm < 0.
            || profile.end()[0] != major
            || visible_start_mm < 0.
            || (!union && visible_start_mm != 0.)
            || thread.spec().length_mm - visible_start_mm <= 4. * profile.rounding_error_mm()
        {
            return Err(ThreadError::Invalid);
        }
        if visible_start_mm > 0.
            && profile.segments().any(|segment| {
                if segment.end[1] <= spec.shank_length_mm {
                    return false;
                }
                match segment.curve {
                    Curve::Line => segment.start[0].min(segment.end[0]) < major,
                    Curve::Arc {
                        center, radius_mm, ..
                    } => {
                        // Error-free subtraction: comparing the rounded result to a
                        // distinct f64 is decisive; equality needs its exact residual.
                        let d = center[0] - radius_mm;
                        let b = d - center[0];
                        let residual = (center[0] - (d - b)) - (radius_mm + b);
                        d < major || (d == major && residual < 0.)
                    }
                }
            })
        {
            return Err(ThreadError::Invalid);
        }
        if spec.socket_across_flats_mm <= 2. * profile.rounding_error_mm()
            || spec.socket_depth_mm <= 2. * profile.rounding_error_mm()
            || !profile.admits_socket(
                spec.socket_across_flats_mm / 3_f64.sqrt(),
                spec.socket_depth_mm,
                spec.minimum_wall_mm,
            )
            || profile.top()[1] + spec.socket_depth_mm >= 0.
        {
            return Err(ThreadError::Invalid);
        }
        Ok(Self {
            thread,
            spec,
            profile,
            visible_start_mm,
        })
    }

    pub fn thread(&self) -> &Thread {
        &self.thread
    }
    /// Original source coordinates, retaining helical phase at the junction.
    pub fn visible_thread_interval_mm(&self) -> [f64; 2] {
        [self.visible_start_mm, self.thread.spec().length_mm]
    }
    pub fn has_head_thread_overlap(&self) -> bool {
        self.visible_start_mm > 0.
    }
    pub fn visible_thread_volume_mm3(&self) -> f64 {
        self.thread.volume_mm3()
            * ((self.thread.spec().length_mm - self.visible_start_mm)
                / self.thread.spec().length_mm)
    }
    pub fn spec(&self) -> SourceSpec {
        self.spec
    }
    pub fn profile(&self) -> &Profile {
        &self.profile
    }
    pub fn head_top_z(&self) -> f64 {
        self.profile.top()[1] - self.spec.shank_length_mm
    }
    pub fn socket_radius(&self, angle: f64) -> f64 {
        hex_radius(self.spec.socket_across_flats_mm, angle)
    }

    pub fn faces(&self) -> Vec<Face> {
        let mut result = vec![Face::ThreadShoulder, Face::ThreadTip];
        result.extend((0..4 * self.thread.spec().starts).map(Face::ThreadWall));
        result.push(Face::HeadTop);
        result.extend((0..self.profile.segments().count() as u32).map(Face::ProfileWall));
        result.push(Face::SocketFloor);
        result.extend((0..6).map(Face::SocketWall));
        result
    }

    pub fn volume_mm3(&self) -> f64 {
        self.visible_thread_volume_mm3() + self.profile.volume_mm3()
            - hex_area(self.spec.socket_across_flats_mm) * self.spec.socket_depth_mm
    }

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

    pub fn contains(&self, point: [f64; 3], tolerance: f64) -> bool {
        if !point.into_iter().all(f64::is_finite) || !tolerance.is_finite() || tolerance < 0. {
            return false;
        }
        if self.thread.contains(point, tolerance) {
            return true;
        }
        let [x, y, z] = point;
        let local_z = z + self.spec.shank_length_mm;
        let [top, end] = [self.profile.top()[1], self.profile.end()[1]];
        if local_z < top - tolerance || local_z > end + tolerance {
            return false;
        }
        let radius = x.hypot(y);
        radius
            <= self
                .profile
                .radius_at(local_z.clamp(top, end))
                .expect("profile domain")
                + tolerance
            && (z >= self.head_top_z() + self.spec.socket_depth_mm - tolerance
                || radius >= self.socket_radius(y.atan2(x)) - tolerance)
    }
}

#[cfg(test)]
mod tests;
