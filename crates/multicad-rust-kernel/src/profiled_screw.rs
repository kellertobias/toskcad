//! Actual helical screws with a retained head and typed physical recess.
use crate::{
    head_profile::{Curve, Head, Profile},
    helical_thread::{Thread, ThreadError},
    hex_drive::HexDrive,
    screw_recess::Recess,
};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Spec {
    pub head: Head,
    pub shank_length_mm: f64,
    pub socket_across_flats_mm: f64,
    pub socket_depth_mm: f64,
    pub minimum_wall_mm: f64,
}

/// Junction data; authored head metadata is absent on a physical import.
/// Outer shape uses the profile, the cavity uses the typed `recess()`;
/// the optional head hint is never geometry authority.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SourceSpec {
    pub head: Option<Head>,
    pub shank_length_mm: f64,
    pub minimum_wall_mm: f64,
}

impl From<Spec> for SourceSpec {
    fn from(spec: Spec) -> Self {
        Self {
            head: Some(spec.head),
            shank_length_mm: spec.shank_length_mm,
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

#[derive(Clone, Debug, PartialEq)]
pub struct Screw {
    thread: Thread,
    spec: SourceSpec,
    profile: Profile,
    visible_start_mm: f64,
    recess: Recess,
}

impl Screw {
    pub fn new(thread: Thread, spec: Spec) -> Result<Self, ThreadError> {
        Self::new_with_shank_radius(thread, thread.spec().major_radius_mm, spec)
    }

    /// The nominal plain shank may exceed the fitted thread's major radius.
    pub fn new_with_shank_radius(
        thread: Thread,
        shank_radius_mm: f64,
        spec: Spec,
    ) -> Result<Self, ThreadError> {
        let profile = Profile::new(spec.head, shank_radius_mm, spec.shank_length_mm)?;
        Self::admit_source(
            thread,
            profile,
            spec.into(),
            false,
            flat_hex(&profile, spec)?,
        )
    }

    /// Physical geometry admission; the recess retains its actual floor and contour.
    pub fn from_source(
        thread: Thread,
        profile: Profile,
        spec: SourceSpec,
        recess: Recess,
    ) -> Result<Self, ThreadError> {
        if spec.head.is_some() {
            return Err(ThreadError::Invalid);
        }
        Self::admit_source(thread, profile, spec, false, recess)
    }

    /// Opt in to the actual head/thread union; the original thread law and
    /// authored length remain unchanged, including the covered prefix.
    pub fn new_with_head_thread_union(thread: Thread, spec: Spec) -> Result<Self, ThreadError> {
        Self::new_with_shank_thread_union(thread, thread.spec().major_radius_mm, spec)
    }

    /// Retain the true shank and the original fitted thread beneath the head.
    pub fn new_with_shank_thread_union(
        thread: Thread,
        shank_radius_mm: f64,
        spec: Spec,
    ) -> Result<Self, ThreadError> {
        let profile =
            Profile::new_for_thread_union(spec.head, shank_radius_mm, spec.shank_length_mm)?;
        Self::admit_source(
            thread,
            profile,
            spec.into(),
            true,
            flat_hex(&profile, spec)?,
        )
    }

    pub fn from_source_with_union(
        thread: Thread,
        profile: Profile,
        spec: SourceSpec,
        recess: Recess,
    ) -> Result<Self, ThreadError> {
        if spec.head.is_some() {
            return Err(ThreadError::Invalid);
        }
        Self::admit_source(thread, profile, spec, true, recess)
    }

    pub fn new_with_recess(
        thread: Thread,
        shank_radius_mm: f64,
        head: Head,
        shank_length_mm: f64,
        minimum_wall_mm: f64,
        recess: Recess,
    ) -> Result<Self, ThreadError> {
        let profile = Profile::new_for_thread_union(head, shank_radius_mm, shank_length_mm)?;
        Self::admit_source(
            thread,
            profile,
            SourceSpec {
                head: Some(head),
                shank_length_mm,
                minimum_wall_mm,
            },
            true,
            recess,
        )
    }

    fn admit_source(
        thread: Thread,
        profile: Profile,
        spec: SourceSpec,
        union: bool,
        recess: Recess,
    ) -> Result<Self, ThreadError> {
        let major = thread.spec().major_radius_mm;
        let visible_start_mm = profile.end()[1] - spec.shank_length_mm;
        if !spec.shank_length_mm.is_finite()
            || spec.shank_length_mm < 0.
            || profile.end()[0] < major
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
        let reserve = profile.rounding_error_mm().max(recess.source_error_mm());
        if recess.maximum_radius_bound_mm() <= 2. * reserve
            || recess.maximum_depth_mm() <= 2. * reserve
            || !profile.admits_socket(
                recess.maximum_radius_bound_mm() + recess.source_error_mm(),
                recess.maximum_depth_mm(),
                spec.minimum_wall_mm,
            )
            || profile.top()[1] + recess.maximum_depth_mm() >= 0.
        {
            return Err(ThreadError::Invalid);
        }
        Ok(Self {
            thread,
            spec,
            profile,
            visible_start_mm,
            recess,
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
    pub fn socket_radius(&self, angle: f64) -> Option<f64> {
        self.recess.radius_at_angle_mm(angle)
    }
    pub fn recess(&self) -> &Recess {
        &self.recess
    }
    pub fn source_error_mm(&self) -> f64 {
        self.profile
            .rounding_error_mm()
            .max(self.recess.source_error_mm())
    }

    pub fn faces(&self) -> Vec<Face> {
        let mut result = vec![Face::ThreadShoulder, Face::ThreadTip];
        result.extend((0..4 * self.thread.spec().starts).map(Face::ThreadWall));
        result.push(Face::HeadTop);
        result.extend((0..self.profile.segments().count() as u32).map(Face::ProfileWall));
        result.push(Face::SocketFloor);
        result.extend((0..self.recess.wall_count()).map(Face::SocketWall));
        result
    }

    /// Nominal source evaluation. Curved recess uncertainty remains available
    /// through the recess's cached volume interval for manufacturing proofs.
    pub fn volume_mm3(&self) -> f64 {
        self.visible_thread_volume_mm3() + self.profile.volume_mm3()
            - self.recess.nominal_volume_mm3()
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
            && self
                .recess
                .contains_void([x, y, z - self.head_top_z()], tolerance)
                == Some(false)
    }
}

fn flat_hex(profile: &Profile, spec: Spec) -> Result<Recess, ThreadError> {
    let contour = HexDrive::canonical(spec.socket_across_flats_mm, profile.rounding_error_mm())?;
    Recess::flat_hex(contour, spec.socket_depth_mm).map_err(|_| ThreadError::Invalid)
}

#[cfg(test)]
mod tests;
