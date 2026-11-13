//! Typed authored composition preceding native and physical-file admission.
use super::*;
mod geometry;
use crate::{
    head_profile::Profile,
    helical_thread::Thread,
    profiled_screw::{SourceSpec, qualify_outer_source},
};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ScrewFace {
    ThreadShoulder,
    ThreadTip,
    ThreadWall(u32),
    HeadTop,
    ProfileWall(u32),
    /// The actual cavity face. The cutter's temporary Mouth cap is absent.
    Cavity(FaceKey),
}
/// Actual authored Phillips screw, preserving its original helical source and
/// finite cutter. This composition proves source containment, not native face
/// properties, head-mouth manufacturing sewing or independent STEP import.
#[derive(Clone, Debug)]
pub struct PhillipsScrew {
    thread: Thread,
    profile: Profile,
    spec: SourceSpec,
    cavity: Cavity,
    host: HostEnvelopeProof,
    visible_start_mm: f64,
}
impl PhillipsScrew {
    pub fn new(
        thread: Thread,
        profile: Profile,
        spec: SourceSpec,
        cavity: Cavity,
        work_limit: usize,
        cancelled: &dyn Fn() -> bool,
    ) -> Result<Self, Error> {
        let visible_start_mm = qualify_outer_source(&thread, &profile, spec.shank_length_mm, true)
            .map_err(|_| Error::Invalid)?;
        let host = cavity.certify_host_envelope(
            &profile,
            &thread,
            spec.shank_length_mm,
            spec.minimum_wall_mm,
            work_limit,
            cancelled,
        )?;
        Ok(Self {
            thread,
            profile,
            spec,
            cavity,
            host,
            visible_start_mm,
        })
    }
    pub fn thread(&self) -> &Thread {
        &self.thread
    }
    pub fn profile(&self) -> &Profile {
        &self.profile
    }
    pub fn spec(&self) -> SourceSpec {
        self.spec
    }
    pub fn cavity(&self) -> &Cavity {
        &self.cavity
    }
    pub fn host_envelope(&self) -> HostEnvelopeProof {
        self.host
    }
    pub fn visible_thread_interval_mm(&self) -> [f64; 2] {
        [self.visible_start_mm, self.thread.spec().length_mm]
    }
    pub fn head_top_z(&self) -> f64 {
        self.profile.top()[1] - self.spec.shank_length_mm
    }
    /// Nominal analytic source evaluation, not a physical-file mass enclosure.
    pub fn nominal_volume_mm3(&self) -> f64 {
        let length = self.thread.spec().length_mm;
        self.thread.volume_mm3() * ((length - self.visible_start_mm) / length)
            + self.profile.volume_mm3()
            - I::new(
                self.cavity.volume.volume_mm3[0],
                self.cavity.volume.volume_mm3[1],
            )
            .midpoint()
    }
    pub fn faces(&self) -> Vec<ScrewFace> {
        let mut faces = vec![ScrewFace::ThreadShoulder, ScrewFace::ThreadTip];
        faces.extend((0..4 * self.thread.spec().starts).map(ScrewFace::ThreadWall));
        faces.push(ScrewFace::HeadTop);
        faces.extend((0..self.profile.segments().count() as u32).map(ScrewFace::ProfileWall));
        faces.extend(
            self.cavity
                .side_face_keys()
                .iter()
                .copied()
                .map(ScrewFace::Cavity),
        );
        faces
    }
}

#[cfg(test)]
mod tests;
