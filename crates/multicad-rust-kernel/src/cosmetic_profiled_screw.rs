//! A distinct smooth-shaft envelope retaining the actual physical head/recess.
//! This is cosmetic geometry; it never earns modeled-thread manufacturing.
use crate::{
    display::{DisplayMesh, MeshError, Tessellation},
    profiled_screw::{Face, Screw},
};
use std::f64::consts::PI;

#[derive(Clone, Debug, PartialEq)]
pub struct Envelope {
    source: Screw,
    faces: Vec<Face>,
}
impl Envelope {
    pub fn new(source: Screw) -> Self {
        let shoulder = source.profile().end()[0] != source.thread().spec().major_radius_mm;
        let faces = source
            .faces()
            .into_iter()
            .filter(|f| match f {
                Face::ThreadWall(index) => *index == 0,
                Face::ThreadShoulder => shoulder,
                _ => true,
            })
            .collect();
        Self { source, faces }
    }
    /// The immutable original source supplies head, recess and junction data.
    /// Its helical wall is not the geometry of this envelope.
    pub fn head_source(&self) -> &Screw {
        &self.source
    }
    pub fn faces(&self) -> &[Face] {
        &self.faces
    }
    pub fn shaft_radius_mm(&self) -> f64 {
        self.source.thread().spec().major_radius_mm
    }
    pub fn shaft_interval_mm(&self) -> [f64; 2] {
        self.source.visible_thread_interval_mm()
    }
    pub fn volume_mm3(&self) -> f64 {
        let [a, b] = self.shaft_interval_mm();
        self.source.profile().volume_mm3() + PI * self.shaft_radius_mm().powi(2) * (b - a)
            - self.source.recess().nominal_volume_mm3()
    }
    pub fn contains(&self, point: [f64; 3], tolerance: f64) -> bool {
        if !point.into_iter().all(f64::is_finite) || !tolerance.is_finite() || tolerance < 0. {
            return false;
        }
        let [x, y, z] = point;
        let r = x.hypot(y);
        let [a, b] = self.shaft_interval_mm();
        let profile = self.source.profile();
        let local = z + self.source.spec().shank_length_mm;
        let inside_shaft =
            z >= a - tolerance && z <= b + tolerance && r <= self.shaft_radius_mm() + tolerance;
        let inside_head = local >= profile.top()[1] - tolerance
            && local <= profile.end()[1] + tolerance
            && r <= profile
                .radius_at(local.clamp(profile.top()[1], profile.end()[1]))
                .expect("retained source profile domain")
                + tolerance;
        (inside_shaft || inside_head)
            && self
                .source
                .recess()
                .contains_void([x, y, z - self.source.head_top_z()], tolerance)
                == Some(false)
    }
    pub fn support(&self, direction: [f64; 3]) -> f64 {
        if !direction.into_iter().all(f64::is_finite) {
            return f64::NAN;
        }
        let [a, b] = self.shaft_interval_mm();
        let cylinder = self.shaft_radius_mm() * direction[0].hypot(direction[1])
            + direction[2] * if direction[2] < 0. { a } else { b };
        cylinder.max(
            self.source.profile().support(direction)
                - direction[2] * self.source.spec().shank_length_mm,
        )
    }
    pub fn mesh(&self, options: Tessellation) -> Result<DisplayMesh, MeshError> {
        let mut mesh =
            crate::helical_thread_mesh::tessellate_profiled_cosmetic_screw(&self.source, options)?;
        let old = self.source.faces();
        for id in &mut mesh.face_indices {
            *id = self
                .faces
                .iter()
                .position(|f| *f == old[*id as usize])
                .expect("cosmetic mesh emits only its retained faces") as u32;
        }
        Ok(mesh)
    }
}
