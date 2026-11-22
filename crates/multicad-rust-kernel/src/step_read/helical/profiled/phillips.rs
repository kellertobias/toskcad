//! Caller-supplied source proof only; never foreign authoring reconstruction.
use super::*;
use crate::step_read::phillips_cutter::{
    PhillipsReadError, PhillipsSourceKnownProof, SourceFaceView,
};
use crate::{phillips_drive::PhillipsScrew, step::helical::HeadSource};
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhillipsScrewSourceKnownProof {
    body_faces: usize,
    cavity: PhillipsSourceKnownProof,
}
impl PhillipsScrewSourceKnownProof {
    pub fn body_face_count(&self) -> usize {
        self.body_faces
    }
    /// The reference cavity includes its virtual Mouth closure. The file's
    /// material head retains only the other28 faces and the shared top hole.
    pub fn reference_cavity(&self) -> PhillipsSourceKnownProof {
        self.cavity
    }
    pub fn material_cavity_face_count(&self) -> usize {
        self.cavity.face_count() - 1
    }
}
fn cavity(
    shell: &mut Shell<'_, '_>,
    source: &PhillipsScrew,
    faces: &[u32],
    placement: [Point; 3],
) -> Result<PhillipsSourceKnownProof> {
    let (surface, forward, bounds) = face(shell, faces[0])?;
    require(forward && bounds.len() == 2)?;
    let uses = ring(shell, bounds[1], false)?
        .into_iter()
        .rev()
        .map(|(id, forward)| (id, !forward))
        .collect();
    let mut views = vec![SourceFaceView {
        surface,
        forward,
        uses,
    }];
    for &id in &faces[source.profile().segments().count() + 1..] {
        let (surface, forward, bounds) = face(shell, id)?;
        let [bound] = bounds.as_slice() else {
            return Err(ReadError::Unsupported);
        };
        let uses = ring(shell, *bound, true)?
            .into_iter()
            .rev()
            .map(|(id, forward)| (id, !forward))
            .collect();
        views.push(SourceFaceView {
            surface,
            forward: !forward,
            uses,
        });
    }
    let frame = Frame::new(placement)?;
    let pose = crate::step::phillips_cutter::support::Pose::new(
        frame.world([0.; 2], source.head_top_z()),
        placement[1],
        placement[2],
    )
    .map_err(|_| ReadError::Invalid)?;
    crate::step_read::phillips_cutter::verify_component(
        shell,
        &views,
        source.cavity().source(),
        &pose,
        shell.tolerance.min(1e-7),
    )
}
pub(in crate::step_read::helical) fn verify(
    proof: Proof,
    shell: &mut Shell<'_, '_>,
    all_faces: &[u32],
    source: &PhillipsScrew,
    rings: &[Vec<u32>; 2],
    corners: &[[[u32; 2]; 2]],
    placement: [Point; 3],
) -> Result<f64> {
    let outer = outer::verify(
        proof,
        shell,
        all_faces,
        HeadSource::phillips(source),
        rings,
        corners,
        placement,
    )?;
    cavity(shell, source, &outer.faces, placement)?;
    // Explicit nominal SOURCE evaluation after actual graph/carrier checks.
    // Public source-known proof discards this internal scalar; it is neither
    // independent decoded-file mass nor an imported body constructor.
    Ok(source.profile().volume_mm3()
        - crate::tapered_thread::interval::I::new(
            source.cavity().volume().volume_mm3[0],
            source.cavity().volume().volume_mm3[1],
        )
        .midpoint())
}
/// Prove the supplied source against this actual file, including full thread
/// carriers, shared finite incidence and real cavity supports/curves. No mass,
/// imported source, original tool Domain or native admission is inferred.
pub fn verify_phillips_screw_source_known(
    input: &str,
    source: &PhillipsScrew,
    origin: Point,
    u: Point,
    v: Point,
    tolerance: f64,
) -> std::result::Result<PhillipsScrewSourceKnownProof, PhillipsReadError> {
    let placement = [origin, u, v];
    let mut results = compound::read_members_by(input, tolerance, true, true, |shell, id, _| {
        let expected = crate::step::helical::Closure::Profiled(HeadSource::phillips(source));
        let _nominal_source_evaluation =
            super::super::verify_member(shell, id, source.thread(), expected, placement, None)?;
        let faces = shell.member_faces(id)?;
        let shoulders = if source.profile().end()[0] == source.thread().spec().major_radius_mm {
            source.thread().spec().starts as usize
        } else {
            1
        };
        let offset = 4 * source.thread().spec().starts as usize + shoulders + 1;
        let cavity = cavity(shell, source, &faces[offset..], placement)?;
        Ok(PhillipsScrewSourceKnownProof {
            body_faces: faces.len(),
            cavity,
        })
    })
    .map_err(PhillipsReadError::from)?;
    results.pop().ok_or(PhillipsReadError::Invalid)
}

#[cfg(test)]
pub(crate) fn verify_cavity_fixture(input: &str, source: &PhillipsScrew) -> Result<()> {
    compound::read_members_by(input, 1e-7, true, true, |shell, id, _| {
        let faces = shell.member_faces(id)?;
        let count = 4 * source.thread().spec().starts as usize + 2;
        cavity(
            shell,
            source,
            &faces[count..],
            [[0.; 3], [1., 0., 0.], [0., 1., 0.]],
        )?;
        Ok(())
    })
    .map(|_| ())
}
