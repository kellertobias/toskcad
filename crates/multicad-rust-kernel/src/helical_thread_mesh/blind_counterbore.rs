//! Shared-index finite thread shoulder, seat cylinder and blind floor.
use super::*;
use crate::thread_counterbore::Counterbore;
pub(crate) fn open_counterbored_cavity(
    seat: &Counterbore,
    options: Tessellation,
) -> Result<(DisplayMesh, Vec<u32>), MeshError> {
    let mut builder = envelope_builder(
        seat.thread(),
        Some(BoreExterior::Cylinder(seat.radius_mm())),
        options,
    )?;
    builder.lower = seat.depth_mm();
    builder.build_wall()?;
    let [shoulder, floor] = std::mem::take(&mut builder.caps);
    builder.add_cap(1, floor)?;
    // Thread wall and floor bound the surrounding stock rather than a male.
    for triangle in builder.mesh.indices.as_chunks_mut::<3>().0 {
        triangle.swap(1, 2);
    }
    let first = builder.mesh.face_indices.len();
    let rim = builder.add_annular_cap(0, shoulder, BoreExterior::Cylinder(seat.radius_mm()))?;
    for face in &mut builder.mesh.face_indices[first..] {
        *face = seat.shoulder_face();
    }
    let bottom = &rim.vertices[..rim.vertices.len() - 1];
    let mut mouth = Vec::with_capacity(bottom.len());
    for &index in bottom {
        let p = builder.mesh.vertices.as_chunks::<3>().0[index as usize];
        mouth.push(builder.mesh.vertices.len() as u32 / 3);
        builder.mesh.vertices.extend([p[0], p[1], 0.]);
    }
    for i in 0..bottom.len() {
        let next = (i + 1) % bottom.len();
        let [a, b, c, d] = [mouth[i], mouth[next], bottom[i], bottom[next]];
        builder.triangle([a, d, b], seat.cylinder_face())?;
        builder.triangle([a, c, d], seat.cylinder_face())?;
    }
    Ok((builder.mesh, mouth))
}

/// The original full-thread bound overestimates its retained suffix. The
/// additional circular wall has radial displacement at most deflection.
pub fn counterbore_volume_error_bound(seat: &Counterbore, options: Tessellation) -> f64 {
    let circle = (((TAU.next_up() * seat.radius_mm()).next_up() * seat.depth_mm()).next_up()
        * options.deflection)
        .next_up();
    (volume_error_bound(seat.thread(), options) + circle).next_up()
}
