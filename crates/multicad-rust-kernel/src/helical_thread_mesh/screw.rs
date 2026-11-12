use super::{Builder, check_options, hex_side, phase_strips, volume_error_bound};
use crate::{
    display::{DisplayMesh, MeshError, Tessellation},
    headed_screw::{Screw, ScrewFace, ScrewHead},
};
use std::f64::consts::TAU;

fn profiled_options(
    screw: &crate::profiled_screw::Screw,
    options: Tessellation,
) -> Result<Tessellation, MeshError> {
    check_options(options)?;
    let reserve = 4. * screw.profile().rounding_error_mm();
    if options.deflection <= reserve {
        return Err(MeshError::UnresolvedBoundary);
    }
    Ok(Tessellation {
        deflection: options.deflection - reserve,
        ..options
    })
}

/// Retain the existing thread's indexed shoulder while sweeping every real
/// head segment. No display envelope or mesh boolean enters the source shape.
pub fn tessellate_profiled_screw(
    screw: &crate::profiled_screw::Screw,
    options: Tessellation,
) -> Result<DisplayMesh, MeshError> {
    use crate::{head_profile::Curve, profiled_screw::Face};
    let options = profiled_options(screw, options)?;
    let radius = screw.profile().support([1., 0., 0.]);
    let max_angle = (options.angle * 0.5).min((options.deflection / radius).sqrt());
    let angular_steps = (TAU / (12. * max_angle)).ceil() * 12.;
    if angular_steps * 2. > options.triangle_limit as f64 {
        return Err(MeshError::BudgetExceeded);
    }
    let thread = screw.thread();
    let strips = phase_strips(thread, TAU / angular_steps, options)?;
    let mut builder = Builder::new(thread, angular_steps as usize, strips, options);
    builder.lower = screw.visible_thread_interval_mm()[0];
    builder.build_wall()?;
    let [mut shoulder, tip] = std::mem::take(&mut builder.caps);
    builder.add_cap(1, tip)?;
    Builder::sort_ring(&mut shoulder)?;
    let faces = screw.faces();
    let face = |f| {
        faces
            .iter()
            .position(|&candidate| candidate == f)
            .expect("profiled screw face") as u32
    };
    let mut shell = HeadMesh {
        builder: &mut builder,
        angles: shoulder.iter().map(|&(angle, _)| angle).collect(),
    };
    let shaft_rim: Vec<_> = shoulder.into_iter().map(|(_, id)| id).collect();
    let shank_radius = screw.profile().end()[0];
    let join = if shank_radius == thread.spec().major_radius_mm {
        shell.shank_join(&shaft_rim, shank_radius)
    } else {
        // Distinct physical circles keep distinct source identities, even
        // when their positive separation is below a display chord request.
        shell.ring(screw.visible_thread_interval_mm()[0], |_| shank_radius)
    };
    shell.annulus(&shaft_rim, &join, face(Face::ThreadShoulder), false)?;
    let top_z = screw.head_top_z();
    let top = shell.ring(top_z, |_| screw.profile().top()[0]);
    let mut previous = top.clone();
    let segments: Vec<_> = screw.profile().segments().collect();
    for (index, segment) in segments.iter().enumerate() {
        let subdivisions = match segment.curve {
            Curve::Line => 1.,
            Curve::Arc {
                radius_mm: r,
                sweep,
                ..
            } => (sweep.abs() / (options.angle * 0.5).min((options.deflection / r).sqrt())).ceil(),
        };
        if subdivisions * shell.angles.len() as f64 * 2. > options.triangle_limit as f64 {
            return Err(MeshError::BudgetExceeded);
        }
        for step in 1..=subdivisions as usize {
            let point = segment.point(step as f64 / subdivisions);
            let next = if index + 1 == segments.len() && step == subdivisions as usize {
                join.clone()
            } else {
                shell.ring(point[1] - screw.spec().shank_length_mm, |_| point[0])
            };
            shell.wall(
                &previous,
                &next,
                |_| face(Face::ProfileWall(index as u32)),
                false,
            )?;
            previous = next;
        }
    }
    let mouth = shell.ring(top_z, |angle| screw.socket_radius(angle));
    let floor_z = top_z + screw.spec().socket_depth_mm;
    let floor = shell.ring(floor_z, |angle| screw.socket_radius(angle));
    shell.annulus(&mouth, &top, face(Face::HeadTop), true)?;
    let first_wall = face(Face::SocketWall(0));
    shell.wall(&mouth, &floor, |angle| first_wall + hex_side(angle), true)?;
    shell.disk(&floor, floor_z, face(Face::SocketFloor))?;
    Ok(builder.mesh)
}

/// Bound the actual frusta and circular chords of the head independently of
/// the thread interpolation bound. The hex drive has exact planar walls.
pub fn profiled_screw_volume_error_bound(
    screw: &crate::profiled_screw::Screw,
    options: Tessellation,
) -> Result<f64, MeshError> {
    use crate::head_profile::{Curve, Segment};
    let options = profiled_options(screw, options)?;
    let radius = screw.profile().support([1., 0., 0.]);
    let angle = (options.angle * 0.5).min((options.deflection / radius).sqrt());
    let steps = (TAU / (12. * angle)).ceil() * 12.;
    if steps * 2. > options.triangle_limit as f64 {
        return Err(MeshError::BudgetExceeded);
    }
    let mut polygonal = 0.;
    for segment in screw.profile().segments() {
        let count = match segment.curve {
            Curve::Line => 1.,
            Curve::Arc {
                radius_mm: r,
                sweep,
                ..
            } => (sweep.abs() / (options.angle * 0.5).min((options.deflection / r).sqrt())).ceil(),
        };
        // Bound work before converting a huge refinement to usize.
        if count * steps * 2. > options.triangle_limit as f64 {
            return Err(MeshError::BudgetExceeded);
        }
        for i in 0..count as usize {
            polygonal += Segment {
                start: segment.point(i as f64 / count),
                end: segment.point((i + 1) as f64 / count),
                curve: Curve::Line,
                role: segment.role,
            }
            .volume_mm3();
        }
    }
    let exact = screw.profile().volume_mm3();
    // Sum (angle - sin(angle)) <= max_angle² * sum(angle) / 6.
    // Extra thread cap stations only subdivide these angular intervals.
    let surfaces = screw
        .profile()
        .segments()
        .map(Segment::area_mm2)
        .sum::<f64>()
        + std::f64::consts::PI
            * (screw.profile().top()[0].powi(2) + screw.profile().end()[0].powi(2));
    let roundoff = 4. * screw.profile().rounding_error_mm() * surfaces
        + 256.
            * f64::EPSILON
            * (exact.abs()
                + polygonal.abs()
                + radius.powi(2) * (screw.profile().end()[1] - screw.profile().top()[1]));
    Ok(volume_error_bound(screw.thread(), options)
        + (polygonal - exact).abs()
        + polygonal * (TAU / steps).powi(2) / 6.
        + roundoff)
}

/// Mesh one connected screw shell using the same trimmed shaft rim at the
/// attachment. Angular columns include the exact corners of the head or socket.
pub fn tessellate_screw(screw: &Screw, options: Tessellation) -> Result<DisplayMesh, MeshError> {
    check_options(options)?;
    let thread = screw.thread();
    let spec = screw.spec();
    let radius = match spec.head {
        ScrewHead::SocketCap { radius_mm, .. } => radius_mm,
        ScrewHead::ExternalHex {
            across_flats_mm, ..
        } => across_flats_mm / 3_f64.sqrt(),
    };
    let max_angle = (options.angle * 0.5).min((2. * options.deflection / radius).sqrt());
    let angular_steps = (TAU / (12. * max_angle)).ceil() * 12.;
    if angular_steps * 2. > options.triangle_limit as f64 {
        return Err(MeshError::BudgetExceeded);
    }
    let steps = angular_steps as usize;
    let strips = phase_strips(thread, TAU / angular_steps, options)?;
    let mut builder = Builder::new(thread, steps, strips, options);
    builder.build_wall()?;
    let [mut shoulder, tip] = std::mem::take(&mut builder.caps);
    builder.add_cap(1, tip)?;
    Builder::sort_ring(&mut shoulder)?;
    let faces = screw.faces();
    let face = |kind| {
        faces
            .iter()
            .position(|&candidate| candidate == kind)
            .expect("screw face role") as u32
    };
    let [seat_z, top_z] = [screw.head_seat_z(), screw.head_top_z()];
    let mut shell = HeadMesh {
        builder: &mut builder,
        angles: shoulder.iter().map(|&(angle, _)| angle).collect(),
    };
    let shaft_rim: Vec<_> = shoulder.into_iter().map(|(_, id)| id).collect();
    let head_seat = shell.ring(seat_z, |angle| screw.head_radius(angle));
    if spec.shank_length_mm > 0. {
        let shank_top = shell.shank_join(&shaft_rim, spec.shank_radius_mm);
        let shank_bottom = shell.ring(seat_z, |_| spec.shank_radius_mm);
        shell.annulus(
            &shaft_rim,
            &shank_top,
            face(ScrewFace::ThreadShoulder),
            false,
        )?;
        let shank_wall = face(ScrewFace::ShankWall);
        shell.wall(&shank_bottom, &shank_top, |_| shank_wall, false)?;
        shell.annulus(&shank_bottom, &head_seat, face(ScrewFace::HeadSeat), false)?;
    } else {
        shell.annulus(
            &shaft_rim,
            &head_seat,
            face(ScrewFace::ThreadShoulder),
            false,
        )?;
    }
    let head_top = shell.ring(top_z, |angle| screw.head_radius(angle));
    let first_wall = face(ScrewFace::HeadWall(0));
    let head_side = match spec.head {
        ScrewHead::SocketCap { .. } => |_| 0,
        ScrewHead::ExternalHex { .. } => hex_side,
    };
    shell.wall(
        &head_top,
        &head_seat,
        |angle| first_wall + head_side(angle),
        false,
    )?;
    if let ScrewHead::SocketCap {
        socket_depth_mm, ..
    } = spec.head
    {
        let floor_z = top_z + socket_depth_mm;
        let socket_radius = |angle| screw.socket_radius(angle).expect("socket head");
        let mouth = shell.ring(top_z, socket_radius);
        let floor = shell.ring(floor_z, socket_radius);
        shell.annulus(&mouth, &head_top, face(ScrewFace::HeadTop), true)?;
        let first_wall = face(ScrewFace::SocketWall(0));
        shell.wall(&mouth, &floor, |angle| first_wall + hex_side(angle), true)?;
        shell.disk(&floor, floor_z, face(ScrewFace::SocketFloor))?;
    } else {
        shell.disk(&head_top, top_z, face(ScrewFace::HeadTop))?;
    }
    Ok(builder.mesh)
}

/// Combine the helical interpolation bound with the area lost by approximating
/// circular walls. The hex faces and recess are planar and use exact corners.
pub fn screw_volume_error_bound(screw: &Screw, options: Tessellation) -> f64 {
    let spec = screw.spec();
    let head_wall = match spec.head {
        ScrewHead::SocketCap {
            radius_mm,
            height_mm,
            ..
        } => radius_mm * height_mm,
        ScrewHead::ExternalHex { .. } => 0.,
    };
    volume_error_bound(screw.thread(), options)
        + TAU * (head_wall + spec.shank_radius_mm * spec.shank_length_mm) * options.deflection
}

struct HeadMesh<'a, 'thread> {
    builder: &'a mut Builder<'thread>,
    angles: Vec<f64>,
}

impl HeadMesh<'_, '_> {
    fn push(&mut self, point: [f64; 3]) -> u32 {
        let id = (self.builder.mesh.vertices.len() / 3) as u32;
        self.builder.mesh.vertices.extend(point);
        id
    }

    fn ring(&mut self, z: f64, radius: impl Fn(f64) -> f64) -> Vec<u32> {
        (0..self.angles.len())
            .map(|index| {
                let angle = self.angles[index];
                let (sin, cos) = angle.sin_cos();
                let radius = radius(angle);
                self.push([radius * cos, radius * sin, z])
            })
            .collect()
    }

    fn shank_join(&mut self, inner: &[u32], radius: f64) -> Vec<u32> {
        inner
            .iter()
            .enumerate()
            .map(|(index, &id)| {
                let [x, y, z] = self.builder.mesh.vertices.as_chunks::<3>().0[id as usize];
                if (radius - x.hypot(y)).abs() <= 64. * f64::EPSILON * radius {
                    // A full crest may coincide with a nominal-diameter shank.
                    // Give those stations one vertex identity on both walls.
                    id
                } else {
                    let (sin, cos) = self.angles[index].sin_cos();
                    self.push([radius * cos, radius * sin, z])
                }
            })
            .collect()
    }

    fn annulus(
        &mut self,
        inner: &[u32],
        outer: &[u32],
        face: u32,
        reverse: bool,
    ) -> Result<(), MeshError> {
        for i in 0..inner.len() {
            let next = (i + 1) % inner.len();
            let [a, b, c, d] = [inner[i], inner[next], outer[i], outer[next]];
            if a != c {
                self.triangle([a, c, d], face, reverse)?;
            }
            if b != d {
                self.triangle([a, d, b], face, reverse)?;
            }
        }
        Ok(())
    }

    fn wall(
        &mut self,
        low: &[u32],
        high: &[u32],
        face: impl Fn(f64) -> u32,
        reverse: bool,
    ) -> Result<(), MeshError> {
        for i in 0..low.len() {
            let next = (i + 1) % low.len();
            let lo = self.angles[i];
            let hi = if next == 0 { TAU } else { self.angles[next] };
            let face = face((lo + hi) * 0.5);
            let [a, b, c, d] = [low[i], low[next], high[i], high[next]];
            self.triangle([a, b, d], face, reverse)?;
            self.triangle([a, d, c], face, reverse)?;
        }
        Ok(())
    }

    fn disk(&mut self, rim: &[u32], z: f64, face: u32) -> Result<(), MeshError> {
        let center = self.push([0., 0., z]);
        for (&a, &b) in rim.iter().zip(rim.iter().cycle().skip(1)) {
            self.builder.triangle([center, b, a], face)?;
        }
        Ok(())
    }

    fn triangle(
        &mut self,
        mut triangle: [u32; 3],
        face: u32,
        reverse: bool,
    ) -> Result<(), MeshError> {
        if reverse {
            triangle.swap(1, 2);
        }
        self.builder.triangle(triangle, face)
    }
}
