use super::{Builder, check_options, hex_side, phase_strips, volume_error_bound};
use crate::{
    display::{DisplayMesh, MeshError, Tessellation},
    headed_screw::{Screw, ScrewFace, ScrewHead},
};
use std::f64::consts::TAU;
#[cfg(test)]
mod contour_tests;
mod torx;

fn profiled_options(
    screw: &crate::profiled_screw::Screw,
    options: Tessellation,
) -> Result<Tessellation, MeshError> {
    check_options(options)?;
    let reserve = 4. * screw.source_error_mm();
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
    use crate::profiled_screw::Face;
    let options = profiled_options(screw, options)?;
    let radius = screw.profile().support([1., 0., 0.]);
    let max_angle = profiled_angle(screw, options, radius);
    let angular_steps = (TAU / (12. * max_angle)).ceil() * 12.;
    if angular_steps * 2. > options.triangle_limit as f64 {
        return Err(MeshError::BudgetExceeded);
    }
    let thread = screw.thread();
    let strips = phase_strips(thread, TAU / angular_steps, options)?;
    let mut builder = Builder::new(thread, angular_steps as usize, strips, options);
    let coupling = screw.source_error_mm() / screw.recess().maximum_radius_bound_mm();
    match screw.recess() {
        crate::screw_recess::Recess::FlatHex(recess) => {
            let (drive, _) = recess.geometry();
            if !drive.is_canonical() {
                let angles = drive.corners_mm().map(|[x, y]| y.atan2(x).rem_euclid(TAU));
                builder.insert_angular_stations(&angles, coupling)?;
            }
        }
        crate::screw_recess::Recess::ConicalHexalobular(recess) => {
            let angles = torx::stations(recess.source(), options)?;
            builder.insert_angular_stations(&angles, coupling)?;
        }
    }
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
    let top_z = screw.head_top_z();
    let top = shell.profile_outer(
        screw.profile(),
        screw.spec().shank_length_mm,
        screw.visible_thread_interval_mm()[0],
        &shaft_rim,
        (face(Face::ThreadShoulder), |index| {
            face(Face::ProfileWall(index))
        }),
        options,
    )?;
    let (mouth, floor, floor_z) = match screw.recess() {
        crate::screw_recess::Recess::FlatHex(recess) => {
            let (drive, depth) = recess.geometry();
            let floor_z = top_z + depth;
            let mouth = if drive.is_canonical() {
                shell.ring(top_z, |angle| drive.radius_at_angle_mm(angle))
            } else {
                shell.drive_ring(top_z, drive, coupling)?
            };
            let floor = if drive.is_canonical() {
                shell.ring(floor_z, |angle| drive.radius_at_angle_mm(angle))
            } else {
                shell.drive_ring(floor_z, drive, coupling)?
            };
            (mouth, floor, floor_z)
        }
        crate::screw_recess::Recess::ConicalHexalobular(recess) => {
            let source = recess.source();
            (
                shell.torx_ring(top_z, source, coupling, false, options)?,
                shell.torx_ring(top_z, source, coupling, true, options)?,
                top_z + source.cone().tip_depth_mm,
            )
        }
    };
    shell.annulus(&mouth, &top, face(Face::HeadTop), true)?;
    let first_wall = face(Face::SocketWall(0));
    shell.wall(
        &mouth,
        &floor,
        |angle| {
            Ok(first_wall
                + match screw.recess() {
                    crate::screw_recess::Recess::FlatHex(recess) => {
                        let (drive, _) = recess.geometry();
                        if drive.is_canonical() {
                            hex_side(angle)
                        } else {
                            drive.side_at_angle(angle).expect("finite drive angle") as u32
                        }
                    }
                    crate::screw_recess::Recess::ConicalHexalobular(recess) => recess
                        .source()
                        .arc_index_at_angle(angle)
                        .ok_or(MeshError::UnresolvedBoundary)?
                        as u32,
                })
        },
        true,
    )?;
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
    let angle = profiled_angle(screw, options, radius);
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
    let error = volume_error_bound(screw.thread(), options)
        + (polygonal - exact).abs()
        + polygonal * (TAU / steps).powi(2) / 6.
        + roundoff;
    if let Some(source) = screw.recess().hexalobular_geometry() {
        let [lo, hi] = screw.recess().volume_bounds_mm3();
        let cavity = torx::volume_error(
            source,
            options,
            4. * screw.source_error_mm(),
            (hi - lo).next_up(),
        )?;
        Ok((crate::tapered_thread::interval::I::point(error)
            + crate::tapered_thread::interval::I::point(cavity))
        .hi)
    } else {
        Ok(error)
    }
}

fn profiled_angle(screw: &crate::profiled_screw::Screw, options: Tessellation, radius: f64) -> f64 {
    let angle = (options.angle * 0.5).min((options.deflection / radius).sqrt());
    screw
        .recess()
        .hexalobular_geometry()
        .map_or(angle, |source| torx::angular_step(source, options, angle))
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
        shell.wall(&shank_bottom, &shank_top, |_| Ok(shank_wall), false)?;
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
        |angle| Ok(first_wall + head_side(angle)),
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
        shell.wall(
            &mouth,
            &floor,
            |angle| Ok(first_wall + hex_side(angle)),
            true,
        )?;
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
    fn profile_outer(
        &mut self,
        profile: &crate::head_profile::Profile,
        plain: f64,
        visible_start: f64,
        shaft_rim: &[u32],
        faces: (u32, impl Fn(u32) -> u32),
        options: Tessellation,
    ) -> Result<Vec<u32>, MeshError> {
        use crate::head_profile::Curve;
        let shank_radius = profile.end()[0];
        let join = if shank_radius == self.builder.thread.spec().major_radius_mm {
            self.shank_join(shaft_rim, shank_radius)
        } else {
            // Distinct physical circles keep distinct source identities, even
            // when their positive separation is below a display chord request.
            self.ring(visible_start, |_| shank_radius)
        };
        self.annulus(shaft_rim, &join, faces.0, false)?;
        let top_z = profile.top()[1] - plain;
        let top = self.ring(top_z, |_| profile.top()[0]);
        let mut previous = top.clone();
        let segments: Vec<_> = profile.segments().collect();
        for (index, segment) in segments.iter().enumerate() {
            let subdivisions = match segment.curve {
                Curve::Line => 1.,
                Curve::Arc {
                    radius_mm: r,
                    sweep,
                    ..
                } => (sweep.abs() / (options.angle * 0.5).min((options.deflection / r).sqrt()))
                    .ceil(),
            };
            if subdivisions * self.angles.len() as f64 * 2. > options.triangle_limit as f64 {
                return Err(MeshError::BudgetExceeded);
            }
            for step in 1..=subdivisions as usize {
                let point = segment.point(step as f64 / subdivisions);
                let next = if index + 1 == segments.len() && step == subdivisions as usize {
                    join.clone()
                } else {
                    self.ring(point[1] - plain, |_| point[0])
                };
                self.wall(&previous, &next, |_| Ok(faces.1(index as u32)), false)?;
                previous = next;
            }
        }
        Ok(top)
    }

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

    fn drive_ring(
        &mut self,
        z: f64,
        drive: &crate::hex_drive::HexDrive,
        coupling: f64,
    ) -> Result<Vec<u32>, MeshError> {
        let corners = self.corner_stations(&drive.corners_mm(), coupling)?;
        Ok((0..self.angles.len())
            .map(|i| {
                let [x, y] = corners[i].unwrap_or_else(|| {
                    let angle = self.angles[i];
                    let radius = drive.radius_at_angle_mm(angle);
                    [radius * angle.cos(), radius * angle.sin()]
                });
                self.push([x, y, z])
            })
            .collect())
    }

    fn corner_stations(
        &self,
        points: &[[f64; 2]],
        coupling: f64,
    ) -> Result<Vec<Option<[f64; 2]>>, MeshError> {
        let mut corners = vec![None; self.angles.len()];
        for &point in points {
            let angle = point[1].atan2(point[0]).rem_euclid(TAU);
            let distance = |a: f64| (a - angle).rem_euclid(TAU).min((angle - a).rem_euclid(TAU));
            let index = (0..self.angles.len())
                .min_by(|&a, &b| distance(self.angles[a]).total_cmp(&distance(self.angles[b])))
                .ok_or(MeshError::UnresolvedBoundary)?;
            if distance(self.angles[index]) > 2. * coupling || corners[index].is_some() {
                return Err(MeshError::UnresolvedBoundary);
            }
            corners[index] = Some(point);
        }
        Ok(corners)
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
        face: impl Fn(f64) -> Result<u32, MeshError>,
        reverse: bool,
    ) -> Result<(), MeshError> {
        for i in 0..low.len() {
            let next = (i + 1) % low.len();
            let lo = self.angles[i];
            let hi = if next == 0 { TAU } else { self.angles[next] };
            let face = face((lo + hi) * 0.5)?;
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
