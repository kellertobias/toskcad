//! Physical finite head patches. No authored head metadata is consulted.
use super::*;
use crate::{
    head_profile::{Curve as HeadCurve, Profile, Role, Segment},
    profiled_screw::Screw,
};
use std::f64::consts::{PI, TAU};
mod decode;
pub(in super::super) use decode::read;

/// Absolute head/drive mass uncertainty from the retained source precision.
/// The separately integrated thread spline error is an additional budget.
pub fn profiled_screw_head_volume_error_bound_mm3(screw: &Screw) -> Option<f64> {
    use crate::tapered_thread::interval::I;
    let profile = screw.profile();
    let error = profile.rounding_error_mm();
    let mut length = profile.top()[0] + profile.end()[0] + profile.end()[1] - profile.top()[1];
    let mut count = 3.;
    for segment in profile.segments() {
        let part = match segment.curve {
            HeadCurve::Line => {
                (segment.end[0] - segment.start[0]).hypot(segment.end[1] - segment.start[1])
            }
            HeadCurve::Arc {
                radius_mm, sweep, ..
            } => {
                if radius_mm <= 32. * error {
                    return None;
                }
                (radius_mm * sweep.abs()).next_up()
            }
        };
        length = (length + part.next_up()).next_up();
        count += 1.;
    }
    // Support-center and trim-angle conversion each spend source error. With
    // r>32e, inverse angular conditioning gives a complete meridian deviation
    // below16e. Lines/caps are bounded by their two endpoints. The area of a
    // radius-d tube around a finite planar curve is at most2Ld+pi*d² per
    // segment; revolving it multiplies by at most2pi*(R+d). This also covers
    // concave necks, without assuming a convex parallel-body formula.
    let d = I::point(16.) * I::point(error);
    let pi = I::new(PI.next_down(), PI.next_up());
    let r = I::point(profile.support([1., 0., 0.]).next_up()) + d;
    let outer = I::point(2.)
        * pi
        * r
        * (I::point(2.) * I::point(length) * d + I::point(count) * pi * d.square());
    let af = screw.spec().socket_across_flats_mm;
    let depth = screw.spec().socket_depth_mm;
    let width = I::new(
        (af - d.hi).max(0.).next_down().max(0.),
        (af + d.hi).next_up(),
    );
    let depth = I::new(
        (depth - d.hi).max(0.).next_down().max(0.),
        (depth + d.hi).next_up(),
    );
    let root3 = I::point(3.).sqrt();
    let void = root3 * I::point(0.5) * width.square() * depth;
    let result = (outer.hi + void.width()).next_up();
    (result.is_finite() && result > 0.).then_some(result)
}

#[derive(Clone, Copy)]
struct Frame {
    origin: Point,
    u: Point,
    v: Point,
    n: Point,
    rows: [Point; 3],
    determinant: f64,
}
impl Frame {
    fn new([origin, u, v]: [Point; 3]) -> Result<Self> {
        let n = cross(u, v);
        let rows = [cross(v, n), cross(n, u), cross(u, v)];
        let determinant = dot(u, rows[0]);
        require(determinant.is_finite() && determinant > 0.)?;
        Ok(Self {
            origin,
            u,
            v,
            n,
            rows,
            determinant,
        })
    }
    fn local(self, p: Point) -> Point {
        self.rows
            .map(|r| dot(sub(p, self.origin), r) / self.determinant)
    }
    fn world(self, [x, y]: [f64; 2], z: f64) -> Point {
        std::array::from_fn(|i| self.origin[i] + x * self.u[i] + y * self.v[i] + z * self.n[i])
    }
}

#[derive(Clone, Copy)]
struct Rim {
    edge: Option<u32>,
    vertex: u32,
    radius: f64,
    z: f64,
    positive: bool,
    center: Point,
}

fn decode_rim(shell: &mut Shell<'_, '_>, edge: u32, frame: Frame) -> Result<Rim> {
    let [a, b] = shell.edge(edge)?;
    require(a == b)?;
    let Curve::Circle(circle) = shell.edges[&edge].curve else {
        return Err(ReadError::Invalid);
    };
    let local = frame.local(circle.center_mm);
    require(
        local[0].hypot(local[1]) <= shell.tolerance / 16.
            && length(cross(circle.normal, frame.n)) < 1e-10,
    )?;
    Ok(Rim {
        edge: Some(edge),
        vertex: a,
        radius: circle.radius_mm,
        z: local[2],
        positive: dot(circle.normal, frame.n) > 0.,
        center: circle.center_mm,
    })
}

struct Junction {
    positive: Vec<EdgeRef>,
    vertex: u32,
}

fn junction(
    proof: Proof,
    shell: &mut Shell<'_, '_>,
    faces: &[u32],
    thread: &Thread,
    rings: &[Vec<u32>; 2],
    corners: &[[[u32; 2]; 2]],
    frame: Frame,
) -> Result<Junction> {
    let spec = thread.spec();
    let starts = spec.starts as usize;
    let shoulders: Vec<_> = std::iter::once(faces[0])
        .chain(faces[2..starts + 1].iter().copied())
        .collect();
    let [_, start, _, _, crest, _] = thread.profile_breaks();
    let rate = thread.helical_rate();
    let arc_axis = frame.n.map(|v| if spec.left_hand { v } else { -v });
    let mut crests = Vec::new();
    let mut complements = Vec::new();
    for (k, face_id) in shoulders.into_iter().enumerate() {
        let (surface, forward, bounds) = face(shell, face_id)?;
        require(forward && bounds.len() == 1)?;
        proof.axis(
            shell,
            hex_head::plane(shell, surface)?,
            frame.origin,
            frame.n,
            frame.u,
        )?;
        let edges = orient(&ring(shell, bounds[0], true)?, !spec.left_hand);
        require(
            edges.len() == 4
                && edges[..3]
                    == [
                        (rings[0][4 * k], true),
                        (rings[0][4 * k + 1], true),
                        (rings[0][4 * k + 2], true),
                    ],
        )?;
        let (complement, sense) = edges[3];
        require(!sense)?;
        let first = corners[4 * k][0][0];
        let last = corners[4 * k + 3][0][0];
        let next = corners[(4 * (k + 1)) % rings[0].len()][0][0];
        let crest_edge = rings[0][4 * k + 3];
        tangent_shank::arc(
            shell,
            crest_edge,
            [last, next],
            frame.origin,
            arc_axis,
            (
                proof,
                spec.major_radius_mm,
                (start + spec.pitch_mm - crest) / rate.abs(),
            ),
        )?;
        tangent_shank::arc(
            shell,
            complement,
            [first, last],
            frame.origin,
            arc_axis,
            (proof, spec.major_radius_mm, (crest - start) / rate.abs()),
        )?;
        crests.push(crest_edge);
        complements.push(complement);
    }
    let mut positive: Vec<_> = complements
        .into_iter()
        .zip(crests)
        .flat_map(|(a, b)| [(a, true), (b, true)])
        .collect();
    if !spec.left_hand {
        positive = positive.into_iter().rev().map(|(e, s)| (e, !s)).collect();
    }
    Ok(Junction {
        positive,
        vertex: corners[0][0][0],
    })
}

fn single_bound_rim(
    shell: &mut Shell<'_, '_>,
    bound: u32,
    outer: bool,
    positive: bool,
    frame: Frame,
) -> Result<Rim> {
    let edges = ring(shell, bound, outer)?;
    let [(edge, sense)] = edges[..] else {
        return Err(ReadError::Invalid);
    };
    let rim = decode_rim(shell, edge, frame)?;
    require(sense == (positive == rim.positive))?;
    Ok(rim)
}

fn plane_at(
    proof: Proof,
    shell: &Shell<'_, '_>,
    surface: u32,
    rim: Rim,
    frame: Frame,
    positive: bool,
) -> Result<()> {
    proof
        .axis(
            shell,
            hex_head::plane(shell, surface)?,
            rim.center,
            frame.n.map(|v| if positive { v } else { -v }),
            frame.u,
        )
        .map(|_| ())
}

fn periodic_end(
    shell: &mut Shell<'_, '_>,
    bound: u32,
    previous: Rim,
    last: bool,
    join: &Junction,
    thread: &Thread,
    frame: Frame,
) -> Result<(Rim, u32)> {
    let mut edges = ring(shell, bound, true)?;
    let start = (previous.edge.ok_or(ReadError::Invalid)?, previous.positive);
    let index = edges
        .iter()
        .position(|&edge| edge == start)
        .ok_or(ReadError::Invalid)?;
    edges.rotate_left(index);
    require(edges.len() >= 4)?;
    let seam = edges[1].0;
    require(edges[1].1 && edges.last() == Some(&(seam, false)))?;
    let middle = &edges[2..edges.len() - 1];
    let end = if last {
        require(middle == orient(&join.positive, false))?;
        Rim {
            edge: None,
            vertex: join.vertex,
            radius: thread.spec().major_radius_mm,
            z: 0.,
            positive: true,
            center: frame.origin,
        }
    } else {
        let [(edge, sense)] = middle else {
            return Err(ReadError::Invalid);
        };
        let rim = decode_rim(shell, *edge, frame)?;
        require(*sense != rim.positive)?;
        rim
    };
    let [Value::String(_), a, b, _, Value::Enum("T")] = shell.args(seam, "EDGE_CURVE")? else {
        return Err(ReadError::Invalid);
    };
    require([reference(a)?, reference(b)?] == [previous.vertex, end.vertex])?;
    Ok((end, seam))
}

pub(super) fn verify(
    proof: Proof,
    shell: &mut Shell<'_, '_>,
    all_faces: &[u32],
    expected: &Screw,
    rings: &[Vec<u32>; 2],
    corners: &[[[u32; 2]; 2]],
    placement: [Point; 3],
) -> Result<f64> {
    closed_incidence(shell, all_faces)?;
    // The generic tensor proof has its own approximation budget. Reconstructing
    // an analytic head may use only its intrinsic source arithmetic reserve.
    let proof = Proof {
        error: proof.error.map(|_| expected.profile().rounding_error_mm()),
        ..proof
    };
    let frame = Frame::new(placement)?;
    let thread = expected.thread();
    let count = 4 * thread.spec().starts as usize;
    let tail = &all_faces[count..];
    let join = junction(proof, shell, tail, thread, rings, corners, frame)?;
    let faces = &tail[thread.spec().starts as usize + 1..];
    if let Some(error) = proof.error {
        head_edge_residuals(shell, faces, frame, error)?;
    }
    let segments: Vec<_> = expected.profile().segments().collect();
    require(faces.len() == segments.len() + 8)?;
    let shift = expected.spec().shank_length_mm;
    let (surface, forward, bounds) = face(shell, faces[0])?;
    require(forward && bounds.len() == 2)?;
    let mut previous = single_bound_rim(shell, bounds[0], true, false, frame)?;
    plane_at(proof, shell, surface, previous, frame, false)?;
    let top = previous;
    let socket_mouth = bounds[1];
    let mut decoded = Vec::new();
    let mut before_z = None;
    for (i, source) in segments.iter().enumerate() {
        require(
            previous.radius == source.start[0]
                && proof.equal(
                    previous.center,
                    frame.world([0.; 2], source.start[1] - shift),
                ),
        )?;
        let (surface, forward, bounds) = face(shell, faces[i + 1])?;
        let last = i + 1 == segments.len();
        let end = if source.role == Role::Bearing {
            require(forward && bounds.len() == 2)?;
            require(
                ring(shell, bounds[0], true)?
                    == [(previous.edge.ok_or(ReadError::Invalid)?, previous.positive)],
            )?;
            plane_at(proof, shell, surface, previous, frame, true)?;
            if last {
                require(ring(shell, bounds[1], false)? == orient(&join.positive, false))?;
                Rim {
                    edge: None,
                    vertex: join.vertex,
                    radius: thread.spec().major_radius_mm,
                    z: 0.,
                    positive: true,
                    center: frame.origin,
                }
            } else {
                single_bound_rim(shell, bounds[1], false, false, frame)?
            }
        } else {
            require(bounds.len() == 1)?;
            let (end, seam) = periodic_end(shell, bounds[0], previous, last, &join, thread, frame)?;
            match source.curve {
                HeadCurve::Line => line_surface(
                    proof,
                    shell,
                    (surface, forward, seam),
                    [previous, end],
                    source,
                    frame,
                )?,
                HeadCurve::Arc {
                    center,
                    radius_mm,
                    sweep,
                } => {
                    let arc = ArcSource {
                        center: [center[0], center[1] - shift],
                        radius: radius_mm,
                        angle: (source.start[1] - center[1]).atan2(source.start[0] - center[0]),
                        sweep,
                        role: source.role,
                    };
                    curved_surface(
                        proof,
                        shell,
                        (surface, forward, seam),
                        [previous, end],
                        arc,
                        frame,
                    )?;
                }
            }
            end
        };
        require(
            end.radius == source.end[0]
                && proof.equal(end.center, frame.world([0.; 2], source.end[1] - shift)),
        )?;
        let start_z = before_z.unwrap_or(previous.z);
        let curve = match source.curve {
            HeadCurve::Line => HeadCurve::Line,
            HeadCurve::Arc {
                center,
                radius_mm,
                sweep,
            } => HeadCurve::Arc {
                radius_mm,
                center: [
                    center[0],
                    frame.local(frame.world([0.; 2], center[1] - shift))[2],
                ],
                sweep,
            },
        };
        decoded.push(Segment {
            start: [previous.radius, start_z],
            end: [end.radius, end.z],
            curve,
            role: source.role,
        });
        before_z = Some(end.z);
        previous = end;
    }
    let (floor_surface, forward, bounds) = face(shell, faces[segments.len() + 1])?;
    require(forward && bounds.len() == 1)?;
    let floor_origin = frame.world(
        [0.; 2],
        expected.head_top_z() + expected.spec().socket_depth_mm,
    );
    let floor = proof.axis(
        shell,
        hex_head::plane(shell, floor_surface)?,
        floor_origin,
        frame.n.map(|v| -v),
        frame.u,
    )?;
    let af = expected.spec().socket_across_flats_mm;
    let points = crate::headed_screw::hex_vertices(af);
    let lower = hex_head::hex_rim(
        proof,
        shell,
        socket_mouth,
        false,
        true,
        points.map(|p| frame.world(p, expected.head_top_z())),
    )?;
    let upper = hex_head::hex_rim(
        proof,
        shell,
        bounds[0],
        true,
        false,
        points.map(|p| frame.world(p, expected.head_top_z() + expected.spec().socket_depth_mm)),
    )?;
    let area = hex_head::walls_area(
        proof,
        shell,
        &faces[segments.len() + 2..],
        [&lower, &upper],
        top.center,
        frame.n,
        true,
    )?;
    let profile =
        Profile::from_transported_segments(&decoded, expected.profile().transport_error())
            .map_err(|_| ReadError::Invalid)?;
    let depth = frame.local(floor.origin)[2] - top.z;
    require(depth > 0.)?;
    Ok(profile.volume_mm3() - area * depth)
}

fn line_surface(
    proof: Proof,
    shell: &mut Shell<'_, '_>,
    (surface, forward, seam): (u32, bool, u32),
    [from, to]: [Rim; 2],
    source: &Segment,
    frame: Frame,
) -> Result<()> {
    require(forward)?;
    let [Value::String(_), _, _, geometry, sense_value] = shell.args(seam, "EDGE_CURVE")? else {
        return Err(ReadError::Invalid);
    };
    let geometry = reference(geometry)?;
    let curve_forward = sense(sense_value)?;
    shell.line_curve(geometry, [from.vertex, to.vertex], curve_forward)?;
    let entity = shell.simple(surface)?;
    let (axis, radius, angle) = match entity.name {
        "CYLINDRICAL_SURFACE" => {
            let [Value::String(_), axis, Value::Number(radius)] = entity.args.as_slice() else {
                return Err(ReadError::Invalid);
            };
            require(from.radius == to.radius && *radius == from.radius)?;
            (reference(axis)?, *radius, 0.)
        }
        "CONICAL_SURFACE" => {
            let [
                Value::String(_),
                axis,
                Value::Number(radius),
                Value::Number(angle),
            ] = entity.args.as_slice()
            else {
                return Err(ReadError::Invalid);
            };
            require(angle.is_finite() && *angle > 0. && *angle < PI / 2.)?;
            (reference(axis)?, *radius, *angle)
        }
        _ => return Err(ReadError::Invalid),
    };
    let dr = source.end[0] - source.start[0];
    let dz = source.end[1] - source.start[1];
    let expected = if dr < 0. { to } else { from };
    let expected_radius = if dr < 0. {
        source.end[0]
    } else {
        source.start[0]
    };
    let source_angle = dr.abs().atan2(dz);
    require(
        radius == expected_radius
            && proof.error.map_or(angle == source_angle, |error| {
                (angle.tan() - source_angle.tan()).abs() * dz <= error
            }),
    )?;
    let axis = shell.placement(axis)?;
    require(
        proof.equal(axis.origin, expected.center)
            && length(sub(
                axis.normal,
                frame.n.map(|x| if dr < 0. { -x } else { x }),
            )) < proof.angular(),
    )?;
    let error = proof.error.unwrap_or(shell.tolerance / 16.);
    for rim in [from, to] {
        let delta = sub(rim.center, axis.origin);
        require(
            length(cross(delta, axis.normal)) <= error
                && (radius + dot(delta, axis.normal) * angle.tan() - rim.radius).abs() <= error,
        )?;
    }
    require(to.z > from.z)
}

struct ArcSource {
    center: [f64; 2],
    radius: f64,
    angle: f64,
    sweep: f64,
    role: Role,
}

fn curved_surface(
    proof: Proof,
    shell: &mut Shell<'_, '_>,
    (surface, forward, seam): (u32, bool, u32),
    [from, to]: [Rim; 2],
    ArcSource {
        center,
        radius,
        angle: start_angle,
        sweep,
        role,
    }: ArcSource,
    frame: Frame,
) -> Result<()> {
    let error = proof.error.unwrap_or(shell.tolerance / 16.);
    let entity = shell.simple(surface)?;
    let axis = if center[0] == 0. {
        let [Value::String(_), axis, Value::Number(r)] =
            shell.args(surface, "SPHERICAL_SURFACE")?
        else {
            return Err(ReadError::Invalid);
        };
        require(*r == radius && forward)?;
        reference(axis)?
    } else {
        let [
            Value::String(_),
            axis,
            Value::Number(major),
            Value::Number(minor),
        ] = entity.args.as_slice()
        else {
            return Err(ReadError::Invalid);
        };
        require(
            entity.name == "TOROIDAL_SURFACE"
                && *major == center[0]
                && *minor == radius
                && forward == (role == Role::Crown),
        )?;
        reference(axis)?
    };
    let support = shell.placement(axis)?;
    require(
        proof.equal(support.origin, frame.world([0.; 2], center[1]))
            && length(sub(support.normal, frame.n)) < proof.angular(),
    )?;
    let [Value::String(_), a, b, geometry, Value::Enum("T")] = shell.args(seam, "EDGE_CURVE")?
    else {
        return Err(ReadError::Invalid);
    };
    require([reference(a)?, reference(b)?] == [from.vertex, to.vertex])?;
    let [
        Value::String(_),
        circle,
        pcurves,
        Value::Enum("PCURVE_S1" | "PCURVE_S2"),
    ] = shell.args(reference(geometry)?, "SEAM_CURVE")?
    else {
        return Err(ReadError::Invalid);
    };
    let circle = reference(circle)?;
    let pcurves = refs(pcurves)?;
    require(pcurves.len() == 2 && pcurves[0] != pcurves[1])?;
    let mut values = Vec::new();
    for id in pcurves {
        let (owner, line) = shell.pcurve_definition(id)?;
        require(owner == surface)?;
        values.push(shell.line(line, true)?);
    }
    values.sort_by(|a, b| a.0[0].total_cmp(&b.0[0]));
    for ((point, direction, scale), u) in values.into_iter().zip([0., TAU]) {
        require(
            proof.error.map_or(point == [u, start_angle, 0.], |_| {
                radius * length(sub(point, [u, start_angle, 0.])) <= error
            }) && direction == [0., sweep.signum(), 0.]
                && scale == 1.,
        )?;
    }
    let [Value::String(_), axis, Value::Number(r)] = shell.args(circle, "CIRCLE")? else {
        return Err(ReadError::Invalid);
    };
    require(*r == radius)?;
    let meridian = shell.placement(reference(axis)?)?;
    let radial = support.u;
    let expected_center: Point = std::array::from_fn(|i| support.origin[i] + center[0] * radial[i]);
    let expected_u: Point = std::array::from_fn(|i| {
        start_angle.cos() * radial[i] + start_angle.sin() * support.normal[i]
    });
    let expected_n = cross(support.normal, radial).map(|v| -sweep.signum() * v);
    require(
        length(sub(meridian.origin, expected_center)) <= error
            && radius * length(sub(meridian.u, expected_u)) <= error
            && radius * length(sub(meridian.normal, expected_n)) <= error,
    )?;
    let end = start_angle + sweep;
    for (rim, angle) in [(from, start_angle), (to, end)] {
        let point: Point = std::array::from_fn(|i| {
            support.origin[i]
                + (center[0] + radius * angle.cos()) * radial[i]
                + radius * angle.sin() * support.normal[i]
        });
        require(length(sub(shell.vertex(rim.vertex)?, point)) <= error)?;
    }
    Ok(())
}

fn closed_incidence(shell: &mut Shell<'_, '_>, faces: &[u32]) -> Result<()> {
    require(faces.iter().collect::<BTreeSet<_>>().len() == faces.len())?;
    let mut directed = BTreeSet::new();
    let mut uses = BTreeMap::<u32, Vec<bool>>::new();
    for &id in faces {
        let (_, _, bounds) = face(shell, id)?;
        let mut outside = 0;
        for bound in bounds {
            let (outer, wire, _) = shell.bound(bound, false)?;
            outside += usize::from(outer);
            let [Value::String(_), oriented] = shell.args(wire, "EDGE_LOOP")? else {
                return Err(ReadError::Invalid);
            };
            for id in refs(oriented)? {
                require(directed.insert(id))?;
            }
            let edges = ring(shell, bound, outer)?;
            require(!edges.is_empty())?;
            let mut first = None;
            let mut previous = None;
            for (edge, forward) in edges {
                let [Value::String(_), a, b, _, sense_value] = shell.args(edge, "EDGE_CURVE")?
                else {
                    return Err(ReadError::Invalid);
                };
                sense(sense_value)?;
                let [a, b] = [reference(a)?, reference(b)?];
                shell.vertex(a)?;
                shell.vertex(b)?;
                let [start, end] = if forward { [a, b] } else { [b, a] };
                require(previous.is_none_or(|last| last == start))?;
                first.get_or_insert(start);
                previous = Some(end);
                uses.entry(edge).or_default().push(forward);
            }
            require(first == previous)?;
        }
        require(outside == 1)?;
    }
    require(uses.values().all(|s| s.len() == 2 && s[0] != s[1]))
}

// A source reconstructed from trims must also retain the literal carriers.
// Bound the whole finite straight edge (by convex endpoint residuals), and
// complete circular support orientation, at the source arithmetic reserve.
fn head_edge_residuals(
    shell: &mut Shell<'_, '_>,
    faces: &[u32],
    frame: Frame,
    error: f64,
) -> Result<()> {
    let mut checked = BTreeSet::new();
    for &face_id in faces {
        let (_, _, bounds) = face(shell, face_id)?;
        for bound in bounds {
            let (outer, _, _) = shell.bound(bound, false)?;
            for (edge, _) in ring(shell, bound, outer)? {
                if !checked.insert(edge) {
                    continue;
                }
                let args = shell.args(edge, "EDGE_CURVE")?;
                let geometry = reference(&args[3])?;
                if shell.simple(geometry)?.name == "SEAM_CURVE" {
                    continue;
                }
                let vertices = shell.edge(edge)?;
                let points = [shell.vertex(vertices[0])?, shell.vertex(vertices[1])?];
                match shell.edges[&edge].curve {
                    Curve::Line => {
                        let (origin, direction, _) = shell.line(geometry, false)?;
                        require(
                            points
                                .iter()
                                .all(|&p| length(cross(sub(p, origin), direction)) <= error),
                        )?;
                    }
                    Curve::Circle(circle) => {
                        require(circle.radius_mm * length(cross(circle.normal, frame.n)) <= error)?;
                        for point in points {
                            let delta = sub(point, circle.center_mm);
                            require(
                                dot(delta, circle.normal).abs() <= error
                                    && (length(delta) - circle.radius_mm).abs() <= error,
                            )?;
                        }
                    }
                    // Curved meridian seams are checked against their literal
                    // circle and both full support PCurves in curved_surface.
                    _ => {}
                }
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        head_profile::Head,
        helical_thread,
        profiled_screw::{Screw, Spec},
        step,
    };

    fn screw(crown: bool, starts: u32, left_hand: bool) -> Screw {
        let thread = Thread::new(helical_thread::Spec {
            major_radius_mm: 1.5,
            pitch_mm: 0.5,
            length_mm: 7.5,
            depth_mm: 0.2,
            crest_width_mm: 0.0625,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            starts,
            left_hand,
            root_radius_mm: None,
        })
        .unwrap();
        Screw::new(
            thread,
            Spec {
                head: if crown {
                    Head::Crowned {
                        radius_mm: 3.,
                        height_mm: 2.5,
                        crown_radius_mm: 4.,
                        top_land_radius_mm: 1.5,
                        neck_radius_mm: 0.1,
                    }
                } else {
                    Head::Countersunk {
                        radius_mm: 3.,
                        height_mm: 1.65,
                        included_angle_deg: 90.,
                        edge_chamfer_mm: 0.1,
                        neck_radius_mm: 0.1,
                    }
                },
                shank_length_mm: 1.,
                socket_across_flats_mm: 2.,
                socket_depth_mm: 0.6,
                minimum_wall_mm: 0.1,
            },
        )
        .unwrap()
    }

    fn topology(input: &str) -> Result<()> {
        compound::read_members_by(input, 1e-7, true, true, |shell, id, _| {
            let faces = shell.member_faces(id)?;
            closed_incidence(shell, &faces)
        })
        .map(|_| ())
    }

    #[test]
    fn writer_sews_actual_heads_to_every_start_without_orphan_or_duplicate_topology() {
        for crown in [false, true] {
            for (starts, hand) in [(1, false), (1, true), (2, false), (2, true)] {
                let body = screw(crown, starts, hand);
                let text =
                    step::write_profiled_screw(&body, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7)
                        .unwrap();
                assert!(text.contains("TOROIDAL_SURFACE("));
                assert_eq!(text.contains("SPHERICAL_SURFACE("), crown);
                topology(&text).unwrap();
                let volume =
                    verify_profiled_screw(&text, &body, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7)
                        .unwrap();
                assert!((volume - body.volume_mm3()).abs() < 1e-6);
                let (_, imported) = crate::step_read::read_import_solid(&text, 1e-7).unwrap();
                assert!(crate::step_read::profiled_screw_import_matches(
                    &imported,
                    &body,
                    [[0.; 3], [1., 0., 0.], [0., 1., 0.]],
                    1e-7
                ));
                if crown && starts == 1 && !hand {
                    let frame = [[0.; 3], [1., 0., 0.], [0., 1., 0.]];
                    for tolerance in [0., f64::NAN] {
                        assert!(!crate::step_read::profiled_screw_import_matches(
                            &imported, &body, frame, tolerance
                        ));
                    }
                    assert!(!crate::step_read::profiled_screw_import_matches(
                        &imported,
                        &screw(false, 1, false),
                        frame,
                        1e-7
                    ));
                    let different_thread = Screw::from_source(
                        Thread::new(helical_thread::Spec {
                            length_mm: 7.6,
                            ..body.thread().spec()
                        })
                        .unwrap(),
                        *body.profile(),
                        crate::profiled_screw::SourceSpec {
                            head: None,
                            ..body.spec()
                        },
                    )
                    .unwrap();
                    assert!(!crate::step_read::profiled_screw_import_matches(
                        &imported,
                        &different_thread,
                        frame,
                        1e-7
                    ));
                    let mut wrong = imported.clone();
                    wrong.profile = ImportProfile::HelicalThread {
                        spec: body.thread().spec(),
                        stock_radius_mm: None,
                    };
                    assert!(!crate::step_read::profiled_screw_import_matches(
                        &wrong, &body, frame, 1e-7
                    ));
                    let no_neck = Screw::new(
                        *body.thread(),
                        Spec {
                            head: Head::Crowned {
                                radius_mm: 3.,
                                height_mm: 2.5,
                                crown_radius_mm: 4.,
                                top_land_radius_mm: 1.5,
                                neck_radius_mm: 0.,
                            },
                            shank_length_mm: 0.,
                            socket_across_flats_mm: 2.,
                            socket_depth_mm: 0.6,
                            minimum_wall_mm: 0.1,
                        },
                    )
                    .unwrap();
                    assert!(!crate::step_read::profiled_screw_import_matches(
                        &imported, &no_neck, frame, 1e-7
                    ));
                }
                assert!(matches!(
                    &imported.profile,
                    ImportProfile::ProfiledScrew { .. }
                ));
                if let ImportProfile::ProfiledScrew {
                    screw,
                    head_volume_error_mm3,
                } = imported.profile
                {
                    assert!(screw.spec().head.is_none());
                    assert_eq!(screw.spec().shank_length_mm, 0.);
                    assert!((screw.volume_mm3() - body.volume_mm3()).abs() < 1e-6);
                    assert!(head_volume_error_mm3 > 0. && head_volume_error_mm3 < 1e-7);
                }
                for (before, after) in [(".T.", ".F."), (".F.", ".T.")] {
                    let line = text
                        .lines()
                        .find(|l| l.contains("ORIENTED_EDGE(") && l.contains(before))
                        .unwrap();
                    assert!(
                        topology(&text.replacen(line, &line.replace(before, after), 1)).is_err()
                    );
                }
            }
        }
    }

    #[test]
    fn physical_heads_survive_rigid_poses_without_radius_fitting() {
        let body = screw(true, 1, false);
        for origin in [
            [30., -20., 10.],
            [1000., -1000., 1000.],
            [100000., -100000., 100000.],
        ] {
            let u = [0.6, 0.8, 0.];
            let v = [0., 0., 1.];
            let text = step::write_profiled_screw(&body, origin, u, v, 1e-7).unwrap();
            let read = crate::step_read::read_import_solid(&text, 1e-7);
            let (proof, imported) = read.unwrap();
            assert!(crate::step_read::profiled_screw_import_matches(
                &imported,
                &body,
                [origin, u, v],
                1e-7
            ));
            assert!(matches!(
                &imported.profile,
                ImportProfile::ProfiledScrew { .. }
            ));
            if let ImportProfile::ProfiledScrew {
                screw,
                head_volume_error_mm3,
            } = imported.profile
            {
                assert!(screw.profile().transport_error_mm() > 0.);
                assert!(head_volume_error_mm3 > 0.);
                assert!(
                    (screw.volume_mm3() - body.volume_mm3()).abs() <= head_volume_error_mm3 + 1e-12
                );
            }
            assert!((proof.volume_mm3 - body.volume_mm3()).abs() < 1e-6);
        }
    }

    #[test]
    fn literal_head_support_changes_below_import_tolerance_are_not_refitted() {
        for crown in [false, true] {
            let old = screw(crown, 1, false);
            let body = Screw::from_source(
                Thread::new(helical_thread::Spec {
                    length_mm: 0.5,
                    ..old.thread().spec()
                })
                .unwrap(),
                *old.profile(),
                crate::profiled_screw::SourceSpec {
                    head: None,
                    ..old.spec()
                },
            )
            .unwrap();
            let text = step::write_profiled_screw(&body, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7)
                .unwrap();
            let data = parse_entities(&text, 1e-7).unwrap();
            let shell = Shell::new(&data, 1e-7, true);
            let faces = shell
                .member_faces(shell.represented_solids(true).unwrap()[0])
                .unwrap();
            let mut changes = Vec::new();
            for &id in &faces[7..faces.len() - 7] {
                let (surface, _, _) = face(&shell, id).unwrap();
                let entity = shell.simple(surface).unwrap();
                let record = match entity.args.as_slice() {
                    [Value::String(_), axis, Value::Number(radius)]
                        if entity.name == "SPHERICAL_SURFACE" =>
                    {
                        Some(format!(
                            "SPHERICAL_SURFACE('',#{},{});",
                            reference(axis).unwrap(),
                            radius + 1e-9
                        ))
                    }
                    [
                        Value::String(_),
                        axis,
                        Value::Number(major),
                        Value::Number(minor),
                    ] if entity.name == "TOROIDAL_SURFACE" => Some(format!(
                        "TOROIDAL_SURFACE('',#{},{},{});",
                        reference(axis).unwrap(),
                        major,
                        minor + 1e-9
                    )),
                    [
                        Value::String(_),
                        axis,
                        Value::Number(radius),
                        Value::Number(angle),
                    ] if entity.name == "CONICAL_SURFACE" => Some(format!(
                        "CONICAL_SURFACE('',#{},{},{});",
                        reference(axis).unwrap(),
                        radius,
                        angle + 1e-9
                    )),
                    _ => None,
                };
                if let Some(record) = record {
                    changes.push((surface, record));
                }
            }
            assert!(changes.len() >= 2);
            for (id, record) in changes {
                let key = format!("#{id}=");
                let line = text.lines().find(|l| l.starts_with(&key)).unwrap();
                let changed = text.replacen(line, &format!("{key}{record}"), 1);
                assert!(
                    crate::step_read::read_import_solid(&changed, 1e-7).is_err(),
                    "support {id} was refitted"
                );
                assert!(
                    verify_profiled_screw(
                        &changed,
                        &body,
                        [0.; 3],
                        [1., 0., 0.],
                        [0., 1., 0.],
                        1e-7
                    )
                    .is_err()
                );
            }
        }
    }

    #[test]
    fn independently_retained_neck_radius_is_literal_even_inside_numeric_reserve() {
        let old = screw(false, 1, false);
        let body = Screw::from_source(
            Thread::new(helical_thread::Spec {
                length_mm: 0.5,
                ..old.thread().spec()
            })
            .unwrap(),
            *old.profile(),
            crate::profiled_screw::SourceSpec {
                head: None,
                ..old.spec()
            },
        )
        .unwrap();
        let text =
            step::write_profiled_screw(&body, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7).unwrap();
        let data = parse_entities(&text, 1e-7).unwrap();
        let mut changed = text.clone();
        let literal = 0.1 + 1e-14;
        for (&id, entities) in &data {
            for entity in entities {
                let record = match entity.args.as_slice() {
                    [
                        Value::String(_),
                        axis,
                        Value::Number(major),
                        Value::Number(minor),
                    ] if entity.name == "TOROIDAL_SURFACE" && *minor == 0.1 => Some(format!(
                        "TOROIDAL_SURFACE('',#{},{},{literal});",
                        reference(axis).unwrap(),
                        major
                    )),
                    [Value::String(_), axis, Value::Number(radius)]
                        if entity.name == "CIRCLE" && *radius == 0.1 =>
                    {
                        Some(format!(
                            "CIRCLE('',#{},{literal});",
                            reference(axis).unwrap()
                        ))
                    }
                    _ => None,
                };
                if let Some(record) = record {
                    let key = format!("#{id}=");
                    let line = text.lines().find(|l| l.starts_with(&key)).unwrap();
                    changed = changed.replacen(line, &format!("{key}{record}"), 1);
                }
            }
        }
        let (_, imported) = crate::step_read::read_import_solid(&changed, 1e-7).unwrap();
        assert!(matches!(
            &imported.profile,
            ImportProfile::ProfiledScrew { .. }
        ));
        if let ImportProfile::ProfiledScrew { screw, .. } = imported.profile {
            let neck = screw
                .profile()
                .segments()
                .find(|s| s.role == Role::Neck)
                .unwrap();
            assert!(matches!(neck.curve,HeadCurve::Arc { radius_mm,.. } if radius_mm==literal));
        }
        assert!(
            verify_profiled_screw(&changed, &body, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7)
                .is_err()
        );
    }

    #[test]
    fn physical_terminal_bearings_and_toroidal_crowns_roundtrip() {
        let old = screw(true, 1, false);
        let thread = Thread::new(helical_thread::Spec {
            length_mm: 0.5,
            ..old.thread().spec()
        })
        .unwrap();
        let plain = Screw::new(
            thread,
            Spec {
                head: Head::Crowned {
                    radius_mm: 3.,
                    height_mm: 2.5,
                    crown_radius_mm: 4.,
                    top_land_radius_mm: 1.5,
                    neck_radius_mm: 0.,
                },
                shank_length_mm: 0.,
                socket_across_flats_mm: 1.5,
                socket_depth_mm: 0.6,
                minimum_wall_mm: 0.1,
            },
        )
        .unwrap();
        let profile = Profile::from_retained_segments(&[
            Segment {
                start: [1., -4.],
                end: [3., -2.],
                curve: HeadCurve::Arc {
                    center: [1., -2.],
                    radius_mm: 2.,
                    sweep: PI / 2.,
                },
                role: Role::Crown,
            },
            Segment {
                start: [3., -2.],
                end: [3., -0.5],
                curve: HeadCurve::Line,
                role: Role::Skirt,
            },
            Segment {
                start: [3., -0.5],
                end: [1.5, -0.5],
                curve: HeadCurve::Line,
                role: Role::Bearing,
            },
            Segment {
                start: [1.5, -0.5],
                end: [1.5, 0.],
                curve: HeadCurve::Line,
                role: Role::Shank,
            },
        ])
        .unwrap();
        let torus = Screw::from_source(
            thread,
            profile,
            crate::profiled_screw::SourceSpec {
                head: None,
                ..plain.spec()
            },
        )
        .unwrap();
        for body in [plain, torus] {
            let text = step::write_profiled_screw(&body, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7)
                .unwrap();
            let (read, _) = crate::step_read::read_import_solid(&text, 1e-7).unwrap();
            assert!((read.volume_mm3 - body.volume_mm3()).abs() < 1e-6);
        }
        let union = Screw::new_with_head_thread_union(
            thread,
            Spec {
                head: Head::Countersunk {
                    radius_mm: 3.,
                    height_mm: 1.65,
                    included_angle_deg: 90.,
                    edge_chamfer_mm: 0.1,
                    neck_radius_mm: 0.1,
                },
                shank_length_mm: 0.,
                socket_across_flats_mm: 1.5,
                socket_depth_mm: 0.6,
                minimum_wall_mm: 0.1,
            },
        )
        .unwrap();
        assert!(union.has_head_thread_overlap());
        assert!(!step::supports_profiled_screw(&union, [0.; 3], 1e-7));
        assert_eq!(
            step::write_profiled_screw(&union, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7),
            Err(step::StepError::Unsupported)
        );
    }

    #[test]
    fn malformed_head_carriers_and_physical_graph_records_are_refused() {
        fn refs_in(value: &Value<'_>, pending: &mut Vec<u32>) {
            match value {
                Value::Reference(id) => pending.push(*id),
                Value::List(values) => {
                    for value in values {
                        refs_in(value, pending);
                    }
                }
                _ => {}
            }
        }
        for crown in [false, true] {
            let old = screw(crown, 1, false);
            let body = Screw::from_source(
                Thread::new(helical_thread::Spec {
                    length_mm: 0.5,
                    ..old.thread().spec()
                })
                .unwrap(),
                *old.profile(),
                crate::profiled_screw::SourceSpec {
                    head: None,
                    ..old.spec()
                },
            )
            .unwrap();
            let text = step::write_profiled_screw(&body, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7)
                .unwrap();
            let data = parse_entities(&text, 1e-7).unwrap();
            let shell = Shell::new(&data, 1e-7, true);
            let faces = shell
                .member_faces(shell.represented_solids(true).unwrap()[0])
                .unwrap();
            let mut pending = faces[6..].to_vec();
            let mut seen = BTreeSet::new();
            while let Some(id) = pending.pop() {
                if !seen.insert(id) {
                    continue;
                }
                for entity in shell.components(id).unwrap() {
                    for arg in &entity.args {
                        refs_in(arg, &mut pending);
                    }
                    if !matches!(
                        entity.name,
                        "EDGE_CURVE"
                            | "EDGE_LOOP"
                            | "CONICAL_SURFACE"
                            | "CYLINDRICAL_SURFACE"
                            | "SPHERICAL_SURFACE"
                            | "TOROIDAL_SURFACE"
                            | "SEAM_CURVE"
                            | "CIRCLE"
                    ) {
                        continue;
                    }
                    let key = format!("#{id}=");
                    let line = text.lines().find(|l| l.starts_with(&key)).unwrap();
                    let changed = text.replacen(line, &format!("{key}{}('');", entity.name), 1);
                    assert!(
                        verify_profiled_screw(
                            &changed,
                            &body,
                            [0.; 3],
                            [1., 0., 0.],
                            [0., 1., 0.],
                            1e-7
                        )
                        .is_err(),
                        "malformed {id}"
                    );
                    assert!(
                        crate::step_read::read_import_solid(&changed, 1e-7).is_err(),
                        "malformed import {id}"
                    );
                }
            }
        }
    }

    #[test]
    fn displaced_physical_head_points_cannot_hide_in_import_tolerance() {
        fn gather(value: &Value<'_>, pending: &mut Vec<u32>) {
            match value {
                Value::Reference(id) => pending.push(*id),
                Value::List(values) => {
                    for value in values {
                        gather(value, pending);
                    }
                }
                _ => {}
            }
        }
        let old = screw(true, 1, false);
        let body = Screw::from_source(
            Thread::new(helical_thread::Spec {
                length_mm: 0.5,
                ..old.thread().spec()
            })
            .unwrap(),
            *old.profile(),
            crate::profiled_screw::SourceSpec {
                head: None,
                ..old.spec()
            },
        )
        .unwrap();
        let text =
            step::write_profiled_screw(&body, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7).unwrap();
        let data = parse_entities(&text, 1e-7).unwrap();
        let shell = Shell::new(&data, 1e-7, true);
        let faces = shell
            .member_faces(shell.represented_solids(true).unwrap()[0])
            .unwrap();
        let mut pending = faces[6..].to_vec();
        let mut seen = BTreeSet::new();
        let mut rejected = 0;
        while let Some(id) = pending.pop() {
            if !seen.insert(id) {
                continue;
            }
            for entity in shell.components(id).unwrap() {
                for value in &entity.args {
                    gather(value, &mut pending);
                }
                if entity.name != "CARTESIAN_POINT"
                    || !matches!(entity.args.as_slice(),[Value::String(_),Value::List(values)] if values.len()==3)
                {
                    continue;
                }
                let point = shell.point(id).unwrap();
                for axis in 0..3 {
                    let mut shifted = point;
                    shifted[axis] += 1e-9;
                    let key = format!("#{id}=");
                    let line = text.lines().find(|l| l.starts_with(&key)).unwrap();
                    let changed = text.replacen(
                        line,
                        &format!(
                            "{key}CARTESIAN_POINT('',({},{},{}));",
                            shifted[0], shifted[1], shifted[2]
                        ),
                        1,
                    );
                    rejected +=
                        usize::from(crate::step_read::read_import_solid(&changed, 1e-7).is_err());
                    let _ = verify_profiled_screw(
                        &changed,
                        &body,
                        [0.; 3],
                        [1., 0., 0.],
                        [0., 1., 0.],
                        1e-7,
                    );
                }
            }
        }
        assert!(rejected > 40);
    }

    #[test]
    fn rejected_head_export_rolls_back_before_retry_and_reserves_precision() {
        let body = screw(false, 1, false);
        assert!(step::supports_profiled_screw(&body, [0.; 3], 1e-7));
        for (origin, tol) in [
            ([0.; 3], 0.),
            ([0.; 3], f64::NAN),
            ([f64::INFINITY, 0., 0.], 1e-7),
            ([1e12; 3], 1e-7),
            ([0.; 3], 1e-16),
        ] {
            assert!(!step::supports_profiled_screw(&body, origin, tol));
            assert!(
                step::write_profiled_screw(&body, origin, [1., 0., 0.], [0., 1., 0.], tol).is_err()
            );
        }
        let mut writer = step::Writer::new();
        assert!(
            writer
                .append_profiled_screw(&body, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-16)
                .is_err()
        );
        writer
            .append_profiled_screw(&body, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7)
            .unwrap();
        let text = writer.finish().unwrap();
        assert_eq!(text.matches("MANIFOLD_SOLID_BREP(").count(), 1);
        crate::step_read::read_import_solid(&text, 1e-7).unwrap();
        assert!(
            step::write_profiled_screw(&body, [0.; 3], [1., 0., 0.], [0., 1.000000001, 0.], 1e-7)
                .is_err()
        );
        let parts: Vec<_> = body.profile().segments().collect();
        let loose = Profile::from_transported_segments(
            &parts,
            crate::head_profile::TransportError {
                numeric_conversion_mm: 0.005,
                arithmetic_mm: 0.,
            },
        )
        .unwrap();
        let loose = Screw::from_source(
            *body.thread(),
            loose,
            crate::profiled_screw::SourceSpec {
                head: None,
                ..body.spec()
            },
        )
        .unwrap();
        assert!(profiled_screw_head_volume_error_bound_mm3(&loose).is_none());
        assert!(!step::supports_profiled_screw(&loose, [0.; 3], 1e-7));
    }
}
