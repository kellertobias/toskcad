//! Actual retained outer meridian extraction; no recess recipe or mass.
use super::*;
pub(super) struct Layout {
    pub frame: Frame,
    pub count: usize,
    pub join: Junction,
    pub head: Vec<u32>,
}
pub(super) struct Meridian {
    pub profile: Profile,
    pub top: Rim,
    pub mouth: u32,
    pub outer_bound: u32,
    pub top_surface: u32,
}
pub(super) fn layout(
    shell: &mut Shell<'_, '_>,
    faces: &[u32],
    thread: &Thread,
    placement: [Point; 3],
    visible_start: f64,
) -> Result<Layout> {
    let frame = Frame::new(placement)?;
    let count = 4 * thread.spec().starts as usize;
    let mut rings = [vec![0; count], vec![0; count]];
    let mut corners = vec![[[0; 2]; 2]; count];
    for (i, &id) in faces[..count].iter().enumerate() {
        let (_, forward, bounds) = face(shell, id)?;
        require(forward == thread.spec().left_hand && bounds.len() == 1)?;
        let edges = orient(&ring(shell, bounds[0], true)?, forward);
        let [(a, true), (right, true), (b, false), (left, false)] = edges[..] else {
            return Err(ReadError::Invalid);
        };
        rings[0][i] = a;
        rings[1][i] = b;
        for (j, rail) in [left, right].into_iter().enumerate() {
            let (vertices, _, forward) = edge_record(shell, rail)?;
            require(forward)?;
            corners[i][j] = vertices;
        }
    }
    let proof = Proof {
        error: Some(shell.tolerance.min(1e-7) / 16.),
        extent: thread.spec().length_mm.max(thread.spec().major_radius_mm),
    };
    let tail = &faces[count..];
    let join = junction(
        proof,
        shell,
        tail,
        thread,
        &rings,
        &corners,
        (frame, visible_start),
    )?;
    let head = tail[join.shoulders + 1..].to_vec();
    Ok(Layout {
        frame,
        count,
        join,
        head,
    })
}
pub(super) fn meridian(
    shell: &mut Shell<'_, '_>,
    layout: &Layout,
    segment_count: usize,
) -> Result<Meridian> {
    let frame = layout.frame;
    let head = &layout.head;
    let join = &layout.join;
    let (top_surface, forward, bounds) = face(shell, head[0])?;
    let bounds_outer = *bounds.first().ok_or(ReadError::Invalid)?;
    require(forward && bounds.len() == 2)?;
    let mouth = bounds[1];
    let mut previous = single_bound_rim(shell, bounds[0], true, false, frame)?;
    let top = previous;
    let mut segments = Vec::new();
    for (i, &id) in head[1..=segment_count].iter().enumerate() {
        let (surface, forward, bounds) = face(shell, id)?;
        let name = shell.simple(surface)?.name;
        let last = i + 1 == segment_count;
        let (end, curve, role) = if name == "PLANE" {
            require(forward && bounds.len() == 2)?;
            require(
                ring(shell, bounds[0], true)?
                    == [(previous.edge.ok_or(ReadError::Invalid)?, previous.positive)],
            )?;
            let end = if last {
                require(ring(shell, bounds[1], false)? == orient(&join.positive, false))?;
                join.end
            } else {
                single_bound_rim(shell, bounds[1], false, false, frame)?
            };
            (end, HeadCurve::Line, Role::Bearing)
        } else {
            require(bounds.len() == 1)?;
            let (end, _) = periodic_end(shell, bounds[0], previous, last, join, frame)?;
            match name {
                "CYLINDRICAL_SURFACE" | "CONICAL_SURFACE" => {
                    let role = if previous.radius == end.radius {
                        if last {
                            Role::Shank
                        } else if segments.iter().any(|s: &Segment| s.role == Role::Crown) {
                            Role::Skirt
                        } else {
                            Role::Rim
                        }
                    } else if end.radius < previous.radius {
                        Role::Cone
                    } else {
                        Role::Chamfer
                    };
                    (end, HeadCurve::Line, role)
                }
                "SPHERICAL_SURFACE" | "TOROIDAL_SURFACE" => {
                    let entity = shell.simple(surface)?;
                    let (axis, major, radius, role) = match entity.args.as_slice() {
                        [Value::String(_), axis, Value::Number(radius)]
                            if name == "SPHERICAL_SURFACE" =>
                        {
                            (reference(axis)?, 0., *radius, Role::Crown)
                        }
                        [
                            Value::String(_),
                            axis,
                            Value::Number(major),
                            Value::Number(radius),
                        ] if name == "TOROIDAL_SURFACE" => (
                            reference(axis)?,
                            *major,
                            *radius,
                            if forward { Role::Crown } else { Role::Neck },
                        ),
                        _ => return Err(ReadError::Invalid),
                    };
                    let support = shell.placement(axis)?;
                    let local = frame.local(support.origin);
                    let center = [major, local[2]];
                    let a = (previous.z - center[1]).atan2(previous.radius - center[0]);
                    let b = (end.z - center[1]).atan2(end.radius - center[0]);
                    let sweep = if role == Role::Crown {
                        (b - a).rem_euclid(TAU)
                    } else {
                        -(a - b).rem_euclid(TAU)
                    };
                    (
                        end,
                        HeadCurve::Arc {
                            center,
                            radius_mm: radius,
                            sweep,
                        },
                        role,
                    )
                }
                _ => return Err(ReadError::Unsupported),
            }
        };
        segments.push(Segment {
            start: [previous.radius, previous.z],
            end: [end.radius, end.z],
            curve,
            role,
        });
        previous = end;
    }
    let error = transport(shell, head, frame)?;
    require(error.total_mm() <= shell.tolerance.min(1e-7) / 32.)?;
    let profile =
        Profile::from_transported_segments(&segments, error).map_err(|_| ReadError::Invalid)?;
    Ok(Meridian {
        profile,
        top,
        mouth,
        outer_bound: bounds_outer,
        top_surface,
    })
}
