//! Retain independently decoded finite carriers; never infer a head recipe.
use super::*;
use crate::profiled_screw::SourceSpec;
use crate::{head_profile::TransportError, tapered_thread::interval::I};

// Decimal STEP numbers are rounded to f64 once by the parser. One adjacent
// representable value on each side encloses that numeric conversion; this is
// independent of the caller's geometry tolerance. Directed intervals then
// enclose every subtraction, cross product, determinant and inverse projection.
struct ChartEnclosure {
    origin: [I; 3],
    rows: [[I; 3]; 3],
    determinant: I,
}
fn icross(a: [I; 3], b: [I; 3]) -> [I; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn idot(a: [I; 3], b: [I; 3]) -> I {
    a.into_iter()
        .zip(b)
        .fold(I::point(0.), |sum, (x, y)| sum + x * y)
}
impl ChartEnclosure {
    fn new(frame: Frame, conversion: bool) -> Result<Self> {
        let enclose = |v: f64| {
            if conversion {
                I::new(v.next_down(), v.next_up())
            } else {
                I::point(v)
            }
        };
        let [u, v] = [frame.u.map(enclose), frame.v.map(enclose)];
        let n = icross(u, v);
        let rows = [icross(v, n), icross(n, u), icross(u, v)];
        let determinant = idot(u, rows[0]);
        require(determinant.valid() && determinant.lo > 0.)?;
        Ok(Self {
            origin: frame.origin.map(enclose),
            rows,
            determinant,
        })
    }
    fn local(&self, point: Point, conversion: bool) -> Result<[I; 3]> {
        let point = point.map(|v| {
            if conversion {
                I::new(v.next_down(), v.next_up())
            } else {
                I::point(v)
            }
        });
        let delta = std::array::from_fn(|i| point[i] - self.origin[i]);
        let local = self.rows.map(|row| idot(row, delta) / self.determinant);
        require(local.into_iter().all(I::valid))?;
        Ok(local)
    }
}
fn references(value: &Value<'_>, out: &mut Vec<u32>) {
    match value {
        Value::Reference(id) => out.push(*id),
        Value::List(values) => {
            for value in values {
                references(value, out);
            }
        }
        _ => {}
    }
}
fn transport(shell: &Shell<'_, '_>, head: &[u32], frame: Frame) -> Result<TransportError> {
    let literal = ChartEnclosure::new(frame, false)?;
    let numeric = ChartEnclosure::new(frame, true)?;
    let mut pending = head.to_vec();
    let mut visited = BTreeSet::new();
    let mut error = TransportError::ZERO;
    while let Some(id) = pending.pop() {
        if !visited.insert(id) {
            continue;
        }
        for entity in shell.components(id)? {
            if entity.name == "CARTESIAN_POINT"
                && matches!(entity.args.as_slice(), [Value::String(_),Value::List(values)] if values.len()==3)
            {
                let point = shell.point(id)?;
                let actual = frame.local(point);
                let distance = |bounds: [I; 3]| {
                    actual.into_iter().zip(bounds).fold(0_f64, |sum, (a, i)| {
                        (sum + (a - i.lo).abs().max((a - i.hi).abs())).next_up()
                    })
                };
                error.arithmetic_mm = error
                    .arithmetic_mm
                    .max(distance(literal.local(point, false)?));
                error.numeric_conversion_mm = error
                    .numeric_conversion_mm
                    .max(distance(numeric.local(point, true)?));
            }
            for value in &entity.args {
                references(value, &mut pending);
            }
        }
    }
    require(error.total_mm().is_finite())?;
    Ok(error)
}

pub(super) fn volume_scale_error(frame: Frame) -> Result<f64> {
    let matrix = ChartEnclosure::new(frame, true)?;
    let error = (matrix.determinant - I::point(1.)).magnitude();
    require(error.is_finite())?;
    Ok(error)
}

pub(in super::super::super) fn read(
    shell: &mut Shell<'_, '_>,
    faces: &[u32],
    thread: &Thread,
    placement: [Point; 3],
) -> Result<Screw> {
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
    let join = junction(proof, shell, tail, thread, &rings, &corners, frame)?;
    let head = &tail[thread.spec().starts as usize + 1..];
    require((9..=13).contains(&head.len()))?;
    let segment_count = head.len() - 8;
    let (_, forward, bounds) = face(shell, head[0])?;
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
            };
            (end, HeadCurve::Line, Role::Bearing)
        } else {
            require(bounds.len() == 1)?;
            let (end, _) = periodic_end(shell, bounds[0], previous, last, &join, thread, frame)?;
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
    let reserve = profile.rounding_error_mm();
    let edges = ring(shell, mouth, false)?;
    require(edges.len() == 6)?;
    let mut half_width = 0_f64;
    let mut mouth_points = Vec::new();
    for (edge, _) in edges {
        for id in shell.edge(edge)? {
            let point = frame.local(shell.vertex(id)?);
            require((point[2] - top.z).abs() <= reserve)?;
            half_width = half_width.max(point[0].abs());
            mouth_points.push(point);
        }
    }
    let across = 2. * half_width;
    let regular = crate::headed_screw::hex_vertices(across);
    require(mouth_points.iter().all(|p| {
        regular
            .iter()
            .any(|q| (p[0] - q[0]).hypot(p[1] - q[1]) <= reserve)
    }))?;
    let (surface, forward, bounds) = face(shell, head[segment_count + 1])?;
    require(forward && bounds.len() == 1)?;
    let [Value::String(_), axis] = shell.args(surface, "PLANE")? else {
        return Err(ReadError::Invalid);
    };
    let floor = shell.placement(reference(axis)?)?;
    let depth = frame.local(floor.origin)[2] - top.z;
    // This is a certified physical wall, not an inferred author's minimum.
    // Keep a real dimensional lower bound for downstream feature policy;
    // using a precision epsilon here would spuriously reject a sound head.
    let corner = across / 3_f64.sqrt();
    let mut wall = 4. * reserve;
    require(profile.admits_socket(corner, depth, wall))?;
    let mut upper = profile.support([1., 0., 0.]);
    for _ in 0..52 {
        let candidate = wall + (upper - wall) * 0.5;
        if profile.admits_socket(corner, depth, candidate) {
            wall = candidate;
        } else {
            upper = candidate;
        }
    }
    Screw::from_source(
        *thread,
        profile,
        SourceSpec {
            head: None,
            shank_length_mm: 0.,
            socket_across_flats_mm: across,
            socket_depth_mm: depth,
            minimum_wall_mm: wall,
        },
    )
    .map_err(|_| ReadError::Invalid)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn transport_encloses_data_and_arithmetic_but_refuses_unresolved_inverse() {
        let frame =
            Frame::new([[100000., -100000., 100000.], [0.6, 0.8, 0.], [0., 0., 1.]]).unwrap();
        let point = [100001.2, -99998.4, 99999.5];
        let literal = ChartEnclosure::new(frame, false)
            .unwrap()
            .local(point, false)
            .unwrap();
        let numeric = ChartEnclosure::new(frame, true)
            .unwrap()
            .local(point, true)
            .unwrap();
        for ((nominal, literal), numeric) in
            frame.local(point).into_iter().zip(literal).zip(numeric)
        {
            assert!(literal.lo <= nominal && nominal <= literal.hi);
            assert!(numeric.lo <= literal.lo && numeric.hi >= literal.hi);
        }
        assert!(
            ChartEnclosure::new(frame, false)
                .unwrap()
                .local([f64::INFINITY, 0., 0.], false)
                .is_err()
        );
        let ill = Frame::new([[0.; 3], [1., 1., 0.], [1., 1. + 1e-15, 0.]]).unwrap();
        assert!(ChartEnclosure::new(ill, true).is_err());
        assert!(Frame::new([[0.; 3], [0.; 3], [0.; 3]]).is_err());
    }
}
