//! Independently retained circular-wall/conical-floor cavity geometry.
use super::*;
use crate::{
    head_profile::TransportError,
    hexalobular_drive::{Arc, Cone, Drive, FiniteBinding, FinitePost, FiniteVertices},
    screw_recess::Recess,
    tapered_thread::interval::I,
};
pub(super) mod canonical;
mod graph;
use canonical::Canonical;

fn source_error(error: crate::hexalobular_drive::Error) -> ReadError {
    match error {
        crate::hexalobular_drive::Error::Invalid => ReadError::Invalid,
        crate::hexalobular_drive::Error::Unresolved => ReadError::Unsupported,
        crate::hexalobular_drive::Error::Limit | crate::hexalobular_drive::Error::Cancelled => {
            ReadError::Limit
        }
    }
}
fn adjacent(value: f64) -> I {
    I::new(value.next_down(), value.next_up())
}
fn zero(values: impl IntoIterator<Item = I>) -> Result<()> {
    require(
        values
            .into_iter()
            .all(|v| v.valid() && v.lo <= 0. && v.hi >= 0.),
    )
}
fn axis_radius(shell: &Shell<'_, '_>, id: u32, name: &str) -> Result<(u32, f64)> {
    let [Value::String(_), axis, Value::Number(radius)] = shell.args(id, name)? else {
        return Err(ReadError::Invalid);
    };
    require(radius.is_finite() && *radius > 0.)?;
    Ok((reference(axis)?, *radius))
}
fn line(
    shell: &Shell<'_, '_>,
    geometry: u32,
    endpoints: [Point; 2],
    frame: Frame,
    top: f64,
    canonical: Canonical,
    vertical: bool,
) -> Result<(FinitePost, TransportError)> {
    let [Value::String(_), origin, vector] = shell.args(geometry, "LINE")? else {
        return Err(ReadError::Invalid);
    };
    let origin = shell.point(reference(origin)?)?;
    let [Value::String(_), direction, Value::Number(magnitude)] =
        shell.args(reference(vector)?, "VECTOR")?
    else {
        return Err(ReadError::Invalid);
    };
    require(magnitude.is_finite() && *magnitude > 0.)?;
    let direction_id = reference(direction)?;
    let [Value::String(_), values] = shell.args(direction_id, "DIRECTION")? else {
        return Err(ReadError::Invalid);
    };
    let raw: Point = numbers(values)?;
    let scale = raw.into_iter().map(f64::abs).fold(0., f64::max);
    require(scale.is_finite() && scale > 0.)?;
    let bounded = raw.map(|v| adjacent(v) / I::point(scale));
    let norm = bounded
        .into_iter()
        .fold(I::point(0.), |s, v| s + v.square())
        .sqrt();
    require(norm.valid() && norm.lo > 0.)?;
    let bounded = bounded.map(|v| v / norm * adjacent(*magnitude));
    let local = canonical.project_vector(bounded, true)?;
    // The represented family has vertical generators. Real tilted support
    // lines cannot be admitted by spending the file's requested tolerance.
    if vertical {
        zero(local[..2].iter().copied())?;
    }
    let vector = shell.direction(direction_id)?.map(|v| v * magnitude);
    let scalar_local = canonical.local_vector(vector, true)?;
    zero((0..3).map(|i| scalar_local[i] - local[i]))?;
    let index = (0..3)
        .max_by(|&a, &b| vector[a].abs().total_cmp(&vector[b].abs()))
        .ok_or(ReadError::Invalid)?;
    require(vector[index] != 0.)?;
    let parameters = endpoints.map(|p| (p[index] - origin[index]) / vector[index]);
    for point in endpoints {
        let parameter = (adjacent(point[index]) - adjacent(origin[index])) / bounded[index];
        zero((0..3).map(|i| adjacent(origin[i]) + bounded[i] * parameter - adjacent(point[i])))?;
    }
    let mut origin_mm = frame.local(origin);
    origin_mm[2] -= top;
    let direction_mm_per_parameter = frame.rows.map(|row| dot(row, vector) / frame.determinant);
    let distance = |bounds: [I; 3]| {
        direction_mm_per_parameter
            .into_iter()
            .zip(bounds)
            .fold(I::point(0.), |sum, (v, b)| {
                sum + I::point((I::point(v) - b).magnitude())
            })
            .hi
    };
    let literal = raw.map(|v| I::point(v) / I::point(scale));
    let norm = literal
        .into_iter()
        .fold(I::point(0.), |s, v| s + v.square())
        .sqrt();
    require(norm.valid() && norm.lo > 0.)?;
    let literal =
        canonical.project_vector(literal.map(|v| v / norm * I::point(*magnitude)), false)?;
    let extent = I::point(parameters[0].abs().max(parameters[1].abs()));
    let error = TransportError {
        arithmetic_mm: (I::point(distance(literal)) * extent).hi,
        numeric_conversion_mm: (I::point(distance(local)) * extent).hi,
    };
    Ok((
        FinitePost {
            origin_mm,
            direction_mm_per_parameter,
            parameters,
        },
        error,
    ))
}

pub(super) fn read(
    shell: &mut Shell<'_, '_>,
    all_faces: &[u32],
    head: &[u32],
    segment_count: usize,
    frame: Frame,
    count: usize,
    top: Rim,
) -> Result<Recess> {
    require(head.len() == segment_count + 14)?;
    let (surface, _, bounds) = face(shell, head[0])?;
    require(bounds.len() == 2)?;
    let [Value::String(_), top_axis] = shell.args(surface, "PLANE")? else {
        return Err(ReadError::Invalid);
    };
    let top_axis = reference(top_axis)?;
    let canonical = Canonical::new(shell, all_faces, frame, count)?;
    let graph = graph::read(
        shell,
        bounds[1],
        head[segment_count + 1],
        &head[segment_count + 2..],
    )?;
    let [
        Value::String(_),
        cone_axis,
        Value::Number(0.),
        Value::Number(angle),
    ] = shell.args(graph.cone, "CONICAL_SURFACE")?
    else {
        return Err(ReadError::Unsupported);
    };
    let (cone_axis, angle) = (reference(cone_axis)?, *angle);
    require(angle > 0. && angle < PI / 2.)?;
    canonical.coaxial_with_thread(shell, cone_axis)?;
    require(dot(shell.placement(cone_axis)?.normal, frame.n) < 0.)?;
    let mut error = canonical.transport(shell, head)?;
    // Every cavity coordinate subtracts the independently transported top
    // ordinate. Summing both point bounds encloses that complete operation.
    error.arithmetic_mm = (2. * error.arithmetic_mm).next_up();
    error.numeric_conversion_mm = (2. * error.numeric_conversion_mm).next_up();
    let relative = |point| {
        let mut p = frame.local(point);
        p[2] -= top.z;
        p
    };
    let apex = finite_graph::vertex(shell, graph.apex)?;
    require(apex == shell.placement(cone_axis)?.origin)?;
    let apex_mm = relative(apex);
    let cone = Cone {
        tip_depth_mm: apex_mm[2],
        depth_per_radius: angle.portable_cos() / angle.portable_sin(),
    };
    let mut arcs = Vec::new();
    let mut vertices = FiniteVertices {
        mouth_mm: [[0.; 3]; 12],
        floor_mm: [[0.; 3]; 12],
        apex_mm,
    };
    let mut posts = Vec::new();
    let mut line_error = TransportError::ZERO;
    let mut carriers = Vec::new();
    for (i, cell) in graph.cells.iter().enumerate() {
        let (circle_axis, radius) = axis_radius(shell, cell.circle, "CIRCLE")?;
        let (cylinder_axis, cylinder_radius) =
            axis_radius(shell, cell.cylinder, "CYLINDRICAL_SURFACE")?;
        require(radius == cylinder_radius)?;
        canonical.coplanar_origins(shell, circle_axis, top_axis)?;
        canonical.same_axis_line(shell, cylinder_axis, circle_axis)?;
        let circle = shell.placement(circle_axis)?;
        let center = relative(circle.origin);
        let points = [shell.vertex(cell.mouth[0])?, shell.vertex(cell.mouth[1])?];
        let local = points.map(relative);
        let mut endpoints = [[0.; 2]; 2];
        for j in 0..2 {
            let p = canonical.local_point(points[j], true)?;
            let c = canonical.local_point(circle.origin, true)?;
            let square = (p[0] - c[0]).square() + (p[1] - c[1]).square();
            zero([square - adjacent(radius).square()])?;
            endpoints[j] = [local[j][0], local[j][1]];
        }
        let start = (endpoints[0][1] - center[1]).portable_atan2(endpoints[0][0] - center[0]);
        let end = (endpoints[1][1] - center[1]).portable_atan2(endpoints[1][0] - center[0]);
        let positive = (dot(circle.normal, frame.n) > 0.) == cell.circle_forward;
        let sweep = if positive {
            (end - start).rem_euclid(TAU)
        } else {
            -(start - end).rem_euclid(TAU)
        };
        require(cell.wall_forward == (sweep < 0.))?;
        arcs.push(Arc {
            center_mm: [center[0], center[1]],
            radius_mm: radius,
            start_mm: endpoints[0],
            end_mm: endpoints[1],
            start_angle_rad: start,
            sweep_rad: sweep,
        });
        let floor = [shell.vertex(cell.floor[0])?, shell.vertex(cell.floor[1])?];
        vertices.mouth_mm[i] = local[0];
        vertices.floor_mm[i] = relative(floor[0]);
        let (post, post_error) = line(
            shell,
            cell.post,
            [points[0], floor[0]],
            frame,
            top.z,
            canonical,
            true,
        )?;
        require((post.parameters[1] > post.parameters[0]) == cell.post_forward)?;
        line_error.arithmetic_mm = line_error.arithmetic_mm.max(post_error.arithmetic_mm);
        line_error.numeric_conversion_mm = line_error
            .numeric_conversion_mm
            .max(post_error.numeric_conversion_mm);
        posts.push(post);
        let mut spans = graph::spline(shell, cell.spline)?;
        require(
            spans.first().ok_or(ReadError::Invalid)?.poles_mm[0] == floor[0]
                && spans.last().ok_or(ReadError::Invalid)?.poles_mm[3] == floor[1],
        )?;
        for span in &mut spans {
            span.poles_mm = span.poles_mm.map(relative);
        }
        carriers.push(spans);
    }
    // Canonical cyclic start stabilizes identities without changing any point,
    // support, knot or carrier. Source construction requires positive even arcs.
    let first = (0..12)
        .filter(|&i| arcs[i].sweep_rad > 0.)
        .max_by(|&i, &j| {
            arcs[i].center_mm[0]
                .total_cmp(&arcs[j].center_mm[0])
                .then(arcs[i].center_mm[1].total_cmp(&arcs[j].center_mm[1]))
        })
        .ok_or(ReadError::Invalid)?;
    arcs.rotate_left(first);
    posts.rotate_left(first);
    carriers.rotate_left(first);
    vertices.mouth_mm.rotate_left(first);
    vertices.floor_mm.rotate_left(first);
    let (seam_line, distance) = graph::seam(shell, &graph)?;
    let (seam_vertices, _, _) = edge_record(shell, graph.seam)?;
    let start = shell.vertex(seam_vertices[0])?;
    let (_, seam_error) = line(
        shell,
        seam_line,
        [start, apex],
        frame,
        top.z,
        canonical,
        false,
    )?;
    error.arithmetic_mm =
        (error.arithmetic_mm + line_error.arithmetic_mm.max(seam_error.arithmetic_mm)).next_up();
    error.numeric_conversion_mm = (error.numeric_conversion_mm
        + line_error
            .numeric_conversion_mm
            .max(seam_error.numeric_conversion_mm))
    .next_up();
    let axis = canonical.local_axis(shell, cone_axis, true)?;
    let raw_apex = canonical::Axis::read(shell, cone_axis, true)?.origin;
    let apex_bounds = canonical.project_point(raw_apex, true)?;
    zero((0..3).map(|i| apex_bounds[i] - axis.origin[i]))?;
    let start_i = canonical.local_point(start, true)?;
    let delta: [I; 3] = std::array::from_fn(|i| start_i[i] - axis.origin[i]);
    let along = decode::idot(delta, axis.axes[2]);
    let radial = decode::idot(delta, axis.axes[0]);
    zero([decode::idot(delta, axis.axes[1])])?;
    require(radial.lo > 0. && along.lo > 0.)?;
    zero([
        radial * adjacent(angle).cos() - along * adjacent(angle).sin(),
        decode::idot(delta, delta) - adjacent(distance).square(),
    ])?;
    let source = Drive::from_retained(
        arcs.try_into().map_err(|_| ReadError::Invalid)?,
        cone,
        error,
    )
    .map_err(source_error)?;
    let binding = FiniteBinding::new_with_posts(
        source,
        vertices,
        posts.try_into().map_err(|_| ReadError::Invalid)?,
        carriers.try_into().map_err(|_| ReadError::Invalid)?,
        shell.tolerance.min(1e-7) / 8.,
    )
    .map_err(source_error)?
    .with_literal_cone_semi_angle(angle)
    .map_err(source_error)?;
    let d = (16. * binding.source().source_error_mm()).next_up();
    let floor = (d * (1. + cone.depth_per_radius)).next_up();
    let source_volume = binding
        .source()
        .boundary_volume_error_bound_mm3(d, floor)
        .ok_or(ReadError::Unsupported)?;
    let trim_volume = binding
        .additional_volume_error_bound_mm3()
        .ok_or(ReadError::Unsupported)?;
    let width = (4. * (source_volume + trim_volume)
        + 1e-8 * (binding.source().maximum_radius_bound_mm().powi(2) + cone.tip_depth_mm.powi(2)))
    .next_up();
    Recess::conical_hexalobular_from_binding(binding, width, 100_000).map_err(source_error)
}

/// Compare finite cavity support geometry in two physical charts. The caller
/// separately proves every serialized carrier and its complete incidence.
pub(in crate::step_read) fn matches(
    a: &Recess,
    b: &Recess,
    frames: [[Point; 3]; 2],
    tops: [f64; 2],
    tolerance: f64,
) -> bool {
    let (Some(a), Some(b)) = (a.hexalobular_geometry(), b.hexalobular_geometry()) else {
        return false;
    };
    let place = |index: usize, p: Point| {
        let [origin, u, v] = frames[index];
        let n = cross(u, v);
        std::array::from_fn(|i| origin[i] + p[0] * u[i] + p[1] * v[i] + (tops[index] + p[2]) * n[i])
    };
    let near = |a: Point, b: Point| length(sub(a, b)) <= tolerance;
    let reserve = a.spatial_source_error_bound_mm() + b.spatial_source_error_bound_mm();
    if !tolerance.is_finite() || tolerance <= 0. || reserve > tolerance / 4. {
        return false;
    }
    let radius = a.maximum_radius_bound_mm().max(b.maximum_radius_bound_mm());
    if (a.cone().depth_per_radius - b.cone().depth_per_radius).abs() * radius + reserve
        > tolerance / 2.
        || !near(
            place(0, [0., 0., a.cone().tip_depth_mm]),
            place(1, [0., 0., b.cone().tip_depth_mm]),
        )
    {
        return false;
    }
    (0..12).any(|shift| {
        a.arcs().iter().enumerate().all(|(i, a)| {
            let b = b.arcs()[(i + shift) % 12];
            let point = |p: [f64; 2]| [p[0], p[1], 0.];
            a.radius_mm == b.radius_mm
                && (a.sweep_rad - b.sweep_rad).abs() * a.radius_mm + reserve <= tolerance / 2.
                && near(place(0, point(a.center_mm)), place(1, point(b.center_mm)))
                && near(place(0, point(a.start_mm)), place(1, point(b.start_mm)))
                && near(place(0, point(a.end_mm)), place(1, point(b.end_mm)))
        })
    })
}
