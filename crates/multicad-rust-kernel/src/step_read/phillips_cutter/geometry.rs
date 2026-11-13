//! Inspect every transmitted support and complete polynomial carrier against
//! the supplied source. No exporter controls or partition are regenerated.
use super::*;
fn axis(shell: &Shell<'_, '_>, id: u32) -> Result<LiteralAxis> {
    let [Value::String(_), origin, normal, reference_axis] =
        shell.args(id, "AXIS2_PLACEMENT_3D")?
    else {
        return Err(ReadError::Invalid);
    };
    let direction = |id| {
        let [Value::String(_), values] = shell.args(id, "DIRECTION")? else {
            return Err(ReadError::Invalid);
        };
        numbers(values)
    };
    Ok(LiteralAxis {
        origin: shell.point(reference(origin)?)?,
        normal: direction(reference(normal)?)?,
        reference: direction(reference(reference_axis)?)?,
    })
}
fn support(shell: &Shell<'_, '_>, id: u32, forward: bool) -> Result<LiteralSupport> {
    let entity = shell.simple(id)?;
    let (axis_id, kind) = match (entity.name, entity.args.as_slice()) {
        ("PLANE", [Value::String(_), a]) => (reference(a)?, Kind::Plane),
        ("CYLINDRICAL_SURFACE", [Value::String(_), a, Value::Number(r)]) => {
            (reference(a)?, Kind::Cylinder(*r))
        }
        (
            "TOROIDAL_SURFACE",
            [
                Value::String(_),
                a,
                Value::Number(major),
                Value::Number(minor),
            ],
        ) => (reference(a)?, Kind::Torus(*major, *minor)),
        ("CONICAL_SURFACE", [Value::String(_), a, Value::Number(r), Value::Number(angle)]) => {
            (reference(a)?, Kind::Cone(*r, *angle))
        }
        _ => return Err(ReadError::Unsupported),
    };
    Ok(LiteralSupport {
        axis: axis(shell, axis_id)?,
        kind,
        forward,
    })
}
fn analytic(shell: &Shell<'_, '_>, edge: graph::Edge) -> Result<Literal> {
    let entity = shell.simple(edge.geometry)?;
    let forward = edge.forward ^ edge.reversed;
    match (entity.name, entity.args.as_slice()) {
        ("LINE", [Value::String(_), origin, vector]) => {
            let [Value::String(_), direction, Value::Number(magnitude)] =
                shell.args(reference(vector)?, "VECTOR")?
            else {
                return Err(ReadError::Invalid);
            };
            require(magnitude.is_finite() && *magnitude > 0.)?;
            let [Value::String(_), values] = shell.args(reference(direction)?, "DIRECTION")? else {
                return Err(ReadError::Invalid);
            };
            Ok(Literal::Line {
                origin: shell.point(reference(origin)?)?,
                direction: numbers(values)?,
                forward,
            })
        }
        ("CIRCLE", [Value::String(_), a, Value::Number(r)]) => Ok(Literal::Circle {
            axis: axis(shell, reference(a)?)?,
            radius: *r,
            forward,
        }),
        ("ELLIPSE", [Value::String(_), a, Value::Number(x), Value::Number(y)]) => {
            Ok(Literal::Ellipse {
                axis: axis(shell, reference(a)?)?,
                semiaxes: [*x, *y],
                forward,
            })
        }
        _ => Err(ReadError::Unsupported),
    }
}
fn spline(
    shell: &Shell<'_, '_>,
    geometry: u32,
    owners: [u32; 2],
    endpoints: [u32; 2],
    work: &mut Work<'_>,
) -> Result<Vec<CubicTrimSpan>> {
    let [Value::String(_), curve, surfaces, Value::Enum("CURVE_3D")] =
        shell.args(geometry, "INTERSECTION_CURVE")?
    else {
        return Err(ReadError::Unsupported);
    };
    let actual = refs(surfaces)?;
    require(
        actual.len() == 2
            && ((actual[0] == owners[0] && actual[1] == owners[1])
                || (actual[0] == owners[1] && actual[1] == owners[0])),
    )?;
    let [
        Value::String(_),
        Value::Number(3.),
        poles,
        Value::Enum("UNSPECIFIED"),
        Value::Enum("F"),
        Value::Enum("F"),
        multiplicities,
        knots,
        Value::Enum("UNSPECIFIED"),
    ] = shell.args(reference(curve)?, "B_SPLINE_CURVE_WITH_KNOTS")?
    else {
        return Err(ReadError::Unsupported);
    };
    let Value::List(knots) = knots else {
        return Err(ReadError::Invalid);
    };
    let Value::List(multiplicities) = multiplicities else {
        return Err(ReadError::Invalid);
    };
    if knots.len() < 2 || knots.len() > 4097 {
        return Err(ReadError::Limit);
    }
    let count = knots.len() - 1;
    work.step(count.checked_mul(256).ok_or(ReadError::Limit)?)?;
    require(knots.len() == multiplicities.len())?;
    let knots = knots
        .iter()
        .map(|v| match v {
            Value::Number(x) if x.is_finite() => Ok(*x),
            _ => Err(ReadError::Invalid),
        })
        .collect::<Result<Vec<_>>>()?;
    require(knots[0] == 0. && knots[count] == 1. && knots.windows(2).all(|w| w[0] < w[1]))?;
    for (i, m) in multiplicities.iter().enumerate() {
        require(matches!(m,Value::Number(x)if *x==if i==0||i==count{4.}else{3.}))?;
    }
    let poles = refs(poles)?;
    require(
        poles.len() == 3 * count + 1
            && poles[0] == endpoints[0]
            && poles[3 * count] == endpoints[1],
    )?;
    let poles = poles
        .into_iter()
        .map(|p| shell.point(p))
        .collect::<Result<Vec<_>>>()?;
    Ok((0..count)
        .map(|i| CubicTrimSpan {
            parameters: [knots[i], knots[i + 1]],
            poles_mm: [
                poles[3 * i],
                poles[3 * i + 1],
                poles[3 * i + 2],
                poles[3 * i + 3],
            ],
        })
        .collect())
}
pub(super) fn verify(
    shell: &Shell<'_, '_>,
    graph: &graph::Graph,
    source: &BoundaryTopology<'_>,
    pose: &Pose,
    tolerance: f64,
    work: &mut Work<'_>,
) -> Result<PhillipsSourceKnownProof> {
    let source_box = source.events().source().bounding_box_mm();
    let mut surface_errors = Vec::new();
    let mut maximum_analytic = 0_f64;
    for (i, face) in graph.faces.iter().enumerate() {
        work.step(1)?;
        let actual = support(shell, face.surface, face.forward)?;
        let expected = support::Support::new(source.face_surface(i).map_err(cutter_error)?)
            .map_err(source_error)?;
        let error = expected
            .distance_bound(pose, &actual, source_box)
            .map_err(source_error)?;
        require(error <= tolerance / 16.)?;
        maximum_analytic = maximum_analytic.max(error);
        surface_errors.push(error);
    }
    let mut maximum_polynomial = 0_f64;
    let mut maximum_support = 0_f64;
    let mut certified_spans = 0;
    for (index, &edge) in graph.edges.iter().enumerate() {
        work.step(1)?;
        if let Some(expected) = source.curve_analytic_support(index).map_err(cutter_error)? {
            let actual = analytic(shell, edge)?;
            let ends = source.edges()[index]
                .vertices()
                .map(|v| graph.vertices[v].2);
            let expected_literal =
                curves::analytic_literal(expected, pose, ends).map_err(source_error)?;
            if let (Literal::Line { forward: a, .. }, Literal::Line { forward: b, .. }) =
                (&expected_literal, &actual)
            {
                require(a == b)?;
            }
            let error = curves::analytic_distance(expected, pose, &actual, source_box)
                .map_err(source_error)?;
            require(error <= tolerance / 16.)?;
            maximum_analytic = maximum_analytic.max(error);
        } else {
            require(edge.forward ^ edge.reversed)?;
            let owners = graph.owners[index].map(|i| graph.faces[i].surface);
            let endpoints = source.edges()[index]
                .vertices()
                .map(|v| graph.vertices[v].1);
            let spans = spline(shell, edge.geometry, owners, endpoints, work)?;
            let (trim, quarter, lower) =
                curves::trim_index(source, source.edges()[index].source_identity())
                    .map_err(source_error)?;
            let axes = curves::orbit_axes(pose, quarter, lower);
            curves::monotone(&spans, axes).map_err(source_error)?;
            let supporting = graph.owners[index]
                .map(|i| surface_errors[i])
                .into_iter()
                .fold(0_f64, f64::max);
            for span in &spans {
                work.step(1)?;
                let bound = source
                    .events()
                    .certify_trim_span_in_frame(trim, span, pose.origin, axes)
                    .map_err(cutter_error)?;
                let conversion = span
                    .poles_mm
                    .into_iter()
                    .map(|p| {
                        support::norm(std::array::from_fn(|i| {
                            support::adjacent(p[i]) - I::point(p[i])
                        }))
                        .hi
                    })
                    .fold(0_f64, f64::max);
                let error =
                    (I::point(bound.source_curve_error_bound_mm()) + I::point(conversion)).hi;
                require(error <= tolerance / 16.)?;
                maximum_polynomial = maximum_polynomial.max(error);
                maximum_support = maximum_support.max((I::point(error) + I::point(supporting)).hi);
                certified_spans += 1;
            }
        }
    }
    Ok(PhillipsSourceKnownProof {
        vertices: graph.vertices.len(),
        edges: graph.edges.len(),
        faces: graph.faces.len(),
        maximum_vertex_error_mm: graph.maximum_vertex_error,
        maximum_analytic_support_error_mm: maximum_analytic,
        maximum_polynomial_error_mm: maximum_polynomial,
        maximum_polynomial_support_error_mm: maximum_support,
        certified_spans,
    })
}
