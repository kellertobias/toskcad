//! Literal analytic supports and finite cubic carriers; no endpoint-fit lines.
use super::*;
use support::{Frame, LiteralAxis, Pose, adjacent, b, map_point, map_vector, norm};
#[derive(Clone, Debug)]
pub(crate) enum Literal {
    Line {
        origin: Point,
        direction: Point,
        forward: bool,
    },
    Circle {
        axis: LiteralAxis,
        radius: f64,
        forward: bool,
    },
    Ellipse {
        axis: LiteralAxis,
        semiaxes: [f64; 2],
        forward: bool,
    },
    Spline {
        spans: Vec<CubicTrimSpan>,
    },
}
pub(super) struct Curve {
    literal: Literal,
}
fn signed(sweep: I) -> Result<bool, Error> {
    if sweep.lo > 0. {
        Ok(true)
    } else if sweep.hi < 0. {
        Ok(false)
    } else {
        Err(Error::Limit)
    }
}
pub(crate) fn orbit_axes(pose: &Pose, quarter: u8, lower: bool) -> [Point; 3] {
    let [x, y, n] = pose.axes;
    let negative = |p: Point| p.map(|v| -v);
    let [u, v] = match quarter {
        0 => [x, y],
        1 => [y, negative(x)],
        2 => [negative(x), negative(y)],
        _ => [negative(y), x],
    };
    [u, if lower { negative(v) } else { v }, n]
}
pub(crate) fn trim_index(
    topology: &BoundaryTopology<'_>,
    key: CurveKey,
) -> Result<(usize, u8, bool), Error> {
    let CurveKey::Trim {
        envelope,
        wall,
        quarter,
        lower,
    } = key
    else {
        return Err(Error::Unsupported);
    };
    topology
        .events()
        .trims()
        .iter()
        .position(|t| t.envelope_role() == envelope && t.wall_kind() == wall)
        .map(|i| (i, quarter, lower))
        .ok_or(Error::Unsupported)
}
fn source_frame(source: AnalyticCurve) -> Result<(Frame, [I; 2], I), Error> {
    let zero = I::point(0.);
    let z = [zero, zero, I::point(1.)];
    let x = [I::point(1.), zero, zero];
    match source {
        AnalyticCurve::CircleArc {
            radius_mm_bounds,
            depth_mm_bounds,
            sweep_rad_bounds,
            ..
        } => Ok((
            Frame::new([zero, zero, b(depth_mm_bounds)], z, x)?,
            [b(radius_mm_bounds); 2],
            b(sweep_rad_bounds),
        )),
        AnalyticCurve::EllipseArc {
            center_mm_bounds,
            semiaxes_mm_bounds,
            quarter_turn,
            sweep_rad_bounds,
            ..
        } => {
            let x = match quarter_turn {
                0 => x,
                1 => [zero, I::point(1.), zero],
                2 => [I::point(-1.), zero, zero],
                _ => [zero, I::point(-1.), zero],
            };
            Ok((
                Frame::new(center_mm_bounds.map(b), z, x)?,
                semiaxes_mm_bounds.map(b),
                b(sweep_rad_bounds),
            ))
        }
        AnalyticCurve::Line(_) => Err(Error::Unsupported),
    }
}
pub(crate) fn analytic_literal(
    source: AnalyticCurve,
    pose: &Pose,
    ends: [Point; 2],
) -> Result<Literal, Error> {
    if let AnalyticCurve::Line(line) = source {
        let d = line.direction_bounds.map(b);
        let direction = pose.vector(d.map(I::midpoint));
        let delta = std::array::from_fn::<_, 3, _>(|i| I::point(ends[1][i]) - I::point(ends[0][i]));
        let forward = signed(
            delta
                .into_iter()
                .zip(map_vector(pose, d))
                .fold(I::point(0.), |s, (a, b)| s + a * b),
        )?;
        Ok(Literal::Line {
            origin: pose.point(line.origin_mm_bounds.map(b).map(I::midpoint)),
            direction,
            forward,
        })
    } else {
        let (frame, radii, sweep) = source_frame(source)?;
        let axis = support::literal_axis(pose, frame);
        let forward = signed(sweep)?;
        Ok(match source {
            AnalyticCurve::CircleArc { .. } => Literal::Circle {
                axis,
                radius: radii[0].midpoint(),
                forward,
            },
            _ => Literal::Ellipse {
                axis,
                semiaxes: radii.map(I::midpoint),
                forward,
            },
        })
    }
}
pub(crate) fn analytic_distance(
    source: AnalyticCurve,
    pose: &Pose,
    raw: &Literal,
    source_box: [[f64; 3]; 2],
) -> Result<f64, Error> {
    let delta = match (source, raw) {
        (
            AnalyticCurve::Line(line),
            Literal::Line {
                origin, direction, ..
            },
        ) => {
            let o = map_point(pose, line.origin_mm_bounds.map(b));
            let d = map_vector(pose, line.direction_bounds.map(b));
            let actual = support::unit(direction.map(adjacent))?;
            let d = std::array::from_fn::<_, 3, _>(|i| d[i] - actual[i]);
            let length = norm(std::array::from_fn(|i| {
                I::new(source_box[0][i], source_box[1][i]) - b(line.origin_mm_bounds[i])
            }));
            norm(std::array::from_fn(|i| o[i] - adjacent(origin[i]))) + length * norm(d)
        }
        (source, raw) => {
            let (frame, radii, sweep) = source_frame(source)?;
            let expected = signed(sweep)?;
            let (axis, actual, forward) = match raw {
                Literal::Circle {
                    axis,
                    radius,
                    forward,
                } => (*axis, [*radius; 2], *forward),
                Literal::Ellipse {
                    axis,
                    semiaxes,
                    forward,
                } => (*axis, *semiaxes, *forward),
                _ => return Err(Error::Unsupported),
            };
            if forward != expected || actual.into_iter().any(|r| !r.is_finite() || r <= 0.) {
                return Err(Error::Unsupported);
            }
            if matches!(source, AnalyticCurve::CircleArc { .. })
                != matches!(raw, Literal::Circle { .. })
            {
                return Err(Error::Unsupported);
            }
            let raw = Frame::literal(&axis)?;
            let o = map_point(pose, frame.origin);
            let error = norm(std::array::from_fn(|i| o[i] - raw.origin[i]));
            (0..2).fold(error, |sum, j| {
                let source = map_vector(pose, frame.axes[j].map(|v| v * radii[j]));
                sum + norm(std::array::from_fn(|i| {
                    source[i] - raw.axes[j][i] * adjacent(actual[j])
                }))
            })
        }
    };
    if !delta.valid() {
        return Err(Error::Unsupported);
    }
    Ok(delta.hi)
}
pub(crate) fn monotone(spans: &[CubicTrimSpan], axes: [Point; 3]) -> Result<(), Error> {
    for span in spans {
        for edge in span.poles_mm.windows(2) {
            let derivative = (0..3).fold(I::point(0.), |sum, i| {
                sum + (adjacent(edge[1][i]) - adjacent(edge[0][i])) * I::point(axes[1][i])
            });
            if !derivative.valid() || derivative.hi >= 0. {
                return Err(Error::Unsupported);
            }
        }
    }
    Ok(())
}
fn placed_spans(
    topology: &BoundaryTopology<'_>,
    index: usize,
    source: &[CubicTrimSpan],
    pose: &Pose,
    axes: [Point; 3],
    ends: [Point; 2],
    tolerance: f64,
) -> Result<Vec<CubicTrimSpan>, Error> {
    if source.is_empty()
        || source.windows(2).any(|w| {
            w[0].parameters[1] != w[1].parameters[0] || w[0].poles_mm[3] != w[1].poles_mm[0]
        })
    {
        return Err(Error::Unsupported);
    }
    let point = |p: Point| {
        std::array::from_fn(|i| {
            pose.origin[i] + axes[0][i] * p[0] + axes[1][i] * p[1] + axes[2][i] * p[2]
        })
    };
    let mut spans = source
        .iter()
        .map(|s| CubicTrimSpan {
            parameters: s.parameters,
            poles_mm: s.poles_mm.map(point),
        })
        .collect::<Vec<_>>();
    spans[0].poles_mm[0] = ends[0];
    let last = spans.len() - 1;
    spans[last].poles_mm[3] = ends[1];
    for span in &spans {
        let bound = topology
            .events()
            .certify_trim_span_in_frame(index, span, pose.origin, axes)?;
        let conversion = span
            .poles_mm
            .into_iter()
            .map(|p| norm(std::array::from_fn(|i| adjacent(p[i]) - I::point(p[i]))).hi)
            .fold(0_f64, f64::max);
        if (I::point(bound.source_curve_error_bound_mm()) + I::point(conversion)).hi
            > tolerance / 16.
        {
            return Err(Error::Limit);
        }
    }
    monotone(&spans, axes)?;
    Ok(spans)
}
impl Curve {
    pub fn new(
        topology: &BoundaryTopology<'_>,
        edge: usize,
        pose: &Pose,
        points: &[Point],
        spans: &[Vec<CubicTrimSpan>; 4],
        tolerance: f64,
    ) -> Result<Self, Error> {
        let ends = topology.edges()[edge].vertices().map(|i| points[i]);
        let literal = if let Some(source) = topology.curve_analytic_support(edge)? {
            let literal = analytic_literal(source, pose, ends)?;
            let error = analytic_distance(
                source,
                pose,
                &literal,
                topology.events().source().bounding_box_mm(),
            )?;
            if error > tolerance / 16. {
                return Err(Error::Limit);
            }
            literal
        } else {
            let (index, quarter, lower) =
                trim_index(topology, topology.edges()[edge].source_identity())?;
            Literal::Spline {
                spans: placed_spans(
                    topology,
                    index,
                    &spans[index],
                    pose,
                    orbit_axes(pose, quarter, lower),
                    ends,
                    tolerance,
                )?,
            }
        };
        Ok(Self { literal })
    }
    pub fn emit(
        &self,
        file: &mut StepFile,
        ends: [Vertex; 2],
        supports: [Entity; 2],
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Entity, Error> {
        let (geometry, forward) = match &self.literal {
            Literal::Line {
                origin,
                direction,
                forward,
            } => {
                let origin = file.point(*origin);
                let direction = file.direction(*direction);
                let vector = file.entity(format_args!("VECTOR('',{direction},1.)"));
                (
                    file.entity(format_args!("LINE('',{origin},{vector})")),
                    *forward,
                )
            }
            Literal::Circle {
                axis,
                radius,
                forward,
            } => {
                let axis = axis.emit(file);
                (
                    file.entity(format_args!("CIRCLE('',{axis},{radius:.16E})")),
                    *forward,
                )
            }
            Literal::Ellipse {
                axis,
                semiaxes: [a, b],
                forward,
            } => {
                let axis = axis.emit(file);
                (
                    file.entity(format_args!("ELLIPSE('',{axis},{a:.16E},{b:.16E})")),
                    *forward,
                )
            }
            Literal::Spline { spans } => {
                let mut poles = vec![ends[0].cartesian];
                for (j, span) in spans.iter().enumerate() {
                    checkpoint(interrupt)?;
                    for k in 1..4 {
                        poles.push(if j + 1 == spans.len() && k == 3 {
                            ends[1].cartesian
                        } else {
                            file.point(span.poles_mm[k])
                        });
                    }
                }
                let mut knots = vec!["0.".to_owned()];
                knots.extend(spans.iter().map(|s| format!("{:.16E}", s.parameters[1])));
                let multiplicities = (0..knots.len())
                    .map(|i| {
                        if i == 0 || i + 1 == knots.len() {
                            "4"
                        } else {
                            "3"
                        }
                    })
                    .collect::<Vec<_>>();
                let spline=file.entity(format_args!("B_SPLINE_CURVE_WITH_KNOTS('',3,{},.UNSPECIFIED.,.F.,.F.,({}),({}),.UNSPECIFIED.)",Entities(&poles),multiplicities.join(","),knots.join(",")));
                let [a, b] = supports;
                (
                    file.entity(format_args!(
                        "INTERSECTION_CURVE('',{spline},({a},{b}),.CURVE_3D.)"
                    )),
                    true,
                )
            }
        };
        let [a, b] = ends.map(|v| v.vertex);
        let sense = if forward { ".T." } else { ".F." };
        Ok(file.entity(format_args!("EDGE_CURVE('',{a},{b},{geometry},{sense})")))
    }
}
