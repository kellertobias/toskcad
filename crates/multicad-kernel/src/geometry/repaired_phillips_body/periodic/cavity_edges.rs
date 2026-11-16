//! Original finite Cutter curves under the actual periodic material map.
use super::*;
use multicad_rust_kernel::phillips_drive::{
    AnalyticCurve, BoundaryTopology, CubicTrimSpan, CurveKey,
};
type B = I;
type Place<'a> = dyn Fn([[f64; 2]; 3], Option<Point>) -> Result<(Point, f64), GeometryError> + 'a;
fn cavity_failure(e: multicad_rust_kernel::phillips_drive::Error) -> GeometryError {
    use multicad_rust_kernel::phillips_drive::Error as E;
    let kind = match e {
        E::Invalid => GeometryErrorKind::Input,
        E::Limit => GeometryErrorKind::Limit,
        E::Cancelled => GeometryErrorKind::Cancelled,
        E::Unresolved => GeometryErrorKind::Unsupported,
    };
    failure(kind, format!("Actual Cutter curve refused: {e:?}"))
}
fn unresolved(message: &str) -> GeometryError {
    error(
        "Rust Phillips display",
        GeometryErrorKind::Unsupported,
        message,
    )
}
fn limit() -> GeometryError {
    error(
        "Rust Phillips display",
        GeometryErrorKind::Limit,
        "Finite source edges exceed the shared display point budget",
    )
}
fn polynomial_second(span: &CubicTrimSpan) -> f64 {
    (0..2)
        .map(|k| {
            (0..3)
                .fold(B::ZERO, |sum, i| {
                    let value = B::point(6.)
                        * (B::point(span.poles_mm[k][i])
                            - B::point(2.) * B::point(span.poles_mm[k + 1][i])
                            + B::point(span.poles_mm[k + 2][i]));
                    sum + B::point(value.magnitude())
                })
                .hi
        })
        .fold(0., f64::max)
}
struct Span {
    domain: [f64; 2],
    second: f64,
    source_error: f64,
}
struct Query<'a, 'b> {
    topology: &'a BoundaryTopology<'b>,
    place: &'a Place<'a>,
    deflection: f64,
    operator: f64,
    endpoints: &'a [[Point; 2]],
}
impl Query<'_, '_> {
    fn point(
        &self,
        edge: usize,
        t: f64,
        existing: Option<Point>,
    ) -> Result<(Point, f64), GeometryError> {
        let point = self
            .topology
            .curve_point_mm_bounds(edge, t)
            .map_err(cavity_failure)?;
        let existing = existing.or(if t == 0. {
            Some(self.endpoints[edge][0])
        } else if t == 1. {
            Some(self.endpoints[edge][1])
        } else {
            None
        });
        let (point, error) = (self.place)(point, existing)?;
        if error > self.deflection / 4. {
            return Err(unresolved(
                "Literal curve placement exceeds the requested display accuracy",
            ));
        }
        Ok((point, error))
    }
    fn append(
        &self,
        edge: usize,
        span: Span,
        output: &mut Vec<Point>,
        remaining: &mut usize,
    ) -> Result<(), GeometryError> {
        let local = (self.deflection / (B::point(4.) * B::point(self.operator)).hi).next_down();
        let capacity = local - 2. * span.source_error;
        if !capacity.is_finite() || capacity <= 0. || !span.second.is_finite() {
            return Err(unresolved(
                "Finite source curve reaches the numeric chord floor",
            ));
        }
        // Reserve half of the local chord budget for arithmetic and actual
        // rounded parameter gaps, checked again below on every emitted chord.
        let count = (B::point(span.second) / (B::point(4.) * B::point(capacity)))
            .hi
            .sqrt()
            .next_up()
            .ceil()
            .max(1.);
        let shared = usize::from(!output.is_empty());
        if !count.is_finite() || count + 1. - shared as f64 > *remaining as f64 {
            return Err(limit());
        }
        let count = count as usize;
        let extent = B::point(span.domain[1]) - B::point(span.domain[0]);
        if !extent.valid() || extent.lo <= 0. {
            return Err(unresolved(
                "Finite source curve parameter span is unresolved",
            ));
        }
        let mut previous = None;
        for i in 0..=count {
            let t = if i == 0 {
                span.domain[0]
            } else if i == count {
                span.domain[1]
            } else {
                (span.domain[1] - span.domain[0]).mul_add(i as f64 / count as f64, span.domain[0])
            };
            let (point, error) =
                self.point(edge, t, (i == 0).then(|| output.last().copied()).flatten())?;
            if let Some((before, before_error)) = previous {
                let gap = (B::point(t) - B::point(before)) / extent;
                let chord = B::point(2.) * B::point(span.source_error)
                    + B::point(span.second) * gap.square() / B::point(8.);
                let world = B::point(self.operator) * chord + B::point(error.max(before_error));
                if !world.valid() || world.hi > self.deflection {
                    return Err(unresolved(
                        "Literal finite curve chord exceeds requested accuracy",
                    ));
                }
            }
            if i != 0 || output.is_empty() {
                output.push(point);
                *remaining -= 1;
            }
            previous = Some((t, error));
        }
        Ok(())
    }
}
fn analytic_span(curve: AnalyticCurve) -> Span {
    let second = match curve {
        AnalyticCurve::Line(_) => 0.,
        AnalyticCurve::CircleArc {
            radius_mm_bounds: [_, r],
            sweep_rad_bounds: [a, b],
            ..
        } => (B::point(r) * B::new(a, b).square()).hi,
        AnalyticCurve::EllipseArc {
            semiaxes_mm_bounds: [a, b],
            sweep_rad_bounds: [lo, hi],
            ..
        } => (B::point(a[1].max(b[1])) * B::new(lo, hi).square()).hi,
    };
    Span {
        domain: [0., 1.],
        second,
        source_error: 0.,
    }
}
pub(super) fn edges(
    cutter: &BoundaryTopology<'_>,
    ordinals: &[u32],
    endpoints: &[[Point; 2]],
    deflection: f64,
    operator: f64,
    remaining: &mut usize,
    place: &Place<'_>,
) -> Result<Vec<DisplayEdge>, GeometryError> {
    if ordinals.len() != cutter.edges().len()
        || endpoints.len() != cutter.edges().len()
        || !deflection.is_finite()
        || deflection <= 0.
        || !operator.is_finite()
        || operator <= 0.
    {
        return Err(input("Invalid actual Cutter edge display request"));
    }
    let query = Query {
        topology: cutter,
        place,
        deflection,
        operator,
        endpoints,
    };
    let local = (deflection / (B::point(4.) * B::point(operator)).hi).next_down();
    if !local.is_finite() || local <= 0. {
        return Err(unresolved(
            "Phillips edge display has no finite source budget",
        ));
    }
    let mut result = Vec::with_capacity(cutter.edges().len());
    for (index, edge) in cutter.edges().iter().enumerate() {
        if *remaining < 2 {
            return Err(limit());
        }
        // Refuse impossible precision before preparing any expensive carrier.
        for t in [0., 0.5, 1.] {
            query.point(index, t, None)?;
        }
        let ordinal = ordinals[index];
        let mut points = Vec::new();
        if let CurveKey::Trim { envelope, wall, .. } = edge.source_identity() {
            let trim = cutter
                .events()
                .trims()
                .iter()
                .position(|t| t.envelope_role() == envelope && t.wall_kind() == wall)
                .expect("retained source trim role");
            let spans = cutter
                .events()
                .trim_carrier(trim, local / 16., (*remaining).min(100_000))
                .map_err(cavity_failure)?;
            for span in spans {
                let error = cutter
                    .events()
                    .certify_trim_span(trim, &span)
                    .map_err(cavity_failure)?
                    .source_curve_error_bound_mm();
                query.append(
                    index,
                    Span {
                        domain: span.parameters,
                        second: polynomial_second(&span),
                        source_error: error,
                    },
                    &mut points,
                    remaining,
                )?;
            }
        } else {
            // The sealed source returns None only for dual-support Trim keys.
            // Every other finite edge retains its actual analytic support.
            let curve = cutter
                .curve_analytic_support(index)
                .map_err(cavity_failure)?
                .expect("nonintersection source edge has retained analytic support");
            query.append(index, analytic_span(curve), &mut points, remaining)?;
        }
        result.push(DisplayEdge {
            index: ordinal,
            points_mm: points,
        });
    }
    Ok(result)
}
