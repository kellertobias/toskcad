//! Original finite master endpoints to the retained reference parameter window.
//! Curve-coefficient transfer alone cannot account for rounded VERTEX_POINTs.
use super::*;
type V = [I; 3];
fn dot(a: V, b: V) -> I {
    (0..3).fold(I::point(0.), |s, i| s + a[i] * b[i])
}
fn vertex(source: &Source<'_>, id: usize, budget: &mut Budget<'_>) -> Result<V> {
    let id = u32::try_from(id).map_err(|_| ReadError::Unsupported)?;
    let [_, point] = source.args(id, "VERTEX_POINT")? else {
        return Err(ReadError::Invalid);
    };
    let mut point = reference(point)?;
    if source.simple(point)?.name == "EVALUATED_DEGENERATE_PCURVE" {
        // Earn the complete original collapsed trace before using its same
        // original location owner. The cached equivalent point is not proof.
        finite_analytic::equivalent_vertex(source, point, budget)?.ok_or(ReadError::Unsupported)?;
        let [_, _, _, equivalent] = source.args(point, "EVALUATED_DEGENERATE_PCURVE")? else {
            return Err(ReadError::Invalid);
        };
        point = reference(equivalent)?;
    }
    precision::point(source, point, 3, budget)
}
pub(super) fn reserve(
    source: &Source<'_>,
    curve_id: u32,
    edge: &Edge,
    budget: &mut Budget<'_>,
) -> Result<f64> {
    // Two endpoint projections/rotated rays, fixed vector folds, finite-span
    // checks/copies and trig interval calls: <1024 local logical actions.
    // Raw numeric/frame leaves have their own unchanged delegated debits.
    budget.charge(4096)?;
    if edge.range.into_iter().any(|t| !t.is_finite()) || edge.range[0] >= edge.range[1] {
        return Err(ReadError::Invalid);
    }
    let id = precision::unwrapped(source, curve_id, budget)?;
    let entity = source.simple(id)?;
    let result = match (entity.name, entity.args.as_slice()) {
        ("LINE", [_, origin, direction]) => {
            let origin = precision::point(source, reference(origin)?, 3, budget)?;
            let d = precision::vector(source, reference(direction)?, 3, budget)?;
            let denominator = dot(d, d);
            if !denominator.valid() || denominator.lo <= 0. {
                return Err(ReadError::Unsupported);
            }
            let mut parameters = [I::point(0.); 2];
            for (i, p) in parameters.iter_mut().enumerate() {
                *p = dot(
                    difference(vertex(source, edge.vertices[i], budget)?, origin),
                    d,
                ) / denominator;
            }
            if parameters.into_iter().any(|p| !p.valid()) || parameters[1].lo <= parameters[0].hi {
                return Err(ReadError::Unsupported);
            }
            let differences = [
                parameters[0] - I::point(edge.range[0]),
                parameters[1] - I::point(edge.range[1]),
            ];
            if differences.into_iter().any(|d| !d.valid()) {
                return Err(ReadError::Unsupported);
            }
            let delta = differences[0].magnitude().max(differences[1].magnitude());
            I::point(delta) * norm(d)
        }
        (name, [_, axis, rest @ ..]) if matches!(name, "CIRCLE" | "ELLIPSE") => {
            let (a, b) = match (name, rest) {
                ("CIRCLE", [r]) => {
                    let r = precision::scalar(source, id, r, budget)?;
                    (r, r)
                }
                ("ELLIPSE", [a, b]) => (
                    precision::scalar(source, id, a, budget)?,
                    precision::scalar(source, id, b, budget)?,
                ),
                _ => return Err(ReadError::Invalid),
            };
            if !a.valid() || !b.valid() || a.lo <= 0. || b.lo <= 0. {
                return Err(ReadError::Unsupported);
            }
            let frame = precision::frame(source, reference(axis)?, 3, budget)?;
            let mut delta = [0.; 2];
            let closed = edge.vertices[0] == edge.vertices[1];
            for (i, item) in delta.iter_mut().enumerate() {
                if closed && i == 1 {
                    continue;
                }
                let p = difference(vertex(source, edge.vertices[i], budget)?, frame.origin);
                let x = dot(p, frame.axes[0]) / a;
                let y = dot(p, frame.axes[1]) / b;
                let theta = I::point(edge.range[i]);
                let (c, s) = (theta.cos(), theta.sin());
                let forward = x * c + y * s;
                let sideways = y * c - x * s;
                if !forward.valid() || !sideways.valid() || forward.lo <= 0. {
                    return Err(ReadError::Unsupported);
                }
                // The actual angle difference is in (-pi/2,pi/2). Since
                // |atan(t)|<=|t|, no unearned floating atan2 inversion enters.
                let angle = I::point(sideways.magnitude()) / I::point(forward.lo);
                if !angle.valid() {
                    return Err(ReadError::Unsupported);
                }
                *item = angle.hi;
            }
            // Independently enclosed true 2*pi; the exact Machin-series test
            // establishes these adjacent binary bounds (not TAU equality).
            let tau = I::new(
                std::f64::consts::TAU.next_down(),
                std::f64::consts::TAU.next_up(),
            );
            if closed {
                if edge.range[1] != edge.range[0] + std::f64::consts::TAU {
                    return Err(ReadError::Unsupported);
                }
                let difference = I::point(edge.range[0]) + tau - I::point(edge.range[1]);
                if !difference.valid() {
                    return Err(ReadError::Unsupported);
                }
                let span_error = difference.magnitude();
                let bound = I::point(delta[0]) + I::point(span_error);
                if !bound.valid() {
                    return Err(ReadError::Unsupported);
                }
                delta[1] = bound.hi;
            } else {
                let span = I::point(edge.range[1]) - I::point(edge.range[0]);
                let uncertainty = I::point(delta[0]) + I::point(delta[1]);
                let lower = span - uncertainty;
                let upper = span + uncertainty;
                if !lower.valid() || !upper.valid() || lower.lo <= 0. || upper.hi >= tau.lo {
                    return Err(ReadError::Unsupported);
                }
            }
            I::point(delta[0].max(delta[1])) * (a + b)
        }
        _ => return Err(ReadError::Unsupported),
    };
    if !result.valid() || result.hi < 0. {
        return Err(ReadError::Unsupported);
    }
    budget.charge(64)?;
    Ok(result.hi)
}
