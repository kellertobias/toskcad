//! Finite original-REAL analytic traces. No bounded UV extent is erased.
use super::*;
use crate::tapered_thread::interval::I;
use num_rational::BigRational as Q;
type V = [I; 3];
fn add(a: V, b: V) -> V {
    std::array::from_fn(|i| a[i] + b[i])
}
fn scale(a: V, b: I) -> V {
    a.map(|x| x * b)
}
fn dot(a: V, b: V) -> I {
    (0..3).fold(I::point(0.), |s, i| s + a[i] * b[i])
}
fn trig(a: I, cosine: bool) -> Result<I> {
    if !a.valid() {
        return Err(ReadError::Unsupported);
    }
    let m = a.midpoint();
    if !m.is_finite() {
        return Err(ReadError::Unsupported);
    }
    let d = (a - I::point(m)).magnitude();
    let v = if cosine {
        I::point(m).cos()
    } else {
        I::point(m).sin()
    };
    let v = I::new(
        (v.lo - d).next_down().max(-1.),
        (v.hi + d).next_up().min(1.),
    );
    if v.valid() {
        Ok(v)
    } else {
        Err(ReadError::Unsupported)
    }
}
fn finish(proof: Proof, target_error: f64, tol: f64) -> Proof {
    if !target_error.is_finite() || target_error < 0. {
        return Proof::Unsupported;
    }
    match proof {
        Proof::Certified(upper) => {
            let value = I::point(upper) + I::point(target_error);
            if value.valid() && value.hi <= tol {
                Proof::Certified(value.hi)
            } else {
                Proof::Unsupported
            }
        }
        Proof::Disagreement(lower) => {
            let value = I::point(lower) - I::point(target_error);
            if value.valid() && value.lo > tol {
                Proof::Disagreement(value.lo)
            } else {
                Proof::Unsupported
            }
        }
        Proof::Unsupported => Proof::Unsupported,
    }
}
fn coordinate<'a, 'i>(source: &'a Source<'i>, point_id: u32, axis: usize) -> Result<&'a Value<'i>> {
    let [_, values] = source.args(point_id, "CARTESIAN_POINT")? else {
        return Err(ReadError::Invalid);
    };
    let values = list(values)?;
    if values.len() != 2 {
        return Err(ReadError::Invalid);
    }
    values.get(axis).ok_or(ReadError::Invalid)
}
fn polyline(source: &Source<'_>, id: u32, budget: &mut Budget<'_>) -> Result<Option<[u32; 2]>> {
    budget.charge(256)?;
    if source.simple(id)?.name != "POLYLINE" {
        return Ok(None);
    }
    let [_, points] = source.args(id, "POLYLINE")? else {
        return Err(ReadError::Invalid);
    };
    let [a, b] = list(points)? else {
        return Err(ReadError::Unsupported);
    };
    Ok(Some([reference(a)?, reference(b)?]))
}
fn raw_line_direction(source: &Source<'_>, id: u32, budget: &mut Budget<'_>) -> Result<V> {
    let id = precision::unwrapped(source, id, budget)?;
    let [_, _, direction] = source.args(id, "LINE")? else {
        return Err(ReadError::Unsupported);
    };
    precision::vector(source, reference(direction)?, 3, budget)
}
fn sphere(
    source: &Source<'_>,
    curve_id: u32,
    definition: &uv::Definition,
    surface: &Surface,
    edge: &Edge,
    tol: f64,
    budget: &mut Budget<'_>,
) -> Result<Option<Proof>> {
    let (Surface::Sphere { .. }, Curve::Line { direction, .. }) = (surface, &definition.geometry)
    else {
        return Ok(None);
    };
    if direction[0] != 0. || direction[1] == 0. {
        return Ok(None);
    }
    if !matches!(edge.curve, Curve::Circle { .. } | Curve::Ellipse { .. }) {
        return Ok(Some(Proof::Unsupported));
    }
    budget.charge(4096)?;
    // This raw coefficient field is an analytic CONTINUATION. Only the
    // comparison map earned below gives authority inside the actual V domain.
    let image = precision::analytic_image(source, definition, surface, budget)?;
    // The actual original vertex-derived finite window need not equal the
    // rounded reference window. Compose its separately earned transport once.
    let target_error = (I::point(precision::analytic_curve_reserve(
        source, curve_id, edge, budget,
    )?) + I::point(finite_axis::endpoint_reserve(
        source, curve_id, edge, budget,
    )?))
    .hi;
    let h = std::f64::consts::FRAC_PI_2.next_down();
    let mut best = f64::INFINITY;
    for shift in [-std::f64::consts::TAU, 0., std::f64::consts::TAU] {
        let Some(map) = analytic::periodic_map(image, edge, shift, budget)? else {
            continue;
        };
        let endpoints = edge.range.map(|t| I::point(t) + I::point(map.phase));
        if endpoints.into_iter().any(|p| !p.valid()) {
            continue;
        }
        let chosen = edge.range.map(|t| (t + map.phase).clamp(-h, h));
        if !chosen.into_iter().all(f64::is_finite) || chosen[0] >= chosen[1] {
            continue;
        }
        // The chosen endpoints define an increasing affine V map whose
        // ENTIRE image lies in [-h,h], strictly inside the true sphere domain.
        // Difference from t+phase is affine too, so endpoint errors suffice.
        let angle = (I::point(chosen[0]) - endpoints[0])
            .magnitude()
            .max((I::point(chosen[1]) - endpoints[1]).magnitude());
        let error = I::point(map.error) + I::point(angle) * I::point(map.speed);
        if error.valid() {
            best = best.min(error.hi)
        }
    }
    if best.is_finite() {
        let proof = finish(Proof::Certified(best), target_error, tol);
        if matches!(proof, Proof::Certified(_)) {
            precision::sense(source, curve_id, definition, surface, edge, budget)?;
            return Ok(Some(proof));
        }
    }
    // A failed domain/upper proof never becomes Invalid. A point-to-ENTIRE
    // continuation lower witness, if earned, is also a lower bound to its
    // actual half-meridian subset. Otherwise refusal remains Unsupported.
    let proof = analytic::certify_image(image, edge, tol, budget)?;
    Ok(Some(if matches!(proof, Proof::Disagreement(_)) {
        finish(proof, target_error, tol)
    } else {
        Proof::Unsupported
    }))
}
fn cone_polyline(
    source: &Source<'_>,
    curve_id: u32,
    definition: &uv::Definition,
    surface: &Surface,
    edge: &Edge,
    tol: f64,
    budget: &mut Budget<'_>,
) -> Result<Option<Proof>> {
    if !matches!(surface, Surface::Cone { .. }) {
        return Ok(None);
    }
    let Some(points) = polyline(source, definition.curve_id, budget)? else {
        return Ok(None);
    };
    if !matches!(edge.curve, Curve::Line { .. }) {
        return Ok(Some(Proof::Unsupported));
    }
    budget.charge(4096)?;
    let ua = precision::rational(source, points[0], coordinate(source, points[0], 0)?, budget)?;
    let ub = precision::rational(source, points[1], coordinate(source, points[1], 0)?, budget)?;
    budget.charge(65_536)?;
    if ua != ub {
        return Ok(Some(Proof::Unsupported));
    }
    let va = precision::rational(source, points[0], coordinate(source, points[0], 1)?, budget)?;
    let vb = precision::rational(source, points[1], coordinate(source, points[1], 1)?, budget)?;
    budget.charge(65_536)?;
    if va == vb {
        return Ok(Some(Proof::Unsupported));
    }
    let u = precision::enclose(&ua, budget)?;
    let v = [
        precision::enclose(&va, budget)?,
        precision::enclose(&vb, budget)?,
    ];
    let [_, axis, radius, angle] = source.args(definition.surface, "CONICAL_SURFACE")? else {
        return Err(ReadError::Invalid);
    };
    let frame = precision::frame(source, reference(axis)?, 3, budget)?;
    let r = precision::scalar(source, definition.surface, radius, budget)?;
    let angle = precision::scalar(source, definition.surface, angle, budget)?;
    if r.lo < 0. || angle.lo <= 0. || angle.hi >= std::f64::consts::FRAC_PI_2 {
        return Ok(Some(Proof::Unsupported));
    }
    let cosine = trig(angle, true)?;
    if cosine.lo <= 0. {
        return Ok(Some(Proof::Unsupported));
    }
    let slope = trig(angle, false)? / cosine;
    let radial = add(
        scale(frame.axes[0], trig(u, true)?),
        scale(frame.axes[1], trig(u, false)?),
    );
    let generator = add(frame.axes[2], scale(radial, slope));
    let origin = add(frame.origin, scale(radial, r));
    let controls = v.map(|v| add(origin, scale(generator, v)));
    // The actual original vertex-derived finite window need not equal the
    // rounded reference window. Compose its separately earned transport once.
    let target_error = (I::point(precision::analytic_curve_reserve(
        source, curve_id, edge, budget,
    )?) + I::point(finite_axis::endpoint_reserve(
        source, curve_id, edge, budget,
    )?))
    .hi;
    // The two actual raw endpoints define a finite segment. The existing
    // finite certificate earns positive coverage of the ACTUAL edge, allowing
    // a legitimate subinterval of an overlong master curve, never extrapolation.
    let proof = affine::nurbs(
        1,
        (&[0., 1.], &[2, 2]),
        None,
        &controls,
        true,
        edge,
        tol,
        budget,
    )?;
    let proof = finish(proof, target_error, tol);
    if matches!(proof, Proof::Certified(_)) {
        let alignment = dot(
            difference(controls[1], controls[0]),
            raw_line_direction(source, curve_id, budget)?,
        );
        if !alignment.valid() || alignment.lo <= 0. {
            return Ok(Some(Proof::Unsupported));
        }
    }
    Ok(Some(proof))
}
/// Final original-REAL correspondence; the coordinator must NOT add another
/// decoded-image reserve. Its SI-radian context gate precedes this dispatch.
#[allow(
    clippy::too_many_arguments,
    reason = "Actual source owner, 3D identity, finite trace/domain and original limits are independent proof inputs."
)]
pub(super) fn certify(
    source: &Source<'_>,
    curve_id: u32,
    definition: &uv::Definition,
    surface: &Surface,
    edge: &Edge,
    tol: f64,
    budget: &mut Budget<'_>,
) -> Result<Option<Proof>> {
    budget.charge(64)?;
    if !tol.is_finite()
        || tol <= 0.
        || edge.range.into_iter().any(|t| !t.is_finite())
        || edge.range[0] >= edge.range[1]
    {
        return Err(ReadError::Invalid);
    }
    if let Some(proof) = sphere(source, curve_id, definition, surface, edge, tol, budget)? {
        return Ok(Some(proof));
    }
    cone_polyline(source, curve_id, definition, surface, edge, tol, budget)
}
/// A collapsed original cone row is an actual geometric POINT. Equality here
/// is earned from the complete finite UV trace and SAME original location
/// owner, not from merely trusting EVALUATED_DEGENERATE_PCURVE's extra point.
pub(super) fn equivalent_vertex(
    source: &Source<'_>,
    point_id: u32,
    budget: &mut Budget<'_>,
) -> Result<Option<Point>> {
    budget.charge(64)?;
    if source.simple(point_id)?.name != "EVALUATED_DEGENERATE_PCURVE" {
        return Ok(None);
    }
    if !si_radians(source, budget)? {
        return Err(ReadError::Unsupported);
    }
    let [_, surface, representation, equivalent] =
        source.args(point_id, "EVALUATED_DEGENERATE_PCURVE")?
    else {
        return Err(ReadError::Invalid);
    };
    let surface = reference(surface)?;
    let equivalent = reference(equivalent)?;
    let [_, axis, radius, angle] = source.args(surface, "CONICAL_SURFACE")? else {
        return Err(ReadError::Unsupported);
    };
    let radius = precision::scalar(source, surface, radius, budget)?;
    let angle = precision::scalar(source, surface, angle, budget)?;
    if radius.lo != 0.
        || radius.hi != 0.
        || angle.lo <= 0.
        || angle.hi >= std::f64::consts::FRAC_PI_2
    {
        return Err(ReadError::Unsupported);
    }
    let axis = reference(axis)?;
    let [_, origin, _, _] = source.args(axis, "AXIS2_PLACEMENT_3D")? else {
        return Err(ReadError::Invalid);
    };
    if reference(origin)? != equivalent {
        return Err(ReadError::Unsupported);
    }
    precision::frame(source, axis, 3, budget)?;
    let [_, items, context] =
        source.args(reference(representation)?, "DEFINITIONAL_REPRESENTATION")?
    else {
        return Err(ReadError::Invalid);
    };
    let context = reference(context)?;
    let Some([dimension]) = source.component(context, "GEOMETRIC_REPRESENTATION_CONTEXT")? else {
        return Err(ReadError::Unsupported);
    };
    let dimension = precision::scalar(source, context, dimension, budget)?;
    if dimension.lo != 2. || dimension.hi != 2. {
        return Err(ReadError::Unsupported);
    }
    let [curve] = list(items)? else {
        return Err(ReadError::Unsupported);
    };
    let points = polyline(source, reference(curve)?, budget)?.ok_or(ReadError::Unsupported)?;
    let mut previous: Option<Q> = None;
    budget.charge(1024)?;
    let upper = Q::from_float(std::f64::consts::TAU).ok_or(ReadError::Unsupported)?;
    for point in points {
        let v = precision::scalar(source, point, coordinate(source, point, 1)?, budget)?;
        if v.lo != 0. || v.hi != 0. {
            return Err(ReadError::Unsupported);
        }
        let u = precision::rational(source, point, coordinate(source, point, 0)?, budget)?;
        // This exact binary64 TAU is below mathematical 2*pi (the independent
        // rational Machin-series control proves that lower-bound constant).
        // Bounded source rationals have at most 4096 bits; the two exact
        // comparisons and optional equality fit this precharged bit-work.
        budget.charge(65_536)?;
        if u < Q::from_integer(0.into()) || u > upper {
            return Err(ReadError::Unsupported);
        }
        if previous.as_ref() == Some(&u) {
            return Err(ReadError::Unsupported);
        }
        previous = Some(u);
    }
    budget.charge(1024)?;
    Ok(Some(source.point(equivalent)?))
}

#[cfg(test)]
mod tests;
