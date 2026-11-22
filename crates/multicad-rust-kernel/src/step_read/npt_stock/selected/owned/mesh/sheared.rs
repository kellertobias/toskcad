//! Whole ORIGINAL-field certificates over sheared UV triangles.
use super::*;
fn cross(a: [I; 3], b: [I; 3]) -> [I; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dot(a: [I; 3], b: [I; 3]) -> I {
    a.into_iter()
        .zip(b)
        .fold(I::point(0.), |sum, (a, b)| sum + a * b)
}
fn length(a: [I; 3]) -> I {
    a.into_iter()
        .fold(I::point(0.), |sum, x| sum + x.square())
        .sqrt()
}
pub(super) struct Quality {
    pub position: f64,
    pub angle: f64,
}
fn aligned_quality(
    root: &normal::Normal,
    edges: [[I; 2]; 2],
    normals: [[I; 3]; 3],
    facet: [I; 3],
    rounding: f64,
    budget: &mut Budget<'_>,
) -> Result<Quality> {
    budget.step(512)?;
    let edges = edges.map(|d| [d[0] + I::point(root.steering[0]) * d[1], d[1]]);
    let weights = hessian::Hessian::weights(edges);
    let combine = |fields: [[I; 3]; 3], weight: [I; 3]| {
        std::array::from_fn(|k| (0..3).fold(I::point(0.), |sum, j| sum + weight[j] * fields[j][k]))
    };
    let remainder = |a: [f64; 3]| {
        ((I::point(a[0]) + I::point(2.) * I::point(a[1]) + I::point(a[2])) / I::point(8.)).hi
    };
    let position = (I::point(remainder(
        weights.map(|w| cells::l1(combine(root.position_aligned, w))),
    )) + I::point(rounding))
    .hi;
    let derivatives = weights.map(|w| combine(root.normal_aligned, w));
    let axial = remainder(derivatives.map(|p| dot(p, facet).magnitude()));
    let perpendicular = remainder(derivatives.map(|p| length(cross(p, facet)).hi));
    let mut angle = 0f64;
    for normal in normals {
        let denominator = dot(normal, facet) - I::point(axial);
        if !denominator.valid() || denominator.lo <= 0. {
            return Ok(Quality {
                position,
                angle: f64::INFINITY,
            });
        }
        let ratio = (I::point(length(cross(normal, facet)).hi) + I::point(perpendicular))
            / I::point(denominator.lo);
        if !ratio.valid() {
            return Ok(Quality {
                position,
                angle: f64::INFINITY,
            });
        }
        angle = angle.max(ratio.hi);
    }
    Ok(Quality { position, angle })
}
fn normal_error(
    root: &normal::Normal,
    weights: [[I; 3]; 3],
    corners: [[I; 3]; 3],
    facet: [I; 3],
    budget: &mut Budget<'_>,
) -> Result<f64> {
    budget.step(root.hessian.coefficients[0].len() * 144)?;
    let mut axial = [0f64; 3];
    let mut perpendicular = [0f64; 3];
    for (direction, weight) in weights.into_iter().enumerate() {
        for index in 0..root.hessian.coefficients[0].len() {
            let derivative = root.hessian.point(index, weight);
            axial[direction] = axial[direction].max(dot(derivative, facet).magnitude());
            perpendicular[direction] =
                perpendicular[direction].max(length(cross(derivative, facet)).hi);
        }
    }
    let remainder = |a: [f64; 3]| {
        ((I::point(a[0]) + I::point(2.) * I::point(a[1]) + I::point(a[2])) / I::point(8.)).hi
    };
    let axial = remainder(axial);
    let perpendicular = remainder(perpendicular);
    let mut angle = 0f64;
    for normal in corners {
        let denominator = dot(normal, facet) - I::point(axial);
        if !denominator.valid() || denominator.lo <= 0. {
            return Ok(f64::INFINITY);
        }
        let numerator = I::point(length(cross(normal, facet)).hi) + I::point(perpendicular);
        let quotient = numerator / I::point(denominator.lo);
        if !quotient.valid() {
            return Ok(f64::INFINITY);
        }
        angle = angle.max(quotient.hi);
    }
    // Taylor's integral remainder is bounded on the convex UV triangle.
    // Linear interpolation has error <=(Mss+2Mst+Mtt)/8: the weighted
    // coordinate variances are <=1/4 and their mixed absolute moment<=1/4.
    // Projection onto/cross the ACTUAL facet normal happens BEFORE each
    // original derivative hull. This preserves helical shear cancellation.
    // Positive corner cones enlarged by that same error imply one cone for
    // every source normal, by barycentric convexity; theta=atan(r)<=r.
    Ok(angle)
}
#[cfg(test)]
pub(super) fn qualify(
    root: &normal::Normal,
    uv: [[f64; 2]; 3],
    points: [[f64; 3]; 3],
    rounding: f64,
    budget: &mut Budget<'_>,
) -> Result<Quality> {
    let mut normals = [[I::point(0.); 3]; 3];
    for (i, at) in uv.into_iter().enumerate() {
        normals[i] = root.point(at, budget)?;
    }
    qualify_cached(root, uv, points, rounding, normals, [0.1, 0.3], budget)
}
pub(super) fn qualify_cached(
    root: &normal::Normal,
    uv: [[f64; 2]; 3],
    points: [[f64; 3]; 3],
    rounding: f64,
    normals: [[I; 3]; 3],
    limits: [f64; 2],
    budget: &mut Budget<'_>,
) -> Result<Quality> {
    if uv
        .into_iter()
        .flatten()
        .any(|x| !x.is_finite() || !(0.0..=1.0).contains(&x))
        || points.into_iter().flatten().any(|x| !x.is_finite())
    {
        return Err(Error::Invalid);
    }
    let edges = [1, 2].map(|i| std::array::from_fn(|k| I::point(uv[i][k]) - I::point(uv[0][k])));
    let determinant = edges[0][0] * edges[1][1] - edges[0][1] * edges[1][0];
    if !determinant.valid() || determinant.hi >= 0. {
        return Err(Error::Unsupported);
    }
    let weights = hessian::Hessian::weights(edges);
    let [a, b, c] = points.map(|p| p.map(I::point));
    let facet = cross(
        std::array::from_fn(|k| b[k] - a[k]),
        std::array::from_fn(|k| c[k] - a[k]),
    );
    let quick = aligned_quality(root, edges, normals, facet, rounding, budget)?;
    if quick.position <= limits[0] && quick.angle <= limits[1] {
        return Ok(quick);
    }
    // These are independent whole-field certificates for the SAME immutable
    // root, UV triangle and literal facet. A failed position bound does not
    // invalidate an already successful normal cone (or vice versa). Compute
    // only the component still needing a certificate; no caller limit changes.
    let position = if quick.position <= limits[0] {
        quick.position
    } else {
        (I::point(root.position_hessian.position_error(weights, budget)?) + I::point(rounding)).hi
    };
    let angle = if quick.angle <= limits[1] {
        quick.angle
    } else {
        normal_error(root, weights, normals, facet, budget)?
    };
    Ok(Quality { position, angle })
}
