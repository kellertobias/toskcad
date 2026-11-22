use super::*;
use crate::tapered_thread::{Scalar, Series};
// Executed-stage ceilings: fallback 32768, two Taylor stages 32768 each,
// then ONE certified Gaussian rule (2:65536 or 4:147456). The full path is
// at most 245760 credits; unavailable stages leave later work unexecuted.
// Credits are upper operation reservations, not measured CPU time.
// Each stage is charged BEFORE running.
const FALLBACK_WORK: usize = 32_768;
// Scalar degree42 basis ~7224 and 1849 coefficient evaluations×4 ~7396,
// actual position jet <9.6k plus finite Series/checks fit 32768.
const TAYLOR_WORK: usize = 32_768;
const GAUSS2_WORK: usize = 65_536;
// 16 actual child first jets+normal² <101k; fixed norm, weights, copies
// and checks fit within 147456. No 43² polynomial evaluation at nodes.
const GAUSS4_WORK: usize = 147_456;
pub(super) struct Enclosure {
    pub value: Values,
    pub split_axis: Option<usize>,
}
#[cfg(test)]
mod tests;
fn intersect(a: I, b: I) -> Result<I, Error> {
    I::new(a.lo.max(b.lo), a.hi.min(b.hi)).checked()
}
fn vector<const N: usize>(p: [[[f64; 2]; 3]; N]) -> [Series<N>; 3] {
    std::array::from_fn(|k| Series {
        c: p.map(|p| I::new(p[k][0], p[k][1])),
    })
}
fn density<const N: usize>(
    p: [Series<N>; 3],
    norm2: Series<N>,
    center: Point,
) -> Option<[Series<N>; 4]> {
    let area = norm2.sqrt().checked().ok()?;
    let result = std::array::from_fn(|k| {
        if k == 0 {
            area
        } else {
            (p[k - 1] - Series::constant(center[k - 1])) * area
        }
    });
    for x in result {
        x.checked().ok()?;
    }
    Some(result)
}
fn fallback(
    geometry: &TensorCell11,
    normal: &SquaredNormal11,
    center: Point,
) -> Result<Values, Error> {
    let domain = geometry.domain();
    let jacobian = (I::point(domain[0][1]) - I::point(domain[0][0]))
        * (I::point(domain[1][1]) - I::point(domain[1][0]));
    let area = normal.bounds(domain)?.sqrt() * jacobian;
    let p = geometry.bounds()?.map(|b| I::new(b[0], b[1]));
    let value = std::array::from_fn(|k| {
        if k == 0 {
            area
        } else {
            (p[k - 1] - I::point(center[k - 1])) * area
        }
    });
    for x in value {
        x.checked()?;
    }
    Ok(value)
}
fn nodes(four: bool) -> Vec<(I, I)> {
    if !four {
        let delta = (I::point(1.) / I::point(3.)).sqrt() / I::point(2.);
        return vec![
            (I::point(0.5) - delta, I::point(0.5)),
            (I::point(0.5) + delta, I::point(0.5)),
        ];
    }
    let root = (I::point(6.) / I::point(5.)).sqrt();
    let outer = ((I::point(3.) + I::point(2.) * root) / I::point(7.)).sqrt() / I::point(2.);
    let inner = ((I::point(3.) - I::point(2.) * root) / I::point(7.)).sqrt() / I::point(2.);
    let wo = (I::point(18.) - I::point(30.).sqrt()) / I::point(72.);
    let wi = (I::point(18.) + I::point(30.).sqrt()) / I::point(72.);
    vec![
        (I::point(0.5) - outer, wo),
        (I::point(0.5) - inner, wi),
        (I::point(0.5) + inner, wi),
        (I::point(0.5) + outer, wo),
    ]
}
fn point_at(geometry: &TensorCell11, uv: [I; 2]) -> Option<([I; 3], I)> {
    let data = geometry.first_derivative_bounds_at(uv.map(pair)).ok()?;
    let p = data.position.map(|b| I::new(b[0], b[1]));
    let u = data.du.map(|b| I::new(b[0], b[1]));
    let v = data.dv.map(|b| I::new(b[0], b[1]));
    let n = std::array::from_fn::<_, 3, _>(|k| {
        u[(k + 1) % 3] * v[(k + 2) % 3] - u[(k + 2) % 3] * v[(k + 1) % 3]
    });
    let q = n
        .into_iter()
        .fold(I::point(0.), |s, n| s + n.square())
        .checked()
        .ok()?;
    Some((p, q))
}
fn gaussian(geometry: &TensorCell11, center: Point, four: bool) -> Option<Values> {
    let nodes = nodes(four);
    let mut sum = [I::point(0.); 4];
    for &(u, wu) in &nodes {
        for &(v, wv) in &nodes {
            let uv = [u, v];
            let (position, q) = point_at(geometry, uv)?;
            let p = position.map(Series::<1>::interval);
            // Actual child tangents already include hu and hv. Their cross
            // squared equals the root scalar field times (hu*hv)^2.
            let at = density(p, Series::<1>::interval(q), center)?;
            for k in 0..4 {
                sum[k] = sum[k] + wu * wv * at[k].c[0];
            }
        }
    }
    Some(sum)
}
pub(super) fn enclosure(
    geometry: &TensorCell11,
    normal: &SquaredNormal11,
    center: Point,
    options: SurfaceOptions11,
    budget: &mut Budget<'_>,
) -> Result<Enclosure, Error> {
    budget.charge(SurfaceStage11::Fallback, FALLBACK_WORK)?;
    budget.evaluated += 1;
    let mut result = fallback(geometry, normal, center)?;
    let mut axis_remainders = [[I::point(0.); 4]; 2];
    let mut remainder2 = [I::point(0.); 4];
    let mut remainder4 = [I::point(0.); 4];
    for along_v in [false, true] {
        budget.charge(SurfaceStage11::Taylor, TAYLOR_WORK)?;
        let p = vector(geometry.position_taylor8_bounds(along_v));
        let n = Series {
            c: normal.taylor8(geometry.domain(), along_v)?,
        };
        let axis = usize::from(along_v);
        let Some(jet) = density(p, n, center) else {
            budget.axis_unresolved[axis] += 1;
            return Ok(Enclosure {
                value: result,
                split_axis: None,
            });
        };
        budget.axis_samples[axis] += 1;
        for k in 0..4 {
            // Tensor-product positive weights: sum of the two whole-axis
            // errors covers sequential quadrature without mixed derivatives.
            // Series stores f^(n)/n!: Gauss2 c4/180, Gauss4 c8/44100.
            remainder2[k] = remainder2[k] + I::point(jet[k].c[4].magnitude()) / I::point(180.);
            remainder4[k] = remainder4[k] + I::point(jet[k].c[8].magnitude()) / I::point(44100.);
            let bound = I::point(jet[k].c[8].magnitude()) / I::point(44100.);
            axis_remainders[axis][k] = bound;
            budget.axis_remainder_max[axis][k] = budget.axis_remainder_max[axis][k].max(bound.hi);
        }
    }
    // Each rule independently includes its complete derivative remainder.
    // Avoid a duplicate certified estimate: use Gauss2 only when its proved
    // radius is no greater for all four quantities; otherwise use Gauss4.
    let four = !remainder2.iter().zip(remainder4).all(|(a, b)| a.hi <= b.hi);
    let remainder = if four { remainder4 } else { remainder2 };
    {
        budget.charge(
            if four {
                SurfaceStage11::Gauss4
            } else {
                SurfaceStage11::Gauss2
            },
            if four { GAUSS4_WORK } else { GAUSS2_WORK },
        )?;
        let Some(estimate) = gaussian(geometry, center, four) else {
            return Ok(Enclosure {
                value: result,
                split_axis: None,
            });
        };
        for k in 0..4 {
            result[k] = intersect(
                result[k],
                estimate[k] + I::new(-remainder[k].hi, remainder[k].hi),
            )?;
        }
    }
    // Refinement priority only. All final same-source interval enclosures
    // and recomposed output widths remain mandatory, independent of scores.
    let scores = axis_remainders.map(|r| {
        let mut score = r[0].hi / options.area_width_mm2;
        for (k, origin) in center.iter().enumerate() {
            score =
                score.max((r[k + 1].hi + origin.abs() * r[0].hi) / options.first_moment_width_mm3);
        }
        score
    });
    let split_axis = if scores[0] == scores[1] {
        None
    } else {
        Some(usize::from(scores[1] > scores[0]))
    };
    Ok(Enclosure {
        value: result,
        split_axis,
    })
}
