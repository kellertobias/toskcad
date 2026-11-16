use super::*;
use crate::tapered_thread::{Scalar, Series};
// Two full21x21 normal Taylor restrictions, eight derivative grids per axis,
// two actual position jets,20 tensor Gaussian evaluations with484 normal
// coefficients each, and fixed Series9 convolutions cost less than2M interval
// primitive operations. This is a delegated upper reservation, not measured
// CPU instructions; no degree9 normal/cost assumption enters this ledger.
const NODE_WORK: usize = 2_000_000;
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
    normal: [Series<N>; 3],
    center: Point,
) -> Result<[Series<N>; 4], Error> {
    let mut norm2 = normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2];
    // The same component multiplied by itself is a true square. Retain the
    // full differentiated convolution, tightening only its zeroth enclosure.
    norm2.c[0] = normal.iter().fold(I::point(0.), |s, n| s + n.c[0].square());
    let area = norm2.sqrt().checked()?;
    let result = std::array::from_fn(|k| {
        if k == 0 {
            area
        } else {
            (p[k - 1] - Series::constant(center[k - 1])) * area
        }
    });
    for x in result {
        x.checked()?;
    }
    Ok(result)
}
fn fallback(
    geometry: &TensorCell11,
    normal: &NormalPolynomial11,
    center: Point,
) -> Result<Values, Error> {
    let domain = geometry.domain();
    let jacobian = (I::point(domain[0][1]) - I::point(domain[0][0]))
        * (I::point(domain[1][1]) - I::point(domain[1][0]));
    let n = normal.lipschitz_bounds(domain)?.map(|b| I::new(b[0], b[1]));
    let area = n
        .into_iter()
        .fold(I::point(0.), |s, x| s + x.square())
        .sqrt()
        * jacobian;
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
fn gaussian(
    geometry: &TensorCell11,
    normal: &NormalPolynomial11,
    center: Point,
    four: bool,
) -> Result<Values, Error> {
    let domain = geometry.domain();
    let widths = domain.map(|d| I::point(d[1]) - I::point(d[0]));
    let jacobian = widths[0] * widths[1];
    let nodes = nodes(four);
    let mut sum = [I::point(0.); 4];
    for &(u, wu) in &nodes {
        for &(v, wv) in &nodes {
            let uv = [u, v];
            let position = geometry.first_derivative_bounds_at(uv.map(pair))?.position;
            let root = std::array::from_fn(|axis| {
                pair(I::point(domain[axis][0]) + widths[axis] * uv[axis])
            });
            let p = position.map(|b| Series::<1>::interval(I::new(b[0], b[1])));
            let n = normal
                .at(root)?
                .map(|b| Series::<1>::interval(I::new(b[0], b[1]) * jacobian));
            let at = density(p, n, center)?;
            for k in 0..4 {
                sum[k] = sum[k] + wu * wv * at[k].c[0];
            }
        }
    }
    Ok(sum)
}
pub(super) fn enclosure(
    geometry: &TensorCell11,
    normal: &NormalPolynomial11,
    center: Point,
    budget: &mut Budget<'_>,
) -> Result<Values, Error> {
    budget.charge(NODE_WORK)?;
    budget.evaluated += 1;
    let mut result = fallback(geometry, normal, center)?;
    let mut remainder2 = [I::point(0.); 4];
    let mut remainder4 = [I::point(0.); 4];
    for along_v in [false, true] {
        let p = vector(geometry.position_taylor8_bounds(along_v));
        let n = vector(normal.taylor8(geometry.domain(), along_v)?);
        let jet = match density(p, n, center) {
            Ok(jet) => jet,
            Err(Error::Unresolved) => return Ok(result),
            Err(e) => return Err(e),
        };
        for k in 0..4 {
            // Tensor-product positive weights: sum of the two whole-axis
            // errors covers sequential quadrature without mixed derivatives.
            // Series stores f^(n)/n!: Gauss2 c4/180, Gauss4 c8/44100.
            remainder2[k] = remainder2[k] + I::point(jet[k].c[4].magnitude()) / I::point(180.);
            remainder4[k] = remainder4[k] + I::point(jet[k].c[8].magnitude()) / I::point(44100.);
        }
    }
    for (four, remainder) in [(false, remainder2), (true, remainder4)] {
        let estimate = gaussian(geometry, normal, center, four)?;
        for k in 0..4 {
            result[k] = intersect(
                result[k],
                estimate[k] + I::new(-remainder[k].hi, remainder[k].hi),
            )?;
        }
    }
    Ok(result)
}
