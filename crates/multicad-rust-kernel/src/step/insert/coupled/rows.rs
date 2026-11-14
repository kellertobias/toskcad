//! Shared angular polynomial times ORIGINAL exact affine radius/height laws.
use super::number::{Grid, charge, rounded_ratio};
use super::*;
const BITS: usize = 96;
type AngularPolynomial = [[N; 2]; 11];
fn split(
    p: &AngularPolynomial,
    t: &N,
    w: &mut Work<'_>,
) -> Result<(AngularPolynomial, AngularPolynomial), Error> {
    let mut p = p.clone();
    let mut left = p.clone();
    let mut right = p.clone();
    for level in 1..11 {
        for i in 0..11 - level {
            let next = p[i + 1].clone();
            for (current, next) in p[i].iter_mut().zip(next) {
                *current = current.add(&next.sub(current, w)?.mul(t, w)?, w)?;
            }
        }
        left[level] = p[0].clone();
        right[10 - level] = p[10 - level].clone();
    }
    Ok((left, right))
}
pub(super) struct Slice<'d, 's> {
    scope: admitted::Scope<'d, 's>,
    poles: [[Grid; 2]; 11],
    errors: [f64; 2],
}
pub(super) fn angular_slice<'d, 's>(
    span: admitted::AngularSpan<'_, 'd, 's>,
    w: &mut Work<'_>,
) -> Result<Slice<'d, 's>, Error> {
    let parent = span.parent();
    let begin = span.begin();
    let end = span.end();
    let p =
        std::array::from_fn(|i| std::array::from_fn(|axis| N(parent[i][axis].actual(0., BITS))));
    let (_, right) = split(&p, begin, w)?;
    let ratio = end.sub(begin, w)?.div(&N::integer(1).sub(begin, w)?, w)?;
    let exact = split(&right, &ratio, w)?.0;
    let mut poles = std::array::from_fn(|_| std::array::from_fn(|_| Grid(BigInt::zero())));
    for i in 0..11 {
        for axis in 0..2 {
            poles[i][axis] = Grid::rational(&exact[i][axis], BITS, w)?;
        }
    }
    let width = span.width().enclosure(w)?;
    let q = I::point(f64::from_bits((1022 - BITS as u64) << 52));
    Ok(Slice {
        scope: span.scope(),
        poles,
        errors: [
            (I::point(span.errors()[0]) + I::point(2.) * q).hi,
            (I::point(span.errors()[1]) + I::point(40.) * q / width).hi,
        ],
    })
}
pub(super) struct Row {
    pub poles: [[Grid; 3]; 12],
    pub source_trim: [N; 2],
    pub errors: [f64; 2],
}
fn affine_value(a: &N, b: &N, t: &N, w: &mut Work<'_>) -> Result<N, Error> {
    if b.0.is_zero() {
        w.charge(1 + (a.0.numer().bits() + a.0.denom().bits()).div_ceil(64) as usize)?;
        return Ok(a.clone());
    }
    if a.0.is_zero() {
        return b.mul(t, w);
    }
    a.add(&b.mul(t, w)?, w)
}
pub(super) fn row(
    angular: &Slice<'_, '_>,
    request: &admitted::OriginalRow<'_, '_>,
    options: CurveOptions,
    w: &mut Work<'_>,
) -> Result<Row, Error> {
    w.charge(1)?;
    if !angular.scope.same(&request.scope) {
        return Err(Error::Invalid);
    }
    let affine = request.affine();
    let t = request.parameters();
    let mut radius = [N::integer(0), N::integer(0)];
    let mut z = radius.clone();
    for i in 0..2 {
        radius[i] = affine_value(&affine[0][0], &affine[1][0], &t[i], w)?;
        z[i] = affine_value(&affine[0][1], &affine[1][1], &t[i], w)?;
    }
    // OriginalRow earns the SAME source's closed finite curve interval. Its
    // affine endpoint radii lie between stored band radii<=majorRadius; heights
    // lie between the original clipped0/L endpoints. Spec already bounds both
    // dimensions by1e7. An outward enclosure may exceed that limit and is not a
    // new source-validity test. Actual numeric/quality checks remain below.
    let r_grid = [
        Grid::rational(&radius[0], BITS, w)?,
        Grid::rational(&radius[1], BITS, w)?,
    ];
    let z_grid = [
        Grid::rational(&z[0], BITS, w)?,
        Grid::rational(&z[1], BITS, w)?,
    ];
    let mut poles = std::array::from_fn(|_| std::array::from_fn(|_| Grid(BigInt::zero())));
    for (i, point) in poles.iter_mut().enumerate() {
        for (axis, coordinate) in point[..2].iter_mut().enumerate() {
            // Two coefficient products, two small weights, sum, quotient and
            // tie correction. Charge all eight bounded limb operations.
            let angular_bits = if i == 0 {
                angular.poles[0][axis].0.bits()
            } else if i == 11 {
                angular.poles[10][axis].0.bits()
            } else {
                angular.poles[i][axis]
                    .0
                    .bits()
                    .max(angular.poles[i - 1][axis].0.bits())
            };
            for _ in 0..8 {
                charge(
                    r_grid[0].0.bits().max(r_grid[1].0.bits()) + angular_bits + 8,
                    w,
                )?;
            }
            let a = if i < 11 {
                &r_grid[0].0 * &angular.poles[i][axis].0 * (11 - i)
            } else {
                BigInt::zero()
            };
            let b = if i > 0 {
                &r_grid[1].0 * &angular.poles[i - 1][axis].0 * i
            } else {
                BigInt::zero()
            };
            *coordinate = Grid(rounded_ratio(a + b, &(BigInt::from(11) << BITS)));
        }
        for _ in 0..5 {
            charge(z_grid[0].0.bits().max(z_grid[1].0.bits()) + 8, w)?;
        }
        point[2] = Grid(rounded_ratio(
            &z_grid[0].0 * (11 - i) + &z_grid[1].0 * i,
            &BigInt::from(11),
        ));
    }
    let r = radius[0]
        .enclosure(w)?
        .magnitude()
        .max(radius[1].enclosure(w)?.magnitude());
    let dr = affine[1][0].enclosure(w)?.magnitude();
    let dv = affine[1][2].enclosure(w)?.magnitude();
    let gain = t[1].sub(&t[0], w)?;
    let gain = N(gain.0.abs()).enclosure(w)?;
    let q = I::point(f64::from_bits((1022 - BITS as u64) << 52));
    let maximum_angle = I::point(2.) + I::point(angular.errors[0]);
    let position = I::point(r) * I::point(angular.errors[0]) + (maximum_angle + I::point(4.)) * q;
    let derivative = I::point(dr) * I::point(angular.errors[0])
        + I::point(r) * I::point(dv) * I::point(angular.errors[1])
        + q * I::point(dv) * (I::point(10.) + I::point(angular.errors[1]))
        + (I::point(2.) * maximum_angle + I::point(68.)) * q / gain;
    if !position.valid()
        || !derivative.valid()
        || position.hi > options.position_error_mm
        || derivative.hi > options.derivative_error_mm
    {
        return Err(Error::Unresolved);
    }
    Ok(Row {
        poles,
        source_trim: t.clone(),
        errors: [position.hi, derivative.hi],
    })
}
