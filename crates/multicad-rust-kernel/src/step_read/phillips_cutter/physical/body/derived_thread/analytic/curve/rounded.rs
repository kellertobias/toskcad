//! Whole-cell actual circular-root radial/trigonometric factor, without a
//! flat-root surrogate or replacing the original cached helical rate.
use super::*;
use crate::tapered_thread::{Scalar, Series};
type S = Series<18>;
pub(super) struct Factor {
    pub poles: Vec<[I; 2]>,
    pub remainder: [f64; 2],
    pub radius: f64,
    pub du: f64,
}
fn value(
    root: crate::helical_thread::RootCircle,
    angle_knots: [I; 2],
    rate: I,
    hand: I,
    t: I,
) -> Result<[S; 2]> {
    let angle =
        S::interval(angle_knots[0]) + S::variable(t) * S::interval(angle_knots[1] - angle_knots[0]);
    let [sin_angle, cos_angle] = angle.sin_cos();
    let radius = S::constant(root.center_radius) - S::constant(root.radius) * cos_angle;
    let axial_phase = S::constant(root.phase) + S::constant(root.radius) * sin_angle;
    let theta = axial_phase * S::interval(-hand / rate);
    let [sin, cos] = theta.sin_cos();
    let a = [radius * cos, radius * sin];
    Ok([
        a[0].checked().map_err(|_| ReadError::Unsupported)?,
        a[1].checked().map_err(|_| ReadError::Unsupported)?,
    ])
}
pub(super) fn bounds(
    root: crate::helical_thread::RootCircle,
    angle_knots: [I; 2],
    rate: I,
    hand: I,
    work: &mut Work<'_>,
) -> Result<Factor> {
    // Two bounded Series18 evaluations, nested trig/multiply recurrences,
    // coefficient conversion and full interval radius/derivative metadata.
    // The real source root uses ANGLE, not linear phase, as its cell chart.
    work.step(131072)?;
    let delta = angle_knots[1] - angle_knots[0];
    let domain = I::new(angle_knots[0].lo, angle_knots[1].hi);
    let cosine = domain.cos();
    let sine = domain.sin();
    let radicand = I::point(root.radius).square() - (I::point(root.radius) * sine).square();
    supported(
        angle_knots
            .into_iter()
            .all(|p| p.valid() && p.lo.is_finite() && p.hi.is_finite())
            && rate.valid()
            && hand.valid()
            && hand.lo == hand.hi
            && (hand.lo == 1. || hand.lo == -1.)
            && delta.lo > 0.
            && delta.hi.is_finite()
            && rate.lo > 0.
            && rate.hi.is_finite()
            && root.radius > 0.
            && root.radius.is_finite()
            && root.center_radius.is_finite()
            && root.phase.is_finite()
            && cosine.lo > 0.
            && cosine.hi.is_finite()
            && radicand.lo > 0.
            && radicand.hi.is_finite(),
    )?;
    let radial = I::point(root.center_radius) - I::point(root.radius) * cosine;
    supported(radial.lo > 0. && radial.hi.is_finite())?;
    let center = value(root, angle_knots, rate, hand, I::point(0.5))?;
    let whole = value(root, angle_knots, rate, hand, I::new(0., 1.))?;
    let mut ordinary = [[I::point(0.); 2]; 17];
    for k in 0..=16 {
        for (i, a) in ordinary.iter_mut().enumerate().take(k + 1) {
            let w = I::point(choose(k, i)) * I::point((-0.5f64).powi((k - i) as i32));
            for axis in 0..2 {
                a[axis] = a[axis] + center[axis].c[k] * w;
            }
        }
    }
    let poles = (0..=16)
        .map(|i| {
            let mut p = [I::point(0.); 2];
            for (k, c) in ordinary.iter().enumerate().take(i + 1) {
                let w = I::point(choose(i, k)) / I::point(choose(16, k));
                for axis in 0..2 {
                    p[axis] = p[axis] + c[axis] * w;
                }
            }
            p
        })
        .collect();
    // Series coefficient17 is f^(17)/17! at EVERY point in the closed
    // original normalized cell. Taylor and derivative-Taylor remainders are
    // respectively c17*h^17 and17*c17*h^16, h<=1/2. Both components enter
    // the Euclidean bound; neither endpoint samples nor a fitted radius suffice.
    let c17 = (whole[0].c[17].square() + whole[1].c[17].square()).sqrt();
    let half16 = I::point(2f64.powi(-16));
    let remainder = [
        (c17 * half16 / I::point(2.)).hi,
        (c17 * I::point(17.) * half16).hi,
    ];
    let dr = delta * I::point(root.radius) * sine;
    let phase_du = delta * I::point(root.radius) * cosine;
    let du = (I::point(dr.magnitude())
        + I::point(radial.hi) * I::point((phase_du / rate).magnitude()))
    .hi;
    supported(
        remainder
            .into_iter()
            .chain([du])
            .all(|x| x.is_finite() && x >= 0.),
    )?;
    Ok(Factor {
        poles,
        remainder,
        radius: radial.hi,
        du,
    })
}

#[cfg(test)]
mod tests;

fn supported(ok: bool) -> Result<()> {
    if ok {
        Ok(())
    } else {
        Err(ReadError::Unsupported)
    }
}
