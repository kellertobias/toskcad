//! Source-neutral reuse of the existing whole-curve Hermite certificate.
//! Parameter t is normalized to [0,1]; physical turns=v0+(v1-v0)*t.
use super::*;
pub(crate) struct AngularCarrier {
    pub knots: Vec<f64>,
    /// x/y approximate the true unit circle; z approximates PHYSICAL turns.
    pub poles: Vec<[f64; 3]>,
    pub position_error: f64,
    /// Derivative error with respect to the WHOLE normalized t, not radians.
    pub derivative_error: f64,
}
pub(crate) fn fit_angular(
    turns: [f64; 2],
    left: bool,
    options: CurveOptions,
    memory: &mut Memory,
    w: &mut Work<'_>,
) -> Result<AngularCarrier, Error> {
    fit(turns, left, options, None, memory, w)
}
/// Only this source-owned unit circle has the proved exact periodic endpoint.
/// No tolerance/proximity equates endpoints of arbitrary curve domains.
pub(crate) fn fit_closed_unit(
    left: bool,
    options: CurveOptions,
    memory: &mut Memory,
    w: &mut Work<'_>,
) -> Result<AngularCarrier, Error> {
    fit(
        [0., 1.],
        left,
        options,
        Some([[1., 0., 0.], [1., 0., 1.]]),
        memory,
        w,
    )
}
fn fit(
    turns: [f64; 2],
    left: bool,
    options: CurveOptions,
    ends: Option<[[f64; 3]; 2]>,
    memory: &mut Memory,
    w: &mut Work<'_>,
) -> Result<AngularCarrier, Error> {
    w.charge(1)?;
    options.validate()?;
    if turns.iter().any(|v| !v.is_finite() || v.abs() > 256.) || turns[0] >= turns[1] {
        return Err(Error::Invalid);
    }
    memory.add(WORKSPACE_BYTES)?;
    let p = I::point;
    let delta = p(turns[1]) - p(turns[0]);
    let law = curves::Law::from_bounds(
        [[p(1.), p(turns[0]), p(turns[0])], [p(0.), delta, delta]],
        left,
    );
    let ends = match ends {
        Some(owned) => owned,
        None => [law.point(0., w)?, law.point(1., w)?].map(|v| v.map(I::midpoint)),
    };
    let basis = curves::derivative_basis(w)?;
    let mut n = 1;
    loop {
        let fit = curves::fixed(&law, n, ends, None, &basis, memory, w)?;
        if fit.position_error_mm <= options.position_error_mm
            && fit.derivative_error_mm <= options.derivative_error_mm
        {
            memory.used -= WORKSPACE_BYTES;
            return Ok(AngularCarrier {
                knots: fit.knots,
                poles: fit.poles,
                position_error: fit.position_error_mm,
                derivative_error: fit.derivative_error_mm,
            });
        }
        memory.release(&fit.knots);
        memory.release(&fit.poles);
        drop(fit);
        if n > options.maximum_spans_per_curve / 2 {
            return Err(Error::Unresolved);
        }
        n *= 2;
    }
}
