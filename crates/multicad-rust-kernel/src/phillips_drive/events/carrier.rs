//! Literal cubic carriers certified against whole finite authored trim curves.
use super::*;
use crate::tapered_thread::Series;
mod correspondence;
pub use correspondence::{
    TRIM_SPAN_COMPARE_WORK, TRIM_SPAN_SOURCE_WORK, TrimSpanComparison, TrimSpanCorrespondence,
};
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CubicTrimSpan {
    /// Normalized parameters on the selected finite source trim, inside[0,1].
    pub parameters: [f64; 2],
    pub poles_mm: [[f64; 3]; 4],
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TrimErrorBound {
    pub hermite_remainder_mm: f64,
    pub control_error_mm: f64,
}
impl TrimErrorBound {
    /// Bound to the exact authored curve, including numeric control evaluation.
    /// This is not an independent emitted-file support/incidence certificate.
    pub fn source_curve_error_bound_mm(self) -> f64 {
        (I::point(self.hermite_remainder_mm) + I::point(self.control_error_mm)).hi
    }
}
fn norm(p: [I; 3]) -> I {
    p.into_iter()
        .fold(I::point(0.), |sum, x| sum + x.square())
        .sqrt()
}
// Taylor4 is the fourth derivative divided by24. The two-end cubic Hermite
// remainder is at most the fourth derivative norm times h^4/384. Literal
// Bernstein-control discrepancy adds a convex-hull distance bound.
fn hermite_error(control_error: f64, fourth: [I; 3], h: I) -> Result<TrimErrorBound, Error> {
    let remainder = norm(fourth) * h.square().square() / I::point(16.);
    let bound = TrimErrorBound {
        hermite_remainder_mm: remainder.hi,
        control_error_mm: control_error,
    };
    if !bound.source_curve_error_bound_mm().is_finite() {
        return Err(Error::Unresolved);
    }
    Ok(bound)
}
impl BoundaryEvents<'_> {
    fn trim_controls(&self, index: usize, [a, b]: [f64; 2]) -> Result<[[I; 3]; 4], Error> {
        let first = evaluate::trim_state::<Series<5>>(self, index, I::point(a))?;
        let last = evaluate::trim_state::<Series<5>>(self, index, I::point(b))?;
        let third = (I::point(b) - I::point(a)) / I::point(3.);
        Ok([
            first.map(|s| s.c[0]),
            first.map(|s| s.c[0] + third * s.c[1]),
            last.map(|s| s.c[0] - third * s.c[1]),
            last.map(|s| s.c[0]),
        ])
    }
    pub fn certify_trim_span(
        &self,
        index: usize,
        span: &CubicTrimSpan,
    ) -> Result<TrimErrorBound, Error> {
        let [a, b] = span.parameters;
        if index >= 4
            || !a.is_finite()
            || !b.is_finite()
            || a < 0.
            || b > 1.
            || a >= b
            || !span.poles_mm.into_iter().flatten().all(f64::is_finite)
        {
            return Err(Error::Invalid);
        }
        if span.poles_mm.into_iter().flatten().any(|v| v.abs() > 1e9) {
            return Err(Error::Limit);
        }
        let ideal = self.trim_controls(index, [a, b])?;
        let control_error = span
            .poles_mm
            .into_iter()
            .zip(ideal)
            .map(|(p, h)| norm(std::array::from_fn(|i| I::point(p[i]) - h[i])).hi)
            .fold(0., f64::max);
        let whole = evaluate::trim_state::<Series<5>>(self, index, I::new(a, b))?;
        let h = I::point(b) - I::point(a);
        hermite_error(control_error, whole.map(|s| s.c[4]), h)
    }
    /// Certify literal controls in an explicit affine source chart. Reflected
    /// corner orbits may preserve depth while reversing one transverse axis.
    /// This evaluates curve distance; it does not admit a frame or solid.
    pub fn certify_trim_span_in_frame(
        &self,
        index: usize,
        span: &CubicTrimSpan,
        origin_mm: [f64; 3],
        axes: [[f64; 3]; 3],
    ) -> Result<TrimErrorBound, Error> {
        let [a, b] = span.parameters;
        if index >= 4
            || !a.is_finite()
            || !b.is_finite()
            || a < 0.
            || b > 1.
            || a >= b
            || !span
                .poles_mm
                .into_iter()
                .flatten()
                .chain(origin_mm)
                .chain(axes.into_iter().flatten())
                .all(f64::is_finite)
        {
            return Err(Error::Invalid);
        }
        if span
            .poles_mm
            .into_iter()
            .flatten()
            .chain(origin_mm)
            .any(|v| v.abs() > 1e12)
            || axes.into_iter().flatten().any(|v| v.abs() > 1e9)
        {
            return Err(Error::Limit);
        }
        let vector = |p: [I; 3]| {
            std::array::from_fn(|i| {
                (0..3).fold(I::point(0.), |sum, j| sum + I::point(axes[j][i]) * p[j])
            })
        };
        let ideal = self.trim_controls(index, [a, b])?.map(|p| {
            let v = vector(p);
            std::array::from_fn::<_, 3, _>(|i| I::point(origin_mm[i]) + v[i])
        });
        let control_error = span
            .poles_mm
            .into_iter()
            .zip(ideal)
            .map(|(p, ideal)| norm(std::array::from_fn(|i| I::point(p[i]) - ideal[i])).hi)
            .fold(0., f64::max);
        let whole = evaluate::trim_state::<Series<5>>(self, index, I::new(a, b))?;
        hermite_error(
            control_error,
            vector(whole.map(|s| s.c[4])),
            I::point(b) - I::point(a),
        )
    }
    pub fn trim_carrier(
        &self,
        index: usize,
        error_mm: f64,
        span_limit: usize,
    ) -> Result<Vec<CubicTrimSpan>, Error> {
        self.trim_carrier_with_cancel(index, error_mm, span_limit, &|| false)
    }
    pub fn trim_carrier_with_cancel(
        &self,
        index: usize,
        error_mm: f64,
        span_limit: usize,
        cancel: &dyn Fn() -> bool,
    ) -> Result<Vec<CubicTrimSpan>, Error> {
        if index >= 4 || !error_mm.is_finite() || error_mm <= 0. {
            return Err(Error::Invalid);
        }
        if span_limit == 0 || span_limit > 1_000_000 {
            return Err(Error::Limit);
        }
        let mut pending = vec![[0., 1.]];
        let mut result = Vec::new();
        while let Some(parameters) = pending.pop() {
            checkpoint(cancel)?;
            if result.len() + pending.len() >= span_limit {
                return Err(Error::Limit);
            }
            let candidate = self.trim_controls(index, parameters).map(|controls| {
                let mut poles = controls.map(|p| p.map(I::midpoint));
                for (endpoint, role) in [
                    (0, self.trims[index].ends[0]),
                    (3, self.trims[index].ends[1]),
                ] {
                    if parameters[endpoint / 3] == (endpoint / 3) as f64 {
                        poles[endpoint] = self.event(role).point.map(I::midpoint)
                    }
                }
                CubicTrimSpan {
                    parameters,
                    poles_mm: poles,
                }
            });
            if let Ok(span) = candidate
                && let Ok(bound) = self.certify_trim_span(index, &span)
                && bound.source_curve_error_bound_mm() <= error_mm
            {
                result.push(span)
            } else {
                let [a, b] = parameters;
                let mid = a + (b - a) / 2.;
                if mid <= a || mid >= b {
                    return Err(Error::Unresolved);
                }
                pending.extend([[mid, b], [a, mid]]);
            }
        }
        Ok(result)
    }
}
#[cfg(test)]
mod tests;
