//! Finite cubic carriers of the actual cylinder/cone intersection. Certification
//! consumes literal controls; an authoring recipe is never a readback proof.
use super::*;
mod polyline;
#[cfg(test)]
mod tests;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CubicIntersectionSpan {
    /// Normalized retained arc parameter; theta=start+sweep*u.
    pub parameters: [f64; 2],
    /// XY contour chart, positive Z depth below the open mouth.
    pub poles_mm: [[f64; 3]; 4],
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IntersectionErrorBound {
    pub nominal_curve_error_mm: f64,
    pub source_error_mm: f64,
}
impl IntersectionErrorBound {
    pub fn total_error_bound_mm(self) -> f64 {
        (I::point(self.nominal_curve_error_mm) + I::point(self.source_error_mm)).hi
    }
}
fn norm(value: [I; 3]) -> I {
    value
        .into_iter()
        .fold(I::point(0.), |sum, v| sum + v.square())
        .sqrt()
}
impl Drive {
    pub(super) fn intersection_state(&self, index: usize, u: f64) -> ([I; 3], [I; 3]) {
        let arc = self.arcs[index];
        let angle = I::point(arc.start_angle_rad) + I::point(arc.sweep_rad) * I::point(u);
        let radial = arc.radial(angle, self.radial_min_squared[index]);
        let r = I::point(arc.radius_mm);
        let k = I::point(self.cone.depth_per_radius);
        let point = [
            I::point(arc.center_mm[0]) + r * angle.cos(),
            I::point(arc.center_mm[1]) + r * angle.sin(),
            I::point(self.cone.tip_depth_mm) - k * radial.0[0],
        ];
        let derivative = [-r * angle.sin(), r * angle.cos(), -k * radial.0[1]]
            .map(|v| v * I::point(arc.sweep_rad));
        (point, derivative)
    }
    fn intersection_controls(&self, index: usize, [a, b]: [f64; 2]) -> [[I; 3]; 4] {
        let (first, first_derivative) = self.intersection_state(index, a);
        let (last, last_derivative) = self.intersection_state(index, b);
        let third = (I::point(b) - I::point(a)) / I::point(3.);
        [
            first,
            std::array::from_fn(|i| first[i] + third * first_derivative[i]),
            std::array::from_fn(|i| last[i] - third * last_derivative[i]),
            last,
        ]
    }
    /// Bound a literal polynomial against the whole finite retained arc lifted
    /// to the actual cone. This proves both support bindings within the bound;
    /// it does not assert that a finite polynomial is an exact circle/cone curve.
    pub fn certify_intersection_span(
        &self,
        index: usize,
        span: &CubicIntersectionSpan,
    ) -> Result<IntersectionErrorBound, Error> {
        let [a, b] = span.parameters;
        if index >= 12
            || ![a, b].into_iter().all(f64::is_finite)
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
        let ideal = self.intersection_controls(index, [a, b]);
        let control_error = span
            .poles_mm
            .into_iter()
            .zip(ideal)
            .map(|(p, h)| norm(std::array::from_fn(|i| I::point(p[i]) - h[i])).hi)
            .fold(0., f64::max);
        let arc = self.arcs[index];
        let angle = I::point(arc.start_angle_rad) + I::point(arc.sweep_rad) * I::new(a, b);
        let radial = arc.radial(angle, self.radial_min_squared[index]);
        let fourth_z = I::point(24.)
            * I::point(self.cone.depth_per_radius)
            * I::point(radial.0[4].magnitude());
        let derivative =
            (I::point(2.) * I::point(arc.radius_mm).square() + fourth_z.square()).sqrt();
        let theta_span = I::point(arc.sweep_rad) * (I::point(b) - I::point(a));
        // Cubic two-end Hermite error <= |f''''| * h^4/384.
        // The coordinate bounds give a conservative Euclidean remainder.
        let remainder = derivative * theta_span.square().square() / I::point(384.);
        let result = I::point(control_error) + remainder;
        // Admission gives supports<=1e7, k<=1e7, rho>32*64eps. The
        // squared-radius coefficients are below4e14 and rho>4e-13. The sqrt
        // recurrence gives |q1..q4|<5e26,4e65,6e104,1e144 respectively;
        // even the squared derivative norm stays below1e307. Literal caller
        // poles<=1e9 also keep the control-error norm finite.
        Ok(IntersectionErrorBound {
            nominal_curve_error_mm: result.hi,
            source_error_mm: self.spatial_source_error_bound_mm(),
        })
    }
    pub fn intersection_carrier(
        &self,
        index: usize,
        error_mm: f64,
        span_limit: usize,
    ) -> Result<Vec<CubicIntersectionSpan>, Error> {
        self.intersection_carrier_with_cancel(index, error_mm, span_limit, &|| false)
    }
    pub fn intersection_carrier_with_cancel(
        &self,
        index: usize,
        error_mm: f64,
        span_limit: usize,
        cancelled: &dyn Fn() -> bool,
    ) -> Result<Vec<CubicIntersectionSpan>, Error> {
        if index >= 12 || !error_mm.is_finite() || error_mm <= 0. {
            return Err(Error::Invalid);
        }
        if span_limit == 0 || span_limit > 1_000_000 {
            return Err(Error::Limit);
        }
        if error_mm <= 4. * self.spatial_source_error_bound_mm() {
            return Err(Error::Unresolved);
        }
        let mut pending = vec![[0., 1.]];
        let mut result = Vec::new();
        while let Some(parameters) = pending.pop() {
            if cancelled() {
                return Err(Error::Cancelled);
            }
            if result.len() + pending.len() >= span_limit {
                return Err(Error::Limit);
            }
            let mut poles = self
                .intersection_controls(index, parameters)
                .map(|p| p.map(I::midpoint));
            let arc = self.arcs[index];
            for (end, point) in [(0, arc.start_mm), (3, arc.end_mm)] {
                if parameters[end / 3] == (end / 3) as f64 {
                    poles[end] = [
                        point[0],
                        point[1],
                        self.cone.tip_depth_mm
                            - self.cone.depth_per_radius * point[0].portable_hypot(point[1]),
                    ];
                }
            }
            let span = CubicIntersectionSpan {
                parameters,
                poles_mm: poles,
            };
            if self
                .certify_intersection_span(index, &span)?
                .total_error_bound_mm()
                <= error_mm
            {
                result.push(span);
            } else {
                let [a, b] = parameters;
                let mid = a + (b - a) * 0.5;
                if mid <= a || mid >= b {
                    return Err(Error::Unresolved);
                }
                pending.extend([[mid, b], [a, mid]]);
            }
        }
        Ok(result)
    }
}
