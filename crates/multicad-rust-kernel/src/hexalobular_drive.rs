//! Retained circular Torx-like recess boundary with its actual conical floor.
//! This source does not substitute a flat bottom or certify an ISO gauge.
use crate::{head_profile::TransportError, tapered_thread::interval::I};
use std::f64::consts::{PI, TAU};

mod authoring;
mod intersection;
mod properties;
#[cfg(test)]
mod tests;
pub use authoring::TangentSpec;
pub use intersection::{CubicIntersectionSpan, IntersectionErrorBound};
pub use properties::{NominalSurfaceEvaluation, NominalSurfaceEvaluations};

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum BoundarySource {
    /// A bounded circular family; not an exact support-intersection receipt.
    RetainedFamily,
    /// Exact tangent-construction authority; arcs are its bounded evaluation.
    AuthoredTangent(TangentSpec),
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Arc {
    pub center_mm: [f64; 2],
    pub radius_mm: f64,
    pub start_mm: [f64; 2],
    pub end_mm: [f64; 2],
    pub start_angle_rad: f64,
    pub sweep_rad: f64,
}

/// In recess coordinates, depth increases into the head and the floor is
/// `depth = tip_depth_mm - radius * depth_per_radius`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Cone {
    pub tip_depth_mm: f64,
    pub depth_per_radius: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    Invalid,
    Unresolved,
    Limit,
    Cancelled,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Drive {
    arcs: [Arc; 12],
    cone: Cone,
    angles: [f64; 13],
    radial_min_squared: [f64; 12],
    radius: f64,
    error: f64,
    boundary_source: BoundarySource,
}

fn dot([ax, ay]: [f64; 2], [bx, by]: [f64; 2]) -> f64 {
    ax * bx + ay * by
}

fn ray_radius(arc: Arc, direction: [f64; 2]) -> Option<f64> {
    let projection = dot(arc.center_mm, direction);
    let perpendicular = arc.center_mm[0] * direction[1] - arc.center_mm[1] * direction[0];
    let discriminant = arc.radius_mm * arc.radius_mm - perpendicular * perpendicular;
    if discriminant <= 0. {
        return None;
    }
    Some(projection + arc.sweep_rad.signum() * discriminant.sqrt())
}

fn polar_span(
    start: [f64; 2],
    end: [f64; 2],
    minimum_radius: f64,
    error: f64,
) -> Result<f64, Error> {
    let a = start[1].atan2(start[0]);
    let b = end[1].atan2(end[0]);
    let span = (b - a).rem_euclid(TAU);
    if span <= 32. * error / minimum_radius || span >= PI {
        return Err(Error::Unresolved);
    }
    Ok(span)
}

fn volume_receipt(result: I, width_mm3: f64) -> Result<[f64; 2], Error> {
    if !result.valid() || result.lo <= 0. || result.width() > width_mm3 {
        return Err(Error::Unresolved);
    }
    Ok([result.lo, result.hi])
}

fn calculation_domain(start: f64, sweep: f64) -> I {
    let end = I::point(start) + I::point(sweep);
    let scale = 1. + start.abs().max(end.magnitude());
    let pad = (I::point(64. * f64::EPSILON) * I::point(scale)).hi;
    // Directed Gauss coordinates use subtraction, binary half scaling, one
    // midpoint addition, sqrt(3), reciprocal/multiplication and one final
    // addition/subtraction. Their accumulated absolute endpoint error is
    // below16eps*(1+max|endpoint|), including subnormals. The64eps proof
    // domain covers every representable subdivision cell and rounded-end
    // tail. It changes no literal support or finite trim.
    I::new(
        (start.min(end.lo) - pad).next_down(),
        (start.max(end.hi) + pad).next_up(),
    )
}

fn gauss_nodes(lo: f64, hi: f64) -> (I, [I; 2]) {
    let half = (I::point(hi) - I::point(lo)) * I::point(0.5);
    let midpoint = I::point(lo) + half;
    let offset = half / I::point(3.).sqrt();
    (half, [midpoint - offset, midpoint + offset])
}

impl Arc {
    fn center_support(self) -> (I, I) {
        let [x, y] = self.center_mm.map(I::point);
        let radius = (x.square() + y.square()).sqrt();
        let angle = self.center_mm[1].atan2(self.center_mm[0]);
        (
            radius,
            I::new(
                (angle - 8. * f64::EPSILON).next_down(),
                (angle + 8. * f64::EPSILON).next_up(),
            ),
        )
    }
    fn center_projection(self, angle: I) -> I {
        let (radius, phase) = self.center_support();
        radius * (angle - phase).cos()
    }
    fn radial_squared(self, angle: I) -> I {
        let [x, y] = self.center_mm.map(I::point);
        let r = I::point(self.radius_mm);
        x.square() + y.square() + r.square() + I::point(2.) * r * self.center_projection(angle)
    }

    fn area(self) -> I {
        let a = I::point(self.start_angle_rad);
        let b = a + I::point(self.sweep_rad);
        let [x, y] = self.center_mm.map(I::point);
        let r = I::point(self.radius_mm);
        (r.square() * I::point(self.sweep_rad)
            + r * (x * (b.sin() - a.sin()) - y * (b.cos() - a.cos())))
            * I::point(0.5)
    }

    fn radial(self, angle: I, radial_min_squared: f64) -> Jet {
        let (center_radius, phase) = self.center_support();
        let [x, y] = self.center_mm;
        let r = self.radius_mm;
        let center_dot = Jet::cos(angle - phase).scale_interval(center_radius);
        let mut squared =
            Jet::constant(I::point(x).square() + I::point(y).square() + I::point(r).square())
                + center_dot.scale(2. * r);
        // Intersect the value enclosure with the whole finite support's
        // positive lower bound proved at admission. This tightens an interval,
        // never the retained circle or its physical radius.
        squared.0[0].lo = squared.0[0].lo.max(radial_min_squared);
        squared.sqrt_positive()
    }
    fn moment(self, angle: I, radial_min_squared: f64) -> Jet {
        let (center_radius, phase) = self.center_support();
        let center_dot = Jet::cos(angle - phase).scale_interval(center_radius);
        let r = self.radius_mm;
        let radial = self.radial(angle, radial_min_squared);
        (radial * (Jet::constant(I::point(r).square()) + center_dot.scale(r)))
            .scale_interval(I::point(1.) / I::point(3.))
    }
}

impl Drive {
    /// Admit the actual finite twelve-arc contour. No radii or points are fitted
    /// from nominal A/B dimensions. Literal supports and shared vertices remain
    /// immutable; the exposed reserve covers their numerical coupling.
    pub fn from_retained(
        arcs: [Arc; 12],
        cone: Cone,
        transport: TransportError,
    ) -> Result<Self, Error> {
        let values = arcs.iter().flat_map(|a| {
            a.center_mm
                .into_iter()
                .chain(a.start_mm)
                .chain(a.end_mm)
                .chain([a.radius_mm, a.start_angle_rad, a.sweep_rad])
        });
        if !values.clone().all(f64::is_finite)
            || !cone.tip_depth_mm.is_finite()
            || cone.tip_depth_mm <= 0.
            || !cone.depth_per_radius.is_finite()
            || cone.depth_per_radius <= 0.
            || !transport.total_mm().is_finite()
            || transport.numeric_conversion_mm < 0.
            || transport.arithmetic_mm < 0.
        {
            return Err(Error::Invalid);
        }
        let scale = values.map(f64::abs).fold(cone.tip_depth_mm, f64::max);
        if scale > 1e7 || cone.depth_per_radius > 1e7 {
            return Err(Error::Limit);
        }
        let error = (64. * f64::EPSILON * scale.max(1.) + transport.total_mm()).next_up();
        let mut angles = [0.; 13];
        angles[0] = arcs[0].start_mm[1].atan2(arcs[0].start_mm[0]);
        let mut radius = 0_f64;
        let mut radial_min_squared = [0.; 12];
        for (i, arc) in arcs.iter().copied().enumerate() {
            if arc.end_mm != arcs[(i + 1) % 12].start_mm
                || arc.radius_mm <= 32. * error
                || arc.sweep_rad == 0.
                || arc.sweep_rad.abs() >= PI
                || arc.start_angle_rad.abs() > TAU
                || arc.radius_mm * arc.sweep_rad.abs() <= 32. * error
                || arc.sweep_rad.is_sign_positive() != i.is_multiple_of(2)
            {
                return Err(Error::Invalid);
            }
            let a = arc.start_angle_rad;
            let b = a + arc.sweep_rad;
            for (angle, point) in [(a, arc.start_mm), (b, arc.end_mm)] {
                let (s, c) = angle.sin_cos();
                if (arc.center_mm[0] + arc.radius_mm * c - point[0])
                    .hypot(arc.center_mm[1] + arc.radius_mm * s - point[1])
                    > error
                {
                    return Err(Error::Unresolved);
                }
            }
            let domain = calculation_domain(a, arc.sweep_rad);
            let r = I::point(arc.radius_mm);
            let center = arc.center_projection(domain);
            let derivative = (r.square() + r * center) * I::point(arc.sweep_rad.signum());
            let squared = arc.radial_squared(domain);
            let radial = squared.sqrt();
            if !derivative.valid()
                || derivative.lo <= 32. * error * arc.radius_mm
                || !radial.valid()
                || radial.lo <= 32. * error
            {
                return Err(Error::Unresolved);
            }
            radius = radius.max(radial.hi);
            radial_min_squared[i] = squared.lo;
            let span = polar_span(arc.start_mm, arc.end_mm, radial.lo, error)?;
            angles[i + 1] = angles[i] + span;
        }
        if ((angles[12] - angles[0]) - TAU).abs() > 128. * f64::EPSILON
            || cone.tip_depth_mm - radius * cone.depth_per_radius
                <= 32. * error * (1. + cone.depth_per_radius)
        {
            return Err(Error::Unresolved);
        }
        Ok(Self {
            arcs,
            cone,
            angles,
            radial_min_squared,
            radius,
            error,
            boundary_source: BoundarySource::RetainedFamily,
        })
    }

    /// Literal descriptors for a retained family, or bounded evaluated
    /// descriptors for an authored tangent source. Shared vertex coordinates
    /// alone do not prove literal circles intersect exactly.
    pub fn arcs(&self) -> &[Arc; 12] {
        &self.arcs
    }
    pub fn cone(&self) -> Cone {
        self.cone
    }
    pub fn boundary_source(&self) -> BoundarySource {
        self.boundary_source
    }
    pub fn source_error_mm(&self) -> f64 {
        self.error
    }
    pub fn maximum_radius_bound_mm(&self) -> f64 {
        self.radius
    }
    /// Spatial reserve after lifting a contour point to the literal cone.
    pub fn spatial_source_error_bound_mm(&self) -> f64 {
        ((I::point(1.) + I::point(self.cone.depth_per_radius)) * I::point(self.error)
            + I::point(64. * f64::EPSILON)
                * (I::point(self.cone.tip_depth_mm)
                    + I::point(self.radius) * I::point(self.cone.depth_per_radius)))
        .hi
    }

    pub fn radius_at_angle_mm(&self, angle: f64) -> Option<f64> {
        let index = self.arc_index_at_angle(angle)?;
        let theta = (angle - self.angles[0]).rem_euclid(TAU) + self.angles[0];
        ray_radius(self.arcs[index], [theta.cos(), theta.sin()])
    }
    pub fn arc_index_at_angle(&self, angle: f64) -> Option<usize> {
        if !angle.is_finite() {
            return None;
        }
        let theta = (angle - self.angles[0]).rem_euclid(TAU) + self.angles[0];
        let index = self.angles[1..12].partition_point(|&a| a < theta);
        Some(index)
    }

    pub fn floor_depth_mm(&self, point: [f64; 2]) -> Option<f64> {
        if !point.into_iter().all(f64::is_finite) {
            return None;
        }
        let radius = point[0].hypot(point[1]);
        if radius > self.radius_at_angle_mm(point[1].atan2(point[0]))? {
            return None;
        }
        Some(self.cone.tip_depth_mm - radius * self.cone.depth_per_radius)
    }

    pub fn contains(&self, point: [f64; 3]) -> bool {
        point[2] >= 0.
            && self
                .floor_depth_mm([point[0], point[1]])
                .is_some_and(|depth| point[2] <= depth)
    }

    fn source_volume_error(&self) -> I {
        let d = I::point(16.) * I::point(self.error);
        let pi = I::new(PI.next_down(), PI.next_up());
        let perimeter = self.arcs.iter().fold(I::point(0.), |sum, a| {
            sum + I::point(a.radius_mm) * I::point(a.sweep_rad.abs())
        });
        let radius = I::point(self.radius) + d;
        let depth = I::point(self.cone.tip_depth_mm)
            + d * (I::point(1.) + I::point(self.cone.depth_per_radius));
        // Planar boundary tubes cover all finite arcs and their coupled joins.
        // The cone is Lipschitz in radius with literal slope depth_per_radius.
        (I::point(2.) * perimeter * d + I::point(12.) * pi * d.square()) * depth
            + pi * radius.square() * d * (I::point(1.) + I::point(self.cone.depth_per_radius))
    }

    /// Bound a boundary approximation using a planar perimeter tube and a
    /// separate floor-height error. The retained cone/curves are unchanged.
    pub fn boundary_volume_error_bound_mm3(&self, xy_mm: f64, floor_mm: f64) -> Option<f64> {
        if ![xy_mm, floor_mm]
            .into_iter()
            .all(|v| v.is_finite() && v >= 0.)
        {
            return None;
        }
        let d = I::point(xy_mm);
        let floor = I::point(floor_mm);
        let pi = I::new(PI.next_down(), PI.next_up());
        let perimeter = self.arcs.iter().fold(I::point(0.), |sum, a| {
            sum + I::point(a.radius_mm) * I::point(a.sweep_rad.abs())
        });
        let radius = I::point(self.radius) + d;
        let result = (I::point(2.) * perimeter * d + I::point(12.) * pi * d.square())
            * (I::point(self.cone.tip_depth_mm) + floor)
            + pi * radius.square() * floor;
        result.valid().then_some(result.hi)
    }

    pub fn volume_bounds_mm3(&self, width_mm3: f64, work_limit: usize) -> Result<[f64; 2], Error> {
        self.volume_bounds_with_cancel(width_mm3, work_limit, &|| false)
    }

    pub fn volume_bounds_with_cancel(
        &self,
        width_mm3: f64,
        work_limit: usize,
        cancelled: &dyn Fn() -> bool,
    ) -> Result<[f64; 2], Error> {
        if !width_mm3.is_finite() || width_mm3 <= 0. {
            return Err(Error::Invalid);
        }
        if work_limit == 0 || work_limit > 1_000_000 {
            return Err(Error::Limit);
        }
        let source = self.source_volume_error();
        let area = self.arcs.iter().fold(I::point(0.), |sum, a| sum + a.area());
        let base = area * I::point(self.cone.tip_depth_mm);
        let available = width_mm3 - 2. * source.hi - base.width();
        if !source.valid() || available <= 0. {
            return Err(Error::Unresolved);
        }
        let mut remaining = work_limit;
        let mut moment = I::point(0.);
        for (index, arc) in self.arcs.into_iter().enumerate() {
            let a = arc.start_angle_rad;
            let b = a + arc.sweep_rad;
            // Adaptive cells have f64 endpoints. Enclose the remaining tail
            // between their rounded end and the literal support's exact sum.
            let end = I::point(a) + I::point(arc.sweep_rad);
            let tail_width = (end - I::point(b)).magnitude();
            let tail = (I::point(tail_width)
                * I::point(
                    arc.moment(
                        I::new(end.lo.min(b), end.hi.max(b)),
                        self.radial_min_squared[index],
                    )
                    .0[0]
                        .magnitude(),
                ))
            .hi;
            moment = moment + I::new(-tail, tail);
            let mut pending = vec![(a.min(b), a.max(b))];
            while let Some((lo, hi)) = pending.pop() {
                if cancelled() {
                    return Err(Error::Cancelled);
                }
                if remaining == 0 {
                    return Err(Error::Limit);
                }
                remaining -= 1;
                let span = I::point(hi) - I::point(lo);
                let (half, [first, last]) = gauss_nodes(lo, hi);
                let values = arc.moment(first, self.radial_min_squared[index]).0[0]
                    + arc.moment(last, self.radial_min_squared[index]).0[0];
                let fourth =
                    arc.moment(I::new(lo, hi), self.radial_min_squared[index]).0[4].magnitude();
                // Two-node Gauss error: |f''''|*h^5/4320. Jet coefficients
                // carry f''''/24, so the denominator is180.
                let remainder =
                    (I::point(fourth) * span.square().square() * span / I::point(180.)).hi;
                let integral = half * values + I::new(-remainder, remainder);
                let budget = available / (24. * self.cone.depth_per_radius)
                    * ((hi - lo) / arc.sweep_rad.abs());
                if integral.valid() && integral.width() <= budget {
                    moment = moment + integral * I::point(arc.sweep_rad.signum());
                } else {
                    let mid = lo + (hi - lo) * 0.5;
                    if mid <= lo || mid >= hi {
                        return Err(Error::Unresolved);
                    }
                    pending.extend([(lo, mid), (mid, hi)]);
                }
            }
        }
        let result =
            base - moment * I::point(self.cone.depth_per_radius) + I::new(-source.hi, source.hi);
        volume_receipt(result, width_mm3)
    }
}

// Normalized Taylor coefficients through order4, with every operation outward.
#[derive(Clone, Copy)]
struct Jet([I; 5]);
impl Jet {
    fn constant(value: I) -> Self {
        Self([
            value,
            I::point(0.),
            I::point(0.),
            I::point(0.),
            I::point(0.),
        ])
    }
    fn cos(t: I) -> Self {
        Self([
            t.cos(),
            -t.sin(),
            -t.cos() * I::point(0.5),
            t.sin() / I::point(6.),
            t.cos() / I::point(24.),
        ])
    }
    fn sin(t: I) -> Self {
        Self([
            t.sin(),
            t.cos(),
            -t.sin() * I::point(0.5),
            -t.cos() / I::point(6.),
            t.sin() / I::point(24.),
        ])
    }
    fn scale(self, value: f64) -> Self {
        self.scale_interval(I::point(value))
    }
    fn scale_interval(self, value: I) -> Self {
        Self(self.0.map(|v| v * value))
    }
    fn sqrt_positive(self) -> Self {
        let mut result = [I::point(0.); 5];
        result[0] = self.0[0].sqrt();
        for n in 1..5 {
            let convolution = (1..n).fold(I::point(0.), |sum, i| sum + result[i] * result[n - i]);
            result[n] = (self.0[n] - convolution) / (I::point(2.) * result[0]);
        }
        Self(result)
    }
}
impl std::ops::Add for Jet {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self(std::array::from_fn(|i| self.0[i] + rhs.0[i]))
    }
}
impl std::ops::Mul for Jet {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        Self(std::array::from_fn(|n| {
            (0..=n).fold(I::point(0.), |sum, i| sum + self.0[i] * rhs.0[n - i])
        }))
    }
}
