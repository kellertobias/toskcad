//! Immutable finite source graph. File identity/incidence proof belongs to the
//! decoder; this type independently certifies the retained chart geometry.
use super::*;
use std::sync::Arc as Shared;
#[cfg(test)]
mod tests;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteVertices {
    pub mouth_mm: [[f64; 3]; 12],
    pub floor_mm: [[f64; 3]; 12],
    pub apex_mm: [f64; 3],
}
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteBinding {
    source: Shared<Drive>,
    vertices: FiniteVertices,
    carriers: [Vec<CubicIntersectionSpan>; 12],
    nominal_binding_error_mm: f64,
}
fn distance(actual: [f64; 3], expected: [I; 3]) -> f64 {
    actual
        .into_iter()
        .zip(expected)
        .fold(I::point(0.), |sum, (p, e)| sum + (I::point(p) - e).square())
        .sqrt()
        .hi
}
fn cross(a: [I; 2], b: [I; 2]) -> I {
    a[0] * b[1] - a[1] * b[0]
}
fn angular_advance(span: &CubicIntersectionSpan) -> Result<bool, Error> {
    let p = span.poles_mm.map(|p| [I::point(p[0]), I::point(p[1])]);
    // A strict minor angular sector keeps the complete Bernstein hull away
    // from the cone axis; a positive angular derivative prevents retracing.
    if cross(p[0], p[3]).lo <= 0.
        || p[1..3]
            .iter()
            .any(|&q| cross(p[0], q).lo <= 0. || cross(q, p[3]).lo <= 0.)
    {
        return Err(Error::Unresolved);
    }
    let choose3 = [1., 3., 3., 1.];
    let choose2 = [1., 2., 1.];
    let choose5 = [1., 5., 10., 10., 5., 1.];
    for (degree, denominator) in choose5.into_iter().enumerate() {
        let mut coefficient = I::point(0.);
        for i in degree.saturating_sub(2)..=degree.min(3) {
            let j = degree - i;
            let d = std::array::from_fn(|k| p[j + 1][k] - p[j][k]);
            let weight = I::point(choose3[i]) * I::point(choose2[j]) / I::point(denominator);
            coefficient = coefficient + cross(p[i], d) * weight;
        }
        if coefficient.lo <= 0. {
            return Err(Error::Unresolved);
        }
    }
    Ok(span.poles_mm[0][1] <= 0. && span.poles_mm[3][1] > 0.)
}
impl FiniteBinding {
    /// Preserve literal finite controls/vertices with a computed support-binding
    /// bound. The requested limit admits/refuses; it never creates uncertainty.
    /// This is a source-graph certificate, not an independently decoded-file claim.
    pub fn new(
        source: Drive,
        vertices: FiniteVertices,
        carriers: [Vec<CubicIntersectionSpan>; 12],
        maximum_error_mm: f64,
    ) -> Result<Self, Error> {
        Self::new_with_cancel(source, vertices, carriers, maximum_error_mm, &|| false)
    }
    pub fn new_with_cancel(
        source: Drive,
        vertices: FiniteVertices,
        carriers: [Vec<CubicIntersectionSpan>; 12],
        maximum_error_mm: f64,
        cancelled: &dyn Fn() -> bool,
    ) -> Result<Self, Error> {
        let values = vertices
            .mouth_mm
            .into_iter()
            .chain(vertices.floor_mm)
            .chain([vertices.apex_mm])
            .flatten();
        if !maximum_error_mm.is_finite()
            || maximum_error_mm <= 0.
            || !values.clone().all(f64::is_finite)
        {
            return Err(Error::Invalid);
        }
        if values.map(f64::abs).any(|v| v > 1e9)
            || carriers.iter().map(Vec::len).sum::<usize>() > 1_000_000
        {
            return Err(Error::Limit);
        }
        let mut crossings = 0;
        let mut error = distance(
            vertices.apex_mm,
            [
                I::point(0.),
                I::point(0.),
                I::point(source.cone().tip_depth_mm),
            ],
        );
        for i in 0..12 {
            if cancelled() {
                return Err(Error::Cancelled);
            }
            let next = (i + 1) % 12;
            let mouth = vertices.mouth_mm[i];
            let floor = vertices.floor_mm[i];
            if mouth[..2] != source.arcs()[i].start_mm
                || mouth == floor
                || floor == vertices.apex_mm
                || carriers[i].is_empty()
            {
                return Err(Error::Invalid);
            }
            let first = &carriers[i][0];
            let last = &carriers[i][carriers[i].len() - 1];
            if first.parameters[0] != 0.
                || last.parameters[1] != 1.
                || first.poles_mm[0] != floor
                || last.poles_mm[3] != vertices.floor_mm[next]
            {
                return Err(Error::Invalid);
            }
            let (contact, _) = source.intersection_state(i, 0.);
            error = error.max(distance(mouth, [contact[0], contact[1], I::point(0.)]));
            for pair in carriers[i].windows(2) {
                if pair[0].parameters[1] != pair[1].parameters[0]
                    || pair[0].poles_mm[3] != pair[1].poles_mm[0]
                {
                    return Err(Error::Invalid);
                }
            }
            for span in &carriers[i] {
                if cancelled() {
                    return Err(Error::Cancelled);
                }
                error = error.max(
                    source
                        .certify_intersection_span(i, span)?
                        .nominal_curve_error_mm,
                );
                crossings += usize::from(angular_advance(span)?);
            }
        }
        if crossings != 1 {
            return Err(Error::Invalid);
        }
        let total = (I::point(error) + I::point(source.spatial_source_error_bound_mm())).hi;
        if total > maximum_error_mm {
            return Err(Error::Unresolved);
        }
        Ok(Self {
            source: Shared::new(source),
            vertices,
            carriers,
            nominal_binding_error_mm: error,
        })
    }
    pub fn source(&self) -> &Drive {
        &self.source
    }
    pub(crate) fn source_handle(&self) -> Shared<Drive> {
        Shared::clone(&self.source)
    }
    pub fn vertices(&self) -> &FiniteVertices {
        &self.vertices
    }
    pub fn carriers(&self) -> &[Vec<CubicIntersectionSpan>; 12] {
        &self.carriers
    }
    /// Additional retained-boundary approximation, separate from transport.
    pub fn nominal_binding_error_mm(&self) -> f64 {
        self.nominal_binding_error_mm
    }
    pub fn spatial_source_error_bound_mm(&self) -> f64 {
        (I::point(self.source.spatial_source_error_bound_mm())
            + I::point(self.nominal_binding_error_mm))
        .hi
    }
    /// Extra trim tube only; the source's volume certificate already includes
    /// its own numerical/transport reserve.
    pub fn additional_volume_error_bound_mm3(&self) -> Option<f64> {
        let d = self.nominal_binding_error_mm;
        let floor =
            ((I::point(1.) + I::point(self.source.cone().depth_per_radius)) * I::point(d)).hi;
        self.source.boundary_volume_error_bound_mm3(d, floor)
    }
}
