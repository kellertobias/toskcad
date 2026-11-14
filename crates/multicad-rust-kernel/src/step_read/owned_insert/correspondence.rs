//! Source-known finite geometry correspondence, never decoded mass/occupancy.
use super::*;
use crate::step::insert::resources::Memory;
use crate::{
    step::insert::{EdgeGeometry, FaceGeometry, TopologyPlan},
    tapered_thread::interval::I,
};
mod coupled;
pub use coupled::CoupledSourceCorrespondence;
#[derive(Clone, Copy, Debug)]
pub struct CorrespondenceOptions {
    pub position_limit_mm: f64,
    /// Absolute derivative in the original whole-row parameter t in[0,1].
    pub derivative_limit_mm: f64,
    /// Aggregate query including both independently earned prerequisites.
    pub maximum_work: usize,
    pub retained_budget_bytes: usize,
}
impl Default for CorrespondenceOptions {
    fn default() -> Self {
        Self {
            position_limit_mm: 1e-7,
            derivative_limit_mm: 1e-6,
            maximum_work: 64_000_000,
            retained_budget_bytes: 64 * 1024 * 1024,
        }
    }
}
pub struct SourceCorrespondence<'a, 'input> {
    actual: &'a SupportBindings<'a, 'input>,
    expected: &'a TopologyPlan<'a, 'a, 'a>,
    position: f64,
    derivative: f64,
    ruling_derivative: f64,
    work: usize,
    bytes: usize,
}
impl<'a, 'i> SourceCorrespondence<'a, 'i> {
    pub fn actual(&self) -> &SupportBindings<'a, 'i> {
        self.actual
    }
    pub fn expected(&self) -> &TopologyPlan<'a, 'a, 'a> {
        self.expected
    }
    pub fn position_error_mm(&self) -> f64 {
        self.position
    }
    pub fn row_derivative_error_mm(&self) -> f64 {
        self.derivative
    }
    /// Absolute derivative error along the normalized affine ruling u.
    pub fn ruling_derivative_error_mm(&self) -> f64 {
        self.ruling_derivative
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.bytes
    }
}
impl<'b, 'i> SupportBindings<'b, 'i> {
    /// Match the actual finite member to this owned source using its canonical
    /// emitted edge/face order, with measured original-decimal coefficient error.
    /// A reordered foreign member can earn MaterialBounds independently while
    /// remaining unsupported here. No source gauge/header produces file mass.
    pub fn source_correspondence<'a>(
        &'a self,
        expected: &'a TopologyPlan<'a, 'a, 'a>,
        options: CorrespondenceOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<SourceCorrespondence<'a, 'i>> {
        if !options.position_limit_mm.is_finite()
            || options.position_limit_mm <= 0.
            || !options.derivative_limit_mm.is_finite()
            || options.derivative_limit_mm <= 0.
            || options.maximum_work == 0
            || options.maximum_work > 64_000_000
            || options.retained_budget_bytes == 0
            || options.retained_budget_bytes > 64 * 1024 * 1024
        {
            return Err(Error::Invalid);
        }
        let mut w = Ledger {
            used: 0,
            remaining: options.maximum_work,
            interrupt,
        };
        w.step(self.work_used() + expected.work_used())?;
        let p = self.source();
        let s = expected.surface();
        if p.faces().len() != expected.faces().len()
            || p.edges().len() != expected.edges().len()
            || p.vertices().len() != expected.used_vertices()
        {
            return Err(Error::Unsupported);
        }
        let fixed = self
            .retained_budget_bytes()
            .saturating_add(expected.retained_bytes())
            .saturating_add(
                size_of::<SourceCorrespondence<'_, '_>>() + 4096 + size_of::<Vec<Option<u32>>>(),
            );
        let mut memory = Memory::new(fixed, options.retained_budget_bytes).map_err(Error::from)?;
        let estimate = s.vertices().len().saturating_mul(size_of::<Option<u32>>());
        memory.add(estimate).map_err(Error::from)?;
        memory.used = fixed;
        let mut vertices = Vec::new();
        vertices
            .try_reserve_exact(s.vertices().len())
            .map_err(|_| Error::Limit)?;
        let mut vertices = memory.account(vertices).map_err(Error::from)?;
        let mut bytes = memory.used;
        vertices.resize(s.vertices().len(), None);
        let mut position = 0_f64;
        for (source, actual) in expected.edges().iter().zip(p.edges()) {
            w.step(64)?;
            let curve = p
                .curve(actual.curve_entity())
                .expect("finite bound actual curve");
            let degree = match source.geometry {
                EdgeGeometry::Row(_) => 11,
                EdgeGeometry::Ruling => 1,
            };
            if curve.degree() != degree {
                return Err(Error::Unsupported);
            }
            for (at, id) in source.ends.into_iter().zip(actual.vertex_entities()) {
                if let Some(old) = vertices[at] {
                    if old != id {
                        return Err(Error::Unsupported);
                    }
                } else {
                    vertices[at] = Some(id);
                    // Measure exactly the used physical endpoints once, at their
                    // first incident edge. Unused source stations are not body.
                    w.step(32)?;
                    let v = p.vertex(id).expect("finite bound actual vertex");
                    let point = p.point(v.point_entity).expect("finite bound actual point");
                    position = position.max(distance(
                        point.coordinates(),
                        s.vertices()[at].point_mm(),
                        &mut w,
                    )?);
                }
            }
        }
        // The map covers every actual edge endpoint, hence every actual vertex.
        // Equal used-vertex counts make it a bijection of the finite boundaries.
        for (source, actual) in expected.faces().iter().zip(p.faces()) {
            w.step(64)?;
            if source.same_sense != actual.same_sense()
                || source.loops.len() != actual.loops().len()
            {
                return Err(Error::Unsupported);
            }
            for (a, b) in source.loops.iter().zip(actual.loops()) {
                if a.len() != b.coedges().len() {
                    return Err(Error::Unsupported);
                }
                for (a, b) in a.iter().zip(b.coedges()) {
                    w.step(1)?;
                    if p.edges()[a.edge].entity() != b.edge || a.forward != b.forward {
                        return Err(Error::Unsupported);
                    }
                }
            }
            // Both sealed protocols use one loop for a ruled side and two
            // for an annular cap. Equal loop arity above therefore already
            // proves support kind, without a second unreachable shape gate.
            if let FaceGeometry::Cap(end) = source.geometry {
                let height = self.cap_heights()[end as usize];
                let intended = if end == 0 {
                    0.
                } else {
                    s.domain().boundary().source().full_length_mm()
                };
                position = position.max(magnitude(enclosure(height) - I::point(intended))?);
            }
        }
        // Full actual UV bindings plus the matched finite face coedges couple
        // each tensor's two literal boundary rows to these owning row curves.
        // Its entire u-degree-one field is their affine interpolation; no
        // independent interior tensor coefficients can escape the row checks.
        let mut derivative = 0_f64;
        for (index, row) in s.rows().iter().enumerate() {
            // TopologyPlan constructs its owned row edges first, in row order.
            let curve = p
                .curve(p.edges()[index].curve_entity())
                .expect("finite actual row");
            if curve.pole_entities().len() != row.poles().len()
                || curve.knots().len() != row.knots().len()
            {
                return Err(Error::Unsupported);
            }
            for (&literal, &wanted) in curve.knots().iter().zip(row.knots()) {
                w.step(literal.spelling().len() + 32)?;
                if !knot(literal, wanted) {
                    return Err(Error::Unsupported);
                }
            }
            let scratch = row.poles().len().saturating_mul(size_of::<[I; 3]>());
            let base = memory.used;
            memory.add(scratch).map_err(Error::from)?;
            memory.used = base;
            let mut delta = Vec::new();
            delta
                .try_reserve_exact(row.poles().len())
                .map_err(|_| Error::Limit)?;
            let mut delta = memory.account(delta).map_err(Error::from)?;
            bytes = bytes.max(memory.used);
            let mut row_position = 0_f64;
            for (&id, &wanted) in curve.pole_entities().iter().zip(row.poles()) {
                let values = p.point(id).expect("finite row pole").coordinates();
                let mut d = [I::point(0.); 3];
                for axis in 0..3 {
                    w.step(values[axis].spelling().len() + 1)?;
                    d[axis] = enclosure(values[axis]) - I::point(wanted[axis]);
                }
                row_position = row_position.max(norm(d)?);
                delta.push(d);
            }
            let mut row_derivative = 0_f64;
            for span in 0..row.knots().len() - 1 {
                let width = I::point(row.knots()[span + 1]) - I::point(row.knots()[span]);
                let gain = I::point(11.) / width;
                for i in 0..11 {
                    w.step(32)?;
                    let d = std::array::from_fn(|axis| {
                        (delta[11 * span + i + 1][axis] - delta[11 * span + i][axis]) * gain
                    });
                    row_derivative = row_derivative.max(norm(d)?);
                }
            }
            position = position.max(sum(row.position_error_mm(), row_position)?);
            derivative = derivative.max(sum(row.derivative_error_mm(), row_derivative)?);
            memory.release(&delta);
        }
        // The exact u ruling interpolates its two actual curves. Its C0 error
        // is their maximum; its derivative error is bounded by their sum.
        if position > options.position_limit_mm || derivative > options.derivative_limit_mm {
            return Err(Error::Unsupported);
        }
        let ruling_derivative = sum(position, position)?;
        w.step(1)?;
        Ok(SourceCorrespondence {
            actual: self,
            expected,
            position,
            derivative,
            ruling_derivative,
            work: w.used,
            bytes,
        })
    }
}
fn enclosure(s: Scalar<'_>) -> I {
    let [lo, hi] = s.bounds();
    I::new(lo, hi)
}
fn magnitude(i: I) -> Result<f64> {
    if !i.valid() {
        return Err(Error::Unsupported);
    }
    Ok(i.lo.abs().max(i.hi.abs()))
}
fn norm(d: [I; 3]) -> Result<f64> {
    let mut total = I::point(0.);
    for i in d {
        total = total + I::point(magnitude(i)?);
    }
    magnitude(total)
}
fn distance(values: &[Scalar<'_>], wanted: [f64; 3], w: &mut Ledger<'_>) -> Result<f64> {
    let mut d = [I::point(0.); 3];
    for i in 0..3 {
        w.step(values[i].spelling().len() + 1)?;
        d[i] = enclosure(values[i]) - I::point(wanted[i]);
    }
    norm(d)
}
fn sum(a: f64, b: f64) -> Result<f64> {
    magnitude(I::point(a) + I::point(b))
}
fn knot(actual: Scalar<'_>, wanted: f64) -> bool {
    if wanted == 0. {
        return actual.is_zero();
    }
    if wanted == 1. {
        return actual.is_one();
    }
    // SurfacePlan's sealed fitter only doubles the uniform span count, <=1024.
    // Every interior i/2^k therefore has <=10 fractional decimal digits; its
    // canonical 17-significant-digit spelling is its EXACT dyadic value.
    actual.spelling() == format!("{wanted:.16E}")
}

#[cfg(test)]
mod tests;
