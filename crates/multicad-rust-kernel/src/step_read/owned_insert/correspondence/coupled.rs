//! Exact source-owned numerator correspondence bound to actual material/map.
//! Decoded volume remains the separate literal-flux authority.
use super::*;
use crate::step::insert::{
    coupled::{Coefficient, topology::CoupledTopology},
    exact_cap,
    resources::Work,
};
pub struct CoupledSourceCorrespondence<'a, 'i> {
    material: &'a MaterialBounds<'a, 'a, 'i>,
    expected: &'a CoupledTopology<'a, 'a>,
    mapping: &'a PhysicalMapping<'a, 'i>,
    position: f64,
    derivative: f64,
    ruling: f64,
    work: usize,
    bytes: usize,
}
impl<'a, 'i> CoupledSourceCorrespondence<'a, 'i> {
    pub fn actual(&self) -> &SupportBindings<'a, 'i> {
        self.material.binding()
    }
    pub fn material(&self) -> &MaterialBounds<'a, 'a, 'i> {
        self.material
    }
    pub fn expected(&self) -> &CoupledTopology<'a, 'a> {
        self.expected
    }
    pub fn mapping(&self) -> &PhysicalMapping<'a, 'i> {
        self.mapping
    }
    pub fn coefficient_unit_mm(&self) -> f64 {
        self.expected.plan().coefficient_unit_mm()
    }
    /// Global local-physical-mm C0, including filled planar annuli. Actual
    /// normalized mapping axes/origin versus Native's frame are a later proof.
    pub fn position_error_mm(&self) -> f64 {
        self.position
    }
    /// Original WHOLE source-curve parameter, not the tiny clipped file span.
    pub fn row_derivative_error_mm(&self) -> f64 {
        self.derivative
    }
    /// Normalized side-chart affine ruling only; no cap-interior C1 is claimed.
    pub fn ruling_derivative_error_mm(&self) -> f64 {
        self.ruling
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.bytes
    }
}
fn equal(
    actual: Scalar<'_>,
    sign: i8,
    numerator: &[u8],
    bits: usize,
    memory: &mut Memory,
    w: &mut Ledger<'_>,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<bool> {
    let remaining = w.remaining.min(16_000_000);
    if remaining == 0 {
        return Err(Error::Limit);
    }
    let mut numeric = Work::new(remaining, interrupt).map_err(Error::from)?;
    let result = exact_cap::equals_dyadic(
        actual.spelling(),
        sign,
        numerator,
        bits,
        memory,
        &mut numeric,
    );
    // Charge the bounded child view immediately; it never renews caller work.
    w.step(numeric.used)?;
    result.map_err(Error::from)
}
fn point_equal(
    actual: &PointLeaf<'_>,
    wanted: &[Coefficient; 3],
    memory: &mut Memory,
    w: &mut Ledger<'_>,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<()> {
    for (actual, wanted) in actual.coordinates().iter().zip(wanted) {
        if !equal(
            *actual,
            wanted.signum(),
            wanted.numerator_magnitude_le(),
            0,
            memory,
            w,
            interrupt,
        )? {
            return Err(Error::Unsupported);
        }
    }
    Ok(())
}
impl<'b, 'c, 'i> MaterialBounds<'b, 'c, 'i> {
    /// Source-known equivalence of this actual CLOSED material, including its
    /// selected representation operator. This constructs no source from a file
    /// header and exposes no inferred/nominal decoded mass.
    pub fn coupled_source_correspondence<'a>(
        &'a self,
        expected: &'a CoupledTopology<'a, 'a>,
        mapping: &'a PhysicalMapping<'a, 'i>,
        options: CorrespondenceOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<CoupledSourceCorrespondence<'a, 'i>> {
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
        let actual = self.binding();
        let p = actual.source();
        if !std::ptr::eq(mapping.parsed(), p) {
            return Err(Error::Invalid);
        }
        if p.mapping_link().map(|m| m.operator_entity()) != Some(mapping.operator_entity())
            || mapping.handedness()? != 1
        {
            return Err(Error::Unsupported);
        }
        let mut w = Ledger {
            used: 0,
            remaining: options.maximum_work,
            interrupt,
        };
        // Material includes parse/C0 once; mapping includes that SAME parse.
        w.step(self.work_used() + expected.work_used() + mapping.work_used() - p.work_used())?;
        let fixed = self
            .retained_budget_bytes()
            .saturating_add(expected.retained_bytes())
            .saturating_add(mapping.query_budget_bytes() - p.retained_budget_bytes())
            .saturating_add(size_of::<CoupledSourceCorrespondence<'_, '_>>());
        let mut memory = Memory::new(fixed, options.retained_budget_bytes).map_err(Error::from)?;
        // Metadata allocation participates in the SAME64M ledger, rather than
        // a fresh16M proof balance.
        let request = expected
            .vertices()
            .len()
            .saturating_mul(size_of::<Option<u32>>());
        memory.add(request).map_err(Error::from)?;
        memory.used = fixed;
        w.step(expected.vertices().len() + 1)?;
        let mut vertices = Vec::<Option<u32>>::new();
        vertices
            .try_reserve_exact(expected.vertices().len())
            .map_err(|_| Error::Limit)?;
        let mut vertices = memory.account(vertices).map_err(Error::from)?;
        vertices.resize(expected.vertices().len(), None);
        let bytes = memory.used.saturating_add(exact_cap::SCRATCH);
        // Preflight numeric scratch before inspecting any original coefficient.
        memory.add(exact_cap::SCRATCH).map_err(Error::from)?;
        memory.used -= exact_cap::SCRATCH;
        let records = &p.entities[&mapping.operator_entity()];
        let base = records
            .iter()
            .find(|e| e.name == "CARTESIAN_TRANSFORMATION_OPERATOR")
            .expect("earned physical operator");
        if !matches!(base.args[3], Value::Number(_)) {
            return Err(Error::Unsupported);
        }
        let mut numeric = Numeric::new(base, &mut w)?;
        let scale = numeric.scalar(&base.args[3])?;
        numeric.finish()?;
        let bits = expected.plan().rows()[0].coefficients()[0][0].fractional_bits();
        if !equal(scale, 1, &[1], bits, &mut memory, &mut w, interrupt)? {
            return Err(Error::Unsupported);
        }
        if p.faces().len() != expected.faces().len()
            || p.edges().len() != expected.edges().len()
            || p.vertices().len() != expected.used_vertices()
        {
            return Err(Error::Unsupported);
        }
        for (source, actual) in expected.edges().iter().zip(p.edges()) {
            w.step(64)?;
            let curve = p.curve(actual.curve_entity()).expect("finite bound edge");
            let degree = if matches!(source.geometry, EdgeGeometry::Row(_)) {
                11
            } else {
                1
            };
            // Parsed earns normalized exact[0,1] knots and
            // poles=(knots.len-1)*degree+1. Bindings earn2-pole rulings.
            // Thus a12-pole row already earns the one fullBezier span.
            if curve.degree() != degree || (degree == 11 && curve.pole_entities().len() != 12) {
                return Err(Error::Unsupported);
            }
            for (at, id) in source.ends.into_iter().zip(actual.vertex_entities()) {
                if let Some(old) = vertices[at] {
                    if old != id {
                        return Err(Error::Unsupported);
                    }
                } else {
                    vertices[at] = Some(id);
                    w.step(32)?;
                    let vertex = p.vertex(id).expect("bound physical vertex");
                    let point = p.point(vertex.point_entity).expect("bound vertex geometry");
                    let [row, control] = expected.vertices()[at].coefficient_reference();
                    point_equal(
                        point,
                        &expected.plan().rows()[row].coefficients()[control],
                        &mut memory,
                        &mut w,
                        interrupt,
                    )?;
                }
            }
        }
        // Actual material admission makes every actual vertex occupied by an
        // edge. Matched complete endpoints+equal counts earn the vertex bijection.
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
        }
        let mut boundary = 0_f64;
        let mut derivative = 0_f64;
        for (index, row) in expected.plan().rows().iter().enumerate() {
            let curve = p
                .curve(p.edges()[index].curve_entity())
                .expect("bound canonical source row");
            for (&id, wanted) in curve.pole_entities().iter().zip(row.coefficients()) {
                w.step(32)?;
                point_equal(
                    p.point(id).expect("finite row pole"),
                    wanted,
                    &mut memory,
                    &mut w,
                    interrupt,
                )?;
            }
            boundary = boundary.max(row.position_error_mm());
            derivative = derivative.max(row.original_curve_derivative_error_mm());
        }
        // Degree-one actual u and exact shared row/PCurve identities bind the
        // entire side field. Caps: winding is unchanged outside theδ boundary
        // tube; nearest boundary transfer gives planar2δ plus plane offsetδ.
        let position = magnitude(I::point(3.) * I::point(boundary))?;
        let ruling = sum(boundary, boundary)?;
        if position > options.position_limit_mm || derivative > options.derivative_limit_mm {
            return Err(Error::Unsupported);
        }
        w.step(1)?;
        Ok(CoupledSourceCorrespondence {
            material: self,
            expected,
            mapping,
            position,
            derivative,
            ruling,
            work: w.used,
            bytes,
        })
    }
}
