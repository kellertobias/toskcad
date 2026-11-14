//! Explicit derived edge expressions, never exact raw-support intersections.
//! This receipt supplies no face filling, embedding, material or mass authority.
use super::*;
mod basis;
mod carrier;
use carrier::Carrier;
#[cfg(test)]
mod tests;
// A clamped cubic has <=4096 spans. Generic raw numeric-leaf lookup can hold
// its token pointers, parsed-number buffer and leaf descriptors together.
// Include geometric Vec growth and a second leaf-buffer reserve, before any
// such query. Original text/AST storage belongs to the separately charged Body.
const MAX_SPLINE_LEAVES: usize = 1 + 2 * 4097;
const QUERY_SCRATCH: usize = 8192
    + MAX_SPLINE_LEAVES
        * (2 * std::mem::size_of::<&str>()
            + 2 * std::mem::size_of::<f64>()
            + 2 * std::mem::size_of::<OriginalNumericLeaf<'static>>());
#[derive(Clone, Copy)]
pub struct SewnOptions {
    pub maximum_position_transfer_mm: f64,
    /// Per unit whole-edge parameter for analytic curves, and per retained
    /// Bezier span for cubics. Parameters follow EDGE_CURVE start->end; a
    /// sameSense=false span reverses the original local parameter and derivative.
    /// This is not an arc-length or original global-knot derivative.
    pub maximum_derivative_transfer_mm: f64,
    pub maximum_work: usize,
    pub maximum_bytes: usize,
}
impl Default for SewnOptions {
    fn default() -> Self {
        Self {
            maximum_position_transfer_mm: 1e-7,
            maximum_derivative_transfer_mm: 1e-6,
            maximum_work: WORK,
            maximum_bytes: BYTES,
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SewnParameter {
    UnitAnalyticEdge,
    UnitBezierSpan {
        original_span: usize,
        raw_curve_forward: bool,
    },
}
#[derive(Clone, Copy)]
struct Bezier {
    poles: [[I; 3]; 4],
    references: [u32; 2],
    original_span: usize,
    raw_curve_forward: bool,
}
enum Expression {
    // Exact mathematical source E(s)=gamma(s)+(1-s)*(V0-gamma(0))+
    // s*(V1-gamma(1)), with Vi the original Body POINT leaves. Numeric fields
    // below are evaluation caches of those coupled expressions, never freely
    // selectable replacement coefficients. Endpoints cancel symbolically.
    EndpointCorrected(Box<Corrected>),
    OriginalCubics(Vec<Bezier>),
}
struct Corrected {
    carrier: Carrier,
    correction: [[I; 3]; 2],
    endpoint_positions: [[I; 3]; 2],
}
fn corrected_fields(value: &Corrected, t: I, envelope: bool) -> ([I; 3], [I; 3]) {
    let (raw, derivative) = if envelope {
        value.carrier.domain_envelope()
    } else {
        value.carrier.evaluate(t)
    };
    (
        std::array::from_fn(|i| {
            raw[i] + (I::point(1.) - t) * value.correction[0][i] + t * value.correction[1][i]
        }),
        std::array::from_fn(|i| derivative[i] + value.correction[1][i] - value.correction[0][i]),
    )
}
fn cubic_derivative(poles: [[I; 3]; 4], t: I) -> [I; 3] {
    casteljau(
        std::array::from_fn::<_, 3, _>(|i| {
            std::array::from_fn(|k| I::point(3.) * (poles[i + 1][k] - poles[i][k]))
        }),
        t,
    )
}
impl Expression {
    fn finite_unit_domain(&self, work: &mut Budget<'_>) -> BranchResult<bool> {
        let unit = I::new(0., 1.);
        let mut valid = true;
        match self {
            Self::EndpointCorrected(value) => {
                work.charge(128)?;
                let (point, derivative) = corrected_fields(value, unit, true);
                valid &= point
                    .into_iter()
                    .chain(derivative)
                    .chain(value.endpoint_positions.into_iter().flatten())
                    .all(I::valid);
            }
            Self::OriginalCubics(spans) => {
                for span in spans {
                    work.charge(128)?;
                    let point = casteljau(span.poles, unit);
                    let derivative = cubic_derivative(span.poles, unit);
                    valid &= point
                        .into_iter()
                        .chain(derivative)
                        .chain(span.poles.into_iter().flatten())
                        .all(I::valid);
                }
            }
        }
        // Runtime operations are only inclusion-isotone directed add/multiply
        // over unit subintervals, plus bounded sin/cos already widened to the
        // complete[-1,1] box above. The retained finite lift bounds also keep
        // theta itself finite. Endpoint overrides use validated original leaves.
        // Thus every runtime intermediate is bounded before source publication;
        // no division or fresh source decoding occurs during evaluation.
        Ok(valid)
    }
}
pub struct SewnEdgeSource {
    edge: usize,
    endpoints: [u32; 2],
    expression: Expression,
    position_transfer: f64,
    derivative_transfer: f64,
    owner_distances: [f64; 2],
}
impl SewnEdgeSource {
    pub fn edge_index(&self) -> usize {
        self.edge
    }
    pub fn source_endpoint_point_references(&self) -> [u32; 2] {
        self.endpoints
    }
    pub fn position_transfer_mm(&self) -> f64 {
        self.position_transfer
    }
    pub fn derivative_transfer_mm(&self) -> f64 {
        self.derivative_transfer
    }
    pub fn raw_owner_distances_mm(&self) -> [f64; 2] {
        self.owner_distances
    }
    pub fn span_count(&self) -> usize {
        match &self.expression {
            Expression::EndpointCorrected(_) => 1,
            Expression::OriginalCubics(spans) => spans.len(),
        }
    }
}
/// All data-derived curves share their ORIGINAL immutable POINT leaf at each
/// declared vertex. Endpoint corrections change the raw carriers explicitly;
/// they do not assert those raw curves or supports already have common roots.
pub struct SewnEdgeSources<'a, 'input> {
    body: &'a Body<'input>,
    lifts: &'a FiniteLifts<'a>,
    cubics: &'a CubicSupportBounds<'a>,
    edges: Vec<SewnEdgeSource>,
    work_used: usize,
}
#[derive(Clone, Copy, Debug)]
pub struct SewnEvaluation {
    pub position_bounds_mm: [[f64; 2]; 3],
    pub derivative_bounds_mm: [[f64; 2]; 3],
    pub parameterization: SewnParameter,
    /// An exact original coordinate entity, not a proximity identity. Only
    /// meaningful together with this receipt's borrowed arena identity.
    pub original_point_reference: Option<u32>,
}
fn norm(p: [I; 3]) -> I {
    p.into_iter()
        .fold(I::point(0.), |s, x| s + x.square())
        .sqrt()
}
fn difference(a: [I; 3], b: [I; 3]) -> [I; 3] {
    std::array::from_fn(|i| a[i] - b[i])
}
fn source_point(body: &Body<'_>, point: u32, work: &mut Budget<'_>) -> BranchResult<[I; 3]> {
    work.charge(1)?;
    let leaf = body.exact_point_leaf(point, work.remaining, work.interrupt)?;
    work.charge(leaf.work_used())?;
    Ok(leaf.coordinate_bounds_mm().map(|x| I::new(x[0], x[1])))
}
fn footprint(
    base: usize,
    edges: usize,
    spans: usize,
    limit: usize,
    work: &mut Budget<'_>,
) -> BranchResult<()> {
    work.charge(1)?;
    let bytes = add(
        base,
        add(
            edges.saturating_mul(std::mem::size_of::<SewnEdgeSource>()),
            spans.saturating_mul(std::mem::size_of::<Bezier>()),
        ),
    );
    if bytes > limit {
        Err(PhillipsReadError::Limit)
    } else {
        Ok(())
    }
}
impl<'source> FiniteLifts<'source> {
    pub fn sewn_edge_sources<'a, 'input>(
        &'a self,
        body: &'a Body<'input>,
        cubics: &'a CubicSupportBounds<'source>,
        options: SewnOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> BranchResult<SewnEdgeSources<'a, 'input>>
    where
        'source: 'a,
    {
        if !options.maximum_position_transfer_mm.is_finite()
            || options.maximum_position_transfer_mm <= 0.
            || !options.maximum_derivative_transfer_mm.is_finite()
            || options.maximum_derivative_transfer_mm <= 0.
            || options.maximum_work == 0
            || options.maximum_work > WORK
            || options.maximum_bytes == 0
            || options.maximum_bytes > BYTES
            || !std::ptr::eq(self.boundary(), cubics.boundary())
            || !std::ptr::eq(self.boundary(), body.cavity())
        {
            return Err(PhillipsReadError::Invalid);
        }
        if interrupt.is_some_and(|f| f()) {
            return Err(PhillipsReadError::Cancelled);
        }
        let inherited = add(self.work_used(), cubics.work_used());
        let mut work = Budget {
            remaining: options
                .maximum_work
                .checked_sub(inherited)
                .ok_or(PhillipsReadError::Limit)?,
            interrupt,
        };
        let b = self.boundary();
        let base = add(
            add(
                body.retained_budget_bytes(),
                self.edges()
                    .len()
                    .saturating_mul(std::mem::size_of::<Corrected>()),
            ),
            add(
                self.retained_capacity_bytes(),
                add(
                    cubics.retained_capacity_bytes(),
                    add(
                        std::mem::size_of::<SewnEdgeSources<'_, '_>>(),
                        QUERY_SCRATCH,
                    ),
                ),
            ),
        );
        let mut remaining_spans = b
            .edges()
            .iter()
            .map(|e| {
                if let Curve::CubicIntersection { spans, .. } = &e.curve {
                    spans.len()
                } else {
                    0
                }
            })
            .sum::<usize>();
        footprint(
            base,
            b.edges().len(),
            remaining_spans,
            options.maximum_bytes,
            &mut work,
        )?;
        let mut edges = Vec::new();
        edges
            .try_reserve_exact(b.edges().len())
            .map_err(|_| PhillipsReadError::Limit)?;
        footprint(
            base,
            edges.capacity(),
            remaining_spans,
            options.maximum_bytes,
            &mut work,
        )?;
        let mut retained_spans = 0usize;
        for (edge_index, edge) in b.edges().iter().enumerate() {
            work.charge(2048)?;
            let endpoints = edge.vertices.map(|i| b.vertices()[i].point_entity);
            let (expression, position_transfer, derivative_transfer, raw) =
                if let Curve::CubicIntersection {
                    spans,
                    pole_entities,
                    spline_entity,
                    ..
                } = &edge.curve
                {
                    basis::prove(body, *spline_entity, spans.len(), &mut work)?;
                    work.charge(spans.len().saturating_mul(16))?;
                    let mut prepared = Vec::new();
                    prepared
                        .try_reserve_exact(spans.len())
                        .map_err(|_| PhillipsReadError::Limit)?;
                    remaining_spans -= spans.len();
                    retained_spans = add(retained_spans, prepared.capacity());
                    footprint(
                        base,
                        edges.capacity(),
                        add(retained_spans, remaining_spans),
                        options.maximum_bytes,
                        &mut work,
                    )?;
                    // The sealed physical parser proves the actual B-spline pole
                    // references alias the declared endpoint POINTs, including
                    // sameSense=false, and adjacent spans share one pole record.
                    // No equality of merely rounded coordinates is used here.
                    for ordinal in 0..spans.len() {
                        work.charge(16)?;
                        let i = if edge.same_sense {
                            ordinal
                        } else {
                            spans.len() - 1 - ordinal
                        };
                        let mut poles = [[I::point(0.); 3]; 4];
                        for (k, p) in poles.iter_mut().enumerate() {
                            *p = source_point(body, pole_entities[3 * i + k], &mut work)?;
                        }
                        let mut references = [pole_entities[3 * i], pole_entities[3 * i + 3]];
                        if !edge.same_sense {
                            poles.reverse();
                            references.reverse();
                        }
                        prepared.push(Bezier {
                            poles,
                            references,
                            original_span: i,
                            raw_curve_forward: edge.same_sense,
                        });
                    }
                    let raw = cubics
                        .edges()
                        .iter()
                        .find(|p| p.edge_index == edge_index)
                        .unwrap()
                        .owner_distances_mm;
                    (
                        Expression::OriginalCubics(prepared),
                        I::point(0.),
                        I::point(0.),
                        raw,
                    )
                } else {
                    let lift = self
                        .edges()
                        .iter()
                        .find(|p| p.edge_index() == edge_index)
                        .unwrap();
                    let carrier = Carrier::new(lift.source);
                    let endpoint_positions = [
                        source_point(body, endpoints[0], &mut work)?,
                        source_point(body, endpoints[1], &mut work)?,
                    ];
                    let correction = [0, 1].map(|i| {
                        difference(
                            endpoint_positions[i],
                            carrier.evaluate(I::point(i as f64)).0,
                        )
                    });
                    let endpoint_errors = correction.map(norm);
                    let position = I::point(endpoint_errors[0].hi.max(endpoint_errors[1].hi));
                    let derivative = norm(difference(correction[1], correction[0]));
                    let transfer = carrier.transfer();
                    let raw = self
                        .analytic
                        .edges()
                        .iter()
                        .find(|p| p.edge_index == edge_index)
                        .unwrap()
                        .owner_distances_mm;
                    work.charge(1)?;
                    (
                        Expression::EndpointCorrected(Box::new(Corrected {
                            carrier,
                            correction,
                            endpoint_positions,
                        })),
                        position + transfer[0],
                        derivative + transfer[1],
                        raw,
                    )
                };
            let owner_distances = raw.map(|x| (I::point(x) + position_transfer).hi);
            if ![position_transfer, derivative_transfer]
                .into_iter()
                .all(I::valid)
                || owner_distances.into_iter().any(|x| !x.is_finite())
                || position_transfer.hi > options.maximum_position_transfer_mm
                || derivative_transfer.hi > options.maximum_derivative_transfer_mm
                || !expression.finite_unit_domain(&mut work)?
            {
                return Err(PhillipsReadError::Unsupported);
            }
            edges.push(SewnEdgeSource {
                edge: edge_index,
                endpoints,
                expression,
                position_transfer: position_transfer.hi,
                derivative_transfer: derivative_transfer.hi,
                owner_distances,
            });
        }
        work.charge(1)?;
        Ok(SewnEdgeSources {
            body,
            lifts: self,
            cubics,
            edges,
            work_used: options.maximum_work - work.remaining,
        })
    }
}
fn casteljau<const N: usize>(mut poles: [[I; 3]; N], t: I) -> [I; 3] {
    for remaining in (1..N).rev() {
        for i in 0..remaining {
            poles[i] =
                std::array::from_fn(|k| (I::point(1.) - t) * poles[i][k] + t * poles[i + 1][k]);
        }
    }
    poles[0]
}
impl<'a, 'input> SewnEdgeSources<'a, 'input> {
    pub fn body(&self) -> &'a Body<'input> {
        self.body
    }
    pub fn boundary(&self) -> &'a Boundary {
        self.lifts.boundary()
    }
    pub fn edges(&self) -> &[SewnEdgeSource] {
        &self.edges
    }
    pub fn work_used(&self) -> usize {
        self.work_used
    }
    pub fn retained_capacity_bytes(&self) -> usize {
        let spans = self
            .edges
            .iter()
            .map(|e| match &e.expression {
                Expression::OriginalCubics(v) => allocation(v),
                Expression::EndpointCorrected(value) => std::mem::size_of_val(value.as_ref()),
            })
            .fold(0, add);
        add(
            std::mem::size_of::<Self>(),
            add(
                allocation(&self.edges),
                add(
                    spans,
                    add(
                        self.lifts.retained_capacity_bytes(),
                        self.cubics.retained_capacity_bytes(),
                    ),
                ),
            ),
        )
    }
    /// Bounded constant work per source query. No triangulation or projection
    /// is performed. The derivative is with respect to the returned unit chart.
    pub fn evaluate(
        &self,
        edge: usize,
        span: usize,
        parameter: [f64; 2],
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> BranchResult<SewnEvaluation> {
        if interrupt.is_some_and(|f| f()) {
            return Err(PhillipsReadError::Cancelled);
        }
        let t = I::new(parameter[0], parameter[1]);
        let record = self.edges.get(edge).ok_or(PhillipsReadError::Invalid)?;
        if !t.valid() || t.lo < 0. || t.hi > 1. || span >= record.span_count() {
            return Err(PhillipsReadError::Invalid);
        }
        let endpoint = if t.lo == t.hi && (t.lo == 0. || t.lo == 1.) {
            Some(t.lo as usize)
        } else {
            None
        };
        let (point, derivative, parameterization, reference) = match &record.expression {
            Expression::EndpointCorrected(value) => {
                let (point, derivative) = corrected_fields(value, t, false);
                let point = if let Some(i) = endpoint {
                    value.endpoint_positions[i]
                } else {
                    point
                };
                (
                    point,
                    derivative,
                    SewnParameter::UnitAnalyticEdge,
                    endpoint.map(|i| record.endpoints[i]),
                )
            }
            Expression::OriginalCubics(spans) => {
                let prepared = spans[span];
                let point = if let Some(i) = endpoint {
                    prepared.poles[3 * i]
                } else {
                    casteljau(prepared.poles, t)
                };
                let derivative = cubic_derivative(prepared.poles, t);
                (
                    point,
                    derivative,
                    SewnParameter::UnitBezierSpan {
                        original_span: prepared.original_span,
                        raw_curve_forward: prepared.raw_curve_forward,
                    },
                    endpoint.map(|i| prepared.references[i]),
                )
            }
        };
        if interrupt.is_some_and(|f| f()) {
            return Err(PhillipsReadError::Cancelled);
        }
        Ok(SewnEvaluation {
            position_bounds_mm: point.map(|x| [x.lo, x.hi]),
            derivative_bounds_mm: derivative.map(|x| [x.lo, x.hi]),
            parameterization,
            original_point_reference: reference,
        })
    }
}
