//! Original finite endpoint windows for an open, unit-weight spline.
//! This descriptor grants neither curve agreement nor geometric admission.
use super::bezier::{Coordinates, CurveCells};
use super::*;

/// Total local scratch INCLUDING selected-record numeric-token workspace
/// and simultaneous stack/small metadata. Source and the borrowed CurveCells remain owned
/// by the caller and must be included in its simultaneous resident ledger.
/// Selected-record preflight bounds numeric_spellings old+new capacities to
/// 64KiB; another 64KiB covers bounded rational temporaries and fixed arrays.
pub(super) const EXTRA_SCRATCH_BYTES: usize = 128 * 1024;

pub(super) struct Window<'s, 'i, 'e> {
    source: &'s Source<'i>,
    edge: &'e Edge,
    definition: &'e super::super::uv::Definition,
    target_id: u32,
    uv_id: u32,
    uv_domain: [I; 2],
    gamma: f64,
    speed: f64,
    acceleration: f64,
    closest: [f64; 2],
    cached: f64,
    numeric: f64,
}
impl Window<'_, '_, '_> {
    pub(super) fn matches(
        &self,
        source: &Source<'_>,
        target_id: u32,
        definition: &super::super::uv::Definition,
        edge: &Edge,
    ) -> bool {
        std::ptr::eq(self.source, source)
            && std::ptr::eq(self.edge, edge)
            && std::ptr::eq(self.definition, definition)
            && self.target_id == target_id
            && self.uv_id == definition.curve_id
    }
    pub(super) fn gamma(&self) -> f64 {
        self.gamma
    }
    pub(super) fn speed(&self) -> f64 {
        self.speed
    }
    pub(super) fn acceleration(&self) -> f64 {
        self.acceleration
    }
    pub(super) fn closest_delta(&self) -> f64 {
        self.closest[0].max(self.closest[1])
    }
    pub(super) fn closest_delta_sum(&self) -> f64 {
        (I::point(self.closest[0]) + I::point(self.closest[1])).hi
    }
    pub(super) fn cached_window_delta(&self) -> f64 {
        self.cached
    }
    pub(super) fn numeric_reserve(&self) -> f64 {
        self.numeric
    }
    pub(super) fn raw_uv_domain(&self) -> [I; 2] {
        self.uv_domain
    }
    pub(super) fn raw_uv_width(&self) -> I {
        self.uv_domain[1] - self.uv_domain[0]
    }
}
// No token Vec is allocated before counting the complete selected record,
// including unknown components and Typed children. At most 2048 numbers
// give <=4096 borrowed-str capacity slots across Vec growth overlap.
pub(super) fn numeric_preflight(
    source: &Source<'_>,
    id: u32,
    budget: &mut Budget<'_>,
) -> Result<()> {
    fn visit(
        items: &[Value<'_>],
        nodes: &mut usize,
        numbers: &mut usize,
        depth: usize,
        budget: &mut Budget<'_>,
    ) -> Result<()> {
        if depth > 64 {
            return Err(ReadError::Limit);
        }
        for value in items {
            budget.charge(16)?;
            *nodes += 1;
            if *nodes > 8192 {
                return Err(ReadError::Limit);
            }
            match value {
                Value::Number(_) => {
                    *numbers += 1;
                    if *numbers > 2048 {
                        return Err(ReadError::Limit);
                    }
                }
                Value::List(v) => visit(v, nodes, numbers, depth + 1, budget)?,
                Value::Typed(e) => visit(&e.args, nodes, numbers, depth + 1, budget)?,
                _ => {}
            }
        }
        Ok(())
    }
    let mut nodes = 0;
    for e in source.components(id)? {
        budget.charge(32)?;
        let mut numbers = 0;
        visit(&e.args, &mut nodes, &mut numbers, 0, budget)?;
    }
    Ok(())
}
fn valid(x: I) -> Result<I> {
    if x.valid() {
        Ok(x)
    } else {
        Err(ReadError::Unsupported)
    }
}
fn magnitude(x: I) -> Result<f64> {
    Ok(valid(x)?.lo.abs().max(x.hi.abs()))
}
fn knot_domain(
    source: &Source<'_>,
    id: u32,
    curve: &crate::brep::NurbsCurve,
    c1: bool,
    budget: &mut Budget<'_>,
) -> Result<[I; 2]> {
    let field = fields(source, id, 1, budget)?;
    let Some((multiplicities, knots)) = field.knots[0] else {
        return Err(ReadError::Unsupported);
    };
    let m = list(multiplicities)?;
    let k = list(knots)?;
    let (cached, counts) = curve.knots();
    if k.len() != cached.len()
        || m.len() != k.len()
        || k.len() < 2
        || counts.first() != Some(&(curve.degree() + 1))
        || counts.last() != Some(&(curve.degree() + 1))
    {
        return Err(ReadError::Unsupported);
    }
    for (i, count) in counts.iter().copied().enumerate() {
        budget.charge(64)?;
        integer(source, id, &m[i], count, budget)?;
        // An interior multiplicity p would allow a derivative jump, for
        // which no global second-derivative Lipschitz bound is earned.
        if c1 && i > 0 && i + 1 < counts.len() && count > curve.degree().saturating_sub(1) {
            return Err(ReadError::Unsupported);
        }
    }
    Ok([
        scalar(source, id, &k[0], budget)?,
        scalar(source, id, &k[k.len() - 1], budget)?,
    ])
}
fn raw_vertex(source: &Source<'_>, id: usize, budget: &mut Budget<'_>) -> Result<V> {
    budget.charge(64)?;
    let id = u32::try_from(id).map_err(|_| ReadError::Unsupported)?;
    let [_, value] = source.args(id, "VERTEX_POINT")? else {
        return Err(ReadError::Invalid);
    };
    // Ordinary original CARTESIAN_POINT only; a collapsed PCurve vertex
    // needs its separate complete-trace proof, which this leaf cannot infer.
    point(source, reference(value)?, 3, budget)
}
fn bounds(data: &CurveCells<'_, '_, '_>, budget: &mut Budget<'_>) -> Result<(f64, f64, f64)> {
    let p = data.degree();
    let mut derivative = [I::new(f64::INFINITY, f64::NEG_INFINITY); 3];
    let mut speed: f64 = 0.;
    let mut acceleration: f64 = 0.;
    for cell in data.cells() {
        budget.charge(1024)?;
        let h = valid(cell.domain()[1] - cell.domain()[0])?;
        if h.lo <= 0. || cell.controls().len() != p + 1 {
            return Err(ReadError::Unsupported);
        }
        for c in cell.controls() {
            budget.charge(64)?;
            if c[3].lo != 1. || c[3].hi != 1. || c.iter().any(|x| !x.valid()) {
                return Err(ReadError::Unsupported);
            }
        }
        for pair in cell.controls().windows(2) {
            budget.charge(512)?;
            let d = std::array::from_fn(|j| (pair[1][j] - pair[0][j]) * I::point(p as f64) / h);
            speed = speed.max(valid(norm(d))?.hi);
            for j in 0..3 {
                valid(d[j])?;
                derivative[j] =
                    I::new(derivative[j].lo.min(d[j].lo), derivative[j].hi.max(d[j].hi));
            }
        }
        for triple in cell.controls().windows(3) {
            budget.charge(768)?;
            let dd = std::array::from_fn(|j| {
                ((triple[2][j] - triple[1][j]) - (triple[1][j] - triple[0][j]))
                    * I::point((p * (p - 1)) as f64)
                    / (h * h)
            });
            acceleration = acceleration.max(valid(norm(dd))?.hi);
        }
    }
    let gamma = derivative
        .into_iter()
        .map(|x| {
            if x.lo > 0. {
                x.lo
            } else if x.hi < 0. {
                -x.hi
            } else {
                0.
            }
        })
        .fold(0., f64::max);
    if gamma <= 0. || !gamma.is_finite() || !speed.is_finite() || !acceleration.is_finite() {
        return Err(ReadError::Unsupported);
    }
    Ok((gamma, speed, acceleration))
}
fn evaluate(data: &CurveCells<'_, '_, '_>, t: I, budget: &mut Budget<'_>) -> Result<V> {
    let p = data.degree();
    let mut out = [I::new(f64::INFINITY, f64::NEG_INFINITY); 3];
    let mut covered = false;
    for cell in data.cells() {
        budget.charge(128)?;
        let [a, b] = cell.domain();
        if t.hi < a.lo || t.lo > b.hi {
            continue;
        }
        budget.charge(512)?;
        let h = valid(b - a)?;
        if h.lo <= 0. {
            return Err(ReadError::Unsupported);
        }
        let u = valid((t - a) / h)?;
        // Only the part of this cell's true domain participating in t is
        // needed. Every possible intersected cell is included in the hull.
        let u = valid(I::new(u.lo.max(0.), u.hi.min(1.)))?;
        let mut work = [[I::point(0.); 3]; 12];
        for (dst, src) in work.iter_mut().zip(cell.controls()) {
            budget.charge(32)?;
            dst.copy_from_slice(&src[..3]);
        }
        for n in (1..=p).rev() {
            for j in 0..n {
                budget.charge(256)?;
                let next = work[j + 1];
                for (value, next) in work[j].iter_mut().zip(next) {
                    *value = valid(*value * (I::point(1.) - u) + next * u)?;
                }
            }
        }
        for (axis, value) in out.iter_mut().enumerate() {
            *value = I::new(
                value.lo.min(work[0][axis].lo),
                value.hi.max(work[0][axis].hi),
            );
        }
        covered = true;
    }
    if !covered || out.iter().any(|x| !x.valid()) {
        Err(ReadError::Unsupported)
    } else {
        Ok(out)
    }
}

pub(super) fn prepare<'s, 'i, 'e>(
    source: &'s Source<'i>,
    target_id: u32,
    definition: &'e super::super::uv::Definition,
    edge: &'e Edge,
    data: &CurveCells<'_, '_, '_>,
    budget: &mut Budget<'_>,
) -> Result<Window<'s, 'i, 'e>> {
    prepare_inner(source, target_id, definition, edge, data, false, budget)
}
/// Original and decoded natural knot domains are related by their positive
/// affine normalization, rather than identifying their endpoint REALs.
/// Its cached-window discrepancy is expressed back in ORIGINAL parameters;
/// the same normalized-domain numeric reserve then transfers to the body.
pub(super) fn prepare_normalized<'s, 'i, 'e>(
    source: &'s Source<'i>,
    target_id: u32,
    definition: &'e super::super::uv::Definition,
    edge: &'e Edge,
    data: &CurveCells<'_, '_, '_>,
    budget: &mut Budget<'_>,
) -> Result<Window<'s, 'i, 'e>> {
    prepare_inner(source, target_id, definition, edge, data, true, budget)
}
fn prepare_inner<'s, 'i, 'e>(
    source: &'s Source<'i>,
    target_id: u32,
    definition: &'e super::super::uv::Definition,
    edge: &'e Edge,
    data: &CurveCells<'_, '_, '_>,
    normalized: bool,
    budget: &mut Budget<'_>,
) -> Result<Window<'s, 'i, 'e>> {
    budget.charge(4096)?;
    let Curve::Nurbs {
        curve: target,
        closed: false,
    } = &edge.curve
    else {
        return Err(ReadError::Unsupported);
    };
    let Curve::Nurbs {
        curve: uv,
        closed: false,
    } = &definition.geometry
    else {
        return Err(ReadError::Unsupported);
    };
    if !data.matches(source, target_id)
        || !matches!(data.coordinates(), Coordinates::Raw(3))
        || data.degree() != target.degree()
        || !(1..=11).contains(&data.degree())
        || data.cells().is_empty()
        || edge.range.into_iter().any(|x| !x.is_finite())
        || edge.range[0] >= edge.range[1]
    {
        return Err(ReadError::Unsupported);
    }
    if data
        .retained_bytes()
        .checked_add(EXTRA_SCRATCH_BYTES)
        .is_none_or(|n| n > 8 * 1024 * 1024)
    {
        return Err(ReadError::Limit);
    }
    numeric_preflight(source, target_id, budget)?;
    numeric_preflight(source, definition.curve_id, budget)?;
    let natural = knot_domain(source, target_id, target, true, budget)?;
    let f = fields(source, target_id, 1, budget)?;
    let Some((_, k)) = f.knots[0] else {
        return Err(ReadError::Unsupported);
    };
    let k = list(k)?;
    let range = target.range();
    // NumericNurbs normalizes its knot domain. Exact equality here earns
    // interpreting its reserve at the SAME original global parameter t.
    if !normalized {
        for (i, index) in [0, k.len() - 1].into_iter().enumerate() {
            if !exact_binary(
                spelling(source, target_id, &k[index], budget)?,
                range[i],
                budget,
            )? {
                return Err(ReadError::Unsupported);
            }
        }
    }
    // The normalized numeric transfer is earned on the entire NATURAL
    // domain. An extrapolated cached parameter has no such authority.
    if normalized && (edge.range[0] < range[0] || edge.range[1] > range[1]) {
        return Err(ReadError::Unsupported);
    }
    budget.charge(256)?;
    let raw_span = valid(natural[1] - natural[0])?;
    let cached_span = valid(I::point(range[1]) - I::point(range[0]))?;
    if raw_span.lo <= 0. || cached_span.lo <= 0. {
        return Err(ReadError::Unsupported);
    }
    let cached_raw = if normalized {
        edge.range
            .map(|t| natural[0] + (I::point(t) - I::point(range[0])) / cached_span * raw_span)
    } else {
        edge.range.map(I::point)
    };
    if cached_raw.iter().any(|v| !v.valid()) {
        return Err(ReadError::Unsupported);
    }
    let clamp = if normalized {
        [natural[0].hi, natural[1].lo]
    } else {
        range
    };
    if clamp.iter().any(|x| !x.is_finite()) || clamp[0] > clamp[1] {
        return Err(ReadError::Unsupported);
    }
    let numeric = curve_nurbs(source, target_id, &edge.curve, 3, budget)?.error;
    curve_nurbs(source, definition.curve_id, &definition.geometry, 2, budget)?;
    let domain = knot_domain(source, definition.curve_id, uv, false, budget)?;
    let width = valid(domain[1] - domain[0])?;
    if width.lo <= 0. {
        return Err(ReadError::Unsupported);
    }
    let (gamma, speed, acceleration) = bounds(data, budget)?;
    let mut closest = [0.; 2];
    let mut cached: f64 = 0.;
    for i in 0..2 {
        budget.charge(1024)?;
        // The legacy branch has exact binary natural endpoints. For the
        // normalized branch [RawA.hi,RawB.lo] is CERTAINLY inside every
        // actual original natural domain. Clamp a candidate into that
        // interior and include its complete distance to the UV endpoint;
        // no original endpoint is silently replaced by its cached value.
        let t = valid(I::new(
            domain[i].lo.clamp(clamp[0], clamp[1]),
            domain[i].hi.clamp(clamp[0], clamp[1]),
        ))?;
        let offset = magnitude(domain[i] - t)?;
        let point = evaluate(data, t, budget)?;
        let vertex = raw_vertex(source, edge.vertices[i], budget)?;
        let r = valid(norm(std::array::from_fn(|j| point[j] - vertex[j])))?;
        closest[i] = valid(I::point(offset) + I::point(2.) * I::point(r.hi) / I::point(gamma))?.hi;
        cached = cached.max(magnitude(domain[i] - cached_raw[i])?);
    }
    if valid(width - I::point(closest[0]) - I::point(closest[1]))?.lo <= 0. {
        return Err(ReadError::Unsupported);
    }
    // For every global distance minimizer t*: d(C(t*),V)<=r. The strict
    // monotone coordinate gives gamma*|t*-t|<=|C(t*)-C(t)|<=2r.
    // Thus ALL closest choices are enclosed, including endpoint minima.
    budget.charge(64)?;
    Ok(Window {
        source,
        edge,
        definition,
        target_id,
        uv_id: definition.curve_id,
        uv_domain: domain,
        gamma,
        speed,
        acceleration,
        closest,
        cached,
        numeric,
    })
}
#[cfg(test)]
mod tests;

#[cfg(test)]
mod normalized_tests;
