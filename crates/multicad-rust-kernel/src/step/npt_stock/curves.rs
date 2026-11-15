use super::*;
pub(super) type Basis = [[I; 12]; 12];
fn p(x: f64) -> I {
    I::point(x)
}
fn power(x: I, n: usize) -> I {
    (0..n).fold(p(1.), |v, _| v * x)
}
fn choose(n: usize, k: usize) -> f64 {
    (0..k.min(n - k)).fold(1., |v, i| v * (n - i) as f64 / (i + 1) as f64)
}
fn norm(v: [I; 3]) -> I {
    v.into_iter().fold(p(0.), |a, x| a + p(x.magnitude()))
}
fn literal(x: f64) -> I {
    I::new(x.next_down(), x.next_up())
}
pub(super) fn basis(work: &mut Work<'_>) -> Result<Basis, Error> {
    let mut basis = [[p(0.); 12]; 12];
    for (i, row) in basis.iter_mut().enumerate() {
        for (k, x) in row.iter_mut().enumerate() {
            for j in 0..=k.min(i) {
                work.charge(16)?;
                *x = *x
                    + p(choose(k, j)) * power(p(-0.5), k - j) * p(choose(i, j)) / p(choose(11, j));
            }
        }
    }
    Ok(basis)
}
pub(super) fn owner(curve: StockClipCurve, index: usize, count: usize) -> CurvePointOwner {
    if index == 0 {
        return CurvePointOwner::Vertex(curve.endpoints()[0]);
    }
    if index == count {
        return CurvePointOwner::Vertex(curve.endpoints()[1]);
    }
    let mut numerator = index;
    let mut denominator = count;
    while numerator.is_multiple_of(2) {
        numerator /= 2;
        denominator /= 2;
    }
    CurvePointOwner::Interior {
        curve,
        numerator: numerator as u16,
        denominator: denominator as u16,
    }
}
fn point(
    law: &StockClipCurveLaw<'_>,
    vertices: &[CarrierVertex; 12],
    index: usize,
    count: usize,
    work: &mut Work<'_>,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<[f64; 3], Error> {
    if let CurvePointOwner::Vertex(key) = owner(law.curve(), index, count) {
        // The sealed source graph contains every endpoint owner exactly once.
        work.charge(12)?;
        return Ok(vertices
            .iter()
            .find(|v| v.owner == key)
            .expect("original source endpoint")
            .position);
    }
    let t = index as f64 / count as f64;
    let result = law
        .plan()
        .curve_bounds(law.curve(), [t, t], source_options(work)?, interrupt)
        .map_err(source_error)?;
    work.charge(result.work_used())?;
    Ok(result.position_mm().map(|b| b[0].midpoint(b[1])))
}
fn span(
    law: &StockClipCurveLaw<'_>,
    domain: [f64; 2],
    ends: [[f64; 3]; 2],
    basis: &Basis,
    work: &mut Work<'_>,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<BezierSpan, Error> {
    let h = domain[1] - domain[0]; // Exact positive dyadic subdivision width.
    let center = domain[0].midpoint(domain[1]);
    let local = law
        .taylor([center, center], source_options(work)?, interrupt)
        .map_err(source_error)?;
    work.charge(local.work_used())?;
    let whole = law
        .taylor(domain, source_options(work)?, interrupt)
        .map_err(source_error)?;
    work.charge(whole.work_used())?;
    let coefficients = local
        .coefficients_mm()
        .map(|v| v.map(|b| I::new(b[0], b[1])));
    let mut ideal = [[p(0.); 3]; 12];
    for (i, pole) in ideal.iter_mut().enumerate() {
        for (k, c) in coefficients[..12].iter().enumerate() {
            work.charge(32)?;
            for j in 0..3 {
                pole[j] = pole[j] + c[j] * power(p(h), k) * basis[i][k];
            }
        }
    }
    let mut poles = ideal.map(|v| v.map(I::midpoint));
    poles[0] = ends[0];
    poles[11] = ends[1];
    if matches!(law.curve(), StockClipCurve::Ring { .. }) {
        // The source station token owns ONE axial coordinate for this whole
        // ring, including the symbolic contact. Choose one literal plane;
        // include the actual resulting transfer in BOTH certificates below.
        // All endpoint owners at that station use the same height expression.
        for pole in &mut poles {
            pole[2] = ends[0][2];
        }
    }
    let errors: [[I; 3]; 12] =
        std::array::from_fn(|i| std::array::from_fn(|j| literal(poles[i][j]) - ideal[i][j]));
    let coefficient_error = errors.into_iter().map(norm).fold(0f64, |a, b| a.max(b.hi));
    let derivative_error = (0..11)
        .map(|i| norm(std::array::from_fn(|j| errors[i + 1][j] - errors[i][j])))
        .fold(0f64, |a, b| a.max(b.hi));
    let remainder = norm(whole.coefficients_mm()[12].map(|b| I::new(b[0], b[1])));
    let position_error = (p(coefficient_error) + remainder * power(p(h / 2.), 12)).hi;
    let derivative_error =
        (p(11.) * p(derivative_error) / p(h) + p(12.) * remainder * power(p(h / 2.), 11)).hi;
    work.charge(1024)?;
    // This helper receives ONLY positive dyadic widths h>=1/256 from fit().
    // Source Taylor coefficients through12 are bounded by1e50 (sealed source
    // angular rate<900, radius/height<=1e7). The centered Bernstein basis
    // magnitude is <(3/2)^11<87, with outward slack<88; h^k<=1. Thus each
    // ideal/literal coefficient and difference is <1e55. The derivative
    // multiplier11/h<=2816 and remainder arithmetic stay below1e60. No
    // computed overflow/NaN case exists on this admitted numeric domain;
    // requested precision is still checked against these actual bounds.
    // Owners are assigned by the exact uniform source station in fit().
    Ok(BezierSpan {
        domain,
        owners: [CurvePointOwner::Vertex(law.curve().endpoints()[0]); 2],
        poles,
        position_error,
        derivative_error,
    })
}
pub(super) fn fit(
    law: &StockClipCurveLaw<'_>,
    vertices: &[CarrierVertex; 12],
    basis: &Basis,
    options: CurveOptions,
    memory: &mut Memory,
    work: &mut Work<'_>,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Vec<BezierSpan>, Error> {
    let mut count = 1usize;
    loop {
        let mut spans = memory.allocate(count, work)?;
        let mut previous = point(law, vertices, 0, count, work, interrupt)?;
        let mut accepted = true;
        for index in 0..count {
            let next = point(law, vertices, index + 1, count, work, interrupt)?;
            let mut candidate = span(
                law,
                [
                    index as f64 / count as f64,
                    (index + 1) as f64 / count as f64,
                ],
                [previous, next],
                basis,
                work,
                interrupt,
            )?;
            if candidate.position_error > options.position_error_mm
                || candidate.derivative_error > options.derivative_error_mm
            {
                accepted = false;
                break;
            }
            candidate.owners = [
                owner(law.curve(), index, count),
                owner(law.curve(), index + 1, count),
            ];
            spans.push(candidate);
            previous = next;
        }
        if accepted {
            return Ok(spans);
        }
        memory.release(&spans);
        drop(spans);
        if count == options.maximum_spans_per_curve {
            return Err(Error::Unresolved);
        }
        count *= 2;
    }
}
