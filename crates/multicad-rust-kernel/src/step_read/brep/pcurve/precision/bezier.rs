//! Bounded extraction of ORIGINAL-REAL open rational curves. These immutable
//! coefficient enclosures are not geometric admission or a decoded proxy.
use super::*;
use crate::brep::NurbsCurve;
type H = [I; 4];
const MAX_DEGREE: usize = 11;
const MAX_CONTROLS: usize = 1024;
const MAX_SPANS: usize = 128;
const SCRATCH: usize = 128 * 1024;
const MAX_BYTES: usize = 8 * 1024 * 1024;
#[derive(Clone, Copy)]
pub(super) enum Coordinates {
    Raw(usize),
}
pub(super) struct Cell {
    domain: [I; 2],
    knots: [usize; 2],
    controls: Vec<H>,
}
impl Cell {
    pub(super) fn domain(&self) -> [I; 2] {
        self.domain
    }
    pub(super) fn knot_indices(&self) -> [usize; 2] {
        self.knots
    }
    pub(super) fn controls(&self) -> &[H] {
        &self.controls
    }
}
pub(super) struct CurveCells<'s, 'i> {
    source: &'s Source<'i>,
    source_id: u32,
    coordinates: Coordinates,
    degree: usize,
    cells: Vec<Cell>,
    peak_bytes: usize,
    retained_bytes: usize,
}
impl CurveCells<'_, '_> {
    /// Source/record association ONLY. A consumer must also inspect the
    /// retained coordinates() chart before interpreting any coefficient.
    pub(super) fn matches(&self, source: &Source<'_>, id: u32) -> bool {
        std::ptr::eq(self.source, source) && id == self.source_id
    }
    pub(super) fn coordinates(&self) -> Coordinates {
        self.coordinates
    }
    pub(super) fn degree(&self) -> usize {
        self.degree
    }
    pub(super) fn cells(&self) -> &[Cell] {
        &self.cells
    }
    /// Fresh helper heap/scratch only. The caller retains the complete Source
    /// arena separately and must compose its own simultaneous resident peak.
    pub(super) fn peak_bytes(&self) -> usize {
        self.peak_bytes
    }
    pub(super) fn retained_bytes(&self) -> usize {
        self.retained_bytes
    }
}
struct Storage {
    current: usize,
    peak: usize,
}
impl Storage {
    fn add(&mut self, bytes: usize) -> Result<()> {
        let next = self.current.checked_add(bytes).ok_or(ReadError::Limit)?;
        if next > MAX_BYTES {
            return Err(ReadError::Limit);
        }
        self.current = next;
        self.peak = self.peak.max(next);
        Ok(())
    }
    fn preflight(&self, bytes: usize) -> Result<()> {
        if self
            .current
            .checked_add(bytes)
            .is_none_or(|v| v > MAX_BYTES)
        {
            Err(ReadError::Limit)
        } else {
            Ok(())
        }
    }
    fn remove(&mut self, bytes: usize) {
        self.current -= bytes;
    }
}
fn reserve<T>(count: usize, storage: &mut Storage, budget: &mut Budget<'_>) -> Result<Vec<T>> {
    budget.charge(64)?;
    storage.preflight(
        count
            .checked_mul(std::mem::size_of::<T>())
            .ok_or(ReadError::Limit)?,
    )?;
    let mut v = Vec::new();
    v.try_reserve_exact(count).map_err(|_| ReadError::Limit)?;
    storage.add(
        v.capacity()
            .checked_mul(std::mem::size_of::<T>())
            .ok_or(ReadError::Limit)?,
    )?;
    Ok(v)
}
fn coordinate(
    source: &Source<'_>,
    id: u32,
    chart: Coordinates,
    budget: &mut Budget<'_>,
) -> Result<V> {
    let Coordinates::Raw(dim) = chart;
    if dim != 2 && dim != 3 {
        return Err(ReadError::Unsupported);
    }
    point(source, id, dim, budget)
}
fn alpha(t: usize, a: usize, b: usize, knots: &[I]) -> Result<I> {
    if t == a {
        return Ok(I::point(0.));
    }
    if t == b {
        return Ok(I::point(1.));
    }
    if !(a < t && t < b) {
        return Err(ReadError::Unsupported);
    }
    let denominator = knots[b] - knots[a];
    if !denominator.valid() || denominator.lo <= 0. {
        return Err(ReadError::Unsupported);
    }
    let ratio = (knots[t] - knots[a]) / denominator;
    if !ratio.valid() {
        return Err(ReadError::Unsupported);
    }
    // Strict original knot order proves exact alpha lies in [0,1]. This
    // intersection narrows an enclosure, never repairs a floating proposal.
    let ratio = I::new(ratio.lo.max(0.), ratio.hi.min(1.));
    if !ratio.valid() {
        Err(ReadError::Unsupported)
    } else {
        Ok(ratio)
    }
}
#[allow(
    clippy::too_many_arguments,
    reason = "Knot identities, literal coefficients, exact constant-weight proof and caller ledgers are distinct proof inputs."
)]
fn insert(
    degree: usize,
    t: usize,
    knots: &[I],
    unit_weights: bool,
    net: &mut Vec<H>,
    tags: &mut Vec<usize>,
    storage: &mut Storage,
    budget: &mut Budget<'_>,
) -> Result<()> {
    budget.charge(tags.len().saturating_mul(4).saturating_add(64))?;
    let k = tags
        .iter()
        .rposition(|&tag| tag == t)
        .ok_or(ReadError::Invalid)?;
    let count = tags.iter().filter(|&&tag| tag == t).count();
    if count >= degree || k < degree || k - count >= net.len() {
        return Err(ReadError::Unsupported);
    }
    let mut new = reserve::<H>(net.len() + 1, storage, budget)?;
    // The precharged initialization writes all new coefficient slots before
    // any indexed assignment. Old and replacement capacities overlap here.
    budget.charge((net.len() + 1).saturating_mul(4))?;
    new.resize(net.len() + 1, [I::point(0.); 4]);
    for i in 0..=k - degree {
        budget.charge(4)?;
        new[i] = net[i];
    }
    for i in k - count..net.len() {
        budget.charge(4)?;
        new[i + 1] = net[i];
    }
    for i in k - degree + 1..=k - count {
        budget.charge(256)?;
        let a = alpha(t, tags[i], tags[i + degree], knots)?;
        new[i] = std::array::from_fn(|j| a * net[i][j] + (I::point(1.) - a) * net[i - 1][j]);
        if unit_weights {
            // Earned from EVERY original raw weight being exactly one.
            // The Bernstein constant identity alpha*1+(1-alpha)*1=1
            // preserves this component, including through all insertions.
            // No cached or merely overlapping-one weight earns this branch.
            budget.charge(4)?;
            new[i][3] = I::point(1.);
        }
        if new[i].into_iter().any(|v| !v.valid()) || new[i][3].lo <= 0. {
            return Err(ReadError::Unsupported);
        }
    }
    let mut new_tags = reserve::<usize>(tags.len() + 1, storage, budget)?;
    for (i, &tag) in tags.iter().enumerate() {
        budget.charge(if i == k { 8 } else { 4 })?;
        new_tags.push(tag);
        if i == k {
            new_tags.push(t)
        }
    }
    storage.remove(net.capacity() * std::mem::size_of::<H>());
    *net = new;
    storage.remove(tags.capacity() * std::mem::size_of::<usize>());
    *tags = new_tags;
    Ok(())
}
/// Original knot ordinals survive every insertion and the resulting Bezier
/// cells keep their original scalar parameter domains. Every coordinate is
/// enclosed from the original raw record; no f64 geometry is substituted.
pub(super) fn original_curve<'s, 'i>(
    source: &'s Source<'i>,
    id: u32,
    decoded: &NurbsCurve,
    coordinates: Coordinates,
    budget: &mut Budget<'_>,
) -> Result<CurveCells<'s, 'i>> {
    budget.charge(4096)?;
    let degree = decoded.degree();
    let count = decoded.points().len();
    if !(1..=MAX_DEGREE).contains(&degree) || !(degree + 1..=MAX_CONTROLS).contains(&count) {
        return Err(ReadError::Unsupported);
    }
    let parts = source.components(id)?;
    let closed = if let [entity] = parts {
        entity.args.get(4)
    } else {
        source
            .component(id, "B_SPLINE_CURVE")?
            .and_then(|base| base.get(3))
    };
    if !matches!(closed, Some(Value::Enum("F"))) {
        return Err(ReadError::Unsupported);
    }
    let field = fields(source, id, 1, budget)?;
    integer(source, id, field.degrees[0], degree, budget)?;
    let controls = list(field.controls)?;
    if controls.len() != count {
        return Err(ReadError::Unsupported);
    }
    let mut storage = Storage {
        current: SCRATCH,
        peak: SCRATCH,
    };
    let mut knots = reserve::<I>(MAX_SPANS + 1, &mut storage, budget)?;
    let mut tags = reserve::<usize>(count + degree + 1, &mut storage, budget)?;
    let mut multiplicities = reserve::<usize>(MAX_SPANS + 1, &mut storage, budget)?;
    if let Some((m, k)) = field.knots[0] {
        let (m, k) = (list(m)?, list(k)?);
        let (_, expected) = decoded.knots();
        if !(2..=MAX_SPANS + 1).contains(&k.len())
            || m.len() != k.len()
            || expected.len() != k.len()
        {
            return Err(ReadError::Unsupported);
        }
        for i in 0..k.len() {
            budget.charge(128)?;
            integer(source, id, &m[i], expected[i], budget)?;
            let n = expected[i];
            if n == 0
                || n > degree + 1
                || ((i == 0 || i + 1 == k.len()) && n != degree + 1)
                || (i > 0 && i + 1 < k.len() && n > degree)
            {
                return Err(ReadError::Unsupported);
            }
            let knot = scalar(source, id, &k[i], budget)?;
            if !knot.valid()
                || knots
                    .last()
                    .is_some_and(|p| !(knot - *p).valid() || (knot - *p).lo <= 0.)
            {
                return Err(ReadError::Unsupported);
            }
            if tags
                .len()
                .checked_add(n)
                .is_none_or(|v| v > count + degree + 1)
            {
                return Err(ReadError::Unsupported);
            }
            knots.push(knot);
            multiplicities.push(n);
            budget.charge(n.saturating_mul(4))?;
            tags.extend(std::iter::repeat_n(i, n));
        }
    } else {
        if count != degree + 1 {
            return Err(ReadError::Unsupported);
        }
        knots.extend([I::point(0.), I::point(1.)]);
        multiplicities.extend([degree + 1, degree + 1]);
        tags.extend(std::iter::repeat_n(0, degree + 1));
        tags.extend(std::iter::repeat_n(1, degree + 1));
    }
    if tags.len() != count + degree + 1 {
        return Err(ReadError::Unsupported);
    }
    let mut net = reserve::<H>(count, &mut storage, budget)?;
    let mut unit_weights = true;
    let raw_weights = field.weights.map(list).transpose()?;
    if raw_weights.is_some_and(|w| w.len() != count) {
        return Err(ReadError::Unsupported);
    }
    for (i, control) in controls.iter().enumerate() {
        budget.charge(256)?;
        let p = coordinate(source, reference(control)?, coordinates, budget)?;
        let w = if let Some(values) = raw_weights {
            scalar(source, id, &values[i], budget)?
        } else {
            I::point(1.)
        };
        if !w.valid() || w.lo <= 0. {
            return Err(ReadError::Unsupported);
        }
        budget.charge(4)?;
        unit_weights &= w.lo == 1. && w.hi == 1.;
        let h = [p[0] * w, p[1] * w, p[2] * w, w];
        if h.into_iter().any(|v| !v.valid()) {
            return Err(ReadError::Unsupported);
        }
        net.push(h);
    }
    for (t, &m) in multiplicities
        .iter()
        .enumerate()
        .take(knots.len() - 1)
        .skip(1)
    {
        budget.charge(8)?;
        for _ in m..degree {
            insert(
                degree,
                t,
                &knots,
                unit_weights,
                &mut net,
                &mut tags,
                &mut storage,
                budget,
            )?;
        }
    }
    let mut cells = reserve::<Cell>(knots.len() - 1, &mut storage, budget)?;
    for k in degree..net.len() {
        budget.charge(16)?;
        if tags[k] == tags[k + 1] {
            continue;
        }
        let mut cell = reserve::<H>(degree + 1, &mut storage, budget)?;
        for h in &net[k - degree..=k] {
            budget.charge(4)?;
            cell.push(*h)
        }
        cells.push(Cell {
            domain: [knots[tags[k]], knots[tags[k + 1]]],
            knots: [tags[k], tags[k + 1]],
            controls: cell,
        });
    }
    if cells.len() != knots.len() - 1 {
        return Err(ReadError::Unsupported);
    }
    budget.charge(cells.len().saturating_mul(4).saturating_add(64))?;
    let retained = std::mem::size_of::<CurveCells<'_, '_>>()
        + cells.capacity() * std::mem::size_of::<Cell>()
        + cells
            .iter()
            .map(|c| c.controls.capacity() * std::mem::size_of::<H>())
            .sum::<usize>();
    Ok(CurveCells {
        source,
        source_id: id,
        coordinates,
        degree,
        cells,
        peak_bytes: storage.peak,
        retained_bytes: retained,
    })
}
#[cfg(test)]
mod tests;
