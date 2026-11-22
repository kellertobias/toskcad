//! Source-bound tensor extraction. Original rational knot ratios and directed
//! homogeneous REAL coefficients enclose Boehm's exact affine identities.
use super::*;
use crate::brep::NurbsSurface;
use crate::step_read::brep::check;
type H = [I; 4];
const MAX_BYTES: usize = 8 * 1024 * 1024;
const SCRATCH: usize = 128 * 1024;
const MAX_CONTROLS: usize = 65_536;
const MAX_KNOTS: usize = 1024;
/// Homogeneous tensor Bernstein controls on local `[0,1]^2`. The actual
/// parameters are each original knot interval's positive affine image;
/// knot ordinals identify those same immutable source records.
pub(super) struct Patch {
    knots: [[usize; 2]; 2],
    shape: [usize; 2],
    controls: Vec<H>,
}
impl Patch {
    pub(super) fn knot_indices(&self) -> [[usize; 2]; 2] {
        self.knots
    }
    pub(super) fn shape(&self) -> [usize; 2] {
        self.shape
    }
    pub(super) fn controls(&self) -> &[H] {
        &self.controls
    }
}
pub(super) struct SurfaceCells<'s, 'i> {
    source: &'s Source<'i>,
    id: u32,
    degrees: [usize; 2],
    counts: [usize; 2],
    patches: Vec<Patch>,
    peak: usize,
    retained: usize,
}
impl SurfaceCells<'_, '_> {
    pub(super) fn matches(&self, source: &Source<'_>, id: u32) -> bool {
        std::ptr::eq(self.source, source) && self.id == id
    }
    pub(super) fn degrees(&self) -> [usize; 2] {
        self.degrees
    }
    pub(super) fn counts(&self) -> [usize; 2] {
        self.counts
    }
    pub(super) fn patches(&self) -> &[Patch] {
        &self.patches
    }
    pub(super) fn peak_bytes(&self) -> usize {
        self.peak
    }
    pub(super) fn retained_bytes(&self) -> usize {
        self.retained
    }
}
struct Storage {
    current: usize,
    peak: usize,
}
impl Storage {
    fn add(&mut self, n: usize) -> Result<()> {
        let next = self.current.checked_add(n).ok_or(ReadError::Limit)?;
        check(next <= MAX_BYTES, ReadError::Limit)?;
        self.current = next;
        self.peak = self.peak.max(next);
        Ok(())
    }
    fn remove(&mut self, n: usize) {
        self.current -= n;
    }
}
fn reserve<T>(n: usize, s: &mut Storage, b: &mut Budget<'_>) -> Result<Vec<T>> {
    b.charge(64)?;
    let bytes = n
        .checked_mul(std::mem::size_of::<T>())
        .ok_or(ReadError::Limit)?;
    check(
        s.current.checked_add(bytes).is_some_and(|v| v <= MAX_BYTES),
        ReadError::Limit,
    )?;
    let mut out = Vec::new();
    out.try_reserve_exact(n).map_err(|_| ReadError::Limit)?;
    s.add(
        out.capacity()
            .checked_mul(std::mem::size_of::<T>())
            .ok_or(ReadError::Limit)?,
    )?;
    Ok(out)
}
fn bytes<T>(v: &Vec<T>) -> usize {
    v.capacity() * std::mem::size_of::<T>()
}
fn bits(q: &BigRational) -> Result<usize> {
    let n = q.numer().bits().max(q.denom().bits());
    check(n <= 4096, ReadError::Limit)?;
    usize::try_from(n).map_err(|_| ReadError::Limit)
}
fn qop(a: &BigRational, c: &BigRational, divide: bool, b: &mut Budget<'_>) -> Result<BigRational> {
    let n = bits(a)?.saturating_add(bits(c)?).saturating_add(2);
    check(n <= 4096, ReadError::Limit)?;
    b.charge(bigint_work(n))?;
    check(!divide || !c.is_zero(), ReadError::Unsupported)?;
    let q = if divide { a / c } else { a - c };
    bits(&q)?;
    Ok(q)
}
struct Axis {
    original: Vec<BigRational>,
    multiplicities: Vec<usize>,
    tags: Vec<usize>,
}
fn spelling_scratch(source: &Source<'_>, id: u32, b: &mut Budget<'_>) -> Result<()> {
    fn visit(v: &Value<'_>, count: &mut usize, b: &mut Budget<'_>) -> Result<()> {
        b.charge(16)?;
        match v {
            Value::Number(_) => {
                *count = count.checked_add(1).ok_or(ReadError::Limit)?;
            }
            Value::List(values) => {
                for value in values {
                    visit(value, count, b)?;
                }
            }
            Value::Typed(entity) => {
                for value in &entity.args {
                    visit(value, count, b)?;
                }
            }
            _ => {}
        }
        Ok(())
    }
    for entity in source.components(id)? {
        b.charge(16)?;
        let mut count = 0;
        for value in &entity.args {
            visit(value, &mut count, b)?;
        }
        // numeric_spellings may geometrically allocate two capacities of
        // borrowed strings. Keep that real delegated temporary inside the
        // scratch allowance, leaving half for bounded rational temporaries.
        check(
            count
                .checked_mul(2 * std::mem::size_of::<&str>())
                .is_some_and(|v| v <= SCRATCH / 2),
            ReadError::Limit,
        )?;
    }
    Ok(())
}
fn axis(
    source: &Source<'_>,
    id: u32,
    field: Option<(&Value<'_>, &Value<'_>)>,
    (degree, count): (usize, usize),
    decoded: (&[f64], &[usize]),
    s: &mut Storage,
    b: &mut Budget<'_>,
) -> Result<Axis> {
    let n = field
        .map(|(_, k)| list(k).map(|v| v.len()))
        .transpose()?
        .unwrap_or(2);
    check((2..=MAX_KNOTS).contains(&n), ReadError::Limit)?;
    let mut out = Axis {
        original: reserve(n, s, b)?,
        multiplicities: reserve(n, s, b)?,
        tags: reserve(count.checked_add(degree + 1).ok_or(ReadError::Limit)?, s, b)?,
    };
    if let Some((m, k)) = field {
        let (m, k) = (list(m)?, list(k)?);
        check(
            m.len() == n && decoded.1.len() == n && decoded.0.len() == n,
            ReadError::Unsupported,
        )?;
        for i in 0..n {
            b.charge(128)?;
            integer(source, id, &m[i], decoded.1[i], b)?;
            let multiplicity = decoded.1[i];
            check(
                multiplicity > 0
                    && multiplicity <= degree + 1
                    && (if i == 0 || i + 1 == n {
                        multiplicity == degree + 1
                    } else {
                        multiplicity <= degree
                    }),
                ReadError::Unsupported,
            )?;
            let q = rational(source, id, &k[i], b)?;
            let size = bits(&q)?;
            // The shared decimal reader caps original inputs at 1023 decimal
            // places/4096 bits. Charge 4096 bytes per retained rational for
            // both BigInt payloads, including normalization spare capacity;
            // reduced bit length alone would undercount retained allocation.
            s.add(4096)?;
            let value = scalar(source, id, &k[i], b)?;
            check(value.valid(), ReadError::Unsupported)?;
            if let Some(previous) = out.original.last() {
                b.charge(bigint_work(
                    bits(previous)?.saturating_add(size).saturating_add(2),
                ))?;
                check(previous < &q, ReadError::Unsupported)?;
            }
            out.original.push(q);
            out.multiplicities.push(multiplicity);
            check(
                out.tags
                    .len()
                    .checked_add(multiplicity)
                    .is_some_and(|v| v <= count + degree + 1),
                ReadError::Unsupported,
            )?;
            b.charge(multiplicity * 4)?;
            out.tags.extend(std::iter::repeat_n(i, multiplicity));
        }
    } else {
        check(count == degree + 1, ReadError::Unsupported)?;
        out.original.extend([
            BigRational::from_integer(0.into()),
            BigRational::from_integer(1.into()),
        ]);
        out.multiplicities.extend([degree + 1; 2]);
        out.tags.extend(std::iter::repeat_n(0, degree + 1));
        out.tags.extend(std::iter::repeat_n(1, degree + 1));
        s.add(64)?;
    }
    check(out.tags.len() == count + degree + 1, ReadError::Unsupported)?;
    Ok(out)
}
fn ratio(axis: &Axis, t: usize, a: usize, c: usize, b: &mut Budget<'_>) -> Result<I> {
    if t == a {
        return Ok(I::point(0.));
    }
    if t == c {
        return Ok(I::point(1.));
    }
    check(a < t && t < c, ReadError::Unsupported)?;
    let numerator = qop(&axis.original[t], &axis.original[a], false, b)?;
    let denominator = qop(&axis.original[c], &axis.original[a], false, b)?;
    let q = qop(&numerator, &denominator, true, b)?;
    let value = enclose(&q, b)?;
    check(value.valid(), ReadError::Unsupported)?;
    let value = I::new(value.lo.max(0.), value.hi.min(1.));
    check(value.valid(), ReadError::Unsupported)?;
    Ok(value)
}
struct Net {
    counts: [usize; 2],
    values: Vec<H>,
    unit_weights: bool,
}
fn insert(
    net: &mut Net,
    axis_index: usize,
    axis: &mut Axis,
    t: usize,
    degree: usize,
    s: &mut Storage,
    b: &mut Budget<'_>,
) -> Result<()> {
    b.charge(axis.tags.len().saturating_mul(4).saturating_add(256))?;
    let k = axis
        .tags
        .iter()
        .rposition(|&v| v == t)
        .ok_or(ReadError::Unsupported)?;
    let multiplicity = axis.tags.iter().filter(|&&v| v == t).count();
    check(k >= degree && multiplicity < degree, ReadError::Unsupported)?;
    let mut counts = net.counts;
    counts[axis_index] = counts[axis_index].checked_add(1).ok_or(ReadError::Limit)?;
    let total = counts[0].checked_mul(counts[1]).ok_or(ReadError::Limit)?;
    check(total <= MAX_CONTROLS, ReadError::Limit)?;
    let mut alpha = reserve::<I>(counts[axis_index], s, b)?;
    for i in 0..counts[axis_index] {
        b.charge(8)?;
        alpha.push(if i > k - degree && i <= k - multiplicity {
            ratio(axis, t, axis.tags[i], axis.tags[i + degree], b)?
        } else {
            I::point(0.)
        });
    }
    let mut out = reserve::<H>(total, s, b)?;
    for u in 0..counts[0] {
        for v in 0..counts[1] {
            b.charge(128)?;
            let (i, j) = if axis_index == 0 { (u, v) } else { (v, u) };
            let get = |i| {
                if axis_index == 0 {
                    net.values[i * net.counts[1] + j]
                } else {
                    net.values[j * net.counts[1] + i]
                }
            };
            let h = if i <= k - degree {
                get(i)
            } else if i > k - multiplicity {
                get(i - 1)
            } else {
                let (a, c) = (get(i - 1), get(i));
                let factor = alpha[i];
                let mut h = std::array::from_fn(|q| a[q] * (I::point(1.) - factor) + c[q] * factor);
                if net.unit_weights {
                    h[3] = I::point(1.)
                }
                h
            };
            check(
                h.iter().all(|v| v.valid()) && h[3].lo > 0.,
                ReadError::Unsupported,
            )?;
            out.push(h);
        }
    }
    let mut tags = reserve(axis.tags.len() + 1, s, b)?;
    b.charge(axis.tags.len().saturating_mul(4).saturating_add(64))?;
    tags.extend_from_slice(&axis.tags[..=k]);
    tags.push(t);
    tags.extend_from_slice(&axis.tags[k + 1..]);
    s.remove(bytes(&alpha));
    drop(alpha);
    s.remove(bytes(&net.values));
    net.values = out;
    s.remove(bytes(&axis.tags));
    axis.tags = tags;
    net.counts = counts;
    Ok(())
}
pub(super) fn original_surface<'s, 'i>(
    source: &'s Source<'i>,
    id: u32,
    decoded: &NurbsSurface,
    b: &mut Budget<'_>,
) -> Result<SurfaceCells<'s, 'i>> {
    b.charge(4096)?;
    let degrees = decoded.degrees();
    let counts = decoded.counts();
    check(
        degrees.iter().all(|&p| (1..=11).contains(&p))
            && counts.iter().zip(degrees).all(|(&n, p)| n > p),
        ReadError::Unsupported,
    )?;
    let total = counts[0].checked_mul(counts[1]).ok_or(ReadError::Limit)?;
    check(total <= MAX_CONTROLS, ReadError::Limit)?;
    spelling_scratch(source, id, b)?;
    let field = fields(source, id, 2, b)?;
    for (i, p) in degrees.iter().enumerate() {
        integer(source, id, field.degrees[i], *p, b)?;
    }
    let mut storage = Storage {
        current: SCRATCH,
        peak: SCRATCH,
    };
    let mut axes = [
        axis(
            source,
            id,
            field.knots[0],
            (degrees[0], counts[0]),
            decoded.knots(0),
            &mut storage,
            b,
        )?,
        axis(
            source,
            id,
            field.knots[1],
            (degrees[1], counts[1]),
            decoded.knots(1),
            &mut storage,
            b,
        )?,
    ];
    let rows = list(field.controls)?;
    check(rows.len() == counts[0], ReadError::Unsupported)?;
    let weight_rows = field.weights.map(list).transpose()?;
    check(
        weight_rows.is_none_or(|r| r.len() == counts[0]),
        ReadError::Unsupported,
    )?;
    let mut net = Net {
        counts,
        values: reserve(total, &mut storage, b)?,
        unit_weights: true,
    };
    for (u, row) in rows.iter().enumerate() {
        let row = list(row)?;
        check(row.len() == counts[1], ReadError::Unsupported)?;
        let weights = weight_rows.map(|r| list(&r[u])).transpose()?;
        check(
            weights.is_none_or(|r| r.len() == counts[1]),
            ReadError::Unsupported,
        )?;
        for (v, p) in row.iter().enumerate() {
            b.charge(128)?;
            let p = point(source, reference(p)?, 3, b)?;
            let w = weights
                .map(|w| scalar(source, id, &w[v], b))
                .transpose()?
                .unwrap_or(I::point(1.));
            check(w.valid() && w.lo > 0., ReadError::Unsupported)?;
            net.unit_weights &= w.lo == 1. && w.hi == 1.;
            let h = [p[0] * w, p[1] * w, p[2] * w, w];
            check(h.iter().all(|v| v.valid()), ReadError::Unsupported)?;
            net.values.push(h);
        }
    }
    for a in 0..2 {
        for t in 1..axes[a].original.len() - 1 {
            let m = axes[a].multiplicities[t];
            for _ in m..degrees[a] {
                insert(&mut net, a, &mut axes[a], t, degrees[a], &mut storage, b)?;
            }
        }
    }
    let patch_counts = axes.each_ref().map(|a| a.original.len() - 1);
    let patch_total = patch_counts[0]
        .checked_mul(patch_counts[1])
        .ok_or(ReadError::Limit)?;
    let mut patches = reserve::<Patch>(patch_total, &mut storage, b)?;
    let mut output_bytes = bytes(&patches);
    for u in degrees[0]..net.counts[0] {
        if axes[0].tags[u] == axes[0].tags[u + 1] {
            continue;
        }
        for v in degrees[1]..net.counts[1] {
            if axes[1].tags[v] == axes[1].tags[v + 1] {
                continue;
            }
            b.charge(256)?;
            let knots = [
                [axes[0].tags[u], axes[0].tags[u + 1]],
                [axes[1].tags[v], axes[1].tags[v + 1]],
            ];
            let shape = degrees.map(|p| p + 1);
            let mut controls = reserve::<H>(shape[0] * shape[1], &mut storage, b)?;
            for i in u - degrees[0]..=u {
                for j in v - degrees[1]..=v {
                    b.charge(8)?;
                    controls.push(net.values[i * net.counts[1] + j]);
                }
            }
            output_bytes = output_bytes
                .checked_add(bytes(&controls))
                .ok_or(ReadError::Limit)?;
            patches.push(Patch {
                knots,
                shape,
                controls,
            });
        }
    }
    check(patches.len() == patch_total, ReadError::Unsupported)?;
    b.charge(64)?;
    Ok(SurfaceCells {
        source,
        id,
        degrees,
        counts: patch_counts,
        patches,
        peak: storage.peak,
        retained: output_bytes,
    })
}

#[cfg(test)]
mod tests;
