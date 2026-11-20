//! Original STEP REALs are enclosed independently of their binary64 caches.
use super::super::{ReadError, Source, Value, list, part21, reference};
use super::{Budget, Result, difference, norm, vector as cached};
use crate::brep::{Curve, Edge, Surface};
use crate::step::phillips_cutter::support::Frame as RawFrame;
use crate::tapered_thread::interval::I;
use multicad_precision::Rational;
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{ToPrimitive, Zero};
mod analytic;
mod analytic_polynomial;
mod bernstein80;
mod bezier;
mod bilinear;
mod c0_curved;
mod c0_trace;
mod cone_polynomial;
mod cone_separation;
mod conic;
mod cylinder_ellipse;
mod cylinder_polynomial;
mod cylinder_trace;
mod finite_boundary;
mod polynomial_axis;
mod polynomial_trace;
mod sense_line;
mod spline_window;
mod surface_bezier;
mod tensor_trace;
pub(super) fn tensor_trace_certificate(
    source: &Source<'_>,
    id: u32,
    definition: &super::uv::Definition,
    surface: &Surface,
    edge: &Edge,
    tolerance: f64,
    budget: &mut Budget<'_>,
) -> Result<Option<super::Proof>> {
    tensor_trace::certify(source, id, definition, surface, edge, tolerance, budget)
}
type V = [I; 3];
/// Original finite tensor-boundary proof, including the complete numeric
/// and finite-window transfer. Its result must not be reserved a second time.
pub(super) fn finite_boundary_certificate(
    source: &Source<'_>,
    curve_id: u32,
    definition: &super::uv::Definition,
    surface: &Surface,
    edge: &Edge,
    tolerance: f64,
    budget: &mut Budget<'_>,
) -> Result<Option<super::Proof>> {
    finite_boundary::certify(
        source, curve_id, definition, surface, edge, tolerance, budget,
    )
}
/// Polynomial traces on tensor surfaces with C0 interior joints.
pub(super) fn c0_trace_certificate(
    source: &Source<'_>,
    curve_id: u32,
    definition: &super::uv::Definition,
    surface: &Surface,
    edge: &Edge,
    tolerance: f64,
    budget: &mut Budget<'_>,
) -> Result<Option<super::Proof>> {
    c0_trace::certify(
        source, curve_id, definition, surface, edge, tolerance, budget,
    )
}
pub(super) fn cone_polynomial_certificate(
    source: &Source<'_>,
    curve_id: u32,
    definition: &super::uv::Definition,
    surface: &Surface,
    edge: &Edge,
    tolerance: f64,
    budget: &mut Budget<'_>,
) -> Result<Option<super::Proof>> {
    cone_polynomial::certify(
        source, curve_id, definition, surface, edge, tolerance, budget,
    )
}
/// This distinct certificate owns the Original-REAL bounded target trace,
/// including its parameter map and whole-curve tangent orientation.
pub(super) fn cylinder_trace_certificate(
    source: &Source<'_>,
    curve_id: u32,
    definition: &super::uv::Definition,
    surface: &Surface,
    edge: &Edge,
    tolerance: f64,
    budget: &mut Budget<'_>,
) -> Result<Option<super::Proof>> {
    cylinder_trace::certify(
        source, curve_id, definition, surface, edge, tolerance, budget,
    )
}
pub(super) fn cylinder_ellipse_certificate(
    source: &Source<'_>,
    curve_id: u32,
    definition: &super::uv::Definition,
    surface: &Surface,
    edge: &Edge,
    tolerance: f64,
    budget: &mut Budget<'_>,
) -> Result<Option<super::Proof>> {
    cylinder_ellipse::certify(
        source, curve_id, definition, surface, edge, tolerance, budget,
    )
}
pub(super) fn bilinear_certificate(
    source: &Source<'_>,
    curve_id: u32,
    definition: &super::uv::Definition,
    surface: &Surface,
    edge: &Edge,
    tolerance: f64,
    budget: &mut Budget<'_>,
) -> Result<Option<super::Proof>> {
    bilinear::certify(
        source, curve_id, definition, surface, edge, tolerance, budget,
    )
}
/// This stage can only prove disjoint original curve images; no positive
/// agreement or material permission is returned by a separating witness.
pub(super) fn cone_separation_certificate(
    source: &Source<'_>,
    curve_id: u32,
    definition: &super::uv::Definition,
    surface: &Surface,
    edge: &Edge,
    tolerance: f64,
    budget: &mut Budget<'_>,
) -> Result<Option<super::Proof>> {
    cone_separation::certify(
        source, curve_id, definition, surface, edge, tolerance, budget,
    )
}
pub(super) fn cylinder_polynomial_certificate(
    source: &Source<'_>,
    curve_id: u32,
    definition: &super::uv::Definition,
    surface: &Surface,
    edge: &Edge,
    tolerance: f64,
    budget: &mut Budget<'_>,
) -> Result<Option<super::Proof>> {
    cylinder_polynomial::certify(
        source, curve_id, definition, surface, edge, tolerance, budget,
    )
}
const MAX_NUMBERS: usize = 262_144;
const TOKEN_WORKSPACE: usize = 8 * 1024 * 1024;

const MAX_SURFACE_CACHE: usize = 4096;
const CACHE_PEAK_BYTES: usize = 256 * 1024;
/// Invocation-local successful WHOLE-surface numeric transfers. The immutable
/// owner and entity ID bind the same Source::surface(id) reference each time.
/// Every queried finite trim/domain/sense check remains outside this cache.
pub(super) struct Cache<'s, 'i> {
    source: &'s Source<'i>,
    surfaces: Vec<(u32, f64)>,
    peak: usize,
}
impl<'s, 'i> Cache<'s, 'i> {
    pub(super) fn new(source: &'s Source<'i>) -> Self {
        Self {
            source,
            surfaces: Vec::new(),
            peak: std::mem::size_of::<Self>(),
        }
    }
    fn owner(&self, source: &Source<'_>) -> Result<()> {
        if !std::ptr::eq(self.source, source) {
            return Err(ReadError::Invalid);
        }
        Ok(())
    }
    fn get(&self, source: &Source<'_>, id: u32, budget: &mut Budget<'_>) -> Result<Option<f64>> {
        self.owner(source)?;
        budget.charge(self.surfaces.len().saturating_mul(8).saturating_add(64))?;
        Ok(self
            .surfaces
            .iter()
            .find_map(|&(key, value)| (key == id).then_some(value)))
    }
    fn remember(
        &mut self,
        source: &Source<'_>,
        id: u32,
        error: f64,
        budget: &mut Budget<'_>,
    ) -> Result<()> {
        self.owner(source)?;
        if !error.is_finite() || error < 0. {
            return Err(ReadError::Unsupported);
        }
        if self.surfaces.len() == MAX_SURFACE_CACHE {
            return Ok(());
        }
        let old = self
            .surfaces
            .capacity()
            .saturating_mul(std::mem::size_of::<(u32, f64)>());
        let grows = self.surfaces.len() == self.surfaces.capacity();
        let proposed = (self.surfaces.len() + 1).saturating_mul(std::mem::size_of::<(u32, f64)>());
        let metadata = std::mem::size_of::<Self>() + 64;
        let pre = if grows {
            old.saturating_add(proposed)
        } else {
            old
        };
        if pre.saturating_add(metadata) > CACHE_PEAK_BYTES {
            return Err(ReadError::Limit);
        }
        // Lookup/copy/comparison and the possible old tuple-buffer relocation
        // are actual new work, not a repeated debit for an unexecuted proof.
        budget.charge(self.surfaces.len().saturating_mul(8).saturating_add(128))?;
        self.surfaces
            .try_reserve_exact(1)
            .map_err(|_| ReadError::Limit)?;
        let retained = self
            .surfaces
            .capacity()
            .saturating_mul(std::mem::size_of::<(u32, f64)>());
        let peak = retained
            .saturating_add(if grows { old } else { 0 })
            .saturating_add(metadata);
        if peak > CACHE_PEAK_BYTES {
            return Err(ReadError::Limit);
        }
        self.peak = self.peak.max(peak);
        self.surfaces.push((id, error));
        Ok(())
    }
    pub(super) fn retained_bytes(&self) -> usize {
        std::mem::size_of::<Self>() + self.surfaces.capacity() * std::mem::size_of::<(u32, f64)>()
    }
    pub(super) fn peak_bytes(&self) -> usize {
        self.peak
    }
}

fn locate(
    items: &[Value<'_>],
    needle: *const Value<'_>,
    count: &mut usize,
    found: &mut Option<usize>,
    budget: &mut Budget<'_>,
) -> Result<()> {
    for value in items {
        budget.charge(4)?;
        match value {
            Value::Number(_) => {
                if std::ptr::eq(value, needle) {
                    *found = Some(*count);
                }
                *count = count.checked_add(1).ok_or(ReadError::Limit)?;
                if *count > MAX_NUMBERS {
                    return Err(ReadError::Limit);
                }
            }
            Value::List(items) => locate(items, needle, count, found, budget)?,
            Value::Typed(entity) => locate(&entity.args, needle, count, found, budget)?,
            _ => {}
        }
    }
    Ok(())
}
/// Records with at most this many numbers keep their spellings cached.
const CACHED_NUMBERS: usize = 1 << 16;
/// How many records keep their spellings at once.
const CACHED_RECORDS: usize = 2;
fn spelling<'i>(
    source: &Source<'i>,
    id: u32,
    value: &Value<'_>,
    budget: &mut Budget<'_>,
) -> Result<&'i str> {
    spelling_in(source, id, None, value, budget)
}
/// [`spelling`] searching only the named component of a complex record. A
/// leaf taken from that component is found there exactly as before; its
/// siblings, which may hold thousands of values but none of its numbers,
/// are not rescanned for every coefficient.
fn spelling_in<'i>(
    source: &Source<'i>,
    id: u32,
    component: Option<&str>,
    value: &Value<'_>,
    budget: &mut Budget<'_>,
) -> Result<&'i str> {
    let Value::Number(cached) = value else {
        return Err(ReadError::Invalid);
    };
    let leaf = std::ptr::from_ref(value) as usize;
    for entity in source.components(id)? {
        budget.charge(16)?;
        if component.is_some_and(|name| entity.name != name) {
            continue;
        }
        let text = entity.source.ok_or(ReadError::Unsupported)?;
        let key = (text.as_ptr() as usize, text.len());
        // A record tokenized before answers from its token table; the
        // numeric leaf is matched by identity as before.
        let caching = source.spellings.borrow().is_some();
        let hit = source.spellings.borrow().as_ref().and_then(|cache| {
            cache
                .iter()
                .find(|spelled| spelled.text == key)
                .map(|spelled| spelled.leaves.get(&leaf).map(|&i| spelled.tokens[i]))
        });
        let range = match hit {
            Some(None) => continue,
            Some(Some(range)) => {
                budget.charge(32)?;
                range
            }
            None => {
                let (mut count, mut found) = (0, None);
                locate(&entity.args, value, &mut count, &mut found, budget)?;
                let Some(index) = found else { continue };
                budget.charge(text.len().saturating_add(count))?;
                // Small geometric records may hide Typed numeric children in
                // their ignored labels. Their schema arity does not bound the
                // token Vec. Reserve at most half of the 128KiB small-helper
                // scratch for tokens; large spline coefficient/knot owners
                // retain their separate bound.
                let token_workspace = if matches!(
                    entity.name,
                    "CARTESIAN_POINT"
                        | "DIRECTION"
                        | "VECTOR"
                        | "CIRCLE"
                        | "ELLIPSE"
                        | "CYLINDRICAL_SURFACE"
                        | "CONICAL_SURFACE"
                        | "SPHERICAL_SURFACE"
                        | "TOROIDAL_SURFACE"
                ) {
                    64 * 1024
                } else {
                    TOKEN_WORKSPACE
                };
                // numeric_spellings grows a Vec geometrically. Two capacities
                // of borrowed str fat pointers fit the helper scratch envelope.
                if count.saturating_mul(2 * std::mem::size_of::<&str>()) > token_workspace {
                    return Err(ReadError::Limit);
                }
                let tokens = part21::numeric_spellings(entity).ok_or(ReadError::Unsupported)?;
                if tokens.len() != count {
                    return Err(ReadError::Unsupported);
                }
                if tokens
                    .capacity()
                    .saturating_mul(std::mem::size_of::<&str>())
                    > token_workspace
                {
                    return Err(ReadError::Limit);
                }
                let offset = |token: &str| token.as_ptr() as usize - text.as_ptr() as usize;
                let range = tokens
                    .get(index)
                    .map(|token| (offset(token), token.len()))
                    .ok_or(ReadError::Unsupported)?;
                if caching && count <= CACHED_NUMBERS {
                    // The table costs one more pass over the leaves and
                    // tokens already visited.
                    budget.charge(count.saturating_mul(8).saturating_add(64))?;
                    let mut leaves = std::collections::HashMap::new();
                    leaves.try_reserve(count).map_err(|_| ReadError::Limit)?;
                    record(&entity.args, &mut leaves);
                    let mut cache = source.spellings.borrow_mut();
                    if let (true, Some(cache)) = (leaves.len() == count, cache.as_mut()) {
                        if cache.len() == CACHED_RECORDS {
                            cache.remove(0);
                        }
                        cache.push(crate::step_read::brep::Spelled {
                            text: key,
                            leaves,
                            tokens: tokens.iter().map(|t| (offset(t), t.len())).collect(),
                        });
                    }
                }
                range
            }
        };
        let token = text
            .get(range.0..range.0 + range.1)
            .ok_or(ReadError::Unsupported)?;
        budget.charge(token.len().saturating_add(64))?;
        if token.parse::<f64>().ok().map(f64::to_bits) != Some(cached.to_bits()) {
            return Err(ReadError::Unsupported);
        }
        return Ok(token);
    }
    // A numerically identical leaf from another arena is not this source leaf.
    Err(ReadError::Invalid)
}
/// The address of every numeric leaf, numbered in traversal order.
fn record(items: &[Value<'_>], leaves: &mut std::collections::HashMap<usize, usize>) {
    for value in items {
        match value {
            Value::Number(_) => {
                let index = leaves.len();
                leaves.insert(std::ptr::from_ref(value) as usize, index);
            }
            Value::List(items) => record(items, leaves),
            Value::Typed(entity) => record(&entity.args, leaves),
            _ => {}
        }
    }
}
fn exact_binary(text: &str, value: f64, budget: &mut Budget<'_>) -> Result<bool> {
    budget.charge(text.len().saturating_add(1024))?;
    let Ok(q) = Rational::parse_decimal(text) else {
        return Ok(false);
    };
    let (n, d) = q.ratio();
    let Some(binary) = BigRational::from_float(value) else {
        return Ok(false);
    };
    // Successful bounded Rational inputs lie in [10^-18,10^18], so this
    // comparison uses <200-bit components, never an unbounded decimal bigint.
    Ok(BigRational::new(BigInt::from(n), BigInt::from(d)) == binary)
}
pub(super) fn scalar(
    source: &Source<'_>,
    id: u32,
    value: &Value<'_>,
    budget: &mut Budget<'_>,
) -> Result<I> {
    scalar_in(source, id, None, value, budget)
}
/// [`scalar`] for a leaf of one named component of a complex record.
fn scalar_in(
    source: &Source<'_>,
    id: u32,
    component: Option<&str>,
    value: &Value<'_>,
    budget: &mut Budget<'_>,
) -> Result<I> {
    let Value::Number(cached) = value else {
        return Err(ReadError::Invalid);
    };
    if !cached.is_finite() {
        return Err(ReadError::Unsupported);
    }
    let token = spelling_in(source, id, component, value, budget)?;
    if *cached == 0. {
        // Only original decimal ZERO can establish a true isocurve. A
        // nonzero underflow such as 1E-1000 must never become a zero direction.
        let mantissa = token.split('E').next().ok_or(ReadError::Unsupported)?;
        if !mantissa
            .bytes()
            .filter(|b| !matches!(b, b'+' | b'-' | b'.'))
            .all(|b| b == b'0')
        {
            return Err(ReadError::Unsupported);
        }
        return Ok(I::point(0.));
    }
    // Ordinary coordinates need no bigint exactness test: adjacent floats
    // already enclose their original REAL. Exact bounded integers are cheap
    // and retain cardinal axes and large integer gauge cancellation.
    let exact_integer = if cached.fract() == 0. && cached.abs() <= 9_007_199_254_740_992. {
        budget.charge(token.len().saturating_add(256))?;
        Rational::parse_decimal(token).is_ok_and(|q| q.ratio() == (*cached as i64, 1))
    } else {
        false
    };
    let value = if exact_integer {
        I::point(*cached)
    } else {
        I::new(cached.next_down(), cached.next_up())
    };
    if !value.valid() {
        return Err(ReadError::Unsupported);
    }
    Ok(value)
}
fn exact_same(
    source: &Source<'_>,
    a: (u32, &Value<'_>),
    b: (u32, &Value<'_>),
    budget: &mut Budget<'_>,
) -> Result<bool> {
    let (left, right) = (
        spelling(source, a.0, a.1, budget)?,
        spelling(source, b.0, b.1, budget)?,
    );
    budget.charge(
        left.len()
            .saturating_add(right.len())
            .saturating_mul(3)
            .saturating_add(1024),
    )?;
    for (value, token) in [(a.1, left), (b.1, right)] {
        let Value::Number(cached) = value else {
            return Err(ReadError::Invalid);
        };
        if !cached.is_finite() {
            return Err(ReadError::Unsupported);
        }
        if *cached == 0.
            && !token
                .split('E')
                .next()
                .ok_or(ReadError::Unsupported)?
                .bytes()
                .filter(|b| !matches!(b, b'+' | b'-' | b'.'))
                .all(|b| b == b'0')
        {
            return Err(ReadError::Unsupported);
        }
    }
    // Identity of actual validated numeric tokens proves exact REAL equality,
    // including tiny normal values beyond the fixed i64 rational range.
    if left == right {
        return Ok(true);
    }
    let left = rational(source, a.0, a.1, budget)?;
    let right = rational(source, b.0, b.1, budget)?;
    let bits = left
        .numer()
        .bits()
        .max(left.denom().bits())
        .saturating_add(right.numer().bits().max(right.denom().bits()))
        .saturating_add(2);
    budget.charge(bigint_work(
        usize::try_from(bits).map_err(|_| ReadError::Limit)?,
    ))?;
    Ok(left == right)
}
fn coordinates(
    source: &Source<'_>,
    id: u32,
    kind: &str,
    dim: usize,
    budget: &mut Budget<'_>,
) -> Result<V> {
    budget.charge(64)?;
    let [_, values] = source.args(id, kind)? else {
        return Err(ReadError::Invalid);
    };
    let values = list(values)?;
    if !(2..=3).contains(&dim) || values.len() != dim {
        return Err(ReadError::Unsupported);
    }
    let mut result = [I::point(0.); 3];
    for i in 0..dim {
        result[i] = scalar(source, id, &values[i], budget)?;
    }
    Ok(result)
}
pub(super) fn point(
    source: &Source<'_>,
    id: u32,
    dim: usize,
    budget: &mut Budget<'_>,
) -> Result<V> {
    coordinates(source, id, "CARTESIAN_POINT", dim, budget)
}
fn unit(p: V, budget: &mut Budget<'_>) -> Result<V> {
    budget.charge(512)?;
    if let Some(axis) = cardinal(p) {
        return Ok(cached(axis));
    }
    let scale = p.into_iter().map(I::magnitude).fold(0_f64, f64::max);
    if !scale.is_finite() || !scale.is_normal() {
        return Err(ReadError::Unsupported);
    }
    let p = p.map(|x| x / I::point(scale));
    let n = norm(p);
    if !n.valid() || n.lo <= 0. {
        return Err(ReadError::Unsupported);
    }
    let result = p.map(|x| x / n);
    if !result.into_iter().all(I::valid) {
        return Err(ReadError::Unsupported);
    }
    Ok(result)
}
pub(super) fn direction(
    source: &Source<'_>,
    id: u32,
    dim: usize,
    budget: &mut Budget<'_>,
) -> Result<V> {
    unit(coordinates(source, id, "DIRECTION", dim, budget)?, budget)
}
pub(super) fn vector(
    source: &Source<'_>,
    id: u32,
    dim: usize,
    budget: &mut Budget<'_>,
) -> Result<V> {
    budget.charge(64)?;
    let [_, axis, magnitude] = source.args(id, "VECTOR")? else {
        return Err(ReadError::Invalid);
    };
    let m = scalar(source, id, magnitude, budget)?;
    if m.lo <= 0. {
        return Err(ReadError::Unsupported);
    }
    Ok(direction(source, reference(axis)?, dim, budget)?.map(|x| x * m))
}
pub(super) fn frame(
    source: &Source<'_>,
    id: u32,
    dim: usize,
    budget: &mut Budget<'_>,
) -> Result<RawFrame> {
    budget.charge(2048)?;
    let (origin, normal, hint) = match dim {
        3 => {
            let [_, origin, normal, hint] = source.args(id, "AXIS2_PLACEMENT_3D")? else {
                return Err(ReadError::Invalid);
            };
            let normal = match normal {
                Value::Omitted => cached([0., 0., 1.]),
                value => coordinates(source, reference(value)?, "DIRECTION", 3, budget)?,
            };
            // Unqualified 3D implicit-reference selection is not guessed.
            if matches!(hint, Value::Omitted) {
                return Err(ReadError::Unsupported);
            }
            (
                point(source, reference(origin)?, 3, budget)?,
                normal,
                coordinates(source, reference(hint)?, "DIRECTION", 3, budget)?,
            )
        }
        2 => {
            let [_, origin, hint] = source.args(id, "AXIS2_PLACEMENT_2D")? else {
                return Err(ReadError::Invalid);
            };
            let hint = match hint {
                Value::Omitted => cached([1., 0., 0.]),
                value => coordinates(source, reference(value)?, "DIRECTION", 2, budget)?,
            };
            (
                point(source, reference(origin)?, 2, budget)?,
                cached([0., 0., 1.]),
                hint,
            )
        }
        _ => return Err(ReadError::Unsupported),
    };
    if let (Some(n), Some(h)) = (cardinal(normal), cardinal(hint)) {
        let cross = |a: [f64; 3], b: [f64; 3]| {
            [
                a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0],
            ]
        };
        let y = cross(n, h);
        if y.iter().filter(|x| **x != 0.).count() == 1 {
            return Ok(RawFrame {
                origin,
                axes: [cached(cross(y, n)), cached(y), cached(n)],
            });
        }
    }
    RawFrame::new(origin, normal, hint).map_err(|_| ReadError::Unsupported)
}
pub(super) fn unwrapped(source: &Source<'_>, mut id: u32, budget: &mut Budget<'_>) -> Result<u32> {
    for _ in 0..=8 {
        budget.charge(64)?;
        let parts = source.components(id)?;
        if let [entity] = parts
            && matches!(
                entity.name,
                "SURFACE_CURVE" | "SEAM_CURVE" | "INTERSECTION_CURVE"
            )
        {
            let [_, curve, ..] = entity.args.as_slice() else {
                return Err(ReadError::Invalid);
            };
            id = reference(curve)?;
        } else {
            return Ok(id);
        }
    }
    Err(ReadError::Unsupported)
}

struct Fields<'a, 'i> {
    degrees: Vec<&'a Value<'i>>,
    controls: &'a Value<'i>,
    knots: Vec<Option<(&'a Value<'i>, &'a Value<'i>)>>,
    weights: Option<&'a Value<'i>>,
}
fn fields<'a, 'i>(
    source: &'a Source<'i>,
    id: u32,
    dimensions: usize,
    budget: &mut Budget<'_>,
) -> Result<Fields<'a, 'i>> {
    budget.charge(256)?;
    let parts = source.components(id)?;
    let (base, knots, weights) = match parts {
        [e] => match (dimensions, e.name, e.args.as_slice()) {
            (1, "B_SPLINE_CURVE_WITH_KNOTS", [_, base @ .., m, k, _]) if base.len() == 5 => {
                (base, vec![Some((m, k))], None)
            }
            (1, "BEZIER_CURVE", [_, base @ ..]) if base.len() == 5 => (base, vec![None], None),
            (2, "B_SPLINE_SURFACE_WITH_KNOTS", [_, base @ .., mu, mv, ku, kv, _])
                if base.len() == 7 =>
            {
                (base, vec![Some((mu, ku)), Some((mv, kv))], None)
            }
            (2, "BEZIER_SURFACE", [_, base @ ..]) if base.len() == 7 => {
                (base, vec![None, None], None)
            }
            _ => return Err(ReadError::Unsupported),
        },
        _ => {
            let name = if dimensions == 1 {
                "B_SPLINE_CURVE"
            } else {
                "B_SPLINE_SURFACE"
            };
            let base = source.component(id, name)?.ok_or(ReadError::Unsupported)?;
            let name = if dimensions == 1 {
                "B_SPLINE_CURVE_WITH_KNOTS"
            } else {
                "B_SPLINE_SURFACE_WITH_KNOTS"
            };
            let knots = match (dimensions, source.component(id, name)?) {
                (1, Some([m, k, _])) => vec![Some((m, k))],
                (2, Some([mu, mv, ku, kv, _])) => vec![Some((mu, ku)), Some((mv, kv))],
                (1, None) if source.component(id, "BEZIER_CURVE")?.is_some() => vec![None],
                (2, None) if source.component(id, "BEZIER_SURFACE")?.is_some() => vec![None, None],
                _ => return Err(ReadError::Unsupported),
            };
            let name = if dimensions == 1 {
                "RATIONAL_B_SPLINE_CURVE"
            } else {
                "RATIONAL_B_SPLINE_SURFACE"
            };
            let weights = match source.component(id, name)? {
                Some([w]) => Some(w),
                None => None,
                _ => return Err(ReadError::Unsupported),
            };
            (base, knots, weights)
        }
    };
    if !(1..=2).contains(&dimensions) || base.len() != dimensions * 2 + 3 {
        return Err(ReadError::Unsupported);
    }
    Ok(Fields {
        degrees: base[..dimensions].iter().collect(),
        controls: &base[dimensions],
        knots,
        weights,
    })
}
fn integer(
    source: &Source<'_>,
    id: u32,
    value: &Value<'_>,
    expected: usize,
    budget: &mut Budget<'_>,
) -> Result<()> {
    let text = spelling(source, id, value, budget)?;
    budget.charge(text.len().saturating_add(1024))?;
    let q = Rational::parse_decimal(text).map_err(|_| ReadError::Unsupported)?;
    if q.ratio() != (expected as i64, 1) {
        return Err(ReadError::Unsupported);
    }
    Ok(())
}
fn visit<'a, 'i>(
    value: &'a Value<'i>,
    dimensions: usize,
    apply: &mut impl FnMut(&'a Value<'i>) -> Result<()>,
) -> Result<()> {
    for value in list(value)? {
        if dimensions == 1 {
            apply(value)?;
        } else {
            for value in list(value)? {
                apply(value)?;
            }
        }
    }
    Ok(())
}
struct NumericNurbs {
    error: f64,
    extent: f64,
}
#[allow(
    clippy::too_many_arguments,
    reason = "Borrowed original tables and reference tables remain explicit in this numerical transfer proof."
)]
fn nurbs_error(
    source: &Source<'_>,
    id: u32,
    field: &Fields<'_, '_>,
    degrees: &[usize],
    knots: &[(&[f64], &[usize])],
    points: &[[f64; 3]],
    weights: Option<&[f64]>,
    dimension: usize,
    mapped: Option<&[V]>,
    budget: &mut Budget<'_>,
) -> Result<NumericNurbs> {
    budget.charge(1024)?;
    if points.is_empty()
        || points.len() > crate::brep::nurbs::MAX_CONTROL_POINTS
        || degrees.len() != field.degrees.len()
    {
        return Err(ReadError::Unsupported);
    }
    let (mut delta_p, mut radius, mut extent) = (0_f64, 0_f64, 0_f64);
    let anchor = cached(points[0]);
    let mut index = 0;
    visit(field.controls, degrees.len(), &mut |value| {
        budget.charge(256)?;
        let decoded = *points.get(index).ok_or(ReadError::Unsupported)?;
        let raw = if let Some(mapped) = mapped {
            *mapped.get(index).ok_or(ReadError::Unsupported)?
        } else {
            point(source, reference(value)?, dimension, budget)?
        };
        let (e, l, w) = (
            norm(difference(raw, cached(decoded))),
            norm(difference(cached(decoded), anchor)),
            norm(cached(decoded)),
        );
        if !e.valid() || !l.valid() || !w.valid() {
            return Err(ReadError::Unsupported);
        }
        delta_p = delta_p.max(e.hi);
        radius = radius.max(l.hi);
        extent = extent.max(w.hi);
        index += 1;
        Ok(())
    })?;
    if index != points.len() {
        return Err(ReadError::Unsupported);
    }
    let (mut wmin, mut wmax, mut delta_w) = (1_f64, 1_f64, 0_f64);
    if let Some(raw_weights) = field.weights {
        let decoded = weights.ok_or(ReadError::Unsupported)?;
        if decoded.len() != points.len() {
            return Err(ReadError::Unsupported);
        }
        wmin = f64::INFINITY;
        wmax = 0.;
        index = 0;
        // Fields reads weights only from this complex-record component.
        let component = if degrees.len() == 1 {
            "RATIONAL_B_SPLINE_CURVE"
        } else {
            "RATIONAL_B_SPLINE_SURFACE"
        };
        visit(raw_weights, degrees.len(), &mut |value| {
            let raw = scalar_in(source, id, Some(component), value, budget)?;
            let decoded = *decoded.get(index).ok_or(ReadError::Unsupported)?;
            let difference = raw - I::point(decoded);
            if !raw.valid()
                || raw.lo <= 0.
                || !decoded.is_finite()
                || decoded <= 0.
                || !difference.valid()
            {
                return Err(ReadError::Unsupported);
            }
            wmin = wmin.min(raw.lo.min(decoded));
            wmax = wmax.max(raw.hi.max(decoded));
            delta_w = delta_w.max(difference.magnitude());
            index += 1;
            Ok(())
        })?;
        if index != points.len() {
            return Err(ReadError::Unsupported);
        }
    } else if weights.is_some() {
        return Err(ReadError::Unsupported);
    }
    let mut knot_term = I::point(0.);
    for axis in 0..degrees.len() {
        integer(source, id, field.degrees[axis], degrees[axis], budget)?;
        let Some((mults, values)) = field.knots[axis] else {
            if knots[axis].0 != [0., 1.] || knots[axis].1 != [degrees[axis] + 1; 2] {
                return Err(ReadError::Unsupported);
            }
            continue;
        };
        let (mults, values) = (list(mults)?, list(values)?);
        let (decoded, counts) = knots[axis];
        if values.len() != decoded.len() || mults.len() != counts.len() || values.len() < 2 {
            return Err(ReadError::Unsupported);
        }
        let (mut delta, mut gap) = (0_f64, f64::INFINITY);
        let mut previous: Option<I> = None;
        for i in 0..values.len() {
            integer(source, id, &mults[i], counts[i], budget)?;
            let mut raw = scalar(source, id, &values[i], budget)?;
            if exact_binary(
                spelling(source, id, &values[i], budget)?,
                decoded[i],
                budget,
            )? {
                raw = I::point(decoded[i]);
            }
            let error = if raw.lo == decoded[i] && raw.hi == decoded[i] {
                0.
            } else {
                (raw - I::point(decoded[i])).magnitude()
            };
            if !error.is_finite() {
                return Err(ReadError::Unsupported);
            }
            delta = delta.max(error);
            if let Some(old) = previous {
                let difference = raw - old;
                if !difference.valid() || difference.lo <= 0. {
                    return Err(ReadError::Unsupported);
                }
                gap = gap.min(
                    difference
                        .lo
                        .min((I::point(decoded[i]) - I::point(decoded[i - 1])).lo),
                );
            }
            previous = Some(raw);
        }
        if !gap.is_finite() || gap <= 0. {
            return Err(ReadError::Unsupported);
        }
        // A displaced discontinuity has no uniform small C0 bound. With
        // continuous internal joins the homotopy below is piecewise smooth.
        if delta > 0. && counts[1..counts.len() - 1].contains(&(degrees[axis] + 1)) {
            return Err(ReadError::Unsupported);
        }
        knot_term = knot_term
            + I::point(degrees[axis] as f64) * I::point(4.) * I::point(delta) / I::point(gap);
    }
    let kappa = I::point(wmax) / I::point(wmin);
    // Homotope normalized-domain t and distinct knots. On every active
    // deBoor blend |alpha'| <= (|t'|+|ki'|+2|kj'|)/g <=4delta/g.
    // Each homogeneous anchored blend contributes <=3*kappa*radius times
    // this perturbation; p levels, continuity across span crossings. Tensor
    // axes telescope. Coordinate and positive-weight perturbations are
    // separate phases, so no unproved dependency cancellation is used.
    let error = kappa * I::point(delta_p)
        + I::point(2.) * I::point(radius) * I::point(delta_w) / I::point(wmin)
        + I::point(3.) * kappa * I::point(radius) * knot_term;
    if !error.valid() || error.hi < 0. {
        return Err(ReadError::Unsupported);
    }
    Ok(NumericNurbs {
        error: error.hi,
        extent,
    })
}
fn curve_nurbs(
    source: &Source<'_>,
    id: u32,
    curve: &Curve,
    dimension: usize,
    budget: &mut Budget<'_>,
) -> Result<NumericNurbs> {
    let Curve::Nurbs {
        curve,
        closed: false,
    } = curve
    else {
        return Err(ReadError::Unsupported);
    };
    let field = fields(source, id, 1, budget)?;
    nurbs_error(
        source,
        id,
        &field,
        &[curve.degree()],
        &[curve.knots()],
        curve.points(),
        curve.weights(),
        dimension,
        None,
        budget,
    )
}

struct Iso {
    moving: usize,
    origin_id: u32,
    endpoint_id: Option<u32>,
    origin: V,
    rate: I,
}
/// Actual finite positive-rational UV segment, in ORIGINAL REAL bounds.
/// This calculation descriptor grants no infinite-support or surface-domain
/// permission; each finite caller must earn the master window separately.
pub(super) struct AxisSegment {
    pub(super) moving: usize,
    pub(super) endpoints: [V; 2],
}
pub(super) fn axis_segment(
    source: &Source<'_>,
    definition: &super::uv::Definition,
    budget: &mut Budget<'_>,
) -> Result<Option<AxisSegment>> {
    if let Curve::Nurbs {
        curve,
        closed: false,
    } = &definition.geometry
        && curve.degree() >= 2
    {
        return polynomial_axis::segment(source, definition, budget);
    }
    if !matches!(definition.geometry, Curve::Nurbs { closed: false, .. }) {
        return Ok(None);
    }
    let iso = isoline(source, definition, budget)?;
    let end = iso.endpoint_id.ok_or(ReadError::Unsupported)?;
    let endpoint = point(source, end, 2, budget)?;
    budget.charge(128)?;
    Ok(Some(AxisSegment {
        moving: iso.moving,
        endpoints: [iso.origin, endpoint],
    }))
}
fn isoline(
    source: &Source<'_>,
    definition: &super::uv::Definition,
    budget: &mut Budget<'_>,
) -> Result<Iso> {
    budget.charge(256)?;
    if let Curve::Nurbs {
        curve,
        closed: false,
    } = &definition.geometry
    {
        if curve.degree() != 1 || curve.points().len() != 2 || curve.knots().1 != [2, 2] {
            return Err(ReadError::Unsupported);
        }
        let a = curve.points()[0];
        let b = curve.points()[1];
        let moving = if a[0] != b[0] && a[1] == b[1] {
            0
        } else if a[0] == b[0] && a[1] != b[1] {
            1
        } else {
            return Err(ReadError::Unsupported);
        };
        // Validate original degree/multiplicity/order and positive weights.
        // A two-pole rational degree-one curve traverses exactly its FINITE
        // endpoint segment, monotonically, whatever those positive weights.
        curve_nurbs(source, definition.curve_id, &definition.geometry, 2, budget)?;
        let field = fields(source, definition.curve_id, 1, budget)?;
        let controls = list(field.controls)?;
        if controls.len() != 2 {
            return Err(ReadError::Unsupported);
        }
        let origin_id = reference(&controls[0])?;
        let endpoint_id = reference(&controls[1])?;
        let origin = point(source, origin_id, 2, budget)?;
        let endpoint = point(source, endpoint_id, 2, budget)?;
        let rate = endpoint[moving] - origin[moving];
        if !rate.valid()
            || !((b[moving] > a[moving] && rate.lo > 0.) || (b[moving] < a[moving] && rate.hi < 0.))
        {
            return Err(ReadError::Unsupported);
        }
        // Original fixed-coordinate equality, not equality of rounded poles.
        let [_, first] = source.args(origin_id, "CARTESIAN_POINT")? else {
            return Err(ReadError::Invalid);
        };
        let [_, last] = source.args(endpoint_id, "CARTESIAN_POINT")? else {
            return Err(ReadError::Invalid);
        };
        if !exact_same(
            source,
            (origin_id, &list(first)?[1 - moving]),
            (endpoint_id, &list(last)?[1 - moving]),
            budget,
        )? {
            return Err(ReadError::Unsupported);
        }
        return Ok(Iso {
            moving,
            origin_id,
            endpoint_id: Some(endpoint_id),
            origin,
            rate,
        });
    }
    let Curve::Line {
        direction: decoded, ..
    } = definition.geometry
    else {
        return Err(ReadError::Unsupported);
    };
    let moving = if decoded[0] != 0. && decoded[1] == 0. {
        0
    } else if decoded[0] == 0. && decoded[1] != 0. {
        1
    } else {
        return Err(ReadError::Unsupported);
    };
    let [_, origin, vector_id] = source.args(definition.curve_id, "LINE")? else {
        return Err(ReadError::Invalid);
    };
    let origin_id = reference(origin)?;
    let [_, direction_id, _] = source.args(reference(vector_id)?, "VECTOR")? else {
        return Err(ReadError::Invalid);
    };
    let direction_id = reference(direction_id)?;
    let [_, values] = source.args(direction_id, "DIRECTION")? else {
        return Err(ReadError::Invalid);
    };
    let values = list(values)?;
    if values.len() != 2 {
        return Err(ReadError::Unsupported);
    }
    let fixed = scalar(source, direction_id, &values[1 - moving], budget)?;
    if fixed.lo != 0. || fixed.hi != 0. {
        return Err(ReadError::Unsupported);
    }
    let rate = vector(source, reference(vector_id)?, 2, budget)?[moving];
    if !rate.valid()
        || !((decoded[moving] > 0. && rate.lo > 0.) || (decoded[moving] < 0. && rate.hi < 0.))
    {
        return Err(ReadError::Unsupported);
    }
    Ok(Iso {
        moving,
        origin_id,
        endpoint_id: None,
        origin: point(source, origin_id, 2, budget)?,
        rate,
    })
}

fn point_boundary(
    source: &Source<'_>,
    id: u32,
    field: &Fields<'_, '_>,
    axis: usize,
    end: usize,
    origin_id: u32,
    budget: &mut Budget<'_>,
) -> Result<()> {
    let [_, coordinates] = source.args(origin_id, "CARTESIAN_POINT")? else {
        return Err(ReadError::Invalid);
    };
    let coordinate = list(coordinates)?.get(axis).ok_or(ReadError::Unsupported)?;
    if let Some((_, knots)) = field.knots[axis] {
        let knots = list(knots)?;
        let knot = if end == 0 {
            knots.first()
        } else {
            knots.last()
        }
        .ok_or(ReadError::Unsupported)?;
        if !exact_same(source, (origin_id, coordinate), (id, knot), budget)? {
            return Err(ReadError::Unsupported);
        }
    } else {
        let value = scalar(source, origin_id, coordinate, budget)?;
        let text = spelling(source, origin_id, coordinate, budget)?;
        if !value.valid() || !exact_binary(text, end as f64, budget)? {
            return Err(ReadError::Unsupported);
        }
    }
    Ok(())
}
fn boundary(
    source: &Source<'_>,
    id: u32,
    field: &Fields<'_, '_>,
    axis: usize,
    end: usize,
    iso: &Iso,
    budget: &mut Budget<'_>,
) -> Result<()> {
    point_boundary(source, id, field, axis, end, iso.origin_id, budget)?;
    if let Some(last) = iso.endpoint_id {
        point_boundary(source, id, field, axis, end, last, budget)?;
    }
    Ok(())
}

fn surface_nurbs(
    source: &Source<'_>,
    definition: &super::uv::Definition,
    surface: &Surface,
    budget: &mut Budget<'_>,
    cache: &mut Cache<'_, '_>,
) -> Result<f64> {
    let Surface::Nurbs { surface, closed } = surface else {
        return Err(ReadError::Unsupported);
    };
    if closed.iter().any(|x| *x) {
        return Err(ReadError::Unsupported);
    }
    let iso = isoline(source, definition, budget)?;
    let fixed = 1 - iso.moving;
    let origin = match &definition.geometry {
        Curve::Line { origin, .. } => *origin,
        Curve::Nurbs { curve, .. } => curve.points()[0],
        _ => return Err(ReadError::Unsupported),
    };
    let range = surface.range(fixed);
    let end = if origin[fixed] == range[0] {
        0
    } else if origin[fixed] == range[1] {
        1
    } else {
        return Err(ReadError::Unsupported);
    };
    let field = fields(source, definition.surface, 2, budget)?;
    boundary(source, definition.surface, &field, fixed, end, &iso, budget)?;
    if let Some(last) = iso.endpoint_id {
        let Curve::Nurbs { curve, .. } = &definition.geometry else {
            return Err(ReadError::Unsupported);
        };
        let moving_range = surface.range(iso.moving);
        let a = curve.points()[0][iso.moving];
        let b = curve.points()[1][iso.moving];
        let forward = if a == moving_range[0] && b == moving_range[1] {
            true
        } else if a == moving_range[1] && b == moving_range[0] {
            false
        } else {
            return Err(ReadError::Unsupported);
        };
        point_boundary(
            source,
            definition.surface,
            &field,
            iso.moving,
            usize::from(!forward),
            iso.origin_id,
            budget,
        )?;
        point_boundary(
            source,
            definition.surface,
            &field,
            iso.moving,
            usize::from(forward),
            last,
            budget,
        )?;
    }
    let degrees = surface.degrees();
    let knots = [surface.knots(0), surface.knots(1)];
    let (_, multiplicities) = knots[fixed];
    if multiplicities[0] != degrees[fixed] + 1
        || multiplicities[multiplicities.len() - 1] != degrees[fixed] + 1
    {
        return Err(ReadError::Unsupported);
    }
    if let Some(error) = cache.get(source, definition.surface, budget)? {
        return Ok(error);
    }
    let error = nurbs_error(
        source,
        definition.surface,
        &field,
        &degrees,
        &knots,
        surface.points(),
        surface.weights(),
        3,
        None,
        budget,
    )?
    .error;
    cache.remember(source, definition.surface, error, budget)?;
    Ok(error)
}
fn plane_nurbs(
    source: &Source<'_>,
    definition: &super::uv::Definition,
    surface: &Surface,
    budget: &mut Budget<'_>,
) -> Result<f64> {
    let Surface::Plane(decoded) = surface else {
        return Err(ReadError::Unsupported);
    };
    let uv = curve_nurbs(source, definition.curve_id, &definition.geometry, 2, budget)?;
    let [_, axis] = source.args(definition.surface, "PLANE")? else {
        return Err(ReadError::Invalid);
    };
    let raw = frame(source, reference(axis)?, 3, budget)?;
    let error = norm(difference(raw.origin, cached(decoded.origin)))
        + I::point(uv.extent)
            * (norm(difference(raw.axes[0], cached(decoded.x)))
                + norm(difference(raw.axes[1], cached(decoded.y))))
        + (norm(raw.axes[0]) + norm(raw.axes[1])) * I::point(uv.error);
    if !error.valid() {
        return Err(ReadError::Unsupported);
    }
    Ok(error.hi)
}
fn extrusion(
    source: &Source<'_>,
    definition: &super::uv::Definition,
    surface: &Surface,
    edge: &Edge,
    budget: &mut Budget<'_>,
) -> Result<f64> {
    let Surface::Extrusion {
        curve,
        direction: decoded_axis,
    } = surface
    else {
        return Err(ReadError::Unsupported);
    };
    let [_, generator, direction_id] =
        source.args(definition.surface, "SURFACE_OF_LINEAR_EXTRUSION")?
    else {
        return Err(ReadError::Invalid);
    };
    let generator = unwrapped(source, reference(generator)?, budget)?;
    // Initial raw-qualified sweep class: actual open retained NURBS. Other
    // generators remain typed refusals until their REAL-domain transfer earns
    // the same independent bound; the decoded sweep helper alone is not it.
    let gen_error = curve_nurbs(source, generator, curve, 3, budget)?;
    let raw_axis = vector(source, reference(direction_id)?, 3, budget)?;
    let iso = isoline(source, definition, budget)?;
    if iso.endpoint_id.is_some() {
        return Err(ReadError::Unsupported);
    }
    let Curve::Line {
        origin: uv_origin,
        direction: uv_direction,
    } = &definition.geometry
    else {
        return Err(ReadError::Unsupported);
    };
    if iso.moving == 0 {
        let shift: V = std::array::from_fn(|i| {
            raw_axis[i] * iso.origin[1] - I::point(decoded_axis[i]) * I::point(uv_origin[1])
        });
        let error = I::point(gen_error.error) + norm(shift);
        if !error.valid() {
            return Err(ReadError::Unsupported);
        }
        return Ok(error.hi);
    }
    let Curve::Nurbs {
        curve: generator_curve,
        closed: false,
    } = curve
    else {
        return Err(ReadError::Unsupported);
    };
    let range = generator_curve.range();
    let end = if uv_origin[0] == range[0] {
        0
    } else if uv_origin[0] == range[1] {
        1
    } else {
        return Err(ReadError::Unsupported);
    };
    let field = fields(source, generator, 1, budget)?;
    let (_, mult) = generator_curve.knots();
    if mult[0] != generator_curve.degree() + 1
        || mult[mult.len() - 1] != generator_curve.degree() + 1
    {
        return Err(ReadError::Unsupported);
    }
    boundary(source, generator, &field, 0, end, &iso, budget)?;
    let poles = list(field.controls)?;
    let pole = if end == 0 {
        poles.first()
    } else {
        poles.last()
    }
    .ok_or(ReadError::Unsupported)?;
    let raw_point = point(source, reference(pole)?, 3, budget)?;
    let cached_point = super::affine::at(curve, uv_origin[0], budget)?;
    let (raw_origin, cached_origin, raw_direction, cached_direction) = (
        std::array::from_fn(|i| raw_point[i] + raw_axis[i] * iso.origin[1]),
        std::array::from_fn(|i| {
            cached_point[i] + I::point(decoded_axis[i]) * I::point(uv_origin[1])
        }),
        raw_axis.map(|x| x * iso.rate),
        cached(*decoded_axis).map(|x| x * I::point(uv_direction[1])),
    );
    analytic::line_image_reserve(
        super::analytic::Image::Line {
            origin: raw_origin,
            direction: raw_direction,
        },
        super::analytic::Image::Line {
            origin: cached_origin,
            direction: cached_direction,
        },
        edge,
        budget,
    )
}
/// Earn a numeric transfer from original REAL leaves, not a tolerance-derived
/// allowance. The caller adds it to a positive decoded certificate or subtracts
/// it from a geometric disagreement lower bound, keeping the original limit.
pub(super) fn reserve(
    source: &Source<'_>,
    curve3d_id: u32,
    definition: &super::uv::Definition,
    surface: &Surface,
    edge: &Edge,
    budget: &mut Budget<'_>,
    cache: &mut Cache<'_, '_>,
) -> Result<f64> {
    budget.charge(256)?;
    let curve3d_id = unwrapped(source, curve3d_id, budget)?;
    let edge_error = edge_curve_error(source, curve3d_id, edge, budget)?;
    let image_error = match (surface, &definition.geometry) {
        (Surface::Plane(_), Curve::Nurbs { .. }) => {
            plane_nurbs(source, definition, surface, budget)?
        }
        (Surface::Nurbs { .. }, _) => surface_nurbs(source, definition, surface, budget, cache)?,
        (Surface::Extrusion { .. }, _) => extrusion(source, definition, surface, edge, budget)?,
        // Raw revolution transfer is deliberately not inferred from the
        // decoded Rodrigues certificate; its independent REAL stage is open.
        (Surface::Revolution { .. }, _) => return Err(ReadError::Unsupported),
        _ => analytic::image_reserve(source, definition, surface, edge, budget)?,
    };
    let total = I::point(edge_error) + I::point(image_error);
    if !total.valid() || total.hi < 0. {
        return Err(ReadError::Unsupported);
    }
    budget.charge(64)?;
    Ok(total.hi)
}

/// Most NURBS edge-curve transfer bounds one pass keeps.
const MAX_CURVE_CACHE: usize = 4096;
/// The numeric transfer bound of the edge's own 3D curve. Each edge is
/// bound to two pcurves; a NURBS curve's whole-control bound depends only on
/// the immutable source and the curve entity, so it is computed once per
/// pass that turned the cache on. Every trim/domain/sense check still runs.
fn edge_curve_error(
    source: &Source<'_>,
    curve3d_id: u32,
    edge: &Edge,
    budget: &mut Budget<'_>,
) -> Result<f64> {
    if !matches!(edge.curve, Curve::Nurbs { .. }) {
        return analytic::curve_reserve(source, curve3d_id, edge, budget);
    }
    if let Some(cache) = source.curve_errors.borrow().as_ref() {
        budget.charge(cache.len().saturating_mul(8).saturating_add(64))?;
        if let Some(&(_, error)) = cache.iter().find(|&&(id, _)| id == curve3d_id) {
            return Ok(error);
        }
    }
    let error = curve_nurbs(source, curve3d_id, &edge.curve, 3, budget)?.error;
    if !error.is_finite() || error < 0. {
        return Err(ReadError::Unsupported);
    }
    if let Some(cache) = source.curve_errors.borrow_mut().as_mut()
        && cache.len() < MAX_CURVE_CACHE
    {
        budget.charge(cache.len().saturating_mul(8).saturating_add(128))?;
        cache.try_reserve(1).map_err(|_| ReadError::Limit)?;
        cache.push((curve3d_id, error));
    }
    Ok(error)
}

#[cfg(test)]
mod tests;

fn cardinal(v: V) -> Option<[f64; 3]> {
    if !v.into_iter().all(I::valid) {
        return None;
    }
    let mut out = [0.; 3];
    let mut count = 0;
    for i in 0..3 {
        if v[i].lo == 0. && v[i].hi == 0. {
            continue;
        }
        out[i] = if v[i].lo > 0. {
            1.
        } else if v[i].hi < 0. {
            -1.
        } else {
            return None;
        };
        count += 1;
    }
    (count == 1).then_some(out)
}
fn bigint_work(bits: usize) -> usize {
    let limbs = bits / 64 + 1;
    limbs
        .saturating_mul(limbs)
        .saturating_mul(bits + 1)
        .saturating_mul(8)
        .saturating_add(4096)
}
/// Bounded exact source decimal, used only for earned affine cancellations.
/// Ordinary coordinate/control enclosures take the cheaper adjacent path.
pub(super) fn rational(
    source: &Source<'_>,
    id: u32,
    value: &Value<'_>,
    budget: &mut Budget<'_>,
) -> Result<BigRational> {
    let text = spelling(source, id, value, budget)?;
    budget.charge(text.len().saturating_add(1024))?;
    if let Ok(q) = Rational::parse_decimal(text) {
        let (n, d) = q.ratio();
        return Ok(BigRational::new(BigInt::from(n), BigInt::from(d)));
    }
    if text.len() > 256 {
        return Err(ReadError::Unsupported);
    }
    let (mantissa, exponent) = text.split_once('E').unwrap_or((text, "0"));
    let exponent = exponent
        .parse::<i32>()
        .map_err(|_| ReadError::Unsupported)?;
    let negative = mantissa.starts_with('-');
    let mantissa = mantissa.trim_start_matches(['+', '-']);
    let fraction = mantissa.split_once('.').map_or(0, |(_, f)| f.len());
    let digits: String = mantissa.chars().filter(|c| *c != '.').collect();
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return Err(ReadError::Unsupported);
    }
    let shift = i64::from(exponent) - fraction as i64;
    let decimal_size = if shift >= 0 {
        (digits.len() as i64).saturating_add(shift)
    } else {
        (digits.len() as i64).max(-shift)
    };
    if decimal_size > 1023 {
        return Err(ReadError::Unsupported);
    }
    let bits = 4 * (decimal_size as usize + 1);
    budget.charge(bigint_work(bits))?;
    let mut numerator = BigInt::parse_bytes(digits.as_bytes(), 10).ok_or(ReadError::Unsupported)?;
    if negative {
        numerator = -numerator;
    }
    let mut denominator = BigInt::from(1);
    if shift >= 0 {
        numerator *= BigInt::from(10).pow(shift as u32);
    } else {
        denominator = BigInt::from(10).pow((-shift) as u32);
    }
    Ok(BigRational::new(numerator, denominator))
}
fn float_bits(value: f64) -> usize {
    if value == 0. {
        return 1;
    }
    let exponent = ((value.to_bits() >> 52) & 2047) as i32;
    let exponent = if exponent == 0 {
        -1022
    } else {
        exponent - 1023
    };
    (53 + (exponent - 52).max(0) as usize).max(1 + (52 - exponent).max(0) as usize)
}
pub(super) fn enclose(q: &BigRational, budget: &mut Budget<'_>) -> Result<I> {
    budget.charge(128)?;
    if q.denom() == &BigInt::from(1)
        && let Some(integer) = q.numer().to_i64()
        && integer.unsigned_abs() <= 1_u64 << 53
    {
        return Ok(I::point(integer as f64));
    }
    let bits = q.numer().bits().max(q.denom().bits()) as usize;
    if bits > 4096 {
        return Err(ReadError::Unsupported);
    }
    budget.charge(bigint_work(bits))?;
    let value = q.to_f64().ok_or(ReadError::Unsupported)?;
    if !value.is_finite() {
        return Err(ReadError::Unsupported);
    }
    budget.charge(bigint_work(bits.max(float_bits(value))))?;
    if BigRational::from_float(value).as_ref() == Some(q) {
        return Ok(I::point(value));
    }
    let (mut lo, mut hi) = (value.next_down(), value.next_up());
    for _ in 0..4 {
        budget.charge(bigint_work(bits.max(float_bits(lo)).max(float_bits(hi))))?;
        let (Some(a), Some(b)) = (BigRational::from_float(lo), BigRational::from_float(hi)) else {
            return Err(ReadError::Unsupported);
        };
        if &a <= q && q <= &b {
            return Ok(I::new(lo, hi));
        }
        if &a > q {
            lo = lo.next_down();
        }
        if &b < q {
            hi = hi.next_up();
        }
    }
    Err(ReadError::Unsupported)
}
/// Exact original decimal plane origin plus UV coordinates for a separately
/// earned cardinal frame. Round only after cancellation, retaining the gauge.
#[cfg(test)]
pub(super) fn translated_point(
    source: &Source<'_>,
    origin_id: u32,
    uv_id: u32,
    axes: [V; 3],
    budget: &mut Budget<'_>,
) -> Result<V> {
    let origin = translation_origin(source, origin_id, budget)?;
    translated_from(source, &origin, uv_id, axes, budget)
}
/// The exact REAL coordinates of a translated frame's origin point.
pub(super) fn translation_origin(
    source: &Source<'_>,
    origin_id: u32,
    budget: &mut Budget<'_>,
) -> Result<[BigRational; 3]> {
    let [_, origin] = source.args(origin_id, "CARTESIAN_POINT")? else {
        return Err(ReadError::Invalid);
    };
    let origin = list(origin)?;
    if origin.len() != 3 {
        return Err(ReadError::Unsupported);
    }
    Ok([
        rational(source, origin_id, &origin[0], budget)?,
        rational(source, origin_id, &origin[1], budget)?,
        rational(source, origin_id, &origin[2], budget)?,
    ])
}
/// [`translated_point`] with its origin's exact coordinates already read.
pub(super) fn translated_from(
    source: &Source<'_>,
    origin: &[BigRational; 3],
    uv_id: u32,
    axes: [V; 3],
    budget: &mut Budget<'_>,
) -> Result<V> {
    let mut coefficients = [[0_i32; 3]; 2];
    for axis in 0..2 {
        for i in 0..3 {
            let c = axes[axis][i];
            if c.lo != c.hi || ![-1., 0., 1.].contains(&c.lo) {
                return Err(ReadError::Unsupported);
            }
            coefficients[axis][i] = c.lo as i32;
        }
    }
    let [_, uv] = source.args(uv_id, "CARTESIAN_POINT")? else {
        return Err(ReadError::Invalid);
    };
    let uv = list(uv)?;
    if uv.len() != 2 {
        return Err(ReadError::Unsupported);
    }
    let q = [
        rational(source, uv_id, &uv[0], budget)?,
        rational(source, uv_id, &uv[1], budget)?,
    ];
    let mut result = [I::point(0.); 3];
    for i in 0..3 {
        budget.charge(64)?;
        let mut sum = origin[i].clone();
        for axis in 0..2 {
            if coefficients[axis][i] == 0 || q[axis].is_zero() {
                continue;
            }
            let bits = sum
                .numer()
                .bits()
                .max(sum.denom().bits())
                .max(q[axis].numer().bits())
                .max(q[axis].denom().bits()) as usize;
            if bits > 4096 {
                return Err(ReadError::Unsupported);
            }
            if sum.is_zero() {
                budget.charge(bits.saturating_mul(2).saturating_add(256))?;
                sum = if coefficients[axis][i] == 1 {
                    q[axis].clone()
                } else {
                    -q[axis].clone()
                };
                continue;
            }
            budget.charge(bigint_work(bits.saturating_mul(2)))?;
            if coefficients[axis][i] == 1 {
                sum += &q[axis];
            } else if coefficients[axis][i] == -1 {
                sum -= &q[axis];
            }
        }
        result[i] = enclose(&sum, budget)?;
    }
    Ok(result)
}

/// Separate from reserve: an uncertain sense may block a positive certificate
/// but must never erase an independently earned point-set disagreement.
pub(super) fn sense(
    source: &Source<'_>,
    curve3d_id: u32,
    definition: &super::uv::Definition,
    surface: &Surface,
    edge: &Edge,
    budget: &mut Budget<'_>,
) -> Result<()> {
    budget.charge(256)?;
    match (surface, &definition.geometry, &edge.curve) {
        (
            Surface::Plane(_),
            Curve::Nurbs { closed: false, .. },
            Curve::Nurbs { closed: false, .. },
        ) => Ok(()),
        (Surface::Nurbs { .. }, _, Curve::Line { .. }) => {
            sense_line::certify(source, curve3d_id, definition, surface, edge, budget)
        }
        (Surface::Nurbs { .. }, _, Curve::Nurbs { closed: false, .. }) => {
            isoline(source, definition, budget)?;
            Ok(())
        }
        (Surface::Extrusion { .. }, _, Curve::Nurbs { closed: false, .. }) => {
            let iso = isoline(source, definition, budget)?;
            if iso.moving != 0 {
                return Err(ReadError::Unsupported);
            }
            Ok(())
        }
        (Surface::Extrusion { .. }, _, Curve::Line { .. }) => {
            let iso = isoline(source, definition, budget)?;
            if iso.moving != 1 {
                return Err(ReadError::Unsupported);
            }
            let [_, _, axis] = source.args(definition.surface, "SURFACE_OF_LINEAR_EXTRUSION")?
            else {
                return Err(ReadError::Invalid);
            };
            let d = vector(source, reference(axis)?, 3, budget)?.map(|x| x * iso.rate);
            let id = unwrapped(source, curve3d_id, budget)?;
            let [_, _, axis] = source.args(id, "LINE")? else {
                return Err(ReadError::Invalid);
            };
            let other = vector(source, reference(axis)?, 3, budget)?;
            let dot = (0..3).fold(I::point(0.), |s, i| s + d[i] * other[i]);
            if !dot.valid() || dot.lo <= 0. {
                return Err(ReadError::Unsupported);
            }
            Ok(())
        }
        _ => analytic::sense(source, curve3d_id, definition, surface, edge, budget),
    }
}

/// Direct ORIGINAL-REAL plane/UV-NURBS certificate. Exact affine cancellation
/// precedes binary enclosure; the large gauge is never a geometric repair.
/// This result already includes the complete numeric transfer. Do not apply
/// reserve() again to it or reinterpret it as a decoded-UV certificate.
pub(super) fn plane_nurbs_certificate(
    source: &Source<'_>,
    curve3d_id: u32,
    definition: &super::uv::Definition,
    surface: &Surface,
    edge: &Edge,
    tol: f64,
    budget: &mut Budget<'_>,
) -> Result<Option<super::Proof>> {
    // An undecided new upper must preserve the existing whole-control
    // separating witness. Unsupported is not evidence of disagreement.
    match conic::certify(source, curve3d_id, definition, surface, edge, tol, budget) {
        Ok(Some(proof @ (super::Proof::Certified(_) | super::Proof::Disagreement(_)))) => {
            return Ok(Some(proof));
        }
        Ok(_) | Err(ReadError::Unsupported) => {}
        Err(error) => return Err(error),
    }
    let (
        Surface::Plane(_),
        Curve::Nurbs {
            curve,
            closed: false,
        },
    ) = (surface, &definition.geometry)
    else {
        return Ok(None);
    };
    budget.charge(512)?;
    let field = fields(source, definition.curve_id, 1, budget)?;
    let controls = list(field.controls)?;
    if controls.len() != curve.points().len() {
        return Err(ReadError::Unsupported);
    }
    let bytes = controls
        .len()
        .saturating_mul(std::mem::size_of::<V>() + std::mem::size_of::<[f64; 3]>());
    if bytes > 8 * 1024 * 1024 {
        return Err(ReadError::Limit);
    }
    let (mut raw, mut centers) = (Vec::new(), Vec::new());
    raw.try_reserve_exact(controls.len())
        .map_err(|_| ReadError::Limit)?;
    centers
        .try_reserve_exact(controls.len())
        .map_err(|_| ReadError::Limit)?;
    if raw
        .capacity()
        .saturating_mul(std::mem::size_of::<V>())
        .saturating_add(
            centers
                .capacity()
                .saturating_mul(std::mem::size_of::<[f64; 3]>()),
        )
        > 8 * 1024 * 1024
    {
        return Err(ReadError::Limit);
    }
    let plane = analytic::plane_lift(source, definition.surface, budget)?;
    for control in controls {
        budget.charge(256)?;
        let p = plane.point(source, reference(control)?, budget)?;
        if !p.into_iter().all(I::valid) {
            return Err(ReadError::Unsupported);
        }
        raw.push(p);
        centers.push(p.map(I::midpoint));
    }
    let numerical = nurbs_error(
        source,
        definition.curve_id,
        &field,
        &[curve.degree()],
        &[curve.knots()],
        &centers,
        curve.weights(),
        2,
        Some(&raw),
        budget,
    )?
    .error;
    let curve3d_id = unwrapped(source, curve3d_id, budget)?;
    let target = edge_curve_error(source, curve3d_id, edge, budget)?;
    let floor = I::point(numerical) + I::point(target);
    if !floor.valid() {
        return Err(ReadError::Unsupported);
    }
    let proof = super::affine::nurbs(
        curve.degree(),
        curve.knots(),
        curve.weights(),
        &raw,
        true,
        edge,
        tol,
        budget,
    )?;
    let proof = match proof {
        super::Proof::Certified(upper) => {
            // Positive natural-domain knot maps preserve the UV/3D parameter
            // sense. The two-control LINE case additionally proves original
            // physical directions agree independently of the C0 tolerance.
            if matches!(edge.curve, Curve::Line { .. }) {
                if curve.degree() != 1 || raw.len() != 2 {
                    return Ok(Some(super::Proof::Unsupported));
                }
                let [_, _, direction] = source.args(curve3d_id, "LINE")? else {
                    return Err(ReadError::Invalid);
                };
                let d = vector(source, reference(direction)?, 3, budget)?;
                let delta = difference(raw[1], raw[0]);
                let dot = (0..3).fold(I::point(0.), |s, i| s + delta[i] * d[i]);
                if !dot.valid() || dot.lo <= 0. {
                    return Ok(Some(super::Proof::Unsupported));
                }
            }
            let total = I::point(upper) + floor;
            if total.valid() && total.hi <= tol {
                super::Proof::Certified(total.hi)
            } else {
                super::Proof::Unsupported
            }
        }
        super::Proof::Disagreement(lower) => {
            let remaining = I::point(lower) - floor;
            if remaining.valid() && remaining.lo > tol {
                super::Proof::Disagreement(remaining.lo)
            } else {
                super::Proof::Unsupported
            }
        }
        super::Proof::Unsupported => super::Proof::Unsupported,
    };
    budget.charge(64)?;
    Ok(Some(proof))
}

/// Direct original-plane analytic certificate. Like the NURBS route this
/// verifies raw REAL coefficients before cancellation and includes numeric
/// transfer exactly once; it never replaces the body's retained geometry.
pub(super) fn plane_analytic_certificate(
    source: &Source<'_>,
    curve3d_id: u32,
    definition: &super::uv::Definition,
    surface: &Surface,
    edge: &Edge,
    tol: f64,
    budget: &mut Budget<'_>,
) -> Result<Option<super::Proof>> {
    if !matches!(surface, Surface::Plane(_)) || matches!(definition.geometry, Curve::Nurbs { .. }) {
        return Ok(None);
    }
    budget.charge(256)?;
    let image = analytic::actual_image(source, definition, surface, budget)?;
    let curve3d_id = unwrapped(source, curve3d_id, budget)?;
    let (target_image, target) = if matches!(edge.curve, Curve::Nurbs { .. }) {
        (
            None,
            curve_nurbs(source, curve3d_id, &edge.curve, 3, budget)?.error,
        )
    } else {
        let (image, reserve) = analytic::curve_packet(source, curve3d_id, edge, budget)?;
        (Some(image), reserve)
    };
    let proof = super::analytic::certify_image(image, edge, tol, budget)?;
    let proof = match proof {
        super::Proof::Certified(upper) => {
            if let Some(target) = target_image {
                analytic::sense_images(target, image, budget)?;
            } else {
                analytic::sense(source, curve3d_id, definition, surface, edge, budget)?;
            }
            let total = I::point(upper) + I::point(target);
            if total.valid() && total.hi <= tol {
                super::Proof::Certified(total.hi)
            } else {
                super::Proof::Unsupported
            }
        }
        super::Proof::Disagreement(lower) => {
            let remaining = I::point(lower) - I::point(target);
            if remaining.valid() && remaining.lo > tol {
                super::Proof::Disagreement(remaining.lo)
            } else {
                super::Proof::Unsupported
            }
        }
        super::Proof::Unsupported => super::Proof::Unsupported,
    };
    budget.charge(64)?;
    Ok(Some(proof))
}

pub(super) fn analytic_image(
    source: &Source<'_>,
    definition: &super::uv::Definition,
    surface: &Surface,
    budget: &mut Budget<'_>,
) -> Result<super::analytic::Image> {
    analytic::actual_image(source, definition, surface, budget)
}
pub(super) fn analytic_curve_reserve(
    source: &Source<'_>,
    curve3d_id: u32,
    edge: &Edge,
    budget: &mut Budget<'_>,
) -> Result<f64> {
    analytic::curve_reserve(source, curve3d_id, edge, budget)
}
