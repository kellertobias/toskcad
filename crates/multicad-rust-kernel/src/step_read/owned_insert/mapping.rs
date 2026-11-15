//! Actual complex CARTESIAN_TRANSFORMATION_OPERATOR_3D semantics.
//! The interval basis follows EXPRESS base_axis, first_proj_axis and
//! second_proj_axis (STEP Tools merged schema). This proves an operator;
//! represented-member ownership and source correspondence remain separate.
//! https://www.steptools.com/stds/stp_aim/html/t_base_axis.html
use super::*;
use crate::tapered_thread::interval::I;
const SCRATCH: usize = 8192;
type V = [I; 3];
struct Direction {
    values: V,
    zero: [bool; 3],
}
pub struct PhysicalMapping<'p, 'input> {
    parsed: &'p Parsed<'input>,
    operator: u32,
    origin: V,
    axes: [V; 3],
    scale: I,
    work: usize,
    bytes: usize,
}
impl<'p, 'input> PhysicalMapping<'p, 'input> {
    pub fn parsed(&self) -> &'p Parsed<'input> {
        self.parsed
    }
    pub fn operator_entity(&self) -> u32 {
        self.operator
    }
    pub fn origin_mm(&self) -> [[f64; 2]; 3] {
        self.origin.map(bounds)
    }
    /// Columns of the actual normalized EXPRESS basis, before uniform scale.
    pub fn axes(&self) -> [[[f64; 2]; 3]; 3] {
        self.axes.map(|v| v.map(bounds))
    }
    /// Oriented material consumers must earn this sign separately; a mapped
    /// operator may reflect space even though its uniform scale is positive.
    pub fn handedness(&self) -> Result<i8> {
        basis_handedness(self.axes)
    }
    pub fn scale_bounds(&self) -> [f64; 2] {
        bounds(self.scale)
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    /// Retained packet, receipt and fixed query scratch; not process RSS.
    pub fn query_budget_bytes(&self) -> usize {
        self.bytes
    }
    /// Box transport only; this does not admit a represented material member.
    pub fn map_box_mm(&self, input: [[f64; 2]; 3]) -> Result<[[f64; 2]; 3]> {
        let input = input.map(|[lo, hi]| I::new(lo, hi));
        for i in input {
            finite(i).map_err(|_| Error::Invalid)?;
        }
        let mut output = self.origin;
        for (i, value) in output.iter_mut().enumerate() {
            for (j, coordinate) in input.into_iter().enumerate() {
                *value = finite(*value + self.scale * self.axes[j][i] * coordinate)?;
            }
        }
        Ok(output.map(bounds))
    }
    /// Projection plus normalization earns a unit orthogonal basis, possibly
    /// reflected. Its absolute determinant is exactly one; volume scales by
    /// scl^3, independently of interval dependency in a determinant expansion.
    pub fn volume_scale_bounds(&self) -> Result<[f64; 2]> {
        let value = finite(self.scale.square() * self.scale)?;
        Ok([value.lo.max(0.), value.hi])
    }
}
fn bounds(i: I) -> [f64; 2] {
    [i.lo, i.hi]
}
fn finite(i: I) -> Result<I> {
    if i.valid() {
        Ok(i)
    } else {
        Err(Error::Unsupported)
    }
}
// Only earned unit-coordinate boxes and their cross products reach this dot
// helper. Components are bounded by two, so all sums/products stay finite.
fn dot(a: V, b: V) -> I {
    a.into_iter()
        .zip(b)
        .fold(I::point(0.), |s, (a, b)| s + a * b)
}
fn basis_handedness([x, y, z]: [V; 3]) -> Result<i8> {
    let cross = [
        z[1] * x[2] - z[2] * x[1],
        z[2] * x[0] - z[0] * x[2],
        z[0] * x[1] - z[1] * x[0],
    ];
    let sign = dot(y, cross);
    if sign.lo > 0. {
        Ok(1)
    } else if sign.hi < 0. {
        Ok(-1)
    } else {
        Err(Error::Unsupported)
    }
}
fn normalize(v: V, work: &mut Ledger<'_>) -> Result<V> {
    work.step(1024)?;
    let magnitude = v.iter().map(|v| v.magnitude()).fold(0_f64, f64::max);
    if !magnitude.is_finite() || magnitude <= 0. {
        return Err(Error::Unsupported);
    }
    let scaled = v.map(|v| v / I::point(magnitude));
    // Scaling bounds each component by one. The sum of three squares is
    // finite; its positive lower bound earns a nonzero actual vector.
    let norm = scaled.into_iter().fold(I::point(0.), |s, x| s + x.square());
    if norm.lo <= 0. {
        return Err(Error::Unsupported);
    }
    let norm = norm.sqrt();
    // sqrt(min_subnormal)>1e-162 bounds the finite quotient. Unit length
    // independently earns this coordinate intersection, even for a broad
    // numerical normalization enclosure.
    let result = scaled.map(|v| {
        let q = v / norm;
        I::new(q.lo.max(-1.), q.hi.min(1.))
    });
    Ok(result)
}
fn remove(v: V, axis: V) -> V {
    let component = dot(v, axis);
    std::array::from_fn(|i| v[i] - component * axis[i])
}
fn component<'a, 'input>(records: &'a [Entity<'input>], name: &str) -> Result<&'a Entity<'input>> {
    records
        .iter()
        .find(|e| e.name == name)
        .ok_or(Error::Unsupported)
}
fn triple(
    records: &BTreeMap<u32, Vec<Entity<'_>>>,
    id: u32,
    name: &str,
    work: &mut Ledger<'_>,
) -> Result<Direction> {
    let e = kind(records, id, name)?;
    let [Value::String(_), values] = e.args.as_slice() else {
        return Err(Error::Invalid);
    };
    let [a, b, c] = list(values)? else {
        return Err(Error::Unsupported);
    };
    let mut numeric = Numeric::new(e, work)?;
    let mut values = [I::point(0.); 3];
    let mut zero = [false; 3];
    for (i, v) in [a, b, c].into_iter().enumerate() {
        let scalar = numeric.scalar(v)?;
        work.step(3 * scalar.spelling().len() + 128)?;
        let [lo, hi] = scalar.bounds();
        values[i] = finite(I::new(lo, hi))?;
        zero[i] = scalar.is_zero();
    }
    numeric.finish()?;
    Ok(Direction { values, zero })
}
fn direction(
    records: &BTreeMap<u32, Vec<Entity<'_>>>,
    value: &Value<'_>,
    default: V,
    work: &mut Ledger<'_>,
) -> Result<Direction> {
    if matches!(value, Value::Omitted) {
        Ok(Direction {
            values: default,
            zero: default.map(|i| i.lo == 0. && i.hi == 0.),
        })
    } else {
        triple(records, reference(value)?, "DIRECTION", work)
    }
}
fn unit(index: usize) -> V {
    std::array::from_fn(|i| I::point(f64::from(i == index)))
}
// One schema evaluator shared by typed source wrappers. The caller owns
// source/member correspondence, inherited work and its own storage scope.
struct MappingValues {
    origin: V,
    axes: [V; 3],
    scale: I,
}
fn evaluate_operator(
    entities: &BTreeMap<u32, Vec<Entity<'_>>>,
    operator: u32,
    work: &mut Ledger<'_>,
) -> Result<MappingValues> {
    let records = entities.get(&operator).ok_or(Error::Invalid)?;
    if records.len() != 5 {
        return Err(Error::Unsupported);
    }
    let base = component(records, "CARTESIAN_TRANSFORMATION_OPERATOR")?;
    let dimensional = component(records, "CARTESIAN_TRANSFORMATION_OPERATOR_3D")?;
    let item = component(records, "REPRESENTATION_ITEM")?;
    let geometry = component(records, "GEOMETRIC_REPRESENTATION_ITEM")?;
    let transformation = component(records, "FUNCTIONALLY_DEFINED_TRANSFORMATION")?;
    if !matches!(item.args.as_slice(), [Value::String(_)])
        || !geometry.args.is_empty()
        || !matches!(
            transformation.args.as_slice(),
            [Value::String(_), Value::String(_)]
        )
    {
        return Err(Error::Unsupported);
    }
    let [axis1, axis2, origin, scale] = base.args.as_slice() else {
        return Err(Error::Unsupported);
    };
    let [axis3] = dimensional.args.as_slice() else {
        return Err(Error::Unsupported);
    };
    let origin = triple(entities, reference(origin)?, "CARTESIAN_POINT", work)?.values;
    let mut numeric = Numeric::new(base, work)?;
    let scale = if matches!(scale, Value::Omitted) {
        I::point(1.)
    } else {
        let scalar = numeric.scalar(scale)?;
        work.step(2 * scalar.spelling().len() + 64)?;
        let [lo, hi] = scalar.bounds();
        finite(I::new(lo, hi))?
    };
    numeric.finish()?;
    if scale.lo <= 0. {
        return Err(Error::Unsupported);
    }
    let z = direction(entities, axis3, unit(2), work)?;
    let default_x = if z.zero[1] && z.zero[2] {
        unit(1)
    } else {
        unit(0)
    };
    let z = normalize(z.values, work)?;
    let x = normalize(direction(entities, axis1, default_x, work)?.values, work)?;
    let x = normalize(remove(x, z), work)?;
    let y = normalize(direction(entities, axis2, unit(1), work)?.values, work)?;
    // second_proj_axis uses the ORIGINAL v for both scalar projections.
    let z_projection = dot(y, z);
    let x_projection = dot(y, x);
    let y = std::array::from_fn(|i| y[i] - z_projection * z[i] - x_projection * x[i]);
    let y = normalize(y, work)?;
    work.step(1)?;
    Ok(MappingValues {
        origin,
        axes: [x, y, z],
        scale,
    })
}
impl<'input> Parsed<'input> {
    /// Read the actual complex operator's original numeric leaves. A caller
    /// must separately prove that this entity maps its represented member.
    pub fn physical_mapping<'p>(
        &'p self,
        operator: u32,
        options: Options,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<PhysicalMapping<'p, 'input>> {
        if options.maximum_work == 0
            || options.maximum_work > 16_000_000
            || options.retained_budget_bytes == 0
            || options.retained_budget_bytes > 64 * 1024 * 1024
        {
            return Err(Error::Invalid);
        }
        let mut work = Ledger {
            used: 0,
            remaining: options.maximum_work,
            interrupt,
        };
        work.step(self.work_used())?;
        work.step(512)?;
        let bytes = self.retained_budget_bytes() + SCRATCH + size_of::<PhysicalMapping<'_, '_>>();
        if bytes > options.retained_budget_bytes {
            return Err(Error::Limit);
        }
        let MappingValues {
            origin,
            axes,
            scale,
        } = evaluate_operator(&self.entities, operator, &mut work)?;
        Ok(PhysicalMapping {
            parsed: self,
            operator,
            origin,
            axes,
            scale,
            work: work.used,
            bytes,
        })
    }
}
#[cfg(test)]
mod tests;
