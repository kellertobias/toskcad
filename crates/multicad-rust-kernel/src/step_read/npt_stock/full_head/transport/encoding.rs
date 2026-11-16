//! Candidate only: replay actual physical records and add a new mapped
//! envelope. The independent reader/source-frame gate must earn admission.
use super::*;
use crate::{
    step::pipe_hardware::numbers,
    step_read::{
        literal_member_binding::{Token, token},
        npt_stock::LiteralOccurrencePlacement,
    },
};
#[cfg(test)]
mod tests;
const MAPPED_HEADER: &str = concat!(
    "ISO-10303-21;\nHEADER;\n",
    "FILE_DESCRIPTION(('MultiCAD actual clipped head mapped candidate'),'2;1');\n",
    "FILE_NAME('MultiCAD','1970-01-01T00:00:00',(''),(''),'MultiCAD','MultiCAD','');\n",
    "FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\nENDSEC;\nDATA;\n"
);
const EDIT_BYTES: usize = 32 * 1024;
#[derive(Clone, Copy, Debug)]
pub struct EncodingOptions {
    pub maximum_work: usize,
    pub maximum_bytes: usize,
    pub additional_resident_bytes: usize,
}
impl Default for EncodingOptions {
    fn default() -> Self {
        Self {
            maximum_work: MAX_WORK,
            maximum_bytes: MAX_BYTES,
            additional_resident_bytes: 0,
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct EncodingReceipt {
    pub inherited_material_work: usize,
    pub inherited_placement_work: usize,
    pub fresh_encoding_work: usize,
    pub aggregate_work: usize,
    pub accounted_peak_bytes: usize,
}
pub struct MappedCandidate {
    text: String,
    receipt: EncodingReceipt,
}
impl MappedCandidate {
    pub fn as_str(&self) -> &str {
        &self.text
    }
    pub fn receipt(&self) -> EncodingReceipt {
        self.receipt
    }
}
struct Root {
    top: u32,
    axis: u32,
    context: u32,
    maximum: u32,
}
fn root(local: &OwnedActualHead, b: &mut Budget<'_>) -> Result<Root> {
    let mut top = None;
    let mut maximum = 0;
    for i in 0..local.original_record_count() {
        b.step(32)?;
        let id = local.original_entity_at(i).ok_or(Error::Invalid)?;
        maximum = maximum.max(id);
        let component = local.original_component_at(i, 0).ok_or(Error::Invalid)?;
        if component
            .trim_start()
            .starts_with("ADVANCED_BREP_SHAPE_REPRESENTATION")
        {
            if top.is_some() {
                return Err(Error::Unsupported);
            }
            b.step(component.len() + 128)?;
            let mut rest = component;
            let mut refs = [0; 3];
            let mut n = 0;
            while let Some(t) = token(&mut rest) {
                b.step(1)?;
                if let Token::Reference(id) = t {
                    *refs.get_mut(n).ok_or(Error::Unsupported)? = id;
                    n += 1;
                }
            }
            if n != 3 {
                return Err(Error::Unsupported);
            }
            let solid = local.physical_roles().solid_entity;
            let axis = if refs[0] == solid {
                refs[1]
            } else if refs[1] == solid {
                refs[0]
            } else {
                return Err(Error::Unsupported);
            };
            top = Some((id, axis, refs[2]));
        }
    }
    let (top, axis, context) = top.ok_or(Error::Unsupported)?;
    Ok(Root {
        top,
        axis,
        context,
        maximum,
    })
}
fn emit(
    out: &mut impl Write,
    local: &OwnedActualHead,
    root: &Root,
    ids: [u32; 8],
    strings: &[[String; 3]; 4],
    scale: &str,
) -> fmt::Result {
    let [rep, map, item, op, origin, x, y, z] = ids;
    out.write_str(MAPPED_HEADER)?;
    for i in 0..local.original_record_count() {
        let id = local.original_entity_at(i).ok_or(fmt::Error)?;
        if id == root.top {
            writeln!(
                out,
                "#{id}=ADVANCED_BREP_SHAPE_REPRESENTATION('',(#{item},#{}),#{});",
                root.axis, root.context
            )?;
        } else {
            write!(out, "#{id}=")?;
            let complex = local.original_component_at(i, 1).is_some();
            if complex {
                out.write_char('(')?;
            }
            let mut part = 0;
            while let Some(s) = local.original_component_at(i, part) {
                out.write_str(s)?;
                part += 1;
            }
            if complex {
                out.write_char(')')?;
            }
            out.write_str(";\n")?;
        }
    }
    writeln!(
        out,
        "#{rep}=SHAPE_REPRESENTATION('',(#{},#{}),#{});",
        local.physical_roles().solid_entity,
        root.axis,
        root.context
    )?;
    writeln!(out, "#{map}=REPRESENTATION_MAP(#{},#{rep});", root.axis)?;
    writeln!(out, "#{item}=MAPPED_ITEM('',#{map},#{op});")?;
    writeln!(
        out,
        "#{op}=(CARTESIAN_TRANSFORMATION_OPERATOR(#{x},#{y},#{origin},{scale})CARTESIAN_TRANSFORMATION_OPERATOR_3D(#{z})FUNCTIONALLY_DEFINED_TRANSFORMATION('','')GEOMETRIC_REPRESENTATION_ITEM()REPRESENTATION_ITEM(''));"
    )?;
    for (i, id) in [origin, x, y, z].into_iter().enumerate() {
        writeln!(
            out,
            "#{id}={}('',({},{},{}));",
            if i == 0 {
                "CARTESIAN_POINT"
            } else {
                "DIRECTION"
            },
            strings[i][0],
            strings[i][1],
            strings[i][2]
        )?;
    }
    out.write_str(TAIL)
}
/// Encode one actual local owner using a rounded normalized output proposal.
/// The original immutable factor chain remains the source-matching target;
/// encoding never equates it with STEP's normalized EXPRESS operator.
pub fn propose_mapped(
    local: &OwnedActualHead,
    placement: &LiteralOccurrencePlacement,
    options: EncodingOptions,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<MappedCandidate> {
    if options.maximum_work == 0
        || options.maximum_work > MAX_WORK
        || options.maximum_bytes == 0
        || options.maximum_bytes > MAX_BYTES
    {
        return Err(Error::Invalid);
    }
    let inherited = local
        .receipt()
        .aggregate_work
        .checked_add(placement.work_used())
        .ok_or(Error::Limit)?;
    let resident = local
        .receipt()
        .retained_bytes
        .checked_add(placement.retained_budget_bytes())
        .and_then(|n| n.checked_add(options.additional_resident_bytes))
        .and_then(|n| n.checked_add(size_of::<MappedCandidate>()))
        .ok_or(Error::Limit)?;
    let initial = resident.checked_add(EDIT_BYTES).ok_or(Error::Limit)?;
    let prior = local
        .receipt()
        .accounted_peak_bytes
        .checked_add(placement.retained_budget_bytes())
        .and_then(|n| n.checked_add(options.additional_resident_bytes))
        .and_then(|n| n.checked_add(size_of::<MappedCandidate>()))
        .ok_or(Error::Limit)?;
    let pose_peak = placement
        .accounted_peak_bytes()
        .checked_add(local.receipt().retained_bytes)
        .and_then(|n| n.checked_add(options.additional_resident_bytes))
        .and_then(|n| n.checked_add(size_of::<MappedCandidate>()))
        .ok_or(Error::Limit)?;
    if initial > options.maximum_bytes
        || prior > options.maximum_bytes
        || pose_peak > options.maximum_bytes
    {
        return Err(Error::Limit);
    }
    let mut b = Budget {
        used: 0,
        remaining: options
            .maximum_work
            .checked_sub(inherited)
            .ok_or(Error::Limit)?
            .min(16_000_000),
        bytes: initial,
        peak: initial.max(prior).max(pose_peak),
        limit: options.maximum_bytes,
        interrupt,
    };
    let root = root(local, &mut b)?;
    let ids = std::array::from_fn::<_, 8, _>(|i| root.maximum.checked_add(i as u32 + 1));
    let [
        Some(rep),
        Some(map),
        Some(item),
        Some(op),
        Some(origin),
        Some(x),
        Some(y),
        Some(z),
    ] = ids
    else {
        return Err(Error::Limit);
    };
    let ids = [rep, map, item, op, origin, x, y, z];
    b.step(1024)?;
    let rows = placement
        .linear_rows_bounds()
        .map(|r| r.map(|x| I::new(x[0], x[1])));
    let scale = (0..3)
        .fold(I::point(0.), |s, i| s + rows[i][0].square())
        .sqrt()
        .midpoint();
    if !scale.is_finite() || scale <= 0. {
        return Err(Error::Unsupported);
    }
    let origin = placement
        .translation_bounds_mm()
        .map(|x| I::new(x[0], x[1]).midpoint());
    let values = std::array::from_fn::<_, 4, _>(|i| {
        if i == 0 {
            origin
        } else {
            std::array::from_fn(|j| rows[j][i - 1].midpoint() / scale)
        }
    });
    let mut strings: [[String; 3]; 4] =
        std::array::from_fn(|_| std::array::from_fn(|_| String::new()));
    for (row, values) in strings.iter_mut().zip(values) {
        for (text, v) in row.iter_mut().zip(values) {
            let (work, bytes) = numbers::cost(v).map_err(geometry_error)?;
            b.step(work.checked_mul(64).ok_or(Error::Limit)?)?;
            if bytes > 4096 {
                return Err(Error::Limit);
            }
            *text = numbers::literal(v).map_err(geometry_error)?;
        }
    }
    let (work, bytes) = numbers::cost(scale).map_err(geometry_error)?;
    b.step(work.checked_mul(64).ok_or(Error::Limit)?)?;
    if bytes > 4096 {
        return Err(Error::Limit);
    }
    let scale = numbers::literal(scale).map_err(geometry_error)?;
    // At most13 finite binary64 leaves, each≤~1100 decimal bytes. The32KiB
    // edit arena includes all String capacities and bigint decimal temporaries.
    // Dense span access is constant-time; count/emit each visit at most four
    // components per admitted record, plus bounded ID formatting. Charge the
    // two structural passes BEFORE either traversal; payload copy is separate.
    b.step(
        local
            .original_record_count()
            .checked_mul(32)
            .ok_or(Error::Limit)?,
    )?;
    let mut count = Count(0);
    emit(&mut count, local, &root, ids, &strings, &scale).map_err(|_| Error::Limit)?;
    b.step(count.0.checked_mul(2).ok_or(Error::Limit)?)?;
    if count.0 > b.limit - b.bytes {
        return Err(Error::Limit);
    }
    let mut text = String::new();
    text.try_reserve_exact(count.0).map_err(|_| Error::Limit)?;
    b.bytes = b.bytes.checked_add(text.capacity()).ok_or(Error::Limit)?;
    b.peak = b.peak.max(b.bytes);
    if b.bytes > b.limit {
        return Err(Error::Limit);
    }
    emit(&mut text, local, &root, ids, &strings, &scale).map_err(|_| Error::Limit)?;
    b.step(0)?;
    let fresh = b.used;
    Ok(MappedCandidate {
        text,
        receipt: EncodingReceipt {
            inherited_material_work: local.receipt().aggregate_work,
            inherited_placement_work: placement.work_used(),
            fresh_encoding_work: fresh,
            aggregate_work: inherited.checked_add(fresh).ok_or(Error::Limit)?,
            accounted_peak_bytes: b.peak,
        },
    })
}
