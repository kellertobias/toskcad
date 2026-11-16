//! Actual one-member mapped envelope over independently earned physical head.
//! The local wrapper is DERIVED; all physical component IDs/scalar slices are
//! copied unchanged. This proves the FILE operator, never an author frame.
use super::*;
use crate::tapered_thread::interval::I;
use std::{
    fmt::{self, Write},
    sync::Arc,
};
mod schema;
#[cfg(test)]
mod tests;
const SCRATCH: usize = 256 * 1024;
const HEADER: &str = concat!(
    "ISO-10303-21;\nHEADER;\n",
    "FILE_DESCRIPTION(('DERIVED identity view of actual mapped clipped head'),'2;1');\n",
    "FILE_NAME('DERIVED','1970-01-01T00:00:00',(''),(''),'MultiCAD','MultiCAD','');\n",
    "FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\nENDSEC;\nDATA;\n"
);
const TAIL: &str = "ENDSEC;\nEND-ISO-10303-21;\n";
#[derive(Clone, Copy, Debug)]
pub struct MappedHeadFrame {
    origin: [[f64; 2]; 3],
    axes: [[[f64; 2]; 3]; 3],
    scale: [f64; 2],
}
impl MappedHeadFrame {
    pub fn origin_mm(self) -> [[f64; 2]; 3] {
        self.origin
    }
    /// Actual normalized EXPRESS basis columns, before positive uniform scale.
    pub fn axes(self) -> [[[f64; 2]; 3]; 3] {
        self.axes
    }
    pub fn scale_bounds(self) -> [f64; 2] {
        self.scale
    }
    pub fn map_box_mm(self, input: [[f64; 2]; 3]) -> Result<[[f64; 2]; 3]> {
        let input = input.map(|v| I::new(v[0], v[1]));
        if input.iter().any(|v| !v.valid()) {
            return Err(Error::Invalid);
        }
        let s = I::new(self.scale[0], self.scale[1]);
        let result = std::array::from_fn::<_, 3, _>(|i| {
            (0..3).fold(I::new(self.origin[i][0], self.origin[i][1]), |v, j| {
                v + s * I::new(self.axes[j][i][0], self.axes[j][i][1]) * input[j]
            })
        });
        if result.iter().any(|v| !v.valid()) {
            return Err(Error::Unsupported);
        }
        Ok(result.map(|v| [v.lo, v.hi]))
    }
}
#[derive(Clone, Copy, Debug)]
pub struct MappedHeadReceipt {
    pub original_cold_work_reserved: usize,
    pub envelope_work: usize,
    pub local: OwnedHeadReceipt,
    pub aggregate_work: usize,
    pub accounted_peak_bytes: usize,
    pub retained_bytes: usize,
}
pub struct MappedActualHead {
    original: Arc<str>,
    local: OwnedActualHead,
    frame: MappedHeadFrame,
    member: u32,
    operator: u32,
    receipt: MappedHeadReceipt,
    world_bounds: [[f64; 2]; 3],
    world_volume: [f64; 2],
}
impl MappedActualHead {
    pub fn actual_file(&self) -> &str {
        &self.original
    }
    pub fn local_geometry(&self) -> &OwnedActualHead {
        &self.local
    }
    pub fn actual_frame(&self) -> MappedHeadFrame {
        self.frame
    }
    pub fn member_entity(&self) -> u32 {
        self.member
    }
    pub fn operator_entity(&self) -> u32 {
        self.operator
    }
    pub fn receipt(&self) -> MappedHeadReceipt {
        self.receipt
    }
    pub fn world_sewn_enclosure_mm(&self) -> [[f64; 2]; 3] {
        self.world_bounds
    }
    /// SAME declared sewn local reference, under the actual proper file map.
    /// Raw physical mass and author-law correspondence remain separate.
    pub fn world_sewn_volume_mm3(&self) -> [f64; 2] {
        self.world_volume
    }
    pub fn read(
        input: Arc<str>,
        options: OwnedHeadOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self> {
        let r = options.resources;
        if r.maximum_work == 0
            || r.maximum_work > MAX_WORK
            || r.maximum_bytes == 0
            || r.maximum_bytes > MAX_BYTES
        {
            return Err(Error::Invalid);
        }
        let header = size_of::<Self>() + 4 * size_of::<usize>();
        let outside = r
            .additional_resident_bytes
            .checked_add(header)
            .ok_or(Error::Limit)?;
        let available = r
            .maximum_bytes
            .checked_sub(outside)
            .filter(|n| *n > 0)
            .ok_or(Error::Limit)?;
        let cold = cold::reserve(&input, r.maximum_work.min(64_000_000), available, interrupt)
            .map_err(geometry_error)?;
        let doc = CompoundDocument::parse(
            &input,
            DocumentOptions {
                maximum_work: cold.delegated_parser_work_reserved,
                retained_budget_bytes: available,
            },
            interrupt,
        )
        .map_err(neutral_error)?;
        let initial = outside
            .checked_add(doc.retained_budget_bytes())
            .and_then(|n| n.checked_add(SCRATCH))
            .ok_or(Error::Limit)?;
        if initial > r.maximum_bytes {
            return Err(Error::Limit);
        }
        let mut b = Budget {
            used: 0,
            remaining: r
                .maximum_work
                .checked_sub(cold.aggregate_work_reserved)
                .ok_or(Error::Limit)?
                .min(64_000_000),
            bytes: initial,
            peak: initial,
            limit: r.maximum_bytes,
            interrupt,
        };
        let envelope = schema::read(&doc, &mut b)?;
        let remaining = b.remaining.min(16_000_000);
        let map = crate::step_read::owned_insert::mapping::evaluate_neutral_operator(
            doc.records(),
            envelope.operator,
            crate::step_read::owned_insert::Options {
                maximum_work: remaining,
                retained_budget_bytes: r
                    .maximum_bytes
                    .checked_sub(b.bytes)
                    .filter(|n| *n > 0)
                    .ok_or(Error::Limit)?,
            },
            interrupt,
        )
        .map_err(neutral_error)?;
        b.step(map.work)?;
        b.peak = b
            .peak
            .max(b.bytes.checked_add(map.bytes).ok_or(Error::Limit)?);
        if b.peak > r.maximum_bytes {
            return Err(Error::Limit);
        }
        if map.handedness != 1 {
            return Err(Error::Unsupported);
        }
        let frame = MappedHeadFrame {
            origin: map.origin,
            axes: map.axes,
            scale: map.scale,
        };
        let text = derive(&doc, &envelope, &mut b)?;
        // Arc slice copy overlaps the original output allocation. Reserve the
        // exact payload and Arc header BEFORE allocating that immutable copy.
        b.step(text.len())?;
        let conversion = b
            .bytes
            .checked_add(text.len())
            .and_then(|n| n.checked_add(2 * size_of::<usize>()))
            .ok_or(Error::Limit)?;
        if conversion > r.maximum_bytes {
            return Err(Error::Limit);
        }
        b.peak = b.peak.max(conversion);
        let local_text: Arc<str> = Arc::from(text.as_str());
        let old_peak = b.peak.max(
            outside
                .checked_add(cold.parser_peak_bytes_reserved)
                .ok_or(Error::Limit)?,
        );
        // Reserve actual fixed interval world-box and scaled-volume work
        // before the child begins. Their numerical guards run before publish.
        b.step(256)?;
        let fresh = b.used;
        let prefix = cold
            .aggregate_work_reserved
            .checked_add(fresh)
            .ok_or(Error::Limit)?;
        let member = envelope.item;
        let operator = envelope.operator;
        drop(text);
        drop(doc);
        let child_outside = outside.checked_add(input.len()).ok_or(Error::Limit)?;
        let local = OwnedActualHead::read(
            local_text,
            OwnedHeadOptions {
                resources: Options {
                    maximum_work: r
                        .maximum_work
                        .checked_sub(prefix)
                        .filter(|n| *n > 0)
                        .ok_or(Error::Limit)?,
                    additional_resident_bytes: child_outside,
                    ..r
                },
                ..options
            },
            interrupt,
        )?;
        let world_bounds = frame.map_box_mm(local.local_sewn_enclosure_mm())?;
        let s = I::new(frame.scale[0], frame.scale[1]);
        let raw_volume = local.local_sewn_volume_mm3();
        let volume = I::new(raw_volume[0], raw_volume[1]) * s.square() * s;
        if !volume.valid() || volume.lo <= 0. {
            return Err(Error::Unsupported);
        }
        let world_volume = [volume.lo, volume.hi];
        let lr = local.receipt();
        let aggregate = prefix.checked_add(lr.aggregate_work).ok_or(Error::Limit)?;
        let peak = old_peak.max(lr.accounted_peak_bytes);
        let retained = header
            .checked_add(input.len())
            .and_then(|n| n.checked_add(lr.retained_bytes))
            .ok_or(Error::Limit)?;
        if peak > r.maximum_bytes || aggregate > r.maximum_work || retained > r.maximum_bytes {
            return Err(Error::Limit);
        }
        if interrupt.is_some_and(|f| f()) {
            return Err(Error::Cancelled);
        }
        Ok(Self {
            original: input,
            local,
            frame,
            member,
            operator,
            world_bounds,
            world_volume,
            receipt: MappedHeadReceipt {
                original_cold_work_reserved: cold.aggregate_work_reserved,
                envelope_work: fresh,
                local: lr,
                aggregate_work: aggregate,
                accounted_peak_bytes: peak,
                retained_bytes: retained,
            },
        })
    }
}
struct Count(usize);
impl Write for Count {
    fn write_str(&mut self, s: &str) -> fmt::Result {
        self.0 = self.0.checked_add(s.len()).ok_or(fmt::Error)?;
        Ok(())
    }
}
fn references(
    values: &[Value<'_>],
    ids: &[u32],
    selected: &mut [bool],
    pending: &mut Vec<usize>,
    b: &mut Budget<'_>,
) -> Result<()> {
    for v in values {
        b.step(32)?;
        match v {
            Value::Reference(id) => {
                let i = ids.binary_search(id).map_err(|_| Error::Invalid)?;
                if !selected[i] {
                    selected[i] = true;
                    pending.push(i);
                }
            }
            Value::List(v) => references(v, ids, selected, pending, b)?,
            Value::Typed(e) => references(&e.args, ids, selected, pending, b)?,
            _ => {}
        }
    }
    Ok(())
}
fn derive(doc: &CompoundDocument<'_>, e: &schema::Envelope, b: &mut Budget<'_>) -> Result<String> {
    let records = doc.records();
    let mut ids = b.allocate(records.len())?;
    ids.extend(records.keys().copied());
    let mut selected = b.allocate(records.len())?;
    selected.resize(records.len(), false);
    let mut pending = b.allocate(records.len())?;
    // Start from actual metadata roots. The original top item's reference is
    // replaced only after actual map membership was proved. All other original
    // references remain literal; typed local metadata rejects smuggled geometry.
    for (i, parts) in records.values().enumerate() {
        b.step(32)?;
        if matches!(parts.as_slice(),[v] if matches!(v.name,"SHAPE_DEFINITION_REPRESENTATION"|"APPLICATION_PROTOCOL_DEFINITION"|"PRODUCT_RELATED_PRODUCT_CATEGORY"))
        {
            selected[i] = true;
            pending.push(i);
        }
    }
    while let Some(i) = pending.pop() {
        b.step(32)?;
        if ids[i] == e.top {
            for id in [e.solid, e.axis, e.context] {
                b.step(32)?;
                let k = ids.binary_search(&id).map_err(|_| Error::Invalid)?;
                if !selected[k] {
                    selected[k] = true;
                    pending.push(k);
                }
            }
        } else {
            for part in &records[&ids[i]] {
                references(&part.args, &ids, &mut selected, &mut pending, b)?;
            }
        }
    }
    for (i, id) in ids.iter().enumerate() {
        b.step(16)?;
        if !selected[i] && !e.removable[..e.count].contains(id) {
            return Err(Error::Unsupported);
        }
    }
    fn emit(
        out: &mut impl Write,
        records: &BTreeMap<u32, Vec<Entity<'_>>>,
        ids: &[u32],
        selected: &[bool],
        e: &schema::Envelope,
    ) -> fmt::Result {
        out.write_str(HEADER)?;
        for (i, id) in ids.iter().enumerate() {
            if !selected[i] {
                continue;
            }
            if *id == e.top {
                writeln!(
                    out,
                    "#{id}=ADVANCED_BREP_SHAPE_REPRESENTATION('',(#{},#{}),#{});",
                    e.solid, e.axis, e.context
                )?;
            } else {
                let parts = &records[id];
                write!(out, "#{id}=")?;
                if parts.len() > 1 {
                    out.write_char('(')?;
                }
                for p in parts {
                    out.write_str(p.source.ok_or(fmt::Error)?)?;
                }
                if parts.len() > 1 {
                    out.write_char(')')?;
                }
                out.write_str(";\n")?;
            }
        }
        out.write_str(TAIL)
    }
    let mut count = Count(0);
    emit(&mut count, records, &ids, &selected, e).map_err(|_| Error::Limit)?;
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
    emit(&mut text, records, &ids, &selected, e).map_err(|_| Error::Limit)?;
    b.release(&ids);
    b.release(&selected);
    b.release(&pending);
    Ok(text)
}
