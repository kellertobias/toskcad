//! Transfer an ACTUAL independently earned local material proof only after
//! neutral OUTPUT parse and complete component/ID/role identity. No writer
//! trace, author recipe or buffer-pointer equality supplies admission.
use super::*;
#[cfg(test)]
mod tests;
fn bind(
    doc: &CompoundDocument<'_>,
    envelope: &schema::Envelope,
    local: &OwnedActualHead,
    b: &mut Budget<'_>,
) -> Result<()> {
    let root = encoding::root(local, b)?;
    if envelope.top != root.top
        || envelope.axis != root.axis
        || envelope.context != root.context
        || envelope.solid != local.physical_roles().solid_entity
    {
        return Err(Error::Unsupported);
    }
    let records = doc.records();
    for i in 0..local.original_record_count() {
        b.step(32)?;
        let id = local.original_entity_at(i).ok_or(Error::Invalid)?;
        let parts = records.get(&id).ok_or(Error::Unsupported)?;
        if id == root.top {
            // The shared envelope schema already earned the exact rewritten
            // membership, identity datum and context. Only its descriptive
            // name and mapped-item replacement differ from local metadata.
            continue;
        }
        let mut part = 0;
        while let Some(expected) = local.original_component_at(i, part) {
            let actual = parts
                .get(part)
                .and_then(|p| p.source)
                .ok_or(Error::Unsupported)?;
            b.step(
                expected
                    .len()
                    .checked_add(actual.len())
                    .and_then(|n| n.checked_add(32))
                    .ok_or(Error::Limit)?,
            )?;
            if expected != actual {
                return Err(Error::Unsupported);
            }
            part += 1;
        }
        if part != parts.len() {
            return Err(Error::Unsupported);
        }
    }
    // NO additional geometry/topology/UV/metadata may be discarded. Only the
    // explicitly typed operator/map/source-datum closure may be new. Every
    // original owner record was already matched above, including scalar refs.
    for &id in records.keys() {
        b.step(64)?;
        if local.original_component(id, 0).is_none()
            && !envelope.removable[..envelope.count].contains(&id)
        {
            return Err(Error::Unsupported);
        }
    }
    Ok(())
}
impl MappedActualHead {
    /// Source-known output path. The output is independently neutrally parsed;
    /// SAME local geometry authority transfers by literal full-record identity.
    /// A mismatch refuses; use the ordinary cold reader for independent import.
    /// `additional_resident_bytes` excludes these two explicit input owners.
    pub fn read_reusing_local(
        input: Arc<str>,
        local: Arc<OwnedActualHead>,
        options: OwnedHeadOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self> {
        if !options.material_tolerance_mm.is_finite()
            || options.material_tolerance_mm <= 0.
            || options.material_tolerance_mm > 1e-4
            || !options.minimum_wall_mm.is_finite()
            || options.minimum_wall_mm <= 0.
        {
            return Err(Error::Invalid);
        }
        if local.material_tolerance_mm() > options.material_tolerance_mm
            || local.minimum_wall_mm() < options.minimum_wall_mm
        {
            return Err(Error::Unsupported);
        }
        let r = options.resources;
        let lr = local.receipt();
        let Prepared {
            doc,
            mut b,
            envelope,
            frame,
            cold_work,
            cold_peak,
            outside,
            header,
        } = prepare_envelope(&input, r, lr.aggregate_work, lr.retained_bytes, interrupt)?;
        bind(&doc, &envelope, &local, &mut b)?;
        b.step(256)?;
        let world_bounds = frame.map_box_mm(local.local_sewn_enclosure_mm())?;
        let s = I::new(frame.scale[0], frame.scale[1]);
        let raw = local.local_sewn_volume_mm3();
        let volume = I::new(raw[0], raw[1]) * s.square() * s;
        if !volume.valid() || volume.lo <= 0. {
            return Err(Error::Unsupported);
        }
        let inherited_peak = lr
            .accounted_peak_bytes
            .checked_sub(lr.additional_resident_bytes)
            .and_then(|n| n.checked_add(r.additional_resident_bytes))
            .and_then(|n| n.checked_add(header))
            .and_then(|n| n.checked_add(input.len()))
            .ok_or(Error::Limit)?;
        let peak = b
            .peak
            .max(inherited_peak)
            .max(outside.checked_add(cold_peak).ok_or(Error::Limit)?);
        let retained = header
            .checked_add(input.len())
            .and_then(|n| n.checked_add(lr.retained_bytes))
            .ok_or(Error::Limit)?;
        let aggregate = lr
            .aggregate_work
            .checked_add(cold_work)
            .and_then(|n| n.checked_add(b.used))
            .ok_or(Error::Limit)?;
        if peak > r.maximum_bytes || aggregate > r.maximum_work || retained > r.maximum_bytes {
            return Err(Error::Limit);
        }
        b.step(0)?;
        Ok(Self {
            original: input,
            local,
            frame,
            member: envelope.item,
            operator: envelope.operator,
            world_bounds,
            world_volume: [volume.lo, volume.hi],
            receipt: MappedHeadReceipt {
                original_cold_work_reserved: cold_work,
                exact_local_material_reused: true,
                envelope_work: b.used,
                local: lr,
                aggregate_work: aggregate,
                accounted_peak_bytes: peak,
                retained_bytes: retained,
                additional_resident_bytes: r.additional_resident_bytes,
            },
        })
    }
}
