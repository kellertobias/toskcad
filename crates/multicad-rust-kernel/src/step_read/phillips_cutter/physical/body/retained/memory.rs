//! Conservative retained-source budget, not exact allocator capacity.
use super::*;
use std::mem::size_of;
fn allocation<T>(values: &Vec<T>) -> usize {
    values.capacity().saturating_mul(size_of::<T>())
}
fn values(values: &Vec<Value<'_>>, work: &mut Work<'_>) -> Result<usize> {
    let mut bytes = allocation(values);
    for value in values {
        work.step(1)?;
        let nested = match value {
            Value::List(items) => self::values(items, work)?,
            Value::Typed(entity) => self::values(&entity.args, work)?,
            _ => 0,
        };
        bytes = bytes.saturating_add(nested);
    }
    Ok(bytes)
}
pub(in super::super) fn budget(body: &Body<'_>, work: &mut Work<'_>) -> Result<usize> {
    // Audited pinned rustc58cf8f9c3 (Rust1.99.1): B=6,11 KV slots,
    //12 child pointers; insertion splits leave at least5 entries per nonroot.
    // https://github.com/rust-lang/rust/blob/58cf8f9c35b355199727627d07ae34a4e2846d67/library/alloc/src/collections/btree/node.rs
    // This parsed map is insertion-only: at most1+N/5 nodes. Eight pointer words
    // conservatively cover parent/index/length and all alignment padding; using
    // internal-node size for every leaf is an intentional additional reserve.
    let node = 20 * size_of::<usize>() + 11 * (size_of::<u32>() + size_of::<Vec<Entity<'_>>>());
    let nodes = 1_usize.saturating_add(body.entities.len() / 5);
    let mut bytes = size_of::<Body<'_>>()
        .saturating_add(body.input.len())
        .saturating_add(nodes.saturating_mul(node));
    for records in body.entities.values() {
        work.step(1)?;
        bytes = bytes.saturating_add(allocation(records));
        for entity in records {
            bytes = bytes.saturating_add(values(&entity.args, work)?);
        }
    }
    for list in [
        &body.material_faces,
        &body.profile_faces,
        &body.cap.removed_cavity_faces,
    ] {
        bytes = bytes.saturating_add(allocation(list));
    }
    let b = &body.cavity;
    bytes = bytes.saturating_add(b.exact_support_keys_capacity_bytes());
    bytes = bytes
        .saturating_add(allocation(&b.vertices))
        .saturating_add(allocation(&b.edges))
        .saturating_add(allocation(&b.faces));
    bytes = bytes
        .saturating_add(allocation(&b.mouth.material_uses))
        .saturating_add(allocation(&b.mouth.material_loop.uses));
    for face in &b.faces {
        work.step(1)?;
        bytes = bytes.saturating_add(allocation(&face.uses));
        if let Some(loop_) = &face.material_loop {
            bytes = bytes.saturating_add(allocation(&loop_.uses));
        }
    }
    for edge in &b.edges {
        work.step(1)?;
        if let Curve::CubicIntersection {
            pole_entities,
            spans,
            ..
        } = &edge.curve
        {
            bytes = bytes
                .saturating_add(allocation(pole_entities))
                .saturating_add(allocation(spans));
        }
    }
    Ok(bytes)
}
