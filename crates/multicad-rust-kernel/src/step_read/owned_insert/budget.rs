//! Retained source accounting, not a measurement of allocator/parser peak.
use super::*;
fn allocation<T>(v: &Vec<T>) -> usize {
    v.capacity().saturating_mul(size_of::<T>())
}
fn values(v: &Vec<Value<'_>>, w: &mut Ledger<'_>) -> Result<usize> {
    let mut bytes = allocation(v);
    for item in v {
        w.step(1)?;
        let n = match item {
            Value::List(v) => values(v, w)?,
            Value::Typed(e) => values(&e.args, w)?,
            _ => 0,
        };
        bytes = bytes.saturating_add(n);
    }
    Ok(bytes)
}
pub(super) fn arena(
    input: &str,
    entities: &BTreeMap<u32, Vec<Entity<'_>>>,
    w: &mut Ledger<'_>,
) -> Result<usize> {
    // Same audited pinned Rust1.99.1 insertion-only BTreeMap bound as the owned
    // Phillips raw-body packet: B=6,11KV slots,12children, >=5keys/nonroot.
    // rust-lang/rust58cf8f9c35b355199727627d07ae34a4e2846d67
    let node = 20 * size_of::<usize>() + 11 * (size_of::<u32>() + size_of::<Vec<Entity<'_>>>());
    let mut bytes = input
        .len()
        .saturating_add((1 + entities.len() / 5).saturating_mul(node));
    for records in entities.values() {
        w.step(1)?;
        bytes = bytes.saturating_add(allocation(records));
        for e in records {
            bytes = bytes.saturating_add(values(&e.args, w)?);
        }
    }
    Ok(bytes)
}
pub(super) fn packet(p: &Parsed<'_>, w: &mut Ledger<'_>) -> Result<usize> {
    let mut b = arena(p.input, &p.entities, w)?.saturating_add(size_of::<Parsed<'_>>());
    b = b.saturating_add(allocation(&p.points));
    for point in &p.points {
        w.step(1)?;
        b = b.saturating_add(allocation(&point.coordinates));
    }
    b = b.saturating_add(allocation(&p.curves));
    for curve in &p.curves {
        w.step(1)?;
        b = b
            .saturating_add(allocation(&curve.poles))
            .saturating_add(allocation(&curve.multiplicities))
            .saturating_add(allocation(&curve.knots));
    }
    b = b.saturating_add(allocation(&p.tensors));
    for s in &p.tensors {
        w.step(1)?;
        b = b
            .saturating_add(allocation(&s.rows[0]))
            .saturating_add(allocation(&s.rows[1]))
            .saturating_add(allocation(&s.multiplicities_v))
            .saturating_add(allocation(&s.knots_v));
    }
    b = b
        .saturating_add(allocation(&p.planes))
        .saturating_add(allocation(&p.vertices))
        .saturating_add(allocation(&p.edges))
        .saturating_add(allocation(&p.faces));
    for f in &p.faces {
        w.step(1)?;
        b = b.saturating_add(allocation(&f.loops));
        for l in &f.loops {
            b = b.saturating_add(allocation(&l.coedges));
        }
    }
    Ok(b)
}
