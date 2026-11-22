//! Cyclic local edge fans of the already closed directed incidence.
use super::*;
pub(super) fn validate(
    faces: &[Face],
    edges: &BTreeMap<u32, Edge>,
    vertices: usize,
    b: &mut Budget<'_>,
) -> Result<(), Error> {
    let corners = edges.len().checked_mul(2).ok_or(Error::Limit)?;
    let reserve = (corners + vertices)
        .checked_mul(512)
        .and_then(|v| v.checked_add(512))
        .ok_or(Error::Limit)?;
    b.retain(reserve)?;
    let tree_cost = 64 * (1 + (edges.len() + 1).ilog2() as usize);
    let mut fans: BTreeMap<u32, BTreeMap<u32, u32>> = BTreeMap::new();
    for face in faces {
        for ring in &face.loops {
            for i in 0..ring.len() {
                b.charge(tree_cost)?;
                let incoming = &ring[i];
                let outgoing = &ring[(i + 1) % ring.len()];
                let edge = edges.get(&incoming.edge).ok_or(Error::Invalid)?;
                let vertex = edge.vertices[usize::from(incoming.forward)];
                if incoming.edge == outgoing.edge {
                    return Err(Error::Invalid);
                }
                if fans
                    .entry(vertex)
                    .or_default()
                    .insert(incoming.edge, outgoing.edge)
                    .is_some()
                {
                    return Err(Error::Invalid);
                }
            }
        }
    }
    if fans.len() != vertices {
        return Err(Error::Invalid);
    }
    for fan in fans.values() {
        let first = *fan.keys().next().ok_or(Error::Invalid)?;
        let mut current = first;
        let mut visited = BTreeSet::new();
        loop {
            b.charge(tree_cost)?;
            if !visited.insert(current) {
                break;
            }
            current = *fan.get(&current).ok_or(Error::Invalid)?;
        }
        // Exactly one successor and predecessor per edge end follows from the
        // opposite two-use check; visiting one cycle must consume the complete fan.
        if current != first || visited.len() != fan.len() || fan.len() < 2 {
            return Err(Error::Invalid);
        }
    }
    b.bytes -= reserve;
    Ok(())
}
