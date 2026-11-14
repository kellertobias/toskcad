//! Actual indexed closure, vertex-link cycles and directed literal volume.
use super::*;
#[derive(Clone, Copy)]
struct Edge {
    key: u64,
    face: usize,
    slot: usize,
    count: u8,
    orientation: i8,
}
pub(super) fn certify(mesh: &DisplayMesh, m: &mut Memory, w: &mut Work<'_>) -> R<()> {
    let vertices = mesh.vertices.len() / 3;
    let triangles = mesh.indices.as_chunks::<3>().0;
    if mesh.vertices.len() % 3 != 0
        || mesh.indices.len() % 3 != 0
        || triangles.len() != mesh.face_indices.len()
        || triangles.is_empty()
    {
        return Err(E::Unresolved);
    }
    let buckets = triangles
        .len()
        .checked_mul(6)
        .and_then(usize::checked_next_power_of_two)
        .ok_or(E::Limit)?;
    let mut edges: Vec<Option<Edge>> = m.allocate(buckets, w)?;
    w.charge(buckets)?;
    edges.resize(buckets, None);
    let mut adjacent = m.allocate(triangles.len(), w)?;
    adjacent.resize(triangles.len(), [usize::MAX; 3]);
    let mut first = m.allocate(vertices, w)?;
    first.resize(vertices, usize::MAX);
    let mut degree = m.allocate(vertices, w)?;
    degree.resize(vertices, 0usize);
    let mut unique = 0usize;
    for (f, tri) in triangles.iter().enumerate() {
        for slot in 0..3 {
            w.step()?;
            let a = tri[slot] as usize;
            let b = tri[(slot + 1) % 3] as usize;
            if a == b || a >= vertices || b >= vertices {
                return Err(E::Unresolved);
            }
            if first[a] == usize::MAX {
                first[a] = f;
            }
            degree[a] += 1;
            let key = ((a.min(b) as u64) << 32) | a.max(b) as u64;
            let orientation = if a < b { 1 } else { -1 };
            let mut index = (key.wrapping_mul(0x9e3779b97f4a7c15) as usize) & (buckets - 1);
            loop {
                w.step()?;
                if let Some(old) = &mut edges[index] {
                    if old.key == key {
                        if old.count != 1 || old.orientation == orientation {
                            return Err(E::Unresolved);
                        }
                        adjacent[f][slot] = old.face;
                        adjacent[old.face][old.slot] = f;
                        old.count = 2;
                        break;
                    }
                } else {
                    edges[index] = Some(Edge {
                        key,
                        face: f,
                        slot,
                        count: 1,
                        orientation,
                    });
                    unique += 1;
                    break;
                }
                index = (index + 1) & (buckets - 1);
            }
        }
    }
    if edges.iter().flatten().any(|e| e.count != 2) || vertices + triangles.len() != unique + 2 {
        return Err(E::Unresolved);
    }
    for v in 0..vertices {
        if first[v] == usize::MAX {
            return Err(E::Unresolved);
        }
        let start = first[v];
        let mut f = start;
        let mut visited = 0;
        loop {
            w.step()?;
            let slot = triangles[f]
                .iter()
                .position(|&x| x as usize == v)
                .ok_or(E::Unresolved)?;
            f = adjacent[f][slot];
            visited += 1;
            if f == start {
                break;
            }
            if f == usize::MAX || visited >= degree[v] {
                return Err(E::Unresolved);
            }
        }
        if visited != degree[v] {
            return Err(E::Unresolved);
        }
    }
    // Every face must belong to one component, in addition to valid links.
    let mut seen = m.allocate(triangles.len(), w)?;
    seen.resize(triangles.len(), false);
    let mut pending = m.allocate(triangles.len(), w)?;
    pending.push(0);
    seen[0] = true;
    let mut count = 0;
    while let Some(f) = pending.pop() {
        w.step()?;
        count += 1;
        for &next in &adjacent[f] {
            if !seen[next] {
                seen[next] = true;
                pending.push(next);
            }
        }
    }
    if count != triangles.len() {
        return Err(E::Unresolved);
    }
    Ok(())
}
#[cfg(test)]
pub(super) fn tetra_volume(mesh: &DisplayMesh, w: &mut Work<'_>) -> R<[f64; 2]> {
    let origin = point(mesh, 0);
    let mut sum = I::point(0.);
    for tri in mesh.indices.as_chunks::<3>().0 {
        w.charge(64)?;
        let p = tri.map(|id| {
            std::array::from_fn(|k| {
                I::point(mesh.vertices[id as usize * 3 + k]) - I::point(origin[k])
            })
        });
        sum = sum + dot(p[0], cross(p[1], p[2])) / I::point(6.);
    }
    if !sum.valid() || sum.lo <= 0. {
        return Err(E::Unresolved);
    }
    Ok([sum.lo, sum.hi])
}

pub(super) fn volume(mesh: &DisplayMesh, w: &mut Work<'_>) -> R<[f64; 2]> {
    // certify() has proved closed oriented incidence. Divergence of(0,0,z-zref)
    // equals1; its exact triangle integral uses ONLY the projected area and
    // affine mean height, rather than three redundant determinant components.
    let zref = I::point(point(mesh, 0)[2]);
    let mut sum = I::point(0.);
    for tri in mesh.indices.as_chunks::<3>().0 {
        w.charge(32)?;
        let [a, b, c] = tri.map(|id| point(mesh, id as usize).map(I::point));
        let nz = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
        let height = (a[2] - zref) + (b[2] - zref) + (c[2] - zref);
        sum = sum + nz * height / I::point(6.);
    }
    if !sum.valid() || sum.lo <= 0. {
        return Err(E::Unresolved);
    }
    Ok([sum.lo, sum.hi])
}
