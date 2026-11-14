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
    if !mesh.vertices.len().is_multiple_of(3)
        || !mesh.indices.len().is_multiple_of(3)
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
    if first.contains(&usize::MAX)
        || edges.iter().flatten().any(|e| e.count != 2)
        || vertices + triangles.len() != unique + 2
    {
        return Err(E::Unresolved);
    }
    for v in 0..vertices {
        vertex_link(v, first[v], degree[v], triangles, &adjacent, w)?;
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
#[cfg(test)]
mod tests {
    use super::*;
    fn tetra() -> DisplayMesh {
        DisplayMesh {
            vertices: vec![0., 0., 0., 1., 0., 0., 0., 1., 0., 0., 0., 1.],
            indices: vec![0, 2, 1, 0, 1, 3, 1, 2, 3, 2, 0, 3],
            face_indices: vec![0; 4],
        }
    }
    fn check(mesh: &DisplayMesh) -> R<()> {
        let mut m = Memory { used: SCRATCH };
        let mut w = Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        certify(mesh, &mut m, &mut w)
    }
    #[test]
    fn actual_closed_oriented_surface_and_volume_are_separate_from_invalid_incidence() {
        let mesh = tetra();
        check(&mesh).unwrap();
        let mut w = Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        let a = volume(&mesh, &mut w).unwrap();
        let b = tetra_volume(&mesh, &mut w).unwrap();
        assert!(a[0] <= 1. / 6. && 1. / 6. <= a[1]);
        assert!(a[0] <= b[1] && b[0] <= a[1]);
        let mut bad = mesh.clone();
        bad.vertices.pop();
        assert!(check(&bad).is_err());
        let mut bad = mesh.clone();
        bad.indices.pop();
        assert!(check(&bad).is_err());
        let mut bad = mesh.clone();
        bad.face_indices.pop();
        assert!(check(&bad).is_err());
        let mut bad = mesh.clone();
        bad.indices[0] = 99;
        assert!(check(&bad).is_err());
        let mut bad = mesh.clone();
        bad.indices[0] = bad.indices[1];
        assert!(check(&bad).is_err());
        let mut bad = mesh.clone();
        bad.indices.swap(0, 1);
        assert!(check(&bad).is_err());
        let mut bad = mesh.clone();
        bad.indices.truncate(9);
        bad.face_indices.truncate(3);
        assert!(check(&bad).is_err());
        let mut bad = mesh.clone();
        bad.vertices.extend([2., 2., 2.]);
        assert!(check(&bad).is_err());
        let mut bad = mesh.clone();
        for tri in bad.indices.as_chunks_mut::<3>().0 {
            tri.swap(1, 2);
        }
        assert!(volume(&bad, &mut w).is_err());
        assert!(tetra_volume(&bad, &mut w).is_err());
        let mut bad = mesh.clone();
        bad.vertices[0] = f64::NAN;
        assert!(volume(&bad, &mut w).is_err());
        assert!(tetra_volume(&bad, &mut w).is_err());
    }
    #[test]
    fn disconnected_or_pinched_actual_closed_surfaces_do_not_become_one_solid() {
        let mut mesh = tetra();
        let other = tetra();
        mesh.vertices.extend(
            other
                .vertices
                .iter()
                .enumerate()
                .map(|(i, &x)| x + if i % 3 == 0 { 3. } else { 0. }),
        );
        mesh.indices.extend(other.indices.iter().map(|&x| x + 4));
        mesh.face_indices.extend(other.face_indices);
        assert!(check(&mesh).is_err());
        // Two genuine closed tetrahedra touching at one shared vertex have a
        // branched vertex link, not a manifold material junction.
        let mut pinched = tetra();
        pinched
            .vertices
            .extend([2., 0., 0., 0., 2., 0., 0., 0., 2.]);
        let map = [0, 4, 5, 6];
        pinched
            .indices
            .extend(tetra().indices.iter().map(|&x| map[x as usize]));
        pinched.face_indices.extend([0; 4]);
        assert!(check(&pinched).is_err());
    }
}

fn vertex_link(
    v: usize,
    start: usize,
    degree: usize,
    triangles: &[[u32; 3]],
    adjacent: &[[usize; 3]],
    w: &mut Work<'_>,
) -> R<()> {
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
        if f == usize::MAX || visited >= degree {
            return Err(E::Unresolved);
        }
    }
    if visited != degree {
        return Err(E::Unresolved);
    }
    Ok(())
}
#[cfg(test)]
mod link_tests {
    use super::*;
    #[test]
    fn real_vertex_links_reject_stale_missing_branch_and_nonclosing_tables() {
        let triangles = [[0, 1, 2], [0, 2, 3], [0, 3, 1]];
        let mut w = Work {
            remaining: 10000,
            used: 0,
            interrupt: None,
        };
        let adjacent = [[1, 0, 0], [2, 0, 0], [0, 0, 0]];
        vertex_link(0, 0, 3, &triangles, &adjacent, &mut w).unwrap();
        assert!(vertex_link(4, 0, 3, &triangles, &adjacent, &mut w).is_err());
        assert!(vertex_link(0, 0, 2, &triangles, &adjacent, &mut w).is_err());
        assert!(vertex_link(0, 0, 4, &triangles, &adjacent, &mut w).is_err());
        let mut missing = adjacent;
        missing[1][0] = usize::MAX;
        assert!(vertex_link(0, 0, 3, &triangles, &missing, &mut w).is_err());
    }
}
#[cfg(test)]
mod component_tests {
    use super::*;
    #[test]
    fn disjoint_sphere_and_torus_pass_local_links_but_not_single_component() {
        let mut mesh = DisplayMesh {
            vertices: vec![0., 0., 0., 1., 0., 0., 0., 1., 0., 0., 0., 1.],
            indices: vec![0, 2, 1, 0, 1, 3, 1, 2, 3, 2, 0, 3],
            face_indices: vec![0; 4],
        };
        for i in 0..3 {
            for j in 0..3 {
                let a = std::f64::consts::TAU * i as f64 / 3.;
                let b = std::f64::consts::TAU * j as f64 / 3.;
                mesh.vertices.extend([
                    10. + (3. + b.cos()) * a.cos(),
                    (3. + b.cos()) * a.sin(),
                    b.sin(),
                ]);
            }
        }
        let id = |i: usize, j: usize| 4 + (i % 3 * 3 + j % 3) as u32;
        for i in 0..3 {
            for j in 0..3 {
                mesh.indices.extend([
                    id(i, j),
                    id(i + 1, j),
                    id(i + 1, j + 1),
                    id(i, j),
                    id(i + 1, j + 1),
                    id(i, j + 1),
                ]);
                mesh.face_indices.extend([1; 2]);
            }
        }
        let mut m = Memory { used: SCRATCH };
        let mut w = Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        assert!(certify(&mesh, &mut m, &mut w).is_err());
    }
}
