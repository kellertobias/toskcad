//! Actual finite cell complex, independent of source event labels or gauges.
use super::*;

fn connected(adjacency: &[Vec<usize>]) -> bool {
    if adjacency.is_empty() {
        return false;
    }
    let mut seen = BTreeSet::new();
    let mut pending = vec![0];
    while let Some(i) = pending.pop() {
        if seen.insert(i) {
            pending.extend(adjacency[i].iter().copied());
        }
    }
    seen.len() == adjacency.len()
}
pub(super) fn read(
    shell: &mut Shell<'_, '_>,
    mouth: Mouth,
    views: &[SourceFaceView],
    entities: &[Option<(u32, bool, MaterialLoop)>],
    work: &mut Work<'_>,
) -> Result<Boundary> {
    require(views.len() == 29 && entities.len() == 29)?;
    let mut raw = BTreeMap::new();
    let mut owners = BTreeMap::<u32, Vec<(usize, bool)>>::new();
    let mut ids = BTreeSet::new();
    for (i, face) in views.iter().enumerate() {
        work.step(1)?;
        for &(id, sense) in &face.uses {
            let [Value::String(_), a, b, geometry, forward] = shell.args(id, "EDGE_CURVE")? else {
                return Err(ReadError::Invalid);
            };
            let vertices = [reference(a)?, reference(b)?];
            require(vertices[0] != vertices[1])?;
            ids.extend(vertices);
            raw.insert(
                id,
                graph::Edge {
                    vertices,
                    geometry: reference(geometry)?,
                    forward: crate::step_read::sense(forward)?,
                    reversed: false,
                },
            );
            owners.entry(id).or_default().push((i, sense));
        }
    }
    require(ids.len() == 49 && raw.len() == 76)?;
    let mut vertices = Vec::new();
    let mut vertex_map = BTreeMap::new();
    for id in ids {
        work.step(1)?;
        let [Value::String(_), point] = shell.args(id, "VERTEX_POINT")? else {
            return Err(ReadError::Invalid);
        };
        let point_entity = reference(point)?;
        let point = shell.point(point_entity)?;
        require(point.into_iter().all(f64::is_finite))?;
        vertex_map.insert(id, vertices.len());
        vertices.push(Vertex {
            entity: id,
            point_entity,
            point,
        });
    }
    let mut dual = vec![Vec::new(); views.len()];
    let mut edges = Vec::new();
    let mut edge_map = BTreeMap::new();
    for (entity, edge) in raw {
        work.step(1)?;
        let [(a, sa), (b, sb)] = owners[&entity].as_slice() else {
            return Err(ReadError::Invalid);
        };
        require(a != b && sa != sb)?;
        dual[*a].push(*b);
        dual[*b].push(*a);
        let owners = [*a, *b];
        let point_ids = edge
            .vertices
            .map(|id| vertices[vertex_map[&id]].point_entity);
        let curve = literal::curve(
            shell,
            edge,
            owners.map(|i| views[i].surface),
            point_ids,
            work,
        )?;
        edge_map.insert(entity, edges.len());
        edges.push(Edge {
            entity,
            vertices: edge.vertices.map(|v| vertex_map[&v]),
            geometry: edge.geometry,
            same_sense: edge.forward,
            owners,
            curve,
        });
    }
    require(connected(&dual))?;
    let mut faces = Vec::new();
    let mut links = vec![BTreeMap::<usize, Vec<usize>>::new(); vertices.len()];
    for (i, view) in views.iter().enumerate() {
        work.step(1)?;
        require(view.uses.len() >= 2)?;
        let uses: Vec<_> = view.uses.iter().map(|&(e, s)| (edge_map[&e], s)).collect();
        let mut seen = BTreeSet::new();
        for (j, &(a, sa)) in uses.iter().enumerate() {
            let (b, sb) = uses[(j + 1) % uses.len()];
            let end = edges[a].vertices[usize::from(sa)];
            let start = edges[b].vertices[usize::from(!sb)];
            require(end == start && seen.insert(end))?;
            links[end].entry(a).or_default().push(b);
            links[end].entry(b).or_default().push(a);
        }
        faces.push(Face {
            material_entity: entities[i].as_ref().map(|e| e.0),
            material_forward: entities[i].as_ref().map_or(mouth.material_forward, |e| e.1),
            material_loop: entities[i].as_ref().map(|e| e.2.clone()),
            cavity_forward: view.forward,
            surface_entity: view.surface,
            support: literal::surface(shell, view.surface, view.forward)?,
            uses,
        });
    }
    for link in links {
        work.step(1)?;
        require(!link.is_empty() && link.values().all(|n| n.len() == 2))?;
        let indexes: BTreeMap<_, _> = link.keys().enumerate().map(|(i, &e)| (e, i)).collect();
        let adjacency: Vec<_> = link
            .values()
            .map(|n| n.iter().map(|e| indexes[e]).collect())
            .collect();
        require(connected(&adjacency))?;
    }
    Ok(Boundary {
        mouth,
        vertices,
        edges,
        faces,
    })
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn an_empty_or_disconnected_link_is_not_a_manifold_cycle() {
        assert!(!connected(&[]));
        assert!(!connected(&[vec![], vec![]]));
    }
}
