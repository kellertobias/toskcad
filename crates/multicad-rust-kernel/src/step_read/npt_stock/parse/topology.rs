use super::*;
fn vertex<'input>(
    records: &BTreeMap<u32, Vec<Entity<'input>>>,
    id: u32,
    points: &[Point<'input>],
    vertices: &mut Vec<Vertex>,
    budget: &mut Budget<'_>,
) -> Result<usize> {
    budget.step(16)?;
    if let Some(index) = vertices.iter().position(|v| v.id == id) {
        return Ok(index);
    }
    if vertices.len() == 12 {
        return Err(Error::Unsupported);
    }
    let entity = kind(records, id, "VERTEX_POINT", budget)?;
    let [Value::String(_), point] = entity.args.as_slice() else {
        return Err(Error::Invalid);
    };
    let point = point_index(points, reference(point)?, budget)?;
    let index = vertices.len();
    vertices.push(Vertex { id, point });
    Ok(index)
}
fn edge<'input>(
    records: &BTreeMap<u32, Vec<Entity<'input>>>,
    id: u32,
    points: &[Point<'input>],
    vertices: &mut Vec<Vertex>,
    curves: &mut Vec<Curve>,
    edges: &mut Vec<Edge>,
    budget: &mut Budget<'_>,
) -> Result<usize> {
    budget.step(32)?;
    if let Some(index) = edges.iter().position(|e| e.id == id) {
        return Ok(index);
    }
    if edges.len() == 20 {
        return Err(Error::Unsupported);
    }
    let entity = kind(records, id, "EDGE_CURVE", budget)?;
    let [Value::String(_), a, b, geometry, Value::Enum("T")] = entity.args.as_slice() else {
        return Err(Error::Unsupported);
    };
    let endpoints = [
        vertex(records, reference(a)?, points, vertices, budget)?,
        vertex(records, reference(b)?, points, vertices, budget)?,
    ];
    if endpoints[0] == endpoints[1] {
        return Err(Error::Unsupported);
    }
    let surface = kind(records, reference(geometry)?, "SURFACE_CURVE", budget)?;
    let [
        Value::String(_),
        curve_id,
        pcurves,
        Value::Enum("PCURVE_S1"),
    ] = surface.args.as_slice()
    else {
        return Err(Error::Unsupported);
    };
    let curve_id = reference(curve_id)?;
    if curves.iter().any(|curve| curve.id == curve_id) {
        return Err(Error::Unsupported);
    }
    let curve = super::curve(records, points, curve_id, budget)?;
    if curve.poles[0] != vertices[endpoints[0]].point
        || *curve.poles.last().expect("nonempty C11 curve") != vertices[endpoints[1]].point
    {
        return Err(Error::Unsupported);
    }
    let pcurves = list(pcurves)?;
    let [p0, p1] = pcurves else {
        return Err(Error::Unsupported);
    };
    let pcurves = [reference(p0)?, reference(p1)?];
    if pcurves[0] == pcurves[1] {
        return Err(Error::Invalid);
    }
    let curve_index = curves.len();
    curves.push(curve);
    let index = edges.len();
    edges.push(Edge {
        id,
        curve: curve_index,
        vertices: endpoints,
        pcurves,
    });
    Ok(index)
}
pub(super) fn wire<'input>(
    records: &BTreeMap<u32, Vec<Entity<'input>>>,
    bound: u32,
    points: &[Point<'input>],
    vertices: &mut Vec<Vertex>,
    curves: &mut Vec<Curve>,
    edges: &mut Vec<Edge>,
    budget: &mut Budget<'_>,
) -> Result<[Coedge; 4]> {
    let bound = kind(records, bound, "FACE_OUTER_BOUND", budget)?;
    let [Value::String(_), wire, Value::Enum("T")] = bound.args.as_slice() else {
        return Err(Error::Unsupported);
    };
    let wire = kind(records, reference(wire)?, "EDGE_LOOP", budget)?;
    let [Value::String(_), values] = wire.args.as_slice() else {
        return Err(Error::Invalid);
    };
    let values = list(values)?;
    if values.len() != 4 {
        return Err(Error::Unsupported);
    }
    let mut coedges = [Coedge {
        edge: 0,
        forward: true,
    }; 4];
    for (index, value) in values.iter().enumerate() {
        let entity = kind(records, reference(value)?, "ORIENTED_EDGE", budget)?;
        let [
            Value::String(_),
            Value::Derived,
            Value::Derived,
            edge_id,
            orientation,
        ] = entity.args.as_slice()
        else {
            return Err(Error::Unsupported);
        };
        let edge = edge(
            records,
            reference(edge_id)?,
            points,
            vertices,
            curves,
            edges,
            budget,
        )?;
        if coedges[..index].iter().any(|c| c.edge == edge) {
            return Err(Error::Invalid);
        }
        coedges[index] = Coedge {
            edge,
            forward: boolean(orientation)?,
        };
    }
    let mut starts = [0usize; 4];
    for i in 0..4 {
        let a = coedges[i];
        let b = coedges[(i + 1) % 4];
        let start = edges[a.edge].vertices[usize::from(!a.forward)];
        let end = edges[a.edge].vertices[usize::from(a.forward)];
        let next = edges[b.edge].vertices[usize::from(!b.forward)];
        if end != next || starts[..i].contains(&start) {
            return Err(Error::Invalid);
        }
        starts[i] = start;
    }
    Ok(coedges)
}
pub(super) fn shell(
    vertices: &[Vertex],
    edges: &[Edge],
    faces: &[Face],
    budget: &mut Budget<'_>,
) -> Result<()> {
    let mut uses = [[(usize::MAX, false); 2]; 20];
    let mut counts = [0usize; 20];
    for (face_index, face) in faces.iter().enumerate() {
        for coedge in face.coedges {
            budget.step(8)?;
            let count = &mut counts[coedge.edge];
            if *count == 2 {
                return Err(Error::Invalid);
            }
            uses[coedge.edge][*count] = (face_index, coedge.forward);
            *count += 1;
        }
    }
    for i in 0..edges.len() {
        if counts[i] != 2 || uses[i][0].0 == uses[i][1].0 || uses[i][0].1 == uses[i][1].1 {
            return Err(Error::Invalid);
        }
    }
    let mut visited = [false; 10];
    let mut pending = [0usize; 10];
    let mut len = 1;
    visited[0] = true;
    while len > 0 {
        len -= 1;
        let face = pending[len];
        for coedge in faces[face].coedges {
            budget.step(8)?;
            for (neighbor, _) in uses[coedge.edge] {
                if !visited[neighbor] {
                    visited[neighbor] = true;
                    pending[len] = neighbor;
                    len += 1;
                }
            }
        }
    }
    if visited[..faces.len()].iter().any(|v| !*v) {
        return Err(Error::Invalid);
    }
    // Every face is a disk with four distinct vertices. Opposite edge uses
    // and one cyclic vertex link establish a closed orientable 2-manifold;
    // geometry, projected embedding and occupied volume remain separate.
    for vertex in 0..vertices.len() {
        let mut incident_edges = [0usize; 4];
        let mut degree = 0;
        for (index, edge) in edges.iter().enumerate() {
            if edge.vertices.contains(&vertex) {
                if degree == 4 {
                    return Err(Error::Unsupported);
                }
                incident_edges[degree] = index;
                degree += 1;
            }
        }
        if degree < 3 {
            return Err(Error::Unsupported);
        }
        let mut incident = [false; 10];
        for edge in &incident_edges[..degree] {
            for (face, _) in uses[*edge] {
                incident[face] = true;
            }
        }
        if incident.iter().filter(|v| **v).count() != degree {
            return Err(Error::Invalid);
        }
        for (face, present) in incident.iter().enumerate() {
            if *present
                && incident_edges[..degree]
                    .iter()
                    .filter(|edge| uses[**edge].iter().any(|(f, _)| *f == face))
                    .count()
                    != 2
            {
                return Err(Error::Invalid);
            }
        }
        let start = incident
            .iter()
            .position(|v| *v)
            .expect("positive vertex degree");
        let mut seen = [false; 10];
        seen[start] = true;
        let mut changed = true;
        while changed {
            changed = false;
            for edge in &incident_edges[..degree] {
                budget.step(8)?;
                let [(a, _), (b, _)] = uses[*edge];
                if seen[a] != seen[b] {
                    seen[a] = true;
                    seen[b] = true;
                    changed = true;
                }
            }
        }
        if incident.iter().zip(seen).any(|(a, b)| *a && !b) {
            return Err(Error::Invalid);
        }
    }
    Ok(())
}
