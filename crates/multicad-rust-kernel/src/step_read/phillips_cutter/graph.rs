//! Match the literal cell complex bijectively, accepting cyclic loop starts and
//! face ordering without changing a vertex, curve, support or topological use.
use super::*;
#[derive(Clone, Copy)]
pub(super) struct Edge {
    pub vertices: [u32; 2],
    pub geometry: u32,
    pub forward: bool,
    pub reversed: bool,
}
pub(super) struct Face {
    pub surface: u32,
    pub forward: bool,
}
pub(super) struct Graph {
    pub vertices: Vec<(u32, u32, Point)>,
    pub edges: Vec<Edge>,
    pub faces: Vec<Face>,
    pub owners: Vec<[usize; 2]>,
    pub maximum_vertex_error: f64,
}
fn edge(shell: &Shell<'_, '_>, id: u32) -> Result<Edge> {
    let [Value::String(_), a, b, geometry, forward] = shell.args(id, "EDGE_CURVE")? else {
        return Err(ReadError::Invalid);
    };
    let vertices = [reference(a)?, reference(b)?];
    require(vertices[0] != vertices[1])?;
    Ok(Edge {
        vertices,
        geometry: reference(geometry)?,
        forward: sense(forward)?,
        reversed: false,
    })
}
fn cyclic(a: &[(usize, bool)], b: &[(usize, bool)]) -> bool {
    a.len() == b.len()
        && (0..a.len()).any(|start| {
            a.iter()
                .enumerate()
                .all(|(i, &e)| e == b[(start + i) % b.len()])
        })
}
pub(super) fn read(
    shell: &mut Shell<'_, '_>,
    face_ids: &[u32],
    source: &BoundaryTopology<'_>,
    pose: &Pose,
    tolerance: f64,
    work: &mut Work<'_>,
) -> Result<Graph> {
    require(face_ids.len() == source.faces().len())?;
    let mut views = Vec::new();
    for &id in face_ids {
        work.step(1)?;
        let (surface, forward, bounds) = helical::face(shell, id)?;
        let [bound] = bounds.as_slice() else {
            return Err(ReadError::Unsupported);
        };
        let uses = helical::ring(shell, *bound, true)?;
        views.push(SourceFaceView {
            surface,
            forward,
            uses,
        });
    }
    read_views(shell, &views, source, pose, tolerance, work)
}
pub(super) fn read_views(
    shell: &mut Shell<'_, '_>,
    views: &[SourceFaceView],
    source: &BoundaryTopology<'_>,
    pose: &Pose,
    tolerance: f64,
    work: &mut Work<'_>,
) -> Result<Graph> {
    require(views.len() == source.faces().len())?;
    let mut faces = Vec::new();
    let mut edges = BTreeMap::new();
    let mut vertex_ids = BTreeSet::new();
    let mut raw_uses = BTreeMap::<u32, Vec<bool>>::new();
    for view in views {
        work.step(1)?;
        let uses = view.uses.clone();
        for &(id, forward) in &uses {
            let e = edge(shell, id)?;
            vertex_ids.extend(e.vertices);
            edges.insert(id, e);
            raw_uses.entry(id).or_default().push(forward);
        }
        faces.push((
            Face {
                surface: view.surface,
                forward: view.forward,
            },
            uses,
        ));
    }
    require(edges.len() == source.edges().len() && vertex_ids.len() == source.vertices().len())?;
    require(raw_uses.values().all(|u| u.len() == 2 && u[0] != u[1]))?;
    let mut vertices = vec![None; source.vertices().len()];
    let mut mapping = BTreeMap::new();
    let mut maximum_vertex_error = 0_f64;
    let expected = source
        .vertices()
        .iter()
        .map(|p| support::map_point(pose, p.point_mm_bounds().map(support::b)))
        .collect::<Vec<_>>();
    for id in vertex_ids {
        work.step(1)?;
        let [Value::String(_), point_ref] = shell.args(id, "VERTEX_POINT")? else {
            return Err(ReadError::Invalid);
        };
        let point_id = reference(point_ref)?;
        let point = shell.point(point_id)?;
        let mut matches = expected.iter().enumerate().filter_map(|(i, p)| {
            let error = support::norm(std::array::from_fn(|j| support::adjacent(point[j]) - p[j]));
            (error.valid() && error.hi <= tolerance / 16.).then_some((i, error.hi))
        });
        let (index, error) = matches.next().ok_or(ReadError::Unsupported)?;
        require(matches.next().is_none() && vertices[index].is_none())?;
        maximum_vertex_error = maximum_vertex_error.max(error);
        vertices[index] = Some((id, point_id, point));
        mapping.insert(id, index);
    }
    let vertices = vertices
        .into_iter()
        .collect::<Option<Vec<_>>>()
        .ok_or(ReadError::Invalid)?;
    let mut ordered_edges = vec![None; source.edges().len()];
    let mut edge_map = BTreeMap::new();
    for (id, mut raw) in edges {
        work.step(1)?;
        let actual = raw.vertices.map(|v| mapping[&v]);
        let mut candidates = source.edges().iter().enumerate().filter(|(_, e)| {
            let ends = e.vertices();
            actual == ends || actual == [ends[1], ends[0]]
        });
        let (index, expected) = candidates.next().ok_or(ReadError::Unsupported)?;
        require(candidates.next().is_none() && ordered_edges[index].is_none())?;
        raw.reversed = actual != expected.vertices();
        if raw.reversed {
            raw.vertices.reverse();
        }
        edge_map.insert(id, (index, raw.reversed));
        ordered_edges[index] = Some(raw);
    }
    let edges = ordered_edges
        .into_iter()
        .collect::<Option<Vec<_>>>()
        .ok_or(ReadError::Invalid)?;
    let mut ordered_faces = (0..source.faces().len()).map(|_| None).collect::<Vec<_>>();
    let mut owners = vec![Vec::new(); source.edges().len()];
    for (raw, uses) in faces {
        work.step(1)?;
        let actual = uses
            .into_iter()
            .map(|(edge, forward)| {
                let (index, reversed) = edge_map[&edge];
                (index, forward ^ reversed)
            })
            .collect::<Vec<_>>();
        let mut candidates = source.faces().iter().enumerate().filter(|(_, f)| {
            let expected = f
                .edge_uses()
                .iter()
                .map(|u| (u.edge, u.forward))
                .collect::<Vec<_>>();
            cyclic(&actual, &expected)
        });
        let (index, _) = candidates.next().ok_or(ReadError::Unsupported)?;
        require(candidates.next().is_none() && ordered_faces[index].is_none())?;
        for &(edge, _) in &actual {
            owners[edge].push(index);
        }
        ordered_faces[index] = Some(raw);
    }
    let faces = ordered_faces
        .into_iter()
        .collect::<Option<Vec<_>>>()
        .ok_or(ReadError::Invalid)?;
    let owners = owners.into_iter().map(|o| [o[0], o[1]]).collect();
    Ok(Graph {
        vertices,
        edges,
        faces,
        owners,
        maximum_vertex_error,
    })
}
