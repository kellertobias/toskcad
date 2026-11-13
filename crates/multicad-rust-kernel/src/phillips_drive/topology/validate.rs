use super::*;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TopologyCertificate {
    pub vertices: usize,
    pub edges: usize,
    pub faces: usize,
    pub euler_characteristic: i64,
}
fn connected(adjacency: &[Vec<usize>]) -> bool {
    if adjacency.is_empty() {
        return false;
    }
    let mut seen = vec![false; adjacency.len()];
    let mut stack = vec![0];
    while let Some(v) = stack.pop() {
        if seen[v] {
            continue;
        }
        seen[v] = true;
        stack.extend(adjacency[v].iter().copied())
    }
    seen.into_iter().all(|v| v)
}
fn vertex_link(incident: &[usize], pairs: &[(usize, usize)]) -> bool {
    if incident.len() < 2 {
        return false;
    }
    let mut adjacency = vec![Vec::new(); incident.len()];
    for &(a, b) in pairs {
        let Some(a) = incident.iter().position(|v| *v == a) else {
            return false;
        };
        let Some(b) = incident.iter().position(|v| *v == b) else {
            return false;
        };
        if a == b {
            return false;
        }
        adjacency[a].push(b);
        adjacency[b].push(a)
    }
    adjacency.iter().all(|a| a.len() == 2) && connected(&adjacency)
}
pub(in crate::phillips_drive) fn certify(
    vertices: &[[I; 3]],
    edges: &[[usize; 2]],
    faces: &[Vec<EdgeUse>],
) -> Result<TopologyCertificate, Error> {
    if vertices.len() < 4 || edges.is_empty() || faces.len() < 4 {
        return Err(Error::Unresolved);
    }
    let mut adjacency = vec![Vec::new(); vertices.len()];
    let mut incident = vec![Vec::new(); vertices.len()];
    let mut links = vec![Vec::new(); vertices.len()];
    let mut uses = vec![Vec::new(); edges.len()];
    for (index, edge) in edges.iter().enumerate() {
        let [a, b] = *edge;
        if a >= vertices.len() || b >= vertices.len() || a == b {
            return Err(Error::Unresolved);
        }
        let distance = std::array::from_fn::<_, 3, _>(|i| vertices[a][i] - vertices[b][i])
            .into_iter()
            .fold(I::point(0.), |sum, v| sum + v.square())
            .sqrt();
        if !distance.valid() || distance.lo <= 0. {
            return Err(Error::Unresolved);
        }
        adjacency[a].push(b);
        adjacency[b].push(a);
        incident[a].push(index);
        incident[b].push(index);
    }
    for (face_index, face) in faces.iter().enumerate() {
        if face.len() < 3 {
            return Err(Error::Unresolved);
        }
        let mut walked = Vec::new();
        for (j, u) in face.iter().enumerate() {
            if u.edge >= edges.len() {
                return Err(Error::Unresolved);
            }
            let [a, b] = edges[u.edge];
            let (start, end) = if u.forward { (a, b) } else { (b, a) };
            if walked.contains(&start) {
                return Err(Error::Unresolved);
            }
            walked.push(start);
            let next = face[(j + 1) % face.len()];
            if next.edge >= edges.len() {
                return Err(Error::Unresolved);
            }
            let next_start = edges[next.edge][usize::from(!next.forward)];
            if end != next_start {
                return Err(Error::Unresolved);
            }
            let previous = face[(j + face.len() - 1) % face.len()];
            links[start].push((previous.edge, u.edge));
            uses[u.edge].push((face_index, u.forward));
        }
    }
    if uses
        .iter()
        .any(|u| u.len() != 2 || u[0].0 == u[1].0 || u[0].1 == u[1].1)
    {
        return Err(Error::Unresolved);
    }
    if !incident
        .iter()
        .zip(&links)
        .all(|(edges, pairs)| vertex_link(edges, pairs))
    {
        return Err(Error::Unresolved);
    }
    if !connected(&adjacency) {
        return Err(Error::Unresolved);
    }
    let euler = vertices.len() as i64 - edges.len() as i64 + faces.len() as i64;
    if euler != 2 {
        return Err(Error::Unresolved);
    }
    Ok(TopologyCertificate {
        vertices: vertices.len(),
        edges: edges.len(),
        faces: faces.len(),
        euler_characteristic: euler,
    })
}
#[cfg(test)]
mod tests;
