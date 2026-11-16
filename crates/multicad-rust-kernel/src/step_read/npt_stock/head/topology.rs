use super::*;
use std::collections::BTreeSet;
fn connected(graph: &[Vec<usize>], work: &mut Work<'_>) -> Result<()> {
    let mut seen = vec![false; graph.len()];
    let mut pending = vec![0];
    while let Some(i) = pending.pop() {
        work.charge(1)?;
        if seen[i] {
            continue;
        }
        seen[i] = true;
        pending.extend(graph[i].iter().copied().filter(|&j| !seen[j]));
    }
    require(seen.into_iter().all(|x| x))
}
pub(super) fn qualify(
    parsed: &Parsed<'_, '_>,
    head: &HeadInput,
    work: &mut Work<'_>,
) -> Result<usize> {
    require(head.vertices.len() == 12 && head.edges.len() == 18 && head.faces.len() == 8)?;
    for p in &head.vertices {
        work.charge(3)?;
        if p.iter().any(|x| !x.is_finite() || x.abs() > 1e7) {
            return Err(Error::Invalid);
        }
    }
    let cap = parsed
        .faces
        .iter()
        .position(|f| f.id == head.auxiliary_cap_entity)
        .ok_or(Error::Invalid)?;
    require(matches!(parsed.faces[cap].support, Support::Plane { .. }))?;
    let annuli = head
        .faces
        .iter()
        .filter(|f| !f.inner.is_empty())
        .collect::<Vec<_>>();
    let [annulus] = annuli.as_slice() else {
        return Err(Error::Unsupported);
    };
    require(annulus.inner.len() == 4)?;
    let original = parsed.faces[cap]
        .coedges
        .map(|e| (parsed.edges[e.edge].id, e.forward));
    require((0..4).any(|shift| (0..4).all(|i| annulus.inner[i] == original[(i + shift) % 4])))?;
    let vertex_count = parsed.vertices.len() + head.vertices.len();
    let edge_count = parsed.edges.len() + head.edges.len();
    let face_count = parsed.faces.len() - 1 + head.faces.len();
    let mut edge_ends = parsed.edges.iter().map(|e| e.vertices).collect::<Vec<_>>();
    for &[a, b] in &head.edges {
        work.charge(1)?;
        require(a < head.vertices.len() && b < head.vertices.len() && a != b)?;
        edge_ends.push([a + parsed.vertices.len(), b + parsed.vertices.len()]);
    }
    let mut loops = Vec::<Vec<Vec<(usize, bool)>>>::new();
    for (i, face) in parsed.faces.iter().enumerate() {
        work.charge(1)?;
        if i != cap {
            loops.push(vec![
                face.coedges.iter().map(|e| (e.edge, e.forward)).collect(),
            ]);
        }
    }
    for face in &head.faces {
        work.charge(1)?;
        require(face.outer.len() >= 3 && face.outer.len() <= 6 && face.inner.len() <= 4)?;
        let outer = face
            .outer
            .iter()
            .map(|&(e, d)| {
                if e < head.edges.len() {
                    Ok((e + parsed.edges.len(), d))
                } else {
                    Err(Error::Invalid)
                }
            })
            .collect::<Result<Vec<_>>>()?;
        let mut rings = vec![outer];
        if !face.inner.is_empty() {
            rings.push(
                face.inner
                    .iter()
                    .map(|&(entity, d)| {
                        parsed
                            .edges
                            .iter()
                            .position(|e| e.id == entity)
                            .map(|i| (i, d))
                            .ok_or(Error::Invalid)
                    })
                    .collect::<Result<Vec<_>>>()?,
            );
        }
        loops.push(rings);
    }
    let mut uses = vec![Vec::<(usize, bool)>::new(); edge_count];
    let mut incident = vec![BTreeSet::new(); vertex_count];
    let mut links = vec![Vec::<(usize, usize)>::new(); vertex_count];
    for (i, &[a, b]) in edge_ends.iter().enumerate() {
        incident[a].insert(i);
        incident[b].insert(i);
    }
    let mut chi = 0isize;
    for (fi, rings) in loops.iter().enumerate() {
        chi += 2 - rings.len() as isize;
        let mut visited = BTreeSet::new();
        for ring in rings {
            for (ri, &(edge, d)) in ring.iter().enumerate() {
                work.charge(4)?;
                let [a, b] = edge_ends[edge];
                let [from, to] = if d { [a, b] } else { [b, a] };
                require(visited.insert(from))?;
                let (next, sense) = ring[(ri + 1) % ring.len()];
                require(to == edge_ends[next][usize::from(!sense)])?;
                uses[edge].push((fi, d));
                links[to].push((edge, next));
            }
        }
    }
    require(vertex_count as isize - edge_count as isize + chi == 2)?;
    let mut neighbors = vec![Vec::new(); face_count];
    for pair in uses {
        work.charge(1)?;
        let [(a, sa), (b, sb)] = pair.as_slice() else {
            return Err(Error::Invalid);
        };
        require(a != b && sa != sb)?;
        neighbors[*a].push(*b);
        neighbors[*b].push(*a);
    }
    connected(&neighbors, work)?;
    for (vi, edges) in incident.iter().enumerate() {
        work.charge(1)?;
        require(edges.len() >= 2)?;
        let edges = edges.iter().copied().collect::<Vec<_>>();
        let mut fan = vec![Vec::new(); edges.len()];
        for &(a, b) in &links[vi] {
            let a = edges.iter().position(|e| *e == a).ok_or(Error::Invalid)?;
            let b = edges.iter().position(|e| *e == b).ok_or(Error::Invalid)?;
            fan[a].push(b);
            fan[b].push(a);
        }
        require(fan.iter().all(|n| n.len() == 2))?;
        connected(&fan, work)?;
    }
    Ok(cap)
}
