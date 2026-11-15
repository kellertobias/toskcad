//! Finite physical source owners. The through bore makes one connected genus1
//! shell: three annular faces, with V-E+sum(face Euler characteristic)=0.
use super::*;
use std::collections::{BTreeMap, BTreeSet};
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Vertex {
    Stock { end: u8, corner: u8 },
    Crown(u8),
    ReliefStart(u8),
    Thread { end: u8, junction: u8 },
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Edge {
    StockEnd { end: u8, side: u8 },
    Rail(u8),
    SplitRail { corner: u8, end: u8 },
    Crown(u8),
    ReliefStart(u8),
    ThreadRing { end: u8, role: u8 },
    Helix(u8),
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Coedge {
    pub edge: Edge,
    pub forward: bool,
}
#[derive(Clone, Debug)]
pub struct Face {
    index: u32,
    loops: Vec<Vec<Coedge>>,
}
impl Face {
    pub fn index(&self) -> u32 {
        self.index
    }
    pub fn loops(&self) -> &[Vec<Coedge>] {
        &self.loops
    }
}
#[derive(Clone, Debug)]
pub struct Boundary {
    faces: Vec<Face>,
    edges: Vec<Edge>,
    vertices: Vec<Vertex>,
    work: usize,
    bytes: usize,
}
impl Boundary {
    pub fn faces(&self) -> &[Face] {
        &self.faces
    }
    pub fn edges(&self) -> &[Edge] {
        &self.edges
    }
    pub fn vertices(&self) -> &[Vertex] {
        &self.vertices
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_bytes(&self) -> usize {
        self.bytes
    }
}
fn stock(end: u8, corner: u8) -> Vertex {
    Vertex::Stock { end, corner }
}
fn root(corner: u8, end: u8) -> u8 {
    match (corner, end) {
        (4, 0) => 1,
        (4, 1) => 0,
        (5, 0) => 2,
        (5, 1) => 3,
        _ => unreachable!("private split rail"),
    }
}
fn junction(role: u8) -> u8 {
    [2, 3, 0, 1][role as usize]
}
impl Edge {
    pub fn endpoints(self) -> Result<[Vertex; 2], Error> {
        Ok(match self {
            Self::StockEnd { end, side } if end < 2 && side < 8 => {
                [stock(end, side), stock(end, (side + 1) % 8)]
            }
            Self::Rail(corner) if corner < 8 && corner != 4 && corner != 5 => {
                [stock(0, corner), stock(1, corner)]
            }
            Self::SplitRail { corner, end } if [4, 5].contains(&corner) && end < 2 => {
                if end == 0 {
                    [stock(0, corner), Vertex::Crown(root(corner, 0))]
                } else {
                    [Vertex::Crown(root(corner, 1)), stock(1, corner)]
                }
            }
            Self::Crown(i) if i < 4 => [Vertex::Crown(i), Vertex::Crown((i + 1) % 4)],
            Self::ReliefStart(i) if i < 4 => {
                [Vertex::ReliefStart(i), Vertex::ReliefStart((i + 1) % 4)]
            }
            Self::ThreadRing { end, role } if end < 2 && role < 4 => {
                let j = junction(role);
                [
                    Vertex::Thread { end, junction: j },
                    Vertex::Thread {
                        end,
                        junction: (j + 1) % 4,
                    },
                ]
            }
            Self::Helix(junction) if junction < 4 => [
                Vertex::Thread { end: 0, junction },
                Vertex::Thread { end: 1, junction },
            ],
            _ => return Err(Error::Invalid),
        })
    }
}
impl Nut {
    pub fn boundary(
        &self,
        maximum_work: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Boundary, Error> {
        checkpoint(interrupt)?;
        if maximum_work < 131072 {
            return Err(Error::Limit);
        }
        if self.thread.spec().starts != 1 {
            return Err(Error::Unsupported);
        }
        let mut edges = Vec::new();
        for end in 0..2 {
            for side in 0..8 {
                edges.push(Edge::StockEnd { end, side })
            }
        }
        for corner in [0, 1, 2, 3, 6, 7] {
            edges.push(Edge::Rail(corner));
        }
        for corner in [4, 5] {
            for end in 0..2 {
                edges.push(Edge::SplitRail { corner, end });
            }
        }
        for i in 0..4 {
            edges.extend([Edge::Crown(i), Edge::ReliefStart(i), Edge::Helix(i)]);
        }
        for end in 0..2 {
            for role in 0..4 {
                edges.push(Edge::ThreadRing { end, role });
            }
        }
        edges.sort();
        let mut lookup = BTreeMap::new();
        for &e in &edges {
            let [a, b] = e.endpoints()?;
            if lookup
                .insert(
                    (a, b),
                    Coedge {
                        edge: e,
                        forward: true,
                    },
                )
                .is_some()
                || lookup
                    .insert(
                        (b, a),
                        Coedge {
                            edge: e,
                            forward: false,
                        },
                    )
                    .is_some()
            {
                return Err(Error::Unresolved);
            }
        }
        let ring = |v: Vec<Vertex>| -> Result<Vec<Coedge>, Error> {
            (0..v.len())
                .map(|i| {
                    lookup
                        .get(&(v[i], v[(i + 1) % v.len()]))
                        .copied()
                        .ok_or(Error::Unresolved)
                })
                .collect()
        };
        let mut faces = (0..17)
            .map(|index| Face {
                index,
                loops: vec![],
            })
            .collect::<Vec<_>>();
        for side in [1u8, 2, 6, 7] {
            faces[side as usize].loops.push(ring(vec![
                stock(0, side),
                stock(1, side),
                stock(1, (side + 1) % 8),
                stock(0, (side + 1) % 8),
            ])?);
        }
        let mut outer = ring(vec![stock(0, 0), stock(1, 0), stock(1, 1), stock(0, 1)])?;
        flip(&mut outer);
        let mut inner = ring(
            (0..4)
                .map(|junction| Vertex::Thread { end: 0, junction })
                .collect(),
        )?;
        if !self.thread.spec().left_hand {
            flip(&mut inner)
        }
        faces[0].loops = vec![outer, inner];
        faces[3].loops.push(ring(vec![
            stock(0, 3),
            stock(1, 3),
            stock(1, 4),
            Vertex::Crown(0),
            Vertex::Crown(1),
            stock(0, 4),
        ])?);
        faces[5].loops.push(ring(vec![
            stock(0, 6),
            stock(1, 6),
            stock(1, 5),
            Vertex::Crown(3),
            Vertex::Crown(2),
            stock(0, 5),
        ])?);
        faces[4].loops.push(ring(vec![
            stock(0, 5),
            Vertex::Crown(2),
            Vertex::Crown(1),
            stock(0, 4),
        ])?);
        faces[10].loops.push(ring(vec![
            stock(1, 5),
            stock(1, 4),
            Vertex::Crown(0),
            Vertex::Crown(3),
        ])?);
        for end in 0..2 {
            faces[8 + end as usize]
                .loops
                .push(ring((0..8).map(|corner| stock(end, corner)).collect())?);
        }
        for role in 0..4 {
            let j = junction(role);
            let next = (j + 1) % 4;
            faces[11 + role as usize].loops.push(ring(vec![
                Vertex::Thread {
                    end: 0,
                    junction: j,
                },
                Vertex::Thread {
                    end: 0,
                    junction: next,
                },
                Vertex::Thread {
                    end: 1,
                    junction: next,
                },
                Vertex::Thread {
                    end: 1,
                    junction: j,
                },
            ])?);
        }
        let outer = ring((0..4).map(Vertex::ReliefStart).collect())?;
        let mut inner = ring(
            (0..4)
                .map(|junction| Vertex::Thread { end: 1, junction })
                .collect(),
        )?;
        if self.thread.spec().left_hand {
            flip(&mut inner)
        }
        faces[15].loops = vec![outer, inner];
        let mut bottom = ring((0..4).map(Vertex::ReliefStart).collect())?;
        flip(&mut bottom);
        faces[16].loops = vec![bottom, ring((0..4).map(Vertex::Crown).collect())?];
        let uses = collect(&faces);
        if uses.values().any(|v| v.len() != 2) {
            return Err(Error::Unresolved);
        }
        // Propagate orientation from the actual outward bottom plane; independently
        // seeded annular loop orientations must agree around the entire face graph.
        let mut orientation = vec![None; 17];
        orientation[0] = Some(false);
        loop {
            let mut changed = false;
            for pair in uses.values() {
                let [(a, da), (b, db)] = [pair[0], pair[1]];
                if let Some(flip_a) = orientation[a] {
                    let needed = flip_a ^ (da == db);
                    if orientation[b].is_some_and(|value| value != needed) {
                        return Err(Error::Unresolved);
                    }
                    if orientation[b].is_none() {
                        orientation[b] = Some(needed);
                        changed = true;
                    }
                } else if let Some(flip_b) = orientation[b] {
                    orientation[a] = Some(flip_b ^ (da == db));
                    changed = true;
                }
            }
            if !changed {
                break;
            }
        }
        if orientation.iter().any(Option::is_none) {
            return Err(Error::Unresolved);
        }
        for (face, flip_) in faces.iter_mut().zip(orientation) {
            if flip_ == Some(true) {
                for ring in &mut face.loops {
                    flip(ring);
                }
            }
        }
        let mut vertices = BTreeSet::new();
        let mut fans = BTreeMap::<Vertex, BTreeMap<Edge, Edge>>::new();
        for f in &faces {
            for ring in &f.loops {
                if ring.len() < 3 {
                    return Err(Error::Unresolved);
                }
                for i in 0..ring.len() {
                    let e = ring[i];
                    let next = ring[(i + 1) % ring.len()];
                    let ends = directed(e)?;
                    if ends[1] != directed(next)?[0] || ends[0] == ends[1] {
                        return Err(Error::Unresolved);
                    }
                    vertices.extend(ends);
                    if fans
                        .entry(ends[1])
                        .or_default()
                        .insert(e.edge, next.edge)
                        .is_some()
                    {
                        return Err(Error::Unresolved);
                    }
                }
            }
        }
        if collect(&faces)
            .values()
            .any(|v| v.len() != 2 || v[0].1 == v[1].1)
        {
            return Err(Error::Unresolved);
        }
        for fan in fans.values() {
            let first = *fan.keys().next().ok_or(Error::Unresolved)?;
            let mut current = first;
            let mut count = 0;
            loop {
                current = *fan.get(&current).ok_or(Error::Unresolved)?;
                count += 1;
                if current == first {
                    break;
                }
                if count > fan.len() {
                    return Err(Error::Unresolved);
                }
            }
            if count != fan.len() {
                return Err(Error::Unresolved);
            }
        }
        let chi = faces
            .iter()
            .map(|f| 2 - f.loops.len() as isize)
            .sum::<isize>();
        if vertices.len() != 32
            || edges.len() != 46
            || vertices.len() as isize - edges.len() as isize + chi != 0
        {
            return Err(Error::Unresolved);
        }
        let vertices = vertices.into_iter().collect::<Vec<_>>();
        let bytes = std::mem::size_of::<Boundary>()
            + vertices.capacity() * std::mem::size_of::<Vertex>()
            + edges.capacity() * std::mem::size_of::<Edge>()
            + faces.capacity() * std::mem::size_of::<Face>()
            + faces
                .iter()
                .map(|f| {
                    f.loops.capacity() * std::mem::size_of::<Vec<Coedge>>()
                        + f.loops
                            .iter()
                            .map(|r| r.capacity() * std::mem::size_of::<Coedge>())
                            .sum::<usize>()
                })
                .sum::<usize>();
        if bytes + 131072 > 64 * 1024 * 1024 {
            return Err(Error::Limit);
        }
        checkpoint(interrupt)?;
        Ok(Boundary {
            faces,
            edges,
            vertices,
            work: 131072,
            bytes,
        })
    }
}
fn flip(ring: &mut [Coedge]) {
    ring.reverse();
    for e in ring {
        e.forward = !e.forward;
    }
}
fn directed(e: Coedge) -> Result<[Vertex; 2], Error> {
    let mut ends = e.edge.endpoints()?;
    if !e.forward {
        ends.reverse()
    }
    Ok(ends)
}
fn collect(faces: &[Face]) -> BTreeMap<Edge, Vec<(usize, bool)>> {
    let mut out = BTreeMap::<Edge, Vec<(usize, bool)>>::new();
    for (i, f) in faces.iter().enumerate() {
        for ring in &f.loops {
            for e in ring {
                out.entry(e.edge).or_default().push((i, e.forward));
            }
        }
    }
    out
}
