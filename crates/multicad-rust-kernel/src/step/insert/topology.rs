//! Actual local tensor/annulus incidence; no geometry or mass admission alone.
use super::*;
use crate::wood_insert::FaceKey;
mod winding;
pub use winding::WindingBounds;
#[cfg(test)]
mod tests;
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum RulingKey {
    Source(BoundaryCurveKey),
    Cell { side: u8, cell: usize, station: u32 },
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EdgeGeometry {
    Row(usize),
    Ruling,
}
#[derive(Clone, Copy, Debug)]
pub struct ShellEdge {
    pub ends: [usize; 2],
    pub geometry: EdgeGeometry,
}
#[derive(Clone, Copy, Debug)]
pub struct ShellCoedge {
    pub edge: usize,
    pub forward: bool,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FaceGeometry {
    Tensor(usize),
    Cap(u8),
}
pub struct ShellFace {
    pub geometry: FaceGeometry,
    pub same_sense: bool,
    pub loops: Vec<Vec<ShellCoedge>>,
}
/// Closed oriented genus-one INCIDENCE. Actual embedding, finite supports,
/// winding, occupied material and file readback require separate proofs.
pub struct TopologyPlan<'s, 'p, 'a> {
    surface: &'s SurfacePlan<'p, 'a>,
    edges: Vec<ShellEdge>,
    faces: Vec<ShellFace>,
    used_vertices: usize,
    work: usize,
    bytes: usize,
}
impl<'s, 'p, 'a> TopologyPlan<'s, 'p, 'a> {
    pub fn surface(&self) -> &'s SurfacePlan<'p, 'a> {
        self.surface
    }
    pub fn edges(&self) -> &[ShellEdge] {
        &self.edges
    }
    pub fn faces(&self) -> &[ShellFace] {
        &self.faces
    }
    pub fn used_vertices(&self) -> usize {
        self.used_vertices
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_bytes(&self) -> usize {
        self.bytes
    }
    pub fn prepare(
        surface: &'s SurfacePlan<'p, 'a>,
        maximum_work: usize,
        maximum_bytes: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        let mut w = Work::new(maximum_work, interrupt)?;
        w.charge(surface.work_used() + 1)?;
        if maximum_bytes == 0 || maximum_bytes > 64 * 1024 * 1024 {
            return Err(Error::Invalid);
        }
        let mut memory = Memory::new(surface.retained_bytes(), maximum_bytes)?;
        memory.add(size_of::<Self>())?;
        let mut edges =
            memory.allocate(surface.rows().len() + 2 * surface.patches().len(), &mut w)?;
        for (i, row) in surface.rows().iter().enumerate() {
            w.charge(1)?;
            edges.push(ShellEdge {
                ends: row.endpoints(),
                geometry: EdgeGeometry::Row(i),
            });
        }
        let mut ruling_keys = memory.allocate(surface.patches().len() * 2, &mut w)?;
        let mut faces = memory.allocate(surface.patches().len() + 2, &mut w)?;
        for (patch_index, patch) in surface.patches().iter().enumerate() {
            let slab = &surface.domain().slabs()[patch.slab];
            let [a, b] = patch.rows.map(|i| &surface.rows()[i]);
            let mut wire = memory.allocate(4, &mut w)?;
            if a.endpoints()[0] != b.endpoints()[0] {
                wire.push(ruling(
                    surface,
                    slab,
                    0,
                    [a.endpoints()[0], b.endpoints()[0]],
                    &mut ruling_keys,
                    &mut edges,
                    &mut w,
                )?);
            }
            wire.push(ShellCoedge {
                edge: patch.rows[1],
                forward: true,
            });
            if a.endpoints()[1] != b.endpoints()[1] {
                wire.push(ruling(
                    surface,
                    slab,
                    1,
                    [b.endpoints()[1], a.endpoints()[1]],
                    &mut ruling_keys,
                    &mut edges,
                    &mut w,
                )?);
            }
            wire.push(ShellCoedge {
                edge: patch.rows[0],
                forward: false,
            });
            if !patch.same_sense {
                wire.reverse();
                for c in &mut wire {
                    c.forward = !c.forward;
                }
            }
            let mut loops = memory.allocate(1, &mut w)?;
            loops.push(wire);
            faces.push(ShellFace {
                geometry: FaceGeometry::Tensor(patch_index),
                same_sense: patch.same_sense,
                loops,
            });
        }
        for (end, key) in [(0, FaceKey::StartCap), (1, FaceKey::EndCap)] {
            let source_loops = surface
                .domain()
                .boundary()
                .face_loops(key)
                .expect("actual paired annular cap");
            let mut loops = memory.allocate(source_loops.len(), &mut w)?;
            for source_loop in source_loops {
                // Each source curve contributes only its retained shared rows;
                // no disk scaffold or circular-envelope substitute is emitted.
                w.charge(source_loop.len() * surface.rows().len())?;
                let count = source_loop
                    .iter()
                    .map(|c| {
                        surface
                            .rows()
                            .iter()
                            .filter(|r| r.key().curve == c.curve)
                            .count()
                    })
                    .sum();
                let mut wire = memory.allocate(count, &mut w)?;
                for c in source_loop {
                    let begin = surface.rows().partition_point(|r| r.key().curve < c.curve);
                    let finish = surface.rows().partition_point(|r| r.key().curve <= c.curve);
                    w.charge(
                        2 * (usize::BITS - surface.rows().len().leading_zeros() + 1) as usize,
                    )?;
                    let rows = &surface.rows()[begin..finish];
                    let original = surface
                        .domain()
                        .boundary()
                        .curve_endpoints(c.curve)
                        .expect("actual cap curve endpoints");
                    let first = surface.vertices()
                        [rows.first().expect("nonempty cap row coverage").endpoints()[0]]
                        .owner();
                    let last = surface.vertices()[rows.last().unwrap().endpoints()[1]].owner();
                    let forward = cap_forward([first, last], original, c.forward)?;
                    for j in 0..rows.len() {
                        w.charge(1)?;
                        wire.push(ShellCoedge {
                            edge: begin + if forward { j } else { rows.len() - 1 - j },
                            forward,
                        });
                    }
                }
                loops.push(wire);
            }
            faces.push(ShellFace {
                geometry: FaceGeometry::Cap(end),
                same_sense: end == 1,
                loops,
            });
        }
        memory.release(&ruling_keys);
        drop(ruling_keys);
        let used_vertices = verify(
            surface.vertices().len(),
            &edges,
            &faces,
            &mut memory,
            &mut w,
        )?;
        w.charge(1)?;
        Ok(Self {
            surface,
            edges,
            faces,
            used_vertices,
            work: w.used,
            bytes: memory.used,
        })
    }
}
#[allow(clippy::too_many_arguments)]
fn ruling(
    surface: &SurfacePlan<'_, '_>,
    slab: &Slab,
    end: usize,
    ordered: [usize; 2],
    keys: &mut Vec<(RulingKey, usize)>,
    edges: &mut Vec<ShellEdge>,
    w: &mut Work<'_>,
) -> Result<ShellCoedge, Error> {
    let mut source = None;
    if let [VertexOwner::Original(a), VertexOwner::Original(b)] =
        ordered.map(|i| surface.vertices()[i].owner())
    {
        let component = surface.domain().boundary().source_boundary(slab.side);
        for c in component
            .cell_coedges(slab.cell)
            .expect("original source cell")
        {
            w.charge(1)?;
            let key = curve_key(slab.side, c.curve as u32);
            let curve = surface.domain().boundary().curve(key).unwrap();
            let v = curve.affine_rzv_bounds()[1][2];
            let endpoints = curve.endpoints().map(|i| vertex_key(slab.side, i));
            if v.lower == 0. && v.upper == 0. && (endpoints == [a, b] || endpoints == [b, a]) {
                source = Some(key);
                break;
            }
        }
    }
    let key = source.map(RulingKey::Source).unwrap_or(RulingKey::Cell {
        side: match slab.side {
            Side::Outer => 0,
            Side::Inner => 1,
        },
        cell: slab.cell,
        station: slab.stations[end],
    });
    let mut low = 0;
    let mut high = keys.len();
    let mut found = None;
    while low < high {
        w.charge(1)?;
        let m = low + (high - low) / 2;
        match keys[m].0.cmp(&key) {
            std::cmp::Ordering::Equal => {
                found = Some(keys[m].1);
                break;
            }
            std::cmp::Ordering::Less => low = m + 1,
            std::cmp::Ordering::Greater => high = m,
        }
    }
    let ends = if ordered[0] < ordered[1] {
        ordered
    } else {
        [ordered[1], ordered[0]]
    };
    let index = if let Some(i) = found {
        if edges[i].ends != ends {
            return Err(Error::Unresolved);
        }
        i
    } else {
        w.charge(keys.len() - low + 1)?;
        let i = edges.len();
        edges.push(ShellEdge {
            ends,
            geometry: EdgeGeometry::Ruling,
        });
        keys.insert(low, (key, i));
        i
    };
    Ok(ShellCoedge {
        edge: index,
        forward: ordered == ends,
    })
}
#[derive(Clone, Copy)]
struct Use {
    face: usize,
    forward: bool,
}
fn verify(
    vertices: usize,
    edges: &[ShellEdge],
    faces: &[ShellFace],
    memory: &mut Memory,
    w: &mut Work<'_>,
) -> Result<usize, Error> {
    let mut uses = memory.allocate(edges.len(), w)?;
    uses.resize(edges.len(), [None::<Use>; 2]);
    let mut links = memory.allocate(edges.len(), w)?;
    links.resize(edges.len(), [[None::<usize>; 2]; 2]);
    let mut first = memory.allocate(vertices, w)?;
    first.resize(vertices, None::<usize>);
    let mut counts = memory.allocate(vertices, w)?;
    counts.resize(vertices, 0usize);
    let mut parents = memory.allocate(faces.len(), w)?;
    parents.extend(0..faces.len());
    for (i, e) in edges.iter().enumerate() {
        w.charge(1)?;
        if e.ends[0] == e.ends[1] || e.ends.iter().any(|&v| v >= vertices) {
            return Err(Error::Unresolved);
        }
        for &v in &e.ends {
            first[v].get_or_insert(i);
            counts[v] += 1;
        }
    }
    let mut euler_faces = 0isize;
    for (face, f) in faces.iter().enumerate() {
        if f.loops.is_empty() {
            return Err(Error::Unresolved);
        }
        euler_faces += 2 - f.loops.len() as isize;
        for wire in &f.loops {
            if wire.len() < 3 {
                return Err(Error::Unresolved);
            }
            for (i, c) in wire.iter().enumerate() {
                w.charge(1)?;
                let Some(edge) = edges.get(c.edge) else {
                    return Err(Error::Unresolved);
                };
                let u = &mut uses[c.edge];
                if u[0].is_none() {
                    u[0] = Some(Use {
                        face,
                        forward: c.forward,
                    });
                } else if u[1].is_none() {
                    u[1] = Some(Use {
                        face,
                        forward: c.forward,
                    });
                } else {
                    return Err(Error::Unresolved);
                }
                let previous = wire[(i + wire.len() - 1) % wire.len()];
                let Some(p) = edges.get(previous.edge) else {
                    return Err(Error::Unresolved);
                };
                let v = edge.ends[usize::from(!c.forward)];
                if p.ends[usize::from(previous.forward)] != v || previous.edge == c.edge {
                    return Err(Error::Unresolved);
                }
                add_link(&mut links[c.edge][usize::from(!c.forward)], previous.edge);
                add_link(
                    &mut links[previous.edge][usize::from(previous.forward)],
                    c.edge,
                );
            }
        }
    }
    for u in &uses {
        w.charge(1)?;
        let [Some(a), Some(b)] = u else {
            return Err(Error::Unresolved);
        };
        if a.forward == b.forward {
            return Err(Error::Unresolved);
        }
        let x = root(&mut parents, a.face, w)?;
        let y = root(&mut parents, b.face, w)?;
        parents[x] = y;
    }
    let connected = root(&mut parents, 0, w)?;
    for i in 1..faces.len() {
        if root(&mut parents, i, w)? != connected {
            return Err(Error::Unresolved);
        }
    }
    let mut used = 0;
    for v in 0..vertices {
        let Some(start) = first[v] else { continue };
        used += 1;
        let mut current = start;
        let mut previous = None;
        let mut visited = 0;
        loop {
            w.charge(1)?;
            let endpoint = usize::from(edges[current].ends[1] == v);
            // Every edge earned exactly two coedge uses above. Each use
            // inserts exactly one link at each endpoint, so both slots exist.
            let [a, b] = links[current][endpoint].map(|x| x.expect("two earned coedge uses"));
            let next = if Some(a) == previous { b } else { a };
            visited += 1;
            if next == start {
                break;
            }
            // A finite degree-two link graph is a union of cycles; following
            // the non-predecessor edge returns to start. The count below
            // independently rejects multiple cycles at a pinched vertex.
            previous = Some(current);
            current = next;
        }
        if visited != counts[v] {
            return Err(Error::Unresolved);
        }
    }
    if used as isize - edges.len() as isize + euler_faces != 0 {
        return Err(Error::Unresolved);
    }
    memory.release(&uses);
    memory.release(&links);
    memory.release(&first);
    memory.release(&counts);
    memory.release(&parents);
    Ok(used)
}
fn add_link(slots: &mut [Option<usize>; 2], edge: usize) {
    // The preceding use-count gate rejects a third edge use before either
    // endpoint is updated. Exactly one link per endpoint is inserted per use.
    slots[usize::from(slots[0].is_some())] = Some(edge);
}
fn root(parents: &mut [usize], mut x: usize, w: &mut Work<'_>) -> Result<usize, Error> {
    while parents[x] != x {
        w.charge(1)?;
        x = parents[x];
    }
    Ok(x)
}

fn cap_forward(
    actual: [VertexOwner; 2],
    original: [BoundaryVertexKey; 2],
    forward: bool,
) -> Result<bool, Error> {
    if actual == original.map(VertexOwner::Original) {
        Ok(forward)
    } else if [actual[1], actual[0]] == original.map(VertexOwner::Original) {
        Ok(!forward)
    } else {
        Err(Error::Unresolved)
    }
}
