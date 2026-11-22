//! Incidence-preserving union of exact translated thread charts. The old
//! cell edge arena supplies ownership; equal coordinates never merge nodes.
use super::*;
use std::collections::{BTreeMap, BTreeSet};
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ChartFaceKind {
    Thread(usize),
    PlainWall { quarter: u8 },
    Plane(FaceKey),
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ChartEdgeKind {
    Original(usize),
    Helix { source_edge: usize },
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ChartCoedge {
    pub edge: usize,
    pub sense: bool,
}
pub struct ChartEdge {
    pub vertices: [usize; 2],
    pub kind: ChartEdgeKind,
}
pub struct ChartFace {
    pub kind: ChartFaceKind,
    pub same_sense: bool,
    loops: Vec<Vec<ChartCoedge>>,
}
impl ChartFace {
    pub fn loops(&self) -> &[Vec<ChartCoedge>] {
        &self.loops
    }
}
pub struct CompactBoundary<'c, 'd, 'f, 's> {
    source: &'c ClosedField<'d, 'f, 's>,
    surfaces: &'c CompactSurfaces<'f, 's>,
    faces: Vec<ChartFace>,
    edges: Vec<ChartEdge>,
    edge_map: Vec<Option<ChartCoedge>>,
    face_map: Vec<usize>,
    work: usize,
    bytes: usize,
}
impl<'c, 'd, 'f, 's> CompactBoundary<'c, 'd, 'f, 's> {
    pub fn source(&self) -> &ClosedField<'d, 'f, 's> {
        self.source
    }
    pub fn surfaces(&self) -> &CompactSurfaces<'f, 's> {
        self.surfaces
    }
    pub fn faces(&self) -> &[ChartFace] {
        &self.faces
    }
    pub fn edges(&self) -> &[ChartEdge] {
        &self.edges
    }
    pub fn original_edges(&self) -> &[Option<ChartCoedge>] {
        &self.edge_map
    }
    pub fn original_faces(&self) -> &[usize] {
        &self.face_map
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_bytes(&self) -> usize {
        self.bytes
    }
    pub fn prepare(
        source: &'c ClosedField<'d, 'f, 's>,
        surfaces: &'c CompactSurfaces<'f, 's>,
        options: DomainOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        if !std::ptr::eq(source.domains().field(), surfaces.field()) {
            return Err(Error::Invalid);
        }
        let inherited = source.retained_bytes() + surfaces.retained_bytes();
        let (mut w, mut m) = crate::step::hanger::domain::ledgers(options, inherited, interrupt)?;
        w.charge(source.work_used() + surfaces.work_used() - surfaces.field().work_used())?;
        let base = m.used;
        let n = source.edges().len();
        // Pinned std BTree capacity B=6: <=11 key/value slots plus12 child
        // pointers per node. These keys/values need <1024 bytes per node.
        // <=10N live map/set nodes and <2048N bytes of transient/final Vec
        // capacity cover this builder's simultaneous arenas, clones and chains.
        let reserve = n.checked_mul(12 * 1024).ok_or(Error::Limit)?;
        m.add(reserve)?;
        let tree_cost = 64 * (1 + n.ilog2() as usize);
        let (faces, edges, edge_map, face_map) = build(source, surfaces, &mut w, tree_cost)?;
        let owned = size_of::<Self>()
            + faces.capacity() * size_of::<ChartFace>()
            + edges.capacity() * size_of::<ChartEdge>()
            + edge_map.capacity() * size_of::<Option<ChartCoedge>>()
            + face_map.capacity() * size_of::<usize>()
            + faces
                .iter()
                .map(|f| {
                    f.loops.capacity() * size_of::<Vec<ChartCoedge>>()
                        + f.loops
                            .iter()
                            .map(|l| l.capacity() * size_of::<ChartCoedge>())
                            .sum::<usize>()
                })
                .sum::<usize>();
        m.used = base;
        m.add(owned)?;
        w.charge(1)?;
        Ok(Self {
            source,
            surfaces,
            faces,
            edges,
            edge_map,
            face_map,
            work: w.used,
            bytes: owned,
        })
    }
}
type Packet = (
    Vec<ChartFace>,
    Vec<ChartEdge>,
    Vec<Option<ChartCoedge>>,
    Vec<usize>,
);
fn build(
    source: &ClosedField<'_, '_, '_>,
    surfaces: &CompactSurfaces<'_, '_>,
    w: &mut Work<'_>,
    cost: usize,
) -> Result<Packet, Error> {
    let mut face_map = Vec::with_capacity(source.faces().len());
    let mut kinds = (0..surfaces.surfaces().len())
        .map(ChartFaceKind::Thread)
        .collect::<Vec<_>>();
    for face in source.faces() {
        w.charge(1)?;
        face_map.push(match face.kind {
            ClosedFaceKind::ThreadCell(i) => {
                surfaces.cells()[source.domains().cells()[i].source_cell].surface()
            }
            ClosedFaceKind::Plane(key) => {
                let id = kinds.len();
                kinds.push(ChartFaceKind::Plane(key));
                id
            }
            ClosedFaceKind::PlainWall { quarter } => {
                let id = kinds.len();
                kinds.push(ChartFaceKind::PlainWall { quarter });
                id
            }
        })
    }
    let mut old_uses = vec![Vec::new(); source.edges().len()];
    let mut face_senses = vec![None; kinds.len()];
    for (i, face) in source.faces().iter().enumerate() {
        w.charge(1)?;
        face_senses[face_map[i]] = Some(face.same_sense);
        for c in face.loops().iter().flatten() {
            w.charge(1)?;
            old_uses[c.edge].push(i)
        }
    }
    let mut groups: BTreeMap<[usize; 2], Vec<usize>> = BTreeMap::new();
    let mut edge_map = vec![None; source.edges().len()];
    let mut edges = Vec::new();
    for (i, edge) in source.edges().iter().enumerate() {
        w.charge(1)?;
        let pair = [face_map[old_uses[i][0]], face_map[old_uses[i][1]]];
        if pair[0] == pair[1] {
            continue;
        }
        let a = old_uses[i][0];
        let b = old_uses[i][1];
        if matches!(source.faces()[a].kind, ClosedFaceKind::ThreadCell(_))
            && matches!(source.faces()[b].kind, ClosedFaceKind::ThreadCell(_))
        {
            w.charge(cost)?;
            groups
                .entry([pair[0].min(pair[1]), pair[0].max(pair[1])])
                .or_default()
                .push(i)
        } else {
            edge_map[i] = Some(ChartCoedge {
                edge: edges.len(),
                sense: true,
            });
            edges.push(ChartEdge {
                vertices: edge.vertices,
                kind: ChartEdgeKind::Original(i),
            })
        }
    }
    for members in groups.values() {
        let mut degrees: BTreeMap<usize, usize> = BTreeMap::new();
        for &i in members {
            for v in source.edges()[i].vertices {
                w.charge(cost)?;
                *degrees.entry(v).or_default() += 1
            }
        }
        w.charge(degrees.len())?;
        let ends = degrees
            .iter()
            .filter(|(_, n)| **n == 1)
            .map(|(&v, _)| v)
            .collect::<Vec<_>>();
        if ends.len() != 2 || degrees.values().any(|&n| n != 1 && n != 2) {
            return Err(Error::Unresolved);
        }
        let target = edges.len();
        edges.push(ChartEdge {
            vertices: [ends[0], ends[1]],
            kind: ChartEdgeKind::Helix {
                source_edge: members[0],
            },
        });
        let mut current = ends[0];
        w.charge(members.len() * cost)?;
        let mut unseen = members.iter().copied().collect::<BTreeSet<_>>();
        while current != ends[1] {
            let mut found = None;
            for &i in &unseen {
                w.charge(8)?;
                if source.edges()[i].vertices.contains(&current) {
                    found = Some(i);
                    break;
                }
            }
            let edge = found.ok_or(Error::Unresolved)?;
            w.charge(cost)?;
            unseen.remove(&edge);
            let old = source.edges()[edge].vertices;
            let forward = old[0] == current;
            edge_map[edge] = Some(ChartCoedge {
                edge: target,
                sense: forward,
            });
            current = if forward { old[1] } else { old[0] };
        }
        if !unseen.is_empty() {
            return Err(Error::Unresolved);
        }
    }
    let mut face_uses = vec![BTreeMap::new(); kinds.len()];
    for (i, face) in source.faces().iter().enumerate() {
        for c in face.loops().iter().flatten() {
            w.charge(1)?;
            if let Some(mapped) = edge_map[c.edge] {
                let sense = c.sense == mapped.sense;
                w.charge(cost)?;
                if let Some(previous) = face_uses[face_map[i]].insert(mapped.edge, sense)
                    && previous != sense
                {
                    return Err(Error::Unresolved);
                }
            }
        }
    }
    let directed = |e: usize, sense: bool| {
        if sense {
            edges[e].vertices
        } else {
            [edges[e].vertices[1], edges[e].vertices[0]]
        }
    };
    let mut faces = Vec::new();
    for (face_id, uses) in face_uses.iter().enumerate() {
        w.charge(64 * uses.len())?;
        let mut remaining = uses.clone();
        let mut loops = Vec::new();
        while let Some((&first, &sense)) = remaining.first_key_value() {
            w.charge(cost)?;
            remaining.remove(&first);
            let [start, mut current] = directed(first, sense);
            let mut cycle = vec![ChartCoedge { edge: first, sense }];
            while current != start {
                let mut found = None;
                for (&e, &s) in &remaining {
                    w.charge(8)?;
                    if directed(e, s)[0] == current {
                        found = Some((e, s));
                        break;
                    }
                }
                let (next, sense) = found.ok_or(Error::Unresolved)?;
                w.charge(cost)?;
                remaining.remove(&next);
                let ends = directed(next, sense);
                cycle.push(ChartCoedge { edge: next, sense });
                current = ends[1];
            }
            loops.push(cycle);
        }
        // Cycle discovery follows edge IDs, not geometric nesting. Shoulder
        // annuli have a unique outer owner because source admission proves
        // metric radius < plain radius < wood radius. Preserve that ownership
        // when choosing the STEP outer bound; never infer it from ID ordering.
        if loops.len() == 2 {
            let outer_plain = match kinds[face_id] {
                ChartFaceKind::Plane(FaceKey::MetricShoulder) => true,
                ChartFaceKind::Plane(FaceKey::WoodShoulder) => false,
                _ => return Err(Error::Unresolved),
            };
            let mut outer = None;
            for (ring, cycle) in loops.iter().enumerate() {
                let edge = cycle.first().ok_or(Error::Unresolved)?.edge;
                for (other, uses) in face_uses.iter().enumerate() {
                    w.charge(cost)?;
                    if other != face_id && uses.contains_key(&edge) {
                        let plain = matches!(kinds[other], ChartFaceKind::PlainWall { .. });
                        if plain == outer_plain && outer.replace(ring).is_some() {
                            return Err(Error::Unresolved);
                        }
                    }
                }
            }
            loops.swap(0, outer.ok_or(Error::Unresolved)?);
        }
        faces.push(ChartFace {
            kind: kinds[face_id],
            same_sense: face_senses[face_id].ok_or(Error::Unresolved)?,
            loops,
        });
    }
    Ok((faces, edges, edge_map, face_map))
}
#[cfg(test)]
mod tests;
