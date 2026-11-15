//! Cap subdivision by exact event ownership. Coordinates never weld vertices.
use super::*;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TrimVertex {
    Graph(usize),
    CapEvent(usize),
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TrimEdgeKind {
    Original(usize),
    CapSegment(usize),
}
pub struct TrimEdge {
    pub vertices: [usize; 2],
    pub kind: TrimEdgeKind,
}
pub struct TrimFace {
    pub kind: ChartFaceKind,
    pub same_sense: bool,
    loops: Vec<Vec<ChartCoedge>>,
}
impl TrimFace {
    pub fn loops(&self) -> &[Vec<ChartCoedge>] {
        &self.loops
    }
}
pub struct CapTopology<'b, 'g, 'v, 'c, 'd, 'f, 's> {
    source: &'b CompactBoundary<'c, 'd, 'f, 's>,
    events: &'g CapEvents<'v, 'c, 'd, 'f, 's>,
    vertices: Vec<TrimVertex>,
    edges: Vec<TrimEdge>,
    faces: Vec<TrimFace>,
    work: usize,
    bytes: usize,
}
impl CapTopology<'_, '_, '_, '_, '_, '_, '_> {
    pub fn source(&self) -> &CompactBoundary<'_, '_, '_, '_> {
        self.source
    }
    pub fn events(&self) -> &CapEvents<'_, '_, '_, '_, '_> {
        self.events
    }
    pub fn vertices(&self) -> &[TrimVertex] {
        &self.vertices
    }
    pub fn edges(&self) -> &[TrimEdge] {
        &self.edges
    }
    pub fn faces(&self) -> &[TrimFace] {
        &self.faces
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_bytes(&self) -> usize {
        self.bytes
    }
    pub fn coordinates(&self, vertex: usize) -> Option<&[EncodedScalar; 3]> {
        match self.vertices.get(vertex)? {
            TrimVertex::Graph(v) => Some(self.events.source().vertices()[*v].coordinates()),
            TrimVertex::CapEvent(v) => self.events.coordinates(*v),
        }
    }
}
impl<'c, 'd, 'f, 's> CompactBoundary<'c, 'd, 'f, 's> {
    pub fn trim_caps<'b, 'g, 'v>(
        &'b self,
        events: &'g CapEvents<'v, 'c, 'd, 'f, 's>,
        options: DomainOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<CapTopology<'b, 'g, 'v, 'c, 'd, 'f, 's>, Error> {
        if !std::ptr::eq(self.source(), events.source().source()) {
            return Err(Error::Invalid);
        }
        let inherited =
            self.retained_bytes() + self.surfaces().retained_bytes() + events.retained_bytes();
        let (mut w, mut m) = crate::step::hanger::domain::ledgers(options, inherited, interrupt)?;
        w.charge(self.work_used() + events.work_used() - self.source().work_used())?;
        // All arenas below have linear bounded capacity. Reserve their complete
        // simultaneous capacity before allocation, including per-face loops and
        // per-original-edge chains and the temporary event-to-vertex lookup.
        let count = self.source().vertices().len()
            + events.events().len()
            + events.segments().len()
            + self.edges().len()
            + self.faces().len();
        let reserve = count.checked_mul(2048).ok_or(Error::Limit)?;
        m.add(reserve)?;
        let mut vertices: Vec<_> = (0..self.source().vertices().len())
            .map(TrimVertex::Graph)
            .collect();
        let mut event_vertices = Vec::with_capacity(events.events().len());
        for (id, event) in events.events().iter().enumerate() {
            w.charge(1)?;
            event_vertices.push(match event.graph_vertex() {
                Some(v) => v,
                None => {
                    let v = vertices.len();
                    vertices.push(TrimVertex::CapEvent(id));
                    v
                }
            });
        }
        let mut edges = Vec::new();
        let mut chains = Vec::with_capacity(self.edges().len());
        for (old, edge) in self.edges().iter().enumerate() {
            let mut pieces = Vec::new();
            if let ChartEdgeKind::Original(original) = edge.kind {
                for (segment, s) in events.segments().iter().enumerate() {
                    w.charge(1)?;
                    if s.original_edge() == original {
                        pieces.push((segment, s.event_ids().map(|e| event_vertices[e])));
                    }
                }
            }
            let mut chain = Vec::new();
            if pieces.is_empty() {
                chain.push(ChartCoedge {
                    edge: edges.len(),
                    sense: true,
                });
                edges.push(TrimEdge {
                    vertices: edge.vertices,
                    kind: TrimEdgeKind::Original(old),
                });
            } else {
                let mut current = edge.vertices[0];
                while !pieces.is_empty() {
                    let mut next = None;
                    for (i, (_, v)) in pieces.iter().enumerate() {
                        w.charge(1)?;
                        if v.contains(&current) {
                            if next.is_some() {
                                return Err(Error::Unresolved);
                            }
                            next = Some(i);
                        }
                    }
                    w.charge(pieces.len() + 1)?;
                    let (segment, v) = pieces.remove(next.ok_or(Error::Unresolved)?);
                    let sense = v[0] == current;
                    chain.push(ChartCoedge {
                        edge: edges.len(),
                        sense,
                    });
                    edges.push(TrimEdge {
                        vertices: v,
                        kind: TrimEdgeKind::CapSegment(segment),
                    });
                    current = if sense { v[1] } else { v[0] };
                }
                if current != edge.vertices[1] {
                    return Err(Error::Unresolved);
                }
            }
            chains.push(chain);
        }
        let mut faces = Vec::with_capacity(self.faces().len());
        for face in self.faces() {
            let mut loops = Vec::with_capacity(face.loops().len());
            for original in face.loops() {
                let mut output = Vec::new();
                for c in original {
                    w.charge(chains[c.edge].len())?;
                    if c.sense {
                        output.extend_from_slice(&chains[c.edge])
                    } else {
                        output.extend(chains[c.edge].iter().rev().map(|c| ChartCoedge {
                            edge: c.edge,
                            sense: !c.sense,
                        }))
                    }
                }
                loops.push(output);
            }
            faces.push(TrimFace {
                kind: face.kind,
                same_sense: face.same_sense,
                loops,
            });
        }
        let bytes = size_of::<CapTopology<'_, '_, '_, '_, '_, '_, '_>>()
            + vertices.capacity() * size_of::<TrimVertex>()
            + edges.capacity() * size_of::<TrimEdge>()
            + faces.capacity() * size_of::<TrimFace>()
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
        if bytes > reserve {
            return Err(Error::Limit);
        }
        w.charge(count)?;
        Ok(CapTopology {
            source: self,
            events,
            vertices,
            edges,
            faces,
            work: w.used,
            bytes,
        })
    }
}
#[cfg(test)]
mod tests;
