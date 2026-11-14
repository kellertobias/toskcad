//! Closed finite field INCIDENCE, with true disk/annular plane domains.
//! This is not a serialized or independently decoded occupied-body receipt.
use super::*;
mod measure;
pub use measure::FieldVolume;
use std::cmp::Ordering;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ClosedVertex {
    Thread(usize),
    Plain { end: End, quarter: u8 },
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ClosedCurve {
    Thread(usize),
    PlainCircle { end: End, quarter: u8 },
    PlainRail { quarter: u8 },
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ClosedEdge {
    pub vertices: [usize; 2],
    pub curve: ClosedCurve,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ClosedCoedge {
    pub edge: usize,
    pub sense: bool,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ClosedFaceKind {
    ThreadCell(usize),
    PlainWall { quarter: u8 },
    Plane(FaceKey),
}
pub struct ClosedFace {
    pub kind: ClosedFaceKind,
    pub same_sense: bool,
    loops: Vec<Vec<ClosedCoedge>>,
}
impl ClosedFace {
    pub fn loops(&self) -> &[Vec<ClosedCoedge>] {
        &self.loops
    }
}
pub struct ClosedField<'d, 'f, 's> {
    domains: &'d FiniteDomains<'f, 's>,
    vertices: Vec<ClosedVertex>,
    edges: Vec<ClosedEdge>,
    faces: Vec<ClosedFace>,
    work: usize,
    bytes: usize,
}
fn ends(edge: &ClosedEdge, sense: bool) -> [usize; 2] {
    if sense {
        edge.vertices
    } else {
        [edge.vertices[1], edge.vertices[0]]
    }
}
fn station(v: &FiniteVertex, field: &EncodedField<'_>, w: &mut Work<'_>) -> Result<N, Error> {
    let s = field.source.source_boundary().source();
    let (pitch, n) = if v.component == Component::Metric {
        let x = s.metric().spec();
        (x.pitch_mm, x.starts)
    } else {
        let x = s.wood().spec();
        (x.pitch_mm, x.starts)
    };
    let pitch = binary(pitch, w)?;
    v.phase
        .number(w)?
        .add(&pitch.mul(&N::integer(i64::from(v.period)), w)?, w)?
        .add(
            &pitch
                .mul(&N::integer(i64::from(n)), w)?
                .mul(&v.turn.number(w)?, w)?,
            w,
        )
}
fn boundary(
    domains: &FiniteDomains<'_, '_>,
    edge: &FiniteEdge,
    w: &mut Work<'_>,
) -> Result<(Component, End), Error> {
    let a = &domains.vertices()[edge.vertices[0]];
    let b = &domains.vertices()[edge.vertices[1]];
    let sa = station(a, domains.field(), w)?;
    let sb = station(b, domains.field(), w)?;
    if sa.cmp(&N::integer(0), w)? == Ordering::Equal
        && sb.cmp(&N::integer(0), w)? == Ordering::Equal
    {
        return Ok((a.component, End::Start));
    }
    let s = domains.field().source.source_boundary().source();
    let length = if a.component == Component::Metric {
        s.metric().full_length_mm()
    } else {
        s.wood().full_length_mm()
    };
    let l = binary(length, w)?;
    if sa.cmp(&l, w)? != Ordering::Equal || sb.cmp(&l, w)? != Ordering::Equal {
        return Err(Error::Unresolved);
    }
    Ok((a.component, End::End))
}
fn ordered(
    mut pending: Vec<ClosedCoedge>,
    edges: &[ClosedEdge],
    m: &mut Memory,
    w: &mut Work<'_>,
) -> Result<Vec<ClosedCoedge>, Error> {
    if pending.is_empty() {
        return Err(Error::Unresolved);
    }
    let mut out = m.allocate(pending.len(), w)?;
    w.charge(pending.len())?;
    let first = pending.remove(0);
    let start = ends(&edges[first.edge], first.sense)[0];
    let mut tail = ends(&edges[first.edge], first.sense)[1];
    out.push(first);
    while !pending.is_empty() {
        if tail == start {
            return Err(Error::Unresolved);
        }
        let mut next = None;
        for (i, c) in pending.iter().enumerate() {
            w.charge(1)?;
            if ends(&edges[c.edge], c.sense)[0] == tail {
                next = Some(i);
                break;
            }
        }
        let i = next.ok_or(Error::Unresolved)?;
        w.charge(pending.len())?;
        let c = pending.remove(i);
        tail = ends(&edges[c.edge], c.sense)[1];
        out.push(c);
    }
    if tail != start {
        return Err(Error::Unresolved);
    }
    m.release(&pending);
    Ok(out)
}
fn cap_ring(
    domains: &FiniteDomains<'_, '_>,
    edges: &[ClosedEdge],
    component: Component,
    end: End,
    m: &mut Memory,
    w: &mut Work<'_>,
) -> Result<Vec<ClosedCoedge>, Error> {
    let source = domains.field().source.source_boundary().source();
    let left = if component == Component::Metric {
        source.metric().spec().left_hand
    } else {
        source.wood().spec().left_hand
    };
    let mut pending = m.allocate(domains.edges().len(), w)?;
    for (i, e) in domains.edges().iter().enumerate() {
        w.charge(1)?;
        if e.faces[1].is_some() {
            continue;
        }
        if boundary(domains, e, w)? == (component, end) {
            let (_, uvsense) = e.faces[0].ok_or(Error::Unresolved)?;
            pending.push(ClosedCoedge {
                edge: i,
                sense: !(uvsense ^ !left),
            });
        }
    }
    ordered(pending, edges, m, w)
}
impl<'d, 'f, 's> ClosedField<'d, 'f, 's> {
    pub fn prepare(
        domains: &'d FiniteDomains<'f, 's>,
        options: DomainOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        let (mut w, mut memory) = super::super::domain::ledgers(
            options,
            domains.retained_bytes() + size_of::<Self>(),
            interrupt,
        )?;
        w.charge(domains.work_used())?;
        memory.add(SCRATCH - super::super::SCRATCH)?;
        let field = domains.field();
        let source = field.source.source_boundary().source();
        let plain = source.spec().plain_radius_mm;
        let mut metric_error = 0_f64;
        let mut wood_error = 0_f64;
        for meridian in field.meridians() {
            w.charge(1)?;
            let e = meridian.bounds.position_error_mm;
            if meridian.component == Component::Metric {
                metric_error = metric_error.max(e);
            } else {
                wood_error = wood_error.max(e);
            }
        }
        // Common actual B(w) cancels in the same-angle radial ordering.
        // These complete source/proxy bounds are conservative and directed.
        if (p(source.metric().maximum_radius_mm()) + p(metric_error)).hi >= plain
            || (p(source.wood().minimum_radius_mm()) - p(wood_error)).lo <= plain
        {
            return Err(Error::Unresolved);
        }
        let mut vertices = memory.allocate(domains.vertices().len() + 8, &mut w)?;
        for i in 0..domains.vertices().len() {
            w.charge(1)?;
            vertices.push(ClosedVertex::Thread(i));
        }
        for end in [End::Start, End::End] {
            for quarter in 0..4 {
                vertices.push(ClosedVertex::Plain { end, quarter });
            }
        }
        let n = domains.vertices().len();
        let e = domains.edges().len();
        let mut edges = memory.allocate(e + 12, &mut w)?;
        for (i, edge) in domains.edges().iter().enumerate() {
            w.charge(1)?;
            edges.push(ClosedEdge {
                vertices: edge.vertices,
                curve: ClosedCurve::Thread(i),
            });
        }
        for end in [End::Start, End::End] {
            let offset = usize::from(end == End::End) * 4;
            for quarter in 0..4 {
                edges.push(ClosedEdge {
                    vertices: [n + offset + quarter, n + offset + (quarter + 1) % 4],
                    curve: ClosedCurve::PlainCircle {
                        end,
                        quarter: quarter as u8,
                    },
                });
            }
        }
        for quarter in 0..4 {
            edges.push(ClosedEdge {
                vertices: [n + quarter, n + 4 + quarter],
                curve: ClosedCurve::PlainRail {
                    quarter: quarter as u8,
                },
            });
        }
        let mut faces = memory.allocate(domains.cells().len() + 8, &mut w)?;
        for (index, cell) in domains.cells().iter().enumerate() {
            let actual = &field.cells[cell.source_cell];
            let profile = &field.meridians[actual.meridian];
            let left = if profile.component == Component::Metric {
                source.metric().spec().left_hand
            } else {
                source.wood().spec().left_hand
            };
            let mut wire = memory.allocate(cell.corners().len(), &mut w)?;
            for i in 0..cell.corners().len() {
                let a = cell.corners()[i].vertex;
                let b = cell.corners()[(i + 1) % cell.corners().len()].vertex;
                let pair = [a.min(b), a.max(b)];
                let mut edge = None;
                for (j, e) in domains.edges().iter().enumerate() {
                    w.charge(1)?;
                    if e.vertices == pair {
                        edge = Some(j);
                        break;
                    }
                }
                wire.push(ClosedCoedge {
                    edge: edge.ok_or(Error::Unresolved)?,
                    sense: (a < b) ^ !left,
                });
            }
            if !left {
                wire.reverse();
            }
            let mut loops = memory.allocate(1, &mut w)?;
            loops.push(wire);
            faces.push(ClosedFace {
                kind: ClosedFaceKind::ThreadCell(index),
                same_sense: left,
                loops,
            });
        }
        for quarter in 0..4 {
            let mut wire = memory.allocate(4, &mut w)?;
            wire.extend_from_slice(&[
                ClosedCoedge {
                    edge: e + quarter,
                    sense: true,
                },
                ClosedCoedge {
                    edge: e + 8 + (quarter + 1) % 4,
                    sense: true,
                },
                ClosedCoedge {
                    edge: e + 4 + quarter,
                    sense: false,
                },
                ClosedCoedge {
                    edge: e + 8 + quarter,
                    sense: false,
                },
            ]);
            let mut loops = memory.allocate(1, &mut w)?;
            loops.push(wire);
            faces.push(ClosedFace {
                kind: ClosedFaceKind::PlainWall {
                    quarter: quarter as u8,
                },
                same_sense: false,
                loops,
            });
        }
        for role in [
            FaceKey::MetricTip,
            FaceKey::MetricShoulder,
            FaceKey::WoodShoulder,
            FaceKey::WoodTip,
        ] {
            let mut loops = memory.allocate(
                if matches!(role, FaceKey::MetricShoulder | FaceKey::WoodShoulder) {
                    2
                } else {
                    1
                },
                &mut w,
            )?;
            if role == FaceKey::MetricShoulder {
                let mut circle = memory.allocate(4, &mut w)?;
                for q in (0..4).rev() {
                    circle.push(ClosedCoedge {
                        edge: e + q,
                        sense: false,
                    });
                }
                loops.push(circle);
            }
            let (component, end) = match role {
                FaceKey::MetricTip => (Component::Metric, End::Start),
                FaceKey::MetricShoulder => (Component::Metric, End::End),
                FaceKey::WoodShoulder => (Component::Wood, End::Start),
                _ => (Component::Wood, End::End),
            };
            loops.push(cap_ring(
                domains,
                &edges,
                component,
                end,
                &mut memory,
                &mut w,
            )?);
            if role == FaceKey::WoodShoulder {
                let mut circle = memory.allocate(4, &mut w)?;
                for q in 0..4 {
                    circle.push(ClosedCoedge {
                        edge: e + 4 + q,
                        sense: true,
                    });
                }
                loops.push(circle);
            }
            faces.push(ClosedFace {
                kind: ClosedFaceKind::Plane(role),
                same_sense: role == FaceKey::WoodTip,
                loops,
            });
        }
        let result = Self {
            domains,
            vertices,
            edges,
            faces,
            work: 0,
            bytes: 0,
        };
        result.certify(&mut memory, &mut w)?;
        w.charge(1)?;
        memory.used -= SCRATCH;
        Ok(Self {
            work: w.used,
            bytes: memory.used,
            ..result
        })
    }
    fn certify(&self, m: &mut Memory, w: &mut Work<'_>) -> Result<(), Error> {
        if self.vertices.len() != self.domains.vertices().len() + 8
            || self.edges.len() != self.domains.edges().len() + 12
            || self.faces.len() != self.domains.cells().len() + 8
        {
            return Err(Error::Unresolved);
        }
        // Bind every produced reference to its immutable source geometry.
        // Vertex order may change if all actual references change with it.
        for (i, v) in self.vertices.iter().enumerate() {
            w.charge(1)?;
            match *v {
                ClosedVertex::Thread(j) => {
                    if j >= self.domains.vertices().len() {
                        return Err(Error::Unresolved);
                    }
                }
                ClosedVertex::Plain { quarter, .. } => {
                    if quarter >= 4 {
                        return Err(Error::Unresolved);
                    }
                }
            }
            for prior in &self.vertices[..i] {
                w.charge(1)?;
                if prior == v {
                    return Err(Error::Unresolved);
                }
            }
        }
        for (index, edge) in self.edges.iter().enumerate() {
            for previous in &self.edges[..index] {
                w.charge(1)?;
                if previous.curve == edge.curve {
                    return Err(Error::Unresolved);
                }
            }
            w.charge(1)?;
            let actual = edge.vertices.map(|i| self.vertices.get(i).copied());
            let expected = match edge.curve {
                ClosedCurve::Thread(i) => self
                    .domains
                    .edges()
                    .get(i)
                    .map(|e| e.vertices.map(ClosedVertex::Thread)),
                ClosedCurve::PlainCircle { end, quarter } if quarter < 4 => Some([
                    ClosedVertex::Plain { end, quarter },
                    ClosedVertex::Plain {
                        end,
                        quarter: (quarter + 1) % 4,
                    },
                ]),
                ClosedCurve::PlainRail { quarter } if quarter < 4 => Some([
                    ClosedVertex::Plain {
                        end: End::Start,
                        quarter,
                    },
                    ClosedVertex::Plain {
                        end: End::End,
                        quarter,
                    },
                ]),
                _ => None,
            }
            .ok_or(Error::Unresolved)?;
            if actual != expected.map(Some) {
                return Err(Error::Unresolved);
            }
        }
        let mut incidence = m.allocate(self.edges.len(), w)?;
        incidence.resize(self.edges.len(), [None; 2]);
        let mut adjacency = m.allocate(self.edges.len(), w)?;
        let mut holes = 0;
        for (face, f) in self.faces.iter().enumerate() {
            if f.loops.is_empty() {
                return Err(Error::Unresolved);
            }
            holes += f.loops.len() - 1;
            for wire in &f.loops {
                if wire.len() < 3 {
                    return Err(Error::Unresolved);
                }
                for i in 0..wire.len() {
                    w.charge(1)?;
                    let c = wire[i];
                    let edge = self.edges.get(c.edge).ok_or(Error::Unresolved)?;
                    let next = wire[(i + 1) % wire.len()];
                    if ends(edge, c.sense)[1]
                        != ends(
                            self.edges.get(next.edge).ok_or(Error::Unresolved)?,
                            next.sense,
                        )[0]
                    {
                        return Err(Error::Unresolved);
                    }
                    let use_ = &mut incidence[c.edge];
                    if let Some((first, sense)) = use_[0] {
                        if use_[1].is_some() || sense == c.sense {
                            return Err(Error::Unresolved);
                        }
                        use_[1] = Some((face, c.sense));
                        adjacency.push([first, face]);
                    } else {
                        use_[0] = Some((face, c.sense));
                    }
                }
            }
        }
        if incidence.iter().any(|v| v[1].is_none())
            || self.vertices.len() + self.faces.len() != self.edges.len() + holes + 2
        {
            return Err(Error::Unresolved);
        }
        let mut seen = m.allocate(self.faces.len(), w)?;
        seen.resize(self.faces.len(), false);
        seen[0] = true;
        for _ in 0..self.faces.len() {
            let mut changed = false;
            for &[a, b] in &adjacency {
                // Recorded exactly when an edge earns its second opposite
                // coface. The completeness guard above proves one/edge.
                w.charge(1)?;
                if seen[a] != seen[b] {
                    seen[a] = true;
                    seen[b] = true;
                    changed = true;
                }
            }
            if !changed {
                break;
            }
        }
        if seen.iter().any(|x| !*x) {
            return Err(Error::Unresolved);
        }
        m.release(&seen);
        m.release(&adjacency);
        m.release(&incidence);
        self.vertex_links(m, w)?;
        self.bind_faces(w)
    }
    fn bind_faces(&self, w: &mut Work<'_>) -> Result<(), Error> {
        let field = self.domains.field();
        let source = field.source.source_boundary().source();
        for (index, face) in self.faces.iter().enumerate() {
            for prior in &self.faces[..index] {
                w.charge(1)?;
                if prior.kind == face.kind {
                    return Err(Error::Unresolved);
                }
            }
            match face.kind {
                ClosedFaceKind::ThreadCell(cell) => {
                    let original = self.domains.cells().get(cell).ok_or(Error::Unresolved)?;
                    let key = &field.cells[original.source_cell];
                    let profile = &field.meridians[key.meridian];
                    let left = if profile.component == Component::Metric {
                        source.metric().spec().left_hand
                    } else {
                        source.wood().spec().left_hand
                    };
                    if face.same_sense != left
                        || face.loops.len() != 1
                        || face.loops[0].len() != original.corners().len()
                    {
                        return Err(Error::Unresolved);
                    }
                    for c in &face.loops[0] {
                        w.charge(original.corners().len())?;
                        let edge = self.edges.get(c.edge).ok_or(Error::Unresolved)?;
                        let actual = ends(edge, c.sense).map(|v| self.vertices.get(v).copied());
                        let valid = (0..original.corners().len()).any(|i| {
                            let pair = [
                                ClosedVertex::Thread(original.corners()[i].vertex),
                                ClosedVertex::Thread(
                                    original.corners()[(i + 1) % original.corners().len()].vertex,
                                ),
                            ];
                            actual
                                == if left {
                                    pair.map(Some)
                                } else {
                                    [Some(pair[1]), Some(pair[0])]
                                }
                        });
                        if !valid {
                            return Err(Error::Unresolved);
                        }
                    }
                }
                ClosedFaceKind::PlainWall { quarter } => {
                    if quarter >= 4
                        || face.same_sense
                        || face.loops.len() != 1
                        || face.loops[0].len() != 4
                    {
                        return Err(Error::Unresolved);
                    }
                    let q = quarter;
                    let next = (q + 1) % 4;
                    let expected = [
                        (
                            ClosedCurve::PlainCircle {
                                end: End::Start,
                                quarter: q,
                            },
                            true,
                        ),
                        (ClosedCurve::PlainRail { quarter: next }, true),
                        (
                            ClosedCurve::PlainCircle {
                                end: End::End,
                                quarter: q,
                            },
                            false,
                        ),
                        (ClosedCurve::PlainRail { quarter: q }, false),
                    ];
                    for c in &face.loops[0] {
                        w.charge(4)?;
                        let e = self.edges.get(c.edge).ok_or(Error::Unresolved)?;
                        if !expected.contains(&(e.curve, c.sense)) {
                            return Err(Error::Unresolved);
                        }
                    }
                }
                ClosedFaceKind::Plane(role) => {
                    let count = match role {
                        FaceKey::MetricTip | FaceKey::WoodTip => 1,
                        FaceKey::MetricShoulder | FaceKey::WoodShoulder => 2,
                        _ => return Err(Error::Unresolved),
                    };
                    if face.loops.len() != count || face.same_sense != (role == FaceKey::WoodTip) {
                        return Err(Error::Unresolved);
                    }
                    for (loop_index, wire) in face.loops.iter().enumerate() {
                        for c in wire {
                            w.charge(1)?;
                            let edge = self.edges.get(c.edge).ok_or(Error::Unresolved)?;
                            let valid = match edge.curve {
                                ClosedCurve::Thread(i) => {
                                    let side =
                                        self.domains.edges().get(i).ok_or(Error::Unresolved)?;
                                    if side.faces[1].is_some() {
                                        false
                                    } else {
                                        let target = boundary(self.domains, side, w)?;
                                        match (role, loop_index) {
                                            (FaceKey::MetricTip, 0) => {
                                                target == (Component::Metric, End::Start)
                                            }
                                            (FaceKey::MetricShoulder, 1) => {
                                                target == (Component::Metric, End::End)
                                            }
                                            (FaceKey::WoodShoulder, 0) => {
                                                target == (Component::Wood, End::Start)
                                            }
                                            (FaceKey::WoodTip, 0) => {
                                                target == (Component::Wood, End::End)
                                            }
                                            _ => false,
                                        }
                                    }
                                }
                                ClosedCurve::PlainCircle { end, .. } => matches!(
                                    (role, loop_index, end),
                                    (FaceKey::MetricShoulder, 0, End::Start)
                                        | (FaceKey::WoodShoulder, 1, End::End)
                                ),
                                _ => false,
                            };
                            if !valid {
                                return Err(Error::Unresolved);
                            }
                        }
                    }
                }
            }
        }
        Ok(())
    }
    fn vertex_links(&self, m: &mut Memory, w: &mut Work<'_>) -> Result<(), Error> {
        for vertex in 0..self.vertices.len() {
            let mut links = m.allocate(self.faces.len(), w)?;
            for face in &self.faces {
                for wire in &face.loops {
                    for i in 0..wire.len() {
                        w.charge(1)?;
                        let here = ends(&self.edges[wire[i].edge], wire[i].sense);
                        if here[0] != vertex {
                            continue;
                        }
                        let previous = wire[(i + wire.len() - 1) % wire.len()];
                        let before = ends(&self.edges[previous.edge], previous.sense);
                        links.push([before[0], here[1]]);
                    }
                }
            }
            check_link(&links, m, w)?;
            m.release(&links);
        }
        Ok(())
    }
    pub fn domains(&self) -> &FiniteDomains<'f, 's> {
        self.domains
    }
    pub fn vertices(&self) -> &[ClosedVertex] {
        &self.vertices
    }
    pub fn edges(&self) -> &[ClosedEdge] {
        &self.edges
    }
    pub fn faces(&self) -> &[ClosedFace] {
        &self.faces
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_bytes(&self) -> usize {
        self.bytes
    }
}
fn check_link(links: &[[usize; 2]], m: &mut Memory, w: &mut Work<'_>) -> Result<(), Error> {
    if links.is_empty() {
        return Err(Error::Unresolved);
    }
    let mut nodes = m.allocate(2 * links.len(), w)?;
    let mut indexed = m.allocate(links.len(), w)?;
    for link in links {
        let mut pair = [0; 2];
        for (k, &v) in link.iter().enumerate() {
            w.charge(nodes.len() + 1)?;
            let index = if let Some(i) = nodes.iter().position(|n: &(usize, usize, bool)| n.0 == v)
            {
                nodes[i].1 += 1;
                i
            } else {
                nodes.push((v, 1, false));
                nodes.len() - 1
            };
            pair[k] = index;
        }
        indexed.push(pair);
    }
    if nodes.iter().any(|n| n.1 != 2) {
        return Err(Error::Unresolved);
    }
    nodes[0].2 = true;
    for _ in 0..nodes.len() {
        let mut changed = false;
        for &[a, b] in &indexed {
            // Both indices were earned by the same completed node pass;
            // nodes are never removed before this adjacency traversal.
            w.charge(1)?;
            if nodes[a].2 != nodes[b].2 {
                nodes[a].2 = true;
                nodes[b].2 = true;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    if nodes.iter().any(|n| !n.2) {
        return Err(Error::Unresolved);
    }
    m.release(&indexed);
    m.release(&nodes);
    Ok(())
}
#[cfg(test)]
mod tests;
