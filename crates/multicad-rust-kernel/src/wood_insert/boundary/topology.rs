//! Independent finite incidence over original namespaced source identities.
use super::*;
#[derive(Clone, Copy)]
pub(super) struct Use {
    pub face: usize,
    pub forward: bool,
}
pub(super) struct Corner {
    pub germs: [usize; 2],
    pub next: usize,
}
impl Boundary {
    fn vertex_index(&self, key: BoundaryVertexKey) -> usize {
        match key {
            BoundaryVertexKey::Outer(i) => i as usize,
            BoundaryVertexKey::Inner(i) => self.outer.points().len() + i as usize,
        }
    }
    fn curve_index(&self, key: BoundaryCurveKey) -> usize {
        match key {
            BoundaryCurveKey::Outer(i) => i as usize,
            BoundaryCurveKey::Inner(i) => self.outer.source_curve_count() + i as usize,
        }
    }
    pub(super) fn qualify(
        &self,
        m: &mut Memory,
        w: &mut Work<'_>,
    ) -> Result<(usize, usize), Error> {
        if self.loops.len() != self.face_count() as usize {
            return Err(Error::Unresolved);
        }
        let points = self.outer.points().len() + self.inner.points().len();
        let curves = self.outer.source_curve_count() + self.inner.source_curve_count();
        let faces = self.face_count() as usize;
        let mut uses = m.allocate(curves, w)?;
        uses.resize(curves, [None::<Use>; 2]);
        let mut heads = m.allocate(points, w)?;
        heads.resize(points, usize::MAX);
        let mut corners = m.allocate(2 * curves, w)?;
        let mut seen = m.allocate(2 * curves, w)?;
        seen.resize(2 * curves, false);
        let mut reached = m.allocate(faces, w)?;
        reached.resize(faces, false);
        let queue = m.allocate(faces, w)?;
        let mut chi = 0isize;
        for (face, loops) in self.loops.iter().enumerate() {
            w.charge(1)?;
            if loops.len() != if face < 2 { 2 } else { 1 } {
                return Err(Error::Unresolved);
            }
            chi += 2 - loops.len() as isize;
            for chain in loops {
                if chain.is_empty() {
                    return Err(Error::Unresolved);
                }
                for (i, edge) in chain.iter().enumerate() {
                    w.charge(1)?;
                    let source = self.curve(edge.curve).ok_or(Error::Unresolved)?;
                    let expected = source
                        .face_uses()
                        .into_iter()
                        .find(|&(f, _)| self.global_face(edge.curve.owner(), f) == face as u32)
                        .ok_or(Error::Unresolved)?;
                    if edge.forward != (expected.1 != (edge.curve.owner() == Side::Inner)) {
                        return Err(Error::Unresolved);
                    }
                    if source.faces()[0] == source.faces()[1] {
                        return Err(Error::Unresolved);
                    }
                    let id = self.curve_index(edge.curve);
                    let pair = uses.get_mut(id).ok_or(Error::Unresolved)?;
                    let at = usize::from(pair[0].is_some());
                    if pair[at].is_some() {
                        return Err(Error::Unresolved);
                    }
                    pair[at] = Some(Use {
                        face,
                        forward: edge.forward,
                    });
                    let next = chain[(i + 1) % chain.len()];
                    let a = self.curve_endpoints(edge.curve).ok_or(Error::Unresolved)?;
                    let b = self.curve_endpoints(next.curve).ok_or(Error::Unresolved)?;
                    if a[usize::from(edge.forward)] != b[usize::from(!next.forward)] {
                        return Err(Error::Unresolved);
                    }
                    let vertex = self.vertex_index(a[usize::from(edge.forward)]);
                    let head = heads.get_mut(vertex).ok_or(Error::Unresolved)?;
                    corners.push(Corner {
                        germs: [
                            2 * id + usize::from(edge.forward),
                            2 * self.curve_index(next.curve) + usize::from(!next.forward),
                        ],
                        next: *head,
                    });
                    *head = corners.len() - 1;
                }
            }
        }
        Incidence {
            uses,
            heads,
            corners,
            seen,
            reached,
            queue,
            chi,
        }
        .finish(
            self.outer
                .source_curves()
                .chain(self.inner.source_curves())
                .map(|(_, c)| c.faces()[0] != c.faces()[1]),
            (
                self.outer.topology_counts().0 + self.inner.topology_counts().0,
                self.outer.topology_counts().1 + self.inner.topology_counts().1,
            ),
            w,
        )
    }
}

/// Collected finite incidence is validated without trusting a geometric source.
/// Keeping this stage separate permits malformed graph controls without forging
/// an admitted thread or insert.
pub(super) struct Incidence {
    pub uses: Vec<[Option<Use>; 2]>,
    pub heads: Vec<usize>,
    pub corners: Vec<Corner>,
    pub seen: Vec<bool>,
    pub reached: Vec<bool>,
    pub queue: Vec<usize>,
    pub chi: isize,
}
impl Incidence {
    pub fn finish(
        mut self,
        physical: impl IntoIterator<Item = bool>,
        expected: (usize, usize),
        w: &mut Work<'_>,
    ) -> Result<(usize, usize), Error> {
        let edges = edge_count(physical.into_iter().zip(self.uses.iter().copied()), w)?;
        vertex_links(&self.heads, &self.corners, &mut self.seen, w)?;
        self.reached[0] = true;
        self.queue.push(0);
        let mut at = 0;
        while at < self.queue.len() {
            let face = self.queue[at];
            at += 1;
            for pair in &self.uses {
                w.charge(1)?;
                if let [Some(a), Some(b)] = *pair {
                    for (a, b) in [(a, b), (b, a)] {
                        if a.face == face && !self.reached[b.face] {
                            self.reached[b.face] = true;
                            self.queue.push(b.face);
                        }
                    }
                }
            }
        }
        let vertices = self.heads.iter().filter(|&&h| h != usize::MAX).count();
        w.charge(1)?;
        shell_receipt((vertices, edges), expected, self.chi, &self.reached)
    }
}

/// Validate and count actual physical edges independently of wire collection.
pub(super) fn edge_count(
    entries: impl IntoIterator<Item = (bool, [Option<Use>; 2])>,
    w: &mut Work<'_>,
) -> Result<usize, Error> {
    let mut count = 0;
    for (physical, uses) in entries {
        w.charge(1)?;
        count += usize::from(edge_receipt(physical, uses)?);
    }
    Ok(count)
}
pub(super) fn edge_receipt(physical: bool, uses: [Option<Use>; 2]) -> Result<bool, Error> {
    match (physical, uses) {
        (true, [Some(a), Some(b)]) if a.face != b.face && a.forward != b.forward => Ok(true),
        (false, [None, None]) => Ok(false),
        _ => Err(Error::Unresolved),
    }
}
pub(super) fn vertex_links(
    heads: &[usize],
    corners: &[Corner],
    seen: &mut [bool],
    w: &mut Work<'_>,
) -> Result<(), Error> {
    for &head in heads {
        if head == usize::MAX {
            continue;
        }
        let start = corners[head].germs[0];
        let mut current = start;
        let mut previous = usize::MAX;
        loop {
            w.charge(1)?;
            let flag = seen.get_mut(current).ok_or(Error::Unresolved)?;
            if *flag {
                return Err(Error::Unresolved);
            }
            *flag = true;
            let mut neighbors = [usize::MAX; 2];
            let mut count = 0;
            let mut at = head;
            while at != usize::MAX {
                w.charge(1)?;
                let c = &corners[at];
                for i in 0..2 {
                    if c.germs[i] == current {
                        if count == 2 {
                            return Err(Error::Unresolved);
                        }
                        neighbors[count] = c.germs[1 - i];
                        count += 1;
                    }
                }
                at = c.next;
            }
            if count != 2 {
                return Err(Error::Unresolved);
            }
            let next = if neighbors[0] == previous {
                neighbors[1]
            } else {
                neighbors[0]
            };
            if next == start {
                break;
            }
            previous = current;
            current = next;
        }
        let mut at = head;
        while at != usize::MAX {
            w.charge(1)?;
            if corners[at].germs.iter().any(|&g| !seen[g]) {
                return Err(Error::Unresolved);
            }
            at = corners[at].next;
        }
    }
    Ok(())
}
pub(super) fn shell_receipt(
    actual: (usize, usize),
    expected: (usize, usize),
    chi: isize,
    reached: &[bool],
) -> Result<(usize, usize), Error> {
    if actual != expected
        || reached.iter().any(|v| !*v)
        || actual.0 as isize - actual.1 as isize + chi != 0
    {
        Err(Error::Unresolved)
    } else {
        Ok(actual)
    }
}
