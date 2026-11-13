//! Independent finite incidence over original namespaced source identities.
use super::*;
#[derive(Clone, Copy)]
struct Use {
    face: usize,
    forward: bool,
}
struct Corner {
    germs: [usize; 2],
    next: usize,
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
        let mut queue = m.allocate(faces, w)?;
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
        let mut edges = 0;
        for (id, curve) in self
            .outer
            .source_curves()
            .map(|(i, c)| (i as usize, c))
            .chain(
                self.inner
                    .source_curves()
                    .map(|(i, c)| (self.outer.source_curve_count() + i as usize, c)),
            )
        {
            w.charge(1)?;
            match (curve.faces()[0] != curve.faces()[1], uses[id]) {
                (true, [Some(a), Some(b)]) if a.face != b.face && a.forward != b.forward => {
                    edges += 1
                }
                (false, [None, None]) => {}
                _ => return Err(Error::Unresolved),
            }
        }
        for &head in &heads {
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
        reached[0] = true;
        queue.push(0);
        let mut at = 0;
        while at < queue.len() {
            let face = queue[at];
            at += 1;
            for pair in &uses {
                w.charge(1)?;
                if let [Some(a), Some(b)] = *pair {
                    for (a, b) in [(a, b), (b, a)] {
                        if a.face == face && !reached[b.face] {
                            reached[b.face] = true;
                            queue.push(b.face);
                        }
                    }
                }
            }
        }
        let vertices = heads.iter().filter(|&&h| h != usize::MAX).count();
        if vertices != self.outer.topology_counts().0 + self.inner.topology_counts().0
            || edges != self.outer.topology_counts().1 + self.inner.topology_counts().1
            || reached.iter().any(|v| !*v)
            || vertices as isize - edges as isize + chi != 0
        {
            return Err(Error::Unresolved);
        }
        w.charge(1)?;
        Ok((vertices, edges))
    }
}
