//! Ordered source coedges; no numerical proximity or emitted triangles.
use super::*;
#[derive(Clone, Copy, Debug)]
pub struct Coedge {
    pub curve: usize,
    pub forward: bool,
}
impl Graph {
    #[cfg(test)]
    pub(super) fn face_loops(
        &self,
        face: u32,
        w: &mut Work<'_>,
    ) -> Result<Vec<Vec<Coedge>>, Error> {
        self.face_loops_reserved(face, 0, w)
    }
    pub(super) fn face_loops_reserved(
        &self,
        face: u32,
        reserved: usize,
        w: &mut Work<'_>,
    ) -> Result<Vec<Vec<Coedge>>, Error> {
        let mut edges = self.buffer(self.curves.len(), reserved, w)?;
        for (index, curve) in self.curves.iter().enumerate() {
            w.step()?;
            if curve.uses[1].is_some_and(|b| b.face == curve.uses[0].unwrap().face) {
                continue;
            }
            for side in curve.uses.into_iter().flatten().filter(|s| s.face == face) {
                edges.push(Coedge {
                    curve: index,
                    forward: side.forward,
                });
            }
        }
        self.trace_reserved(&edges, reserved + edges.capacity() * size_of::<Coedge>(), w)
    }
    pub(super) fn rim(&self, end: usize, w: &mut Work<'_>) -> Result<Vec<Coedge>, Error> {
        let mut edges = self.buffer(self.curves.len(), 0, w)?;
        for (index, curve) in self.curves.iter().enumerate() {
            w.step()?;
            if curve.uses[1].is_none() && curve.station == Some(if end == 0 { 0 } else { 2 }) {
                edges.push(Coedge {
                    curve: index,
                    forward: curve.uses[0].unwrap().forward,
                });
            }
        }
        let mut loops = self.trace_reserved(&edges, edges.capacity() * size_of::<Coedge>(), w)?;
        if loops.len() != 1 {
            return Err(Error::Unresolved);
        }
        let rim = loops.pop().unwrap();
        // Canonical shared endpoint angles telescope exactly. Only explicit
        // true-periodic integer jumps survive, so this proves winding without
        // testing whether a rounded angular sum happens to enclose2pi.
        let mut winding = 0i64;
        let positive = (end == 0) != self.source.thread.spec().left_hand;
        for coedge in &rim {
            let curve = &self.curves[coedge.curve];
            let direction = if coedge.forward {
                curve.v_direction
            } else {
                curve.v_direction.reverse()
            };
            if direction
                != if positive {
                    Ordering::Greater
                } else {
                    Ordering::Less
                }
            {
                return Err(Error::Unresolved);
            }
            winding += i64::from(curve.wrap_jump) * if coedge.forward { 1 } else { -1 };
        }
        if winding != if positive { 1 } else { -1 } {
            return Err(Error::Unresolved);
        }
        Ok(rim)
    }
    #[cfg(test)]
    pub(super) fn trace(
        &self,
        edges: &[Coedge],
        w: &mut Work<'_>,
    ) -> Result<Vec<Vec<Coedge>>, Error> {
        self.trace_reserved(edges, std::mem::size_of_val(edges), w)
    }
    fn trace_reserved(
        &self,
        edges: &[Coedge],
        reserved: usize,
        w: &mut Work<'_>,
    ) -> Result<Vec<Vec<Coedge>>, Error> {
        w.step()?;
        if edges.is_empty() {
            return Err(Error::Unresolved);
        }
        let count = self.vertices.len();
        let seen_request = edges.len() * size_of::<bool>();
        let mut incoming = self.buffer::<usize>(
            count,
            reserved + count * size_of::<usize>() + seen_request,
            w,
        )?;
        incoming.resize(count, usize::MAX);
        let incoming_bytes = incoming.capacity() * size_of::<usize>();
        let mut outgoing = self.buffer(count, reserved + incoming_bytes + seen_request, w)?;
        outgoing.resize(count, usize::MAX);
        let mut seen = self.buffer(
            edges.len(),
            reserved + incoming_bytes + outgoing.capacity() * size_of::<usize>(),
            w,
        )?;
        seen.resize(edges.len(), false);
        let bytes = reserved
            + incoming_bytes
            + outgoing.capacity() * size_of::<usize>()
            + seen.capacity() * size_of::<bool>();
        for (i, edge) in edges.iter().enumerate() {
            w.step()?;
            let curve = self.curves.get(edge.curve).ok_or(Error::Unresolved)?;
            let [a, b] = if edge.forward {
                curve.endpoints
            } else {
                [curve.endpoints[1], curve.endpoints[0]]
            };
            let slot = outgoing.get_mut(a as usize).ok_or(Error::Unresolved)?;
            if *slot != usize::MAX {
                return Err(Error::Unresolved);
            }
            *slot = i;
            let slot = incoming.get_mut(b as usize).ok_or(Error::Unresolved)?;
            if *slot != usize::MAX {
                return Err(Error::Unresolved);
            }
            *slot = i;
        }
        if incoming
            .iter()
            .zip(&outgoing)
            .any(|(a, b)| (*a == usize::MAX) != (*b == usize::MAX))
        {
            return Err(Error::Unresolved);
        }
        let mut loops = self.buffer::<Vec<Coedge>>(2, bytes, w)?;
        let mut retained = loops.capacity() * size_of::<Vec<Coedge>>();
        for start in 0..edges.len() {
            if seen[start] {
                continue;
            }
            if loops.len() == 2 {
                return Err(Error::Unresolved);
            }
            let mut chain = self.buffer(edges.len(), bytes + retained, w)?;
            retained += chain.capacity() * size_of::<Coedge>();
            self.check(bytes + retained)?;
            let mut current = start;
            loop {
                w.step()?;
                // One incoming and outgoing edge per used vertex was proved
                // above. A finite bijection is disjoint cycles; this walk can
                // only revisit its starting edge (handled below).
                seen[current] = true;
                let coedge = edges[current];
                chain.push(coedge);
                let curve = &self.curves[coedge.curve];
                let end = curve.endpoints[usize::from(coedge.forward)];
                current = outgoing[end as usize];
                if current == start {
                    break;
                }
            }
            loops.push(chain);
        }
        Ok(loops)
    }
}
