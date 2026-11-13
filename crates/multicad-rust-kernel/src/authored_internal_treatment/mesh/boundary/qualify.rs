//! Independent incidence and vertex-germ links: connected cavity sphere or
//! through-cavity material torus, with face Euler contribution 2 minus loops.
use super::*;
struct Corner {
    germs: [usize; 2],
    next: usize,
}
impl Graph {
    pub(super) fn qualify(&self, w: &mut Work<'_>) -> Result<(usize, usize), Error> {
        let faces = (2
            + 4 * self.source.thread.spec().starts
            + u32::from(self.source.entry.is_some())
            + 2 * u32::from(self.source.tail.is_some())
            + u32::from(self.material.is_some())) as usize;
        let requests = [
            self.vertices.len() * size_of::<usize>(),
            2 * self.curves.len() * size_of::<Corner>(),
            2 * self.curves.len() * size_of::<bool>(),
            faces * size_of::<bool>(),
            faces * size_of::<usize>(),
        ];
        let mut remaining: usize = requests.iter().sum();
        let mut workspace = 0;
        remaining -= requests[0];
        let mut heads = self.buffer(self.vertices.len(), workspace + remaining)?;
        heads.resize(self.vertices.len(), usize::MAX);
        workspace += heads.capacity() * size_of::<usize>();
        remaining -= requests[1];
        let mut corners = self.buffer(2 * self.curves.len(), workspace + remaining)?;
        workspace += corners.capacity() * size_of::<Corner>();
        remaining -= requests[2];
        let mut seen = self.buffer(2 * self.curves.len(), workspace + remaining)?;
        seen.resize(2 * self.curves.len(), false);
        workspace += seen.capacity() * size_of::<bool>();
        remaining -= requests[3];
        let mut reached = self.buffer(faces, workspace + remaining)?;
        reached.resize(faces, false);
        workspace += reached.capacity() * size_of::<bool>();
        let mut queue = self.buffer(faces, workspace)?;
        workspace += queue.capacity() * size_of::<usize>();
        let base =
            2 + 4 * self.source.thread.spec().starts + u32::from(self.source.entry.is_some());
        let mut chi_faces = 0isize;
        for face in 0..faces {
            let loops = self.face_loops_reserved(face as u32, workspace, w)?;
            let expected = if (self.material.is_some() && (face < 2 || face + 1 == faces))
                || (self.source.entry.is_some()
                    && face as u32 == 2 + 4 * self.source.thread.spec().starts)
                || (self.source.tail.is_some() && face as u32 >= base)
            {
                2
            } else {
                1
            };
            if loops.len() != expected {
                return Err(Error::Unresolved);
            }
            let mut held = loops.capacity() * size_of::<Vec<loops::Coedge>>();
            for chain in &loops {
                held += chain.capacity() * size_of::<loops::Coedge>();
            }
            self.check(workspace + held)?;
            chi_faces += 2 - loops.len() as isize;
            for chain in &loops {
                for (i, incoming) in chain.iter().enumerate() {
                    w.step()?;
                    let outgoing = chain[(i + 1) % chain.len()];
                    let first = &self.curves[incoming.curve];
                    let second = &self.curves[outgoing.curve];
                    let a = usize::from(incoming.forward);
                    let b = usize::from(!outgoing.forward);
                    let vertex = first.endpoints[a] as usize;
                    if second.endpoints[b] as usize != vertex {
                        return Err(Error::Unresolved);
                    }
                    let germs = [2 * incoming.curve + a, 2 * outgoing.curve + b];
                    if germs[0] == germs[1] {
                        return Err(Error::Unresolved);
                    }
                    corners.push(Corner {
                        germs,
                        next: heads[vertex],
                    });
                    heads[vertex] = corners.len() - 1;
                }
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
                w.step()?;
                if seen[current] {
                    return Err(Error::Unresolved);
                }
                seen[current] = true;
                let mut neighbors = [usize::MAX; 2];
                let mut count = 0;
                let mut at = head;
                while at != usize::MAX {
                    w.step()?;
                    let link = &corners[at];
                    for i in 0..2 {
                        if link.germs[i] == current {
                            if count == 2 {
                                return Err(Error::Unresolved);
                            }
                            neighbors[count] = link.germs[1 - i];
                            count += 1;
                        }
                    }
                    at = link.next;
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
                w.step()?;
                let link = &corners[at];
                if link.germs.iter().any(|&g| !seen[g]) {
                    return Err(Error::Unresolved);
                }
                at = link.next;
            }
        }
        let mut edges = 0usize;
        for curve in &self.curves {
            w.step()?;
            let [Some(a), Some(b)] = curve.uses else {
                return Err(Error::Unresolved);
            };
            if a.forward == b.forward || a.face as usize >= faces || b.face as usize >= faces {
                return Err(Error::Unresolved);
            }
            if a.face != b.face {
                edges += 1;
            }
        }
        reached[0] = true;
        queue.push(0);
        let mut at = 0;
        while at < queue.len() {
            let face = queue[at];
            at += 1;
            for curve in &self.curves {
                w.step()?;
                let [Some(a), Some(b)] = curve.uses else {
                    return Err(Error::Unresolved);
                };
                for (from, to) in [(a.face, b.face), (b.face, a.face)] {
                    if from as usize == face && !reached[to as usize] {
                        reached[to as usize] = true;
                        queue.push(to as usize);
                    }
                }
            }
        }
        if reached.iter().any(|v| !v) {
            return Err(Error::Unresolved);
        }
        let vertices = heads.iter().filter(|&&v| v != usize::MAX).count();
        if vertices as isize - edges as isize + chi_faces
            != if self.material.is_some() { 0 } else { 2 }
        {
            return Err(Error::Unresolved);
        }
        w.step()?;
        Ok((vertices, edges))
    }
}
