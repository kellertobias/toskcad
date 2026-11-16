//! Actual oriented loops, connected face adjacency and every vertex link.
//! The three annular faces contribute0 to Euler characteristic, not1.
use super::*;
impl Boundary {
    fn vertex_ports(&self, e: Coedge) -> [usize; 2] {
        let i = self.curves.iter().position(|c| c.key == e.curve).unwrap();
        if e.forward {
            [2 * i, 2 * i + 1]
        } else {
            [2 * i + 1, 2 * i]
        }
    }
    fn port_vertex(&self, p: usize) -> VertexKey {
        self.curves[p / 2].vertices[p % 2]
    }
    pub(super) fn validate(
        &self,
        nv: usize,
        ne: usize,
        nf: usize,
        w: &mut Work<'_>,
    ) -> Result<(), Error> {
        w.charge(nv * nv + ne * ne + nf * nf + 2 * ne * ne)?;
        if self.vertices.len() != nv
            || self.curves.len() != ne
            || self.faces.len() != nf
            || self.coedges.len() != 2 * ne
            || self
                .vertices
                .iter()
                .enumerate()
                .any(|(i, v)| self.vertices[..i].contains(v))
            || self
                .curves
                .iter()
                .enumerate()
                .any(|(i, c)| self.curves[..i].iter().any(|p| p.key == c.key))
            || self.faces.iter().enumerate().any(|(i, f)| {
                self.faces[..i].iter().any(|p| p.key == f.key)
                    || self.source.face_index(f.key).is_none()
            })
            || self.faces.iter().any(|f| {
                f.spans[0].0 >= f.spans[0].1
                    || f.spans
                        .iter()
                        .any(|&(a, b)| a > b || b > self.coedges.len())
            })
            || self
                .coedges
                .iter()
                .any(|e| !self.curves.iter().any(|c| c.key == e.curve))
        {
            return Err(Error::Unresolved);
        }
        for c in &self.curves {
            w.charge(self.coedges.len() + self.vertices.len() + self.faces.len())?;
            if c.vertices.iter().any(|v| !self.vertices.contains(v)) {
                return Err(Error::Unresolved);
            }
            let mut uses = [None; 2];
            let mut count = 0;
            for f in &self.faces {
                for (a, b) in f.spans {
                    for e in &self.coedges[a..b] {
                        if e.curve == c.key {
                            if count == 2 {
                                return Err(Error::Unresolved);
                            }
                            uses[count] = Some((f.key, e.forward));
                            count += 1;
                        }
                    }
                }
            }
            if count != 2
                || uses[0].unwrap().1 == uses[1].unwrap().1
                || uses[0].unwrap().0 == uses[1].unwrap().0
                || !uses.iter().flatten().all(|(f, _)| c.faces.contains(f))
            {
                return Err(Error::Unresolved);
            }
        }
        // Each endpoint germ is a distinct port, including both endpoints of a
        // closed plain circle. Compute its two face-corner neighbors only once.
        let mut adjacent = vec![[usize::MAX; 2]; 2 * ne];
        let mut degree = vec![0; 2 * ne];
        for f in &self.faces {
            for (a, b) in f.spans {
                let edges = &self.coedges[a..b];
                for i in 0..edges.len() {
                    w.charge(2 * ne)?;
                    let a = self.vertex_ports(edges[i])[1];
                    let b = self.vertex_ports(edges[(i + 1) % edges.len()])[0];
                    if self.port_vertex(a) != self.port_vertex(b)
                        || a == b
                        || degree[a] == 2
                        || degree[b] == 2
                    {
                        return Err(Error::Unresolved);
                    }
                    adjacent[a][degree[a]] = b;
                    degree[a] += 1;
                    adjacent[b][degree[b]] = a;
                    degree[b] += 1;
                }
            }
        }
        // Exactly two opposite uses of every curve, and every closed face
        // loop, give each endpoint port exactly two corner neighbors.
        debug_assert!(degree.iter().all(|&d| d == 2));
        for &v in &self.vertices {
            w.charge(2 * ne + 36)?;
            let mut ports = [usize::MAX; 6];
            let mut n = 0;
            for (p, _) in adjacent.iter().enumerate() {
                if self.port_vertex(p) == v {
                    if n == ports.len() {
                        return Err(Error::Unresolved);
                    }
                    ports[n] = p;
                    n += 1;
                }
            }
            link(&ports[..n], &adjacent)?;
        }
        let mut reached = [false; 69];
        reached[0] = true;
        for _ in 0..nf {
            for c in &self.curves {
                w.charge(2 * nf)?;
                let a = self.faces.iter().position(|f| f.key == c.faces[0]).unwrap();
                let b = self.faces.iter().position(|f| f.key == c.faces[1]).unwrap();
                if reached[a] || reached[b] {
                    reached[a] = true;
                    reached[b] = true;
                }
            }
        }
        let chi_faces = self
            .faces
            .iter()
            .map(|f| {
                if f.spans[1].0 == f.spans[1].1 {
                    1isize
                } else {
                    0
                }
            })
            .sum::<isize>();
        finish_surface(&reached[..nf], nv, ne, chi_faces)
    }
}
fn link(ports: &[usize], adjacent: &[[usize; 2]]) -> Result<(), Error> {
    if ports.len() < 2 {
        return Err(Error::Unresolved);
    }
    let mut seen = [usize::MAX; 6];
    let mut previous = usize::MAX;
    let mut current = ports[0];
    for i in 0..ports.len() {
        if !ports.contains(&current) || seen[..i].contains(&current) {
            return Err(Error::Unresolved);
        }
        seen[i] = current;
        let [a, b] = adjacent[current];
        let next = if a != previous { a } else { b };
        previous = current;
        current = next;
    }
    if current != ports[0] {
        return Err(Error::Unresolved);
    }
    Ok(())
}

fn finish_surface(reached: &[bool], nv: usize, ne: usize, chi_faces: isize) -> Result<(), Error> {
    if reached.iter().any(|x| !*x) || nv as isize - ne as isize + chi_faces != 2 {
        Err(Error::Unresolved)
    } else {
        Ok(())
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn actual_endpoint_link_and_annular_euler_rules_refuse_branches_and_wrong_disk_count() {
        let cycle = [[1, 3], [0, 2], [1, 3], [2, 0]];
        assert!(link(&[0, 1, 2, 3], &cycle).is_ok());
        assert!(link(&[], &cycle).is_err());
        let disconnected = [[1, 1], [0, 0], [3, 3], [2, 2]];
        assert!(link(&[0, 1, 2, 3], &disconnected).is_err());
        let leaves = [[1, 3], [0, 2], [1, 3], [2, 4], [3, 3]];
        assert!(link(&[0, 1, 2, 3], &leaves).is_err());
        let skips = [[4, 3], [0, 2], [1, 3], [2, 0], [0, 0]];
        assert!(link(&[0, 1, 2, 3], &skips).is_err());
        // Real RH1 graph:18vertices26edges8side disks+2tip disks.
        assert!(finish_surface(&[true; 13], 18, 26, 10).is_ok());
        assert!(finish_surface(&[true; 13], 18, 26, 13).is_err());
        assert!(finish_surface(&[true, false], 18, 26, 10).is_err());
    }
}
