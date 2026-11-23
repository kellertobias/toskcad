//! Remove straight-edge subdivisions that only record how one source face was
//! clipped into convex fragments.
//!
//! Clipping a face against a convex cutter leaves convex fragments whose seams
//! can end in the middle of an original straight edge. That vertex has no
//! geometric meaning: both sides of the edge belong to the same two source
//! faces. Leaving it splits one selectable edge into collinear pieces and
//! makes a whole-edge reference ambiguous.
//!
//! A vertex is removed only when it is not an input vertex, its seams separate
//! coplanar fragments of one source face, and once those are merged exactly
//! two collinear edges remain at it. A vertex where any other edge starts
//! (a genuine Boolean split) and every authored input vertex stay explicit.
use super::*;
use std::collections::BTreeSet;

impl Boundary {
    /// Merge the fragments that meet at an irrelevant straight-edge vertex and
    /// drop the vertex. `protected` lists every input vertex. Must run after
    /// [`Self::split_edges`], when every incidence is installed in both loops.
    pub(in crate::planar_boolean) fn remove_collinear_splits(
        &mut self,
        protected: impl Iterator<Item = Point>,
    ) -> Result<(), Error> {
        // Bounded separately from the Boolean that produced this boundary:
        // every step is linear in the incidences, which are already bounded.
        let mut work = MAX_WORK;
        let budget = &mut work;
        let protected = self.protected(protected, budget)?;
        let edges = self.edges(budget)?;
        let mut merge = UnionFind::new(self.faces.len());
        for (vertex, _) in protected.iter().enumerate().filter(|(_, kept)| !**kept) {
            if self.straight(&edges, vertex).is_some() {
                for (_, uses) in edges.at(vertex) {
                    if let [(a, _), (b, _)] = uses[..]
                        && self.seam(uses)
                    {
                        merge.join(a, b);
                    }
                }
            }
        }
        self.merge_fragments(&mut merge, budget)?;
        loop {
            let edges = self.edges(budget)?;
            // Remove an independent set per pass, so each removal sees the
            // final neighbors of its straight run and never doubles an edge.
            let mut blocked = vec![false; self.vertices.len()];
            let mut chosen = Vec::new();
            for vertex in 0..self.vertices.len() {
                if protected[vertex] || blocked[vertex] {
                    continue;
                }
                let Some((ends, 0)) = self.straight(&edges, vertex) else {
                    continue;
                };
                if edges
                    .uses
                    .contains_key(&[ends[0].min(ends[1]), ends[0].max(ends[1])])
                {
                    continue;
                }
                blocked[vertex] = true;
                for end in ends {
                    blocked[end] = true;
                }
                chosen.push(vertex);
            }
            if chosen.is_empty() {
                break;
            }
            for face in &mut self.faces {
                charge(budget, face.vertices.len())?;
                face.vertices
                    .retain(|vertex| chosen.binary_search(vertex).is_err());
            }
        }
        self.compact(budget)
    }

    fn protected(
        &self,
        points: impl Iterator<Item = Point>,
        budget: &mut usize,
    ) -> Result<Vec<bool>, Error> {
        let mut points: Vec<Point> = points.collect();
        charge(budget, points.len().saturating_mul(2))?;
        points.sort_unstable_by(|a, b| a[0].total_cmp(&b[0]));
        let mut flags = Vec::with_capacity(self.vertices.len());
        for &vertex in &self.vertices {
            let first = points.partition_point(|p| p[0] < vertex[0] - self.tolerance);
            let near = points[first..]
                .iter()
                .take_while(|p| p[0] <= vertex[0] + self.tolerance);
            let mut found = false;
            for &point in near {
                charge(budget, 1)?;
                if norm(sub(point, vertex)) <= self.tolerance {
                    found = true;
                    break;
                }
            }
            flags.push(found);
        }
        Ok(flags)
    }

    fn edges(&self, budget: &mut usize) -> Result<Edges, Error> {
        let mut uses: BTreeMap<[usize; 2], Vec<(usize, bool)>> = BTreeMap::new();
        for (index, face) in self.faces.iter().enumerate() {
            charge(budget, face.vertices.len())?;
            for (&a, &b) in face
                .vertices
                .iter()
                .zip(face.vertices.iter().cycle().skip(1))
            {
                uses.entry([a.min(b), a.max(b)])
                    .or_default()
                    .push((index, a < b));
            }
        }
        let mut star = vec![Vec::new(); self.vertices.len()];
        for &[a, b] in uses.keys() {
            star[a].push([a, b]);
            star[b].push([a, b]);
        }
        Ok(Edges { uses, star })
    }

    /// Two opposite uses by coplanar, same-facing fragments of one source face.
    fn seam(&self, uses: &[(usize, bool)]) -> bool {
        let [(a, forward_a), (b, forward_b)] = uses[..] else {
            return false;
        };
        let (a, b) = (&self.faces[a], &self.faces[b]);
        forward_a != forward_b
            && a.source == b.source
            && a.source_face == b.source_face
            && dot(a.normal, b.normal) > 1. - 1e-12
    }

    /// The two collinear non-seam edges through `vertex` and its seam count,
    /// when `vertex` lies strictly between their far ends.
    fn straight(&self, edges: &Edges, vertex: usize) -> Option<([usize; 2], usize)> {
        let mut ends = Vec::with_capacity(2);
        let mut seams = 0;
        for (key, uses) in edges.at(vertex) {
            if self.seam(uses) {
                seams += 1;
            } else {
                if ends.len() == 2 {
                    return None;
                }
                ends.push(key[0] + key[1] - vertex);
            }
        }
        let [u, w] = ends[..] else {
            return None;
        };
        let start = self.vertices[u];
        let direction = sub(self.vertices[w], start);
        let length_squared = dot(direction, direction);
        let margin = self.tolerance / length_squared.sqrt();
        let relative = sub(self.vertices[vertex], start);
        let fraction = dot(relative, direction) / length_squared;
        (fraction > margin
            && fraction < 1. - margin
            && norm(sub(relative, scale(direction, fraction))) <= self.tolerance)
            .then_some(([u, w], seams))
    }

    /// Replace each group of fragments by its outer loop when that is one
    /// simple loop whose interior vertices no other face uses.
    fn merge_fragments(&mut self, merge: &mut UnionFind, budget: &mut usize) -> Result<(), Error> {
        let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for face in 0..self.faces.len() {
            groups.entry(merge.root(face)).or_default().push(face);
        }
        let mut removed = vec![false; self.faces.len()];
        for members in groups.into_values().filter(|members| members.len() > 1) {
            let Some(ring) = self.outer_ring(&members, budget)? else {
                continue;
            };
            let first = members[0];
            self.faces[first].vertices = ring;
            for &member in &members[1..] {
                removed[member] = true;
            }
        }
        let mut index = 0;
        self.faces.retain(|_| {
            index += 1;
            !removed[index - 1]
        });
        Ok(())
    }

    fn outer_ring(
        &self,
        members: &[usize],
        budget: &mut usize,
    ) -> Result<Option<Vec<usize>>, Error> {
        let mut directed = BTreeMap::new();
        for &member in members {
            let ring = &self.faces[member].vertices;
            charge(budget, ring.len())?;
            for (&a, &b) in ring.iter().zip(ring.iter().cycle().skip(1)) {
                if directed.remove(&(b, a)).is_none() && directed.insert((a, b), ()).is_some() {
                    return Ok(None);
                }
            }
        }
        let mut next = BTreeMap::new();
        for &(a, b) in directed.keys() {
            if next.insert(a, b).is_some() {
                return Ok(None);
            }
        }
        let Some((&start, _)) = next.iter().next() else {
            return Ok(None);
        };
        let mut ring = vec![start];
        let mut current = next[&start];
        while current != start {
            charge(budget, 1)?;
            if ring.len() > next.len() {
                return Ok(None);
            }
            ring.push(current);
            current = *next.get(&current).ok_or(Error::Unresolved)?;
        }
        if ring.len() != next.len() || ring.len() < 3 {
            return Ok(None);
        }
        // A vertex dropped from the group's loops must not support any other
        // face; otherwise merging would detach it from the shell.
        let mut outside = BTreeSet::new();
        for (index, face) in self.faces.iter().enumerate() {
            if !members.contains(&index) {
                charge(budget, face.vertices.len())?;
                outside.extend(face.vertices.iter().copied());
            }
        }
        let kept: BTreeSet<_> = ring.iter().copied().collect();
        if members
            .iter()
            .flat_map(|&member| &self.faces[member].vertices)
            .any(|vertex| !kept.contains(vertex) && outside.contains(vertex))
        {
            return Ok(None);
        }
        Ok(Some(ring))
    }

    fn compact(&mut self, budget: &mut usize) -> Result<(), Error> {
        let mut map = vec![usize::MAX; self.vertices.len()];
        for face in &self.faces {
            charge(budget, face.vertices.len())?;
            for &vertex in &face.vertices {
                map[vertex] = 0;
            }
        }
        let mut vertices = Vec::with_capacity(self.vertices.len());
        for (old, slot) in map.iter_mut().enumerate() {
            if *slot == 0 {
                *slot = vertices.len();
                vertices.push(self.vertices[old]);
            }
        }
        for face in &mut self.faces {
            for vertex in &mut face.vertices {
                *vertex = map[*vertex];
            }
        }
        self.vertices = vertices;
        Ok(())
    }
}

struct Edges {
    uses: BTreeMap<[usize; 2], Vec<(usize, bool)>>,
    star: Vec<Vec<[usize; 2]>>,
}

impl Edges {
    fn at(&self, vertex: usize) -> impl Iterator<Item = ([usize; 2], &[(usize, bool)])> {
        self.star[vertex]
            .iter()
            .map(|key| (*key, self.uses[key].as_slice()))
    }
}

struct UnionFind(Vec<usize>);

impl UnionFind {
    fn new(len: usize) -> Self {
        Self((0..len).collect())
    }

    fn root(&mut self, mut item: usize) -> usize {
        while self.0[item] != item {
            self.0[item] = self.0[self.0[item]];
            item = self.0[item];
        }
        item
    }

    fn join(&mut self, a: usize, b: usize) {
        let (a, b) = (self.root(a), self.root(b));
        self.0[a.max(b)] = a.min(b);
    }
}
