//! Complete paired-source candidate; Native literal embedding is mandatory.
use super::{chord::CapProof, *};
use crate::{
    authored_flat_thread::{Error as SourceError, ShellEnd, ThreadShell},
    display::{DisplayMesh, Tessellation},
};
use std::mem::size_of;
#[cfg(test)]
mod tests;
const WORK: usize = 16_000_000;
const BYTES: usize = 128 * 1024 * 1024;
const SCRATCH: usize = 4 * 1024 * 1024;
type R<T> = Result<T, Error>;
pub struct Candidate {
    source: Insert,
    mesh: DisplayMesh,
    chord: f64,
    normal: f64,
    work: usize,
    source_volume: [f64; 2],
    literal_volume: [f64; 2],
}
impl Candidate {
    pub fn source(&self) -> &Insert {
        &self.source
    }
    pub fn mesh(&self) -> &DisplayMesh {
        &self.mesh
    }
    pub fn into_mesh(self) -> DisplayMesh {
        self.mesh
    }
    pub fn source_chord_bound_mm(&self) -> f64 {
        self.chord
    }
    pub fn source_normal_angle_bound_rad(&self) -> f64 {
        self.normal
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn source_volume_mm3(&self) -> [f64; 2] {
        self.source_volume
    }
    /// Signed volume of these actual literal facets, independently enclosed.
    /// This is not exact source mass and not independently decoded-file mass.
    pub fn literal_volume_mm3(&self) -> [f64; 2] {
        self.literal_volume
    }
}
struct Memory {
    used: usize,
}
impl Memory {
    fn check(&self, extra: usize) -> R<()> {
        if extra > BYTES - self.used {
            Err(Error::Limit)
        } else {
            Ok(())
        }
    }
    fn add(&mut self, n: usize) -> R<()> {
        self.check(n)?;
        self.used += n;
        Ok(())
    }
    fn remove(&mut self, n: usize) {
        self.used -= n;
    }
}
fn allocate<T>(count: usize, m: &mut Memory, w: &mut Work<'_>) -> R<Vec<T>> {
    w.charge(1)?;
    let bytes = count.checked_mul(size_of::<T>()).ok_or(Error::Limit)?;
    m.check(bytes)?;
    let mut v = Vec::new();
    v.try_reserve_exact(count).map_err(|_| Error::Limit)?;
    m.add(v.capacity() * size_of::<T>())?;
    Ok(v)
}
fn failure(e: SourceError) -> Error {
    match e {
        SourceError::Invalid => Error::Invalid,
        SourceError::Unsupported => Error::Unsupported,
        SourceError::Unresolved => Error::Unresolved,
        SourceError::Limit => Error::Limit,
        SourceError::Cancelled => Error::Cancelled,
    }
}
struct Node {
    neighbors: [usize; 2],
    degree: usize,
    angle: Angle,
}
struct Vertex {
    index: usize,
    angle: Angle,
}
fn encoded(q: &Q, m: &mut Memory, w: &mut Work<'_>) -> R<Angle> {
    let maximum =
        q.numer().bits().div_ceil(8).max(1) as usize + q.denom().bits().div_ceil(8).max(1) as usize;
    m.check(maximum)?;
    let a = Angle::encode(q, w)?;
    m.add(a.bytes())?;
    Ok(a)
}
fn ring(
    shell: &ThreadShell,
    end: ShellEnd,
    offset: usize,
    m: &mut Memory,
    w: &mut Work<'_>,
) -> R<Vec<Vertex>> {
    let count = shell.mesh().vertices.len() / 3;
    let mut nodes: Vec<Option<Node>> = allocate(count, m, w)?;
    w.charge(count)?;
    nodes.resize_with(count, || None);
    let headers = nodes.capacity() * size_of::<Option<Node>>();
    let mut unique = 0;
    for edge in 0..shell.rim_count(end) {
        w.charge(1)?;
        let e = shell.rim_edge(end, edge).ok_or(Error::Unresolved)?;
        let [a, b] = e.vertex_indices;
        if a == b || a >= count || b >= count {
            return Err(Error::Unresolved);
        }
        for (vertex, neighbor) in [(a, b), (b, a)] {
            if nodes[vertex].is_none() {
                let (q, cost) = shell
                    .exact_rim_turns(vertex, w.limit - w.used, w.interrupt)
                    .map_err(failure)?;
                w.charge(cost)?;
                nodes[vertex] = Some(Node {
                    neighbors: [usize::MAX; 2],
                    degree: 0,
                    angle: encoded(&q, m, w)?,
                });
                unique += 1;
            }
            let n = nodes[vertex].as_mut().unwrap();
            if n.degree >= 2 || n.neighbors.contains(&neighbor) {
                return Err(Error::Unresolved);
            }
            n.neighbors[n.degree] = neighbor;
            n.degree += 1;
        }
    }
    if unique < 3 {
        return Err(Error::Unresolved);
    }
    let mut start: Option<usize> = None;
    for (i, node) in nodes.iter().enumerate() {
        w.charge(1)?;
        if let Some(node) = node {
            if node.degree != 2 {
                return Err(Error::Unresolved);
            }
            if let Some(old) = start {
                let old: &Node = nodes[old].as_ref().unwrap();
                if cmp(&node.angle.decode(w)?, &old.angle.decode(w)?, w)? == Ordering::Less {
                    start = Some(i);
                }
            } else {
                start = Some(i);
            }
        }
    }
    let start = start.unwrap();
    let mut result = allocate(unique, m, w)?;
    let first = nodes[start].take().unwrap();
    let [a, b] = first.neighbors;
    let next = if cmp(
        &nodes[a]
            .as_ref()
            .ok_or(Error::Unresolved)?
            .angle
            .decode(w)?,
        &nodes[b]
            .as_ref()
            .ok_or(Error::Unresolved)?
            .angle
            .decode(w)?,
        w,
    )? == Ordering::Less
    {
        a
    } else {
        b
    };
    result.push(Vertex {
        index: start + offset,
        angle: first.angle,
    });
    let mut previous = start;
    let mut current = next;
    while current != start {
        w.charge(1)?;
        let node = nodes
            .get_mut(current)
            .and_then(Option::take)
            .ok_or(Error::Unresolved)?;
        let next = if node.neighbors[0] == previous {
            node.neighbors[1]
        } else if node.neighbors[1] == previous {
            node.neighbors[0]
        } else {
            return Err(Error::Unresolved);
        };
        if cmp(
            &result.last().unwrap().angle.decode(w)?,
            &node.angle.decode(w)?,
            w,
        )? != Ordering::Less
        {
            return Err(Error::Unresolved);
        }
        result.push(Vertex {
            index: current + offset,
            angle: node.angle,
        });
        previous = current;
        current = next;
    }
    if result.len() != unique {
        return Err(Error::Unresolved);
    }
    m.remove(headers);
    Ok(result)
}
fn point(ring: &[Vertex], index: usize, t: i32, w: &mut Work<'_>) -> R<(usize, [Q; 2])> {
    let v = &ring[index % ring.len()];
    let mut q = v.angle.decode(w)?;
    if index == ring.len() {
        q = op(&q, &Q::one(), Op::Add, w)?;
    }
    Ok((v.index, [q, Q::from_integer(t.into())]))
}
fn cap(
    mesh: &mut DisplayMesh,
    source: &CapDomains,
    outer: &[Vertex],
    inner: &[Vertex],
    options: Tessellation,
    w: &mut Work<'_>,
) -> R<f64> {
    let proof = CapProof::new(source, w)?;
    let mut o = 0;
    let mut i = 0;
    let mut maximum = 0_f64;
    let z = if source.end() == CapEnd::Start {
        0.
    } else {
        source.owner().full_length_mm()
    };
    let mut area = Q::zero();
    while o < outer.len() || i < inner.len() {
        w.charge(1)?;
        let advance_outer = if i == inner.len() {
            true
        } else if o == outer.len() {
            false
        } else {
            cmp(
                &point(outer, o + 1, 1, w)?.1[0],
                &point(inner, i + 1, 0, w)?.1[0],
                w,
            )? != Ordering::Greater
        };
        let vertices = if advance_outer {
            let v = [
                point(outer, o, 1, w)?,
                point(inner, i, 0, w)?,
                point(outer, o + 1, 1, w)?,
            ];
            o += 1;
            v
        } else {
            let v = [
                point(inner, i, 0, w)?,
                point(inner, i + 1, 0, w)?,
                point(outer, o, 1, w)?,
            ];
            i += 1;
            v
        };
        let mut ids = vertices.each_ref().map(|v| v.0);
        let mut uv = vertices.map(|v| v.1);
        let du1 = op(&uv[1][0], &uv[0][0], Op::Sub, w)?;
        let du2 = op(&uv[2][0], &uv[0][0], Op::Sub, w)?;
        let dt1 = op(&uv[1][1], &uv[0][1], Op::Sub, w)?;
        let dt2 = op(&uv[2][1], &uv[0][1], Op::Sub, w)?;
        let a = op(&du1, &dt2, Op::Mul, w)?;
        let b = op(&du2, &dt1, Op::Mul, w)?;
        let determinant = op(&a, &b, Op::Sub, w)?;
        if determinant <= Q::zero() {
            return Err(Error::Unresolved);
        }
        area = op(&area, &determinant, Op::Add, w)?;
        let mut points = ids.map(|id| mesh.vertices[id * 3..id * 3 + 3].try_into().unwrap());
        if points.iter().any(|p: &[f64; 3]| p[2] != z) {
            return Err(Error::Unresolved);
        }
        // Exact literal planar determinant: no tolerance or rounded normal is used.
        let x1 = op(
            &super::chord::binary_work(points[1][0], w)?,
            &super::chord::binary_work(points[0][0], w)?,
            Op::Sub,
            w,
        )?;
        let x2 = op(
            &super::chord::binary_work(points[2][0], w)?,
            &super::chord::binary_work(points[0][0], w)?,
            Op::Sub,
            w,
        )?;
        let y1 = op(
            &super::chord::binary_work(points[1][1], w)?,
            &super::chord::binary_work(points[0][1], w)?,
            Op::Sub,
            w,
        )?;
        let y2 = op(
            &super::chord::binary_work(points[2][1], w)?,
            &super::chord::binary_work(points[0][1], w)?,
            Op::Sub,
            w,
        )?;
        let a = op(&x1, &y2, Op::Mul, w)?;
        let b = op(&x2, &y1, Op::Mul, w)?;
        let det = op(&a, &b, Op::Sub, w)?;
        if det.is_zero() {
            return Err(Error::Unresolved);
        }
        let expected_positive = source.end() == CapEnd::End;
        if (det > Q::zero()) != expected_positive {
            ids.swap(1, 2);
            uv.swap(1, 2);
            points.swap(1, 2);
        }
        let bound = proof.certify(&uv, points, options.deflection, w)?;
        maximum = maximum.max(bound);
        if mesh.face_indices.len() >= options.triangle_limit {
            return Err(Error::Limit);
        }
        mesh.indices.extend(ids.map(|x| x as u32));
        mesh.face_indices
            .push(if source.end() == CapEnd::Start { 0 } else { 1 });
    }
    // Sum of signed double areas is exactly2 for the complete unit annulus
    // parallelogram, regardless of the two independently chosen seam origins.
    if area != Q::from_integer(2.into()) {
        return Err(Error::Unresolved);
    }
    Ok(maximum)
}
fn literal_volume(mesh: &DisplayMesh, w: &mut Work<'_>) -> R<[f64; 2]> {
    let origin: [f64; 3] = mesh.vertices[..3].try_into().unwrap();
    let mut sum = I::point(0.);
    for tri in mesh.indices.as_chunks::<3>().0 {
        w.charge(64)?;
        let p = tri.map(|id| {
            std::array::from_fn::<_, 3, _>(|i| {
                I::point(mesh.vertices[id as usize * 3 + i]) - I::point(origin[i])
            })
        });
        let c = [
            p[1][1] * p[2][2] - p[1][2] * p[2][1],
            p[1][2] * p[2][0] - p[1][0] * p[2][2],
            p[1][0] * p[2][1] - p[1][1] * p[2][0],
        ];
        sum = sum + (p[0][0] * c[0] + p[0][1] * c[1] + p[0][2] * c[2]) / I::point(6.);
    }
    if !sum.valid() || sum.lo <= 0. {
        return Err(Error::Unresolved);
    }
    Ok([sum.lo, sum.hi])
}
impl Insert {
    pub fn tessellate_candidate(
        &self,
        options: Tessellation,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> R<Candidate> {
        options.check().map_err(|_| Error::Invalid)?;
        let mut w = Work {
            used: 0,
            limit: WORK,
            bytes: BYTES,
            interrupt,
        };
        w.charge(1)?;
        let mut memory = Memory {
            used: SCRATCH + size_of::<Candidate>(),
        };
        let cap_options = |w: &Work<'_>| CapOptions {
            maximum_work: w.limit - w.used,
            maximum_bytes: 1024 * 1024,
        };
        let bottom = self.cap_domains(CapEnd::Start, cap_options(&w), interrupt)?;
        w.charge(bottom.construction_work())?;
        memory.add(bottom.retained_bytes())?;
        let top = self.cap_domains(CapEnd::End, cap_options(&w), interrupt)?;
        w.charge(top.construction_work())?;
        memory.add(top.retained_bytes())?;
        let outer = self
            .outer()
            .tessellate_shell(options, w.limit - w.used, interrupt)
            .map_err(failure)?;
        w.charge(outer.work_used())?;
        memory.add(outer.retained_bytes())?;
        let inner = self
            .inner()
            .tessellate_shell(options, w.limit - w.used, interrupt)
            .map_err(failure)?;
        w.charge(inner.work_used())?;
        memory.add(inner.retained_bytes())?;
        let vertices = outer.mesh().vertices.len() + inner.mesh().vertices.len();
        let side_triangles = outer.mesh().face_indices.len() + inner.mesh().face_indices.len();
        let caps = outer.rim_count(ShellEnd::Start)
            + outer.rim_count(ShellEnd::End)
            + inner.rim_count(ShellEnd::Start)
            + inner.rim_count(ShellEnd::End);
        if side_triangles + caps > options.triangle_limit {
            return Err(Error::Limit);
        }
        let mut mesh = DisplayMesh {
            vertices: allocate(vertices, &mut memory, &mut w)?,
            indices: allocate((side_triangles + caps) * 3, &mut memory, &mut w)?,
            face_indices: allocate(side_triangles + caps, &mut memory, &mut w)?,
        };
        w.charge(vertices + side_triangles * 4)?;
        mesh.vertices.extend(&outer.mesh().vertices);
        mesh.vertices.extend(&inner.mesh().vertices);
        mesh.indices.extend(&outer.mesh().indices);
        mesh.face_indices.extend(&outer.mesh().face_indices);
        let offset = outer.mesh().vertices.len() / 3;
        let inner_face = 2 + 4 * self.outer().spec().starts;
        for (tri, face) in inner
            .mesh()
            .indices
            .as_chunks::<3>()
            .0
            .iter()
            .zip(&inner.mesh().face_indices)
        {
            mesh.indices.extend([
                tri[0] + offset as u32,
                tri[2] + offset as u32,
                tri[1] + offset as u32,
            ]);
            mesh.face_indices.push(inner_face + face - 2);
        }
        let mut chord = outer
            .source_chord_bound_mm()
            .max(inner.source_chord_bound_mm());
        for (end, source) in [(ShellEnd::Start, &bottom), (ShellEnd::End, &top)] {
            let outer_ring = ring(&outer, end, 0, &mut memory, &mut w)?;
            let inner_ring = ring(&inner, end, offset, &mut memory, &mut w)?;
            chord = chord.max(cap(
                &mut mesh,
                source,
                &outer_ring,
                &inner_ring,
                options,
                &mut w,
            )?);
            // Ring arenas are small and only two are live. Their retained byte charges
            // remain conservative until request completion; no forgotten live allocation.
        }
        let source_volume = self
            .volume(
                QueryOptions {
                    volume_error_mm3: 1e-7,
                    maximum_cells: 10,
                },
                interrupt,
            )?
            .material_mm3;
        w.charge(512)?;
        let literal_volume = literal_volume(&mesh, &mut w)?;
        w.charge(1)?;
        Ok(Candidate {
            source: *self,
            mesh,
            chord,
            normal: outer
                .source_normal_angle_bound_rad()
                .max(inner.source_normal_angle_bound_rad()),
            work: w.used,
            source_volume,
            literal_volume,
        })
    }
}
