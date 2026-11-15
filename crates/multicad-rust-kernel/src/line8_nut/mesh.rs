//! One indexed shell, retaining original thread rim owners and sharing every
//! circle/roof intersection with its incident planar face. No proximity joins.
use super::*;
use crate::{
    authored_internal_treatment::{ShellEnd, ThreadShell},
    display::{DisplayMesh, Tessellation},
    manufacturing::{self, MeshProof},
};
use spade::Point2;
use std::collections::BTreeMap;
#[derive(Debug)]
pub struct Candidate {
    source: Nut,
    mesh: DisplayMesh,
    proof: MeshProof,
    chord: f64,
    normal: f64,
    work: usize,
    bytes: usize,
}
impl Candidate {
    pub fn source(&self) -> &Nut {
        &self.source
    }
    pub fn mesh(&self) -> &DisplayMesh {
        &self.mesh
    }
    pub fn into_mesh(self) -> DisplayMesh {
        self.mesh
    }
    pub fn connectivity(&self) -> MeshProof {
        self.proof
    }
    pub fn source_distance_mm(&self) -> f64 {
        self.chord
    }
    pub fn normal_chord(&self) -> f64 {
        self.normal
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_bytes(&self) -> usize {
        self.bytes
    }
}
struct Build<'a> {
    mesh: DisplayMesh,
    work: usize,
    maximum: usize,
    bytes: usize,
    source_bytes: usize,
    chord: f64,
    normal: f64,
    q: Tessellation,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
fn norm(a: [I; 3]) -> I {
    a.into_iter()
        .fold(point(0.), |s, v| {
            let lo = if v.lo <= 0. && v.hi >= 0. {
                0.
            } else {
                v.lo.abs().min(v.hi.abs()).powi(2).next_down().max(0.)
            };
            s + I::new(lo, v.magnitude().powi(2).next_up())
        })
        .sqrt()
}
fn unit(a: [I; 3]) -> Result<[I; 3], Error> {
    let length = norm(a);
    if !length.lo.is_finite() || length.lo <= 0. {
        return Err(Error::Unresolved);
    }
    Ok(a.map(|v| v / length))
}
fn at(v: &[f64], i: usize) -> [I; 3] {
    std::array::from_fn(|k| point(v[3 * i + k]))
}
fn facet(v: &[f64], ids: [usize; 3]) -> Result<[I; 3], Error> {
    let points = ids.map(|id| std::array::from_fn(|i| v[id * 3 + i]));
    crate::literal_triangle::normal_bounds(points)
        .map(|n| n.map(|v| I::new(v[0], v[1])))
        .map_err(|_| Error::Unresolved)
}

fn allowed(q: Tessellation) -> f64 {
    let a = point(q.angle) / point(2.);
    (point(2.) * (a - a * a * a / point(6.))).lo
}
impl Build<'_> {
    fn charge(&mut self, n: usize) -> Result<(), Error> {
        checkpoint(self.interrupt)?;
        self.work = self.work.checked_add(n).ok_or(Error::Limit)?;
        if self.work > self.maximum {
            return Err(Error::Limit);
        }
        Ok(())
    }
    fn observe(&mut self, extra: usize) -> Result<(), Error> {
        let bytes = self.source_bytes
            + self.mesh.vertices.capacity() * 8
            + self.mesh.indices.capacity() * 4
            + self.mesh.face_indices.capacity() * 4
            + extra
            + 16384;
        self.bytes = self.bytes.max(bytes);
        if bytes > 64 * 1024 * 1024 {
            return Err(Error::Limit);
        }
        Ok(())
    }
    fn vertex(&mut self, exact: [I; 3]) -> Result<usize, Error> {
        self.charge(128)?;
        let p = exact.map(I::midpoint);
        if p.iter().any(|v| !v.is_finite() || v.abs() > 1e7) {
            return Err(Error::Limit);
        }
        let error = norm(std::array::from_fn(|i| point(p[i]) - exact[i])).hi;
        self.chord = self.chord.max(error);
        let id = self.mesh.vertices.len() / 3;
        self.mesh
            .vertices
            .try_reserve(3)
            .map_err(|_| Error::Limit)?;
        self.mesh.vertices.extend(p);
        self.observe(0)?;
        Ok(id)
    }
    fn triangle(&mut self, mut ids: [usize; 3], face: u32, normal: [I; 3]) -> Result<(), Error> {
        self.charge(crate::literal_triangle::NORMAL_WORK + 256)?;
        if self.mesh.face_indices.len() >= self.q.triangle_limit {
            return Err(Error::Limit);
        }
        let mut actual = facet(&self.mesh.vertices, ids)?;
        let dot = actual
            .into_iter()
            .zip(normal)
            .fold(point(0.), |sum, (a, b)| sum + a * b);
        if dot.hi < 0. {
            ids.swap(0, 1);
            actual = facet(&self.mesh.vertices, ids)?;
        } else if dot.lo <= 0. {
            return Err(Error::Unresolved);
        }
        let error = norm(std::array::from_fn(|i| actual[i] - normal[i])).hi;
        if error > allowed(self.q) {
            return Err(Error::Unresolved);
        }
        self.normal = self.normal.max(error);
        self.mesh.indices.try_reserve(3).map_err(|_| Error::Limit)?;
        self.mesh
            .face_indices
            .try_reserve(1)
            .map_err(|_| Error::Limit)?;
        self.mesh.indices.extend(ids.map(|i| i as u32));
        self.mesh.face_indices.push(face);
        self.observe(0)
    }
    fn cap(
        &mut self,
        mut rings: Vec<Vec<usize>>,
        projection: [usize; 2],
        face: u32,
        normal: [I; 3],
    ) -> Result<(), Error> {
        let n = rings.iter().map(Vec::len).sum::<usize>();
        self.charge(
            n.checked_mul(n)
                .and_then(|v| v.checked_mul(16))
                .ok_or(Error::Limit)?,
        )?;
        self.observe(n * 1024)?;
        for (i, ring) in rings.iter_mut().enumerate() {
            if ring.len() < 3 {
                return Err(Error::Unresolved);
            }
            let mut area = point(0.);
            for j in 0..ring.len() {
                let a = at(&self.mesh.vertices, ring[j]);
                let b = at(&self.mesh.vertices, ring[(j + 1) % ring.len()]);
                area = area + a[projection[0]] * b[projection[1]]
                    - a[projection[1]] * b[projection[0]];
            }
            if area.lo <= 0. && area.hi >= 0. {
                return Err(Error::Unresolved);
            }
            if (area.hi < 0.) != (i > 0) {
                ring.reverse();
            }
        }
        let mut ids = Vec::with_capacity(n);
        let mut points = Vec::with_capacity(n);
        let mut edges = Vec::with_capacity(n);
        for ring in &rings {
            let offset = ids.len();
            for (i, &id) in ring.iter().enumerate() {
                let p = at(&self.mesh.vertices, id);
                ids.push(id);
                points.push(Point2::new(p[projection[0]].lo, p[projection[1]].lo));
                edges.push([offset + i, offset + (i + 1) % ring.len()]);
            }
        }
        let triangles = crate::display::triangulate_region(points, &edges, rings.len() - 1)
            .map_err(|_| Error::Unresolved)?;
        for t in triangles.triangles {
            self.triangle(t.map(|i| ids[i]), face, normal)?;
        }
        Ok(())
    }
}
fn rim(shell: &ThreadShell, end: ShellEnd) -> Result<Vec<usize>, Error> {
    let mut adjacent = BTreeMap::<usize, Vec<usize>>::new();
    for i in 0..shell.rim_count(end) {
        let [a, b] = shell
            .rim_edge(end, i)
            .ok_or(Error::Unresolved)?
            .vertex_indices;
        adjacent.entry(a).or_default().push(b);
        adjacent.entry(b).or_default().push(a);
    }
    if adjacent.values().any(|v| v.len() != 2) {
        return Err(Error::Unresolved);
    }
    let start = *adjacent.keys().next().ok_or(Error::Unresolved)?;
    let mut out = vec![start];
    let mut previous = usize::MAX;
    let mut current = start;
    loop {
        let pair = &adjacent[&current];
        let next = if pair[0] != previous {
            pair[0]
        } else {
            pair[1]
        };
        if next == start {
            break;
        }
        if out.len() >= adjacent.len() {
            return Err(Error::Unresolved);
        }
        out.push(next);
        previous = current;
        current = next;
    }
    if out.len() != adjacent.len() {
        return Err(Error::Unresolved);
    }
    Ok(out)
}
struct Circle {
    angle: I,
    xy: [I; 2],
    top: I,
}
impl Nut {
    /// Candidate establishes original source correspondence and one closed
    /// indexed component. Whole literal embedding remains a separate obligation.
    pub fn mesh_candidate(
        &self,
        q: Tessellation,
        maximum_work: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Candidate, Error> {
        checkpoint(interrupt)?;
        q.check().map_err(|_| Error::Invalid)?;
        if maximum_work == 0 || maximum_work > 64_000_000 {
            return Err(Error::Limit);
        }
        if self.thread.spec().starts != 1 {
            return Err(Error::Unsupported);
        }
        let inner_q = Tessellation {
            deflection: q.deflection * 0.9,
            angle: (q.angle * 0.9).max(0.01),
            ..q
        };
        let shell = self
            .thread
            .tessellate_prefix_shell_shifted(
                self.relief.start_mm,
                inner_q,
                maximum_work.min(16_000_000),
                interrupt,
            )
            .map_err(|e| match e {
                crate::authored_flat_thread::Error::Cancelled => Error::Cancelled,
                crate::authored_flat_thread::Error::Limit => Error::Limit,
                _ => Error::Unresolved,
            })?;
        let lower = rim(&shell, ShellEnd::Start)?;
        let upper = rim(&shell, ShellEnd::End)?;
        let work = shell.work_used();
        let retained = shell.retained_bytes();
        let old_chord = shell.source_chord_bound_mm();
        let old_normal = shell.source_normal_angle_bound_rad();
        let mut literal = shell.into_mesh();
        let mut build = Build {
            mesh: DisplayMesh {
                vertices: vec![],
                indices: vec![],
                face_indices: vec![],
            },
            work,
            maximum: maximum_work,
            bytes: 0,
            source_bytes: 2 * self.retained_bytes()
                + retained
                + (lower.capacity() + upper.capacity()) * 8,
            chord: 0.,
            normal: 0.,
            q,
            interrupt,
        };
        build.observe(0)?;
        let m = self.phase_rotation.map(|r| r.map(point));
        let cofactor = [
            [m[1][1], -m[1][0], point(0.)],
            [-m[0][1], m[0][0], point(0.)],
            [point(0.), point(0.), self.determinant],
        ];
        let mut max_error = 0_f64;
        for &p in literal.vertices.as_chunks::<3>().0 {
            let exact = [
                m[0][0] * point(p[0]) + m[0][1] * point(p[1]),
                m[1][0] * point(p[0]) + m[1][1] * point(p[1]),
                point(p[2]),
            ];
            let id = build.vertex(exact)?;
            let actual = at(&build.mesh.vertices, id);
            max_error = max_error.max(norm(std::array::from_fn(|i| actual[i] - exact[i])).hi);
        }
        let scale = (self.maximum_thread / point(self.thread.maximum_radius_mm())).hi;
        build.chord = (point(scale.max(1.)) * point(old_chord) + point(max_error)).hi;
        let minimum = (self.minimum_thread / point(self.thread.minimum_radius_mm()))
            .lo
            .min(1.);
        let condition = (point(scale.max(1.)) / point(minimum)).hi;
        for (t, face) in literal
            .indices
            .as_chunks::<3>()
            .0
            .iter()
            .zip(&literal.face_indices)
        {
            build.charge(2 * crate::literal_triangle::NORMAL_WORK + 256)?;
            let ids = t.map(|v| v as usize);
            let old = facet(&literal.vertices, ids)?;
            let exact = unit(cofactor.map(|r| {
                r.into_iter()
                    .zip(old)
                    .fold(point(0.), |s, (a, b)| s + a * b)
            }))?;
            let actual = facet(&build.mesh.vertices, ids)?;
            let error = (point(condition) * point(old_normal)
                + norm(std::array::from_fn(|i| actual[i] - exact[i])))
            .hi;
            if error > allowed(q) {
                return Err(Error::Unresolved);
            }
            build.normal = build.normal.max(error);
            build.mesh.indices.extend([t[1], t[0], t[2]]);
            build.mesh.face_indices.push(face + 9);
        }
        literal.vertices.clear();
        literal.vertices.shrink_to_fit();
        drop(literal);
        build.source_bytes -= retained;
        let r = point(self.relief.radius_mm);
        let a = point(self.width[3]);
        let root = (r * r - a * a).sqrt();
        let alpha = (a / root).atan();
        let pi = I::new(
            std::f64::consts::PI.next_down(),
            std::f64::consts::PI.next_up(),
        );
        let half = pi / point(2.);
        let tau = pi * point(2.);
        let cuts = [
            point(0.),
            alpha,
            half,
            pi - alpha,
            pi,
            pi + alpha,
            pi + half,
            tau - alpha,
            tau,
        ];
        let slope = (point(self.height[3]) - point(self.height[2]))
            / (point(self.width[2]) - point(self.width[3]));
        let mut n = 16usize;
        loop {
            let step = tau / point(n as f64);
            let error = (r * (point(1.) + slope) * step * step / point(8.)).hi;
            if error < q.deflection / 4. && step.hi < q.angle / 4. {
                break;
            }
            n = n.checked_mul(2).ok_or(Error::Limit)?;
            if n > 4096 {
                return Err(Error::Limit);
            }
        }
        // Keep all live circle and rim buffers in the same capacity ledger.
        // Eight ceilings add at most16 entries to the proposed full grid.
        let circle_limit = n.checked_add(16).ok_or(Error::Limit)?;
        build.observe(
            circle_limit
                .checked_mul(std::mem::size_of::<Circle>() + 16)
                .ok_or(Error::Limit)?,
        )?;
        let mut circle = Vec::new();
        circle
            .try_reserve_exact(circle_limit)
            .map_err(|_| Error::Limit)?;
        let mut critical = [0usize; 9];
        for segment in 0..8 {
            critical[segment] = circle.len();
            let count = (((cuts[segment + 1].hi - cuts[segment].lo)
                / (std::f64::consts::TAU / n as f64))
                .ceil() as usize)
                .max(1);
            for j in 0..count {
                let angle = cuts[segment]
                    + (cuts[segment + 1] - cuts[segment]) * point(j as f64) / point(count as f64);
                let xy = if j == 0 {
                    match segment {
                        0 => [r, point(0.)],
                        1 => [root, a],
                        2 => [point(0.), r],
                        3 => [-root, a],
                        4 => [-r, point(0.)],
                        5 => [-root, -a],
                        6 => [point(0.), -r],
                        7 => [root, -a],
                        _ => unreachable!(),
                    }
                } else {
                    [r * angle.cos(), r * angle.sin()]
                };
                let region = if segment == 1 || segment == 2 {
                    1
                } else if segment == 5 || segment == 6 {
                    2
                } else {
                    0
                };
                let top = if region == 0 || (j == 0 && [1, 3, 5, 7].contains(&segment)) {
                    point(self.height[3])
                } else {
                    let y = if region == 1 { xy[1] } else { -xy[1] };
                    point(self.height[3]) - (y - a) * slope
                };
                circle.push(Circle { angle, xy, top });
            }
        }
        critical[8] = circle.len();
        let mut bottom_circle = Vec::new();
        let mut top_circle = Vec::new();
        bottom_circle
            .try_reserve_exact(circle.len())
            .map_err(|_| Error::Limit)?;
        top_circle
            .try_reserve_exact(circle.len())
            .map_err(|_| Error::Limit)?;
        build.source_bytes += circle.capacity() * std::mem::size_of::<Circle>()
            + (bottom_circle.capacity() + top_circle.capacity()) * 8;
        build.observe(0)?;
        for node in &circle {
            bottom_circle.push(build.vertex([
                node.xy[0],
                node.xy[1],
                point(self.relief.start_mm),
            ])?);
            top_circle.push(build.vertex([node.xy[0], node.xy[1], node.top])?);
        }
        for i in 0..circle.len() {
            let next = (i + 1) % circle.len();
            let end = if next == 0 { tau } else { circle[next].angle };
            let delta = end - circle[i].angle;
            let angular = I::new(circle[i].angle.lo, end.hi);
            let normal = [-angular.cos(), -angular.sin(), point(0.)];
            let chord = (r * (point(1.) + slope) * delta * delta / point(8.)).hi;
            build.chord = build.chord.max(chord);
            for ids in [
                [bottom_circle[i], bottom_circle[next], top_circle[next]],
                [bottom_circle[i], top_circle[next], top_circle[i]],
            ] {
                build.triangle(ids, 16, normal)?;
            }
        }
        let x = [
            point(-self.bore_from_start_mm),
            point(self.stock.height_mm()) - point(self.bore_from_start_mm),
        ];
        let mut corners = [[0usize; 8]; 2];
        for end in 0..2 {
            for (i, &p) in self.stock.base()[0].iter().enumerate() {
                corners[end][i] = build.vertex([x[end], point(p[0]), point(p[1])])?;
            }
        }
        build.cap(
            vec![
                vec![corners[0][0], corners[1][0], corners[1][1], corners[0][1]],
                lower,
            ],
            [0, 1],
            0,
            [point(0.), point(0.), point(-1.)],
        )?;
        build.cap(
            vec![bottom_circle.clone(), upper],
            [0, 1],
            15,
            [point(0.), point(0.), point(1.)],
        )?;
        for side in [1usize, 2, 6, 7] {
            let a = self.stock.base()[0][side];
            let b = self.stock.base()[0][(side + 1) % 8];
            let normal = unit([
                point(0.),
                point(b[1]) - point(a[1]),
                point(a[0]) - point(b[0]),
            ])?;
            build.cap(
                vec![vec![
                    corners[0][side],
                    corners[1][side],
                    corners[1][(side + 1) % 8],
                    corners[0][(side + 1) % 8],
                ]],
                [0, 2],
                side as u32,
                normal,
            )?;
        }
        for (end, corner) in corners.iter().enumerate() {
            build.cap(
                vec![corner.to_vec()],
                [1, 2],
                8 + end as u32,
                [point(if end == 0 { -1. } else { 1. }), point(0.), point(0.)],
            )?;
        }
        // Source-owned cylinder/roof intersections are reused by every incident
        // face. Critical corners have exactly original Y=+/-topWidth and Z=top.
        let positive = (critical[1]..=critical[3])
            .map(|i| top_circle[i])
            .collect::<Vec<_>>();
        let negative = (critical[5]..=critical[7])
            .rev()
            .map(|i| top_circle[i])
            .collect::<Vec<_>>();
        for (side, arc) in [(3usize, positive), (5usize, negative)] {
            let mut ring = vec![corners[0][3], corners[1][3], corners[1][4]];
            if side == 5 {
                ring = vec![corners[0][6], corners[1][6], corners[1][5]];
            }
            ring.extend(arc);
            ring.push(corners[0][if side == 3 { 4 } else { 5 }]);
            let a = self.stock.base()[0][side];
            let b = self.stock.base()[0][side + 1];
            let normal = unit([
                point(0.),
                point(b[1]) - point(a[1]),
                point(a[0]) - point(b[0]),
            ])?;
            build.cap(vec![ring], [0, 2], side as u32, normal)?;
        }
        let mut right = vec![corners[1][5], corners[1][4]];
        right.extend((0..=critical[1]).rev().map(|i| top_circle[i]));
        right.extend((critical[7]..circle.len()).rev().map(|i| top_circle[i]));
        build.cap(vec![right], [0, 1], 10, [point(0.), point(0.), point(1.)])?;
        let mut left = vec![corners[0][5]];
        left.extend((critical[3]..=critical[5]).rev().map(|i| top_circle[i]));
        left.push(corners[0][4]);
        build.cap(vec![left], [0, 1], 4, [point(0.), point(0.), point(1.)])?;
        if build.chord > q.deflection {
            return Err(Error::Unresolved);
        }
        let vertices = build.mesh.vertices.len() / 3;
        let triangles = build.mesh.face_indices.len();
        build.charge(triangles.checked_mul(256).ok_or(Error::Limit)?)?;
        build.observe(
            triangles * 128
                + vertices * 32
                + circle.capacity() * std::mem::size_of::<Circle>()
                + top_circle.capacity() * 8
                + bottom_circle.capacity() * 8,
        )?;
        let proof = manufacturing::certify(&build.mesh).map_err(|_| Error::Unresolved)?;
        let mut parent = Vec::<u32>::new();
        parent
            .try_reserve_exact(vertices)
            .map_err(|_| Error::Limit)?;
        parent.extend(0..vertices as u32);
        let mut rank = Vec::<u8>::new();
        rank.try_reserve_exact(vertices).map_err(|_| Error::Limit)?;
        rank.resize(vertices, 0);
        build.observe(parent.capacity() * 4 + rank.capacity())?;
        for t in build.mesh.indices.as_chunks::<3>().0 {
            for [a, b] in [[t[0], t[1]], [t[1], t[2]], [t[2], t[0]]] {
                let mut a = a;
                let mut b = b;
                while parent[a as usize] != a {
                    a = parent[a as usize]
                }
                while parent[b as usize] != b {
                    b = parent[b as usize]
                }
                if a != b {
                    if rank[a as usize] < rank[b as usize] {
                        std::mem::swap(&mut a, &mut b)
                    }
                    parent[b as usize] = a;
                    if rank[a as usize] == rank[b as usize] {
                        rank[a as usize] += 1;
                    }
                }
            }
        }
        // Union by rank bounds every traversal by log2(V). The fixed charge
        // covers all six finds per triangle, not merely successful joins.
        build.charge(triangles.checked_mul(6 * 32 * 8).ok_or(Error::Limit)?)?;
        let mut component = None;
        for i in 0..vertices {
            let mut root = i as u32;
            while parent[root as usize] != root {
                root = parent[root as usize]
            }
            if component.is_some_and(|first| first != root) {
                return Err(Error::Unresolved);
            }
            component = Some(root);
        }
        checkpoint(interrupt)?;
        Ok(Candidate {
            source: self.clone(),
            mesh: build.mesh,
            proof,
            chord: build.chord,
            normal: build.normal,
            work: build.work,
            bytes: build.bytes,
        })
    }
}
