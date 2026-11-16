//! Exact original owner/parameter sharing; angles steer sampling only.
use super::*;
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(super) enum SegmentKey {
    External(usize, usize),
    U(usize, usize, usize),
    V(usize, usize, usize),
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(super) enum VertexKey {
    Original(usize),
    EdgeKnot(usize, usize),
    TensorKnot(usize, usize, usize),
    Sample(SegmentKey, u64),
    Interior(usize, usize, usize, [u64; 2]),
    Cap(usize),
}
#[derive(Clone, Copy)]
pub(super) struct Sample {
    pub t: f64,
    pub vertex: u32,
    pub ray: Option<u64>,
}
pub(super) struct Segment {
    pub points: Vec<Sample>,
}
pub(super) struct Vertex {
    pub local: [f64; 3],
    pub point: [f64; 3],
}
struct Rays {
    values: [u64; 8],
    len: usize,
}
pub(super) struct Registry {
    pub mesh: DisplayMesh,
    pub vertices: Vec<Vertex>,
    pub segments: BTreeMap<SegmentKey, Segment>,
    ids: BTreeMap<VertexKey, u32>,
    extras: BTreeMap<SegmentKey, Rays>,
    pub rays: usize,
    pub triangle_limit: usize,
    deflection: f64,
}
pub(super) fn reserve<T>(
    values: &mut Vec<T>,
    additional: usize,
    budget: &mut Budget<'_>,
) -> Result<()> {
    let count = values.len().checked_add(additional).ok_or(Error::Limit)?;
    if count <= values.capacity() {
        return Ok(());
    }
    let old = values.capacity();
    let next = old.max(16).checked_mul(2).ok_or(Error::Limit)?.max(count);
    hold(
        budget,
        next.checked_mul(size_of::<T>()).ok_or(Error::Limit)?,
    )?;
    values
        .try_reserve_exact(next - values.len())
        .map_err(|_| Error::Limit)?;
    if values.capacity() > next {
        hold(budget, (values.capacity() - next) * size_of::<T>())?;
    }
    budget.bytes -= old * size_of::<T>();
    Ok(())
}
fn entry<K: Ord, V>(map: &BTreeMap<K, V>, budget: &mut Budget<'_>) -> Result<()> {
    if map.is_empty() {
        hold(budget, 4096)?;
    } // first root may have only one occupied slot
    // Non-root pinned BTree nodes contain at least5 of11 slots. Three inline
    // entry sizes plus64 bytes dominate slot and child-pointer occupancy;
    // the separately charged4096 covers the exceptional sparse root.
    hold(budget, 3 * (size_of::<K>() + size_of::<V>()) + 64)
}
fn curve(
    source: &OwnedMaterial,
    controls: [usize; 12],
    t: f64,
    budget: &mut Budget<'_>,
) -> Result<[I; 3]> {
    budget.step(1024)?;
    let weights = basis::<12>(t);
    let anchor = interval_point(source, controls[0]);
    let mut point = anchor;
    for (i, w) in weights.into_iter().enumerate() {
        let p = interval_point(source, controls[i]);
        for k in 0..3 {
            point[k] = point[k] + w * (p[k] - anchor[k]);
        }
    }
    Ok(point)
}
pub(super) fn root_point(
    source: &OwnedMaterial,
    t: usize,
    root: [usize; 2],
    uv: [f64; 2],
    budget: &mut Budget<'_>,
) -> Result<[I; 3]> {
    budget.step(4096)?;
    let tensor = &source.tensors()[t];
    let height = 11 * tensor.spans[1] + 1;
    let weights = uv.map(basis::<12>);
    let anchor = interval_point(source, tensor.poles[11 * root[0] * height + 11 * root[1]]);
    let mut point = anchor;
    for i in 0..12 {
        for j in 0..12 {
            let w = weights[0][i] * weights[1][j];
            let p = interval_point(
                source,
                tensor.poles[(11 * root[0] + i) * height + 11 * root[1] + j],
            );
            for k in 0..3 {
                point[k] = point[k] + w * (p[k] - anchor[k]);
            }
        }
    }
    Ok(point)
}
fn float_curve(source: &OwnedMaterial, controls: [usize; 12], t: f64) -> [f64; 3] {
    let mut points = controls.map(|i| source.0.points[i].coordinates.each_ref().map(|s| s.value));
    for length in (1..12).rev() {
        for i in 0..length {
            points[i] = std::array::from_fn(|k| (1. - t) * points[i][k] + t * points[i + 1][k]);
        }
    }
    points[0]
}
pub(super) fn float_root(
    source: &OwnedMaterial,
    t: usize,
    root: [usize; 2],
    uv: [f64; 2],
) -> [f64; 3] {
    let weights = uv.map(|value| {
        let mut out = [0.; 12];
        out[0] = 1.;
        for degree in 1..12 {
            for i in (0..=degree).rev() {
                out[i] = if i == 0 {
                    (1. - value) * out[i]
                } else if i == degree {
                    value * out[i - 1]
                } else {
                    (1. - value) * out[i] + value * out[i - 1]
                };
            }
        }
        out
    });
    let tensor = &source.tensors()[t];
    let height = 11 * tensor.spans[1] + 1;
    let mut point = [0.; 3];
    for i in 0..12 {
        for j in 0..12 {
            let p = &source.0.points[tensor.poles[(11 * root[0] + i) * height + 11 * root[1] + j]];
            for (k, value) in point.iter_mut().enumerate() {
                *value += weights[0][i] * weights[1][j] * p.coordinates[k].value;
            }
        }
    }
    point
}
fn ray_angle(angle: f64) -> f64 {
    let tau = std::f64::consts::TAU;
    let value = angle.rem_euclid(tau);
    let rounded = (value * 2f64.powi(40)).round() * 2f64.powi(-40);
    if rounded < 1e-11 || tau - rounded < 1e-11 {
        0.
    } else {
        rounded
    }
}
pub(super) fn snap(t: f64) -> f64 {
    (t * 2f64.powi(26)).round() * 2f64.powi(-26)
}
impl Registry {
    pub fn new(rays: usize, triangle_limit: usize, deflection: f64) -> Self {
        Self {
            mesh: DisplayMesh {
                vertices: Vec::new(),
                indices: Vec::new(),
                face_indices: Vec::new(),
            },
            vertices: Vec::new(),
            segments: BTreeMap::new(),
            ids: BTreeMap::new(),
            extras: BTreeMap::new(),
            rays,
            triangle_limit,
            deflection,
        }
    }
    pub fn vertex(
        &mut self,
        source: &OwnedMaterial,
        key: VertexKey,
        local: [I; 3],
        budget: &mut Budget<'_>,
    ) -> Result<u32> {
        budget.step(32)?;
        if let Some(&id) = self.ids.get(&key) {
            return Ok(id);
        }
        let world = source.map_box_mm(local.map(bounds))?.map(interval);
        let vertex = Vertex {
            local: local.map(I::midpoint),
            point: world.map(I::midpoint),
        };
        let id = self.vertices.len();
        if id >= u32::MAX as usize {
            return Err(Error::Limit);
        }
        entry(&self.ids, budget)?;
        reserve(&mut self.vertices, 1, budget)?;
        reserve(&mut self.mesh.vertices, 3, budget)?;
        self.mesh.vertices.extend(vertex.point);
        self.vertices.push(vertex);
        self.ids.insert(key, id as u32);
        Ok(id as u32)
    }
    pub fn triangle(
        &mut self,
        vertices: [u32; 3],
        face: usize,
        budget: &mut Budget<'_>,
    ) -> Result<()> {
        budget.step(32)?;
        if self.mesh.face_indices.len() == self.triangle_limit {
            return Err(Error::Limit);
        }
        reserve(&mut self.mesh.indices, 3, budget)?;
        reserve(&mut self.mesh.face_indices, 1, budget)?;
        self.mesh.indices.extend(vertices);
        self.mesh.face_indices.push(face as u32);
        Ok(())
    }
    fn endpoint(
        &mut self,
        source: &OwnedMaterial,
        key: VertexKey,
        budget: &mut Budget<'_>,
    ) -> Result<u32> {
        let point = match key {
            VertexKey::Original(v) => source.vertices()[v].point,
            VertexKey::EdgeKnot(edge, knot) => {
                source.curves()[source.edges()[edge].curve].poles[11 * knot]
            }
            VertexKey::TensorKnot(t, u, v) => {
                let tensor = &source.tensors()[t];
                tensor.poles[11 * u * (11 * tensor.spans[1] + 1) + 11 * v]
            }
            _ => unreachable!("source endpoint owner"),
        };
        self.vertex(source, key, interval_point(source, point), budget)
    }
    pub fn prepare_rays(&mut self, source: &OwnedMaterial, budget: &mut Budget<'_>) -> Result<()> {
        for (t, tensor) in source.tensors().iter().enumerate() {
            for u in 0..tensor.spans[0] {
                for v in 0..tensor.spans[1] {
                    let mut rays = [0u64; 4];
                    for (i, [a, b]) in [[u, v], [u + 1, v], [u + 1, v + 1], [u, v + 1]]
                        .into_iter()
                        .enumerate()
                    {
                        let key = tensor_vertex(source, t, a, b);
                        let vertex = self.endpoint(source, key, budget)?;
                        let p = self.vertices[vertex as usize].local;
                        rays[i] = ray_angle(p[1].atan2(p[0])).to_bits();
                    }
                    for segment in root_sides(source, t, [u, v]) {
                        if !self.extras.contains_key(&segment) {
                            entry(&self.extras, budget)?;
                            self.extras.insert(
                                segment,
                                Rays {
                                    values: [0; 8],
                                    len: 0,
                                },
                            );
                        }
                        let values = self
                            .extras
                            .get_mut(&segment)
                            .expect("inserted original segment");
                        for ray in rays {
                            if !values.values[..values.len].contains(&ray) {
                                assert!(
                                    values.len < 8,
                                    "a manifold original segment meets at most two root cells"
                                );
                                values.values[values.len] = ray;
                                values.len += 1;
                            }
                        }
                    }
                }
            }
        }
        Ok(())
    }
    pub fn ensure(
        &mut self,
        source: &OwnedMaterial,
        key: SegmentKey,
        budget: &mut Budget<'_>,
    ) -> Result<()> {
        budget.step(32)?;
        if self.segments.contains_key(&key) {
            return Ok(());
        }
        let (controls, endpoints) = segment(source, key);
        let mut candidates = budget.allocate(
            (self.rays + 8)
                .checked_mul(32)
                .and_then(|n| n.checked_add(35))
                .ok_or(Error::Limit)?,
        )?;
        candidates.extend([(0., None), (0.5, None), (1., None)]);
        let step = std::f64::consts::TAU / self.rays as f64;
        let mut directions = budget.allocate(self.rays + 8)?;
        directions.extend((0..self.rays).map(|i| ray_angle(i as f64 * step).to_bits()));
        if let Some(extra) = self.extras.get(&key) {
            directions.extend_from_slice(&extra.values[..extra.len]);
        }
        directions.sort_unstable();
        directions.dedup();
        let mut previous = 0.;
        let first = float_curve(source, controls, 0.);
        let mut angle = first[1].atan2(first[0]);
        for part in 1..=32 {
            budget.step(directions.len() * 16 + 1024)?;
            let next = part as f64 / 32.;
            let point = float_curve(source, controls, next);
            let raw = point[1].atan2(point[0]);
            let delta = (raw - angle).sin().atan2((raw - angle).cos());
            let end = angle + delta;
            for &bits in &directions {
                let theta = f64::from_bits(bits);
                let lo = ((angle.min(end) - theta) / std::f64::consts::TAU).ceil() as i64;
                let hi = ((angle.max(end) - theta) / std::f64::consts::TAU).floor() as i64;
                if lo > hi {
                    continue;
                }
                let direction = [theta.cos(), theta.sin()];
                let f = |t| {
                    let p = float_curve(source, controls, t);
                    direction[0] * p[1] - direction[1] * p[0]
                };
                let (mut a, mut b) = (previous, next);
                let mut fa = f(a);
                let fb = f(b);
                if fa != 0. && fb != 0. && fa.is_sign_positive() == fb.is_sign_positive() {
                    continue;
                }
                for _ in 0..42 {
                    budget.step(1024)?;
                    let mid = a.midpoint(b);
                    let fm = f(mid);
                    if fa == 0. || fm == 0. || fa.is_sign_positive() != fm.is_sign_positive() {
                        b = mid;
                    } else {
                        a = mid;
                        fa = fm;
                    }
                }
                let t = snap(a.midpoint(b));
                let p = float_curve(source, controls, t);
                // Endpoint angular rounding must not create an almost coincident
                // steering vertex. Original endpoints remain literal, and the
                // complete triangle proof still controls admission.
                if t > 1e-10 && t < 1. - 1e-10 && direction[0] * p[0] + direction[1] * p[1] > 0. {
                    candidates.push((t, Some(bits)));
                }
            }
            previous = next;
            angle = end;
        }
        budget.release(&directions);
        candidates.sort_by(|a, b| a.0.total_cmp(&b.0));
        candidates.dedup_by(|a, b| {
            if a.0 == b.0 {
                if b.1.is_none() {
                    b.1 = a.1;
                }
                true
            } else {
                false
            }
        });
        budget.step(2048)?;
        let mut curvature = 0f64;
        for i in 0..10 {
            let p = [i, i + 1, i + 2].map(|j| interval_point(source, controls[j]));
            let d: [I; 3] = std::array::from_fn(|k| {
                (p[2][k] - I::point(2.) * p[1][k] + p[0][k]) * I::point(110.)
            });
            let world = std::array::from_fn(|k| {
                (0..3).fold(I::point(0.), |sum, j| sum + source.0.linear[j][k] * d[j])
            });
            curvature = curvature.max(cells::l1(world));
        }
        let mut index = 0;
        while index + 1 < candidates.len() {
            budget.step(32)?;
            let width = I::point(candidates[index + 1].0) - I::point(candidates[index].0);
            let error = I::point(curvature) * width.square() / I::point(8.);
            if error.hi > self.deflection * 0.25 {
                let mid = snap(candidates[index].0.midpoint(candidates[index + 1].0));
                if mid <= candidates[index].0 || mid >= candidates[index + 1].0 {
                    return Err(Error::Unsupported);
                }
                if candidates.len() > self.triangle_limit {
                    return Err(Error::Limit);
                }
                reserve(&mut candidates, 1, budget)?;
                candidates.insert(index + 1, (mid, None));
            } else {
                index += 1;
            }
        }
        let mut selected = budget.allocate(candidates.len())?;
        for &(t, ray) in &candidates {
            let vertex = if t == 0. {
                self.endpoint(source, endpoints[0], budget)?
            } else if t == 1. {
                self.endpoint(source, endpoints[1], budget)?
            } else {
                self.vertex(
                    source,
                    VertexKey::Sample(key, t.to_bits()),
                    curve(source, controls, t, budget)?,
                    budget,
                )?
            };
            let ray = ray.or_else(|| {
                let p = self.vertices[vertex as usize].local;
                let actual = ray_angle(p[1].atan2(p[0]));
                let uniform = ray_angle((actual / step).round() * step);
                if (actual - uniform).abs() < 1e-11 {
                    return Some(uniform.to_bits());
                }
                self.extras.get(&key).and_then(|extra| {
                    extra.values[..extra.len]
                        .iter()
                        .copied()
                        .find(|bits| (f64::from_bits(*bits) - actual).abs() < 1e-11)
                })
            });
            selected.push(Sample { t, vertex, ray });
        }
        budget.release(&candidates);
        entry(&self.segments, budget)?;
        self.segments.insert(key, Segment { points: selected });
        Ok(())
    }
}
fn face_edges(source: &OwnedMaterial, t: usize) -> [usize; 4] {
    let face = source
        .faces()
        .iter()
        .find(|f| matches!(f.support,Support::Tensor(i)if i==t))
        .expect("sealed unique tensor face");
    side_edges(face)
}
fn edge_vertex(source: &OwnedMaterial, edge: usize, knot: usize) -> VertexKey {
    let e = &source.edges()[edge];
    let count = source.curves()[e.curve].spans;
    if knot == 0 {
        VertexKey::Original(e.vertices[0])
    } else if knot == count {
        VertexKey::Original(e.vertices[1])
    } else {
        VertexKey::EdgeKnot(edge, knot)
    }
}
fn tensor_vertex(source: &OwnedMaterial, t: usize, u: usize, v: usize) -> VertexKey {
    let tensor = &source.tensors()[t];
    let edges = face_edges(source, t);
    if u == 0 {
        edge_vertex(source, edges[0], v)
    } else if u == tensor.spans[0] {
        edge_vertex(source, edges[2], v)
    } else if v == 0 {
        edge_vertex(source, edges[3], u)
    } else if v == tensor.spans[1] {
        edge_vertex(source, edges[1], u)
    } else {
        VertexKey::TensorKnot(t, u, v)
    }
}
fn segment(source: &OwnedMaterial, key: SegmentKey) -> ([usize; 12], [VertexKey; 2]) {
    match key {
        SegmentKey::External(edge, span) => {
            let curve = &source.curves()[source.edges()[edge].curve];
            (
                std::array::from_fn(|i| curve.poles[11 * span + i]),
                [
                    edge_vertex(source, edge, span),
                    edge_vertex(source, edge, span + 1),
                ],
            )
        }
        SegmentKey::U(t, knot, span) => {
            let tensor = &source.tensors()[t];
            (
                std::array::from_fn(|i| {
                    tensor.poles[11 * knot * (11 * tensor.spans[1] + 1) + 11 * span + i]
                }),
                [
                    tensor_vertex(source, t, knot, span),
                    tensor_vertex(source, t, knot, span + 1),
                ],
            )
        }
        SegmentKey::V(t, span, knot) => {
            let tensor = &source.tensors()[t];
            (
                std::array::from_fn(|i| {
                    tensor.poles[(11 * span + i) * (11 * tensor.spans[1] + 1) + 11 * knot]
                }),
                [
                    tensor_vertex(source, t, span, knot),
                    tensor_vertex(source, t, span + 1, knot),
                ],
            )
        }
    }
}
pub(super) fn root_sides(source: &OwnedMaterial, t: usize, [u, v]: [usize; 2]) -> [SegmentKey; 4] {
    let spans = source.tensors()[t].spans;
    let e = face_edges(source, t);
    [
        if u == 0 {
            SegmentKey::External(e[0], v)
        } else {
            SegmentKey::U(t, u, v)
        },
        if v + 1 == spans[1] {
            SegmentKey::External(e[1], u)
        } else {
            SegmentKey::V(t, u, v + 1)
        },
        if u + 1 == spans[0] {
            SegmentKey::External(e[2], v)
        } else {
            SegmentKey::U(t, u + 1, v)
        },
        if v == 0 {
            SegmentKey::External(e[3], u)
        } else {
            SegmentKey::V(t, u, v)
        },
    ]
}
