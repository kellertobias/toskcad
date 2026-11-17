//! Angular constraints steer a conforming UV triangulation. Admission uses
//! only original-field full-triangle certificates, never sampled angles.
use super::*;
use spade::{
    ConstrainedDelaunayTriangulation, HasPosition, Point2, Triangulation,
    handles::FixedVertexHandle,
};
pub(super) mod registry;
use registry::{Registry, SegmentKey, VertexKey};
const ROOT_NODES: usize = 4096;
#[derive(Clone, Copy)]
struct Node {
    uv: [f64; 2],
    steering: [f64; 2],
    vertex: u32,
    rounding: f64,
    normal: [I; 3],
    ray: Option<(i64, u64)>,
    generation: u8,
}
impl HasPosition for Node {
    type Scalar = f64;
    fn position(&self) -> Point2<f64> {
        Point2::new(self.steering[0], self.steering[1])
    }
}
type Cdt = ConstrainedDelaunayTriangulation<Node>;
pub(super) struct Output {
    pub mesh: DisplayMesh,
    pub connectivity: crate::manufacturing::MeshProof,
    pub position: f64,
    pub angle: f64,
    pub root_triangles: Vec<Vec<usize>>,
    pub caps: [embedding::Cap; 2],
}
fn side_uv(side: usize, t: f64) -> [f64; 2] {
    match side {
        0 => [0., t],
        1 => [t, 1.],
        2 => [1., t],
        _ => [t, 0.],
    }
}
#[derive(Clone, Copy)]
struct Patch<'a> {
    plan: &'a LiteralGeometryPlan,
    t: usize,
    root: [usize; 2],
    normal: &'a normal::Normal,
}
fn root_node(
    patch: &Patch<'_>,
    registry: &Registry,
    uv: [f64; 2],
    vertex: u32,
    ray: Option<(i64, u64)>,
    generation: u8,
    budget: &mut Budget<'_>,
) -> Result<Node> {
    let Patch {
        plan,
        t,
        root,
        normal,
    } = *patch;
    let source = plan.source();
    let bounds = source
        .map_box_mm(registry::root_point(source, t, root, uv, budget)?.map(bounds))?
        .map(interval);
    let p = registry.vertices[vertex as usize].point;
    let rounding = cells::l1(std::array::from_fn(|k| bounds[k] - I::point(p[k])));
    // Exact dyadic affine steering preserves original boundary collinearity:
    // 26 initial UV bits +12 midpoint levels +at most13 shear-product bits.
    // The metric is prepared once per original root, never per vertex.
    let [shear, scale] = normal.steering;
    let steering = [uv[0] + shear * uv[1], scale * uv[1]];
    Ok(Node {
        uv,
        steering,
        vertex,
        rounding,
        normal: plan.cell_normal_point(t, root, uv, budget)?,
        ray,
        generation,
    })
}
fn boundary(
    patch: &Patch<'_>,
    registry: &mut Registry,
    cdt: &mut Cdt,
    budget: &mut Budget<'_>,
) -> Result<()> {
    let Patch { plan, t, root, .. } = *patch;
    let source = plan.source();
    let sides = registry::root_sides(source, t, root);
    for key in sides {
        registry.ensure(source, key, budget)?;
    }
    let mut nodes = Vec::new();
    let mut levels: BTreeMap<(i64, u64), Vec<FixedVertexHandle>> = BTreeMap::new();
    let mut phase = 0.;
    let mut previous = 0.;
    let mut first = true;
    for (side, forward) in [(3, true), (2, true), (1, false), (0, false)] {
        let samples = registry.segments[&sides[side]].points.clone();
        let count = samples.len();
        for index in 0..count - 1 {
            let sample = samples[if forward { index } else { count - index - 1 }];
            let local = registry.vertices[sample.vertex as usize].local;
            let angle = local[1].atan2(local[0]);
            if first {
                phase = angle;
                first = false;
            } else {
                let delta = angle - previous;
                phase += delta.sin().atan2(delta.cos());
            }
            previous = angle;
            let ray = sample.ray.map(|bits| {
                (
                    ((phase - f64::from_bits(bits)) / std::f64::consts::TAU).round() as i64,
                    bits,
                )
            });
            let node = root_node(
                patch,
                registry,
                side_uv(side, sample.t),
                sample.vertex,
                ray,
                0,
                budget,
            )?;
            let handle = cdt.insert(node).map_err(|_| Error::Unsupported)?;
            if cdt.vertex(handle).data().vertex != sample.vertex {
                return Err(Error::Unsupported);
            }
            nodes.push(handle);
            if let Some(level) = ray {
                levels.entry(level).or_default().push(handle);
            }
        }
    }
    for i in 0..nodes.len() {
        let a = nodes[i];
        let b = nodes[(i + 1) % nodes.len()];
        if a == b {
            return Err(Error::Unsupported);
        }
        if !cdt.exists_constraint(a, b) {
            if !cdt.can_add_constraint(a, b) {
                return Err(Error::Unsupported);
            }
            cdt.add_constraint(a, b);
        }
    }
    for handles in levels.values_mut() {
        handles.sort_unstable();
        handles.dedup();
        if let [a, b] = handles.as_slice()
            && a != b
            && !cdt.exists_constraint(*a, *b)
            && cdt.can_add_constraint(*a, *b)
        {
            cdt.add_constraint(*a, *b);
        }
    }
    Ok(())
}
fn on_boundary(a: [f64; 2], b: [f64; 2]) -> bool {
    (0..2).any(|k| (a[k] == 0. && b[k] == 0.) || (a[k] == 1. && b[k] == 1.))
}
fn refine(
    patch: &Patch<'_>,
    registry: &mut Registry,
    cdt: &mut Cdt,
    vertices: [FixedVertexHandle; 3],
    budget: &mut Budget<'_>,
) -> Result<()> {
    let Patch { plan, t, root, .. } = *patch;
    let source = plan.source();
    let mut selected = None;
    for i in 0..3 {
        let a = vertices[i];
        let b = vertices[(i + 1) % 3];
        let na = *cdt.vertex(a).data();
        let nb = *cdt.vertex(b).data();
        if on_boundary(na.uv, nb.uv) {
            continue;
        }
        let distance =
            (na.steering[0] - nb.steering[0]).powi(2) + (na.steering[1] - nb.steering[1]).powi(2);
        if selected.is_none_or(|(_, _, d)| distance > d) {
            selected = Some((a, b, distance));
        }
    }
    let (a, b, _) = selected.ok_or(Error::Unsupported)?;
    let na = *cdt.vertex(a).data();
    let nb = *cdt.vertex(b).data();
    let generation = na.generation.max(nb.generation) + 1;
    if generation > 12 {
        return Err(Error::Unsupported);
    }
    let mut uv = std::array::from_fn(|k| (na.uv[k] + nb.uv[k]) * 0.5);
    let ray = na.ray.filter(|ray| Some(*ray) == nb.ray);
    let constraint = cdt
        .get_edge_from_neighbors(a, b)
        .filter(|e| e.is_constraint_edge())
        .map(|e| e.fix().as_undirected());
    // Only a constrained contour owns an angular level. A free CDT chord
    // may connect two same-level nodes across a curved contour; projecting
    // that chord would repeat an existing contour vertex.
    let mut ray = ray.filter(|_| constraint.is_some());
    if let Some((_, bits)) = ray {
        let theta = f64::from_bits(bits);
        let direction = [theta.cos(), theta.sin()];
        let axis = usize::from((na.uv[0] - nb.uv[0]).abs() >= (na.uv[1] - nb.uv[1]).abs());
        let evaluate = |value| {
            let mut at = uv;
            at[axis] = value;
            let p = registry::float_root(source, t, root, at);
            direction[0] * p[1] - direction[1] * p[0]
        };
        budget.step(4096 * 42)?;
        let (mut lo, mut hi) = (0., 1.);
        let mut flo = evaluate(lo);
        let fhi = evaluate(hi);
        if flo.is_sign_positive() != fhi.is_sign_positive() {
            for _ in 0..40 {
                let middle = lo.midpoint(hi);
                let value = evaluate(middle);
                if value.is_sign_positive() == flo.is_sign_positive() {
                    lo = middle;
                    flo = value;
                } else {
                    hi = middle;
                }
            }
            uv[axis] = registry::snap(lo.midpoint(hi));
        }
    }
    if let Some(edge) = constraint {
        cdt.remove_constraint_edge(edge);
        // A curved angular contour can bend across a different steering
        // constraint. This is only a triangulation proposal: if its two
        // proposed segments conflict, split the original straight chord
        // instead and drop the unearned angular level. The complete original
        // field C0/C1 certificate below remains the sole admission authority.
        if ray.is_some() {
            let [shear, scale] = patch.normal.steering;
            let point = Point2::new(uv[0] + shear * uv[1], scale * uv[1]);
            budget.step(cdt.num_directed_edges() * 128 + 128)?;
            if [na, nb].into_iter().any(|n| {
                cdt.intersects_constraint(Point2::new(n.steering[0], n.steering[1]), point)
            }) {
                uv = std::array::from_fn(|k| (na.uv[k] + nb.uv[k]) * 0.5);
                ray = None;
            }
        }
    }
    if uv == na.uv || uv == nb.uv || uv.into_iter().any(|x| x <= 0. || x >= 1.) {
        #[cfg(test)]
        eprintln!(
            "NPT refine point refused uv={uv:?} a={:?} b={:?} ray={ray:?}",
            na.uv, nb.uv
        );
        return Err(Error::Unsupported);
    }
    let key = VertexKey::Interior(t, root[0], root[1], uv.map(f64::to_bits));
    let local = registry::root_point(source, t, root, uv, budget)?;
    let vertex = registry.vertex(source, key, local, budget)?;
    let node = root_node(patch, registry, uv, vertex, ray, generation, budget)?;
    let old = cdt.num_vertices();
    let added = cdt.insert(node).map_err(|_| Error::Unsupported)?;
    if constraint.is_some() {
        for endpoint in [a, b] {
            if !cdt.can_add_constraint(endpoint, added) {
                #[cfg(test)]
                eprintln!(
                    "NPT refine constraint refused uv={uv:?} a={:?} b={:?} ray={ray:?}",
                    na.uv, nb.uv
                );
                return Err(Error::Unsupported);
            }
            cdt.add_constraint(endpoint, added);
        }
    }
    if cdt.num_vertices() == old {
        #[cfg(test)]
        eprintln!("NPT refine duplicate uv={uv:?}");
        return Err(Error::Unsupported);
    }
    Ok(())
}
fn root_mesh(
    patch: &Patch<'_>,
    registry: &mut Registry,
    face: usize,
    options: LiteralMeshOptions,
    budget: &mut Budget<'_>,
) -> Result<(f64, f64, usize)> {
    let normal = patch.normal;
    #[cfg(test)]
    let (t, root) = (patch.t, patch.root);
    // Bound all delegated CDT vertices/halfedges/faces, their growth peaks,
    // local source caches and certificate map within one live root workspace.
    let workspace = ROOT_NODES * 4096 + 65536;
    hold(budget, workspace)?;
    let mut cdt = Cdt::new();
    boundary(patch, registry, &mut cdt, budget).inspect_err(|error| {
        #[cfg(test)]
        eprintln!("NPT angular boundary error t={t} root={root:?}: {error:?}");
        #[cfg(not(test))]
        let _ = error;
    })?;
    let mut proofs: BTreeMap<[usize; 3], (f64, f64)> = BTreeMap::new();
    loop {
        budget.step(cdt.num_inner_faces() * 32 + 1)?;
        if cdt.num_vertices() > ROOT_NODES {
            return Err(Error::Limit);
        }
        let mut failed = None;
        let mut position = 0f64;
        let mut angle = 0f64;
        for triangle in cdt.inner_faces() {
            let handles = triangle.vertices().map(|v| v.fix());
            let nodes = handles.map(|v| *cdt.vertex(v).data());
            let order = [0usize, 2, 1];
            let mut key = handles.map(|v| v.index());
            key.sort_unstable();
            let value = if let Some(&value) = proofs.get(&key) {
                value
            } else {
                let nodes = order.map(|i| nodes[i]);
                let points = nodes.map(|n| registry.vertices[n.vertex as usize].point);
                let rounding = nodes.into_iter().map(|n| n.rounding).fold(0f64, f64::max);
                let quality = sheared::qualify_cached(
                    normal,
                    nodes.map(|n| n.uv),
                    points,
                    rounding,
                    nodes.map(|n| n.normal),
                    [options.tessellation.deflection, options.tessellation.angle],
                    budget,
                )
                .inspect_err(|error| {
                    #[cfg(test)]
                    eprintln!(
                        "NPT angular triangle error t={t} root={root:?} uv={:?}: {error:?}",
                        nodes.map(|n| n.uv)
                    );
                    #[cfg(not(test))]
                    let _ = error;
                })?;
                let value = (quality.position, quality.angle);
                if value.0 <= options.tessellation.deflection
                    && value.1 <= options.tessellation.angle
                    && proofs.len() < ROOT_NODES * 8
                {
                    proofs.insert(key, value);
                }
                value
            };
            if !(value.0 <= options.tessellation.deflection
                && value.1 <= options.tessellation.angle)
            {
                #[cfg(test)]
                if nodes.iter().any(|n| n.generation > 5) {
                    eprintln!(
                        "NPT angular difficult t={t} root={root:?} count={} uv={:?} C0={} angle={}",
                        cdt.num_vertices(),
                        nodes.map(|n| n.uv),
                        value.0,
                        value.1
                    );
                }
                failed = Some(handles);
                break;
            }
            position = position.max(value.0);
            angle = angle.max(value.1);
        }
        if let Some(vertices) = failed {
            refine(patch, registry, &mut cdt, vertices, budget)?;
            continue;
        }
        let count = cdt.num_inner_faces();
        for face_handle in cdt.inner_faces() {
            let nodes = face_handle.vertices().map(|v| *v.data());
            registry.triangle(
                [nodes[0].vertex, nodes[2].vertex, nodes[1].vertex],
                face,
                budget,
            )?;
        }
        budget.bytes -= workspace;
        return Ok((position, angle, count));
    }
}
pub(super) fn build(
    plan: &LiteralGeometryPlan,
    options: LiteralMeshOptions,
    budget: &mut Budget<'_>,
) -> Result<Output> {
    let source = plan.source();
    let cache = &plan.normals;
    let local = source.local_enclosure_mm();
    let radius = (local[0][0].abs().max(local[0][1].abs()).powi(2)
        + local[1][0].abs().max(local[1][1].abs()).powi(2))
    .sqrt();
    let scale = source.0.linear[0]
        .into_iter()
        .fold(I::point(0.), |sum, x| sum + x.square())
        .sqrt()
        .hi;
    let step = (options.tessellation.angle * 0.125)
        .min((options.tessellation.deflection / (radius * scale)).sqrt());
    let rays = (std::f64::consts::TAU / step).ceil();
    if !rays.is_finite() || rays > options.tessellation.triangle_limit as f64 {
        return Err(Error::Limit);
    }
    let mut registry = Registry::new(
        rays as usize,
        options.tessellation.triangle_limit,
        options.tessellation.deflection,
    );
    registry.prepare_rays(source, budget)?;
    let mut root_triangles = budget.allocate(source.tensors().len())?;
    for tensor in source.tensors() {
        let mut counts = budget.allocate(tensor.spans[0] * tensor.spans[1])?;
        counts.resize(tensor.spans[0] * tensor.spans[1], 0);
        root_triangles.push(counts);
    }
    let mut position = 0f64;
    let mut angle = 0f64;
    for (face_index, face) in source.faces().iter().enumerate() {
        let Support::Tensor(t) = face.support else {
            continue;
        };
        let tensor = &source.tensors()[t];
        for u in 0..tensor.spans[0] {
            for v in 0..tensor.spans[1] {
                let p = root_mesh(
                    &Patch {
                        plan,
                        t,
                        root: [u, v],
                        normal: &cache[t][u * tensor.spans[1] + v],
                    },
                    &mut registry,
                    face_index,
                    options,
                    budget,
                )
                .inspect_err(|error| {
                    #[cfg(test)]
                    eprintln!(
                        "NPT angular root error t={t} root={u},{v}: {error:?} work={} bytes={}",
                        budget.used, budget.peak
                    );
                    #[cfg(not(test))]
                    let _ = error;
                })?;
                position = position.max(p.0);
                angle = angle.max(p.1);
                root_triangles[t][u * tensor.spans[1] + v] = p.2;
            }
        }
    }
    let mut caps = [None, None];
    for (face_index, face) in source.faces().iter().enumerate() {
        let Support::Plane { height } = face.support else {
            continue;
        };
        let original_axis = [
            I::point(0.),
            I::point(0.),
            interval_point(source, height)[2],
        ];
        // The small, explicitly charged outward apex gives a literal fan
        // strictly outside the rounded rim. Source cap correspondence uses
        // the ORIGINAL axis anchor, retaining the complete radial-disk bound.
        let column = cells::l1(source.0.linear[2]);
        let offset = options.tessellation.deflection * 0.001 / column;
        let mut axis = original_axis;
        axis[2] = axis[2] + I::point(if face.forward { offset } else { -offset });
        let anchor = registry.vertex(source, VertexKey::Cap(face_index), axis, budget)?;
        let actual = source.map_box_mm(original_axis.map(bounds))?.map(interval);
        let mut rim = Vec::new();
        let mut cap_edges = Vec::new();
        let a = registry.vertices[anchor as usize].point;
        position = position.max(cells::l1(std::array::from_fn(|k| {
            actual[k] - I::point(a[k])
        })));
        for coedge in &face.coedges {
            let spans = source.curves()[source.edges()[coedge.edge].curve].spans;
            for span_index in 0..spans {
                let span = if coedge.forward {
                    span_index
                } else {
                    spans - span_index - 1
                };
                let key = SegmentKey::External(coedge.edge, span);
                registry.ensure(source, key, budget)?;
                let sample_count = registry.segments[&key].points.len();
                for index in 0..sample_count - 1 {
                    // Copy only two fixed samples. Do not allocate a duplicate
                    // span buffer while extending the retained cap trace.
                    budget.step(64)?;
                    let pair = {
                        let samples = &registry.segments[&key].points;
                        if coedge.forward {
                            [samples[index], samples[index + 1]]
                        } else {
                            [
                                samples[sample_count - index - 1],
                                samples[sample_count - index - 2],
                            ]
                        }
                    };
                    let [a, b] = pair.map(|s| s.vertex);
                    registry::reserve(&mut cap_edges, 1, budget)?;
                    cap_edges.push(LiteralCapEdge {
                        edge: coedge.edge,
                        span,
                        parameters: pair.map(|s| s.t),
                        vertices: [a, b],
                    });
                    registry::reserve(&mut rim, 1, budget)?;
                    rim.push(a);
                    let triangle = [anchor, a, b];
                    let error = cells::cap_angle(
                        source,
                        face.forward,
                        triangle.map(|i| registry.vertices[i as usize].point),
                    );
                    if error
                        .partial_cmp(&options.tessellation.angle)
                        .is_none_or(|o| o.is_gt())
                    {
                        return Err(Error::Unsupported);
                    }
                    angle = angle.max(error);
                    registry.triangle(triangle, face_index, budget)?;
                }
            }
        }
        caps[usize::from(face.forward)] = Some(embedding::Cap {
            face: face_index,
            apex: anchor,
            rim,
            edges: cap_edges,
        });
    }
    if position
        .partial_cmp(&options.tessellation.deflection)
        .is_none_or(|o| o.is_gt())
    {
        return Err(Error::Unsupported);
    }
    let certification_base = budget.bytes;
    let connectivity = crate::manufacturing::certify_counted(&registry.mesh, &mut |n, scratch| {
        budget.step(n)?;
        if let Some(bytes) = scratch {
            let total = certification_base.checked_add(bytes).ok_or(Error::Limit)?;
            if total > budget.limit {
                return Err(Error::Limit);
            }
            budget.peak = budget.peak.max(total);
        }
        Ok(())
    })
    .map_err(|e| match e {
        crate::manufacturing::CountedError::Mesh(_) => Error::Unsupported,
        crate::manufacturing::CountedError::Interrupted(e) => e,
    })?;
    budget.step(1)?;
    Ok(Output {
        mesh: registry.mesh,
        connectivity,
        position,
        angle,
        root_triangles,
        caps: caps.map(|cap| cap.expect("sealed material has one cap per orientation")),
    })
}

#[cfg(test)]
mod probes {
    use super::*;
    #[test]
    fn actual_clipped_root_mesh_probe() {
        let source = super::super::tests::source(true, false);
        let options = LiteralMeshOptions::default();
        let mut budget = Budget {
            used: 0,
            remaining: options.maximum_work,
            bytes: source.retained_budget_bytes() + SCRATCH,
            peak: source.retained_budget_bytes() + SCRATCH,
            limit: options.maximum_bytes,
            interrupt: None,
        };
        let plan = source.geometry_plan(Default::default(), None).unwrap();
        let normal = &plan.normals[5][3];
        let mut registry = Registry::new(
            168,
            options.tessellation.triangle_limit,
            options.tessellation.deflection,
        );
        registry.prepare_rays(&source, &mut budget).unwrap();
        let face = source
            .faces()
            .iter()
            .position(|f| matches!(f.support, Support::Tensor(5)))
            .unwrap();
        let start = std::time::Instant::now();
        let result = root_mesh(
            &Patch {
                plan: &plan,
                t: 5,
                root: [0, 3],
                normal,
            },
            &mut registry,
            face,
            options,
            &mut budget,
        );
        eprintln!(
            "NPT END ROOT {result:?} {:?} work={} peak={}",
            start.elapsed(),
            budget.used,
            budget.peak
        );
        result.unwrap();
    }
}
