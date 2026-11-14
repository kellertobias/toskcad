//! Actual literal manifold and radial degree-one material proof. No authored
//! source/header, fitted thread, display mesh or nominal mass is consulted.
use super::*;
use crate::{
    step::insert::{
        projection::{self, polynomial::P},
        resources::{Memory, Work},
        topology::{self, EdgeGeometry, FaceGeometry, ShellCoedge, ShellEdge, ShellFace},
    },
    tapered_thread::interval::I,
};
const SCRATCH: usize = 256 * 1024;
pub struct MaterialBounds<'a, 'b, 'i> {
    binding: &'a SupportBindings<'b, 'i>,
    outer_component: u8,
    components: Vec<u8>,
    bounds: [[f64; 2]; 3],
    radii: [[f64; 2]; 2],
    minimum_jacobian: f64,
    minimum_cap_derivative: f64,
    work: usize,
    bytes: usize,
}
impl<'a, 'b, 'i> MaterialBounds<'a, 'b, 'i> {
    pub fn binding(&self) -> &'a SupportBindings<'b, 'i> {
        self.binding
    }
    pub fn bounds_mm(&self) -> [[f64; 2]; 3] {
        self.bounds
    }
    pub fn component_radius_mm(&self) -> [[f64; 2]; 2] {
        self.radii
    }
    pub fn minimum_signed_jacobian(&self) -> f64 {
        self.minimum_jacobian
    }
    pub fn minimum_signed_cap_derivative(&self) -> f64 {
        self.minimum_cap_derivative
    }
    pub fn outer_component(&self) -> u8 {
        self.outer_component
    }
    pub fn tensor_component(&self, face: usize) -> Option<u8> {
        self.components.get(face).copied().filter(|&i| i < 2)
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.bytes
    }
}
// Global argument: literal incidence earns one oriented manifold with cyclic
// vertex links and chi=0. Removing its two disjoint annular caps leaves two
// connected sides, each with two Jordan boundary loops; hence both are
// cylinders. Positive radius and actual signed (theta,z) Jacobians (including
// true owned-vertex corner factors) make each side a proper positive covering.
// Every corner has a positive projected wedge, so its cyclic fan has positive
// local degree. The actual cap degree is one; consequently every fiber has one
// point. Strictly separated radial ranges then enclose one embedded material
// region. These hypotheses are earned below, never inferred from source names.
impl<'b, 'i> SupportBindings<'b, 'i> {
    pub fn material_bounds<'a>(
        &'a self,
        options: Options,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<MaterialBounds<'a, 'b, 'i>> {
        if options.retained_budget_bytes == 0 || options.retained_budget_bytes > 64 * 1024 * 1024 {
            return Err(Error::Invalid);
        }
        let mut w = Work::new(options.maximum_work, interrupt).map_err(Error::from)?;
        w.charge(self.work_used()).map_err(Error::from)?;
        let mut memory = Memory::new(self.retained_budget_bytes(), options.retained_budget_bytes)
            .map_err(Error::from)?;
        memory
            .add(SCRATCH + size_of::<MaterialBounds<'_, '_, '_>>())
            .map_err(Error::from)?;
        let p = self.source();
        let heights = self.cap_heights();
        let mut cache = memory
            .allocate(p.points().len(), &mut w)
            .map_err(Error::from)?;
        for point in p.points() {
            let c = point.coordinates();
            if c.len() == 3 {
                let mut xyz = [I::point(0.); 3];
                for i in 0..3 {
                    w.charge(c[i].bounds_work()).map_err(Error::from)?;
                    let [lo, hi] = c[i].bounds();
                    let x = I::new(lo, hi);
                    if !x.valid() {
                        return Err(Error::Unsupported);
                    }
                    xyz[i] = x;
                }
                cache.push(Some(xyz));
            } else {
                cache.push(None);
            }
        }
        manifold(self, &mut memory, &mut w)?;
        let components = components(self, &mut memory, &mut w)?;
        let mut stack = memory.allocate(3 * 32 + 1, &mut w).map_err(Error::from)?;
        let mut height_checked = memory
            .allocate(p.points().len(), &mut w)
            .map_err(Error::from)?;
        height_checked.resize(p.points().len(), false);
        let span_count = self
            .tensors()
            .iter()
            .map(|patch| p.tensors()[patch.tensor_index()].knots_v().len() - 1)
            .sum();
        // The radius and Jacobian proofs use the SAME literal polynomial nets.
        // Decode each net once; cache coefficients without rounding or fitting.
        let mut nets = memory.allocate(span_count, &mut w).map_err(Error::from)?;
        let mut radii = [[f64::INFINITY, 0_f64]; 2];
        let mut bounds = [[f64::INFINITY, f64::NEG_INFINITY]; 3];
        // Radius comes from each actual literal tensor, BEFORE assigning outer/inner.
        for patch in self.tensors() {
            let t = &p.tensors()[patch.tensor_index()];
            let component = components[patch.face_index()] as usize;
            for span in 0..t.knots_v().len() - 1 {
                let mut coordinates = ControlCoordinates {
                    p,
                    cache: &cache,
                    heights,
                    checked: &mut height_checked,
                    bounds: &mut bounds,
                };
                let net = controls(&mut coordinates, t, span, &mut memory, &mut w)?;
                nets.push(net);
                let r = convex_radius(net.data, &mut w)?;
                radii[component][0] = radii[component][0].min(r.lo);
                radii[component][1] = radii[component][1].max(r.hi);
            }
        }
        memory.release(&height_checked);
        drop(height_checked);
        if !(radii[0][0] > radii[1][1] || radii[1][0] > radii[0][1]) {
            // The cheap convex enclosure is sufficient only when the actual
            // components separate. Preserve the stronger polynomial proof for
            // literal nets whose control hulls overlap.
            radii = [[f64::INFINITY, 0.]; 2];
            for (patch, data) in self
                .tensors()
                .iter()
                .flat_map(|patch| {
                    let count = p.tensors()[patch.tensor_index()].knots_v().len() - 1;
                    std::iter::repeat_n(patch, count)
                })
                .zip(&nets)
            {
                let r = projection::radius_squared(data.data, &mut w).map_err(Error::from)?;
                let (r, _) = projection::positive(r, &mut stack, &mut w).map_err(Error::from)?;
                let r = r.sqrt();
                let component = components[patch.face_index()] as usize;
                radii[component][0] = radii[component][0].min(r.lo);
                radii[component][1] = radii[component][1].max(r.hi);
            }
        }
        let outer = if radii[0][0] > radii[1][1] {
            0
        } else if radii[1][0] > radii[0][1] {
            1
        } else {
            return Err(Error::Unsupported);
        };
        let mut minimum = f64::INFINITY;
        let mut net_index = 0;
        let mut decimal_memo = None;
        for patch in self.tensors() {
            let t = &p.tensors()[patch.tensor_index()];
            let face = &p.faces()[patch.face_index()];
            let side = if components[patch.face_index()] == outer {
                1.
            } else {
                -1.
            };
            let same = if face.same_sense() { 1. } else { -1. };
            let collapsed = patch.collapsed_v_edges();
            for span in 0..t.knots_v().len() - 1 {
                let reverse = collapsed[1] && span + 2 == t.knots_v().len();
                let factor = (collapsed[0] && span == 0) || reverse;
                let mut net = nets[net_index];
                net_index += 1;
                if reverse {
                    w.charge(6 * 12).map_err(Error::from)?;
                    for row in &mut net.data {
                        for axis in row {
                            axis.c[..=axis.degree].reverse();
                        }
                    }
                    for row in &mut net.points {
                        row.reverse();
                    }
                }
                let data = net.data;
                let j = if let Some(j) = strip_jacobian(p, &net, factor, &mut w)? {
                    j
                } else {
                    projection::jacobian(data, factor, &mut w).map_err(Error::from)?
                };
                let sign = side * same * if reverse { -1. } else { 1. };
                let j = [
                    j[0].scale(sign, &mut w).map_err(Error::from)?,
                    j[1].scale(sign, &mut w).map_err(Error::from)?,
                    j[2].scale(sign, &mut w).map_err(Error::from)?,
                ];
                // Select arithmetic, never acceptance: almost-stationary angular
                // rows need literal differences before interval conversion.
                let tiny_v = data.iter().all(|row| {
                    let h = [row[0].hull(), row[1].hull()];
                    let scale = h[0].magnitude().max(h[1].magnitude()).max(1.);
                    h.iter()
                        .all(|axis| axis.hi - axis.lo <= 64. * f64::EPSILON * scale)
                });
                let tiny_u = [0, 11].into_iter().any(|i| {
                    (!factor || i != 0)
                        && (0..3).all(|axis| {
                            let a = data[0][axis].c[i];
                            let b = data[1][axis].c[i];
                            let scale = (0..3)
                                .map(|a| {
                                    data[0][a].c[i].magnitude().max(data[1][a].c[i].magnitude())
                                })
                                .fold(1_f64, f64::max);
                            (b - a).magnitude() <= 64. * f64::EPSILON * scale
                        })
                });
                let tiny = tiny_u || tiny_v;
                w.charge(84).map_err(Error::from)?;
                let proof = if tiny {
                    None
                } else {
                    Some(projection::positive(j, &mut stack, &mut w))
                };
                let proof = if matches!(
                    proof,
                    None | Some(Err(crate::step::insert::Error::Unresolved))
                ) {
                    stack.clear();
                    let mut tokens = [[crate::step::insert::exact_cap::LiteralPoint {
                        index: 0,
                        coordinates: [""; 3],
                    }; 12]; 2];
                    for (r, row) in tokens.iter_mut().enumerate() {
                        for (i, point) in row.iter_mut().enumerate() {
                            let index = net.points[r][i];
                            point.index = index;
                            for (axis, token) in point.coordinates.iter_mut().enumerate() {
                                let scalar = p.points()[index].coordinates()[axis];
                                *token = scalar.spelling();
                            }
                        }
                    }
                    if decimal_memo.is_none() {
                        decimal_memo = Some(
                            crate::step::insert::exact_cap::DecimalMemo::new(
                                p.points().len(),
                                &mut memory,
                                &mut w,
                            )
                            .map_err(Error::from)?,
                        );
                    }
                    let [d, da, dd] = crate::step::insert::exact_cap::tensor_differentials(
                        tokens,
                        factor,
                        decimal_memo.as_mut().expect("memo initialized above"),
                        &mut memory,
                        &mut w,
                    )
                    .map_err(Error::from)?;
                    let exact =
                        projection::jacobian_differentials(data[0], d, da, dd, factor, &mut w)
                            .map_err(Error::from)?;
                    let signed = [
                        exact[0].scale(sign, &mut w).map_err(Error::from)?,
                        exact[1].scale(sign, &mut w).map_err(Error::from)?,
                        exact[2].scale(sign, &mut w).map_err(Error::from)?,
                    ];
                    projection::positive(signed, &mut stack, &mut w)
                } else {
                    proof.expect("ordinary interval proof was selected")
                };
                let (j, _) = proof.map_err(Error::from)?;
                minimum = minimum.min(j.lo);
            }
        }
        memory.release(&nets);
        drop(nets);
        if let Some(memo) = decimal_memo {
            memo.release(&mut memory);
        }
        let cap_min = caps(
            self,
            &components,
            outer,
            &cache,
            &mut stack,
            &mut memory,
            &mut w,
        )?;
        memory.release(&cache);
        drop(cache);
        memory.release(&stack);
        drop(stack);
        memory.used -= SCRATCH;
        bounds[2] = [heights[0].bounds()[0], heights[1].bounds()[1]];
        w.charge(1).map_err(Error::from)?;
        Ok(MaterialBounds {
            binding: self,
            outer_component: outer,
            components,
            bounds,
            radii: [radii[outer as usize], radii[1 - outer as usize]],
            minimum_jacobian: minimum,
            minimum_cap_derivative: cap_min,
            work: w.used,
            bytes: memory.used,
        })
    }
}
/// A Bézier tensor is a convex combination of its literal controls. The norm
/// is convex, and every positive separating projection bounds it from below.
fn convex_radius(data: [[P; 3]; 2], w: &mut Work<'_>) -> Result<I> {
    let mut lower = 0_f64;
    let mut upper = 0_f64;
    for row in &data {
        for i in 0..=row[0].degree {
            w.charge(5).map_err(Error::from)?;
            let x = row[0].c[i];
            let y = row[1].c[i];
            let r = (x * x + y * y).sqrt();
            if !r.valid() {
                return Err(Error::Unsupported);
            }
            upper = upper.max(r.hi);
        }
    }
    for [dx, dy] in [
        [1., 0.],
        [1., 1.],
        [0., 1.],
        [-1., 1.],
        [-1., 0.],
        [-1., -1.],
        [0., -1.],
        [1., -1.],
    ] {
        let mut projected = f64::INFINITY;
        for row in &data {
            for i in 0..=row[0].degree {
                w.charge(3).map_err(Error::from)?;
                projected =
                    projected.min((row[0].c[i] * I::point(dx) + row[1].c[i] * I::point(dy)).lo);
            }
        }
        let scale = (I::point(dx * dx) + I::point(dy * dy)).sqrt();
        lower = lower.max((I::point(projected) / scale).lo);
    }
    Ok(I::new(lower, upper))
}
fn point_index(p: &Parsed<'_>, id: u32, w: &mut Work<'_>) -> Result<usize> {
    let (mut lo, mut hi) = (0, p.points().len());
    while lo < hi {
        w.charge(1).map_err(Error::from)?;
        let mid = lo + (hi - lo) / 2;
        match p.points()[mid].entity().cmp(&id) {
            std::cmp::Ordering::Equal => return Ok(mid),
            std::cmp::Ordering::Less => lo = mid + 1,
            std::cmp::Ordering::Greater => hi = mid,
        }
    }
    Err(Error::Invalid)
}
fn cached(p: &Parsed<'_>, cache: &[Option<[I; 3]>], id: u32, w: &mut Work<'_>) -> Result<[I; 3]> {
    let i = point_index(p, id, w)?;
    cache[i].ok_or(Error::Invalid)
}
struct ControlCoordinates<'a, 'i> {
    p: &'a Parsed<'i>,
    cache: &'a [Option<[I; 3]>],
    heights: [Scalar<'i>; 2],
    checked: &'a mut [bool],
    bounds: &'a mut [[f64; 2]; 3],
}
#[derive(Clone, Copy)]
struct LiteralNet {
    data: [[P; 3]; 2],
    points: [[usize; 12]; 2],
}
fn controls(
    coordinates: &mut ControlCoordinates<'_, '_>,
    t: &TensorSurface<'_>,
    span: usize,
    memory: &mut Memory,
    w: &mut Work<'_>,
) -> Result<LiteralNet> {
    let mut rows = [[P::controls(&[I::point(0.)]); 3]; 2];
    let mut indices = [[0; 12]; 2];
    for (r, ids) in t.rows().iter().enumerate() {
        let mut points = [[I::point(0.); 3]; 12];
        for (i, point) in points.iter_mut().enumerate() {
            let index = point_index(coordinates.p, ids[span * 11 + i], w)?;
            indices[r][i] = index;
            if !coordinates.checked[index] {
                height_bound(
                    &coordinates.p.points()[index],
                    coordinates.heights,
                    memory,
                    w,
                )?;
                coordinates.checked[index] = true;
            }
            *point = coordinates.cache[index].ok_or(Error::Invalid)?;
            for (axis, value) in point.iter().enumerate() {
                coordinates.bounds[axis][0] = coordinates.bounds[axis][0].min(value.lo);
                coordinates.bounds[axis][1] = coordinates.bounds[axis][1].max(value.hi);
            }
        }
        for a in 0..3 {
            rows[r][a] = P::controls(&std::array::from_fn::<_, 12, _>(|i| points[i][a]));
        }
    }
    Ok(LiteralNet {
        data: rows,
        points: indices,
    })
}
fn strip_jacobian(
    p: &Parsed<'_>,
    net: &LiteralNet,
    factor: bool,
    w: &mut Work<'_>,
) -> Result<Option<[P; 3]>> {
    for axis in 0..2 {
        for i in 0..12 {
            w.charge(1).map_err(Error::from)?;
            let (a, b) = (net.data[0][axis].c[i], net.data[1][axis].c[i]);
            if a.lo != b.lo || a.hi != b.hi {
                return Ok(None);
            }
            let a = p.points()[net.points[0][i]].coordinates()[axis].spelling();
            let b = p.points()[net.points[1][i]].coordinates()[axis].spelling();
            w.charge(a.len() + b.len() + 1).map_err(Error::from)?;
            if a != b {
                return Ok(None);
            }
        }
    }
    // Actual identicalXY coefficients give Dxy=0 on the WHOLE tensor.
    // Consequently Jac=-Dz*(Ax Ay'-Ay Ax'), independent ofu andZ'.
    let a = net.data[0];
    let k = a[0]
        .product(a[1].derivative(w).map_err(Error::from)?, w)
        .map_err(Error::from)?
        .sub(
            a[1].product(a[0].derivative(w).map_err(Error::from)?, w)
                .map_err(Error::from)?,
            w,
        )
        .map_err(Error::from)?;
    let mut dz = net.data[1][2].sub(a[2], w).map_err(Error::from)?;
    if factor {
        dz.c[0] = I::point(0.);
        dz = dz.remove_t(w).map_err(Error::from)?;
    }
    let j = k
        .product(dz, w)
        .map_err(Error::from)?
        .scale(-1., w)
        .map_err(Error::from)?;
    Ok(Some([j; 3]))
}
fn height_bound(
    point: &PointLeaf<'_>,
    heights: [Scalar<'_>; 2],
    memory: &mut Memory,
    w: &mut Work<'_>,
) -> Result<()> {
    let z = point.coordinates()[2];
    let equal = |h: Scalar<'_>| {
        z.spelling() == h.spelling() || (z.is_zero() && h.is_zero()) || (z.is_one() && h.is_one())
    };
    // Strict cache order implies strict order of the encoded reals. Equality
    // requires original-leaf equality, never equality of rounded caches alone.
    for (i, height) in heights.into_iter().enumerate() {
        w.charge(1).map_err(Error::from)?;
        if if i == 0 {
            z.value() > height.value()
        } else {
            z.value() < height.value()
        } {
            continue;
        }
        // Only an unresolved cache comparison scans original numeric leaves.
        w.charge(z.spelling().len() + height.spelling().len() + 1)
            .map_err(Error::from)?;
        if equal(height) {
            continue;
        }
        let order = crate::step::insert::exact_cap::compare_decimal(
            z.spelling(),
            height.spelling(),
            memory,
            w,
        )
        .map_err(Error::from)?;
        if (i == 0 && order.is_lt()) || (i == 1 && order.is_gt()) {
            return Err(Error::Unsupported);
        }
    }
    Ok(())
}

fn root(parents: &mut [usize], mut i: usize, w: &mut Work<'_>) -> Result<usize> {
    let mut cursor = i;
    while parents[i] != i {
        w.charge(1).map_err(Error::from)?;
        i = parents[i];
    }
    while parents[cursor] != i {
        w.charge(1).map_err(Error::from)?;
        let next = parents[cursor];
        parents[cursor] = i;
        cursor = next;
    }
    Ok(i)
}

fn components(
    b: &SupportBindings<'_, '_>,
    memory: &mut Memory,
    w: &mut Work<'_>,
) -> Result<Vec<u8>> {
    let p = b.source();
    partition_sides(
        p.faces().len(),
        &|i| matches!(p.faces()[i].surface(), SurfaceIndex::Tensor(_)),
        b.edge_face_indices(),
        memory,
        w,
    )
}
// This graph check is independent of the geometric charts: exactly two side
// components are required before radial separation can assign their ownership.
fn partition_sides(
    count: usize,
    is_side: &dyn Fn(usize) -> bool,
    owners: &[[usize; 2]],
    memory: &mut Memory,
    w: &mut Work<'_>,
) -> Result<Vec<u8>> {
    let mut parents = memory.allocate(count, w).map_err(Error::from)?;
    parents.extend(0..count);
    for &[a, c] in owners {
        w.charge(1).map_err(Error::from)?;
        if is_side(a) && is_side(c) {
            let a = root(&mut parents, a, w)?;
            let c = root(&mut parents, c, w)?;
            parents[a] = c;
        }
    }
    let mut roots = [None, None];
    let mut out = memory.allocate(count, w).map_err(Error::from)?;
    for i in 0..count {
        if !is_side(i) {
            out.push(2);
            continue;
        }
        let r = root(&mut parents, i, w)?;
        let component = if roots[0] == Some(r) {
            0
        } else if roots[1] == Some(r) {
            1
        } else if roots[0].is_none() {
            roots[0] = Some(r);
            0
        } else if roots[1].is_none() {
            roots[1] = Some(r);
            1
        } else {
            return Err(Error::Unsupported);
        };
        out.push(component);
    }
    if roots[1].is_none() {
        return Err(Error::Unsupported);
    }
    memory.release(&parents);
    drop(parents);
    Ok(out)
}
fn manifold(b: &SupportBindings<'_, '_>, memory: &mut Memory, w: &mut Work<'_>) -> Result<()> {
    let p = b.source();
    let mut edges = memory.allocate(p.edges().len(), w).map_err(Error::from)?;
    for edge in p.edges() {
        w.charge(64).map_err(Error::from)?;
        let [a, c] = edge.vertex_entities();
        let a = p
            .vertices()
            .binary_search_by_key(&a, |v| v.entity)
            .map_err(|_| Error::Invalid)?;
        let c = p
            .vertices()
            .binary_search_by_key(&c, |v| v.entity)
            .map_err(|_| Error::Invalid)?;
        edges.push(ShellEdge {
            ends: [a, c],
            geometry: EdgeGeometry::Ruling,
        });
    }
    let mut faces = memory.allocate(p.faces().len(), w).map_err(Error::from)?;
    for (i, face) in p.faces().iter().enumerate() {
        let mut loops = memory
            .allocate(face.loops().len(), w)
            .map_err(Error::from)?;
        for wire in face.loops() {
            let mut uses = memory
                .allocate(wire.coedges().len(), w)
                .map_err(Error::from)?;
            for c in wire.coedges() {
                w.charge(32).map_err(Error::from)?;
                let edge = p
                    .edges()
                    .binary_search_by_key(&c.edge, Edge::entity)
                    .map_err(|_| Error::Invalid)?;
                uses.push(ShellCoedge {
                    edge,
                    forward: c.forward,
                });
            }
            loops.push(uses);
        }
        faces.push(ShellFace {
            geometry: match face.surface() {
                SurfaceIndex::Tensor(t) => FaceGeometry::Tensor(t),
                SurfaceIndex::Plane(_) => FaceGeometry::Cap(u8::from(i == b.cap_face_indices()[1])),
            },
            same_sense: face.same_sense(),
            loops,
        });
    }
    topology::verify(p.vertices().len(), &edges, &faces, memory, w).map_err(Error::from)?;
    for face in &faces {
        for wire in &face.loops {
            memory.release(wire);
        }
        memory.release(&face.loops);
    }
    memory.release(&faces);
    memory.release(&edges);
    drop(faces);
    drop(edges);
    Ok(())
}
fn caps(
    b: &SupportBindings<'_, '_>,
    components: &[u8],
    outer: u8,
    cache: &[Option<[I; 3]>],
    stack: &mut Vec<projection::Cell>,
    memory: &mut Memory,
    w: &mut Work<'_>,
) -> Result<f64> {
    let p = b.source();
    let mut minimum = f64::INFINITY;
    for &face_index in &b.cap_face_indices() {
        let face = &p.faces()[face_index];
        for wire in face.loops() {
            let desired =
                if face.same_sense() { 1 } else { -1 } * if wire.is_outer() { 1 } else { -1 };
            let mut winding = 0;
            for co in wire.coedges() {
                w.charge(64).map_err(Error::from)?;
                let edge_index = p
                    .edges()
                    .binary_search_by_key(&co.edge, Edge::entity)
                    .map_err(|_| Error::Invalid)?;
                let edge = &p.edges()[edge_index];
                let other = b.edge_face_indices()[edge_index]
                    .iter()
                    .copied()
                    .find(|&i| i != face_index)
                    .ok_or(Error::Invalid)?;
                let expected = if wire.is_outer() { outer } else { 1 - outer };
                if components[other] != expected {
                    return Err(Error::Unsupported);
                }
                let curve = p.curve(edge.curve_entity()).ok_or(Error::Invalid)?;
                let spans = curve.knots().len() - 1;
                for j in 0..spans {
                    let span = if co.forward { j } else { spans - 1 - j };
                    let mut points = [[I::point(0.); 3]; 12];
                    for (i, point) in points.iter_mut().enumerate() {
                        *point = cached(p, cache, curve.pole_entities()[span * 11 + i], w)?;
                    }
                    halfplane(&points, w)?;
                    let [a, c] = if co.forward {
                        [points[0], points[11]]
                    } else {
                        [points[11], points[0]]
                    };
                    winding += crossing(a, c, w)?;
                    let x = P::controls(&points.map(|v| v[0]));
                    let y = P::controls(&points.map(|v| v[1]));
                    let sign = f64::from(desired) * if co.forward { 1. } else { -1. };
                    let h = [x.hull(), y.hull()];
                    let scale = h[0].magnitude().max(h[1].magnitude()).max(1.);
                    w.charge(24).map_err(Error::from)?;
                    let tiny = h.iter().all(|a| a.hi - a.lo <= 64. * f64::EPSILON * scale);
                    let proof = if tiny {
                        None
                    } else {
                        let k = x
                            .product(y.derivative(w).map_err(Error::from)?, w)
                            .map_err(Error::from)?
                            .sub(
                                y.product(x.derivative(w).map_err(Error::from)?, w)
                                    .map_err(Error::from)?,
                                w,
                            )
                            .map_err(Error::from)?
                            .scale(sign, w)
                            .map_err(Error::from)?;
                        Some(projection::positive([k; 3], stack, w))
                    };
                    let lower = if matches!(
                        proof,
                        None | Some(Err(crate::step::insert::Error::Unresolved))
                    ) {
                        stack.clear();
                        let mut tokens = [[""; 2]; 12];
                        for (i, xy) in tokens.iter_mut().enumerate() {
                            let index = point_index(p, curve.pole_entities()[span * 11 + i], w)?;
                            for (axis, token) in xy.iter_mut().enumerate() {
                                *token = p.points()[index].coordinates()[axis].spelling();
                            }
                        }
                        crate::step::insert::exact_cap::angular_lower(tokens, sign < 0., memory, w)
                            .map_err(Error::from)?
                    } else {
                        proof
                            .expect("ordinary cap proof selected")
                            .map_err(Error::from)?
                            .0
                            .lo
                    };
                    minimum = minimum.min(lower);
                }
            }
            if winding != desired {
                return Err(Error::Unsupported);
            }
        }
    }
    Ok(minimum)
}
fn halfplane(points: &[[I; 3]; 12], w: &mut Work<'_>) -> Result<()> {
    for [x, y] in [
        [1., 0.],
        [1., 1.],
        [0., 1.],
        [-1., 1.],
        [-1., 0.],
        [-1., -1.],
        [0., -1.],
        [1., -1.],
    ] {
        let mut fits = true;
        for p in points {
            w.charge(1).map_err(Error::from)?;
            let h = p[0] * I::point(x) + p[1] * I::point(y);
            if !h.valid() || h.lo <= 0. {
                fits = false;
                break;
            }
        }
        if fits {
            return Ok(());
        }
    }
    Err(Error::Unsupported)
}
fn crossing(a: [I; 3], b: [I; 3], w: &mut Work<'_>) -> Result<i32> {
    w.charge(1).map_err(Error::from)?;
    let up = a[1].hi <= 0. && b[1].lo > 0.;
    let down = b[1].hi <= 0. && a[1].lo > 0.;
    if !up && !down {
        if (a[1].lo > 0. && b[1].lo > 0.) || (a[1].hi <= 0. && b[1].hi <= 0.) {
            return Ok(0);
        }
        return Err(Error::Unsupported);
    }
    let cross = a[0] * b[1] - a[1] * b[0];
    if !cross.valid() || (cross.lo <= 0. && cross.hi >= 0.) {
        return Err(Error::Unsupported);
    }
    Ok(if up && cross.lo > 0. {
        1
    } else if down && cross.hi < 0. {
        -1
    } else {
        0
    })
}

#[cfg(test)]
mod predicate_tests {
    use super::*;
    #[test]
    fn metered_point_lookup_matches_actual_unique_arena_and_charges_real_comparisons() {
        let text = crate::step_read::owned_insert::tests::file(false);
        let parsed = parse(&text, Options::default(), None).unwrap();
        let mut total = 0;
        for (expected, point) in parsed.points().iter().enumerate() {
            let mut w = Work::new(16_000_000, None).unwrap();
            assert_eq!(
                point_index(&parsed, point.entity(), &mut w).unwrap(),
                expected
            );
            assert_eq!(
                parsed
                    .points()
                    .binary_search_by_key(&point.entity(), PointLeaf::entity),
                Ok(expected)
            );
            total += w.used;
            if w.used > 1 {
                let mut short = Work::new(w.used - 1, None).unwrap();
                assert!(matches!(
                    point_index(&parsed, point.entity(), &mut short),
                    Err(Error::Limit)
                ));
            }
        }
        assert!(total < parsed.points().len() * 32);
        for absent in [0, u32::MAX] {
            let mut w = Work::new(16_000_000, None).unwrap();
            assert!(matches!(
                point_index(&parsed, absent, &mut w),
                Err(Error::Invalid)
            ));
        }
        let middle = parsed.points()[parsed.points().len() / 2].entity();
        let mut one = Work::new(1, None).unwrap();
        assert_eq!(
            point_index(&parsed, middle, &mut one).unwrap(),
            parsed.points().len() / 2
        );
        let mut cancelled = Work::new(16_000_000, Some(&|| true)).unwrap();
        assert!(matches!(
            point_index(&parsed, middle, &mut cancelled),
            Err(Error::Cancelled)
        ));
    }
    #[test]
    fn actual_side_graph_partition_requires_two_connected_components() {
        let mut w = Work::new(16_000_000, None).unwrap();
        let mut memory = Memory::new(0, 64 * 1024 * 1024).unwrap();
        assert!(matches!(
            partition_sides(3, &|_| true, &[[0, 1], [1, 2]], &mut memory, &mut w),
            Err(Error::Unsupported)
        ));
        assert!(matches!(
            partition_sides(3, &|_| true, &[], &mut memory, &mut w),
            Err(Error::Unsupported)
        ));
        let ids = partition_sides(
            5,
            &|i| i != 4,
            &[[0, 1], [2, 3], [0, 4], [2, 4]],
            &mut memory,
            &mut w,
        )
        .unwrap();
        assert_eq!(ids, vec![0, 0, 1, 1, 2]);
    }
    #[test]
    fn finite_convex_control_hulls_refuse_overflowing_radial_norms() {
        let mut w = Work::new(16_000_000, None).unwrap();
        let row = [
            P::controls(&[I::point(1e308)]),
            P::controls(&[I::point(1e308)]),
            P::controls(&[I::point(0.)]),
        ];
        assert!(matches!(
            convex_radius([row; 2], &mut w),
            Err(Error::Unsupported)
        ));
        let ordinary = [
            P::controls(&[I::point(3.)]),
            P::controls(&[I::point(4.)]),
            P::controls(&[I::point(0.)]),
        ];
        let bounds = convex_radius([ordinary; 2], &mut w).unwrap();
        assert!(bounds.lo <= 5. && bounds.hi >= 5.);
    }
    #[test]
    fn ray_and_homotopy_proofs_refuse_ambiguous_axis_crossings_and_origin_hulls() {
        let mut w = Work::new(16_000_000, None).unwrap();
        let p = |x, y| [I::point(x), I::point(y), I::point(0.)];
        assert!(matches!(
            halfplane(&[p(0., 0.); 12], &mut w),
            Err(Error::Unsupported)
        ));
        assert!(matches!(
            crossing(
                [I::point(1.), I::new(-1e-16, 1e-16), I::point(0.)],
                p(1., 1.),
                &mut w
            ),
            Err(Error::Unsupported)
        ));
        assert!(matches!(
            crossing(p(-1., -1.), p(1., 1.), &mut w),
            Err(Error::Unsupported)
        ));
        assert_eq!(crossing(p(1., -1.), p(1., 1.), &mut w).unwrap(), 1);
        assert_eq!(crossing(p(1., 1.), p(1., -1.), &mut w).unwrap(), -1);
        assert_eq!(crossing(p(-1., -1.), p(-1., 1.), &mut w).unwrap(), 0);
    }
}
