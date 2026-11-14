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
                    w.charge(c[i].spelling().len() + 1).map_err(Error::from)?;
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
            for row in t.rows() {
                for &id in row {
                    // One real POINT lookup serves bounds and memoized height
                    // validation. The old loop performed two separate lookups.
                    w.charge(32).map_err(Error::from)?;
                    let index = p
                        .points()
                        .binary_search_by_key(&id, PointLeaf::entity)
                        .map_err(|_| Error::Invalid)?;
                    if !height_checked[index] {
                        height_bound(&p.points()[index], heights, &mut w)?;
                        height_checked[index] = true;
                    }
                    let point = cache[index].ok_or(Error::Invalid)?;
                    for a in 0..3 {
                        bounds[a][0] = bounds[a][0].min(point[a].lo);
                        bounds[a][1] = bounds[a][1].max(point[a].hi);
                    }
                }
            }
            for span in 0..t.knots_v().len() - 1 {
                let data = controls(p, &cache, t, span, false, &mut w)?;
                nets.push(data);
                let r = projection::radius_squared(data, &mut w).map_err(Error::from)?;
                let (r, _) = projection::positive(r, &mut stack, &mut w).map_err(Error::from)?;
                let r = r.sqrt();
                radii[component][0] = radii[component][0].min(r.lo);
                radii[component][1] = radii[component][1].max(r.hi);
            }
        }
        memory.release(&height_checked);
        drop(height_checked);
        let outer = if radii[0][0] > radii[1][1] {
            0
        } else if radii[1][0] > radii[0][1] {
            1
        } else {
            return Err(Error::Unsupported);
        };
        let mut minimum = f64::INFINITY;
        let mut net_index = 0;
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
                let mut data = nets[net_index];
                net_index += 1;
                if reverse {
                    w.charge(6 * 12).map_err(Error::from)?;
                    for row in &mut data {
                        for axis in row {
                            axis.c[..=axis.degree].reverse();
                        }
                    }
                }
                let j = projection::jacobian(data, factor, &mut w).map_err(Error::from)?;
                let sign = side * same * if reverse { -1. } else { 1. };
                let j = [
                    j[0].scale(sign, &mut w).map_err(Error::from)?,
                    j[1].scale(sign, &mut w).map_err(Error::from)?,
                    j[2].scale(sign, &mut w).map_err(Error::from)?,
                ];
                let (j, _) = projection::positive(j, &mut stack, &mut w).map_err(Error::from)?;
                minimum = minimum.min(j.lo);
            }
        }
        memory.release(&nets);
        drop(nets);
        let cap_min = caps(self, &components, outer, p, &cache, &mut stack, &mut w)?;
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
fn cached(p: &Parsed<'_>, cache: &[Option<[I; 3]>], id: u32) -> Result<[I; 3]> {
    let i = p
        .points()
        .binary_search_by_key(&id, PointLeaf::entity)
        .map_err(|_| Error::Invalid)?;
    cache[i].ok_or(Error::Invalid)
}
fn controls(
    p: &Parsed<'_>,
    cache: &[Option<[I; 3]>],
    t: &TensorSurface<'_>,
    span: usize,
    reverse: bool,
    w: &mut Work<'_>,
) -> Result<[[P; 3]; 2]> {
    let mut rows = [[P::controls(&[I::point(0.)]); 3]; 2];
    for (r, ids) in t.rows().iter().enumerate() {
        let mut points = [[I::point(0.); 3]; 12];
        for (i, point) in points.iter_mut().enumerate() {
            w.charge(32).map_err(Error::from)?;
            *point = cached(p, cache, ids[span * 11 + if reverse { 11 - i } else { i }])?;
        }
        for a in 0..3 {
            rows[r][a] = P::controls(&std::array::from_fn::<_, 12, _>(|i| points[i][a]));
        }
    }
    Ok(rows)
}
fn height_bound(point: &PointLeaf<'_>, heights: [Scalar<'_>; 2], w: &mut Work<'_>) -> Result<()> {
    let z = point.coordinates()[2];
    w.charge(z.spelling().len() + heights[0].spelling().len() + heights[1].spelling().len() + 1)
        .map_err(Error::from)?;
    let equal = |h: Scalar<'_>| {
        z.spelling() == h.spelling() || (z.is_zero() && h.is_zero()) || (z.is_one() && h.is_one())
    };
    // Strict cache order implies strict order of the encoded reals. Equality
    // requires original-leaf equality, never equality of rounded caches alone.
    if !(equal(heights[0]) || z.value() > heights[0].value())
        || !(equal(heights[1]) || z.value() < heights[1].value())
    {
        return Err(Error::Unsupported);
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
    p: &Parsed<'_>,
    cache: &[Option<[I; 3]>],
    stack: &mut Vec<projection::Cell>,
    w: &mut Work<'_>,
) -> Result<f64> {
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
                        w.charge(32).map_err(Error::from)?;
                        *point = cached(p, cache, curve.pole_entities()[span * 11 + i])?;
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
                    let (bound, _) = projection::positive([k; 3], stack, w).map_err(Error::from)?;
                    minimum = minimum.min(bound.lo);
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
