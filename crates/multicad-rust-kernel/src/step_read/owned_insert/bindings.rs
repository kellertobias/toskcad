//! Exact literal coefficient incidence and full parameter-domain boundaries.
//! This receipt does not certify global embedding or material volume.
use super::*;
#[derive(Clone, Copy, Debug)]
pub struct TensorBinding {
    face: usize,
    tensor: usize,
    corners: [u32; 4],
    collapsed: [bool; 2],
}
impl TensorBinding {
    pub fn face_index(&self) -> usize {
        self.face
    }
    pub fn tensor_index(&self) -> usize {
        self.tensor
    }
    pub fn corner_vertex_entities(&self) -> [u32; 4] {
        self.corners
    }
    pub fn collapsed_v_edges(&self) -> [bool; 2] {
        self.collapsed
    }
}
pub struct SupportBindings<'a, 'i> {
    source: &'a Parsed<'i>,
    tensors: Vec<TensorBinding>,
    caps: [usize; 2],
    heights: [Scalar<'i>; 2],
    owners: Vec<[usize; 2]>,
    work: usize,
    bytes: usize,
}
impl<'a, 'i> SupportBindings<'a, 'i> {
    pub fn source(&self) -> &'a Parsed<'i> {
        self.source
    }
    pub fn tensors(&self) -> &[TensorBinding] {
        &self.tensors
    }
    pub fn cap_face_indices(&self) -> [usize; 2] {
        self.caps
    }
    pub fn cap_heights(&self) -> [Scalar<'i>; 2] {
        self.heights
    }
    pub fn edge_face_indices(&self) -> &[[usize; 2]] {
        &self.owners
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.bytes
    }
}
#[derive(Clone, Copy)]
struct Use {
    face: usize,
    forward: bool,
}
impl<'i> Parsed<'i> {
    pub fn support_bindings<'a>(
        &'a self,
        options: Options,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<SupportBindings<'a, 'i>> {
        if options.maximum_work == 0
            || options.maximum_work > 16_000_000
            || options.retained_budget_bytes == 0
            || options.retained_budget_bytes > 64 * 1024 * 1024
        {
            return Err(Error::Invalid);
        }
        let mut w = Ledger {
            used: 0,
            remaining: options.maximum_work,
            interrupt,
        };
        w.step(self.work)?;
        let mut uses = vec![[None::<Use>; 2]; self.edges.len()];
        let mut caps = [None, None];
        let mut tensors = Vec::new();
        let heights = self
            .planes
            .iter()
            .map(|plane| plane_height(self, plane, &mut w))
            .collect::<Result<Vec<_>>>()?;
        // Different correctly-rounded caches prove strict order of the actual
        // decimal heights. Equal caches are ambiguous; never invent a datum.
        let start = if heights[0].value() < heights[1].value() {
            0
        } else if heights[1].value() < heights[0].value() {
            1
        } else {
            return Err(Error::Unsupported);
        };
        let end = 1 - start;
        let mut support_ids = self
            .faces
            .iter()
            .map(|f| f.surface_entity)
            .collect::<Vec<_>>();
        w.step(support_ids.len() * 128)?;
        support_ids.sort_unstable();
        if support_ids.windows(2).any(|p| p[0] == p[1]) {
            return Err(Error::Unsupported);
        }
        for (face_index, face) in self.faces.iter().enumerate() {
            w.step(128)?;
            let mut borders = [false; 4];
            let mut corners = [None; 4];
            for wire in &face.loops {
                for (i, coedge) in wire.coedges.iter().enumerate() {
                    w.step(256)?;
                    let edge_index = self
                        .edges
                        .binary_search_by_key(&coedge.edge, |e| e.entity)
                        .map_err(|_| Error::Invalid)?;
                    let edge = &self.edges[edge_index];
                    let slots = &mut uses[edge_index];
                    let at = usize::from(slots[0].is_some());
                    if slots[at].is_some() {
                        return Err(Error::Invalid);
                    }
                    slots[at] = Some(Use {
                        face: face_index,
                        forward: coedge.forward,
                    });
                    let previous = wire.coedges[(i + wire.coedges.len() - 1) % wire.coedges.len()];
                    let previous = &self.edges[self
                        .edges
                        .binary_search_by_key(&previous.edge, |e| e.entity)
                        .map_err(|_| Error::Invalid)?];
                    let previous_forward =
                        wire.coedges[(i + wire.coedges.len() - 1) % wire.coedges.len()].forward;
                    if previous.ends[usize::from(previous_forward)]
                        != edge.ends[usize::from(!coedge.forward)]
                    {
                        return Err(Error::Invalid);
                    }
                    let [a, b] = edge.pcurves;
                    let pc = match (
                        a.support_entity == face.surface_entity,
                        b.support_entity == face.surface_entity,
                    ) {
                        (true, false) => a,
                        (false, true) => b,
                        _ => return Err(Error::Unsupported),
                    };
                    let curve = self.curve(edge.curve_entity).ok_or(Error::Invalid)?;
                    let uv = self.curve(pc.curve_entity).ok_or(Error::Invalid)?;
                    match face.surface {
                        SurfaceIndex::Plane(plane) => {
                            cap_curve(self, curve, uv, heights[plane], &mut w)?
                        }
                        SurfaceIndex::Tensor(tensor) => {
                            let tensor = &self.tensors[tensor];
                            let coordinates = tensor_curve(self, tensor, curve, uv, &mut w)?;
                            let oriented = if coedge.forward {
                                coordinates
                            } else {
                                [coordinates[1], coordinates[0]]
                            };
                            let border = border(oriented, face.same_sense)?;
                            if borders[border] {
                                return Err(Error::Invalid);
                            }
                            borders[border] = true;
                            for (i, uv) in coordinates.into_iter().enumerate() {
                                let corner = corner(uv);
                                if corners[corner].is_some_and(|old| old != edge.ends[i]) {
                                    return Err(Error::Invalid);
                                }
                                corners[corner] = Some(edge.ends[i]);
                            }
                        }
                    }
                }
            }
            match face.surface {
                SurfaceIndex::Plane(plane) => {
                    let end_index = usize::from(plane != start);
                    if face.same_sense != (end_index == 1) || caps[end_index].is_some() {
                        return Err(Error::Invalid);
                    }
                    caps[end_index] = Some(face_index);
                }
                SurfaceIndex::Tensor(tensor) => {
                    let net = &self.tensors[tensor];
                    let last = net.rows[0].len() - 1;
                    let collapsed = [
                        net.rows[0][0] == net.rows[1][0],
                        net.rows[0][last] == net.rows[1][last],
                    ];
                    if collapsed == [true, true]
                        || borders != [!collapsed[0], true, !collapsed[1], true]
                    {
                        return Err(Error::Unsupported);
                    }
                    let [Some(a), Some(b), Some(c), Some(d)] = corners else {
                        return Err(Error::Invalid);
                    };
                    if (collapsed[0] && a != b) || (collapsed[1] && c != d) {
                        return Err(Error::Invalid);
                    }
                    tensors.push(TensorBinding {
                        face: face_index,
                        tensor,
                        corners: [a, b, c, d],
                        collapsed,
                    });
                }
            }
        }
        let [Some(a), Some(b)] = caps else {
            return Err(Error::Unsupported);
        };
        if tensors.len() != self.tensors.len() {
            return Err(Error::Unsupported);
        }
        let mut owners = Vec::with_capacity(uses.len());
        for u in &uses {
            w.step(1)?;
            let [Some(a), Some(b)] = u else {
                return Err(Error::Invalid);
            };
            if a.face == b.face || a.forward == b.forward {
                return Err(Error::Invalid);
            }
            owners.push([a.face, b.face]);
        }
        let bytes = self
            .bytes
            .saturating_add(size_of::<SupportBindings<'_, '_>>())
            .saturating_add(tensors.capacity() * size_of::<TensorBinding>())
            .saturating_add(owners.capacity() * size_of::<[usize; 2]>());
        if bytes > options.retained_budget_bytes {
            return Err(Error::Limit);
        }
        w.step(1)?;
        Ok(SupportBindings {
            source: self,
            tensors,
            caps: [a, b],
            heights: [heights[start], heights[end]],
            owners,
            work: w.used,
            bytes,
        })
    }
}
fn same(a: Scalar<'_>, b: Scalar<'_>, w: &mut Ledger<'_>) -> Result<bool> {
    w.step(a.spelling().len() + b.spelling().len() + 1)?;
    Ok(a.spelling() == b.spelling() || (a.is_zero() && b.is_zero()) || (a.is_one() && b.is_one()))
}
fn same_knots(
    a: &SplineCurve<'_>,
    mult: &[u32],
    knots: &[Scalar<'_>],
    w: &mut Ledger<'_>,
) -> Result<bool> {
    w.step(a.multiplicities.len() + 1)?;
    if a.multiplicities != mult || a.knots.len() != knots.len() {
        return Ok(false);
    }
    for (&a, &b) in a.knots.iter().zip(knots) {
        if !same(a, b, w)? {
            return Ok(false);
        }
    }
    Ok(true)
}
fn coords<'a, 'i>(p: &'a Parsed<'i>, id: u32, n: usize) -> Result<&'a [Scalar<'i>]> {
    let value = p.point(id).ok_or(Error::Invalid)?.coordinates();
    if value.len() != n {
        return Err(Error::Invalid);
    }
    Ok(value)
}
fn cap_curve(
    p: &Parsed<'_>,
    curve: &SplineCurve<'_>,
    uv: &SplineCurve<'_>,
    height: Scalar<'_>,
    w: &mut Ledger<'_>,
) -> Result<()> {
    if curve.degree != 11
        || uv.degree != 11
        || !same_knots(curve, &uv.multiplicities, &uv.knots, w)?
        || curve.poles.len() != uv.poles.len()
    {
        return Err(Error::Unsupported);
    }
    for (&a, &b) in curve.poles.iter().zip(&uv.poles) {
        w.step(64)?;
        let xyz = coords(p, a, 3)?;
        let xy = coords(p, b, 2)?;
        if !same(xyz[2], height, w)? || !same(xyz[0], xy[0], w)? || !same(xyz[1], xy[1], w)? {
            return Err(Error::Unsupported);
        }
    }
    Ok(())
}
fn bit(s: Scalar<'_>, w: &mut Ledger<'_>) -> Result<u8> {
    w.step(s.spelling().len() + 1)?;
    if s.is_zero() {
        Ok(0)
    } else if s.is_one() {
        Ok(1)
    } else {
        Err(Error::Unsupported)
    }
}
fn tensor_curve(
    p: &Parsed<'_>,
    tensor: &TensorSurface<'_>,
    curve: &SplineCurve<'_>,
    uv: &SplineCurve<'_>,
    w: &mut Ledger<'_>,
) -> Result<[[u8; 2]; 2]> {
    if uv.degree != 1 || uv.poles.len() != 2 {
        return Err(Error::Unsupported);
    }
    let a = coords(p, uv.poles[0], 2)?;
    let b = coords(p, uv.poles[1], 2)?;
    let values = [
        [bit(a[0], w)?, bit(a[1], w)?],
        [bit(b[0], w)?, bit(b[1], w)?],
    ];
    if curve.degree == 11 {
        w.step(curve.poles.len() * 4)?;
        let row = tensor
            .rows
            .iter()
            .position(|row| row == &curve.poles)
            .ok_or(Error::Unsupported)?;
        if !same_knots(curve, &tensor.multiplicities_v, &tensor.knots_v, w)?
            || values != [[row as u8, 0], [row as u8, 1]]
        {
            return Err(Error::Unsupported);
        }
    } else {
        if curve.poles.len() != 2 || values[0][1] != values[1][1] || values[0][0] == values[1][0] {
            return Err(Error::Unsupported);
        }
        let last = tensor.rows[0].len() - 1;
        let expected = values.map(|[u, v]| tensor.rows[u as usize][if v == 0 { 0 } else { last }]);
        if curve.poles.as_slice() != expected {
            return Err(Error::Unsupported);
        }
    }
    Ok(values)
}
fn corner([u, v]: [u8; 2]) -> usize {
    match (u, v) {
        (0, 0) => 0,
        (1, 0) => 1,
        (1, 1) => 2,
        _ => 3,
    }
}
fn border(mut endpoints: [[u8; 2]; 2], same_sense: bool) -> Result<usize> {
    if !same_sense {
        endpoints.reverse();
    }
    [
        [[0, 0], [1, 0]],
        [[1, 0], [1, 1]],
        [[1, 1], [0, 1]],
        [[0, 1], [0, 0]],
    ]
    .iter()
    .position(|&edge| edge == endpoints)
    .ok_or(Error::Unsupported)
}
fn direction<'i>(p: &Parsed<'i>, id: u32, w: &mut Ledger<'_>) -> Result<[Scalar<'i>; 3]> {
    let e = kind(&p.entities, id, "DIRECTION")?;
    let [Value::String(_), values] = e.args.as_slice() else {
        return Err(Error::Invalid);
    };
    let [a, b, c] = list(values)? else {
        return Err(Error::Unsupported);
    };
    let mut n = Numeric::new(e, w)?;
    let out = [n.scalar(a)?, n.scalar(b)?, n.scalar(c)?];
    n.finish()?;
    Ok(out)
}
fn plane_height<'i>(
    p: &Parsed<'i>,
    plane: &PlaneSurface,
    w: &mut Ledger<'_>,
) -> Result<Scalar<'i>> {
    let origin = coords(p, plane.origin, 3)?;
    let normal = direction(p, plane.normal, w)?;
    let reference = direction(p, plane.reference, w)?;
    // Positive axial ratios normalize to EXACTLY these axes. No proximity,
    // caller tolerance, fitting or inferred authoring frame is used.
    if !origin[0].is_zero()
        || !origin[1].is_zero()
        || !normal[0].is_zero()
        || !normal[1].is_zero()
        || normal[2].value() <= 0.
        || reference[0].value() <= 0.
        || !reference[1].is_zero()
        || !reference[2].is_zero()
    {
        return Err(Error::Unsupported);
    }
    Ok(origin[2])
}
