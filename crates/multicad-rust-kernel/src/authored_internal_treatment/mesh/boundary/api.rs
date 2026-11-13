//! Immutable finite source graph, with construction caches discarded on success.
use super::*;
use crate::authored_flat_thread::PointSample;
pub struct Boundary {
    pub(super) graph: Graph,
    work_used: usize,
    vertices: usize,
    edges: usize,
}
pub use super::loops::Coedge as BoundaryCoedge;
impl BoundaryOptions {
    pub(super) fn validate(self) -> Result<(), Error> {
        if self.cells == 0
            || self.cells > 2048
            || self.vertices == 0
            || self.vertices > 8192
            || self.curves == 0
            || self.curves > 8192
            || self.work == 0
            || self.work > 16_000_000
            || self.bytes == 0
            || self.bytes > BYTES
        {
            Err(Error::Invalid)
        } else {
            Ok(())
        }
    }
}
impl Cavity {
    pub fn boundary_with_interrupt(
        &self,
        options: BoundaryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Boundary, Error> {
        build(*self, None, options, interrupt)
    }
}
impl Material {
    /// Geometry and coordinates are stock-relative, as in stock_chart_origin().
    pub fn boundary_with_interrupt(
        &self,
        options: BoundaryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Boundary, Error> {
        build(*self.cavity(), Some(*self), options, interrupt)
    }
}
fn build(
    source: Cavity,
    material: Option<Material>,
    options: BoundaryOptions,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Boundary, Error> {
    if interrupt.is_some_and(|f| f()) {
        return Err(Error::Cancelled);
    }
    options.validate()?;
    let mut work = Work {
        remaining: options.work,
        used: 0,
        interrupt,
    };
    let mut graph = Graph::raw_in(source, options, &mut work)?.closed(&mut work)?;
    if let Some(material) = material {
        graph = graph.material(material, &mut work)?;
    }
    let (vertices, edges) = graph.qualify(&mut work)?;
    work.charge(graph.curves.len() + graph.vertices.len() + graph.cells.len())?;
    for curve in &graph.curves {
        if curve.physical() {
            for &id in &curve.endpoints {
                graph.vertices[id as usize].boundary = true;
            }
        }
    }
    work.step()?;
    graph.point_owners = None;
    graph.curve_owners = None;
    // Copy only descriptors, MOVING exact boxed source values. Charge both old
    // and new descriptor buffers before allocation/publication; no BigInt or
    // oversized construction hash table survives in the immutable boundary.
    let requested = graph.vertices.len() * size_of::<BoundaryPoint>()
        + graph.curves.len() * size_of::<BoundaryCurve>()
        + graph.cells.len() * size_of::<BoundaryCell>();
    graph.check(requested)?;
    let header = size_of::<Boundary>() - size_of::<Graph>();
    let points = graph.buffer::<BoundaryPoint>(
        graph.vertices.len(),
        header
            + graph.curves.len() * size_of::<BoundaryCurve>()
            + graph.cells.len() * size_of::<BoundaryCell>(),
        &mut work,
    )?;
    let curves = graph.buffer::<BoundaryCurve>(
        graph.curves.len(),
        header
            + points.capacity() * size_of::<BoundaryPoint>()
            + graph.cells.len() * size_of::<BoundaryCell>(),
        &mut work,
    )?;
    let cells = graph.buffer::<BoundaryCell>(
        graph.cells.len(),
        header
            + points.capacity() * size_of::<BoundaryPoint>()
            + curves.capacity() * size_of::<BoundaryCurve>(),
        &mut work,
    )?;
    let (mut points, mut curves, mut cells) = (points, curves, cells);
    // The final buffer reservation included both earlier actual capacities,
    // this buffer's actual capacity and the Boundary header while all original
    // graph buffers remain live. Rechecking the same immutable sum adds no
    // allocation or publication protection.
    points.append(&mut graph.vertices);
    curves.append(&mut graph.curves);
    cells.append(&mut graph.cells);
    graph.vertices = points;
    graph.curves = curves;
    graph.cells = cells;
    graph.temporary = 0;
    work.step()?;
    Ok(Boundary {
        graph,
        work_used: work.used,
        vertices,
        edges,
    })
}
impl Boundary {
    // Public queries retain the Boundary extension beyond Graph. At most eight
    // dynamic Vec descriptors are live while tracing nested wires (including
    // construction/return descriptors); reserve them even when moves can share
    // storage. Their allocated capacities are charged separately by Graph.
    pub(super) fn query_headers(&self) -> usize {
        size_of::<Self>() - size_of::<Graph>() + 8 * size_of::<Vec<usize>>()
    }
    pub fn cavity(&self) -> &Cavity {
        &self.graph.source
    }
    pub fn material(&self) -> Option<&Material> {
        self.graph.material.as_ref()
    }
    pub fn face_count(&self) -> u32 {
        2 + 4 * self.graph.source.thread.spec().starts
            + u32::from(self.graph.source.entry.is_some())
            + 2 * u32::from(self.graph.source.tail.is_some())
            + u32::from(self.graph.material.is_some())
    }
    pub fn topology_counts(&self) -> (usize, usize) {
        (self.vertices, self.edges)
    }
    pub fn work_used(&self) -> usize {
        self.work_used
    }
    pub fn retained_bytes(&self) -> usize {
        self.graph.resident() + size_of::<Self>() - size_of::<Graph>()
    }
    pub fn points(&self) -> &[BoundaryPoint] {
        &self.graph.vertices
    }
    pub fn vertices(&self) -> impl Iterator<Item = (u32, &BoundaryPoint)> {
        self.graph
            .vertices
            .iter()
            .enumerate()
            .filter(|(_, p)| p.boundary)
            .map(|(i, p)| (i as u32, p))
    }
    pub fn curves(&self) -> impl Iterator<Item = (u32, &BoundaryCurve)> {
        self.graph
            .curves
            .iter()
            .enumerate()
            .filter(|(_, c)| c.physical())
            .map(|(i, c)| (i as u32, c))
    }
    pub fn cells(&self) -> &[BoundaryCell] {
        &self.graph.cells
    }
    pub fn face_loops(
        &self,
        face: u32,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Vec<Vec<BoundaryCoedge>>, Error> {
        let mut work = Work {
            remaining: self.graph.limits.work,
            used: 0,
            interrupt,
        };
        work.step()?;
        if face >= self.face_count() {
            return Err(Error::Invalid);
        }
        let result = self
            .graph
            .face_loops_reserved(face, self.query_headers(), &mut work)?;
        work.step()?;
        Ok(result)
    }
    /// Original source-cell q/v footpoint enclosure, retained from the actual
    /// exact values before shared geometry undergoes periodic canonicalization.
    pub fn cell_parameter_bounds(
        &self,
        cell: usize,
        corner: usize,
    ) -> Result<[Enclosure; 2], Error> {
        let cell = self.graph.cells.get(cell).ok_or(Error::Invalid)?;
        if corner >= cell.count {
            return Err(Error::Invalid);
        }
        Ok(cell.parameters[corner].map(interval))
    }
    pub fn curve_point(
        &self,
        index: u32,
        t: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<PointSample, Error> {
        if interrupt.is_some_and(|f| f()) {
            return Err(Error::Cancelled);
        }
        let curve = self
            .graph
            .curves
            .get(index as usize)
            .filter(|c| c.physical())
            .ok_or(Error::Invalid)?;
        if !t.is_finite() || !(0. ..=1.).contains(&t) {
            return Err(Error::Invalid);
        }
        let point = if t == 0. || t == 1. {
            self.graph.vertices[curve.endpoints[usize::from(t == 1.)] as usize].rzv
        } else {
            std::array::from_fn(|i| curve.base[i] + curve.delta[i] * I::point(t))
        };
        let result = sample(position(point, self.graph.source.thread.spec().left_hand));
        if interrupt.is_some_and(|f| f()) {
            return Err(Error::Cancelled);
        }
        Ok(result)
    }
}
impl BoundaryPoint {
    pub fn role(&self) -> BoundaryPointRole {
        match self.origin {
            Origin::Parameter { .. } => BoundaryPointRole::Source(self.source_key),
            Origin::Tail { end, .. } => BoundaryPointRole::Tail {
                original: self.source_key,
                end,
            },
            Origin::Stock { end } => BoundaryPointRole::StockRim { end },
        }
    }
    /// Opaque exact source parameter token; not a persistent import schema.
    pub fn exact_parameter_key(&self) -> &[u8] {
        &self.exact
    }
    pub fn position_bounds(&self) -> [Enclosure; 3] {
        position(self.rzv, self.left_hand).map(interval)
    }
    pub fn is_vertex(&self) -> bool {
        self.boundary
    }
    pub fn source_cell(&self) -> Option<(usize, usize)> {
        match self.origin {
            Origin::Parameter { cell, corner } => Some((cell, corner)),
            _ => None,
        }
    }
}
impl BoundaryCurve {
    /// Opaque exact affine (r,z,v) base/delta law. Orientation follows endpoints(),
    /// and is not claimed invariant under reversal or changes of source trims.
    pub fn exact_affine_key(&self) -> &[u8] {
        &self.exact
    }
    pub(super) fn physical(&self) -> bool {
        self.uses[1].is_some_and(|b| b.face != self.uses[0].unwrap().face)
    }
    pub fn endpoints(&self) -> [u32; 2] {
        self.endpoints
    }
    pub fn faces(&self) -> [u32; 2] {
        [self.uses[0].unwrap().face, self.uses[1].unwrap().face]
    }
    pub fn affine_rzv_bounds(&self) -> [[Enclosure; 3]; 2] {
        [self.base.map(interval), self.delta.map(interval)]
    }
}
impl BoundaryCell {
    pub fn period(&self) -> i64 {
        self.period
    }
    pub fn band(&self) -> usize {
        self.band
    }
    pub fn is_entry_cone(&self) -> bool {
        matches!(self.kind, Kind::Cone)
    }
    pub fn face(&self) -> u32 {
        self.face
    }
    pub fn points(&self) -> &[u32] {
        &self.vertices[..self.count]
    }
}
fn interval(v: I) -> Enclosure {
    Enclosure {
        lower: v.lo,
        upper: v.hi,
    }
}
fn position([r, z, v]: [I; 3], left: bool) -> [I; 3] {
    let tau = I::new(
        std::f64::consts::TAU.next_down(),
        std::f64::consts::TAU.next_up(),
    );
    let angle = v * if left { -tau } else { tau };
    [r * angle.cos(), r * angle.sin(), z]
}
fn sample(bounds: [I; 3]) -> PointSample {
    let point_mm = bounds.map(|v| v.lo.midpoint(v.hi));
    let error_mm = bounds.into_iter().zip(point_mm).fold(0., |s, (v, p)| {
        (s + (p - v.lo).next_up().max((v.hi - p).next_up())).next_up()
    });
    PointSample { point_mm, error_mm }
}
