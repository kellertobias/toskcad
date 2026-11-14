//! Actual coupled ruled-shell STEP data. No authoring recipe, inferred thread,
//! material mass or native Solid admission is produced by structural parsing.
use super::{
    Shell,
    part21::{self, Entity, Value},
};
use std::{collections::BTreeMap, mem::size_of};
mod bindings;
mod flux;
pub use flux::{FluxOptions, OccupiedVolume, OrientedFlux};
mod geometry;
mod numbers;
pub use bindings::{SupportBindings, TensorBinding};
pub use geometry::MaterialBounds;
use numbers::Numeric;
pub use numbers::Scalar;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    Invalid,
    Unsupported,
    Limit,
    Cancelled,
}
type Result<T> = std::result::Result<T, Error>;
impl From<crate::step::insert::Error> for Error {
    fn from(e: crate::step::insert::Error) -> Self {
        match e {
            crate::step::insert::Error::Invalid => Self::Invalid,
            crate::step::insert::Error::Unresolved => Self::Unsupported,
            crate::step::insert::Error::Limit => Self::Limit,
            crate::step::insert::Error::Cancelled => Self::Cancelled,
        }
    }
}
impl From<super::ReadError> for Error {
    fn from(e: super::ReadError) -> Self {
        match e {
            super::ReadError::Invalid => Self::Invalid,
            super::ReadError::Unsupported => Self::Unsupported,
            super::ReadError::Limit => Self::Limit,
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Options {
    pub maximum_work: usize,
    /// Retained input/AST/cache budget. The common bounded Part21 parser keeps
    /// its own fixed parsing limits; this is not an allocator peak measurement.
    pub retained_budget_bytes: usize,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            maximum_work: 16_000_000,
            retained_budget_bytes: 64 * 1024 * 1024,
        }
    }
}
struct Ledger<'a> {
    used: usize,
    remaining: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Ledger<'_> {
    fn step(&mut self, n: usize) -> Result<()> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(Error::Cancelled);
        }
        self.remaining = self.remaining.checked_sub(n).ok_or(Error::Limit)?;
        self.used = self.used.checked_add(n).ok_or(Error::Limit)?;
        Ok(())
    }
}
#[derive(Debug)]
pub struct PointLeaf<'i> {
    entity: u32,
    coordinates: Vec<Scalar<'i>>,
}
impl<'i> PointLeaf<'i> {
    pub fn entity(&self) -> u32 {
        self.entity
    }
    pub fn coordinates(&self) -> &[Scalar<'i>] {
        &self.coordinates
    }
}
#[derive(Debug)]
pub struct SplineCurve<'i> {
    entity: u32,
    degree: u32,
    poles: Vec<u32>,
    multiplicities: Vec<u32>,
    knots: Vec<Scalar<'i>>,
}
impl<'i> SplineCurve<'i> {
    pub fn entity(&self) -> u32 {
        self.entity
    }
    pub fn degree(&self) -> u32 {
        self.degree
    }
    pub fn pole_entities(&self) -> &[u32] {
        &self.poles
    }
    pub fn multiplicities(&self) -> &[u32] {
        &self.multiplicities
    }
    pub fn knots(&self) -> &[Scalar<'i>] {
        &self.knots
    }
}
#[derive(Debug)]
pub struct TensorSurface<'i> {
    entity: u32,
    rows: [Vec<u32>; 2],
    multiplicities_v: Vec<u32>,
    knots_v: Vec<Scalar<'i>>,
}
impl<'i> TensorSurface<'i> {
    pub fn entity(&self) -> u32 {
        self.entity
    }
    pub fn rows(&self) -> [&[u32]; 2] {
        self.rows.each_ref().map(Vec::as_slice)
    }
    pub fn multiplicities_v(&self) -> &[u32] {
        &self.multiplicities_v
    }
    pub fn knots_v(&self) -> &[Scalar<'i>] {
        &self.knots_v
    }
}
#[derive(Clone, Copy, Debug)]
pub struct PlaneSurface {
    entity: u32,
    axis: u32,
    origin: u32,
    normal: u32,
    reference: u32,
}
impl PlaneSurface {
    pub fn entity(&self) -> u32 {
        self.entity
    }
    pub fn axis_entity(&self) -> u32 {
        self.axis
    }
    pub fn origin_entity(&self) -> u32 {
        self.origin
    }
    pub fn direction_entities(&self) -> [u32; 2] {
        [self.normal, self.reference]
    }
}
#[derive(Clone, Copy, Debug)]
pub enum SurfaceIndex {
    Tensor(usize),
    Plane(usize),
}
#[derive(Clone, Copy, Debug)]
pub struct Coedge {
    pub edge: u32,
    pub forward: bool,
}
#[derive(Debug)]
pub struct Loop {
    bound: u32,
    entity: u32,
    outer: bool,
    coedges: Vec<Coedge>,
}
impl Loop {
    pub fn bound_entity(&self) -> u32 {
        self.bound
    }
    pub fn entity(&self) -> u32 {
        self.entity
    }
    pub fn is_outer(&self) -> bool {
        self.outer
    }
    pub fn coedges(&self) -> &[Coedge] {
        &self.coedges
    }
}
#[derive(Debug)]
pub struct Face {
    entity: u32,
    surface_entity: u32,
    surface: SurfaceIndex,
    same_sense: bool,
    loops: Vec<Loop>,
}
impl Face {
    pub fn entity(&self) -> u32 {
        self.entity
    }
    pub fn surface_entity(&self) -> u32 {
        self.surface_entity
    }
    pub fn surface(&self) -> SurfaceIndex {
        self.surface
    }
    pub fn same_sense(&self) -> bool {
        self.same_sense
    }
    pub fn loops(&self) -> &[Loop] {
        &self.loops
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Vertex {
    pub entity: u32,
    pub point_entity: u32,
}
#[derive(Clone, Copy, Debug)]
pub struct Pcurve {
    pub entity: u32,
    pub support_entity: u32,
    pub curve_entity: u32,
    pub representation: u32,
    pub context: u32,
}
#[derive(Debug)]
pub struct Edge {
    entity: u32,
    ends: [u32; 2],
    curve_entity: u32,
    surface_curve: u32,
    pcurves: [Pcurve; 2],
}
impl Edge {
    pub fn entity(&self) -> u32 {
        self.entity
    }
    pub fn vertex_entities(&self) -> [u32; 2] {
        self.ends
    }
    pub fn curve_entity(&self) -> u32 {
        self.curve_entity
    }
    pub fn surface_curve_entity(&self) -> u32 {
        self.surface_curve
    }
    pub fn pcurves(&self) -> &[Pcurve; 2] {
        &self.pcurves
    }
}
pub struct Parsed<'input> {
    input: &'input str,
    entities: BTreeMap<u32, Vec<Entity<'input>>>,
    points: Vec<PointLeaf<'input>>,
    curves: Vec<SplineCurve<'input>>,
    tensors: Vec<TensorSurface<'input>>,
    planes: Vec<PlaneSurface>,
    vertices: Vec<Vertex>,
    edges: Vec<Edge>,
    faces: Vec<Face>,
    work: usize,
    bytes: usize,
}
impl<'i> Parsed<'i> {
    pub fn literal_input(&self) -> &'i str {
        self.input
    }
    pub fn literal_entity_count(&self) -> usize {
        self.entities.len()
    }
    pub fn points(&self) -> &[PointLeaf<'i>] {
        &self.points
    }
    pub fn curves(&self) -> &[SplineCurve<'i>] {
        &self.curves
    }
    pub fn tensors(&self) -> &[TensorSurface<'i>] {
        &self.tensors
    }
    pub fn planes(&self) -> &[PlaneSurface] {
        &self.planes
    }
    pub fn vertices(&self) -> &[Vertex] {
        &self.vertices
    }
    pub fn edges(&self) -> &[Edge] {
        &self.edges
    }
    pub fn faces(&self) -> &[Face] {
        &self.faces
    }
    pub fn point(&self, id: u32) -> Option<&PointLeaf<'i>> {
        self.points
            .binary_search_by_key(&id, |p| p.entity)
            .ok()
            .map(|i| &self.points[i])
    }
    pub fn curve(&self, id: u32) -> Option<&SplineCurve<'i>> {
        self.curves
            .binary_search_by_key(&id, |p| p.entity)
            .ok()
            .map(|i| &self.curves[i])
    }
    pub fn edge(&self, id: u32) -> Option<&Edge> {
        self.edges
            .binary_search_by_key(&id, |e| e.entity)
            .ok()
            .map(|i| &self.edges[i])
    }
    pub fn vertex(&self, id: u32) -> Option<&Vertex> {
        self.vertices
            .binary_search_by_key(&id, |v| v.entity)
            .ok()
            .map(|i| &self.vertices[i])
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.bytes
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
}
fn simple<'a, 'i>(entities: &'a BTreeMap<u32, Vec<Entity<'i>>>, id: u32) -> Result<&'a Entity<'i>> {
    let [e] = entities.get(&id).ok_or(Error::Invalid)?.as_slice() else {
        return Err(Error::Unsupported);
    };
    Ok(e)
}
fn kind<'a, 'i>(
    entities: &'a BTreeMap<u32, Vec<Entity<'i>>>,
    id: u32,
    name: &str,
) -> Result<&'a Entity<'i>> {
    let e = simple(entities, id)?;
    if e.name != name {
        return Err(Error::Unsupported);
    }
    Ok(e)
}
fn reference(v: &Value<'_>) -> Result<u32> {
    super::reference(v).map_err(Error::from)
}
fn refs(v: &Value<'_>) -> Result<Vec<u32>> {
    super::refs(v).map_err(Error::from)
}
fn sense(v: &Value<'_>) -> Result<bool> {
    super::sense(v).map_err(Error::from)
}
fn list<'a, 'i>(v: &'a Value<'i>) -> Result<&'a [Value<'i>]> {
    if let Value::List(v) = v {
        Ok(v)
    } else {
        Err(Error::Invalid)
    }
}
mod budget;
mod schema;
/// Parse the actual represented graph and original carriers. This is an
/// unqualified structural packet, never a Solid, geometry-recipe match or
/// decoded material/volume receipt. Support/trim and embedding are separate.
pub fn parse<'input>(
    input: &'input str,
    options: Options,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Parsed<'input>> {
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
    w.step(input.len())?;
    let entities = part21::parse(input).map_err(Error::from)?;
    let retained = budget::arena(input, &entities, &mut w)?;
    if retained > options.retained_budget_bytes {
        return Err(Error::Limit);
    }
    let mut points = Vec::new();
    let mut curves = Vec::new();
    let mut tensors = Vec::new();
    let mut planes = Vec::new();
    for (&id, records) in &entities {
        w.step(1)?;
        for e in records {
            if e.name == "GEOMETRIC_REPRESENTATION_CONTEXT" {
                let [dimension] = e.args.as_slice() else {
                    return Err(Error::Invalid);
                };
                let mut n = Numeric::new(e, &mut w)?;
                if !matches!(n.integer(dimension)?, 2 | 3) {
                    return Err(Error::Unsupported);
                }
                n.finish()?;
            }
        }
        let [e] = records.as_slice() else { continue };
        match e.name {
            "CARTESIAN_POINT" => points.push(schema::point(id, e, &mut w)?),
            "B_SPLINE_CURVE_WITH_KNOTS" => curves.push(schema::curve(id, e, &mut w)?),
            "B_SPLINE_SURFACE_WITH_KNOTS" => tensors.push(schema::tensor(id, e, &mut w)?),
            "PLANE" => planes.push(schema::plane(&entities, id, e)?),
            _ => {}
        }
    }
    if planes.len() != 2 || tensors.is_empty() {
        return Err(Error::Unsupported);
    }
    // No source recipe or product label selects the body. Follow the actual
    // represented solid and retain its complete unchanged material face graph.
    let shell = Shell::new(&entities, 1e-7, false);
    let solids = shell.represented_solids(true).map_err(Error::from)?;
    if entities
        .values()
        .filter(|r| matches!(r.as_slice(),[e] if e.name=="MANIFOLD_SOLID_BREP"))
        .count()
        != 1
    {
        return Err(Error::Unsupported);
    }
    let mut faces = Vec::new();
    for face in shell.member_faces(solids[0]).map_err(Error::from)? {
        faces.push(schema::face(&entities, face, &tensors, &planes, &mut w)?);
    }
    let mut edge_ids = faces
        .iter()
        .flat_map(|f| f.loops.iter())
        .flat_map(|l| l.coedges.iter().map(|c| c.edge))
        .collect::<Vec<_>>();
    w.step(edge_ids.len() * 16)?;
    edge_ids.sort_unstable();
    edge_ids.dedup();
    let mut edges = Vec::new();
    for id in edge_ids {
        edges.push(schema::edge(&entities, id, &mut w)?);
    }
    let mut vertex_ids = edges.iter().flat_map(|e| e.ends).collect::<Vec<_>>();
    w.step(vertex_ids.len() * 16)?;
    vertex_ids.sort_unstable();
    vertex_ids.dedup();
    let mut vertices = Vec::new();
    for id in vertex_ids {
        w.step(1)?;
        let [Value::String(_), point] = kind(&entities, id, "VERTEX_POINT")?.args.as_slice() else {
            return Err(Error::Invalid);
        };
        vertices.push(Vertex {
            entity: id,
            point_entity: reference(point)?,
        });
    }
    let mut p = Parsed {
        input,
        entities,
        points,
        curves,
        tensors,
        planes,
        vertices,
        edges,
        faces,
        work: 0,
        bytes: 0,
    };
    validate_references(&p, &mut w)?;
    p.bytes = budget::packet(&p, &mut w)?;
    if p.bytes > options.retained_budget_bytes {
        return Err(Error::Limit);
    }
    w.step(1)?;
    p.work = w.used;
    Ok(p)
}

fn validate_references(p: &Parsed<'_>, w: &mut Ledger<'_>) -> Result<()> {
    for curve in &p.curves {
        let mut dimension = None;
        for &id in &curve.poles {
            w.step(32)?;
            let n = p.point(id).ok_or(Error::Invalid)?.coordinates.len();
            if dimension.is_some_and(|old| old != n) {
                return Err(Error::Invalid);
            }
            dimension = Some(n);
        }
    }
    for tensor in &p.tensors {
        for &id in tensor.rows.iter().flatten() {
            w.step(32)?;
            if p.point(id).ok_or(Error::Invalid)?.coordinates.len() != 3 {
                return Err(Error::Invalid);
            }
        }
    }
    for plane in &p.planes {
        if p.point(plane.origin)
            .ok_or(Error::Invalid)?
            .coordinates
            .len()
            != 3
        {
            return Err(Error::Invalid);
        }
        for id in [plane.normal, plane.reference] {
            let e = kind(&p.entities, id, "DIRECTION")?;
            let [Value::String(_), v] = e.args.as_slice() else {
                return Err(Error::Invalid);
            };
            let mut n = Numeric::new(e, w)?;
            if list(v)?.len() != 3 {
                return Err(Error::Unsupported);
            }
            for v in list(v)? {
                n.scalar(v)?;
            }
            n.finish()?;
        }
    }
    for vertex in &p.vertices {
        w.step(32)?;
        if p.point(vertex.point_entity)
            .ok_or(Error::Invalid)?
            .coordinates
            .len()
            != 3
        {
            return Err(Error::Invalid);
        }
    }
    for edge in &p.edges {
        w.step(256)?;
        let curve = p.curve(edge.curve_entity).ok_or(Error::Invalid)?;
        if p.point(curve.poles[0])
            .ok_or(Error::Invalid)?
            .coordinates
            .len()
            != 3
        {
            return Err(Error::Invalid);
        }
        for (i, id) in edge.ends.iter().enumerate() {
            let vertex = p
                .vertices
                .binary_search_by_key(id, |v| v.entity)
                .map_err(|_| Error::Invalid)?;
            if p.vertices[vertex].point_entity
                != curve.poles[if i == 0 { 0 } else { curve.poles.len() - 1 }]
            {
                return Err(Error::Unsupported);
            }
        }
        for pc in &edge.pcurves {
            let curve = p.curve(pc.curve_entity).ok_or(Error::Invalid)?;
            if p.point(curve.poles[0])
                .ok_or(Error::Invalid)?
                .coordinates
                .len()
                != 2
            {
                return Err(Error::Invalid);
            }
            let context = p.entities.get(&pc.context).ok_or(Error::Invalid)?;
            let e = context
                .iter()
                .find(|e| e.name == "GEOMETRIC_REPRESENTATION_CONTEXT")
                .ok_or(Error::Unsupported)?;
            let [dimension] = e.args.as_slice() else {
                return Err(Error::Invalid);
            };
            let mut n = Numeric::new(e, w)?;
            if n.integer(dimension)? != 2 {
                return Err(Error::Unsupported);
            }
            n.finish()?;
        }
    }
    Ok(())
}
#[cfg(test)]
mod packet_tests;
#[cfg(test)]
mod schema_tests;
#[cfg(test)]
mod tests;
