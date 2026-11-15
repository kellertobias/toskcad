//! Actual literal degree11 rectangular shells. Structural and numeric binding
//! does not infer an authoring recipe or admit occupied material by itself.
use super::{
    owned_insert::CompoundDocument,
    part21::{self, Entity, Value},
};
use std::{collections::BTreeMap, mem::size_of};
mod material;
mod moments;
mod numeric;
mod selected;
pub use material::{MaterialOptions, Occupied, Projection};
pub use numeric::Scalar;
pub use selected::{
    FaceMeasure, FaceMeasureOptions, LiteralEdgePolyline, LiteralEmbeddingProof,
    LiteralGeometryOptions, LiteralGeometryPlan, LiteralMesh, LiteralMeshCandidate,
    LiteralMeshOptions, LiteralOccurrencePlacement, LiteralTransform, LiteralVertexSample,
    Location, MaterialQuery, OwnedMaterial, PlacedLiteralMesh, Placement, QueryOptions, Selected,
};
mod parse;
mod records;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    Invalid,
    Unsupported,
    Limit,
    Cancelled,
}
type Result<T> = std::result::Result<T, Error>;
#[derive(Clone, Copy, Debug)]
pub struct Options {
    pub maximum_work: usize,
    pub maximum_bytes: usize,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            maximum_work: 16_000_000,
            maximum_bytes: 64 * 1024 * 1024,
        }
    }
}
const SCRATCH: usize = 256 * 1024;
struct Budget<'a> {
    used: usize,
    remaining: usize,
    bytes: usize,
    peak: usize,
    limit: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl<'a> Budget<'a> {
    fn new(options: Options, interrupt: Option<&'a dyn Fn() -> bool>) -> Result<Self> {
        if options.maximum_work == 0
            || options.maximum_work > 16_000_000
            || options.maximum_bytes == 0
            || options.maximum_bytes > 64 * 1024 * 1024
        {
            return Err(Error::Invalid);
        }
        if SCRATCH > options.maximum_bytes {
            return Err(Error::Limit);
        }
        Ok(Self {
            used: 0,
            remaining: options.maximum_work,
            bytes: SCRATCH,
            peak: SCRATCH,
            limit: options.maximum_bytes,
            interrupt,
        })
    }
    fn step(&mut self, n: usize) -> Result<()> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(Error::Cancelled);
        }
        self.remaining = self.remaining.checked_sub(n).ok_or(Error::Limit)?;
        self.used += n;
        Ok(())
    }
    fn allocate<T>(&mut self, count: usize) -> Result<Vec<T>> {
        self.step(count + 1)?;
        let requested = count.checked_mul(size_of::<T>()).ok_or(Error::Limit)?;
        if requested > self.limit - self.bytes {
            return Err(Error::Limit);
        }
        let mut out = Vec::new();
        out.try_reserve_exact(count).map_err(|_| Error::Limit)?;
        self.bytes = self
            .bytes
            .checked_add(out.capacity() * size_of::<T>())
            .ok_or(Error::Limit)?;
        self.peak = self.peak.max(self.bytes);
        if self.bytes > self.limit {
            return Err(Error::Limit);
        }
        Ok(out)
    }
    fn release<T>(&mut self, values: &Vec<T>) {
        self.bytes -= values.capacity() * size_of::<T>();
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Point<'input> {
    id: u32,
    coordinates: [Scalar<'input>; 3],
}
impl<'input> Point<'input> {
    pub fn entity(&self) -> u32 {
        self.id
    }
    pub fn coordinates(&self) -> [Scalar<'input>; 3] {
        self.coordinates
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Vertex {
    id: u32,
    point: usize,
}
impl Vertex {
    pub fn entity(&self) -> u32 {
        self.id
    }
    pub fn point_index(&self) -> usize {
        self.point
    }
}
#[derive(Debug)]
pub struct Curve {
    id: u32,
    spans: usize,
    poles: Vec<usize>,
}
impl Curve {
    pub fn entity(&self) -> u32 {
        self.id
    }
    pub fn spans(&self) -> usize {
        self.spans
    }
    pub fn point_indices(&self) -> &[usize] {
        &self.poles
    }
}
#[derive(Debug)]
pub struct Tensor {
    id: u32,
    spans: [usize; 2],
    poles: Vec<usize>,
}
impl Tensor {
    pub fn entity(&self) -> u32 {
        self.id
    }
    pub fn spans(&self) -> [usize; 2] {
        self.spans
    }
    pub fn point_indices(&self) -> &[usize] {
        &self.poles
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Edge {
    id: u32,
    curve: usize,
    vertices: [usize; 2],
    pcurves: [u32; 2],
}
impl Edge {
    pub fn entity(&self) -> u32 {
        self.id
    }
    pub fn curve_index(&self) -> usize {
        self.curve
    }
    pub fn vertex_indices(&self) -> [usize; 2] {
        self.vertices
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Coedge {
    edge: usize,
    forward: bool,
}
impl Coedge {
    pub fn edge_index(&self) -> usize {
        self.edge
    }
    pub fn forward(&self) -> bool {
        self.forward
    }
}
#[derive(Clone, Copy, Debug)]
pub enum Support {
    Tensor(usize),
    Plane { height: usize },
}
#[derive(Debug)]
pub struct Face {
    id: u32,
    support: Support,
    coedges: [Coedge; 4],
    tensor_sides: Option<[usize; 4]>,
    forward: bool,
}
impl Face {
    pub fn entity(&self) -> u32 {
        self.id
    }
    pub fn support(&self) -> Support {
        self.support
    }
    pub fn coedges(&self) -> &[Coedge; 4] {
        &self.coedges
    }
    /// Original PCurve-to-tensor boundary binding in coedge order. Sides are
    /// 0: u=0, 1: v=1, 2: u=1, 3: v=0. Plane faces return None.
    /// A cyclic wire rotation preserves these validated owners.
    pub fn tensor_sides(&self) -> Option<[usize; 4]> {
        self.tensor_sides
    }
    pub fn forward(&self) -> bool {
        self.forward
    }
}
pub struct Parsed<'document, 'input> {
    document: &'document CompoundDocument<'input>,
    solid: u32,
    dependencies: Vec<u32>,
    points: Vec<Point<'input>>,
    vertices: Vec<Vertex>,
    curves: Vec<Curve>,
    tensors: Vec<Tensor>,
    edges: Vec<Edge>,
    faces: Vec<Face>,
    work: usize,
    bytes: usize,
}
impl<'document, 'input> Parsed<'document, 'input> {
    pub fn parse(
        document: &'document CompoundDocument<'input>,
        solid: u32,
        options: Options,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self> {
        parse::body(document, solid, options, interrupt)
    }
    pub fn document(&self) -> &'document CompoundDocument<'input> {
        self.document
    }
    pub fn solid_entity(&self) -> u32 {
        self.solid
    }
    pub fn points(&self) -> &[Point<'input>] {
        &self.points
    }
    pub fn vertices(&self) -> &[Vertex] {
        &self.vertices
    }
    pub fn curves(&self) -> &[Curve] {
        &self.curves
    }
    pub fn tensors(&self) -> &[Tensor] {
        &self.tensors
    }
    pub fn edges(&self) -> &[Edge] {
        &self.edges
    }
    pub fn faces(&self) -> &[Face] {
        &self.faces
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    /// Query arena only; the borrowed neutral document retains its additional
    /// original input/AST work and storage receipt independently.
    pub fn query_budget_bytes(&self) -> usize {
        self.bytes
    }
}

#[cfg(test)]
mod tests;
