//! Original finite thread carriers compose a paired material boundary.
//! Component cap disks are source scaffolding only: the occupied material has
//! exactly two plane annuli with the actual outer/inner source rim wires.
use super::*;
use crate::authored_internal_treatment::{
    BoundaryCell, BoundaryCurve, BoundaryPoint, ThreadBoundary,
};
use resources::{Memory, Work, failure, loop_bytes};
use std::mem::size_of;
mod resources;
#[cfg(test)]
mod tests;
mod topology;
/// Packet-local identity, qualified by the actual component source. The point's
/// original constraint role remains available through point().role().
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BoundaryVertexKey {
    Outer(u32),
    Inner(u32),
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BoundaryCurveKey {
    Outer(u32),
    Inner(u32),
}
impl BoundaryCurveKey {
    pub fn owner(self) -> Side {
        match self {
            Self::Outer(_) => Side::Outer,
            Self::Inner(_) => Side::Inner,
        }
    }
    pub fn source_index(self) -> u32 {
        match self {
            Self::Outer(i) | Self::Inner(i) => i,
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BoundaryCoedge {
    pub curve: BoundaryCurveKey,
    pub forward: bool,
}
pub struct Boundary {
    source: Insert,
    outer: ThreadBoundary,
    inner: ThreadBoundary,
    loops: Vec<Vec<Vec<BoundaryCoedge>>>,
    limits: BoundaryOptions,
    work_used: usize,
    vertices: usize,
    edges: usize,
}
impl Insert {
    /// Source incidence only. No mesh, file import, mass or Native Solid
    /// authority is inferred from this independently qualified finite graph.
    pub fn boundary_with_interrupt(
        &self,
        options: BoundaryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Boundary, Error> {
        checkpoint(interrupt)?;
        options.validate().map_err(failure)?;
        let mut w = Work {
            remaining: options.work,
            used: 0,
            interrupt,
        };
        w.charge(1)?;
        let header = size_of::<Boundary>() - 2 * size_of::<ThreadBoundary>();
        let mut memory = Memory {
            used: header,
            limit: options.bytes,
        };
        if w.remaining == 0 {
            return Err(Error::Limit);
        }
        let outer = self
            .outer()
            .boundary_with_interrupt(
                BoundaryOptions {
                    bytes: options
                        .bytes
                        .checked_sub(header)
                        .filter(|&v| v > 0)
                        .ok_or(Error::Limit)?,
                    work: w.remaining,
                    ..options
                },
                interrupt,
            )
            .map_err(failure)?;
        w.charge(outer.work_used())?;
        memory.add(outer.retained_bytes())?;
        let inner_options = BoundaryOptions {
            cells: options
                .cells
                .checked_sub(outer.cells().len())
                .filter(|&v| v > 0)
                .ok_or(Error::Limit)?,
            vertices: options
                .vertices
                .checked_sub(outer.points().len())
                .filter(|&v| v > 0)
                .ok_or(Error::Limit)?,
            curves: options
                .curves
                .checked_sub(outer.source_curve_count())
                .filter(|&v| v > 0)
                .ok_or(Error::Limit)?,
            bytes: options
                .bytes
                .checked_sub(memory.used)
                .filter(|&v| v > 0)
                .ok_or(Error::Limit)?,
            work: w.remaining,
        };
        if inner_options.work == 0 {
            return Err(Error::Limit);
        }
        let inner = self
            .inner()
            .boundary_with_interrupt(inner_options, interrupt)
            .map_err(failure)?;
        w.charge(inner.work_used())?;
        memory.add(inner.retained_bytes())?;
        let mut result = Boundary {
            source: *self,
            outer,
            inner,
            loops: Vec::new(),
            limits: options,
            work_used: 0,
            vertices: 0,
            edges: 0,
        };
        result.build_loops(&mut memory, &mut w)?;
        let (vertices, edges) = result.qualify(&mut memory, &mut w)?;
        result.vertices = vertices;
        result.edges = edges;
        w.charge(1)?;
        result.work_used = w.used;
        Ok(result)
    }
}
impl Boundary {
    pub fn source(&self) -> &Insert {
        &self.source
    }
    /// Original closed component source. Its disk caps are not occupied
    /// material faces; use this Boundary's face_loops for the paired annuli.
    pub fn source_boundary(&self, owner: Side) -> &ThreadBoundary {
        match owner {
            Side::Outer => &self.outer,
            Side::Inner => &self.inner,
        }
    }
    pub fn face_count(&self) -> u32 {
        self.source.face_count()
    }
    pub fn work_used(&self) -> usize {
        self.work_used
    }
    pub fn topology_counts(&self) -> (usize, usize) {
        (self.vertices, self.edges)
    }
    pub fn retained_bytes(&self) -> usize {
        size_of::<Self>() - 2 * size_of::<ThreadBoundary>()
            + self.outer.retained_bytes()
            + self.inner.retained_bytes()
            + self.loops.capacity() * size_of::<Vec<Vec<BoundaryCoedge>>>()
            + self.loops.iter().map(loop_bytes).sum::<usize>()
    }
    pub fn face_loops(&self, face: FaceKey) -> Option<&[Vec<BoundaryCoedge>]> {
        self.source
            .face_index(face)
            .and_then(|i| self.loops.get(i as usize))
            .map(Vec::as_slice)
    }
    pub fn curve(&self, key: BoundaryCurveKey) -> Option<&BoundaryCurve> {
        self.source_boundary(key.owner())
            .source_curve(key.source_index())
    }
    pub fn point(&self, key: BoundaryVertexKey) -> Option<&BoundaryPoint> {
        match key {
            BoundaryVertexKey::Outer(i) => self.outer.points().get(i as usize),
            BoundaryVertexKey::Inner(i) => self.inner.points().get(i as usize),
        }
    }
    pub fn vertices(&self) -> impl Iterator<Item = (BoundaryVertexKey, &BoundaryPoint)> {
        self.outer
            .vertices()
            .map(|(i, p)| (BoundaryVertexKey::Outer(i), p))
            .chain(
                self.inner
                    .vertices()
                    .map(|(i, p)| (BoundaryVertexKey::Inner(i), p)),
            )
    }
    pub fn curves(&self) -> impl Iterator<Item = (BoundaryCurveKey, &BoundaryCurve)> {
        self.outer
            .curves()
            .map(|(i, c)| (BoundaryCurveKey::Outer(i), c))
            .chain(
                self.inner
                    .curves()
                    .map(|(i, c)| (BoundaryCurveKey::Inner(i), c)),
            )
    }
    pub fn curve_endpoints(&self, key: BoundaryCurveKey) -> Option<[BoundaryVertexKey; 2]> {
        self.curve(key)
            .map(|c| c.endpoints().map(|i| vertex_key(key.owner(), i)))
    }
    pub fn curve_faces(&self, key: BoundaryCurveKey) -> Option<[FaceKey; 2]> {
        self.curve(key).map(|c| {
            c.faces().map(|f| {
                self.source
                    .face_key(self.global_face(key.owner(), f))
                    .expect("sealed source face")
            })
        })
    }
    pub fn curve_point(
        &self,
        key: BoundaryCurveKey,
        t: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<crate::authored_flat_thread::PointSample, Error> {
        self.source_boundary(key.owner())
            .curve_point(key.source_index(), t, interrupt)
            .map_err(failure)
    }
    pub fn cells(&self, owner: Side) -> &[BoundaryCell] {
        self.source_boundary(owner).cells()
    }
    pub fn cell_surface_same_sense(&self, owner: Side) -> bool {
        self.source_boundary(owner).cell_surface_same_sense() != (owner == Side::Inner)
    }
    pub fn cell_face(&self, owner: Side, cell: usize) -> Option<FaceKey> {
        self.cells(owner)
            .get(cell)
            .and_then(|c| self.source.face_key(self.global_face(owner, c.face())))
    }
    pub fn cell_coedges(
        &self,
        owner: Side,
        cell: usize,
    ) -> Result<impl Iterator<Item = BoundaryCoedge> + '_, Error> {
        let mut edges = self
            .source_boundary(owner)
            .cell_coedges(cell)
            .map_err(failure)?;
        Ok(std::iter::from_fn(move || {
            let edge = if owner == Side::Inner {
                edges.next_back()
            } else {
                edges.next()
            }?;
            Some(BoundaryCoedge {
                curve: curve_key(owner, edge.curve as u32),
                forward: edge.forward != (owner == Side::Inner),
            })
        }))
    }
    fn global_face(&self, owner: Side, face: u32) -> u32 {
        if face < 2 || owner == Side::Outer {
            face
        } else {
            face + 4 * self.source.outer().spec().starts
        }
    }
    fn build_loops(&mut self, memory: &mut Memory, w: &mut Work<'_>) -> Result<(), Error> {
        self.loops = memory.allocate(self.face_count() as usize, w)?;
        for face in 0..self.face_count() {
            let cap = face < 2;
            let mut loops = memory.allocate(if cap { 2 } else { 1 }, w)?;
            let (owner, source_face) = if cap || face < 2 + 4 * self.source.outer().spec().starts {
                (Side::Outer, face)
            } else {
                (Side::Inner, face - 4 * self.source.outer().spec().starts)
            };
            self.append_loops(owner, source_face, &mut loops, memory, w)?;
            if cap {
                self.append_loops(Side::Inner, face, &mut loops, memory, w)?;
            }
            self.loops.push(loops);
        }
        Ok(())
    }
    fn append_loops(
        &self,
        owner: Side,
        face: u32,
        out: &mut Vec<Vec<BoundaryCoedge>>,
        memory: &mut Memory,
        w: &mut Work<'_>,
    ) -> Result<(), Error> {
        let source = self.source_boundary(owner);
        // Each component constructor already reserved the then-live packet.
        // Add only later live objects here; do not double-charge that reservation.
        let construction_reserved = self.limits.bytes - source.query_byte_limit();
        let external = memory
            .used
            .checked_sub(source.retained_bytes())
            .and_then(|v| v.checked_sub(construction_reserved))
            .ok_or(Error::Limit)?;
        let (loops, used) = source
            .face_loops_accounted(face, w.remaining, external, w.interrupt)
            .map_err(failure)?;
        w.charge(used)?;
        let temporary = loop_bytes(&loops);
        memory.add(temporary)?;
        for chain in &loops {
            let mut owned = memory.allocate(chain.len(), w)?;
            for index in 0..chain.len() {
                w.charge(1)?;
                let edge = chain[if owner == Side::Inner {
                    chain.len() - 1 - index
                } else {
                    index
                }];
                owned.push(BoundaryCoedge {
                    curve: curve_key(owner, edge.curve as u32),
                    forward: edge.forward != (owner == Side::Inner),
                });
            }
            out.push(owned);
        }
        drop(loops);
        memory.used -= temporary;
        Ok(())
    }
}
fn curve_key(owner: Side, id: u32) -> BoundaryCurveKey {
    match owner {
        Side::Outer => BoundaryCurveKey::Outer(id),
        Side::Inner => BoundaryCurveKey::Inner(id),
    }
}
fn vertex_key(owner: Side, id: u32) -> BoundaryVertexKey {
    match owner {
        Side::Outer => BoundaryVertexKey::Outer(id),
        Side::Inner => BoundaryVertexKey::Inner(id),
    }
}
