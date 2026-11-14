//! Owned insert export candidates. These source-known carrier certificates do
//! not admit STEP solids, foreign imports, occupied material or decoded mass.
pub(crate) mod exact_cap;
use crate::{
    tapered_thread::interval::I,
    wood_insert::{Boundary, BoundaryCurveKey, BoundaryVertexKey, Side},
};
use std::mem::size_of;
mod curves;
mod domains;
pub(crate) mod projection;
pub(crate) mod topology;
mod writer;
pub use topology::{
    CapRegularity, EdgeGeometry, FaceGeometry, ShellCoedge, ShellEdge, ShellFace, TopologyPlan,
    WindingBounds,
};
pub use writer::{write_local_candidate, write_mapped_candidate};
mod angular;
mod surfaces;
pub(crate) use angular::fit_closed_unit;
pub(crate) use angular::{AngularCarrier, fit_angular};
pub use domains::{DomainOptions, DomainPlan, Slab, Station};
pub use projection::ProjectionBounds;
pub use surfaces::{RowKey, SurfacePatch, SurfacePlan, SurfaceRow, SurfaceVertex, VertexOwner};
pub(crate) mod resources;
use resources::{Memory, Work};
#[cfg(test)]
mod tests;
// Conservative peak stack-array charge: basis + endpoint/midpoint jets,
// ideal/literal controls and their nested interval temporaries total <16 KiB.
const WORKSPACE_BYTES: usize = 16 * 1024;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    Invalid,
    Unresolved,
    Limit,
    Cancelled,
}
#[derive(Clone, Copy, Debug)]
pub struct CurveOptions {
    pub position_error_mm: f64,
    /// Absolute first-derivative error in the ORIGINAL whole-edge parameter t.
    pub derivative_error_mm: f64,
    pub maximum_spans_per_curve: usize,
    pub work: usize,
    /// Original borrowed graph + owned output + live dynamic workspace.
    pub bytes: usize,
}
impl CurveOptions {
    pub(crate) fn validate(self) -> Result<(), Error> {
        if !self.position_error_mm.is_finite()
            || self.position_error_mm <= 0.
            || !self.derivative_error_mm.is_finite()
            || self.derivative_error_mm <= 0.
            || self.maximum_spans_per_curve == 0
            || self.maximum_spans_per_curve > 1024
            || self.bytes == 0
            || self.bytes > 64 * 1024 * 1024
            || self.work == 0
            || self.work > 16_000_000
        {
            Err(Error::Invalid)
        } else {
            Ok(())
        }
    }
}
impl Default for CurveOptions {
    fn default() -> Self {
        Self {
            position_error_mm: 1e-7,
            derivative_error_mm: 1e-6,
            maximum_spans_per_curve: 1024,
            work: 16_000_000,
            bytes: 64 * 1024 * 1024,
        }
    }
}
#[derive(Clone, Debug)]
pub struct CurveCarrier {
    key: BoundaryCurveKey,
    endpoints: [BoundaryVertexKey; 2],
    degree: usize,
    knots: Vec<f64>,
    poles: Vec<[f64; 3]>,
    position_error_mm: f64,
    derivative_error_mm: f64,
}
impl CurveCarrier {
    pub fn source_key(&self) -> BoundaryCurveKey {
        self.key
    }
    pub fn endpoints(&self) -> [BoundaryVertexKey; 2] {
        self.endpoints
    }
    pub fn degree(&self) -> usize {
        self.degree
    }
    pub fn knots(&self) -> &[f64] {
        &self.knots
    }
    pub fn poles(&self) -> &[[f64; 3]] {
        &self.poles
    }
    pub fn position_error_mm(&self) -> f64 {
        self.position_error_mm
    }
    pub fn derivative_error_mm(&self) -> f64 {
        self.derivative_error_mm
    }
}
pub struct CurvePlan<'a> {
    boundary: &'a Boundary,
    vertices: Vec<(BoundaryVertexKey, [f64; 3])>,
    curves: Vec<CurveCarrier>,
    work_used: usize,
    retained_bytes: usize,
}
impl<'a> CurvePlan<'a> {
    pub fn prepare(
        boundary: &'a Boundary,
        options: CurveOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        let mut work = Work::new(options.work, interrupt)?;
        work.charge(1)?;
        options.validate()?;
        let mut memory = Memory::new(boundary.retained_bytes(), options.bytes)?;
        memory.add(size_of::<Self>())?;
        memory.add(WORKSPACE_BYTES)?;
        let count = [Side::Outer, Side::Inner]
            .into_iter()
            .map(|s| boundary.source_boundary(s).points().len())
            .sum();
        let mut vertices = memory.allocate(count, &mut work)?;
        for side in [Side::Outer, Side::Inner] {
            for (i, p) in boundary.source_boundary(side).points().iter().enumerate() {
                work.charge(1)?;
                let bounds = p.position_bounds().map(|v| I::new(v.lower, v.upper));
                let point = bounds.map(|v| v.midpoint());
                vertices.push((vertex_key(side, i as u32), point));
            }
        }
        // Cap ownership, rather than a rounded coordinate comparison, selects
        // shared literal planar stations. Every incident carrier reuses them.
        for side in [Side::Outer, Side::Inner] {
            for (_, source) in boundary.source_boundary(side).source_curves() {
                work.charge(1)?;
                if let Some(cap) = source.faces().into_iter().find(|&f| f < 2) {
                    for endpoint in source.endpoints() {
                        let index = vertex_index(boundary, vertex_key(side, endpoint));
                        vertices[index].1[2] = if cap == 0 {
                            0.
                        } else {
                            boundary.source().full_length_mm()
                        };
                    }
                }
            }
        }
        let count = [Side::Outer, Side::Inner]
            .into_iter()
            .map(|s| boundary.source_boundary(s).source_curve_count())
            .sum();
        let mut carriers = memory.allocate(count, &mut work)?;
        let basis = curves::derivative_basis(&mut work)?;
        for side in [Side::Outer, Side::Inner] {
            for (i, source) in boundary.source_boundary(side).source_curves() {
                work.charge(1)?;
                let key = curve_key(side, i);
                let endpoints = source.endpoints().map(|v| vertex_key(side, v));
                let literal = endpoints.map(|key| vertex(&vertices, boundary, key));
                let cap = source.faces().into_iter().find(|&f| f < 2).map(|f| {
                    if f == 0 {
                        0.
                    } else {
                        boundary.source().full_length_mm()
                    }
                });
                let carrier = curves::prepare(
                    key,
                    endpoints,
                    source,
                    side == Side::Inner && boundary.source().inner().spec().left_hand
                        || side == Side::Outer && boundary.source().outer().spec().left_hand,
                    literal,
                    cap,
                    options,
                    &basis,
                    &mut memory,
                    &mut work,
                )?;
                carriers.push(carrier);
            }
        }
        work.charge(1)?;
        memory.used -= WORKSPACE_BYTES;
        Ok(Self {
            boundary,
            vertices,
            curves: carriers,
            work_used: work.used,
            retained_bytes: memory.used,
        })
    }
    /// Re-evaluate supplied literal controls against this SAME original source.
    /// This is not a foreign-file or material certificate.
    pub fn certify_controls(
        &self,
        index: usize,
        poles: &[[f64; 3]],
        maximum_work: usize,
        maximum_bytes: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<CurveCertificate, Error> {
        let mut work = Work::new(maximum_work, interrupt)?;
        work.charge(1)?;
        if maximum_bytes == 0 || maximum_bytes > 64 * 1024 * 1024 {
            return Err(Error::Invalid);
        }
        let curve = self.curves.get(index).ok_or(Error::Invalid)?;
        work.charge(poles.len().saturating_mul(3))?;
        if poles.len() != curve.poles.len() || poles.iter().flatten().any(|x| !x.is_finite()) {
            return Err(Error::Invalid);
        }
        if poles.first() != curve.poles.first() || poles.last() != curve.poles.last() {
            return Err(Error::Unresolved);
        }
        let mut memory = Memory::new(self.retained_bytes, maximum_bytes)?;
        memory.add(WORKSPACE_BYTES)?;
        if poles.as_ptr() != curve.poles.as_ptr() {
            memory.add(std::mem::size_of_val(poles))?;
        }
        let basis = curves::derivative_basis(&mut work)?;
        let side = curve.key.owner();
        let source = self.boundary.curve(curve.key).ok_or(Error::Invalid)?;
        let (position_error_mm, derivative_error_mm) = curves::recheck(
            source,
            self.boundary
                .source_boundary(side)
                .source()
                .spec()
                .left_hand,
            curve,
            poles,
            &basis,
            &mut work,
        )?;
        work.charge(1)?;
        Ok(CurveCertificate {
            position_error_mm,
            derivative_error_mm,
            work_used: work.used,
        })
    }
    pub fn boundary(&self) -> &'a Boundary {
        self.boundary
    }
    pub fn vertices(&self) -> &[(BoundaryVertexKey, [f64; 3])] {
        &self.vertices
    }
    /// Includes retained cell seams as well as semantic physical edges.
    pub fn curves(&self) -> &[CurveCarrier] {
        &self.curves
    }
    pub fn work_used(&self) -> usize {
        self.work_used
    }
    pub fn retained_bytes(&self) -> usize {
        self.retained_bytes
    }
}
fn vertex_key(side: Side, i: u32) -> BoundaryVertexKey {
    match side {
        Side::Outer => BoundaryVertexKey::Outer(i),
        Side::Inner => BoundaryVertexKey::Inner(i),
    }
}
fn curve_key(side: Side, i: u32) -> BoundaryCurveKey {
    match side {
        Side::Outer => BoundaryCurveKey::Outer(i),
        Side::Inner => BoundaryCurveKey::Inner(i),
    }
}
fn vertex(
    points: &[(BoundaryVertexKey, [f64; 3])],
    body: &Boundary,
    key: BoundaryVertexKey,
) -> [f64; 3] {
    points[vertex_index(body, key)].1
}
fn vertex_index(body: &Boundary, key: BoundaryVertexKey) -> usize {
    match key {
        BoundaryVertexKey::Outer(i) => i as usize,
        BoundaryVertexKey::Inner(i) => {
            body.source_boundary(Side::Outer).points().len() + i as usize
        }
    }
}

/// Same-source numerical carrier transfer, with no occupied-material authority.
#[derive(Clone, Copy, Debug)]
pub struct CurveCertificate {
    pub position_error_mm: f64,
    pub derivative_error_mm: f64,
    pub work_used: usize,
}
