//! Literal shared rows for local insert tensors. No closed/material/file claim.
use super::*;
use super::{
    curves::{Basis, Fixed, Law},
    domains::N,
};
mod source;
#[cfg(test)]
mod tests;
const SCRATCH: usize = 16 * 1024 * 1024;
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct RowKey {
    pub curve: BoundaryCurveKey,
    pub interval: u32,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VertexOwner {
    Original(BoundaryVertexKey),
    Station {
        curve: BoundaryCurveKey,
        station: u32,
    },
}
pub struct SurfaceVertex {
    owner: VertexOwner,
    point: [f64; 3],
}
impl SurfaceVertex {
    pub fn owner(&self) -> VertexOwner {
        self.owner
    }
    pub fn point_mm(&self) -> [f64; 3] {
        self.point
    }
}
pub struct SurfaceRow {
    key: RowKey,
    ends: [usize; 2],
    fit: Fixed,
}
impl SurfaceRow {
    pub fn key(&self) -> RowKey {
        self.key
    }
    pub fn endpoints(&self) -> [usize; 2] {
        self.ends
    }
    pub fn poles(&self) -> &[[f64; 3]] {
        &self.fit.poles
    }
    pub fn knots(&self) -> &[f64] {
        &self.fit.knots
    }
    pub fn position_error_mm(&self) -> f64 {
        self.fit.position_error_mm
    }
    /// Derivative in the WHOLE original angular-interval coordinate t in[0,1].
    pub fn derivative_error_mm(&self) -> f64 {
        self.fit.derivative_error_mm
    }
}
#[derive(Clone, Copy, Debug)]
pub struct SurfacePatch {
    pub slab: usize,
    pub rows: [usize; 2],
    pub same_sense: bool,
}
/// Every tensor is u-degree1, v-degree11, with its actual two owning rows.
/// Adjacent patches borrow the SAME row object; no independently fitted seam.
pub struct SurfacePlan<'p, 'a> {
    domain: &'p DomainPlan<'a>,
    vertices: Vec<SurfaceVertex>,
    rows: Vec<SurfaceRow>,
    patches: Vec<SurfacePatch>,
    work: usize,
    bytes: usize,
}
impl<'p, 'a> SurfacePlan<'p, 'a> {
    pub fn domain(&self) -> &'p DomainPlan<'a> {
        self.domain
    }
    pub fn vertices(&self) -> &[SurfaceVertex] {
        &self.vertices
    }
    pub fn rows(&self) -> &[SurfaceRow] {
        &self.rows
    }
    pub fn patches(&self) -> &[SurfacePatch] {
        &self.patches
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_bytes(&self) -> usize {
        self.bytes
    }
    /// Aggregate work/bytes include the borrowed DomainPlan (and original graph)
    /// once. This admits row correspondence only, never material or STEP.
    pub fn prepare(
        domain: &'p DomainPlan<'a>,
        options: CurveOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        let mut w = Work::new(options.work, interrupt)?;
        w.charge(1)?;
        options.validate()?;
        w.charge(domain.work_used())?;
        let mut memory = Memory::new(domain.retained_bytes(), options.bytes)?;
        memory.add(size_of::<Self>() + SCRATCH + WORKSPACE_BYTES)?;
        let mut keys: Vec<RowKey> = memory.allocate(domain.slabs().len() * 2, &mut w)?;
        for slab in domain.slabs() {
            for curve in slab.rows {
                let key = RowKey {
                    curve,
                    interval: slab.stations[0],
                };
                let mut low = 0;
                let mut high = keys.len();
                let mut found = false;
                while low < high {
                    w.charge(1)?;
                    let mid = low + (high - low) / 2;
                    match keys[mid].cmp(&key) {
                        std::cmp::Ordering::Less => low = mid + 1,
                        std::cmp::Ordering::Greater => high = mid,
                        std::cmp::Ordering::Equal => {
                            found = true;
                            break;
                        }
                    }
                }
                if !found {
                    w.charge(1 + keys.len() - low)?;
                    keys.insert(low, key);
                }
            }
        }
        let count = [Side::Outer, Side::Inner]
            .into_iter()
            .map(|s| domain.boundary().source_boundary(s).points().len())
            .sum::<usize>();
        let mut vertices = memory.allocate(count + keys.len() * 2, &mut w)?;
        source::original_vertices(domain, &mut vertices, &mut w)?;
        let mut requests = memory.allocate(keys.len(), &mut w)?;
        for key in &keys {
            let previous = requests.last();
            let request = source::request(
                domain,
                *key,
                previous,
                &mut vertices,
                &mut memory,
                &mut w,
                interrupt,
            )?;
            requests.push(request);
        }
        let mut spans = [
            memory.allocate(domain.stations(Side::Outer).len() - 1, &mut w)?,
            memory.allocate(domain.stations(Side::Inner).len() - 1, &mut w)?,
        ];
        for (i, side) in [Side::Outer, Side::Inner].into_iter().enumerate() {
            spans[i].resize(domain.stations(side).len() - 1, 1usize);
        }
        let basis = curves::derivative_basis(&mut w)?;
        let rows = fit_rows(
            &requests,
            &vertices,
            &mut spans,
            options,
            &basis,
            &mut memory,
            &mut w,
        )?;
        let mut patches = memory.allocate(domain.slabs().len(), &mut w)?;
        for (slab_index, slab) in domain.slabs().iter().enumerate() {
            w.charge(6 * (usize::BITS - keys.len().leading_zeros() + 1) as usize)?;
            let row_ids = slab.rows.map(|curve| {
                keys.binary_search(&RowKey {
                    curve,
                    interval: slab.stations[0],
                })
                .expect("registered exact source row")
            });
            let [a, b] = row_ids.map(|i| &rows[i]);
            verify_boundaries(slab, a, b, &vertices)?;
            let left = domain
                .boundary()
                .source_boundary(slab.side)
                .source()
                .spec()
                .left_hand;
            patches.push(SurfacePatch {
                slab: slab_index,
                rows: row_ids,
                same_sense: left != (slab.side == Side::Inner),
            });
        }
        memory.release(&requests);
        memory.release(&keys);
        for s in &spans {
            memory.release(s);
        }
        memory.used -= SCRATCH + WORKSPACE_BYTES;
        w.charge(1)?;
        Ok(Self {
            domain,
            vertices,
            rows,
            patches,
            work: w.used,
            bytes: memory.used,
        })
    }
}
struct Request {
    key: RowKey,
    ends: [usize; 2],
    law: Law,
    cap: Option<f64>,
}
fn side_index(side: Side) -> usize {
    match side {
        Side::Outer => 0,
        Side::Inner => 1,
    }
}
fn fit_rows(
    requests: &[Request],
    vertices: &[SurfaceVertex],
    spans: &mut [Vec<usize>; 2],
    options: CurveOptions,
    basis: &Basis,
    memory: &mut Memory,
    w: &mut Work<'_>,
) -> Result<Vec<SurfaceRow>, Error> {
    let mut rows = memory.allocate(requests.len(), w)?;
    loop {
        let mut failed = [false; 2];
        for r in requests {
            let i = side_index(r.key.curve.owner());
            let n = spans[i][r.key.interval as usize];
            let ends = r.ends.map(|v| vertices[v].point);
            let fit = curves::fixed(&r.law, n, ends, r.cap, basis, memory, w)?;
            if fit.position_error_mm > options.position_error_mm
                || fit.derivative_error_mm > options.derivative_error_mm
            {
                if n > options.maximum_spans_per_curve / 2 {
                    return Err(Error::Unresolved);
                }
                // Refine only after ALL rows use the same current interval grid.
                failed[i] = true;
            }
            rows.push(SurfaceRow {
                key: r.key,
                ends: r.ends,
                fit,
            });
        }
        if failed == [false, false] {
            return Ok(rows);
        }
        // Only intervals whose actual rows failed are refined. A successful
        // neighboring interval is not a certificate for this one.
        for row in &rows {
            if row.fit.position_error_mm > options.position_error_mm
                || row.fit.derivative_error_mm > options.derivative_error_mm
            {
                let i = side_index(row.key.curve.owner());
                let n = &mut spans[i][row.key.interval as usize];
                let original = row.fit.knots.len() - 1;
                if *n == original {
                    *n *= 2;
                }
            }
        }
        for row in rows.drain(..) {
            memory.release(&row.fit.knots);
            memory.release(&row.fit.poles);
        }
    }
}

pub(super) fn verify_boundaries(
    slab: &Slab,
    a: &SurfaceRow,
    b: &SurfaceRow,
    vertices: &[SurfaceVertex],
) -> Result<(), Error> {
    if a.knots() != b.knots() {
        return Err(Error::Unresolved);
    }
    for row in [a, b] {
        for end in 0..2 {
            let pole = if end == 0 {
                row.poles().first()
            } else {
                row.poles().last()
            };
            if pole != vertices.get(row.ends[end]).map(|v| &v.point) {
                return Err(Error::Unresolved);
            }
        }
    }
    // U-degree1: its COMPLETE boundary row consists of precisely these two
    // owned controls. Equal spatial coordinates alone never erase an edge.
    for end in 0..2 {
        if slab.collapsed[end] && a.ends[end] != b.ends[end] {
            return Err(Error::Unresolved);
        }
    }
    Ok(())
}
