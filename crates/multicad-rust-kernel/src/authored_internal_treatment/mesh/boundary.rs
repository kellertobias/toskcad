//! Actual unsampled source cells and source-owned finite boundary curves.
//! Display/manufacturing triangles never supply geometry or incidence.
use super::{
    clip::Point,
    exact::Q,
    owners::{Constraint, Key, Owners},
    regular::Patch,
    *,
};
use std::{cmp::Ordering, mem::size_of};
mod api;
mod close;
mod loops;
mod material;
mod polyline;
mod qualify;
#[cfg(test)]
mod tests;
const BYTES: usize = 64 * 1024 * 1024;
// Construction keeps only one Patch/Polygon and one exact registration alive.
// These contain fewer than 128 distinct Q values (shared Rc clones do not add
// limbs). Each Q operand is <=8192 bits; arithmetic preflights products before
// allocation. The pinned bigint algorithms use bounded school/Karatsuba scratch
// at this size, including <=16384-bit comparison/reduction temporaries. Reserve
// 8 MiB for these <=128 values, arithmetic scratch, encoded staging and the
// <=6-line/6-point polygon buffers, rather than charging reduced bit length as
// historical Vec capacity. No BigInt survives the construction/query phase.
const EXACT_SCRATCH: usize = 8 * 1024 * 1024;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BoundaryOptions {
    pub cells: usize,
    pub vertices: usize,
    pub curves: usize,
    pub work: usize,
    /// Aggregate owned graph, query output and dynamic workspace bytes,
    /// including container descriptors. Fixed call-stack frames are not an
    /// allocator budget; exact-arithmetic dynamic scratch has its own reserve.
    pub bytes: usize,
}
impl Default for BoundaryOptions {
    fn default() -> Self {
        Self {
            cells: 2048,
            vertices: 8192,
            curves: 8192,
            work: 16_000_000,
            bytes: BYTES,
        }
    }
}
#[derive(Clone, Copy)]
enum Origin {
    Parameter { cell: usize, corner: usize },
    Tail { source: u32, end: u8 },
    Stock { end: u8 },
}
/// Exact unsampled authored constraints. Rail n is the global phase knot
/// period*5+band (including the next period at band5); Row0 is the angular cut.
/// Station codes: 0 mouth, 1 entry/thread seam, 2 thread end. A nonzero
/// Station domain identifies the actual phase band crossing the angular cut.
/// These identify source constraints, not immutable fragments across all edits.
pub use super::owners::{Constraint as BoundaryConstraint, Key as BoundaryPointKey};
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BoundaryPointRole {
    Source(BoundaryPointKey),
    Tail { original: BoundaryPointKey, end: u8 },
    StockRim { end: u8 },
}
pub struct BoundaryPoint {
    exact: Box<[u8]>,
    source_key: Key,
    rzv: [I; 3],
    origin: Origin,
    boundary: bool,
    left_hand: bool,
}
#[derive(Clone, Copy)]
struct Use {
    face: u32,
    forward: bool,
}
pub struct BoundaryCurve {
    exact: Box<[u8]>,
    endpoints: [u32; 2],
    base: [I; 3],
    delta: [I; 3],
    v_direction: Ordering,
    wrap_jump: i8,
    station: Option<u8>,
    uses: [Option<Use>; 2],
}
/// Outward coefficients of the actual normalized original source chart:
/// Q=q_origin+q_span*u, r=r0+ru*u+rv*v, z=z0+zu*u+zv*v.
/// Angular v spans one exact true revolution; physical trims remain separate.
#[derive(Clone, Copy, Debug)]
pub struct BoundaryCellSource {
    q: [I; 2],
    radial: [I; 3],
    axial: [I; 3],
}
pub struct BoundaryCell {
    period: i64,
    band: usize,
    kind: Kind,
    face: u32,
    vertices: [u32; 6],
    count: usize,
    parameters: [[I; 2]; 6],
    source: BoundaryCellSource,
    coedges: [BoundaryCoedge; 6],
}
struct Graph {
    source: Cavity,
    material: Option<Material>,
    vertices: Vec<BoundaryPoint>,
    curves: Vec<BoundaryCurve>,
    cells: Vec<BoundaryCell>,
    payload: usize,
    point_owners: Option<Owners>,
    curve_owners: Option<Owners>,
    limits: BoundaryOptions,
    temporary: usize,
}
impl Graph {
    fn point_owner_bytes(&self) -> usize {
        self.point_owners.as_ref().map_or(0, Owners::memory_bytes)
    }
    fn curve_owner_bytes(&self) -> usize {
        self.curve_owners.as_ref().map_or(0, Owners::memory_bytes)
    }
    fn resident(&self) -> usize {
        size_of::<Self>()
            + self.temporary
            + self.payload
            + self.vertices.capacity() * size_of::<BoundaryPoint>()
            + self.curves.capacity() * size_of::<BoundaryCurve>()
            + self.cells.capacity() * size_of::<BoundaryCell>()
    }
    fn check(&self, additional: usize) -> Result<(), Error> {
        storage::bytes(
            &[(self.point_owner_bytes(), 1), (self.curve_owner_bytes(), 1)],
            self.resident()
                .checked_add(additional)
                .ok_or(Error::Limit)?,
            self.limits.bytes,
        )
        .map(|_| ())
    }
    fn buffer<T>(&self, count: usize, reserved: usize, w: &mut Work<'_>) -> Result<Vec<T>, Error> {
        w.step()?;
        let requested = reserved
            .checked_add(count.checked_mul(size_of::<T>()).ok_or(Error::Limit)?)
            .ok_or(Error::Limit)?;
        self.check(requested)?;
        let mut result = Vec::new();
        result.try_reserve_exact(count).map_err(|_| Error::Limit)?;
        let actual = reserved
            .checked_add(result.capacity() * size_of::<T>())
            .ok_or(Error::Limit)?;
        self.check(actual)?;
        w.step()?;
        Ok(result)
    }
    fn new(source: Cavity, limits: BoundaryOptions, w: &mut Work<'_>) -> Result<Self, Error> {
        w.step()?;
        // Both entry points validate immutable BoundaryOptions before work or
        // allocation; no second, unreachable options gate is needed here.
        let point_bytes = Owners::requested_bytes(limits.vertices)?;
        let curve_bytes = Owners::requested_bytes(limits.curves)?;
        let mut graph = Self {
            source,
            material: None,
            vertices: Vec::new(),
            curves: Vec::new(),
            cells: Vec::new(),
            payload: 0,
            point_owners: None,
            curve_owners: None,
            limits,
            temporary: EXACT_SCRATCH,
        };
        let vp = limits.vertices * size_of::<BoundaryPoint>();
        let cp = limits.curves * size_of::<BoundaryCurve>();
        let sp = limits.cells * size_of::<BoundaryCell>();
        graph.check(vp + cp + sp + point_bytes + curve_bytes)?;
        graph
            .vertices
            .try_reserve_exact(limits.vertices)
            .map_err(|_| Error::Limit)?;
        graph.check(cp + sp + point_bytes + curve_bytes)?;
        graph
            .curves
            .try_reserve_exact(limits.curves)
            .map_err(|_| Error::Limit)?;
        graph.check(sp + point_bytes + curve_bytes)?;
        graph
            .cells
            .try_reserve_exact(limits.cells)
            .map_err(|_| Error::Limit)?;
        graph.check(point_bytes + curve_bytes)?;
        graph.point_owners = Some(Owners::new_reserved(
            limits.vertices,
            graph.resident() + curve_bytes,
            limits.bytes,
            w,
        )?);
        graph.curve_owners = Some(Owners::new_reserved(
            limits.curves,
            graph.resident() + graph.point_owner_bytes(),
            limits.bytes,
            w,
        )?);
        graph.check(0)?;
        Ok(graph)
    }
    fn point(
        &mut self,
        patch: &Patch,
        polygon: &regular::Polygon,
        p: &Point,
        cell: usize,
        corner: usize,
        w: &mut Work<'_>,
    ) -> Result<u32, Error> {
        let (key, encoded, rzv) = patch.owned(p, polygon, self.source.thread.spec().starts, w)?;
        self.register_point(key, encoded, rzv, Origin::Parameter { cell, corner }, w)
    }
    fn register_point(
        &mut self,
        key: Key,
        encoded: Box<[u8]>,
        rzv: [I; 3],
        origin: Origin,
        w: &mut Work<'_>,
    ) -> Result<u32, Error> {
        self.check(encoded.len().checked_mul(2).ok_or(Error::Limit)?)?;
        let reserved = self.resident() + self.curve_owner_bytes() + encoded.len();
        let proposed = self.vertices.len() as u32;
        let copy = encoded.clone();
        let id = self
            .point_owners
            .as_mut()
            .expect("construction point table")
            .vertex_reserved(key, encoded, proposed, reserved, w)?;
        if id == proposed {
            w.step()?;
            self.payload += copy.len();
            self.vertices.push(BoundaryPoint {
                exact: copy,
                source_key: match origin {
                    Origin::Tail { source, .. } => self.vertices[source as usize].source_key,
                    _ => key,
                },
                rzv,
                origin,
                boundary: false,
                left_hand: self.source.thread.spec().left_hand,
            });
        }
        self.check(0)?;
        Ok(id)
    }
    fn edge(
        &mut self,
        patch: &Patch,
        polygon: &regular::Polygon,
        endpoints: [u32; 2],
        points: [&Point; 2],
        face: u32,
        w: &mut Work<'_>,
    ) -> Result<BoundaryCoedge, Error> {
        w.step()?;
        let shared = points[0].active & points[1].active;
        let station = (4..6)
            .find(|&i| shared & (1 << i) != 0)
            .map(|i| polygon.stations[i - 4]);
        let values = [
            coordinates(patch, polygon.kind, points[0], w)?,
            coordinates(patch, polygon.kind, points[1], w)?,
        ];
        let forward = (endpoints[0] <= endpoints[1]) == self.source.thread.spec().left_hand;
        self.register(
            endpoints,
            values,
            station,
            Use {
                face,
                forward: self.source.thread.spec().left_hand,
            },
            w,
        )
        .map(|curve| BoundaryCoedge {
            curve: curve as usize,
            forward,
        })
    }
    fn register(
        &mut self,
        endpoints: [u32; 2],
        mut values: [[Q; 3]; 2],
        station: Option<u8>,
        side: Use,
        w: &mut Work<'_>,
    ) -> Result<u32, Error> {
        let one = Q::integer(1, w)?;
        let mut wraps = [
            values[0][2].cmp(&one, w)? == Ordering::Equal,
            values[1][2].cmp(&one, w)? == Ordering::Equal,
        ];
        // A whole edge on the angular cut has the SAME exact physical law at
        // v0 and v1. For other edges retain the actual one-turn interval.
        if wraps == [true; 2] {
            values[0][2] = Q::integer(0, w)?;
            values[1][2] = Q::integer(0, w)?;
            wraps = [false; 2];
        }
        let forward = endpoints[0] <= endpoints[1];
        let ids = if forward {
            endpoints
        } else {
            values.swap(0, 1);
            wraps.swap(0, 1);
            [endpoints[1], endpoints[0]]
        };
        let delta = [
            values[1][0].sub(&values[0][0], w)?,
            values[1][1].sub(&values[0][1], w)?,
            values[1][2].sub(&values[0][2], w)?,
        ];
        let encoded = Q::encode(
            &[
                &values[0][0],
                &values[0][1],
                &values[0][2],
                &delta[0],
                &delta[1],
                &delta[2],
            ],
            w,
        )?;
        self.check(encoded.len().checked_mul(2).ok_or(Error::Limit)?)?;
        let reserved = self.resident() + self.point_owner_bytes() + encoded.len();
        let proposed = self.curves.len() as u32;
        let copy = encoded.clone();
        let index = self
            .curve_owners
            .as_mut()
            .expect("construction curve table")
            .vertex_reserved(
                Key([
                    Constraint::Rail(i64::from(ids[0])),
                    Constraint::Rail(i64::from(ids[1])),
                ]),
                encoded,
                proposed,
                reserved,
                w,
            )?;
        let side = Use {
            face: side.face,
            forward: if forward { side.forward } else { !side.forward },
        };
        if index == proposed {
            let base = [
                values[0][0].enclosure(w)?.1,
                values[0][1].enclosure(w)?.1,
                values[0][2].enclosure(w)?.1,
            ];
            let bounds = [
                delta[0].enclosure(w)?.1,
                delta[1].enclosure(w)?.1,
                delta[2].enclosure(w)?.1,
            ];
            self.payload += copy.len();
            self.curves.push(BoundaryCurve {
                exact: copy,
                endpoints: ids,
                base,
                delta: bounds,
                v_direction: delta[2].sign(),
                wrap_jump: i8::from(wraps[1]) - i8::from(wraps[0]),
                station,
                uses: [Some(side), None],
            });
        } else {
            let curve = &mut self.curves[index as usize];
            if curve.uses[1].is_some()
                || curve.uses[0].unwrap().forward == side.forward
                || curve.station != station
            {
                return Err(Error::Unresolved);
            }
            curve.uses[1] = Some(side);
        }
        self.check(0)?;
        Ok(index)
    }
    #[cfg(test)]
    fn raw(
        source: Cavity,
        limits: BoundaryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        limits.validate()?;
        let mut work = Work {
            remaining: limits.work,
            used: 0,
            interrupt,
        };
        Self::raw_in(source, limits, &mut work)
    }
    fn raw_in(source: Cavity, limits: BoundaryOptions, work: &mut Work<'_>) -> Result<Self, Error> {
        work.step()?;
        if source
            .entry
            .is_some_and(|e| e.radius <= source.thread.maximum_radius_mm())
        {
            return Err(Error::Unresolved);
        }
        let starts = source.thread.spec().starts;
        let periods = source.tail.map_or(source.full_length_mm(), |t| t.station)
            / source.thread.spec().pitch_mm;
        // AuthoredThread construction qualifies its immutable legacy Spec with
        // length/pitch + 2*starts <=128. The retained prefix is no longer than
        // that source, so this cast and the <=137*5*2 count product are bounded.
        let count = periods.ceil() as usize + starts as usize + 1;
        if count * 5 * (1 + usize::from(source.entry.is_some())) > limits.cells {
            return Err(Error::Limit);
        }
        let mut graph = Self::new(source, limits, work)?;
        for period in -i64::from(starts)..=periods.ceil() as i64 {
            for band in 0..5 {
                let patch = Patch::new(&source, period, band, 1, work)?;
                for kind in [Kind::Cone, Kind::Wall] {
                    if matches!(kind, Kind::Cone) && source.entry.is_none() {
                        continue;
                    }
                    let polygon = patch.cell(0, 0, 1, kind, work)?;
                    if polygon.points.is_empty() {
                        continue;
                    }
                    let index = graph.cells.len();
                    let face = face(
                        &source.thread,
                        Node {
                            band,
                            start: period.rem_euclid(i64::from(starts)) as u32,
                        },
                        kind,
                    );
                    let mut cell = BoundaryCell {
                        period,
                        band,
                        kind,
                        face,
                        vertices: [0; 6],
                        count: polygon.points.len(),
                        parameters: [[I::point(0.); 2]; 6],
                        source: cell_source(&source, &patch, kind, work)?,
                        coedges: [BoundaryCoedge {
                            curve: 0,
                            forward: false,
                        }; 6],
                    };
                    for (corner, point) in polygon.points.iter().enumerate() {
                        cell.vertices[corner] =
                            graph.point(&patch, &polygon, point, index, corner, work)?;
                        // Retain the original source chart BEFORE ownership's
                        // true-periodic canonicalization. Queries need no later
                        // decoding or recomputation of the exact q/v transfer.
                        cell.parameters[corner] =
                            [point.q.enclosure(work)?.1, point.v.enclosure(work)?.1];
                    }
                    for corner in 0..cell.count {
                        let next = (corner + 1) % cell.count;
                        cell.coedges[corner] = graph.edge(
                            &patch,
                            &polygon,
                            [cell.vertices[corner], cell.vertices[next]],
                            [&polygon.points[corner], &polygon.points[next]],
                            face,
                            work,
                        )?;
                    }
                    graph.cells.push(cell);
                }
            }
        }
        graph.check(0)?;
        work.step()?;
        Ok(graph)
    }
}
fn cell_source(
    body: &Cavity,
    patch: &Patch,
    kind: Kind,
    w: &mut Work<'_>,
) -> Result<BoundaryCellSource, Error> {
    let band = body.thread.bands()[patch.band];
    let q0 = patch.q0.enclosure(w)?.1;
    // The exact shared period*p cancels from q1-q0. Retain the original band
    // endpoints instead of subtracting independently rounded global positions.
    let dq = I::point(band.phase_mm[1]) - I::point(band.phase_mm[0]);
    let lead = patch.lead.enclosure(w)?.1; // exact p*starts, not rounded lead_mm
    let radial = match kind {
        Kind::Wall => [
            I::point(band.radius_mm[0]),
            if band.radius_mm[0] == band.radius_mm[1] {
                I::point(0.)
            } else {
                I::point(band.radius_mm[1]) - I::point(band.radius_mm[0])
            },
            I::point(0.),
        ],
        Kind::Cone => {
            let e = body.entry.expect("cone source requires actual entry");
            let k = I::new(e.slope[0], e.slope[1]);
            [I::point(e.radius) - k * q0, -k * dq, -k * lead]
        }
    };
    Ok(BoundaryCellSource {
        q: [q0, dq],
        radial,
        axial: [q0, dq, lead],
    })
}

fn coordinates(patch: &Patch, kind: Kind, p: &Point, w: &mut Work<'_>) -> Result<[Q; 3], Error> {
    let z = p.q.add(&patch.lead.mul(&p.v, w)?, w)?;
    let r = match kind {
        Kind::Wall => patch.slope.mul(&p.q, w)?.add(&patch.intercept, w)?,
        Kind::Cone => {
            let (a, k) = patch.entry.as_ref().ok_or(Error::Invalid)?;
            a.sub(&k.mul(&z, w)?, w)?
        }
    };
    Ok([r, z, p.v.clone()])
}

pub use api::{Boundary, BoundaryCoedge};
pub use polyline::BoundaryPolyline;
