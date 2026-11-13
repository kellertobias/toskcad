//! Circle-aware actual-source shell with exact rational q/v ownership.
use super::{
    clip::{self, Line, Point},
    exact::Q,
    jet::Jet,
    owners::{Constraint, Key, Owners},
    *,
};
use crate::authored_metric_thread::{AuthoredMetricThread, Knot, phase_coefficients};
use std::{cmp::Ordering, mem::size_of};
mod patch;
use patch::Patch;
#[cfg(test)]
mod tests;
const BYTES: usize = 64 * 1024 * 1024;
// Plan/cell keeps <128 bounded Q values live (<=8192-bit operands), including
// clipping lines and encoded staging. Reserve for library arithmetic scratch
// and transient containers separately from retained compact ownership bytes.
const SCRATCH: usize = 8 * 1024 * 1024;
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MetricRimEdge {
    pub vertex_indices: [usize; 2],
    /// Unsigned source v enclosures; these do not prove equality or ordering.
    pub angular_turns: [[f64; 2]; 2],
    pub source_period: i64,
    pub source_band: usize,
}
pub struct MetricShell {
    source: AuthoredMetricThread,
    mesh: DisplayMesh,
    rims: [Vec<MetricRimEdge>; 2],
    trace: trace::Trace,
    chord: f64,
    normal: f64,
    work: usize,
}
impl MetricShell {
    pub fn source(&self) -> &AuthoredMetricThread {
        &self.source
    }
    pub fn mesh(&self) -> &DisplayMesh {
        &self.mesh
    }
    pub fn source_chord_bound_mm(&self) -> f64 {
        self.chord
    }
    pub fn source_normal_angle_bound_rad(&self) -> f64 {
        self.normal
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_bytes(&self) -> usize {
        size_of::<Self>()
            + self.mesh.vertices.capacity() * size_of::<f64>()
            + self.mesh.indices.capacity() * size_of::<u32>()
            + self.mesh.face_indices.capacity() * size_of::<u32>()
            + self
                .rims
                .iter()
                .map(|r| r.capacity() * size_of::<MetricRimEdge>())
                .sum::<usize>()
            + self.trace.memory_bytes()
    }
    pub fn rim_count(&self, end: ShellEnd) -> usize {
        self.rims[match end {
            ShellEnd::Start => 0,
            ShellEnd::End => 1,
        }]
        .len()
    }
    pub fn rim_edge(&self, end: ShellEnd, index: usize) -> Option<MetricRimEdge> {
        self.rims[match end {
            ShellEnd::Start => 0,
            ShellEnd::End => 1,
        }]
        .get(index)
        .copied()
    }
    /// Exact physical turn in[0,1), tied to an actual retained rim vertex.
    /// The returned rational is ephemeral; it is not an independently owned
    /// geometry descriptor or a tolerance-based point-equality test.
    pub fn exact_rim_turns(
        &self,
        vertex: usize,
        maximum_work: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<(num_rational::BigRational, usize), Error> {
        let mut w = Work {
            remaining: maximum_work.min(16_000_000),
            used: 0,
            interrupt,
        };
        w.step()?;
        let original = self.trace.parameter(vertex, &mut w)?;
        let zero = Q::integer(0, &mut w)?;
        let one = Q::integer(1, &mut w)?;
        let mut physical = if original.cmp(&one, &mut w)? == Ordering::Less {
            original
        } else {
            original.sub(&one, &mut w)?
        };
        if self.source.spec().left_hand && physical.cmp(&zero, &mut w)? != Ordering::Equal {
            physical = one.sub(&physical, &mut w)?;
        }
        let value = physical.rational(&mut w)?;
        w.step()?;
        Ok((value, w.used))
    }
}
#[derive(Clone, Copy)]
struct Grid {
    counts: [usize; 5],
    rows: u32,
    odd: bool,
    shifted: bool,
}
impl Grid {
    fn new(source: &AuthoredMetricThread, options: Tessellation) -> Result<Self, Error> {
        let root = (2. * std::f64::consts::PI / 3. / options.angle).ceil() as usize;
        let angle = (options.angle * 0.5)
            .min((2. * options.deflection / source.maximum_radius_mm()).sqrt());
        let rows = (std::f64::consts::TAU / angle).ceil();
        if !rows.is_finite() || rows > f64::from(sampling::MAX_ROWS) {
            return Err(Error::Limit);
        }
        Ok(Self {
            counts: [1, 1, root.next_power_of_two().max(1), 1, 1],
            rows: (rows as u32).next_power_of_two().max(8),
            odd: false,
            shifted: false,
        })
    }
    fn effective_rows(self) -> u32 {
        self.rows + if self.odd { 3 } else { 0 }
    }
    fn total(self) -> usize {
        self.counts.into_iter().sum()
    }
    fn advance(&mut self, reason: Reason) -> Result<(), Error> {
        match reason {
            Reason::Repartition if !self.odd && self.rows <= sampling::MAX_ROWS - 3 => {
                self.odd = true
            }
            Reason::Repartition if !self.shifted && self.rows < sampling::MAX_ROWS => {
                self.shifted = true;
                self.odd = false;
            }
            Reason::Refine(band) => {
                if self.counts[band] >= sampling::MAX_SUBDIVISIONS {
                    return Err(Error::Limit);
                }
                self.counts[band] *= 2;
            }
            _ => {
                if self.rows >= sampling::MAX_ROWS {
                    return Err(Error::Limit);
                }
                self.rows *= 2;
                self.odd = false;
            }
        }
        Ok(())
    }
}
enum Reason {
    Refine(usize),
    Repartition,
}
pub(crate) fn metric_shell(
    source: &AuthoredMetricThread,
    options: Tessellation,
    maximum_work: usize,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<MetricShell, Error> {
    options.check().map_err(|_| Error::Invalid)?;
    let mut w = Work {
        remaining: maximum_work.min(16_000_000),
        used: 0,
        interrupt,
    };
    w.step()?;
    let mut grid = Grid::new(source, options)?;
    loop {
        let attempted = attempt(source, options, grid, &mut w);
        match attempted? {
            Ok(mut s) => {
                s.work = w.used;
                w.step()?;
                s.work = w.used;
                return Ok(s);
            }
            Err(reason) => grid.advance(reason)?,
        }
    }
}
fn source_error(e: crate::authored_metric_thread::Error) -> Error {
    match e {
        crate::authored_metric_thread::Error::Invalid => Error::Invalid,
        crate::authored_metric_thread::Error::Limit => Error::Limit,
        crate::authored_metric_thread::Error::Unresolved => Error::Unresolved,
        crate::authored_metric_thread::Error::Cancelled => Error::Cancelled,
    }
}
fn reserve<T>(count: usize) -> Result<Vec<T>, Error> {
    let mut v = Vec::new();
    v.try_reserve_exact(count).map_err(|_| Error::Limit)?;
    Ok(v)
}
fn attempt(
    source: &AuthoredMetricThread,
    options: Tessellation,
    grid: Grid,
    w: &mut Work<'_>,
) -> Result<Result<MetricShell, Reason>, Error> {
    w.step()?;
    let spec = source.spec();
    let rows = grid.effective_rows();
    let total = grid.total();
    if rows > sampling::MAX_ROWS
        || grid
            .counts
            .into_iter()
            .any(|n| n == 0 || n > sampling::MAX_SUBDIVISIONS)
    {
        return Err(Error::Limit);
    }
    let first = -i64::from(spec.starts) - i64::from(grid.shifted);
    let last = (spec.length_mm / spec.pitch_mm).ceil() as i64;
    let periods = (last - first + 1) as usize;
    let rails = periods.checked_mul(total).ok_or(Error::Limit)?;
    let cells = rails.checked_mul(rows as usize).ok_or(Error::Limit)?;
    let maximum = periods
        .checked_mul(total + 1)
        .and_then(|n| n.checked_mul(rows as usize + 1))
        .and_then(|n| n.checked_add(4 * (rails + rows as usize + 1)))
        .ok_or(Error::Limit)?
        .min(options.triangle_limit.saturating_mul(3));
    if maximum == 0 || maximum > 131_072 {
        return Err(Error::Limit);
    }
    let triangles = cells
        .checked_mul(4)
        .ok_or(Error::Limit)?
        .min(options.triangle_limit);
    let rim_limit = 2 * (rails + rows as usize + 1);
    let planned = [
        (maximum * 3, size_of::<f64>()),
        (maximum, size_of::<f64>()),
        (triangles * 3, size_of::<u32>()),
        (triangles, size_of::<u32>()),
        (2 * rim_limit, size_of::<MetricRimEdge>()),
    ];
    let base = SCRATCH + size_of::<MetricShell>();
    let preflight = storage::bytes(&planned, base, BYTES)?;
    let ownership = [
        (Owners::requested_bytes(maximum)?, 1),
        (trace::Trace::requested_bytes(maximum)?, 1),
    ];
    storage::bytes(&ownership, preflight, BYTES)?;
    w.charge(maximum + triangles)?;
    let mut mesh = DisplayMesh {
        vertices: reserve(maximum * 3)?,
        indices: reserve(triangles * 3)?,
        face_indices: reserve(triangles)?,
    };
    let mut rims = [reserve(rim_limit)?, reserve(rim_limit)?];
    let mut point_errors: Vec<f64> = reserve(maximum)?;
    let retained = [
        (mesh.vertices.capacity(), size_of::<f64>()),
        (point_errors.capacity(), size_of::<f64>()),
        (mesh.indices.capacity(), size_of::<u32>()),
        (mesh.face_indices.capacity(), size_of::<u32>()),
        (
            rims[0].capacity() + rims[1].capacity(),
            size_of::<MetricRimEdge>(),
        ),
    ];
    let other = storage::bytes(&retained, base, BYTES)?;
    let origin = if grid.shifted {
        sampling::TURN_OFFSET
    } else {
        0.
    };
    let mut trace = trace::Trace::at_origin(maximum, origin, w)?;
    let mut owners = Owners::new_reserved(maximum, other + trace.memory_bytes(), BYTES, w)?;
    let mut chord = 0_f64;
    let mut normal = 0_f64;
    for period in first..=last {
        for band in 0..5 {
            let patch = Patch::new(source, period, band, grid, w)?;
            for index in 0..grid.counts[band] {
                let strip = patch.strip(index, w)?;
                for row in 0..rows {
                    let polygon = patch.cell(&strip, row, w)?;
                    if polygon.len() < 3 {
                        continue;
                    }
                    let field = patch.field(&polygon, &strip, w)?;
                    let mut ids = Vec::with_capacity(6);
                    let mut errors = Vec::with_capacity(6);
                    for p in &polygon {
                        let (key, token) = patch.identity(p, &strip, row, w)?;
                        let proposed = (mesh.vertices.len() / 3) as u32;
                        let vertex = owners.vertex_reserved(
                            key,
                            token,
                            proposed,
                            other + trace.memory_bytes(),
                            w,
                        )? as usize;
                        if vertex == proposed as usize {
                            let point = patch.point(p, &strip, w)?;
                            room(mesh.vertices.len(), 3, mesh.vertices.capacity())?;
                            room(point_errors.len(), 1, point_errors.capacity())?;
                            let (literal, error) = literal_point(point, options.deflection)?;
                            mesh.vertices.extend(literal);
                            point_errors.push(error);
                        }
                        // Exact key+canonical coordinate equality was checked by
                        // Owners before this cache access. The immutable continuous
                        // source determines one point, so its EARNED literal error
                        // remains valid for every incident facet at fixed quality.
                        w.step()?;
                        let error = point_errors[vertex];
                        ids.push(vertex);
                        errors.push(error);
                    }
                    for (cap, rim) in rims.iter_mut().enumerate() {
                        let positions: Vec<_> = polygon
                            .iter()
                            .enumerate()
                            .filter_map(|(i, p)| (p.active & (1 << (4 + cap)) != 0).then_some(i))
                            .collect();
                        if positions.len() < 2 {
                            continue;
                        }
                        let [a, b] = [positions[0], positions[1]];
                        for pos in [a, b] {
                            let occupied = owners.memory_bytes() + other;
                            trace.remember(ids[pos], &polygon[pos].v, occupied, w)?;
                        }
                        room(rim.len(), 1, rim.capacity())?;
                        rim.push(MetricRimEdge {
                            vertex_indices: [ids[a], ids[b]],
                            angular_turns: [
                                pair(polygon[a].v.enclosure(w)?.1),
                                pair(polygon[b].v.enclosure(w)?.1),
                            ],
                            source_period: period,
                            source_band: band,
                        });
                    }
                    let points: Vec<[f64; 3]> = ids
                        .iter()
                        .map(|&id| mesh.vertices[id * 3..id * 3 + 3].try_into().unwrap())
                        .collect();
                    let face = patch.face();
                    let proofs = match qualify_fan(&field, &points, &errors, options, w)? {
                        Attempt::Complete(proofs) => proofs,
                        Attempt::Refine => return Ok(Err(Reason::Refine(band))),
                        Attempt::Repartition => return Ok(Err(Reason::Repartition)),
                    };
                    for proof in proofs.into_iter().flatten() {
                        room(mesh.face_indices.len(), 1, mesh.face_indices.capacity())?;
                        chord = chord.max(proof.chord);
                        normal = normal.max(proof.angle);
                        append_qualified(&mut mesh, &ids, proof, face, options)?;
                    }
                }
            }
        }
    }
    Ok(Ok(MetricShell {
        source: *source,
        mesh,
        rims,
        trace,
        chord,
        normal,
        work: 0,
    }))
}
fn pair(v: I) -> [f64; 2] {
    [v.lo, v.hi]
}

/// A convex exact source polygon has at most six vertices. Every cyclic fan
/// covers that SAME polygon once; only complete whole-source facet proofs may
/// select a fan. Failed proposals never publish a partial polygon.
fn qualify_fan(
    field: &source::Field,
    points: &[[f64; 3]],
    errors: &[f64],
    options: Tessellation,
    w: &mut Work<'_>,
) -> Result<Attempt<[Option<Qualified>; 4]>, Error> {
    let mut repartition = false;
    for anchor in 0..points.len() {
        w.step()?;
        let mut proofs = [None; 4];
        let mut complete = true;
        for offset in 1..points.len() - 1 {
            let ids = [
                anchor,
                (anchor + offset) % points.len(),
                (anchor + offset + 1) % points.len(),
            ];
            match qualify_attempt(field, points, errors, ids, options, w)? {
                Attempt::Complete(proof) => proofs[offset - 1] = Some(proof),
                failure => {
                    repartition |= matches!(failure, Attempt::Repartition);
                    complete = false;
                    break;
                }
            }
        }
        if complete {
            return Ok(Attempt::Complete(proofs));
        }
    }
    Ok(if repartition {
        Attempt::Repartition
    } else {
        Attempt::Refine
    })
}

fn room(length: usize, additional: usize, capacity: usize) -> Result<(), Error> {
    if additional > capacity.saturating_sub(length) {
        Err(Error::Limit)
    } else {
        Ok(())
    }
}
fn literal_point(point: [I; 3], deflection: f64) -> Result<([f64; 3], f64), Error> {
    let literal = point.map(I::midpoint);
    let error = point
        .into_iter()
        .zip(literal)
        .fold(I::point(0.), |sum, (p, l)| {
            sum + I::point((p - I::point(l)).magnitude())
        })
        .hi;
    if !error.is_finite() || error > deflection / 4. {
        Err(Error::Unresolved)
    } else {
        Ok((literal, error))
    }
}
