//! Immutable source-owned D11 material, extracted only from earned full volume.
//! No parser/head/receipt borrow survives and no manufacturing admission follows.
use super::*;
use crate::head_profile::OwnedRepairedHead;
use std::sync::Arc;
const SCRATCH: usize = 16 * 1024;
const ARC_BYTES: usize = 64;
#[derive(Clone, Copy)]
pub struct Degree11OwnedMaterialOptions {
    /// Full inherited volume/material aggregate plus fresh copy work.
    pub maximum_work: usize,
    pub maximum_bytes: usize,
}
impl Default for Degree11OwnedMaterialOptions {
    fn default() -> Self {
        Self {
            maximum_work: WORK_LIMIT,
            maximum_bytes: MAX_BYTES,
        }
    }
}
pub struct OwnedDegree11Material {
    raw: Arc<str>,
    head: OwnedRepairedHead,
    thread: Vec<Degree11ThreadFace>,
    inventory: Vec<Degree11MaterialFace>,
    rim: Vec<Degree11RimUse>,
    cavity: Cutter,
    chart: CommonMouthChart,
    determinant: [f64; 2],
    stations: [f64; 2],
    collar: Degree11JoinCollar,
    join_margin: f64,
    containment: [f64; 2],
    enclosure: [[f64; 2]; 3],
    core: f64,
    winding: [[f64; 2]; 2],
    local_volume: [f64; 2],
    components: [[f64; 2]; 3],
    work: usize,
    peak: usize,
}
impl OwnedDegree11Material {
    pub fn original_input(&self) -> &str {
        &self.raw
    }
    pub fn head(&self) -> &OwnedRepairedHead {
        &self.head
    }
    pub fn thread_faces(&self) -> &[Degree11ThreadFace] {
        &self.thread
    }
    pub fn cavity(&self) -> &Cutter {
        &self.cavity
    }
    pub fn chart(&self) -> CommonMouthChart {
        self.chart
    }
    pub fn determinant_bounds(&self) -> [f64; 2] {
        self.determinant
    }
    pub fn cap_stations_mm(&self) -> [f64; 2] {
        self.stations
    }
    pub fn transfer_collar(&self) -> Degree11JoinCollar {
        self.collar
    }
    pub fn join_radial_margin_mm(&self) -> f64 {
        self.join_margin
    }
    pub fn containment_clearances_mm(&self) -> [f64; 2] {
        self.containment
    }
    pub fn local_enclosure_bounds_mm(&self) -> [[f64; 2]; 3] {
        self.enclosure
    }
    /// Full occupied THREAD cylinder only between the retained cap stations;
    /// the head cavity is a separate part of this material source.
    pub fn occupied_core_radius_mm(&self) -> f64 {
        self.core
    }
    pub fn cap_winding_angle_bounds_rad(&self) -> [[f64; 2]; 2] {
        self.winding
    }
    pub fn local_volume_mm3(&self) -> [f64; 2] {
        self.local_volume
    }
    pub fn world_volume_mm3(&self) -> [f64; 2] {
        let v = interval(self.local_volume) * interval(self.determinant);
        [v.lo, v.hi]
    }
    pub fn component_volumes_mm3(&self) -> [[f64; 2]; 3] {
        self.components
    }
    pub fn face_count(&self) -> usize {
        self.inventory.len()
    }
    pub fn face_source(&self, index: usize) -> Option<Degree11MaterialFace> {
        self.inventory.get(index).copied()
    }
    pub fn join_inner_uses(&self) -> &[Degree11RimUse] {
        &self.rim
    }
    pub fn entry_rim_poles(
        &self,
        index: usize,
    ) -> Option<impl ExactSizeIterator<Item = Point> + '_> {
        let source = self.thread.get(index)?.source();
        let [width, height] = source.dimensions();
        Some((0..width).map(move |u| source.poles()[u * height]))
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn accounted_peak_bytes(&self) -> usize {
        self.peak
    }
    /// Complete logical live source capacities after all borrowed owners drop.
    /// Shared immutable Arc poles are counted in full, never treated as free.
    pub fn retained_budget_bytes(&self) -> usize {
        std::mem::size_of::<Self>()
            + self.raw.len()
            + ARC_BYTES
            + self.thread.capacity() * std::mem::size_of::<Degree11ThreadFace>()
            + self.inventory.capacity() * std::mem::size_of::<Degree11MaterialFace>()
            + self.rim.capacity() * std::mem::size_of::<Degree11RimUse>()
            + self
                .thread
                .iter()
                .map(|f| {
                    std::mem::size_of_val(f.source.poles())
                        + ARC_BYTES
                        + f.original_controls.capacity() * std::mem::size_of::<u32>()
                })
                .sum::<usize>()
    }
}
fn payload_bound(source: &Degree11ClosedMaterial<'_, '_, '_>) -> usize {
    let faces = source.thread().carrier().faces();
    // Twice all immutable pole payloads, descriptor/ID buffers and metadata
    // conservatively cover clone capacities and old/new construction overlap.
    // Pole Arcs actually share their exact allocations; they are still counted.
    let thread = faces
        .iter()
        .map(|f| {
            std::mem::size_of::<Degree11ThreadFace>()
                + std::mem::size_of_val(f.source().poles())
                + ARC_BYTES
                + std::mem::size_of_val(f.original_control_entities())
        })
        .sum::<usize>();
    std::mem::size_of::<OwnedDegree11Material>()
        + 2 * thread
        + 2 * source.face_count() * std::mem::size_of::<Degree11MaterialFace>()
        + 2 * faces.len() * std::mem::size_of::<Degree11RimUse>()
        + source.body().literal_input().len()
        + ARC_BYTES
}
fn reserve<T>(count: usize) -> std::result::Result<Vec<T>, PhillipsReadError> {
    let mut v = Vec::new();
    v.try_reserve_exact(count)
        .map_err(|_| PhillipsReadError::Limit)?;
    Ok(v)
}
fn check_peak(
    prefix: usize,
    payload: usize,
    maximum: usize,
) -> std::result::Result<usize, PhillipsReadError> {
    let peak = prefix
        .checked_add(payload)
        .ok_or(PhillipsReadError::Limit)?;
    if peak > maximum {
        Err(PhillipsReadError::Limit)
    } else {
        Ok(peak)
    }
}
impl Degree11MaterialVolume<'_, '_, '_, '_> {
    /// Complete old source + earned volume receipt + copy workspace + new
    /// owner capacity envelope, reserved before any new geometry allocation.
    pub fn owned_peak_budget_bytes(&self) -> usize {
        self.source.retained_budget_bytes()
            + std::mem::size_of::<Self>()
            + std::mem::size_of::<OwnedDegree11Material>()
            + SCRATCH
            + payload_bound(self.source)
    }
    /// Copy from this same immutable volume proof. Refusal/cancellation keeps
    /// the original receipt/source usable; numerical prerequisites are not rerun.
    pub fn owned_snapshot(
        &self,
        options: Degree11OwnedMaterialOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> std::result::Result<OwnedDegree11Material, PhillipsReadError> {
        let mut work = Work {
            remaining: options.maximum_work,
            interrupt,
            cancelled: Cell::new(false),
        };
        work.step(1).map_err(|e| bound_error(e, &work))?;
        if options.maximum_work > WORK_LIMIT || options.maximum_bytes > MAX_BYTES {
            return Err(PhillipsReadError::Invalid);
        }
        let source = self.source;
        work.step(self.work_used())
            .and_then(|_| work.step(256 + 16 * source.thread().carrier().faces().len()))
            .map_err(|e| bound_error(e, &work))?;
        let payload = payload_bound(source);
        let prefix = source.retained_budget_bytes()
            + std::mem::size_of::<Self>()
            + std::mem::size_of::<OwnedDegree11Material>()
            + SCRATCH;
        let peak = check_peak(prefix, payload, options.maximum_bytes)?;
        // Byte-linear credit covers actual descriptor copies, immutable raw
        // text, original IDs, source snapshots and all fixed metadata queries.
        work.step(payload).map_err(|e| bound_error(e, &work))?;
        let mut enclosure = source.thread().local_enclosure_bounds_mm();
        for i in 0..source.head().segment_count() {
            work.step(257).map_err(|e| bound_error(e, &work))?;
            let p = source
                .head()
                .certified_interval_bounds(i, [0., 1.])
                .position_mm;
            let radius = p[0][1];
            for axis in enclosure.iter_mut().take(2) {
                axis[0] = axis[0].min(-radius);
                axis[1] = axis[1].max(radius);
            }
            enclosure[2][0] = enclosure[2][0].min(p[1][0]);
            enclosure[2][1] = enclosure[2][1].max(p[1][1]);
        }
        let mut thread = reserve(source.thread().carrier().faces().len())?;
        for f in source.thread().carrier().faces() {
            work.step(1).map_err(|e| bound_error(e, &work))?;
            let mut ids = reserve(f.original_controls.len())?;
            ids.extend_from_slice(&f.original_controls);
            thread.push(Degree11ThreadFace {
                source: f.source.clone(),
                original_controls: ids,
                forward: f.forward,
                transfer: f.transfer,
            });
        }
        let mut inventory = reserve(source.face_count())?;
        for i in 0..source.face_count() {
            inventory.push(
                source
                    .face_source(i)
                    .expect("sealed complete material face inventory"),
            );
        }
        let mut rim = reserve(thread.len())?;
        rim.extend(source.join_inner_uses());
        work.step(1).map_err(|e| bound_error(e, &work))?;
        let result = OwnedDegree11Material {
            raw: Arc::from(source.body().literal_input()),
            head: source.head().owned_snapshot(),
            thread,
            inventory,
            rim,
            cavity: source.cavity().clone(),
            chart: source.chart(),
            determinant: source.determinant_bounds(),
            stations: source.thread().cap_stations_mm(),
            collar: source.transfer_collar(),
            join_margin: source.join_radial_margin_mm(),
            containment: [
                source.containment().minimum_radial_clearance_mm(),
                source.containment().axial_floor_clearance_mm(),
            ],
            enclosure,
            core: source.thread().occupied_core_radius_mm(),
            winding: source.thread().cap_winding_angle_bounds_rad(),
            local_volume: self.local_volume_mm3(),
            components: self.component_volumes_mm3(),
            work: 0,
            peak,
        };
        publish(result, prefix, options, &mut work)
    }
}
// Validate actual post-allocation capacities before ownership escapes. The
// constructor's conservative reservation is distinct from this publication
// check; cancellation drops the complete ready value without exposing it.
fn publish(
    mut result: OwnedDegree11Material,
    prefix: usize,
    options: Degree11OwnedMaterialOptions,
    work: &mut Work<'_>,
) -> std::result::Result<OwnedDegree11Material, PhillipsReadError> {
    check_peak(
        prefix,
        result.retained_budget_bytes(),
        options.maximum_bytes,
    )?;
    work.step(1).map_err(|e| bound_error(e, work))?;
    result.work = options.maximum_work - work.remaining;
    Ok(result)
}
#[cfg(test)]
mod tests;

mod query;
pub use query::{Degree11MaterialQuery, Degree11MaterialQueryError, Degree11MaterialQueryOptions};

mod properties;
pub use properties::{Degree11EndDiskOptions, Degree11EndDiskProperties, Degree11EndDisks};
