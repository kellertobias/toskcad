//! Independent ownership of the actual qualified periodic repair source.
//! No parser arena, cached-rate thread replacement or borrowed proof survives.
//! This does not qualify a manufacturing approximation or raw STEP material.
use super::*;
use crate::head_profile::OwnedRepairedHead;
use std::sync::Arc;
#[cfg(test)]
mod tests;
const SCRATCH: usize = 16 * 1024;
#[derive(Clone, Copy)]
pub struct PeriodicOwnedOptions {
    pub maximum_work: usize,
    pub maximum_total_work: usize,
    pub maximum_bytes: usize,
}
impl Default for PeriodicOwnedOptions {
    fn default() -> Self {
        Self {
            maximum_work: WORK,
            maximum_total_work: TOTAL,
            maximum_bytes: BYTES,
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct PeriodicCurveOwner {
    pub face_entity: u32,
    pub edge_entity: u32,
    pub curve_entity: u32,
    pub vertex_entities: [u32; 2],
    pub maximum_distance_mm: f64,
    pub maximum_cell_derivative_mm: f64,
    pub endpoint_incidence_reserve_mm: f64,
}
impl PeriodicCurveOwner {
    fn from_curve(c: &PeriodicFiniteCurve) -> Self {
        Self {
            face_entity: c.face_entity,
            edge_entity: c.edge_entity,
            curve_entity: c.curve_entity,
            vertex_entities: c.vertex_entities,
            maximum_distance_mm: c.distance_mm,
            maximum_cell_derivative_mm: c.cell_derivative_mm,
            endpoint_incidence_reserve_mm: c.endpoint_incidence_reserve_mm,
        }
    }
}
pub struct OwnedPeriodicMaterial {
    raw: Arc<str>,
    head: OwnedRepairedHead,
    thread: PeriodicThreadSource,
    cavity: Cutter,
    chart: CommonMouthChart,
    determinant: [f64; 2],
    stations: [f64; 2],
    bottom: [f64; 2],
    containment: [f64; 2],
    owners: [PeriodicCurveOwner; 17],
    thread_transfers: [AnalyticThreadTransfer; 4],
    cavity_faces: [CavityFace; 88],
    cavity_count: usize,
    local_volume: [f64; 2],
    components: [[f64; 2]; 3],
    enclosure: [[f64; 2]; 3],
    inherited: usize,
    new_work: usize,
}
impl OwnedPeriodicMaterial {
    pub fn original_input(&self) -> &str {
        &self.raw
    }
    pub fn head(&self) -> &OwnedRepairedHead {
        &self.head
    }
    pub fn thread(&self) -> &PeriodicThreadSource {
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
    pub fn head_bottom_station_mm(&self) -> [f64; 2] {
        self.bottom
    }
    pub fn containment_clearances_mm(&self) -> [f64; 2] {
        self.containment
    }
    pub fn original_curve_owners(&self) -> &[PeriodicCurveOwner; 17] {
        &self.owners
    }
    pub fn thread_transfers(&self) -> &[AnalyticThreadTransfer; 4] {
        &self.thread_transfers
    }
    pub fn local_volume_mm3(&self) -> [f64; 2] {
        self.local_volume
    }
    pub fn component_volumes_mm3(&self) -> [[f64; 2]; 3] {
        self.components
    }
    pub fn world_volume_mm3(&self) -> [f64; 2] {
        let v = interval(self.local_volume) * interval(self.determinant);
        [v.lo, v.hi]
    }
    pub fn local_enclosure_bounds(&self) -> [[f64; 2]; 3] {
        self.enclosure
    }
    pub fn new_work_used(&self) -> usize {
        self.new_work
    }
    pub fn inherited_work_used(&self) -> usize {
        self.inherited
    }
    pub fn total_work_used(&self) -> usize {
        self.inherited + self.new_work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        std::mem::size_of::<Self>()
            .saturating_add(self.raw.len())
            .saturating_add(64)
    }
    pub fn face_count(&self) -> usize {
        self.head.view().segment_count() + 4 + 3 + self.cavity_count
    }
    /// The three shoulder patches retain one original planar face identity.
    pub fn face_source(&self, mut index: usize) -> Option<MaterialFaceSource> {
        let h = self.head.view().segment_count();
        if index < h {
            return Some(MaterialFaceSource::HeadMeridian { segment: index });
        }
        index -= h;
        if index < 4 {
            return Some(MaterialFaceSource::ThreadSide { face: index });
        }
        index -= 4;
        match index {
            0 => Some(MaterialFaceSource::ThreadTip),
            1 => Some(MaterialFaceSource::JoinAnnulus),
            2 => Some(MaterialFaceSource::TopAnnulus),
            _ => self.cavity_faces[..self.cavity_count]
                .get(index - 3)
                .copied()
                .map(|face| MaterialFaceSource::ReversedCavity { face }),
        }
    }
}
impl PeriodicMaterialVolume<'_, '_, '_> {
    pub fn owned_peak_budget_bytes(&self) -> usize {
        self.retained_budget_bytes()
            .saturating_add(SCRATCH)
            .saturating_add(std::mem::size_of::<OwnedPeriodicMaterial>())
            .saturating_add(
                self.source
                    .join
                    .finite_curves()
                    .thread()
                    .body()
                    .literal_input()
                    .len(),
            )
            .saturating_add(64)
    }
    pub fn owned_source(
        &self,
        o: PeriodicOwnedOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> FResult<OwnedPeriodicMaterial> {
        let mut work = Work {
            remaining: o.maximum_work,
            interrupt,
            cancelled: Cell::new(false),
        };
        let result = (|| {
            work.step(1)?;
            if o.maximum_work > WORK || o.maximum_total_work > TOTAL || o.maximum_bytes > BYTES {
                return Err(PhillipsReadError::Invalid);
            }
            let allowance = o.maximum_work.min(
                o.maximum_total_work
                    .checked_sub(self.total_work_used())
                    .ok_or(PhillipsReadError::Limit)?,
            );
            work.remaining = allowance.checked_sub(1).ok_or(PhillipsReadError::Limit)?;
            if self.owned_peak_budget_bytes() > o.maximum_bytes {
                return Err(PhillipsReadError::Limit);
            }
            let source = self.source;
            let join = source.join;
            let body = join.finite_curves().thread().body();
            let raw = body.literal_input();
            let payload = std::mem::size_of::<OwnedPeriodicMaterial>() + raw.len() + 64;
            // All geometry is fixed inline arrays/caches. Only Arc<str> allocates; its
            // exact slice storage and conservative header/alignment are preflighted.
            work.step(payload)?;
            let radius = join.periodic_source().radius_nodes_mm()[0];
            let stations = join.finite_curves().cap_stations_mm();
            let mut enclosure = [[-radius, radius], [-radius, radius], stations];
            for part in 0..join.head_source().segment_count() {
                work.step(257)?;
                let field = join.head_source().certified_interval_bounds(part, [0., 1.]);
                let r = field.position_mm[0][1];
                for axis in &mut enclosure[..2] {
                    axis[0] = axis[0].min(-r);
                    axis[1] = axis[1].max(r);
                }
                enclosure[2][0] = enclosure[2][0].min(field.position_mm[1][0]);
                enclosure[2][1] = enclosure[2][1].max(field.position_mm[1][1]);
            }
            let owners = std::array::from_fn(|i| {
                PeriodicCurveOwner::from_curve(if i == 16 {
                    join.original_complement()
                } else {
                    &join.finite_curves().curves()[i]
                })
            });
            let thread_transfers =
                std::array::from_fn(|i| join.finite_curves().thread().transfers()[i]);
            let mut owned = OwnedPeriodicMaterial {
                raw: Arc::<str>::from(raw),
                head: join.head_source().owned_snapshot(),
                thread: *join.periodic_source(),
                cavity: source.containment.source().clone(),
                chart: source.chart,
                determinant: source.determinant,
                stations,
                bottom: join.head_bottom_station_mm(),
                containment: [
                    source.containment.minimum_radial_clearance_mm(),
                    source.containment.axial_floor_clearance_mm(),
                ],
                owners,
                thread_transfers,
                cavity_faces: source.cavity_faces,
                cavity_count: source.count,
                local_volume: self.local,
                components: self.components,
                enclosure,
                inherited: self.total_work_used(),
                new_work: 0,
            };
            // Arc slices have exact length, and all other payloads are inline. Reread
            // the real retained footprint before publication while both sources live.
            assert!(
                self.retained_budget_bytes()
                    .saturating_add(SCRATCH)
                    .saturating_add(owned.retained_budget_bytes())
                    <= o.maximum_bytes,
                "preflighted inline payload and exact-length Arc slice"
            );
            work.step(1)?;
            owned.new_work = allowance - work.remaining;
            Ok(owned)
        })();
        if work.cancelled.get() {
            Err(PhillipsReadError::Cancelled)
        } else {
            result
        }
    }
}
