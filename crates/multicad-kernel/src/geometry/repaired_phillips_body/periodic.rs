//! Explicit Native ownership and whole-box queries of the actual periodic
//! repair source. Retained polynomial getters and receipts remain distinct.
use super::query::{exhausted, query_failure};
use super::*;
use multicad_rust_kernel::{
    phillips_drive::Classification, step_read::PhillipsOwnedMaterialQueryOptions as SourceOptions,
};
pub(super) mod candidate;
mod cavity_edges;
mod display_edges;
mod properties;
mod topology;
mod topology_planes;
const SOURCE_WORK: usize = 18_008_192;
const INVERSE_WORK: usize = 1024;
const MAX: usize = SOURCE_WORK + INVERSE_WORK + 1;
#[cfg(test)]
mod tests;
pub struct PeriodicRepairedPhillipsQuery<'s> {
    source: &'s OwnedPeriodicMaterial,
    world: [[f64; 2]; 3],
    local: [[f64; 2]; 3],
    occurrence: PlaneFrame,
    location: Classification,
    used: usize,
    remaining: usize,
}
impl PeriodicRepairedPhillipsQuery<'_> {
    pub fn source(&self) -> &OwnedPeriodicMaterial {
        self.source
    }
    pub fn world_bounds_mm(&self) -> [[f64; 2]; 3] {
        self.world
    }
    pub fn local_bounds_mm(&self) -> [[f64; 2]; 3] {
        self.local
    }
    pub fn occurrence(&self) -> PlaneFrame {
        self.occurrence
    }
    pub fn location(&self) -> Classification {
        self.location
    }
    pub fn work_used(&self) -> usize {
        self.used
    }
    pub fn remaining_work(&self) -> usize {
        self.remaining
    }
}
impl RepairedPhillipsBody {
    pub(super) fn classify_periodic(
        &self,
        frame: PlaneFrame,
        world: [[f64; 2]; 3],
        options: RepairedPhillipsQueryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<PeriodicRepairedPhillipsQuery<'_>, GeometryError> {
        checkpoint(interrupt)?;
        let source = self
            .source
            .periodic()
            .ok_or_else(|| unsupported("Solid has no periodic repaired Phillips source"))?;
        if options.maximum_work > MAX
            || options.maximum_head_cells == 0
            || options.maximum_head_cells > 1_000_000
            || world.iter().any(|b| !I::new(b[0], b[1]).valid())
        {
            return Err(input("Invalid periodic repaired material box query"));
        }
        let mut remaining = options
            .maximum_work
            .checked_sub(INVERSE_WORK)
            .ok_or_else(exhausted)?;
        let local = self.placement(frame)?.inverse_bounds(world)?;
        let query = source
            .classify_local_box(
                local,
                SourceOptions {
                    maximum_work: remaining.min(SOURCE_WORK),
                    maximum_head_cells: options.maximum_head_cells,
                },
                interrupt,
            )
            .map_err(query_failure)?;
        remaining -= query.work_used();
        checkpoint(interrupt)?;
        remaining = remaining.checked_sub(1).ok_or_else(exhausted)?;
        Ok(PeriodicRepairedPhillipsQuery {
            source,
            world,
            local,
            occurrence: frame,
            location: query.location(),
            used: options.maximum_work - remaining,
            remaining,
        })
    }
}
impl Solid {
    pub fn rust_periodic_repaired_phillips_material(
        source: Rc<OwnedPeriodicMaterial>,
        feature: EntityId,
        occurrence: PlaneFrame,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, GeometryError> {
        checkpoint(interrupt)?;
        policy.validate()?;
        // Every public caller must establish its requested policy against
        // the actual owned material. An authored seed's older wall request,
        // like a decoded proposal-only margin, is not policy authority.
        if source
            .containment_clearances_mm()
            .into_iter()
            .any(|v| v < policy.minimum_feature_mm)
        {
            return Err(unsupported(
                "Actual cavity clearances do not establish the policy minimum feature",
            ));
        }
        let body =
            RepairedPhillipsBody::new_source(Source::Periodic(source), feature, policy, interrupt)?;
        body.validate_frame(occurrence)?;
        checkpoint(interrupt)?;
        Ok(Self {
            body: SolidBody::RepairedPhillips(Rc::new(body), occurrence),
            policy,
            provenance: vec![],
            rust_origins: vec![],
            retained_placement: None,
        })
    }
    pub fn rust_periodic_repaired_phillips_source(
        &self,
    ) -> Option<(&OwnedPeriodicMaterial, PlaneFrame)> {
        self.repaired_phillips_body()
            .and_then(|(body, frame)| body.source.periodic().map(|source| (source, frame)))
    }
    pub fn rust_periodic_repaired_phillips_classify_world_box(
        &self,
        bounds: [[f64; 2]; 3],
        options: RepairedPhillipsQueryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<PeriodicRepairedPhillipsQuery<'_>, GeometryError> {
        let (body, frame) = self
            .repaired_phillips_body()
            .ok_or_else(|| unsupported("Solid has no periodic repaired Phillips source"))?;
        body.classify_periodic(frame, bounds, options, interrupt)
    }
}
