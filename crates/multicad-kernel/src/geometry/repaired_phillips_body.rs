//! Native ownership of an explicitly repaired, actual-data-derived material.
//! Raw STEP provenance is retained, but is not exact raw-material authority.
use super::internal_thread_treatment_body::interval::Interval as I;
use super::*;
use multicad_rust_kernel::step_read::PhillipsOwnedDerivedMaterial as OwnedDerivedMaterial;
use std::cell::RefCell;
mod affine;
mod mesh;
mod workspace;
use affine::Placement;
pub use mesh::{RepairedPhillipsMeshBudget, RepairedPhillipsMeshReceipt};
type Point = [f64; 3];
const OPERATION: &str = "repaired Phillips material";
const WORK: usize = 512_000_000;
const BYTES: usize = 64 * 1024 * 1024;
#[cfg(test)]
mod tests;

pub(crate) struct RepairedPhillipsBody {
    source: Rc<OwnedDerivedMaterial>,
    pub(super) feature: EntityId,
    policy: GeometryPolicy,
    local_bounds: Bounds,
    display: RefCell<Option<mesh::Cached>>,
    query_depth: std::cell::Cell<usize>,
    #[cfg(test)]
    builds: std::cell::Cell<usize>,
}
fn failure(kind: GeometryErrorKind, message: impl Into<String>) -> GeometryError {
    error(OPERATION, kind, message)
}
fn unsupported(message: impl Into<String>) -> GeometryError {
    failure(GeometryErrorKind::Unsupported, message)
}
fn checkpoint(interrupt: Option<&dyn Fn() -> bool>) -> Result<(), GeometryError> {
    if interrupt.is_some_and(|f| f()) {
        Err(failure(
            GeometryErrorKind::Cancelled,
            "Repaired source query cancelled",
        ))
    } else {
        Ok(())
    }
}
fn source_failure(cause: multicad_rust_kernel::step_read::PhillipsReadError) -> GeometryError {
    use multicad_rust_kernel::step_read::PhillipsReadError as E;
    failure(
        match cause {
            E::Cancelled => GeometryErrorKind::Cancelled,
            E::Limit => GeometryErrorKind::Limit,
            _ => GeometryErrorKind::Unsupported,
        },
        format!("Repaired source query refused: {cause:?}"),
    )
}
fn check_source_storage(
    source: &OwnedDerivedMaterial,
    maximum: usize,
) -> Result<(), GeometryError> {
    let needed = source.retained_budget_bytes().saturating_add(
        std::mem::size_of::<RepairedPhillipsBody>() + 8 * std::mem::size_of::<usize>(),
    );
    if needed > maximum {
        return Err(failure(
            GeometryErrorKind::Limit,
            "Repaired source ownership exceeds native storage",
        ));
    }
    Ok(())
}
impl RepairedPhillipsBody {
    fn new(
        source: Rc<OwnedDerivedMaterial>,
        feature: EntityId,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, GeometryError> {
        checkpoint(interrupt)?;
        policy.validate()?;
        if feature.is_nil() {
            return Err(input(
                "Repaired material needs its original feature identity",
            ));
        }
        if policy.boolean_fuzzy_mm != 0. {
            return Err(unsupported(
                "Repaired material requires its original feature and no Boolean fuzz",
            ));
        }
        check_source_storage(&source, BYTES)?;
        // The sealed source earned this from actual typed thread cells/caps
        // and actual repaired head fields. Its cavity only removes material.
        let envelope = source.local_enclosure_bounds();
        let bounds = Bounds {
            min: envelope.map(|b| b[0]),
            max: envelope.map(|b| b[1]),
        };
        let result = Self {
            source,
            feature,
            policy,
            local_bounds: bounds,
            display: RefCell::new(None),
            query_depth: std::cell::Cell::new(0),
            #[cfg(test)]
            builds: std::cell::Cell::new(0),
        };
        result.placement(PlaneFrame::default())?;
        checkpoint(interrupt)?;
        Ok(result)
    }
    fn placement(&self, occurrence: PlaneFrame) -> Result<Placement, GeometryError> {
        let placed = Placement::new(self.source.chart(), occurrence)?;
        occurrence
            .checked_metric_bounds(placed.constructed_bounds(self.local_bounds)?, self.policy)?;
        let bounds = placed.bounds(self.local_bounds)?;
        self.policy.point(bounds.min)?;
        self.policy.point(bounds.max)?;
        Ok(placed)
    }
    pub(super) fn validate_frame(&self, frame: PlaneFrame) -> Result<(), GeometryError> {
        self.placement(frame).map(|_| ())
    }
    pub(super) fn bounds(&self, frame: PlaneFrame) -> Result<Bounds, GeometryError> {
        self.placement(frame)?.bounds(self.local_bounds)
    }
    pub(super) fn volume(&self, frame: PlaneFrame) -> Result<[f64; 2], GeometryError> {
        self.placement(frame)?
            .volume_interval(self.source.local_volume_mm3())
    }
}
impl Solid {
    /// Own the qualified repaired geometry. This does not admit the original
    /// raw STEP shell, face properties, finite topology or physical STEP output.
    /// Policy minimum feature controls authoring/placement requests; this
    /// retained-source handoff does not earn a new global separation reserve
    /// for an arbitrary stricter policy. Original source proofs stay intact.
    pub fn rust_repaired_phillips_material(
        source: Rc<OwnedDerivedMaterial>,
        feature: EntityId,
        occurrence: PlaneFrame,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, GeometryError> {
        let body = RepairedPhillipsBody::new(source, feature, policy, interrupt)?;
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
    pub(crate) fn repaired_phillips_body(&self) -> Option<(&RepairedPhillipsBody, PlaneFrame)> {
        match &self.body {
            SolidBody::RepairedPhillips(body, frame) => Some((body, *frame)),
            _ => None,
        }
    }
    pub fn rust_repaired_phillips_source(&self) -> Option<(&OwnedDerivedMaterial, PlaneFrame)> {
        self.repaired_phillips_body()
            .map(|(body, frame)| (body.source.as_ref(), frame))
    }
    /// Interval for the NEW repaired material under both actual affine maps.
    /// It is neither a decoded raw-material scalar nor a literal-mesh volume.
    pub fn rust_repaired_phillips_volume_bounds(&self) -> Result<[f64; 2], GeometryError> {
        let (body, frame) = self
            .repaired_phillips_body()
            .ok_or_else(|| unsupported("Solid has no repaired Phillips source"))?;
        body.volume(frame)
    }
}
