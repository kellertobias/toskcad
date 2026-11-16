//! Actual owned degree-eleven material through the native affine occurrence.
use super::*;
mod planes;
mod properties;
mod topology;
use multicad_rust_kernel::{
    head_profile::RepairedHeadError as H,
    phillips_drive::{Classification, Error as C},
    step_read::{
        PhillipsDegree11MaterialQueryError as E,
        PhillipsDegree11MaterialQueryOptions as SourceOptions,
        PhillipsOwnedDegree11Material as Owned,
    },
    tapered_thread::Error as T,
};
const INVERSE_WORK: usize = 1024;
const SOURCE_WORK: usize = 80_016_384;
const QUERY_WORK: usize = SOURCE_WORK + INVERSE_WORK + 1;
#[derive(Clone, Copy, Debug)]
pub struct Degree11RepairedPhillipsQueryOptions {
    pub maximum_work: usize,
    pub maximum_bytes: usize,
    pub maximum_head_cells: usize,
    pub maximum_thread_cells: usize,
    pub maximum_thread_depth: usize,
}
impl Default for Degree11RepairedPhillipsQueryOptions {
    fn default() -> Self {
        Self {
            maximum_work: QUERY_WORK,
            maximum_bytes: BYTES,
            maximum_head_cells: 4096,
            maximum_thread_cells: 4096,
            maximum_thread_depth: 48,
        }
    }
}
pub struct Degree11RepairedPhillipsQuery<'s> {
    source: &'s Owned,
    world: [[f64; 2]; 3],
    local: [[f64; 2]; 3],
    occurrence: PlaneFrame,
    location: Classification,
    used: usize,
    remaining: usize,
    cells: usize,
    peak: usize,
}
impl Degree11RepairedPhillipsQuery<'_> {
    pub fn source(&self) -> &Owned {
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
    pub fn visited_thread_cells(&self) -> usize {
        self.cells
    }
    /// Logical complete reachable source/cache plus conservative scratch.
    pub fn accounted_peak_bytes(&self) -> usize {
        self.peak
    }
}
fn query_error(e: E) -> GeometryError {
    let kind = match e {
        E::Cancelled
        | E::Head(H::Cancelled)
        | E::Thread(T::Cancelled)
        | E::Cavity(C::Cancelled) => GeometryErrorKind::Cancelled,
        E::Limit | E::Head(H::Limit) | E::Thread(T::Limit) | E::Cavity(C::Limit) => {
            GeometryErrorKind::Limit
        }
        E::Invalid => GeometryErrorKind::Input,
        _ => GeometryErrorKind::Unsupported,
    };
    failure(kind, format!("Degree-eleven material query refused: {e:?}"))
}
impl RepairedPhillipsBody {
    pub(crate) fn has_degree11_source(&self) -> bool {
        self.source.degree11().is_some()
    }
    pub(super) fn classify_degree11(
        &self,
        frame: PlaneFrame,
        world: [[f64; 2]; 3],
        options: Degree11RepairedPhillipsQueryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Degree11RepairedPhillipsQuery<'_>, GeometryError> {
        checkpoint(interrupt)?;
        let source = self
            .source
            .degree11()
            .ok_or_else(|| unsupported("Solid has no degree-eleven repaired material"))?;
        if options.maximum_work > QUERY_WORK
            || options.maximum_bytes > BYTES
            || world
                .iter()
                .any(|b| !b[0].is_finite() || !b[1].is_finite() || b[0] > b[1])
        {
            return Err(input("Invalid degree-eleven material box query"));
        }
        let mut remaining = options
            .maximum_work
            .checked_sub(INVERSE_WORK)
            .ok_or_else(query::exhausted)?;
        // Both the source chart and actual occurrence are inverted as whole
        // boxes. No nominal centre replaces an uncertain world selection.
        let local = self.placement(frame)?.inverse_bounds(world)?;
        let outside = self
            .resident_bytes()
            .checked_sub(source.retained_budget_bytes())
            .ok_or_else(query::exhausted)?;
        let available_bytes = options
            .maximum_bytes
            .checked_sub(outside)
            .ok_or_else(query::exhausted)?;
        let q = source
            .classify_local_box(
                local,
                SourceOptions {
                    maximum_work: remaining.min(SOURCE_WORK),
                    maximum_bytes: available_bytes,
                    maximum_head_cells: options.maximum_head_cells,
                    maximum_thread_cells: options.maximum_thread_cells,
                    maximum_thread_depth: options.maximum_thread_depth,
                },
                interrupt,
            )
            .map_err(query_error)?;
        remaining = remaining
            .checked_sub(q.work_used())
            .ok_or_else(query::exhausted)?;
        let peak = outside
            .checked_add(q.accounted_peak_bytes())
            .ok_or_else(query::exhausted)?;
        if peak > options.maximum_bytes {
            return Err(query::exhausted());
        }
        checkpoint(interrupt)?;
        remaining = remaining.checked_sub(1).ok_or_else(query::exhausted)?;
        Ok(Degree11RepairedPhillipsQuery {
            source,
            world,
            local,
            occurrence: frame,
            location: q.location(),
            used: options.maximum_work - remaining,
            remaining,
            cells: q.visited_thread_cells(),
            peak,
        })
    }
}
impl Solid {
    pub fn rust_degree11_repaired_phillips_material(
        source: Rc<Owned>,
        feature: EntityId,
        occurrence: PlaneFrame,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, GeometryError> {
        checkpoint(interrupt)?;
        policy.validate()?;
        if source
            .containment_clearances_mm()
            .into_iter()
            .chain([source.join_radial_margin_mm()])
            .any(|v| v < policy.minimum_feature_mm)
        {
            return Err(unsupported(
                "Actual cavity and join clearances do not establish the policy minimum feature",
            ));
        }
        let body =
            RepairedPhillipsBody::new_source(Source::Degree11(source), feature, policy, interrupt)?;
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
    pub fn rust_degree11_repaired_phillips_source(&self) -> Option<(&Owned, PlaneFrame)> {
        self.repaired_phillips_body()
            .and_then(|(body, frame)| body.source.degree11().map(|s| (s, frame)))
    }
    pub fn rust_degree11_repaired_phillips_classify_world_box(
        &self,
        bounds: [[f64; 2]; 3],
        options: Degree11RepairedPhillipsQueryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Degree11RepairedPhillipsQuery<'_>, GeometryError> {
        let (body, frame) = self
            .repaired_phillips_body()
            .ok_or_else(|| unsupported("Solid has no degree-eleven repaired material"))?;
        body.classify_degree11(frame, bounds, options, interrupt)
    }
}
