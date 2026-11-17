//! Whole-box queries of the same repaired sources through actual affine maps.
use super::*;
use multicad_rust_kernel::{
    head_profile::RepairedHeadError as H,
    phillips_drive::{Classification, Error as C},
    step_read::{
        PhillipsOwnedMaterialQueryError as E, PhillipsOwnedMaterialQueryOptions as SourceOptions,
    },
    tapered_thread::Error as T,
};
const SOURCE_WORK: usize = 18_008_192;
const INVERSE_WORK: usize = 1024;
const MAXIMUM_WORK: usize = SOURCE_WORK + INVERSE_WORK + 1;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RepairedPhillipsQueryOptions {
    pub maximum_work: usize,
    pub maximum_head_cells: usize,
}
impl Default for RepairedPhillipsQueryOptions {
    fn default() -> Self {
        Self {
            maximum_work: MAXIMUM_WORK,
            maximum_head_cells: 4096,
        }
    }
}
/// Same-owner whole-set classification. UnresolvedBoundary includes virtual
/// join disks and is not a certificate that an exposed boundary was found.
pub struct RepairedPhillipsQuery<'s> {
    source: &'s OwnedDerivedMaterial,
    world: [[f64; 2]; 3],
    local: [[f64; 2]; 3],
    occurrence: PlaneFrame,
    location: Classification,
    used: usize,
    remaining: usize,
}
impl RepairedPhillipsQuery<'_> {
    pub fn source(&self) -> &OwnedDerivedMaterial {
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
pub(super) fn query_failure(cause: E) -> GeometryError {
    let kind = match cause {
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
    failure(
        kind,
        format!("Retained repaired material source query refused: {cause:?}"),
    )
}
pub(super) fn exhausted() -> GeometryError {
    failure(
        GeometryErrorKind::Limit,
        "Repaired material query work exhausted",
    )
}
impl RepairedPhillipsBody {
    fn classify(
        &self,
        frame: PlaneFrame,
        world: [[f64; 2]; 3],
        options: RepairedPhillipsQueryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<RepairedPhillipsQuery<'_>, GeometryError> {
        let source = self
            .source
            .retained()
            .ok_or_else(|| unsupported("Use the explicit periodic source query API"))?;
        checkpoint(interrupt)?;
        if options.maximum_work > MAXIMUM_WORK
            || options.maximum_head_cells == 0
            || options.maximum_head_cells > 1_000_000
            || world.iter().any(|b| !I::new(b[0], b[1]).valid())
        {
            return Err(input("Invalid repaired material box query"));
        }
        let mut remaining = options
            .maximum_work
            .checked_sub(INVERSE_WORK)
            .ok_or_else(exhausted)?;
        // These are numerical query enclosures, not new authored coordinates:
        // inverse uncertainty may extend beyond model coordinate input limits.
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
        remaining -= query.work_used(); // Source earned used <= the supplied balance.
        checkpoint(interrupt)?;
        remaining = remaining.checked_sub(1).ok_or_else(exhausted)?;
        Ok(RepairedPhillipsQuery {
            source,
            world,
            local,
            occurrence: frame,
            location: query.location(),
            used: options.maximum_work - remaining,
            remaining,
        })
    }
    pub(in crate::geometry) fn contains(
        &self,
        frame: PlaneFrame,
        point: Point,
    ) -> Result<bool, GeometryError> {
        // Preserve the entire world tolerance box through BOTH inverses;
        // taking its nominal centre could erase a real source crossing.
        let tolerance = self.policy.linear_tolerance_mm;
        let world = point.map(|v| {
            let b = I::point(v) + I::new(-tolerance, tolerance);
            [b.lo, b.hi]
        });
        let location = if self.source.degree11().is_some() {
            self.classify_degree11(frame, world, Default::default(), None)?
                .location()
        } else if self.source.periodic().is_some() {
            self.classify_periodic(frame, world, Default::default(), None)?
                .location()
        } else {
            self.classify(frame, world, Default::default(), None)?
                .location()
        };
        match location {
            Classification::Inside => Ok(true),
            Classification::Outside => Ok(false),
            Classification::UnresolvedBoundary => Err(unsupported(
                "Actual repaired material tolerance box crosses an unresolved surface or virtual join",
            )),
        }
    }
}
impl Solid {
    /// Every point of this finite world box is classified against the same
    /// owned repaired head/thread/cavity, never a nominal recipe or mesh.
    pub fn rust_repaired_phillips_classify_world_box(
        &self,
        bounds: [[f64; 2]; 3],
        options: RepairedPhillipsQueryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<RepairedPhillipsQuery<'_>, GeometryError> {
        let (body, frame) = self
            .repaired_phillips_body()
            .ok_or_else(|| unsupported("Solid has no repaired Phillips source"))?;
        body.classify(frame, bounds, options, interrupt)
    }
}
#[cfg(test)]
mod tests;
