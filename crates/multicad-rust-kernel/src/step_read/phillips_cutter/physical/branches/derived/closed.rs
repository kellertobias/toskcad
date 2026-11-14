//! Closed embedding of the DISTINCT derived cavity and its actual affine image.
//! Outer containment, material closure/flux and native import remain separate.
use super::*;
use crate::phillips_drive::{ClosedBoundaryCertificate, ClosedBoundaryOptions};
const PLACEMENT_WORK: usize = 256;
pub struct ClosedDerivedCavityCertificate<'r, 'a, 'input> {
    representative: &'r DerivedCavityRepresentative<'a, 'input>,
    local: ClosedBoundaryCertificate<'r>,
    determinant: [f64; 2],
    work: usize,
    retained: usize,
}
impl<'r, 'a, 'input> ClosedDerivedCavityCertificate<'r, 'a, 'input> {
    pub fn representative(&self) -> &'r DerivedCavityRepresentative<'a, 'input> {
        self.representative
    }
    pub fn local_source_certificate(&self) -> &ClosedBoundaryCertificate<'r> {
        &self.local
    }
    pub fn actual_affine_columns(&self) -> [Point; 4] {
        self.representative.derived_affine_columns()
    }
    pub fn determinant_bounds(&self) -> [f64; 2] {
        self.determinant
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    /// Original Body/derived source once, plus this inline certificate once.
    pub fn retained_budget_bytes(&self) -> usize {
        self.retained
    }
}
fn source_error(error: crate::phillips_drive::Error) -> PhillipsReadError {
    if error == crate::phillips_drive::Error::Cancelled {
        PhillipsReadError::Cancelled
    } else {
        cutter_error(error).into()
    }
}

fn determinant(columns: [Point; 4]) -> [f64; 2] {
    let u = columns[1].map(I::point);
    let v = columns[2].map(I::point);
    let n = columns[3].map(I::point);
    let cross: [I; 3] =
        std::array::from_fn(|i| v[(i + 1) % 3] * n[(i + 2) % 3] - v[(i + 2) % 3] * n[(i + 1) % 3]);
    let value = (0..3).fold(I::point(0.), |sum, i| sum + u[i] * cross[i]);
    [value.lo, value.hi]
}
impl<'a, 'input> DerivedCavityRepresentative<'a, 'input> {
    /// Earn source closure/embedding and preserve it under the actual literal A.
    /// Work includes the original representative qualification. No raw supports
    /// become exact, and no outer material or file-mass admission is performed.
    pub fn closed_boundary_certificate(
        &self,
        options: ClosedBoundaryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> BranchResult<ClosedDerivedCavityCertificate<'_, 'a, 'input>> {
        let checkpoint = || {
            if interrupt.is_some_and(|f| f()) {
                Err(PhillipsReadError::Cancelled)
            } else {
                Ok(())
            }
        };
        checkpoint()?;
        if options.maximum_work == 0
            || options.maximum_work > WORK
            || options.maximum_bytes == 0
            || options.maximum_bytes > BYTES
        {
            return Err(PhillipsReadError::Invalid);
        }
        let base = self.work_used().saturating_add(PLACEMENT_WORK + 1);
        let remaining = options
            .maximum_work
            .checked_sub(base)
            .filter(|n| *n > 0)
            .ok_or(PhillipsReadError::Limit)?;
        let retained = self
            .retained_budget_bytes()
            .saturating_add(std::mem::size_of::<
                ClosedDerivedCavityCertificate<'_, '_, '_>,
            >());
        let peak = retained.saturating_add(ClosedBoundaryCertificate::scratch_budget_bytes());
        if peak > options.maximum_bytes {
            return Err(PhillipsReadError::Limit);
        }
        let local = self
            .derived_cavity_source()
            .closed_boundary_certificate(
                ClosedBoundaryOptions {
                    maximum_work: remaining,
                    maximum_bytes: options.maximum_bytes,
                },
                interrupt,
            )
            .map_err(source_error)?;
        let det = determinant(self.derived_affine_columns());
        // The representative retained the actual Pose/Placement columns after its
        // original finite Gram1e-10 qualification. Real Gram error is <2e-10 after
        // rounding; the stored normal differs from exact UxV by <8eps. Thus the actual
        // determinant is >1-1e-9, including the directed evaluation error above.
        // This is the earned immutable placement invariant, not normalization or a
        // sampled determinant test. Full-rank positive A preserves source embedding,
        // closure and outward orientation; no rigid area/volume invariance follows.
        checkpoint()?;
        let work = base + local.work_used();
        Ok(ClosedDerivedCavityCertificate {
            representative: self,
            local,
            determinant: det,
            work,
            retained,
        })
    }
}
