//! Actual analytic raw-support maps, distinct from finite trim/face filling.
use super::*;

fn step_error(e: crate::step::phillips_cutter::Error) -> PhillipsReadError {
    use crate::step::phillips_cutter::Error;
    match e {
        Error::Cancelled => PhillipsReadError::Cancelled,
        Error::Limit => PhillipsReadError::Limit,
        Error::InvalidGeometry | Error::InvalidPlacement => PhillipsReadError::Invalid,
        Error::Unsupported => PhillipsReadError::Unsupported,
    }
}
const ATLAS_WORK: usize = 1_000_000;
const PATCH_WORK: usize = 16_384;
const EXTRA_BYTES: usize = 2 * 1024 * 1024;

#[derive(Clone, Copy)]
pub struct RawCavitySupportOptions {
    /// Fresh work only, in addition to the retained representative's prefix.
    pub maximum_fresh_work: usize,
    /// Includes the original representative and all live new storage.
    pub maximum_bytes: usize,
}
impl Default for RawCavitySupportOptions {
    fn default() -> Self {
        Self {
            maximum_fresh_work: WORK,
            maximum_bytes: BYTES,
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct RawCavitySupportFace {
    pub source_role: FaceKey,
    pub raw_surface_entity: u32,
    pub maximum_position_mm: f64,
    /// Original atlas unit-square chart derivatives, NOT leaf-scaled or normals.
    pub maximum_chart_derivative_mm: [f64; 2],
    pub source_patches: usize,
}
pub struct RawCavitySupportCorrespondence<'r, 'a, 'input> {
    source: &'r DerivedCavityRepresentative<'a, 'input>,
    faces: [RawCavitySupportFace; 29],
    fresh_work: usize,
    peak: usize,
}
impl<'r, 'a, 'input> RawCavitySupportCorrespondence<'r, 'a, 'input> {
    pub fn source(&self) -> &'r DerivedCavityRepresentative<'a, 'input> {
        self.source
    }
    pub fn faces(&self) -> &[RawCavitySupportFace; 29] {
        &self.faces
    }
    pub fn new_work_used(&self) -> usize {
        self.fresh_work
    }
    pub fn work_used(&self) -> usize {
        self.source.work_used().saturating_add(self.fresh_work)
    }
    pub fn accounted_peak_bytes(&self) -> usize {
        self.peak
    }
}
impl<'a, 'input> DerivedCavityRepresentative<'a, 'input> {
    /// Compare every actual atlas patch with its original raw analytic support.
    /// This earns neither finite raw trim filling nor exact raw occupied mass.
    pub fn raw_support_correspondence(
        &self,
        options: RawCavitySupportOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> BranchResult<RawCavitySupportCorrespondence<'_, 'a, 'input>> {
        let checkpoint = || {
            if interrupt.is_some_and(|f| f()) {
                Err(PhillipsReadError::Cancelled)
            } else {
                Ok(())
            }
        };
        checkpoint()?;
        if options.maximum_fresh_work == 0
            || options.maximum_fresh_work > WORK
            || options.maximum_bytes == 0
            || options.maximum_bytes > BYTES
        {
            return Err(PhillipsReadError::Invalid);
        }
        let peak = self
            .retained_budget_bytes()
            .saturating_add(std::mem::size_of::<
                RawCavitySupportCorrespondence<'_, '_, '_>,
            >())
            .saturating_add(EXTRA_BYTES);
        if peak > options.maximum_bytes {
            return Err(PhillipsReadError::Limit);
        }
        let mut budget = Budget {
            remaining: options.maximum_fresh_work,
            interrupt,
        };
        budget.charge(ATLAS_WORK)?;
        let atlas = self
            .source
            .surface_patches_with_cancel(&|| interrupt.is_some_and(|f| f()))
            .map_err(|e| {
                if e == crate::phillips_drive::Error::Cancelled {
                    PhillipsReadError::Cancelled
                } else {
                    cutter_error(e).into()
                }
            })?;
        // Atlas construction is bounded by the already-qualified fixed rounded
        // family. Include its ACTUAL capacities before any following proof.
        if atlas.retained_budget_bytes() > EXTRA_BYTES / 2 {
            return Err(PhillipsReadError::Limit);
        }
        let shell = Shell::new(self.body.source_entities(), 1e-7, true);
        let pose = Pose {
            origin: self.frame[0],
            axes: [self.frame[1], self.frame[2], self.frame[3]],
        };
        let mut faces = self.faces.map(|f| RawCavitySupportFace {
            source_role: f.source_role,
            raw_surface_entity: f.raw_surface_entity,
            maximum_position_mm: 0.,
            maximum_chart_derivative_mm: [0.; 2],
            source_patches: 0,
        });
        for (index, patch) in atlas.patches().iter().enumerate() {
            budget.charge(PATCH_WORK)?;
            let face = patch.face_index();
            let actual = self.faces.get(face).ok_or(PhillipsReadError::Unsupported)?;
            if atlas.topology().faces()[face].source_identity() != actual.source_role {
                return Err(PhillipsReadError::Unsupported);
            }
            let source =
                support::Support::new(atlas.topology().face_surface(face).map_err(cutter_error)?)
                    .map_err(step_error)?;
            let raw =
                geometry::support(&shell, actual.raw_surface_entity, actual.raw_cavity_forward)?;
            let bounds = atlas
                .patch_bounds(index, [[0., 1.]; 2])
                .map_err(cutter_error)?;
            let source_box = [bounds.point_mm.map(|x| x[0]), bounds.point_mm.map(|x| x[1])];
            let position = source
                .distance_bound(&pose, &raw, source_box)
                .map_err(step_error)?;
            let speed = [bounds.du_mm, bounds.dv_mm].map(|v| support::norm(v.map(support::b)).hi);
            let derivative = source
                .chart_derivative_bounds(&pose, &raw, speed)
                .map_err(step_error)?;
            let result = &mut faces[face];
            result.maximum_position_mm = result.maximum_position_mm.max(position);
            for (old, new) in result
                .maximum_chart_derivative_mm
                .iter_mut()
                .zip(derivative)
            {
                *old = old.max(new);
            }
            result.source_patches += 1;
        }
        if faces.iter().any(|f| f.source_patches == 0) {
            return Err(PhillipsReadError::Unsupported);
        }
        budget.charge(1)?;
        checkpoint()?;
        Ok(RawCavitySupportCorrespondence {
            source: self,
            faces,
            fresh_work: options.maximum_fresh_work - budget.remaining,
            peak,
        })
    }
}
