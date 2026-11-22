//! Immutable owned normal/curvature preparation. The literal arena is shared,
//! never reparsed or reconstructed as an authored thread recipe.
use super::*;
const PREPARATION_WORK: usize = 512_000_000;
#[derive(Clone, Copy, Debug)]
pub struct LiteralGeometryOptions {
    pub maximum_work: usize,
    pub maximum_bytes: usize,
}
impl Default for LiteralGeometryOptions {
    fn default() -> Self {
        Self {
            maximum_work: PREPARATION_WORK,
            maximum_bytes: BYTES,
        }
    }
}
pub struct LiteralGeometryPlan {
    pub(super) source: OwnedMaterial,
    pub(super) normals: Vec<Vec<normal::Normal>>,
    preparation: usize,
    retained: usize,
    construction_peak: usize,
}
impl LiteralGeometryPlan {
    pub(super) fn cell_normal_point(
        &self,
        tensor: usize,
        index: [usize; 2],
        uv: [f64; 2],
        budget: &mut Budget<'_>,
    ) -> Result<[I; 3]> {
        normal::point_values(
            self.original_normal_coefficients(tensor, index)?,
            uv,
            budget,
        )
    }

    /// Outward Sv×Su,22×22 u-major, ORIGINAL unit-root derivatives,
    /// already transported by the actual STEP operator cofactor.
    pub(crate) fn original_normal_coefficients(
        &self,
        tensor: usize,
        index: [usize; 2],
    ) -> Result<&[[I; 3]]> {
        let t = self.source.tensors().get(tensor).ok_or(Error::Invalid)?;
        if index[0] >= t.spans[0] || index[1] >= t.spans[1] {
            return Err(Error::Invalid);
        }
        Ok(self.normals[tensor][index[0] * t.spans[1] + index[1]].coefficients())
    }
    pub fn source(&self) -> &OwnedMaterial {
        &self.source
    }
    pub fn preparation_work_used(&self) -> usize {
        self.preparation
    }
    /// Includes the one shared original literal arena, original normal nets,
    /// fixed descriptors and cache capacities. A mesh receipt already includes
    /// this same arena/cache; their reported peaks must not be added.
    pub fn retained_budget_bytes(&self) -> usize {
        self.retained
    }
    pub fn construction_peak_bytes(&self) -> usize {
        self.construction_peak
    }
}
impl OwnedMaterial {
    pub fn geometry_plan(
        &self,
        options: LiteralGeometryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<LiteralGeometryPlan> {
        if options.maximum_work == 0
            || options.maximum_work > PREPARATION_WORK
            || options.maximum_bytes == 0
            || options.maximum_bytes > BYTES
        {
            return Err(Error::Invalid);
        }
        let bytes = self.retained_budget_bytes() + SCRATCH + size_of::<LiteralGeometryPlan>();
        if bytes > options.maximum_bytes {
            return Err(Error::Limit);
        }
        let mut budget = Budget {
            used: 0,
            remaining: options.maximum_work,
            bytes,
            peak: bytes,
            limit: options.maximum_bytes,
            interrupt,
        };
        budget.step(1)?;
        let normals = normal::prepare(self, &mut budget)?;
        budget.step(1)?;
        Ok(LiteralGeometryPlan {
            source: self.clone(),
            normals,
            preparation: budget.used,
            retained: budget.bytes - SCRATCH,
            construction_peak: budget.peak,
        })
    }
}
