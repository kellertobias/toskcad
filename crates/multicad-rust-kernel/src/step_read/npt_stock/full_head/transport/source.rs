//! Supplied author/source matching follows actual file material admission.
//! World C0 includes the exact retained original factor chain; D1 remains
//! explicitly ORIGINAL raw thread parameters, not a sewn/head derivative.
use super::*;
use crate::{
    step::npt_stock::{CorrespondenceOptions, CorrespondenceReceipt, SurfacePlan},
    step_read::npt_stock::LiteralOccurrencePlacement,
};
#[cfg(test)]
mod tests;
#[derive(Clone, Copy, Debug)]
pub struct MappedSourceReceipt {
    pub actual_file: MappedHeadReceipt,
    /// Original thread law versus RAW retained local STEP coefficients only.
    pub original_raw_thread: CorrespondenceReceipt,
    pub head_author_corner_error_mm: f64,
    pub thread_axial_transfer_displacement_mm: f64,
    pub local_sewn_reference_position_error_mm: f64,
    /// SAME actual sewn local enclosure, original chain versus file operator.
    pub factor_chain_discrepancy_mm: f64,
    pub original_linear_l2_norm_upper: f64,
    pub world_sewn_reference_position_error_mm: f64,
    pub head_and_frame_work: usize,
    pub fresh_thread_correspondence_work: usize,
    pub aggregate_work: usize,
    pub accounted_peak_bytes: usize,
}
pub struct MappedSourceQualification<'actual, 'plan, 'source, 'pose> {
    actual: &'actual MappedActualHead,
    source: &'plan SurfacePlan<'source>,
    original: &'pose LiteralOccurrencePlacement,
    receipt: MappedSourceReceipt,
}
impl<'actual, 'plan, 'source, 'pose> MappedSourceQualification<'actual, 'plan, 'source, 'pose> {
    pub fn actual(&self) -> &'actual MappedActualHead {
        self.actual
    }
    pub fn source_plan(&self) -> &'plan SurfacePlan<'source> {
        self.source
    }
    pub fn original_placement(&self) -> &'pose LiteralOccurrencePlacement {
        self.original
    }
    pub fn receipt(&self) -> MappedSourceReceipt {
        self.receipt
    }
}
impl MappedActualHead {
    pub fn qualify_original_source<'actual, 'plan, 'source, 'pose>(
        &'actual self,
        source: &'plan SurfacePlan<'source>,
        head: HexHeadIntent,
        original: &'pose LiteralOccurrencePlacement,
        options: SourceOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<MappedSourceQualification<'actual, 'plan, 'source, 'pose>> {
        super::super::source::validate_source_options(head, options)?;
        let r = options.resources;
        let owned = &self.local;
        if owned.material_tolerance_mm() > options.material_tolerance_mm
            || owned.minimum_wall_mm() < options.minimum_wall_mm
        {
            return Err(Error::Unsupported);
        }
        if original.determinant_bounds()[0] <= 0. {
            return Err(Error::Unsupported);
        }
        let source_work = source.work_used();
        let pose_work = original.work_used();
        let inherited = self
            .receipt
            .aggregate_work
            .checked_add(source_work)
            .and_then(|n| n.checked_add(pose_work))
            .ok_or(Error::Limit)?;
        let residents = r
            .additional_resident_bytes
            .checked_add(self.receipt.retained_bytes)
            .and_then(|n| n.checked_add(source.retained_budget_bytes()))
            .and_then(|n| n.checked_add(original.retained_budget_bytes()))
            .and_then(|n| n.checked_add(size_of::<MappedSourceQualification<'_, '_, '_, '_>>()))
            .ok_or(Error::Limit)?;
        let initial = residents.checked_add(16 * 1024).ok_or(Error::Limit)?;
        // Previously earned actual-file construction remains an inherited
        // phase peak; current source/frame allocations are live additions.
        let inherited_peak = self
            .receipt
            .accounted_peak_bytes
            .checked_sub(self.receipt.additional_resident_bytes)
            .and_then(|n| n.checked_add(source.retained_budget_bytes()))
            .and_then(|n| n.checked_add(original.retained_budget_bytes()))
            .and_then(|n| n.checked_add(r.additional_resident_bytes))
            .and_then(|n| n.checked_add(size_of::<MappedSourceQualification<'_, '_, '_, '_>>()))
            .ok_or(Error::Limit)?;
        let pose_peak = original
            .accounted_peak_bytes()
            .checked_add(self.receipt.retained_bytes)
            .and_then(|n| n.checked_add(source.retained_budget_bytes()))
            .and_then(|n| n.checked_add(r.additional_resident_bytes))
            .and_then(|n| n.checked_add(size_of::<MappedSourceQualification<'_, '_, '_, '_>>()))
            .ok_or(Error::Limit)?;
        let inherited_peak = inherited_peak.max(pose_peak);
        if initial > r.maximum_bytes || inherited_peak > r.maximum_bytes {
            return Err(Error::Limit);
        }
        let mut b = Budget {
            used: 0,
            remaining: r.maximum_work.checked_sub(inherited).ok_or(Error::Limit)?,
            bytes: initial,
            peak: initial,
            limit: r.maximum_bytes,
            interrupt,
        };
        b.step(1024)?;
        let a = original
            .linear_rows_bounds()
            .map(|r| r.map(|x| I::new(x[0], x[1])));
        let origin = original.translation_bounds_mm().map(|x| I::new(x[0], x[1]));
        let frame = self.actual_frame();
        let scale = I::new(frame.scale[0], frame.scale[1]);
        let bound = owned.local_sewn_enclosure_mm().map(|x| I::new(x[0], x[1]));
        let delta = std::array::from_fn::<_, 3, _>(|i| {
            let translation = origin[i] - I::new(frame.origin[i][0], frame.origin[i][1]);
            (0..3).fold(translation, |sum, j| {
                sum + (a[i][j] - scale * I::new(frame.axes[j][i][0], frame.axes[j][i][1]))
                    * bound[j]
            })
        });
        if delta.iter().any(|x| !x.valid()) {
            return Err(Error::Unsupported);
        }
        let discrepancy = delta
            .into_iter()
            .fold(I::point(0.), |s, x| s + I::point(x.magnitude()))
            .hi;
        let norm = original.linear_l2_norm_upper();
        let local_tolerance =
            ((I::point(options.position_error_mm) - I::point(discrepancy)) / I::point(norm)).lo;
        if !local_tolerance.is_finite() || local_tolerance <= 0. {
            return Err(Error::Unsupported);
        }
        let corner = super::super::source::actual_head_error(
            owned.head(),
            owned.head_polygon_mm(),
            owned.head_stations_mm(),
            source,
            head,
            local_tolerance,
            &mut b,
        )?;
        let head_error = (I::point(corner) + I::point(owned.planar_sewing_mm())).hi;
        let displacement = owned.axial_sewing_transfer().domain_displacement_mm();
        let raw_tolerance = (I::point(local_tolerance) - I::point(displacement)).lo;
        if !raw_tolerance.is_finite() || raw_tolerance <= 0. || head_error > local_tolerance {
            return Err(Error::Unsupported);
        }
        let auxiliary = owned.auxiliary_thread();
        let shared_prefix = source_work
            .checked_add(auxiliary.work_used())
            .ok_or(Error::Limit)?;
        let available = r
            .maximum_work
            .checked_sub(inherited)
            .and_then(|n| n.checked_sub(b.used))
            .ok_or(Error::Limit)?;
        let child_work = shared_prefix
            .checked_add(available)
            .ok_or(Error::Limit)?
            .min(96_000_000);
        // The correspondence child includes SAME source+auxiliary allocations.
        // Remove only those exact retained owners from its additional scope.
        let outside = residents
            .checked_sub(source.retained_budget_bytes())
            .and_then(|n| n.checked_sub(auxiliary.retained_budget_bytes()))
            .ok_or(Error::Limit)?;
        let proof = source
            .correspond_literal(
                auxiliary,
                CorrespondenceOptions {
                    position_error_mm: raw_tolerance,
                    original_derivative_error_mm: options.original_thread_derivative_error_mm,
                    maximum_work: child_work,
                    maximum_bytes: r
                        .maximum_bytes
                        .checked_sub(outside)
                        .filter(|n| *n > 0)
                        .ok_or(Error::Limit)?,
                },
                interrupt,
            )
            .map_err(super::super::source::source_error)?;
        let fresh = proof
            .aggregate_work_used()
            .checked_sub(shared_prefix)
            .ok_or(Error::Limit)?;
        b.step(fresh)?;
        b.step(64)?;
        let local = (I::point(proof.position_error_mm()) + I::point(displacement))
            .hi
            .max(head_error);
        let world = (I::point(norm) * I::point(local) + I::point(discrepancy)).hi;
        if !world.is_finite() || world > options.position_error_mm {
            return Err(Error::Unsupported);
        }
        let aggregate = inherited.checked_add(b.used).ok_or(Error::Limit)?;
        let peak = inherited_peak.max(b.peak).max(
            outside
                .checked_add(proof.accounted_peak_bytes())
                .ok_or(Error::Limit)?,
        );
        if peak > r.maximum_bytes {
            return Err(Error::Limit);
        }
        b.step(0)?;
        Ok(MappedSourceQualification {
            actual: self,
            source,
            original,
            receipt: MappedSourceReceipt {
                actual_file: self.receipt,
                original_raw_thread: proof,
                head_author_corner_error_mm: corner,
                thread_axial_transfer_displacement_mm: displacement,
                local_sewn_reference_position_error_mm: local,
                factor_chain_discrepancy_mm: discrepancy,
                original_linear_l2_norm_upper: norm,
                world_sewn_reference_position_error_mm: world,
                head_and_frame_work: b.used - fresh,
                fresh_thread_correspondence_work: fresh,
                aggregate_work: aggregate,
                accounted_peak_bytes: peak,
            },
        })
    }
}
