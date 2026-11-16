//! Previously unsupported combined cuts use the authored flat-band material.
use super::*;
#[cfg(test)]
mod tests;
use multicad_kernel::{
    GeometryErrorKind, InternalThreadTreatmentCut,
    rust_kernel::authored_internal_treatment as authored,
};
pub(super) fn build(
    feature: EntityId,
    thread: &ResolvedThread,
    extent: multicad_kernel::rust_kernel::Rational,
    target: &Solid,
    policy: GeometryPolicy,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Solid, ModelError> {
    let (profile, root, frame) = rust_profile_geometry(feature, thread)?;
    if root.is_some() {
        return Err(error(
            Some(feature),
            "Combined internal treatments require the actual flat-band source; rounded roots are not supported",
        ));
    }
    let ends = thread
        .end_treatments
        .as_ref()
        .expect("combined treatment route");
    target
        .rust_cut_authored_internal_thread_treatments_with_exact_extent_and_interrupt(
            InternalThreadTreatmentCut {
                profile,
                intent: authored::Intent {
                    entry_depth_mm: Some(ends.start_chamfer_radial_depth_mm),
                    terminal: Some(authored::Terminal {
                        width_mm: ends.end_relief_length_mm,
                        radial_depth_mm: ends.end_relief_radial_depth_mm,
                    }),
                },
                minimum_wall_mm: thread.minimum_wall_mm,
                frame,
            },
            extent,
            feature,
            policy,
            interrupt,
        )
        .map_err(|cause| {
            if cause.kind == GeometryErrorKind::Cancelled {
                error(Some(feature), "Model evaluation superseded")
            } else {
                error(Some(feature), cause)
            }
        })
}
