//! Canonical metric/imperial dimensions feed the same owned physical source.
use super::*;
use multicad_kernel::rust_kernel::{
    hexalobular_drive::{Cone, Drive, TangentSpec},
    profiled_screw::SourceSpec,
    screw_recess::Recess,
};

pub(super) fn build(
    feature: EntityId,
    underhead_length_mm: f64,
    shank_radius_mm: f64,
    head: &ResolvedScrewHead,
    drive: &ResolvedScrewDrive,
    thread: &ResolvedThread,
    policy: GeometryPolicy,
) -> Result<Option<Solid>, ModelError> {
    let head = match *head {
        ResolvedScrewHead::SocketCap {
            radius_mm,
            height_mm,
        } => Head::Cylindrical {
            radius_mm,
            height_mm,
        },
        ResolvedScrewHead::Countersunk {
            radius_mm,
            height_mm,
            included_angle_deg,
            edge_chamfer_mm,
            neck_radius_mm,
        } => Head::Countersunk {
            radius_mm,
            height_mm,
            included_angle_deg,
            edge_chamfer_mm,
            neck_radius_mm,
        },
        ResolvedScrewHead::Pan {
            radius_mm,
            height_mm,
            crown_radius_mm,
            top_land_radius_mm,
            neck_radius_mm,
        }
        | ResolvedScrewHead::Button {
            radius_mm,
            height_mm,
            crown_radius_mm,
            top_land_radius_mm,
            neck_radius_mm,
        } => Head::Crowned {
            radius_mm,
            height_mm,
            crown_radius_mm,
            top_land_radius_mm,
            neck_radius_mm,
        },
        _ => return Ok(None),
    };
    let ResolvedScrewDrive::Hexalobular {
        major_radius_mm,
        minor_radius_mm,
        outer_lobe_radius_mm,
        depth_mm,
        bottom_angle_deg,
    } = *drive
    else {
        unreachable!("typed curved-drive dispatcher");
    };
    // Canonical depth measures the straight wall at the major radius. The
    // conical floor continues to its real apex, preserving the original tool.
    let slope = 1. / (bottom_angle_deg / 2.).to_radians().tan();
    let source = Drive::from_authored_tangent(
        TangentSpec {
            major_diameter_mm: 2. * major_radius_mm,
            minor_diameter_mm: 2. * minor_radius_mm,
            outer_lobe_radius_mm,
        },
        Cone {
            tip_depth_mm: depth_mm + major_radius_mm * slope,
            depth_per_radius: slope,
        },
    )
    .map_err(|cause| {
        error(
            Some(feature),
            format!("Invalid physical hexalobular drive: {cause:?}"),
        )
    })?;
    // A geometric tolerance times a squared source extent is a volume budget.
    // Scaling the budget avoids oversolving larger screws to an arbitrary
    // fixed absolute volume precision.
    let extent = source
        .maximum_radius_bound_mm()
        .max(source.cone().tip_depth_mm);
    let width = policy.linear_tolerance_mm / 32. * extent * extent;
    let recess = Recess::conical_hexalobular(source, width, 100_000).map_err(|cause| {
        error(
            Some(feature),
            format!("Unresolved physical hexalobular drive: {cause:?}"),
        )
    })?;
    let (profile, root, frame) = crate::threads::rust_profile(feature, thread)?;
    let solid = Solid::rust_profiled_machine_screw_with_recess(
        profile,
        root,
        shank_radius_mm,
        SourceSpec {
            head: Some(head),
            shank_length_mm: underhead_length_mm - thread.length_mm,
            minimum_wall_mm: policy.minimum_feature_mm,
        },
        recess,
        feature,
        frame,
        policy,
    )
    .map_err(|cause| error(Some(feature), cause))?;
    if !solid.supports_step_export() {
        return Err(error(
            Some(feature),
            "This curved screw configuration requires qualified physical STEP exchange",
        ));
    }
    Ok(Some(solid))
}
