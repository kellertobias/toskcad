//! Retained modeled heads share the canonical thread and physical shank.
use super::*;
use multicad_kernel::rust_kernel::{head_profile::Head, profiled_screw::Spec};
mod hexalobular;

pub(super) fn build(
    feature: EntityId,
    underhead_length_mm: f64,
    shank_radius_mm: f64,
    head: &ResolvedScrewHead,
    drive: &ResolvedScrewDrive,
    thread: &ResolvedThread,
    policy: GeometryPolicy,
) -> Result<Option<Solid>, ModelError> {
    if matches!(drive, ResolvedScrewDrive::Hexalobular { .. }) {
        return hexalobular::build(
            feature,
            underhead_length_mm,
            shank_radius_mm,
            head,
            drive,
            thread,
            policy,
        );
    }
    let ResolvedScrewDrive::HexSocket {
        across_flats_mm,
        depth_mm,
    } = *drive
    else {
        return Ok(None);
    };
    let head = match *head {
        ResolvedScrewHead::SocketCap {
            radius_mm,
            height_mm,
        } if depth_mm >= height_mm => Head::Cylindrical {
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
    let (profile, root_radius, frame) = crate::threads::rust_profile(feature, thread)?;
    let solid = Solid::rust_profiled_machine_screw_with_shank(
        profile,
        root_radius,
        shank_radius_mm,
        Spec {
            head,
            shank_length_mm: underhead_length_mm - thread.length_mm,
            socket_across_flats_mm: across_flats_mm,
            socket_depth_mm: depth_mm,
            minimum_wall_mm: policy.minimum_feature_mm,
        },
        feature,
        frame,
        policy,
    )
    .map_err(|cause| error(Some(feature), cause))?;
    if !solid.supports_step_export() {
        return Err(error(
            Some(feature),
            "This profiled screw configuration still requires qualified shoulder and head/thread union STEP exchange",
        ));
    }
    Ok(Some(solid))
}
