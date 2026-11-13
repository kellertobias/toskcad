//! Canonical authored Phillips parameters retain the actual finite cutter.
use super::*;
use multicad_core::model::ResolvedPhillipsInnerForm;
use multicad_kernel::rust_kernel::{
    head_profile::Profile,
    helical_thread::{self, Thread},
    phillips_drive::{Cavity, Cutter, Domain, Inner, PhillipsScrew, Spec as DriveSpec},
    profiled_screw::SourceSpec,
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
    let ResolvedScrewDrive::PhillipsH { profile: drive } = drive else {
        return Ok(None);
    };
    let head = match *head {
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
        } => Head::Crowned {
            radius_mm,
            height_mm,
            crown_radius_mm,
            top_land_radius_mm,
            neck_radius_mm,
        },
        _ => return Ok(None),
    };
    let fail = |cause| {
        error(
            Some(feature),
            format!("Unresolved authored Phillips screw: {cause:?}"),
        )
    };
    let (canonical, root_radius_mm, underhead) = crate::threads::rust_profile(feature, thread)?;
    let thread = Thread::new(helical_thread::Spec {
        major_radius_mm: canonical.major_radius_mm,
        pitch_mm: canonical.pitch_mm,
        length_mm: canonical.length_mm,
        depth_mm: canonical.depth_mm,
        crest_width_mm: canonical.crest_width_mm,
        lower_flank_deg: canonical.lower_flank_deg,
        upper_flank_deg: canonical.upper_flank_deg,
        starts: canonical.starts,
        left_hand: canonical.left_hand,
        root_radius_mm,
    })
    .map_err(|cause| {
        error(
            Some(feature),
            format!("Invalid retained Phillips thread: {cause:?}"),
        )
    })?;
    let plain = underhead_length_mm - thread.spec().length_mm;
    let profile = Profile::new_for_thread_union(head, shank_radius_mm, plain).map_err(|cause| {
        error(
            Some(feature),
            format!("Invalid retained Phillips head: {cause:?}"),
        )
    })?;
    let spec = DriveSpec {
        mouth_diameter_mm: drive.mouth_diameter_mm,
        mouth_radius_mm: drive.mouth_radius_mm,
        b_width_mm: drive.b_width_mm,
        g_diameter_mm: drive.g_diameter_mm,
        tilt_angle_deg: drive.tilt_angle_deg,
        outer_flank_angle_deg: drive.outer_flank_angle_deg,
        floor_angle_deg: drive.floor_angle_deg,
        outer_corner_angle_deg: drive.outer_corner_angle_deg,
        inner: match drive.inner {
            ResolvedPhillipsInnerForm::Rounded { radius_mm } => Inner::Rounded { radius_mm },
            ResolvedPhillipsInnerForm::ShallowV {
                angle_deg,
                e_width_mm,
            } => Inner::ShallowV {
                angle_deg,
                e_width_mm,
            },
        },
    };
    let tip = (spec.mouth_diameter_mm - spec.g_diameter_mm)
        / (2. * spec.outer_flank_angle_deg.to_radians().tan())
        + spec.g_diameter_mm / 2. * spec.floor_angle_deg.to_radians().tan();
    let extent = profile
        .support([1., 0., 0.])
        .max(spec.mouth_diameter_mm)
        .max(tip.abs());
    // These finite construction candidates are not a geometry certificate.
    // Cutter/Cavity must independently prove every corner closure lies beyond
    // the real envelope. Inadequate or nonfinite candidates refuse unchanged.
    let cutter = Cutter::new(
        spec,
        Domain {
            corner_start_depth_mm: -2. * spec.mouth_diameter_mm,
            corner_end_depth_mm: tip + spec.mouth_diameter_mm,
            corner_reach_mm: 8. * extent,
        },
    )
    .map_err(fail)?;
    let cavity = Cavity::new(
        cutter,
        policy.linear_tolerance_mm / 32. * extent * extent,
        100_000,
    )
    .map_err(fail)?;
    let source = PhillipsScrew::new(
        thread,
        profile,
        SourceSpec {
            head: Some(head),
            shank_length_mm: plain,
            minimum_wall_mm: policy.minimum_feature_mm,
        },
        cavity,
        100_000,
        &|| false,
    )
    .map_err(fail)?;
    let frame = PlaneFrame {
        origin: std::array::from_fn(|i| underhead.origin[i] + underhead.normal()[i] * plain),
        ..underhead
    };
    // Manufacturing uses the qualified retained source. Physical STEP remains
    // explicitly unavailable until its independent import contract is earned.
    Solid::rust_phillips_screw_from_source(source, feature, frame, policy)
        .map(Some)
        .map_err(|cause| error(Some(feature), cause))
}
