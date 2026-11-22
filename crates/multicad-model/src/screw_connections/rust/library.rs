//! Pinned identity is checked separately from actual material and placement.
use super::*;
use multicad_core::{document::EntityData, model::ThreadFamily, symbols};

pub(super) fn qualify(
    input: &connection_receivers::Inputs<'_>,
    evaluated: &Model,
    row: &model::ScrewConnectionHardwareRequest,
    actual: &Hardware<'_>,
    id: EntityId,
) -> Result<(), ModelError> {
    let fail = |message| error(Some(id), message);
    let ScrewConnectionHardwareSource::Library {
        definition,
        placement,
    } = &row.source
    else {
        unreachable!("library dispatch")
    };
    let pin = symbols::read_instance(input.document, placement.record)
        .map_err(|cause| error(Some(id), cause))?;
    if placement.instance != row.instance
        || pin.definition != placement.definition
        || pin.wrapper_metadata != placement.wrapper_metadata
        || pin.member != placement.member
        || pin.instance != row.instance
        || pin.wrapper != placement.wrapper
        || pin.ids != placement.ids
        || pin.pin != definition.pin().map_err(|cause| error(Some(id), cause))?
    {
        return Err(fail(
            "Pinned hardware ownership differs from its current Library source",
        ));
    }
    let root = placement
        .ids
        .get(&definition.root)
        .and_then(|mapped| input.document.get(*mapped))
        .ok_or_else(|| fail("Pinned hardware lost its mapped root"))?;
    let EntityData::Part { features, .. } = &root.data else {
        return Err(fail("Pinned hardware must own one mechanical Part"));
    };
    let producer = features
        .iter()
        .rev()
        .find(|feature| {
            matches!(
                input.document.get(**feature).map(|entity| &entity.data),
                Some(EntityData::Feature {
                    suppressed: false,
                    ..
                })
            )
        })
        .ok_or_else(|| fail("Pinned hardware has no current mechanical producer"))?;
    if !evaluated
        .snapshot()
        .leaf_occurrences
        .values()
        .any(|leaf| leaf.path.root == row.instance && leaf.producer == *producer)
    {
        return Err(fail(
            "Pinned hardware must retain its actual mapped mechanical output",
        ));
    }
    let mapped = model::read_mechanical_part(input.document, *producer)
        .map_err(|cause| error(Some(id), cause))?;
    let metric_nut = matches!(mapped.parameters, model::MechanicalPartParameters::HexNut {thread,..} if thread.family == ThreadFamily::IsoMetric);
    match (row.role, &actual.part) {
        (Role::HangerBolt, ResolvedMechanicalPart::CustomConnection { parameters }) => {
            let model::ResolvedCustomConnectionHardware::HangerBolt {
                wood_thread,
                metric_thread,
                plain_length_mm,
                plain_radius_mm,
                ..
            } = parameters.as_ref()
            else {
                return Err(fail(
                    "Pinned hanger requires its actual connected hanger source",
                ));
            };
            let data = actual
                .solid
                .rust_hanger_data()
                .ok_or_else(|| fail("Pinned hanger lost its actual connected material"))?;
            let (wood, root, _) = crate::threads::rust_profile(id, wood_thread)?;
            let (metric, metric_root, _) = crate::threads::rust_profile(id, metric_thread)?;
            let expected_wood = multicad_kernel::rust_kernel::helical_thread::Spec {
                major_radius_mm: wood.major_radius_mm,
                pitch_mm: wood.pitch_mm,
                length_mm: wood.length_mm,
                depth_mm: wood.depth_mm,
                crest_width_mm: wood.crest_width_mm,
                lower_flank_deg: wood.lower_flank_deg,
                upper_flank_deg: wood.upper_flank_deg,
                starts: wood.starts,
                left_hand: wood.left_hand,
                root_radius_mm: root,
            };
            let expected_metric = multicad_kernel::rust_kernel::authored_metric_thread::Spec {
                major_radius_mm: metric.major_radius_mm,
                pitch_mm: metric.pitch_mm,
                length_mm: metric.length_mm,
                crest_width_mm: metric.crest_width_mm,
                starts: metric.starts,
                left_hand: metric.left_hand,
            };
            if root.is_some()
                || metric_root.is_none()
                || data.source.wood().spec() != expected_wood
                || data.source.metric().spec() != expected_metric
                || data.source.spec().plain_length_mm != *plain_length_mm
                || data.source.spec().plain_radius_mm != *plain_radius_mm
            {
                return Err(fail(
                    "Actual pinned hanger differs from its current resolved material",
                ));
            }
            // Current mapped producer owns the actual source. Complete stack
            // certification separately checks its raw frame, phase and every
            // contact; the Library pin grants no placement/material exemption.
            Ok(())
        }
        (Role::HeadWasher | Role::NutWasher, ResolvedMechanicalPart::Washer { .. }) => {
            // The common bearing proof checks actual ring dimensions and the
            // complete placed frame; a pin cannot waive that physical check.
            Ok(())
        }
        (
            Role::Nut,
            ResolvedMechanicalPart::HexNut {
                across_flats_mm,
                thread,
            },
        ) if metric_nut => {
            let (nut, frame) = actual
                .solid
                .rust_hex_nut_geometry()
                .ok_or_else(|| fail("Pinned nut has no actual retained hex-nut material"))?;
            let (profile, root, _) = crate::threads::rust_profile(id, thread)?;
            let expected = multicad_kernel::rust_kernel::helical_thread::Spec {
                major_radius_mm: profile.major_radius_mm,
                pitch_mm: profile.pitch_mm,
                length_mm: profile.length_mm,
                depth_mm: profile.depth_mm,
                crest_width_mm: profile.crest_width_mm,
                lower_flank_deg: profile.lower_flank_deg,
                upper_flank_deg: profile.upper_flank_deg,
                starts: profile.starts,
                left_hand: profile.left_hand,
                root_radius_mm: root,
            };
            if nut.across_flats_mm() != *across_flats_mm || nut.thread().spec() != expected {
                return Err(fail(
                    "Actual pinned nut material differs from its current resolved source",
                ));
            }
            let expected = placed_frame(&actual.part, actual.pose);
            let distance = |a: [f64; 3], b: [f64; 3]| {
                (a[0] - b[0])
                    .portable_hypot(a[1] - b[1])
                    .portable_hypot(a[2] - b[2])
            };
            let residual = distance(frame.origin, expected.origin_mm)
                + nut.across_flats_mm() / 3_f64.sqrt()
                    * (distance(frame.u, expected.u) + distance(frame.v, expected.v))
                + nut.thread().spec().length_mm * distance(frame.normal(), expected.normal());
            if !residual.is_finite() || residual > input.policy.linear_tolerance_mm / 8. {
                return Err(fail(
                    "Actual pinned nut placement differs from its source mount",
                ));
            }
            Ok(())
        }
        _ => Err(fail(
            "Pinned connection hardware requires an actual hanger, ISO metric nut or flat washer",
        )),
    }
}
