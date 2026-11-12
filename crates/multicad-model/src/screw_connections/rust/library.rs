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
            let distance =
                |a: [f64; 3], b: [f64; 3]| (a[0] - b[0]).hypot(a[1] - b[1]).hypot(a[2] - b[2]);
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
            "Pinned connection hardware requires an ISO metric nut or a flat washer",
        )),
    }
}
