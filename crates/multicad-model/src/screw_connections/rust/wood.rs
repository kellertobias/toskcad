//! Source-owned intentional wood bite, independently from through-nut mating.
use super::*;
use multicad_core::model::{ConnectionReceiverExtent, ResolvedCustomConnectionHardware as Custom};
use multicad_kernel::rust_kernel::wood_thread_contact::Requirements;

pub(super) fn check(
    input: &connection_receivers::Inputs<'_>,
    source: &model::ScrewConnectionRequest,
    evaluated: &Model,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<ScrewConnectionCheck, ModelError> {
    let id = source.connection;
    let fail = |message: &str| error(Some(id), message);
    let unsupported = || {
        fail(
            "Rust direct wood contact currently requires one actual generated wood screw, one blind pilot, and no clamp or washer stack",
        )
    };
    let [row] = source.hardware.as_slice() else {
        return Err(unsupported());
    };
    let [receiver] = source.receivers.as_slice() else {
        return Err(unsupported());
    };
    if source.mode != Mode::DirectWood || row.role != Role::Screw {
        return Err(unsupported());
    }
    let model::ScrewConnectionHardwareSource::Generated { part } = &row.source else {
        return Err(unsupported());
    };
    let h = hardware(input.document, input.parameters, evaluated, row, id)?;
    let ResolvedMechanicalPart::CustomConnection { parameters } = &h.part else {
        return Err(unsupported());
    };
    let Custom::WoodScrew { .. } = parameters.as_ref() else {
        return Err(unsupported());
    };
    let path = model::OccurrencePath {
        root: row.instance,
        members: vec![],
    };
    if evaluated
        .snapshot()
        .leaf_occurrences
        .get(&path)
        .is_none_or(|leaf| leaf.producer != part.feature)
    {
        return Err(fail(
            "The wood screw must retain its actual generated source",
        ));
    }
    let physical = evaluated
        .snapshot()
        .leaf_occurrences
        .get(&receiver.source.occurrence)
        .ok_or_else(|| fail("Pilot must be a current physical leaf"))?;
    if physical.producer != receiver.feature {
        return Err(fail("Pilot must retain the current drilled source"));
    }
    let grip = model::resolve_screw_connection_grip(
        input.document,
        input.parameters,
        source,
        input.revision,
        input.snapshots,
    )
    .map_err(|e| error(Some(id), e))?;
    let resolved = receiver
        .source
        .resolve_current(
            input.document,
            input.parameters,
            input.revision,
            input.snapshots,
        )
        .map_err(|e| error(Some(id), e))?;
    if grip != 0.
        || resolved.station_mm != 0.
        || resolved.thread.is_some()
        || resolved.seat.is_some()
        || !matches!(
            receiver.source.extent,
            ConnectionReceiverExtent::Blind { .. }
        )
    {
        return Err(unsupported());
    }
    let floor = resolved.floor_mm.ok_or_else(unsupported)?;
    let axis = source
        .axis
        .source
        .resolve(
            input.document,
            input.parameters,
            &source.axis.face,
            input.revision,
            input.snapshots,
        )
        .map_err(|e| error(Some(id), e))?
        .frame(ReferenceFrame {
            origin_mm: [0.; 3],
            u: [1., 0., 0.],
            v: [0., 1., 0.],
        });
    at(
        axis,
        placed_frame(&h.part, h.pose),
        0.,
        input.policy.linear_tolerance_mm,
        id,
    )?;
    let minimum = mm(&source.parameters.minimum_engagement, input.parameters, id)?;
    let clearance = mm(&source.parameters.blind_tip_clearance, input.parameters, id)?;
    let recipient = evaluated
        .leaf_solid(&receiver.source.occurrence)
        .ok_or_else(|| fail("The current pilot has no actual material source"))?;
    checkpoint(interrupt)?;
    let contact = recipient
        .rust_direct_wood_contact(
            h.solid,
            Requirements {
                minimum_wall_mm: resolved.minimum_wall_mm,
                minimum_floor_mm: floor,
                minimum_engagement_mm: minimum,
                tip_clearance_mm: clearance,
                linear_tolerance_mm: input.policy.linear_tolerance_mm,
                minimum_feature_mm: input.policy.minimum_feature_mm,
            },
            interrupt,
        )
        .map_err(|e| error(Some(id), e))?;
    let mut queries = 0;
    for other_path in evaluated.snapshot().leaf_occurrences.keys() {
        if *other_path == path || *other_path == receiver.source.occurrence {
            continue;
        }
        checkpoint(interrupt)?;
        let other = evaluated
            .leaf_solid(other_path)
            .ok_or_else(|| fail("Physical component has no actual solid"))?;
        if queries >= 4096 {
            return Err(fail(
                "Wood collision validation exceeds the 4096 candidate budget",
            ));
        }
        queries += 1;
        if h.solid
            .rust_wood_separated_from(other, interrupt)
            .map_err(|e| error(Some(id), e))?
            == Some(true)
        {
            continue;
        }

        if h.solid
            .rust_wood_obstruction_witness(other, interrupt)
            .map_err(|e| error(Some(id), e))?
            .is_some()
        {
            return Err(fail(
                "Actual wood hardware has unintended component interference",
            ));
        }
        return Err(fail(
            "Actual source separation from this component could not be certified; no nominal bounding-box or Boolean fallback is admitted",
        ));
    }
    checkpoint(interrupt)?;
    let engagement = contact.engagement_mm[0].midpoint(contact.engagement_mm[1]);
    Ok(ScrewConnectionCheck{connection:id,mode:source.mode,hardware_instances:vec![row.instance],
        receiver_features:vec![receiver.feature],grip_mm:grip,minimum_engagement_mm:minimum,
        blind_tip_clearance_mm:clearance,actual_engagement_mm:Some(engagement),actual_wood_engagement_mm:Some(engagement),
        native_component_collision_checked:true,
        scope:"Actual retained flat-root wood screw and blind pilot: complete stock wall/floor, exact shared bearing plane, positive open material bite and all other component collision checks; no clamp/washer, installation or strength certification".into()})
}
