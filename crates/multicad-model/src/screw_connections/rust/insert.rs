//! Original generated insert, real metric mating, and one retained blind pilot.
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
            "Rust insert connections require one actual generated insert, one generated metric screw, one blind pilot, and no clamp or washer stack",
        )
    };
    let [receiver] = source.receivers.as_slice() else {
        return Err(unsupported());
    };
    if source.mode != Mode::WoodInsert
        || source.hardware.len() != 2
        || [Role::Screw, Role::Insert]
            .into_iter()
            .any(|role| source.hardware.iter().filter(|r| r.role == role).count() != 1)
    {
        return Err(unsupported());
    }
    let mut placed = Vec::with_capacity(2);
    for row in &source.hardware {
        checkpoint(interrupt)?;
        let ScrewConnectionHardwareSource::Generated { part } = &row.source else {
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
                "Insert hardware must retain its original generated physical source",
            ));
        }
        placed.push(hardware(
            input.document,
            input.parameters,
            evaluated,
            row,
            id,
        )?);
    }
    let screw = placed
        .iter()
        .find(|row| row.role == Role::Screw)
        .expect("one screw checked");
    let insert = placed
        .iter()
        .find(|row| row.role == Role::Insert)
        .expect("one insert checked");
    if !matches!(screw.part, ResolvedMechanicalPart::MachineScrew { .. })
        || !matches!(&insert.part, ResolvedMechanicalPart::CustomConnection { parameters } if matches!(parameters.as_ref(),Custom::WoodInsert{..}))
    {
        return Err(unsupported());
    }
    let physical = evaluated
        .snapshot()
        .leaf_occurrences
        .get(&receiver.source.occurrence)
        .filter(|leaf| leaf.producer == receiver.feature)
        .ok_or_else(|| fail("Insert pilot must retain the current physical drilled source"))?;
    let recipient = evaluated
        .leaf_solid(&physical.path)
        .ok_or_else(|| fail("Insert pilot has no actual material source"))?;
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
        placed_frame(&insert.part, insert.pose),
        0.,
        input.policy.linear_tolerance_mm,
        id,
    )?;
    let minimum = mm(&source.parameters.minimum_engagement, input.parameters, id)?;
    let clearance = mm(&source.parameters.blind_tip_clearance, input.parameters, id)?;
    let contact = recipient
        .rust_wood_insert_contact(
            insert.solid,
            Requirements {
                minimum_wall_mm: resolved.minimum_wall_mm,
                minimum_floor_mm: resolved.floor_mm.ok_or_else(unsupported)?,
                minimum_engagement_mm: minimum,
                tip_clearance_mm: clearance,
                linear_tolerance_mm: input.policy.linear_tolerance_mm,
                minimum_feature_mm: input.policy.minimum_feature_mm,
            },
            interrupt,
        )
        .map_err(|e| error(Some(id), e))?;
    recipient
        .rust_wood_insert_screw_clearance(screw.solid, clearance, interrupt)
        .map_err(|e| error(Some(id), e))?;
    let mating = screw
        .solid
        .rust_wood_insert_mating(insert.solid, minimum, interrupt)
        .map_err(|e| error(Some(id), e))?;
    let mut queries = 0;
    for path in evaluated.snapshot().leaf_occurrences.keys() {
        if *path == receiver.source.occurrence || placed.iter().any(|h| path.root == h.instance) {
            continue;
        }
        checkpoint(interrupt)?;
        let other = evaluated
            .leaf_solid(path)
            .ok_or_else(|| fail("Physical component has no actual solid"))?;
        for row in &placed {
            collision(row.solid, other, evaluated, id, &mut queries)?;
        }
    }
    checkpoint(interrupt)?;
    Ok(ScrewConnectionCheck{
        connection:id,mode:source.mode,hardware_instances:placed.iter().map(|h|h.instance).collect(),
        receiver_features:vec![receiver.feature],grip_mm:grip,minimum_engagement_mm:minimum,
        blind_tip_clearance_mm:clearance,
        actual_engagement_mm:Some(mating.engagement_mm[0].midpoint(mating.engagement_mm[1])),
        actual_wood_engagement_mm:Some(contact.engagement_mm[0].midpoint(contact.engagement_mm[1])),
        native_component_collision_checked:true,
        scope:"Actual generated insert and metric screw: retained finite thread mating, blind pilot wall/floor and positive wood contact, whole screw separation from receiving stock, and all unrelated component collision checks; no clamp/washer, strength or installation certification".into(),
    })
}
