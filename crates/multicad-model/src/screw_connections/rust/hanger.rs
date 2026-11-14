//! Actual connected source, machined stock, mating and all related materials.
use super::*;
use multicad_core::model::{ConnectionReceiverExtent, ResolvedCustomConnectionHardware as Custom};
use multicad_kernel::rust_kernel::wood_thread_contact::{
    HangerMatingRequirements, PassageRequirements, Requirements, StackRequirements,
};

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
            "Hanger connections require one blind pilot, one through bore, an authored hanger, a hex nut and optional flat nut washers",
        )
    };
    if source.mode != Mode::HangerBolt
        || !(2..=64).contains(&source.hardware.len())
        || source
            .hardware
            .iter()
            .any(|row| !matches!(row.role, Role::HangerBolt | Role::Nut | Role::NutWasher))
        || [Role::HangerBolt, Role::Nut].into_iter().any(|role| {
            source
                .hardware
                .iter()
                .filter(|row| row.role == role)
                .count()
                != 1
        })
    {
        return Err(unsupported());
    }
    let [a, b] = source.receivers.as_slice() else {
        return Err(unsupported());
    };
    let (wood, clamp) = match (&a.source.extent, &b.source.extent) {
        (ConnectionReceiverExtent::Blind { .. }, ConnectionReceiverExtent::Through {}) => (a, b),
        (ConnectionReceiverExtent::Through {}, ConnectionReceiverExtent::Blind { .. }) => (b, a),
        _ => return Err(unsupported()),
    };
    let mut placed = Vec::with_capacity(source.hardware.len());
    let mut related = BTreeSet::new();
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
                "Hanger hardware must retain its actual generated physical source",
            ));
        }
        related.insert(path);
        placed.push(hardware(
            input.document,
            input.parameters,
            evaluated,
            row,
            id,
        )?);
    }
    let hanger = placed
        .iter()
        .find(|h| h.role == Role::HangerBolt)
        .expect("one hanger checked");
    let nut = placed
        .iter()
        .find(|h| h.role == Role::Nut)
        .expect("one nut checked");
    if !matches!(&hanger.part,ResolvedMechanicalPart::CustomConnection {parameters}if matches!(parameters.as_ref(),Custom::HangerBolt {..}))
        || !matches!(nut.part, ResolvedMechanicalPart::HexNut { .. })
        || placed
            .iter()
            .filter(|h| h.role == Role::NutWasher)
            .any(|h| !matches!(h.part, ResolvedMechanicalPart::Washer { .. }))
    {
        return Err(unsupported());
    }
    let grip = model::resolve_screw_connection_grip(
        input.document,
        input.parameters,
        source,
        input.revision,
        input.snapshots,
    )
    .map_err(|e| error(Some(id), e))?;
    let minimum = mm(&source.parameters.minimum_engagement, input.parameters, id)?;
    let tip = mm(&source.parameters.blind_tip_clearance, input.parameters, id)?;
    let wood_data = wood
        .source
        .resolve_current(
            input.document,
            input.parameters,
            input.revision,
            input.snapshots,
        )
        .map_err(|e| error(Some(id), e))?;
    let clamp_data = clamp
        .source
        .resolve_current(
            input.document,
            input.parameters,
            input.revision,
            input.snapshots,
        )
        .map_err(|e| error(Some(id), e))?;
    let tolerance = input.policy.linear_tolerance_mm;
    if !grip.is_finite()
        || grip <= 0.
        || wood_data.thread.is_some()
        || clamp_data.thread.is_some()
        || wood_data.seat.is_some()
        || clamp_data.seat.is_some()
        || (wood_data.station_mm - grip).abs() > tolerance
        || clamp_data.station_mm.abs() > tolerance
    {
        return Err(unsupported());
    }
    let floor = wood_data.floor_mm.ok_or_else(unsupported)?;
    let mut materials = Vec::with_capacity(placed.len() + 2);
    for receiver in [wood, clamp] {
        let physical = evaluated
            .snapshot()
            .leaf_occurrences
            .get(&receiver.source.occurrence)
            .filter(|leaf| leaf.producer == receiver.feature)
            .ok_or_else(|| {
                fail("Hanger receivers must retain their current drilled physical sources")
            })?;
        related.insert(physical.path.clone());
        materials.push(
            evaluated
                .leaf_solid(&physical.path)
                .ok_or_else(|| fail("A current receiving part has no actual material source"))?,
        );
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
    // This is an intent gate only. Actual translated/tilted material proofs
    // below use original raw source frames, never this approximate comparison.
    at(
        axis,
        placed_frame(&hanger.part, hanger.pose),
        grip,
        tolerance,
        id,
    )?;
    let washers = placed
        .iter()
        .filter(|h| h.role == Role::NutWasher)
        .rev()
        .map(|h| h.solid)
        .collect::<Vec<_>>();
    checkpoint(interrupt)?;
    let proof = hanger
        .solid
        .rust_hanger_stack(
            materials[0],
            materials[1],
            nut.solid,
            &washers,
            StackRequirements {
                wood: Requirements {
                    minimum_wall_mm: wood_data.minimum_wall_mm,
                    minimum_floor_mm: floor,
                    minimum_engagement_mm: minimum,
                    tip_clearance_mm: tip,
                    linear_tolerance_mm: tolerance,
                    minimum_feature_mm: input.policy.minimum_feature_mm,
                },
                minimum_clamp_wall_mm: clamp_data.minimum_wall_mm,
                passage: PassageRequirements {
                    minimum_radial_clearance_mm: input.policy.minimum_feature_mm,
                    linear_tolerance_mm: tolerance,
                    minimum_feature_mm: input.policy.minimum_feature_mm,
                },
                mating: HangerMatingRequirements {
                    minimum_engagement_mm: minimum,
                    minimum_radial_clearance_mm: input.policy.minimum_feature_mm,
                    linear_tolerance_mm: tolerance,
                    ..Default::default()
                },
                minimum_bearing_land_mm: input.policy.minimum_feature_mm,
                maximum_bearing_gap_mm: tolerance,
            },
            interrupt,
        )
        .map_err(|e| error(Some(id), e))?;
    materials.extend(placed.iter().map(|h| h.solid));
    let mut pairs = 0;
    for path in evaluated.snapshot().leaf_occurrences.keys() {
        if related.contains(path) {
            continue;
        }
        let other = evaluated
            .leaf_solid(path)
            .ok_or_else(|| fail("A physical scene component has no actual material source"))?;
        for material in &materials {
            checkpoint(interrupt)?;
            if pairs >= 4096 {
                return Err(fail(
                    "Hanger scene clearance exceeds the 4096 pair query budget",
                ));
            }
            pairs += 1;
            if material
                .rust_hanger_component_separation(other, interrupt)
                .map_err(|e| error(Some(id), e))?
                != Some(true)
            {
                return Err(fail(
                    "Clearance from another component could not be verified for this geometry",
                ));
            }
        }
    }
    checkpoint(interrupt)?;
    Ok(ScrewConnectionCheck {
        connection:id,mode:source.mode,hardware_instances:placed.iter().map(|h|h.instance).collect(),
        receiver_features:source.receivers.iter().map(|r|r.feature).collect(),grip_mm:grip,
        minimum_engagement_mm:minimum,blind_tip_clearance_mm:tip,
        actual_engagement_mm:Some(proof.nut.engagement_mm[0].midpoint(proof.nut.engagement_mm[1])),
        actual_wood_engagement_mm:Some(proof.wood.engagement_mm[0].midpoint(proof.wood.engagement_mm[1])),
        native_component_collision_checked:true,
        scope:"Actual authored hanger and finite nut: pilot wall/floor and wood contact, complete clamp/washer passage, nonnegative bearing gaps, retained metric phase and strict source-owned separation from unrelated components. No strength or installation certification".into(),
    })
}
