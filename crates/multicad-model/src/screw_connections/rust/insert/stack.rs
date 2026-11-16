//! Separate actual clamp and wood receivers; no combined-envelope substitute.
use super::*;
fn unsupported(id: EntityId) -> ModelError {
    error(
        Some(id),
        "Rust clamped inserts require two actual generated roles, a flat counterbored through clamp, and one separate blind wood pilot",
    )
}
fn datum(
    axis: ReferenceFrame,
    bounds: [[f64; 2]; 3],
    station: f64,
    tolerance: f64,
    id: EntityId,
) -> Result<(), ModelError> {
    let n = axis.normal();
    let mut reach = 0_f64;
    let mut magnitude = 0_f64;
    for i in 0..3 {
        let expected = axis.origin_mm[i] + n[i] * station;
        reach = (reach
            + (bounds[i][0] - expected)
                .abs()
                .max((bounds[i][1] - expected).abs()))
        .next_up();
        magnitude = (magnitude + axis.origin_mm[i].abs() + n[i].abs() * station.abs()).next_up();
    }
    // Six finite multiply/add operations plus the outward L1 accumulation;
    // 64eps*max(1,sum|terms|) also dominates all subnormal absolute error.
    let reserve = (64. * f64::EPSILON * magnitude.max(1.)).next_up();
    if !reach.is_finite() || (reach + reserve).next_up() > tolerance {
        return Err(error(
            Some(id),
            "Actual receiving source datum disagrees with the current connection axis",
        ));
    }
    Ok(())
}
pub(super) fn check(
    input: &connection_receivers::Inputs<'_>,
    source: &model::ScrewConnectionRequest,
    evaluated: &Model,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<ScrewConnectionCheck, ModelError> {
    let id = source.connection;
    let fail = |s: &str| error(Some(id), s);
    let [clamp, wood] = source.receivers.as_slice() else {
        return Err(unsupported(id));
    };
    if source.mode != Mode::WoodInsert
        || source.hardware.len() != 2
        || [Role::Screw, Role::Insert].into_iter().any(|role| {
            source
                .hardware
                .iter()
                .filter(|row| row.role == role)
                .count()
                != 1
        })
        || !matches!(clamp.source.extent, ConnectionReceiverExtent::Through {})
        || !matches!(wood.source.extent, ConnectionReceiverExtent::Blind { .. })
        || !matches!(clamp.source.cut, model::ConnectionReceiverCut::Round { .. })
        || !matches!(wood.source.cut, model::ConnectionReceiverCut::Round { .. })
    {
        return Err(unsupported(id));
    }
    let mut placed = Vec::with_capacity(2);
    for row in &source.hardware {
        checkpoint(interrupt)?;
        let ScrewConnectionHardwareSource::Generated { part } = &row.source else {
            return Err(unsupported(id));
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
                "Clamped insert hardware must retain its current generated source",
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
        .expect("unique screw checked");
    let insert = placed
        .iter()
        .find(|row| row.role == Role::Insert)
        .expect("unique insert checked");
    if !matches!(screw.part, ResolvedMechanicalPart::MachineScrew { .. })
        || !matches!(&insert.part,ResolvedMechanicalPart::CustomConnection {parameters} if matches!(parameters.as_ref(),Custom::WoodInsert {..}))
    {
        return Err(unsupported(id));
    }
    for receiver in [clamp, wood] {
        if evaluated
            .snapshot()
            .leaf_occurrences
            .get(&receiver.source.occurrence)
            .is_none_or(|leaf| leaf.producer != receiver.feature)
        {
            return Err(fail(
                "Both receiving components must retain their actual current cut sources",
            ));
        }
    }
    let grip = model::resolve_screw_connection_grip(
        input.document,
        input.parameters,
        source,
        input.revision,
        input.snapshots,
    )
    .map_err(|e| error(Some(id), e))?;
    let c = clamp
        .source
        .resolve_current(
            input.document,
            input.parameters,
            input.revision,
            input.snapshots,
        )
        .map_err(|e| error(Some(id), e))?;
    let w = wood
        .source
        .resolve_current(
            input.document,
            input.parameters,
            input.revision,
            input.snapshots,
        )
        .map_err(|e| error(Some(id), e))?;
    let Some((false, diameter, depth)) = c.seat else {
        return Err(unsupported(id));
    };
    if grip <= 0.
        || c.station_mm != 0.
        || w.station_mm != grip
        || w.seat.is_some()
        || c.thread.is_some()
        || w.thread.is_some()
    {
        return Err(unsupported(id));
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
        placed_frame(&screw.part, screw.pose),
        depth,
        input.policy.linear_tolerance_mm,
        id,
    )?;
    at(
        axis,
        placed_frame(&insert.part, insert.pose),
        grip,
        input.policy.linear_tolerance_mm,
        id,
    )?;
    let clamp_solid = evaluated
        .leaf_solid(&clamp.source.occurrence)
        .ok_or_else(|| fail("Actual clamp material is missing"))?;
    let wood_solid = evaluated
        .leaf_solid(&wood.source.occurrence)
        .ok_or_else(|| fail("Actual wood material is missing"))?;
    let minimum = mm(&source.parameters.minimum_engagement, input.parameters, id)?;
    let clearance = mm(&source.parameters.blind_tip_clearance, input.parameters, id)?;
    let bearing = clamp_solid
        .rust_wood_insert_clamp_contact(screw.solid, insert.solid, c.minimum_wall_mm, interrupt)
        .map_err(|e| error(Some(id), e))?;
    if bearing.bore_radius_mm != c.diameter_mm / 2.
        || bearing.seat_radius_mm != diameter / 2.
        || bearing.seat_depth_mm != depth
        || bearing.source_height_mm != grip
    {
        return Err(fail(
            "Actual clamp source differs from the current receiver",
        ));
    }
    for (mouth, station) in bearing.bore_mouths_mm.into_iter().zip([0., grip]) {
        datum(axis, mouth, station, input.policy.linear_tolerance_mm, id)?;
    }
    let wood_contact = wood_solid
        .rust_inset_wood_insert_contact(
            insert.solid,
            Requirements {
                minimum_wall_mm: w.minimum_wall_mm,
                minimum_floor_mm: w.floor_mm.ok_or_else(|| unsupported(id))?,
                minimum_engagement_mm: minimum,
                tip_clearance_mm: clearance,
                linear_tolerance_mm: input.policy.linear_tolerance_mm,
                minimum_feature_mm: input.policy.minimum_feature_mm,
            },
            input.policy.linear_tolerance_mm,
            interrupt,
        )
        .map_err(|e| error(Some(id), e))?;
    if wood_contact.pilot_radius_mm != w.diameter_mm / 2.
        || Some(wood_contact.pilot_depth_mm) != w.depth_mm
    {
        return Err(fail("Actual wood pilot differs from the current receiver"));
    }
    datum(
        axis,
        wood_contact.pilot_mouth_mm,
        grip,
        input.policy.linear_tolerance_mm,
        id,
    )?;
    // Retain original successful providers first. Only their Unsupported scope
    // selects the additive actual-clock proof; all other errors propagate.
    match wood_solid.rust_wood_insert_screw_clearance(screw.solid, clearance, interrupt) {
        Ok(_) => {}
        Err(e) if e.kind == multicad_kernel::GeometryErrorKind::Unsupported => {
            wood_solid
                .rust_clocked_wood_insert_screw_clearance(screw.solid, clearance, interrupt)
                .map_err(|e| error(Some(id), e))?;
        }
        Err(e) => return Err(error(Some(id), e)),
    }
    let mating = match screw
        .solid
        .rust_wood_insert_mating(insert.solid, minimum, interrupt)
    {
        Ok(v) => v,
        Err(e) if e.kind == multicad_kernel::GeometryErrorKind::Unsupported => {
            screw
                .solid
                .rust_clocked_wood_insert_mating(insert.solid, minimum, interrupt)
                .map_err(|e| error(Some(id), e))?
                .mating
        }
        Err(e) => return Err(error(Some(id), e)),
    };
    let mut queries = 0;
    for path in evaluated.snapshot().leaf_occurrences.keys() {
        if *path == clamp.source.occurrence
            || *path == wood.source.occurrence
            || placed.iter().any(|row| path.root == row.instance)
        {
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
    Ok(ScrewConnectionCheck {
        connection: id,
        mode: source.mode,
        hardware_instances: placed.iter().map(|row| row.instance).collect(),
        receiver_features: vec![clamp.feature, wood.feature],
        grip_mm: grip,
        minimum_engagement_mm: minimum,
        blind_tip_clearance_mm: clearance,
        actual_engagement_mm: Some(mating.engagement_mm[0].midpoint(mating.engagement_mm[1])),
        actual_wood_engagement_mm: Some(
            wood_contact.contact.engagement_mm[0].midpoint(wood_contact.contact.engagement_mm[1]),
        ),
        actual_head_bearing_gap_mm: Some(bearing.bearing_face_gap_mm),
        actual_insert_inset_mm: Some(wood_contact.inset_mm),
        native_component_collision_checked: true,
        scope: format!(
            "Actual separate counterbored clamp and blind wood pilot: original insert wood bite and metric mating, whole hardware clearance, positive actual head land; measured head gap <= {:.6e} mm and insert inset [{:.6e}, {:.6e}] mm. All unrelated hardware-component collisions checked; receiving stock-stock disjointness is not certified; no washers or strength certification",
            bearing.bearing_face_gap_mm[1], wood_contact.inset_mm[0], wood_contact.inset_mm[1]
        ),
    })
}
