use super::*;
mod bearing;

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
            "Rust blind tapped connections require one generated metric screw, optional flat head washers, no grip stack and one blind tap",
        )
    };
    let [receiver] = source.receivers.as_slice() else {
        return Err(unsupported());
    };
    if !(1..=64).contains(&source.hardware.len())
        || source
            .hardware
            .iter()
            .filter(|r| r.role == Role::Screw)
            .count()
            != 1
        || source
            .hardware
            .iter()
            .any(|r| !matches!(r.role, Role::Screw | Role::HeadWasher))
    {
        return Err(unsupported());
    }
    let mut placed = Vec::with_capacity(source.hardware.len());
    for row in &source.hardware {
        checkpoint(interrupt)?;
        let ScrewConnectionHardwareSource::Generated { part } = &row.source else {
            return Err(unsupported());
        };
        if !evaluated
            .snapshot()
            .leaf_occurrences
            .values()
            .any(|leaf| leaf.path.root == row.instance && leaf.producer == part.feature)
        {
            return Err(unsupported());
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
    let ResolvedMechanicalPart::MachineScrew {
        thread: external, ..
    } = &screw.part
    else {
        return Err(unsupported());
    };

    let physical = evaluated
        .snapshot()
        .leaf_occurrences
        .get(&receiver.source.occurrence)
        .filter(|leaf| leaf.producer == receiver.feature)
        .ok_or_else(|| fail("Blind tap must be the current physical receiving output"))?;
    let receiving = evaluated
        .leaf_solid(&physical.path)
        .ok_or_else(|| fail("Blind tap has no current physical solid"))?;
    let proof = screw
        .solid
        .rust_screw_blind_thread(receiving)
        .ok_or_else(unsupported)?;
    if !proof.entry_is_support_plane {
        return Err(fail(
            "Rust cannot prove the entire screw head clears the blind tap entry plane",
        ));
    }
    let grip = model::resolve_screw_connection_grip(
        input.document,
        input.parameters,
        source,
        input.revision,
        input.snapshots,
    )
    .map_err(|cause| error(Some(id), cause))?;
    let tolerance = input.policy.linear_tolerance_mm;
    if !grip.is_finite() || grip.abs() > tolerance {
        return Err(unsupported());
    }
    let minimum = mm(&source.parameters.minimum_engagement, input.parameters, id)?;
    let blind_clearance = mm(&source.parameters.blind_tip_clearance, input.parameters, id)?;
    if !blind_clearance.is_finite() || blind_clearance < 0. {
        return Err(fail("Blind tip clearance must be finite and nonnegative"));
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
        .map_err(|cause| error(Some(id), cause))?
        .frame(ReferenceFrame {
            origin_mm: [0.; 3],
            u: [1., 0., 0.],
            v: [0., 1., 0.],
        });
    let [(screw_thread, screw_frame), (tap_thread, tap_frame)] = proof.threads;
    let reference = |frame: PlaneFrame| ReferenceFrame {
        origin_mm: frame.origin,
        u: frame.u,
        v: frame.v,
    };
    let screw_frame = reference(screw_frame);
    let tap_frame = reference(tap_frame);
    let shaft = bearing::prove(input, &placed, screw, axis, &proof, screw_frame, id)?;
    at(axis, tap_frame, 0., tolerance, id)?;
    let (receiver_station, _) = station(axis, tap_frame.origin_mm);
    let (screw_station, _) = station(axis, screw_frame.origin_mm);
    if screw_station + screw_thread.spec().length_mm
        > receiver_station + tap_thread.spec().length_mm - blind_clearance + tolerance
    {
        return Err(fail("Screw tip violates the declared blind-hole clearance"));
    }
    let resolved = receiver
        .source
        .resolve_current(
            input.document,
            input.parameters,
            input.revision,
            input.snapshots,
        )
        .map_err(|cause| error(Some(id), cause))?;
    let expected_seat = resolved
        .seat
        .map(|(conical, diameter, depth)| (conical, diameter / 2., depth));
    let actual_seat = proof
        .counterbore
        .map(|seat| (false, seat.radius_mm, seat.depth_mm));
    if expected_seat != actual_seat {
        return Err(fail(
            "Actual blind entry seat differs from its current definition",
        ));
    }
    let internal = resolved.thread.ok_or_else(unsupported)?;
    let mating = model::resolve_metric_thread_mating_clipped(
        external,
        screw_frame,
        &internal,
        tap_frame,
        minimum,
        tolerance,
        model::MetricThreadMatingIntervals {
            external: None,
            internal: proof
                .counterbore
                .map(|seat| [seat.depth_mm, tap_thread.spec().length_mm]),
        },
    )
    .map_err(|cause| error(Some(id), cause))?;
    no_phase_correction(&mating, id)?;
    check_radial_clearance(
        id,
        [(screw_thread, screw_frame), (tap_thread, tap_frame)],
        mating.phase.axial_station_mm,
        shaft,
        tolerance,
    )?;
    for hardware in &placed {
        let bounds = hardware
            .solid
            .bounds()
            .map_err(|cause| error(Some(id), cause))?;
        for leaf in evaluated.snapshot().leaf_occurrences.values() {
            checkpoint(interrupt)?;
            if leaf.path == receiver.source.occurrence
                || placed.iter().any(|row| leaf.path.root == row.instance)
            {
                continue;
            }
            if !(0..3).any(|i| {
                bounds.max[i] + tolerance < leaf.bounds.min[i]
                    || leaf.bounds.max[i] + tolerance < bounds.min[i]
            }) {
                return Err(fail(
                    "Rust cannot certify clearance to another nearby physical component",
                ));
            }
        }
    }
    Ok(ScrewConnectionCheck {
        connection:id,
        mode:source.mode,
        hardware_instances:placed.iter().map(|row| row.instance).collect(),
        receiver_features:vec![receiver.feature],
        grip_mm:grip,
        minimum_engagement_mm:minimum,
        blind_tip_clearance_mm:blind_clearance,
        actual_engagement_mm:Some(mating.overlap_length_mm),
        actual_wood_engagement_mm:None,
        actual_head_bearing_gap_mm: None,
        actual_insert_inset_mm: None,
        native_component_collision_checked:true,
        scope:"Exact placed metric profiles, blind-floor and tip clearance, flat seating and separation from other components".into(),
    })
}
