//! Complete flat bearing chains inside an actual retained blind counterbore.
use super::*;
use multicad_kernel::{BlindThreadMating, rust_kernel::headed_screw::ScrewHead};

pub(super) fn prove(
    input: &connection_receivers::Inputs<'_>,
    placed: &[Hardware<'_>],
    screw: &Hardware<'_>,
    axis: ReferenceFrame,
    proof: &BlindThreadMating,
    thread_frame: ReferenceFrame,
    id: EntityId,
) -> Result<f64, ModelError> {
    let fail = |message| error(Some(id), message);
    let tolerance = input.policy.linear_tolerance_mm;
    let margin = input.policy.minimum_feature_mm.max(16. * tolerance);
    let (actual, frame) = screw
        .solid
        .rust_machine_screw_geometry()
        .ok_or_else(|| fail("Blind bearing requires actual regular screw source geometry"))?;
    let spec = actual.spec();
    let (minimum_radius, maximum_radius) = match spec.head {
        ScrewHead::SocketCap { radius_mm, .. } => (radius_mm, radius_mm),
        ScrewHead::ExternalHex {
            across_flats_mm, ..
        } => (across_flats_mm / 2., across_flats_mm / 3_f64.sqrt()),
    };
    let shaft = spec
        .shank_radius_mm
        .max(actual.thread().spec().major_radius_mm);
    let mut sum = 0.;
    for row in placed.iter().filter(|row| row.role == Role::HeadWasher) {
        sum += washers::actual_dimensions(row, axis, input.policy, id)?[2];
    }
    if !sum.is_finite() {
        return Err(fail("Washer stack exceeds the numerical budget"));
    }
    let seat_depth = proof.counterbore.map_or(0., |seat| seat.depth_mm);
    let head_station = seat_depth - sum;
    let head_frame = ReferenceFrame {
        origin_mm: std::array::from_fn(|i| {
            frame.origin[i] + frame.normal()[i] * actual.head_seat_z()
        }),
        u: frame.u,
        v: frame.v,
    };
    let expected = placed_frame(&screw.part, screw.pose);
    let distance = |a: [f64; 3], b: [f64; 3]| {
        (a[0] - b[0])
            .portable_hypot(a[1] - b[1])
            .portable_hypot(a[2] - b[2])
    };
    let extent = spec.head.height_mm() + spec.shank_length_mm + actual.thread().spec().length_mm;
    let residual = distance(head_frame.origin_mm, expected.origin_mm)
        + maximum_radius
            * (distance(head_frame.u, expected.u) + distance(head_frame.v, expected.v))
        + extent * distance(head_frame.normal(), expected.normal());
    let tilt = cross(head_frame.normal(), axis.normal());
    if !residual.is_finite()
        || residual > tolerance / 8.
        || dot(tilt, tilt).sqrt() * (maximum_radius + extent) > tolerance / 8.
    {
        return Err(fail(
            "Actual screw placement cannot prove a flat bearing plane",
        ));
    }
    at(axis, head_frame, head_station, tolerance / 8., id)?;
    at(
        axis,
        thread_frame,
        head_station + spec.shank_length_mm,
        tolerance / 8.,
        id,
    )?;
    if spec.shank_length_mm > 0.
        && proof
            .counterbore
            .is_none_or(|seat| head_station + spec.shank_length_mm + margin >= seat.depth_mm)
    {
        return Err(fail(
            "The actual plain shank must end inside the free counterbore",
        ));
    }
    if proof
        .counterbore
        .is_some_and(|seat| maximum_radius + 5. * margin >= seat.radius_mm)
    {
        return Err(fail(
            "The actual head envelope cannot clear the counterbore wall",
        ));
    }
    let mut previous = [shaft, minimum_radius];
    let mut position = head_station;
    for row in placed.iter().filter(|row| row.role == Role::HeadWasher) {
        let [inner, outer, height] = washers::actual_dimensions(row, axis, input.policy, id)?;
        if inner - shaft < margin {
            return Err(fail("Washer bore interferes with the shaft"));
        }
        if proof
            .counterbore
            .is_some_and(|seat| outer + 5. * margin >= seat.radius_mm)
        {
            return Err(fail(
                "Actual washer material cannot clear the counterbore wall",
            ));
        }
        let (_, frame) = row.solid.rust_annulus_data().expect("checked actual ring");
        at(
            axis,
            ReferenceFrame {
                origin_mm: frame.origin,
                u: frame.u,
                v: frame.v,
            },
            position,
            tolerance / 8.,
            id,
        )?;
        washers::annular_land(previous, [inner, outer], input.policy, id)?;
        previous = [inner, outer];
        position += height;
    }
    washers::annular_land(
        previous,
        [
            proof.threads[1].0.spec().major_radius_mm,
            proof
                .counterbore
                .map_or(proof.protected_radius_mm, |seat| seat.radius_mm),
        ],
        input.policy,
        id,
    )?;
    Ok(shaft)
}
