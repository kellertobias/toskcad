//! Physical flat washer stacks on the two certified through-bore mouths.
use super::*;
use multicad_kernel::GeometryPolicy;

pub(super) struct Inputs<'a, 'solid> {
    pub placed: &'a [Hardware<'solid>],
    pub axis: ReferenceFrame,
    pub bore: rust_through::ThroughReceiverProof,
    pub grip: f64,
    pub head: [f64; 2],
    pub head_extent: f64,
    pub conical_bearing: Option<&'a super::conical::Proof>,
    pub nut: [f64; 2],
    pub shaft: f64,
    pub policy: GeometryPolicy,
    pub connection: EntityId,
}

fn fail(input: &Inputs<'_, '_>, message: &str) -> ModelError {
    error(Some(input.connection), message)
}

pub(super) fn actual_dimensions(
    row: &Hardware<'_>,
    axis: ReferenceFrame,
    policy: GeometryPolicy,
    id: EntityId,
) -> Result<[f64; 3], ModelError> {
    let fail = |message| error(Some(id), message);
    let ResolvedMechanicalPart::Washer {
        bore_radius_mm,
        outside_radius_mm,
        thickness_mm,
        ..
    } = row.part
    else {
        return Err(fail(
            "A washer role must contain actual flat washer geometry",
        ));
    };
    let (ring, actual) = row
        .solid
        .rust_annulus_data()
        .ok_or_else(|| fail("A washer role has no retained circular ring material"))?;
    let dimensions = [bore_radius_mm, outside_radius_mm, thickness_mm];
    if dimensions
        != [
            ring.inner_radius_mm(),
            ring.outer_radius_mm(),
            ring.height_mm(),
        ]
    {
        return Err(fail(
            "Actual washer material differs from its current definition",
        ));
    }
    // A generated Part can still have a placement. Check its actual retained
    // chart against the current source mount, at the full physical extent.
    let expected = placed_frame(&row.part, row.pose);
    let distance = |a: [f64; 3], b: [f64; 3]| (a[0] - b[0]).hypot(a[1] - b[1]).hypot(a[2] - b[2]);
    let normal = actual.normal();
    let residual = distance(actual.origin, expected.origin_mm)
        + outside_radius_mm * (distance(actual.u, expected.u) + distance(actual.v, expected.v))
        + thickness_mm * distance(normal, expected.normal());
    let tilt = cross(normal, axis.normal());
    let tilt = dot(tilt, tilt).sqrt();
    if !residual.is_finite()
        || residual > policy.linear_tolerance_mm / 8.
        || tilt * outside_radius_mm > policy.linear_tolerance_mm / 8.
        || dot(normal, axis.normal()) <= 0.
    {
        return Err(fail(
            "Actual washer placement cannot prove a flat bearing plane",
        ));
    }
    Ok(dimensions)
}

pub(super) fn annular_land(
    a: [f64; 2],
    b: [f64; 2],
    policy: GeometryPolicy,
    id: EntityId,
) -> Result<(), ModelError> {
    let fail = |message| error(Some(id), message);
    let margin = policy
        .minimum_feature_mm
        .max(16. * policy.linear_tolerance_mm);
    // This complete annular land fits both adjacent materials. The reserve
    // exceeds the admitted frame and reconstruction errors on both sides.
    if a[1].min(b[1]) - a[0].max(b[0]) < 5. * margin {
        return Err(fail(
            "Bearing components have no complete overlapping flat annular land",
        ));
    }
    Ok(())
}

pub(super) fn prove(input: Inputs<'_, '_>) -> Result<[f64; 2], ModelError> {
    if input.bore.head_seat.is_some_and(|s| s.conical) && input.conical_bearing.is_none() {
        return Err(fail(
            &input,
            "Conical head bearing requires a retained matching screw profile",
        ));
    }
    let mut sums = [0.; 2];
    for row in input.placed {
        let side = match row.role {
            Role::HeadWasher => 0,
            Role::NutWasher => 1,
            _ => continue,
        };
        sums[side] += actual_dimensions(row, input.axis, input.policy, input.connection)?[2];
    }
    if sums.into_iter().any(|sum| !sum.is_finite()) {
        return Err(fail(
            &input,
            "Washer stack dimensions exceed the numerical budget",
        ));
    }
    if input.conical_bearing.is_some() && sums[0] != 0. {
        return Err(fail(
            &input,
            "Flat head washers cannot replace certified conical bearing",
        ));
    }
    let stock = [input.bore.radius_mm, input.bore.protected_radius_mm];
    let margin = input
        .policy
        .minimum_feature_mm
        .max(16. * input.policy.linear_tolerance_mm);
    let head_station = input.bore.head_seat.map_or(0., |s| s.depth_mm);
    let head_stock = input
        .bore
        .head_seat
        .map_or(stock, |s| [input.bore.radius_mm, s.radius_mm]);
    if let Some(seat) = input.bore.head_seat
        && !seat.conical
        && input.head_extent + 5. * margin >= seat.radius_mm
    {
        return Err(fail(
            &input,
            "The actual head envelope cannot clear the counterbore wall",
        ));
    }
    for (role, mut station, mut previous) in [
        (Role::HeadWasher, head_station - sums[0], input.head),
        (Role::NutWasher, input.grip, stock),
    ] {
        if role == Role::HeadWasher && input.conical_bearing.is_some() {
            continue; // The actual finite cone proof replaces a planar head land.
        }
        for row in input.placed.iter().filter(|row| row.role == role) {
            let [inner, outer, height] =
                actual_dimensions(row, input.axis, input.policy, input.connection)?;
            if role == Role::HeadWasher
                && input
                    .bore
                    .head_seat
                    .is_some_and(|s| outer + 5. * margin >= s.radius_mm)
            {
                return Err(fail(
                    &input,
                    "Actual washer material cannot clear the counterbore wall",
                ));
            }
            let (_, actual) = row
                .solid
                .rust_annulus_data()
                .expect("checked retained ring");
            at(
                input.axis,
                ReferenceFrame {
                    origin_mm: actual.origin,
                    u: actual.u,
                    v: actual.v,
                },
                station,
                input.policy.linear_tolerance_mm / 8.,
                input.connection,
            )?;
            // A head washer must clear the neck blend as well as the plain
            // shank and thread crests. Nut washers only surround the shaft.
            let clearance = if role == Role::HeadWasher {
                input.head[0]
            } else {
                input.shaft
            };
            if inner - clearance < margin {
                return Err(fail(
                    &input,
                    "Washer bore interferes with the shaft or head neck",
                ));
            }
            annular_land(previous, [inner, outer], input.policy, input.connection)?;
            previous = [inner, outer];
            station += height;
        }
        annular_land(
            previous,
            if role == Role::HeadWasher {
                head_stock
            } else {
                input.nut
            },
            input.policy,
            input.connection,
        )?;
    }
    Ok(sums)
}
