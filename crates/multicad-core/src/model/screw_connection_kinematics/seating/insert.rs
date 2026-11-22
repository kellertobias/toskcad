//! Bounded numerical proposals for managed generated clamped inserts.
//! Source geometry/dimensions stay intact; final Native proofs own acceptance.
use super::*;
fn extent(part: &ResolvedMechanicalPart) -> Result<(f64, f64), DocumentError> {
    match part {
        ResolvedMechanicalPart::MachineScrew {
            head,
            underhead_length_mm,
            ..
        } => Ok((
            head.maximum_radius_mm(),
            I::p(*underhead_length_mm).add(I::p(head.height_mm())).hi,
        )),
        ResolvedMechanicalPart::CustomConnection { parameters } => {
            if let model::ResolvedCustomConnectionHardware::WoodInsert { outer_thread, .. } =
                parameters.as_ref()
            {
                Ok((outer_thread.major_radius_mm, outer_thread.length_mm))
            } else {
                Err(invalid(
                    "Managed insert seating requires an actual insert definition",
                ))
            }
        }
        _ => Err(invalid(
            "Managed insert seating requires screw and insert definitions",
        )),
    }
}
pub(in crate::model::screw_connection_kinematics) fn pose(
    request: &ScrewConnectionRequest,
    parts: &[ResolvedMechanicalPart],
    role: EntityId,
    grip: f64,
    axis: RigidPose,
    p: &EvaluatedParameters,
    options: Options,
) -> Result<(f64, RigidPose), DocumentError> {
    axis.validate()?;
    if parts.len() != request.hardware.len()
        || parts.len() != 2
        || !grip.is_finite()
        || grip <= 0.
        || !options.maximum_gap_mm.is_finite()
        || options.maximum_gap_mm <= 0.
        || options.maximum_attempts == 0
        || options.maximum_attempts > 4
    {
        return Err(invalid(
            "Managed insert numerical seating inputs exceed their contract",
        ));
    }
    let bounds = parts.iter().map(extent).collect::<Result<Vec<_>, _>>()?;
    let radius = bounds.iter().map(|x| x.0).fold(0., f64::max);
    let height = bounds.iter().map(|x| x.1).fold(0., f64::max);
    let primary = request
        .hardware
        .iter()
        .zip(parts)
        .find(|(row, _)| row.role == Role::Screw)
        .ok_or_else(|| invalid("Managed insert seating requires one screw role"))?
        .1;
    let bearing = model::resolve_screw_connection_head_station(request, p, primary, 0.)?;
    let step = (grip.next_up() - grip).max(bearing.next_up() - bearing);
    let mut gap = 0.;
    let mut inset = 0.;
    for _ in 0..options.maximum_attempts {
        let head = I::p(bearing).sub(I::p(gap));
        let tail = I::p(grip).add(I::p(inset));
        let head_station = head.midpoint();
        let insert_station = tail.midpoint();
        let (start, lead, hand) =
            super::super::phase_reference(request, parts, head_station, grip)?
                .ok_or_else(|| invalid("Managed insert seating requires a metric phase source"))?;
        let roll = super::super::signed_turn(hand, start - insert_station, lead)?;
        let mut error = 0_f64;
        let mut poses = Vec::with_capacity(2);
        for (row, part) in request.hardware.iter().zip(parts) {
            let (target, angle_rad) = match row.role {
                Role::Screw => (head, roll),
                Role::Insert => (tail, 0.),
                _ => return Err(invalid("Managed insert seating has an incompatible role")),
            };
            let control = model::NumericTransform {
                translation: [length(0.), length(0.), length(target.midpoint())],
                angle: angle(angle_rad.to_degrees()),
                ..Default::default()
            }
            .resolve(p)?;
            let mount = part.mount_pose();
            let pose = mount.inverse()?.then(control)?.then(axis)?;
            error = error.max(residual(axis, pose, mount, target, radius, height));
            poses.push((row.role_id, target.midpoint(), pose));
        }
        if I::p(gap).add(I::p(2.).mul(I::p(error))).hi > options.maximum_gap_mm
            || I::p(inset).add(I::p(2.).mul(I::p(error))).hi > options.maximum_gap_mm
        {
            return Err(invalid(
                "Managed insert contact tolerance cannot contain actual pose rounding",
            ));
        }
        // Measured actual pose arithmetic plus the source station's own ULP.
        // These are positive proposals only; no stock or hardware is welded.
        let needed = I::p(4.).mul(I::p(error.max(step))).hi;
        if gap >= needed && inset >= needed && insert_station > grip && head_station < bearing {
            return poses
                .into_iter()
                .find(|x| x.0 == role)
                .map(|x| (x.1, x.2))
                .ok_or_else(|| invalid("Managed insert seating role is absent"));
        }
        gap = I::p(gap.max(needed)).mul(I::p(1.25)).hi;
        inset = I::p(inset.max(needed)).mul(I::p(1.25)).hi;
    }
    Err(invalid(
        "Managed insert numerical seating did not converge within four attempts",
    ))
}
