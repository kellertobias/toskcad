//! Plain blind round receivers retain actual pocket material and source charts.
use super::*;
use multicad_kernel::CircularPocketTool;
use multicad_precision::PortableMath as _;
pub(super) fn build(
    inputs: Inputs<'_>,
    id: EntityId,
    source: &ConnectionReceiverDefinition,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Solid, ModelError> {
    let fail = |message| error(Some(id), message);
    let resolved = source
        .resolve_current(
            inputs.document,
            inputs.parameters,
            inputs.revision,
            inputs.snapshots,
        )
        .map_err(|cause| error(Some(id), cause))?;
    let (Some(depth), Some(floor)) = (resolved.depth_mm, resolved.floor_mm) else {
        return Err(fail(
            "Plain blind receivers require an explicit depth and floor",
        ));
    };
    if resolved.thread.is_some() || resolved.seat.is_some_and(|s| s.0) {
        return Err(fail(
            "Rust plain blind receivers do not yet support conical seats or threads",
        ));
    }
    let pose = source
        .local_pose(
            inputs.document,
            inputs.parameters,
            inputs.revision,
            inputs.snapshots,
        )
        .map_err(|cause| error(Some(id), cause))?;
    let local = pose.frame(ReferenceFrame {
        origin_mm: [0.; 3],
        u: [1., 0., 0.],
        v: [0., 1., 0.],
    });
    let frame = PlaneFrame {
        origin: local.origin_mm,
        u: local.u,
        v: local.v,
    };
    let overrun = inputs
        .policy
        .minimum_feature_mm
        .max(8. * inputs.policy.linear_tolerance_mm);
    let n = frame.normal();
    let tool = CircularPocketTool {
        frame: PlaneFrame {
            origin: std::array::from_fn(|i| frame.origin[i] - n[i] * overrun),
            ..frame
        },
        center_uv: [0.; 2],
        radius_mm: resolved.diameter_mm / 2.,
        displacement_mm: n.map(|x| x * (depth + overrun)),
        origin: origin(id, "bore"),
    };
    let stock = inputs
        .solids
        .get(&source.input.producer)
        .ok_or_else(|| fail("Blind receiving predecessor has no current Rust solid"))?;
    if stock.rust_cylindrical_stock_available() {
        return stock
            .rust_cut_cylindrical_blind_pocket(
                multicad_kernel::CylindricalBlindPocketSpec {
                    radius_mm: resolved.diameter_mm / 2.,
                    depth_mm: depth,
                    minimum_wall_mm: resolved.minimum_wall_mm,
                    minimum_floor_mm: floor,
                    seat: resolved.seat.map(|(_, diameter, depth)| {
                        multicad_kernel::rust_kernel::cylindrical_blind_pocket::Seat {
                            radius_mm: diameter / 2.,
                            depth_mm: depth,
                        }
                    }),
                },
                frame,
                id,
                inputs.policy,
                interrupt,
            )
            .map_err(|cause| error(Some(id), cause));
    }
    if let Some((_, diameter, seat_depth)) = resolved.seat {
        return stock
            .rust_cut_rectangular_blind_counterbore(
                multicad_kernel::CylindricalBlindPocketSpec {
                    radius_mm: resolved.diameter_mm / 2.,
                    depth_mm: depth,
                    minimum_wall_mm: resolved.minimum_wall_mm,
                    minimum_floor_mm: floor,
                    seat: Some(
                        multicad_kernel::rust_kernel::cylindrical_blind_pocket::Seat {
                            radius_mm: diameter / 2.,
                            depth_mm: seat_depth,
                        },
                    ),
                },
                frame,
                id,
                inputs.policy,
                interrupt,
            )
            .map_err(|cause| error(Some(id), cause));
    }
    let (result, mouth) = stock
        .rust_connection_blind_bore(&tool, resolved.minimum_wall_mm, floor, inputs.policy)
        .map_err(|cause| error(Some(id), cause))?;
    let distance = (mouth[0] - frame.origin[0])
        .portable_hypot(mouth[1] - frame.origin[1])
        .portable_hypot(mouth[2] - frame.origin[2]);
    if distance > inputs.policy.linear_tolerance_mm {
        return Err(fail(
            "Blind receiving station does not match the actual stock entry face",
        ));
    }
    Ok(result)
}
