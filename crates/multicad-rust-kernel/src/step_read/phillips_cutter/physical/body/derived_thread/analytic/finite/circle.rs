//! Bounded finite circle-to-periodic-crest transfer. Original circle endpoint
//! vertices need not be exactly incident on the decimal circle: their radial
//! projection is an explicit derived trim, with a separately measured reserve.
//! A chord bound earns each endpoint angle error without atan. The requested
//! directed span stays strictly between0 and2pi even after those errors, so raw
//! normal and EDGE_CURVE sense select the same finite arc, not its complement.
use super::*;
pub(super) fn direction(
    body: &Body<'_>,
    shell: &Shell<'_, '_>,
    id: u32,
    work: &mut Work<'_>,
) -> FResult<[I; 3]> {
    let a = shell.args(id, "DIRECTION")?;
    numbers(&a[1])?;
    let n = numeric(body, id, work)?;
    Ok(n.try_into()
        .expect("validated three-coordinate DIRECTION schema"))
}

mod core;
pub(super) mod target;
use target::ArcTarget;
pub(super) fn qualify(
    body: &Body<'_>,
    shell: &Shell<'_, '_>,
    thread: &CommonChartPeriodicThreadRepair<'_, '_>,
    band: usize,
    stations: [f64; 2],
    record: &mut PeriodicFiniteCurve,
    work: &mut Work<'_>,
) -> FResult<()> {
    // Same actual edge/geometry ID proved by verify_prefix's skipped entry
    // crest and profiled::junction -> tangent_shank::arc. This producer accepts
    // one-start, non-rounded sources only; exit_crests=false, hence role0/band3.
    assert_eq!(record.role, 0);
    assert_eq!(band, 3);
    let radii = thread.source().radius_nodes_mm();
    assert_eq!(radii[band], radii[band + 1]);
    let target = ArcTarget::new(
        *thread.source(),
        [band, band + 1],
        radii[band],
        stations[0],
        work,
    )?;
    qualify_target(body, shell, thread.chart(), target, record, work)
}
pub(super) fn qualify_target(
    body: &Body<'_>,
    shell: &Shell<'_, '_>,
    chart: CommonMouthChart,
    target: ArcTarget,
    record: &mut PeriodicFiniteCurve,
    work: &mut Work<'_>,
) -> FResult<()> {
    let args = shell.args(record.curve_entity, "CIRCLE")?;
    let axis = shell.args(reference(&args[1])?, "AXIS2_PLACEMENT_3D")?;
    let origin = point(body, reference(&axis[1])?, work)?;
    let normal = direction(body, shell, reference(&axis[2])?, work)?;
    let reference_direction = direction(body, shell, reference(&axis[3])?, work)?;
    let radius = numeric(body, record.curve_entity, work)?[0];
    let raw = core::RawCircle::new(origin, normal, reference_direction, radius, work)?;
    let vertices = [
        point(body, record.vertex_point_entities[0], work)?,
        point(body, record.vertex_point_entities[1], work)?,
    ];
    let result = core::certify(
        &raw,
        core::Input {
            target,
            chart,
            vertices,
            same_sense: record.curve_same_sense,
        },
        work,
    )?;
    record.kind = PeriodicFiniteKind::Circle;
    record.distance_mm = result.distance;
    record.cell_derivative_mm = result.derivative;
    record.endpoint_incidence_reserve_mm = result.incidence;
    Ok(())
}
