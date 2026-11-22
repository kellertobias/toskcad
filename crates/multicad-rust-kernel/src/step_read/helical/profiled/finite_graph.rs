//! Finite topological uses and genuine collapsed cone parameter vertices.
//! Full carrier/support qualification remains the owning head proof's job.
use super::*;

pub(super) struct Wire {
    pub(super) edges: Vec<EdgeRef>,
    pub(super) uses: Vec<u32>,
    pub(super) seam_pcurves: Vec<Option<u32>>,
}
pub(super) fn wire(shell: &Shell<'_, '_>, bound: u32, outside: bool) -> Result<Wire> {
    let (outer, loop_id, forward) = shell.bound(bound, false)?;
    require(outer == outside)?;
    let [Value::String(_), edges] = shell.args(loop_id, "EDGE_LOOP")? else {
        return Err(ReadError::Invalid);
    };
    let mut uses = refs(edges)?;
    let mut result = Vec::new();
    let mut seam_pcurves = Vec::new();
    for &id in &uses {
        let entity = shell.simple(id)?;
        let (edge, sense, pcurve) = match (entity.name, entity.args.as_slice()) {
            (
                "ORIENTED_EDGE",
                [
                    Value::String(_),
                    Value::Derived,
                    Value::Derived,
                    edge,
                    forward,
                ],
            ) => (reference(edge)?, sense(forward)?, None),
            (
                "SEAM_EDGE",
                [
                    Value::String(_),
                    Value::Derived,
                    Value::Derived,
                    edge,
                    forward,
                    pcurve,
                ],
            ) => {
                let edge = reference(edge)?;
                let (_, geometry, _) = edge_record(shell, edge)?;
                let [
                    Value::String(_),
                    _,
                    associated,
                    Value::Enum("PCURVE_S1" | "PCURVE_S2"),
                ] = shell.args(geometry, "SEAM_CURVE")?
                else {
                    return Err(ReadError::Invalid);
                };
                let curves = refs(associated)?;
                let pcurve = reference(pcurve)?;
                require(
                    curves.len() == 2
                        && curves[0] != curves[1]
                        && curves.contains(&pcurve)
                        && shell.simple(pcurve)?.name == "PCURVE",
                )?;
                (edge, sense(forward)?, Some(pcurve))
            }
            _ => return Err(ReadError::Unsupported),
        };
        result.push((edge, sense));
        seam_pcurves.push(pcurve);
    }
    if !forward {
        uses.reverse();
        seam_pcurves.reverse();
    }
    Ok(Wire {
        edges: orient(&result, forward),
        uses,
        seam_pcurves,
    })
}

pub(super) fn uv_polyline(shell: &Shell<'_, '_>, id: u32) -> Result<[[f64; 2]; 2]> {
    let [Value::String(_), items, context] = shell.args(id, "DEFINITIONAL_REPRESENTATION")? else {
        return Err(ReadError::Invalid);
    };
    require(
        shell
            .components(reference(context)?)?
            .iter()
            .any(|c| c.name == "GEOMETRIC_REPRESENTATION_CONTEXT" && c.args == [Value::Number(2.)]),
    )?;
    let [curve] = refs(items)?.try_into().map_err(|_| ReadError::Invalid)?;
    let [Value::String(_), points] = shell.args(curve, "POLYLINE")? else {
        return Err(ReadError::Invalid);
    };
    let points: [u32; 2] = refs(points)?.try_into().map_err(|_| ReadError::Invalid)?;
    let mut output = [[0.; 2]; 2];
    for (i, id) in points.into_iter().enumerate() {
        let [Value::String(_), Value::List(coordinates)] = shell.args(id, "CARTESIAN_POINT")?
        else {
            return Err(ReadError::Invalid);
        };
        let [Value::Number(u), Value::Number(v)] = coordinates.as_slice() else {
            return Err(ReadError::Invalid);
        };
        require(
            [*u, *v]
                .into_iter()
                .all(|n| n.is_finite() && n.abs() <= 1e12),
        )?;
        output[i] = [*u, *v];
    }
    Ok(output)
}

pub(super) fn vertex(shell: &mut Shell<'_, '_>, id: u32) -> Result<Point> {
    let [Value::String(_), geometry] = shell.args(id, "VERTEX_POINT")? else {
        return Err(ReadError::Invalid);
    };
    let geometry = reference(geometry)?;
    if shell.simple(geometry)?.name != "EVALUATED_DEGENERATE_PCURVE" {
        return shell.vertex(id);
    }
    let [Value::String(_), surface, definition, equivalent] =
        shell.args(geometry, "EVALUATED_DEGENERATE_PCURVE")?
    else {
        return Err(ReadError::Invalid);
    };
    let [
        Value::String(_),
        axis,
        Value::Number(0.),
        Value::Number(angle),
    ] = shell.args(reference(surface)?, "CONICAL_SURFACE")?
    else {
        return Err(ReadError::Unsupported);
    };
    require(*angle > 0. && *angle < PI / 2.)?;
    let apex = shell.placement(reference(axis)?)?.origin;
    let point = shell.point(reference(equivalent)?)?;
    require(point == apex && uv_polyline(shell, reference(definition)?)? == [[0., 0.], [TAU, 0.]])?;
    shell.points.insert(id, point);
    Ok(point)
}
