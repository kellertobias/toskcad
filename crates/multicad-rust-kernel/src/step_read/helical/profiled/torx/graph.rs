//! Literal cavity incidence, without a source recipe or coordinate fitting.
use super::*;

pub(super) struct Cell {
    pub mouth: [u32; 2],
    pub floor: [u32; 2],
    pub circle: u32,
    pub circle_forward: bool,
    pub cylinder: u32,
    pub wall_forward: bool,
    pub post: u32,
    pub post_forward: bool,
    pub spline: u32,
    pub lower: u32,
}
pub(super) struct Graph {
    pub cells: Vec<Cell>,
    pub cone: u32,
    pub apex: u32,
    pub seam: u32,
    pub seam_pcurves: [(bool, u32); 2],
}

pub(super) fn seam(shell: &Shell<'_, '_>, graph: &Graph) -> Result<(u32, f64)> {
    let (vertices, geometry, _) = edge_record(shell, graph.seam)?;
    require(vertices[1] == graph.apex)?;
    let [Value::String(_), line, associated, Value::Enum("PCURVE_S1")] =
        shell.args(geometry, "SEAM_CURVE")?
    else {
        return Err(ReadError::Unsupported);
    };
    let pcurves = refs(associated)?;
    require(pcurves.len() == 2)?;
    let mut length = None;
    for (i, &pcurve) in pcurves.iter().enumerate() {
        let [Value::String(_), surface, definition] = shell.args(pcurve, "PCURVE")? else {
            return Err(ReadError::Invalid);
        };
        require(reference(surface)? == graph.cone)?;
        let uv = finite_graph::uv_polyline(shell, reference(definition)?)?;
        let u = if i == 0 { 0. } else { TAU };
        require(
            uv[0][0] == u
                && uv[1] == [u, 0.]
                && uv[0][1] > 0.
                && length.is_none_or(|previous| previous == uv[0][1]),
        )?;
        length = Some(uv[0][1]);
    }
    for (forward, pcurve) in graph.seam_pcurves {
        require(pcurve == pcurves[usize::from(forward)])?;
    }
    Ok((reference(line)?, length.ok_or(ReadError::Invalid)?))
}

pub(super) fn read(
    shell: &mut Shell<'_, '_>,
    mouth_bound: u32,
    floor_face: u32,
    wall_faces: &[u32],
) -> Result<Graph> {
    require(wall_faces.len() == 12)?;
    let mouth = finite_graph::wire(shell, mouth_bound, false)?;
    require(mouth.edges.len() == 12)?;
    let (cone, forward, floor_bounds) = face(shell, floor_face)?;
    require(!forward && floor_bounds.len() == 1)?;
    require(shell.simple(cone)?.name == "CONICAL_SURFACE")?;
    let floor = finite_graph::wire(shell, floor_bounds[0], true)?;
    require(floor.edges.len() == 14)?;
    let mut used_walls = BTreeSet::new();
    let mut cells = Vec::with_capacity(12);
    for &(mouth_edge, mouth_forward) in &mouth.edges {
        let (mut vertices, circle, circle_forward) = edge_record(shell, mouth_edge)?;
        require(shell.simple(circle)?.name == "CIRCLE")?;
        if !mouth_forward {
            vertices.reverse();
        }
        let mut selected = None;
        for &wall_face in wall_faces {
            let (cylinder, wall_forward, bounds) = face(shell, wall_face)?;
            require(bounds.len() == 1)?;
            let wire = finite_graph::wire(shell, bounds[0], true)?;
            if let Some(index) = wire.edges.iter().position(|&(id, _)| id == mouth_edge) {
                require(selected.is_none() && used_walls.insert(wall_face))?;
                require(wire.edges.len() == 4 && wire.edges[index].1 != mouth_forward)?;
                require(shell.simple(cylinder)?.name == "CYLINDRICAL_SURFACE")?;
                let [post, lower, next_post] =
                    std::array::from_fn(|j| wire.edges[(index + j + 1) % 4]);
                require(post.1 && lower.1 && !next_post.1)?;
                let (post_vertices, post_curve, post_forward) = edge_record(shell, post.0)?;
                let (next_vertices, _, _) = edge_record(shell, next_post.0)?;
                let (lower_vertices, intersection, lower_forward) = edge_record(shell, lower.0)?;
                require(
                    post_vertices[0] == vertices[0]
                        && next_vertices[0] == vertices[1]
                        && lower_vertices == [post_vertices[1], next_vertices[1]],
                )?;
                require(shell.simple(post_curve)?.name == "LINE" && lower_forward)?;
                let [Value::String(_), spline, supports, Value::Enum("CURVE_3D")] =
                    shell.args(intersection, "INTERSECTION_CURVE")?
                else {
                    return Err(ReadError::Unsupported);
                };
                let supports = refs(supports)?;
                require(
                    supports.len() == 2
                        && supports.contains(&cone)
                        && supports.contains(&cylinder)
                        && cone != cylinder,
                )?;
                require(
                    floor
                        .edges
                        .iter()
                        .filter(|&&e| e == (lower.0, false))
                        .count()
                        == 1,
                )?;
                selected = Some(Cell {
                    mouth: vertices,
                    floor: lower_vertices,
                    circle,
                    circle_forward: circle_forward == mouth_forward,
                    cylinder,
                    wall_forward,
                    post: post_curve,
                    post_forward,
                    spline: reference(spline)?,
                    lower: lower.0,
                });
            }
        }
        cells.push(selected.ok_or(ReadError::Invalid)?);
    }
    require(used_walls.len() == 12)?;
    for i in 0..12 {
        require(
            cells[i].mouth[1] == cells[(i + 1) % 12].mouth[0]
                && cells[i].floor[1] == cells[(i + 1) % 12].floor[0],
        )?;
    }
    let lower_ids: BTreeSet<_> = cells.iter().map(|cell| cell.lower).collect();
    require(lower_ids.len() == 12)?;
    let seam: Vec<_> = floor
        .edges
        .iter()
        .zip(&floor.seam_pcurves)
        .filter(|(edge, _)| !lower_ids.contains(&edge.0))
        .collect();
    require(seam.len() == 2 && seam[0].0.0 == seam[1].0.0 && seam[0].0.1 != seam[1].0.1)?;
    let seam_id = seam[0].0.0;
    let (vertices, geometry, same_sense) = edge_record(shell, seam_id)?;
    require(same_sense && cells.iter().any(|cell| cell.floor[0] == vertices[0]))?;
    require(shell.simple(geometry)?.name == "SEAM_CURVE")?;
    Ok(Graph {
        cells,
        cone,
        apex: vertices[1],
        seam: seam_id,
        seam_pcurves: [
            (seam[0].0.1, seam[0].1.ok_or(ReadError::Invalid)?),
            (seam[1].0.1, seam[1].1.ok_or(ReadError::Invalid)?),
        ],
    })
}

pub(super) fn spline(
    shell: &Shell<'_, '_>,
    id: u32,
) -> Result<Vec<crate::hexalobular_drive::CubicIntersectionSpan>> {
    let [
        Value::String(_),
        Value::Number(3.),
        points,
        Value::Enum("UNSPECIFIED"),
        Value::Enum("F"),
        Value::Enum("F"),
        Value::List(multiplicities),
        Value::List(knots),
        Value::Enum("UNSPECIFIED"),
    ] = shell.args(id, "B_SPLINE_CURVE_WITH_KNOTS")?
    else {
        return Err(ReadError::Unsupported);
    };
    require(knots.len() >= 2 && knots.len() <= 100_001 && multiplicities.len() == knots.len())?;
    let mut parameters = Vec::with_capacity(knots.len());
    for (i, (knot, multiplicity)) in knots.iter().zip(multiplicities).enumerate() {
        let Value::Number(knot) = knot else {
            return Err(ReadError::Invalid);
        };
        require(
            knot.is_finite()
                && *knot >= 0.
                && *knot <= 1.
                && multiplicity
                    == &Value::Number(if i == 0 || i + 1 == knots.len() {
                        4.
                    } else {
                        3.
                    }),
        )?;
        parameters.push(*knot);
    }
    require(
        parameters[0] == 0.
            && *parameters.last().ok_or(ReadError::Invalid)? == 1.
            && parameters.windows(2).all(|pair| pair[0] < pair[1]),
    )?;
    let points = refs(points)?;
    require(points.len() == 3 * (parameters.len() - 1) + 1)?;
    let points = points
        .into_iter()
        .map(|id| shell.point(id))
        .collect::<Result<Vec<_>>>()?;
    Ok(parameters
        .windows(2)
        .enumerate()
        .map(
            |(i, pair)| crate::hexalobular_drive::CubicIntersectionSpan {
                parameters: [pair[0], pair[1]],
                poles_mm: std::array::from_fn(|j| points[3 * i + j]),
            },
        )
        .collect())
}
