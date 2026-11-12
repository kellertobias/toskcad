//! Curved cavity patches; publication stays gated until physical readback.
use super::*;
use crate::hexalobular_drive::{Drive, Error};
#[cfg(test)]
mod tests;

fn source_error(error: Error) -> StepError {
    match error {
        Error::Invalid => StepError::InvalidGeometry,
        Error::Unresolved => StepError::Unsupported,
        Error::Limit | Error::Cancelled => StepError::Limit,
    }
}
pub(super) fn encode(
    file: &mut StepFile,
    p: &Placement,
    source: &Drive,
    top: f64,
    tolerance: f64,
) -> Result<(Entity, Vec<Entity>), StepError> {
    if !tolerance.is_finite() || tolerance <= 0. || !top.is_finite() {
        return Err(StepError::InvalidGeometry);
    }
    let tolerance = tolerance.min(MAX_ERROR_MM);
    let half_angle = source.cone().depth_per_radius.recip().atan();
    if half_angle <= 0. || half_angle >= std::f64::consts::FRAC_PI_2 {
        return Err(StepError::Unsupported);
    }
    let curves = (0..12)
        .map(|i| {
            source
                .intersection_carrier(i, tolerance / 16., 100_000)
                .map_err(source_error)
        })
        .collect::<Result<Vec<_>, _>>()?;
    let contacts = source.arcs().map(|a| a.start_mm);
    let mouth = contacts.map(|xy| file.vertex(p.point_f64(xy, top)));
    let floor = contacts.map(|xy| {
        file.vertex(p.point_f64(
            xy,
            top + source.cone().tip_depth_mm - source.cone().depth_per_radius * xy[0].hypot(xy[1]),
        ))
    });
    let posts = (0..12)
        .map(|i| line_edge(file, mouth[i], floor[i]))
        .collect::<Result<Vec<_>, _>>()?;
    let apex_point = p.point_f64([0.; 2], top + source.cone().tip_depth_mm);
    let apex_cartesian = file.point(apex_point);
    let [x, y] = contacts[0];
    let radius = x.hypot(y);
    let axis = file.axis(
        apex_cartesian,
        p.normal.map(|x| -x),
        p.vector([x / radius, y / radius, 0.]),
    );
    let cone = file.entity(format_args!(
        "CONICAL_SURFACE('',{axis},0.,{half_angle:.16E})"
    ));
    let context=file.entity(format_args!("(GEOMETRIC_REPRESENTATION_CONTEXT(2)PARAMETRIC_REPRESENTATION_CONTEXT()REPRESENTATION_CONTEXT('2D SPACE',''))"));
    let [seam, first_pcurve, last_pcurve] =
        apex_seam(file, context, cone, floor[0], apex_cartesian, apex_point)?;
    let mut top_edges = Vec::with_capacity(12);
    let mut lower_edges = Vec::with_capacity(12);
    let mut sides = Vec::with_capacity(12);
    for (i, arc) in source.arcs().iter().enumerate() {
        let next = (i + 1) % 12;
        let center = p.point_f64(arc.center_mm, top);
        let (top_edge, _) = circular_arc_edge(
            file,
            [mouth[i], mouth[next]],
            center,
            p.normal.map(|x| x * arc.sweep_rad.signum()),
            p.vector([arc.start_angle_rad.cos(), arc.start_angle_rad.sin(), 0.]),
            arc.radius_mm,
            arc.sweep_rad.abs(),
        )?;
        let point = file.point(center);
        let cylinder_axis = file.axis(point, p.normal, p.u);
        let r = arc.radius_mm;
        let cylinder = file.entity(format_args!(
            "CYLINDRICAL_SURFACE('',{cylinder_axis},{r:.16E})"
        ));
        let mut poles = vec![floor[i].cartesian];
        for (j, span) in curves[i].iter().enumerate() {
            for index in 1..4 {
                let entity = if j + 1 == curves[i].len() && index == 3 {
                    floor[next].cartesian
                } else {
                    let [x, y, z] = span.poles_mm[index];
                    file.point(p.point_f64([x, y], top + z))
                };
                poles.push(entity);
            }
        }
        let mut knots = vec!["0.".to_owned()];
        knots.extend(
            curves[i]
                .iter()
                .map(|s| format!("{:.16E}", s.parameters[1])),
        );
        let multiplicities = (0..knots.len())
            .map(|j| {
                if j == 0 || j + 1 == knots.len() {
                    "4"
                } else {
                    "3"
                }
            })
            .collect::<Vec<_>>();
        let spline = file.entity(format_args!(
            "B_SPLINE_CURVE_WITH_KNOTS('',3,{},.UNSPECIFIED.,.F.,.F.,({}),({}),.UNSPECIFIED.)",
            Entities(&poles),
            multiplicities.join(","),
            knots.join(",")
        ));
        let intersection = file.entity(format_args!(
            "INTERSECTION_CURVE('',{spline},({cylinder},{cone}),.CURVE_3D.)"
        ));
        let (a, b) = (floor[i].vertex, floor[next].vertex);
        let lower = file.entity(format_args!("EDGE_CURVE('',{a},{b},{intersection},.T.)"));
        // Face.same_sense selects the support normal, not topological uses.
        let wire = [
            (posts[i], true),
            (lower, true),
            (posts[next], false),
            (top_edge, false),
        ];
        let wire = wire
            .into_iter()
            .map(|(e, s)| file.oriented_edge(e, s))
            .collect::<Vec<_>>();
        let bound = file.bound(&wire, true);
        sides.push(file.face_on(&[bound], cylinder, arc.sweep_rad < 0.));
        top_edges.push(top_edge);
        lower_edges.push(lower);
    }
    let mouth_wire = top_edges
        .into_iter()
        .map(|e| file.oriented_edge(e, true))
        .collect::<Vec<_>>();
    let mouth_bound = file.bound(&mouth_wire, false);
    let mut floor_wire = lower_edges
        .into_iter()
        .rev()
        .map(|e| file.oriented_edge(e, false))
        .collect::<Vec<_>>();
    floor_wire.push(file.entity(format_args!("SEAM_EDGE('',*,*,{seam},.T.,{last_pcurve})")));
    floor_wire.push(file.entity(format_args!("SEAM_EDGE('',*,*,{seam},.F.,{first_pcurve})")));
    let bound = file.bound(&floor_wire, true);
    let floor_face = file.face_on(&[bound], cone, false);
    let mut faces = vec![floor_face];
    faces.extend(sides);
    Ok((mouth_bound, faces))
}

fn uv_definition(file: &mut StepFile, context: Entity, points: [[f64; 2]; 2]) -> Entity {
    let points =
        points.map(|[u, v]| file.entity(format_args!("CARTESIAN_POINT('',({u:.16E},{v:.16E}))")));
    let curve = file.entity(format_args!("POLYLINE('',{})", Entities(&points)));
    file.entity(format_args!(
        "DEFINITIONAL_REPRESENTATION('',({curve}),{context})"
    ))
}
fn apex_seam(
    file: &mut StepFile,
    context: Entity,
    cone: Entity,
    rim: Vertex,
    cartesian: Entity,
    point: Point,
) -> Result<[Entity; 3], StepError> {
    // A collapsed cone parameter row is a point, not a zero-length 3D edge.
    // STEP geometry_schema evaluated_degenerate_pcurve preserves that authority.
    let definition = uv_definition(file, context, [[0., 0.], [TAU, 0.]]);
    let evaluated = file.entity(format_args!(
        "EVALUATED_DEGENERATE_PCURVE('',{cone},{definition},{cartesian})"
    ));
    let vertex = file.entity(format_args!("VERTEX_POINT('',{evaluated})"));
    let apex = Vertex {
        point,
        cartesian,
        vertex,
    };
    let distance = (rim
        .point
        .into_iter()
        .zip(point)
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>())
    .sqrt();
    let direction = file.direction(direction(rim.point, apex.point)?);
    let vector = file.entity(format_args!("VECTOR('',{direction},1.)"));
    let start = rim.cartesian;
    let line = file.entity(format_args!("LINE('',{start},{vector})"));
    let pcurves = [0., TAU].map(|u| {
        let definition = uv_definition(file, context, [[u, distance], [u, 0.]]);
        file.entity(format_args!("PCURVE('',{cone},{definition})"))
    });
    let [first, last] = pcurves;
    let curve = file.entity(format_args!(
        "SEAM_CURVE('',{line},({first},{last}),.PCURVE_S1.)"
    ));
    let from = rim.vertex;
    let to = apex.vertex;
    let edge = file.entity(format_args!("EDGE_CURVE('',{from},{to},{curve},.T.)"));
    Ok([edge, first, last])
}
