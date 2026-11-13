//! Curved cavity patches; publication stays gated until physical readback.
use super::*;
use crate::{
    hexalobular_drive::{Error, FinitePost, FiniteVertices},
    screw_recess::ConicalHexalobular,
};
use std::borrow::Cow;
mod preflight;
#[cfg(test)]
mod tests;
pub(super) use preflight::{admission, placement};

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
    recess: &ConicalHexalobular,
    top: f64,
    tolerance: f64,
) -> Result<(Entity, Vec<Entity>), StepError> {
    if !tolerance.is_finite() || tolerance <= 0. || !top.is_finite() {
        return Err(StepError::InvalidGeometry);
    }
    let source = recess.source();
    let binding = recess.finite_binding();
    let tolerance = tolerance.min(MAX_ERROR_MM);
    let half_angle = binding
        .and_then(|b| b.literal_cone_semi_angle_rad())
        .unwrap_or_else(|| source.cone().depth_per_radius.recip().atan());
    if half_angle <= 0. || half_angle >= std::f64::consts::FRAC_PI_2 {
        return Err(StepError::Unsupported);
    }
    let (curves, vertices) = if let Some(binding) = binding {
        if binding.spatial_source_error_bound_mm() > tolerance / 8. {
            return Err(StepError::Limit);
        }
        (
            Cow::Borrowed(binding.carriers().as_slice()),
            *binding.vertices(),
        )
    } else {
        let curves = (0..12)
            .map(|i| {
                source
                    .intersection_carrier(i, tolerance / 16., 100_000)
                    .map_err(source_error)
            })
            .collect::<Result<Vec<_>, _>>()?;
        let mouth_mm = source.arcs().map(|a| [a.start_mm[0], a.start_mm[1], 0.]);
        let floor_mm = mouth_mm.map(|[x, y, _]| {
            [
                x,
                y,
                source.cone().tip_depth_mm - source.cone().depth_per_radius * x.hypot(y),
            ]
        });
        (
            Cow::Owned(curves),
            FiniteVertices {
                mouth_mm,
                floor_mm,
                apex_mm: [0., 0., source.cone().tip_depth_mm],
            },
        )
    };
    let placed = |[x, y, z]: Point| p.point_f64([x, y], top + z);
    let mouth = vertices.mouth_mm.map(|point| file.vertex(placed(point)));
    let floor = vertices.floor_mm.map(|point| file.vertex(placed(point)));
    let posts = (0..12)
        .map(|i| {
            if let Some(binding) = binding {
                post_edge(file, p, binding.post_lines()[i], [mouth[i], floor[i]], top)
            } else {
                // Authored generators share the cylinder's actual support
                // axis. Subtracting rounded world vertices would introduce a
                // real tilt at far placements, unrelated to STEP precision.
                post_edge(
                    file,
                    p,
                    FinitePost {
                        origin_mm: vertices.mouth_mm[i],
                        direction_mm_per_parameter: [0., 0., 1.],
                        parameters: [0., vertices.floor_mm[i][2] - vertices.mouth_mm[i][2]],
                    },
                    [mouth[i], floor[i]],
                    top,
                )
            }
        })
        .collect::<Result<Vec<_>, _>>()?;
    let apex_point = placed(vertices.apex_mm);
    let apex_cartesian = file.point(apex_point);
    let [x, y] = std::array::from_fn(|i| vertices.floor_mm[0][i] - vertices.apex_mm[i]);
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

fn post_edge(
    file: &mut StepFile,
    p: &Placement,
    post: FinitePost,
    [from, to]: [Vertex; 2],
    top: f64,
) -> Result<Entity, StepError> {
    let vector = p.vector(post.direction_mm_per_parameter);
    let magnitude = vector[0].hypot(vector[1]).hypot(vector[2]);
    if from.point == to.point || !magnitude.is_finite() || magnitude <= 0. {
        return Err(StepError::InvalidGeometry);
    }
    let [x, y, z] = post.origin_mm;
    let origin = file.point(p.point_f64([x, y], top + z));
    let direction = file.direction(vector);
    let vector = file.entity(format_args!("VECTOR('',{direction},{magnitude:.16E})"));
    let line = file.entity(format_args!("LINE('',{origin},{vector})"));
    let (a, b) = (from.vertex, to.vertex);
    let sense = if post.parameters[1] > post.parameters[0] {
        ".T."
    } else {
        ".F."
    };
    Ok(file.entity(format_args!("EDGE_CURVE('',{a},{b},{line},{sense})")))
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
