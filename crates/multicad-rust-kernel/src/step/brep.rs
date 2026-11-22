//! Retained imported boundary written back with its own supports and trims.
use super::*;
use crate::brep::{Curve, Loop, NurbsCurve, NurbsSurface, Solid, Surface, surface::Frame};

struct Reals<'a>(&'a [f64]);

impl fmt::Display for Reals<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("(")?;
        for (i, value) in self.0.iter().enumerate() {
            write!(f, "{}{value:.16E}", if i == 0 { "" } else { "," })?;
        }
        f.write_str(")")
    }
}

struct Counts<'a>(&'a [usize]);

impl fmt::Display for Counts<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("(")?;
        for (i, value) in self.0.iter().enumerate() {
            write!(f, "{}{value}", if i == 0 { "" } else { "," })?;
        }
        f.write_str(")")
    }
}

fn frame(file: &mut StepFile, frame: &Frame) -> Entity {
    let origin = file.point(frame.origin);
    file.axis(origin, frame.z, frame.x)
}

/// A direction scaled by its length, as `VECTOR`.
fn vector(file: &mut StepFile, vector: Point) -> Entity {
    let length = vector[0].hypot(vector[1]).hypot(vector[2]);
    let direction = file.direction(vector.map(|x| x / length));
    file.entity(format_args!("VECTOR('',{direction},{length:.16E})"))
}

fn points(file: &mut StepFile, points: &[Point]) -> Vec<Entity> {
    points.iter().map(|&point| file.point(point)).collect()
}

fn spline_curve(file: &mut StepFile, curve: &NurbsCurve) -> Entity {
    let poles = points(file, curve.points());
    let (degree, poles) = (curve.degree(), Entities(&poles));
    let (knots, multiplicities) = curve.knots();
    let (knots, multiplicities) = (Reals(knots), Counts(multiplicities));
    match curve.weights() {
        None => file.entity(format_args!(
            "B_SPLINE_CURVE_WITH_KNOTS('',{degree},{poles},.UNSPECIFIED.,.F.,.F.,{multiplicities},{knots},.UNSPECIFIED.)"
        )),
        Some(weights) => file.entity(format_args!(
            "(BOUNDED_CURVE()B_SPLINE_CURVE({degree},{poles},.UNSPECIFIED.,.F.,.F.)B_SPLINE_CURVE_WITH_KNOTS({multiplicities},{knots},.UNSPECIFIED.)CURVE()GEOMETRIC_REPRESENTATION_ITEM()RATIONAL_B_SPLINE_CURVE({})REPRESENTATION_ITEM(''))",
            Reals(weights)
        )),
    }
}

fn curve(file: &mut StepFile, curve: &Curve) -> Entity {
    match curve {
        Curve::Line { origin, direction } => {
            let origin = file.point(*origin);
            let direction = vector(file, *direction);
            file.entity(format_args!("LINE('',{origin},{direction})"))
        }
        Curve::Circle {
            frame: placement,
            radius,
        } => {
            let placement = frame(file, placement);
            file.entity(format_args!("CIRCLE('',{placement},{radius:.16E})"))
        }
        Curve::Ellipse {
            frame: placement,
            semiaxes: [a, b],
        } => {
            let placement = frame(file, placement);
            file.entity(format_args!("ELLIPSE('',{placement},{a:.16E},{b:.16E})"))
        }
        Curve::Hyperbola {
            frame: placement,
            semiaxes: [a, b],
        } => {
            let placement = frame(file, placement);
            file.entity(format_args!("HYPERBOLA('',{placement},{a:.16E},{b:.16E})"))
        }
        Curve::Parabola {
            frame: placement,
            focal,
        } => {
            let placement = frame(file, placement);
            file.entity(format_args!("PARABOLA('',{placement},{focal:.16E})"))
        }
        Curve::Nurbs { curve, .. } => spline_curve(file, curve),
    }
}

fn spline_surface(file: &mut StepFile, surface: &NurbsSurface) -> Entity {
    let [rows, columns] = surface.counts();
    let poles = points(file, surface.points());
    let mut net = String::from("(");
    let mut weights = String::from("(");
    for row in 0..rows {
        let separator = if row == 0 { "" } else { "," };
        let range = row * columns..(row + 1) * columns;
        write!(net, "{separator}{}", Entities(&poles[range.clone()])).expect("string write");
        if let Some(all) = surface.weights() {
            write!(weights, "{separator}{}", Reals(&all[range])).expect("string write");
        }
    }
    net.push(')');
    weights.push(')');
    let [du, dv] = surface.degrees();
    let ((ku, mu), (kv, mv)) = (surface.knots(0), surface.knots(1));
    let knots = format!("{},{},{},{}", Counts(mu), Counts(mv), Reals(ku), Reals(kv));
    match surface.weights() {
        None => file.entity(format_args!(
            "B_SPLINE_SURFACE_WITH_KNOTS('',{du},{dv},{net},.UNSPECIFIED.,.F.,.F.,.F.,{knots},.UNSPECIFIED.)"
        )),
        Some(_) => file.entity(format_args!(
            "(BOUNDED_SURFACE()B_SPLINE_SURFACE({du},{dv},{net},.UNSPECIFIED.,.F.,.F.,.F.)B_SPLINE_SURFACE_WITH_KNOTS({knots},.UNSPECIFIED.)GEOMETRIC_REPRESENTATION_ITEM()RATIONAL_B_SPLINE_SURFACE({weights})REPRESENTATION_ITEM('')SURFACE())"
        )),
    }
}

fn surface(file: &mut StepFile, surface: &Surface) -> Entity {
    match surface {
        Surface::Plane(placement) => {
            let placement = frame(file, placement);
            file.entity(format_args!("PLANE('',{placement})"))
        }
        Surface::Cylinder {
            frame: placement,
            radius,
        } => {
            let placement = frame(file, placement);
            file.entity(format_args!(
                "CYLINDRICAL_SURFACE('',{placement},{radius:.16E})"
            ))
        }
        Surface::Cone {
            frame: placement,
            radius,
            half_angle,
        } => {
            let placement = frame(file, placement);
            file.entity(format_args!(
                "CONICAL_SURFACE('',{placement},{radius:.16E},{half_angle:.16E})"
            ))
        }
        Surface::Sphere {
            frame: placement,
            radius,
        } => {
            let placement = frame(file, placement);
            file.entity(format_args!(
                "SPHERICAL_SURFACE('',{placement},{radius:.16E})"
            ))
        }
        Surface::Torus {
            frame: placement,
            major,
            minor,
        } => {
            let placement = frame(file, placement);
            file.entity(format_args!(
                "TOROIDAL_SURFACE('',{placement},{major:.16E},{minor:.16E})"
            ))
        }
        Surface::Extrusion {
            curve: profile,
            direction,
        } => {
            let profile = curve(file, profile);
            let direction = vector(file, *direction);
            file.entity(format_args!(
                "SURFACE_OF_LINEAR_EXTRUSION('',{profile},{direction})"
            ))
        }
        Surface::Revolution {
            curve: profile,
            origin,
            axis,
        } => {
            let profile = curve(file, profile);
            let origin = file.point(*origin);
            let axis = file.direction(*axis);
            let axis = file.entity(format_args!("AXIS1_PLACEMENT('',{origin},{axis})"));
            file.entity(format_args!("SURFACE_OF_REVOLUTION('',{profile},{axis})"))
        }
        Surface::Nurbs { surface, .. } => spline_surface(file, surface),
    }
}

fn sense(forward: bool) -> &'static str {
    if forward { ".T." } else { ".F." }
}

pub(super) fn encode(file: &mut StepFile, solid: &Solid) -> Result<Entity, StepError> {
    let boundary = solid.boundary();
    let finite = |point: &Point| point.iter().all(|x| x.is_finite() && x.abs() <= 1e12);
    if !boundary.vertices.iter().all(finite) {
        return Err(StepError::InvalidPlacement);
    }
    let vertices: Vec<Entity> = boundary
        .vertices
        .iter()
        .map(|&point| file.vertex(point).vertex)
        .collect();
    let mut edges = Vec::with_capacity(boundary.edges.len());
    for edge in &boundary.edges {
        let geometry = curve(file, &edge.curve);
        let [start, end] = edge.vertices.map(|vertex| vertices[vertex]);
        edges.push(file.entity(format_args!("EDGE_CURVE('',{start},{end},{geometry},.T.)")));
    }
    let mut shells = Vec::with_capacity(boundary.shells.len());
    for shell in &boundary.shells {
        let mut faces = Vec::with_capacity(shell.faces.len());
        for &index in &shell.faces {
            let face = &boundary.faces[index];
            let mut bounds = Vec::with_capacity(face.loops.len());
            for Loop { uses, vertex } in &face.loops {
                let trim = match vertex {
                    Some(vertex) => {
                        file.entity(format_args!("VERTEX_LOOP('',{})", vertices[*vertex]))
                    }
                    None => {
                        let oriented: Vec<Entity> = uses
                            .iter()
                            .map(|&(edge, forward)| file.oriented_edge(edges[edge], forward))
                            .collect();
                        file.entity(format_args!("EDGE_LOOP('',{})", Entities(&oriented)))
                    }
                };
                bounds.push(file.entity(format_args!("FACE_BOUND('',{trim},.T.)")));
            }
            let support = surface(file, &face.surface);
            faces.push(file.face_on(&bounds, support, face.same_sense));
        }
        let closed = file.entity(format_args!("CLOSED_SHELL('',{})", Entities(&faces)));
        shells.push((closed, shell.reversed));
    }
    let outer = shells[0].0;
    if shells.len() == 1 {
        return Ok(file.entity(format_args!("MANIFOLD_SOLID_BREP('',{outer})")));
    }
    let voids: Vec<Entity> = shells[1..]
        .iter()
        .map(|&(shell, reversed)| {
            file.entity(format_args!(
                "ORIENTED_CLOSED_SHELL('',*,{shell},{})",
                sense(!reversed)
            ))
        })
        .collect();
    Ok(file.entity(format_args!(
        "BREP_WITH_VOIDS('',{outer},{})",
        Entities(&voids)
    )))
}
