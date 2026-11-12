//! Exact finite surfaces for a retained head, sharing every thread start rim.
use super::*;
use crate::{
    head_profile::{Curve, Role},
    profiled_screw::Screw,
};
use std::f64::consts::TAU;
#[cfg(test)]
mod tests;

#[derive(Clone, Copy)]
struct Station {
    center: Entity,
    axis: Option<Entity>,
    rim: Option<Entity>,
    vertex: Vertex,
    radius: f64,
}
impl From<CircularRim> for Station {
    fn from(rim: CircularRim) -> Self {
        Self {
            center: rim.center,
            axis: Some(rim.axis),
            rim: Some(rim.edge),
            vertex: rim.vertex,
            radius: rim.radius,
        }
    }
}

pub(super) struct Junction {
    pub(super) crests: Vec<Entity>,
    pub(super) complements: Vec<Entity>,
    vertex: Vertex,
    reference: Point,
    positive: Vec<(Entity, bool)>,
    pub(super) outer: Option<CircularRim>,
}

pub(super) fn admission(screw: &Screw, origin: Point, tolerance: f64) -> Result<(), StepError> {
    let (drive, socket_depth) = screw
        .recess()
        .flat_hex_geometry()
        .ok_or(StepError::Unsupported)?;
    if !tolerance.is_finite() || tolerance <= 0. || !origin.into_iter().all(f64::is_finite) {
        return Err(StepError::InvalidGeometry);
    }
    let tolerance = tolerance.min(MAX_ERROR_MM);
    let step = screw.profile().end()[0] - screw.thread().spec().major_radius_mm;
    if step > 0. && step <= 10. * tolerance {
        return Err(StepError::Limit);
    }
    if step > 0.
        && screw.profile().segments().count() == 1
        && screw.profile().top()[0] == screw.profile().end()[0]
    {
        return Err(StepError::Unsupported);
    }
    let mut scale =
        origin.into_iter().map(f64::abs).fold(0., f64::max) + screw.thread().spec().length_mm + 1.;
    for segment in screw.profile().segments() {
        for p in [segment.start, segment.end] {
            scale = scale.max(p[0].abs() + p[1].abs());
        }
        let dr = (segment.end[0] - segment.start[0]).abs();
        let dz = segment.end[1] - segment.start[1];
        if [dr, dz]
            .into_iter()
            .any(|x| x != 0. && x <= 10. * tolerance)
            || segment.start[0].min(segment.end[0]) <= 10. * tolerance
        {
            return Err(StepError::Limit);
        }
        if let Curve::Arc {
            center,
            radius_mm: r,
            sweep,
        } = segment.curve
        {
            scale += center[0].abs() + center[1].abs() + r;
            if r <= 10. * tolerance {
                return Err(StepError::Limit);
            }
            let a = (segment.start[1] - center[1]).atan2(segment.start[0] - center[0]);
            let b = a + sweep;
            let error = (center[0] + r * b.cos() - segment.end[0])
                .hypot(center[1] + r * b.sin() - segment.end[1]);
            if error + screw.profile().rounding_error_mm() > tolerance / 8. {
                return Err(StepError::Limit);
            }
        }
    }
    if 128. * f64::EPSILON * scale + 4. * screw.profile().rounding_error_mm() > tolerance / 4.
        || socket_depth <= 10. * tolerance
        || drive.across_flats_mm() <= 10. * tolerance
        || !screw.profile().admits_socket(
            drive.across_flats_mm() / 3_f64.sqrt(),
            socket_depth,
            10. * tolerance,
        )
    {
        return Err(StepError::Limit);
    }
    Ok(())
}

pub(super) fn placement(screw: &Screw, p: &Placement, tolerance: f64) -> Result<(), StepError> {
    let norm = |v: Point| v[0].hypot(v[1]).hypot(v[2]);
    let n = p.normal.map(|v| v / norm(p.normal));
    let projection: f64 = p.u.into_iter().zip(n).map(|(a, b)| a * b).sum();
    let u: Point = std::array::from_fn(|i| p.u[i] - projection * n[i]);
    let u = u.map(|v| v / norm(u));
    let v = cross(n, u);
    let error = |a: Point, b: Point| norm(std::array::from_fn(|i| a[i] - b[i]));
    let mut radial = screw.thread().spec().major_radius_mm;
    let mut axial = screw.thread().spec().length_mm;
    for s in screw.profile().segments() {
        radial = radial.max(s.start[0]).max(s.end[0]);
        axial = axial
            .max((s.start[1] - screw.spec().shank_length_mm).abs())
            .max((s.end[1] - screw.spec().shank_length_mm).abs());
        if let Curve::Arc {
            center,
            radius_mm: r,
            ..
        } = s.curve
        {
            radial = radial.max(center[0] + r);
            axial = axial.max((center[1] - screw.spec().shank_length_mm).abs() + r);
        }
    }
    let residual = (error(p.u, u) + error(p.v, v)) * radial + error(p.normal, n) * axial;
    let scale = p.origin.into_iter().map(f64::abs).fold(0., f64::max) + radial + axial + 1.;
    if !residual.is_finite()
        || residual > 4. * screw.profile().rounding_error_mm()
        || residual + 128. * f64::EPSILON * scale + 4. * screw.profile().rounding_error_mm()
            > tolerance.min(MAX_ERROR_MM) / 4.
    {
        return Err(StepError::InvalidPlacement);
    }
    Ok(())
}

pub(super) fn junction(
    file: &mut StepFile,
    p: &Placement,
    thread: &Thread,
    shank_radius: f64,
    z: f64,
    corners: &[[Vertex; 2]],
) -> Result<Junction, StepError> {
    let spec = thread.spec();
    let [_, start, _, _, crest, _] = thread.profile_breaks();
    let rate = thread.helical_rate();
    let radial = |phase: f64| {
        let (sin, cos) = ((z - phase) / rate).sin_cos();
        p.vector([cos, sin, 0.])
    };
    let center = p.point_f64([0.; 2], z);
    let axis = p.normal.map(|x| if spec.left_hand { x } else { -x });
    let mut crests = Vec::new();
    let mut complements = Vec::new();
    for (k, &[first, last]) in corners.iter().enumerate() {
        let shift = k as f64 * spec.pitch_mm;
        let next = corners[(k + 1) % corners.len()][0];
        crests.push(
            circular_arc_edge(
                file,
                [last, next],
                center,
                axis,
                radial(crest + shift),
                spec.major_radius_mm,
                (start + spec.pitch_mm - crest) / rate.abs(),
            )?
            .0,
        );
        if shank_radius == spec.major_radius_mm {
            complements.push(
                circular_arc_edge(
                    file,
                    [first, last],
                    center,
                    axis,
                    radial(start + shift),
                    spec.major_radius_mm,
                    (crest - start) / rate.abs(),
                )?
                .0,
            );
        }
    }
    let mut positive: Vec<_> = complements
        .iter()
        .zip(&crests)
        .flat_map(|(&a, &b)| [(a, true), (b, true)])
        .collect();
    if !spec.left_hand {
        positive = positive.into_iter().rev().map(|(e, s)| (e, !s)).collect();
    }
    let outer = if shank_radius > spec.major_radius_mm {
        Some(CircularRim::new(
            file,
            center,
            p.normal,
            radial(start),
            shank_radius,
        )?)
    } else {
        None
    };
    if let Some(outer) = &outer {
        positive = vec![(outer.edge, true)];
    }
    Ok(Junction {
        crests,
        complements,
        vertex: corners[0][0],
        reference: radial(start),
        positive,
        outer,
    })
}

pub(super) fn head(
    file: &mut StepFile,
    p: &Placement,
    screw: &Screw,
    join: Junction,
) -> Result<Vec<Entity>, StepError> {
    let (drive, socket_depth) = screw
        .recess()
        .flat_hex_geometry()
        .ok_or(StepError::Unsupported)?;
    let segments: Vec<_> = screw.profile().segments().collect();
    let zshift = screw.spec().shank_length_mm;
    let rim = |file: &mut StepFile, [r, z]: [f64; 2]| {
        CircularRim::new(
            file,
            p.point_f64([0.; 2], z - zshift),
            p.normal,
            join.reference,
            r,
        )
        .map(Station::from)
    };
    let first = rim(file, screw.profile().top())?;
    let top_edge = first.rim.expect("initial full circle");
    let top_center = first.center;
    let terminal = join.outer.map(Station::from);
    let mut previous = first;
    let mut walls = Vec::new();
    for (index, segment) in segments.iter().enumerate() {
        let last = index + 1 == segments.len();
        let next = if last {
            // The nominal circle is split into actual crest/complement edges.
            terminal.unwrap_or(Station {
                center: file.point(p.point_f64([0.; 2], screw.visible_thread_interval_mm()[0])),
                axis: None,
                rim: None,
                vertex: join.vertex,
                radius: screw.thread().spec().major_radius_mm,
            })
        } else {
            rim(file, segment.end)?
        };
        let boundary: Vec<_> = if last {
            join.positive.clone()
        } else {
            vec![(next.rim.expect("nonterminal circle"), true)]
        };
        let face = match segment.curve {
            Curve::Line if segment.start[1] == segment.end[1] => {
                let outer =
                    circular_bound(file, previous.rim.expect("nonterminal circle"), true, true);
                let edges: Vec<_> = boundary
                    .into_iter()
                    .rev()
                    .map(|(edge, sense)| file.oriented_edge(edge, !sense))
                    .collect();
                let inner = file.bound(&edges, false);
                file.face(&[outer, inner], previous.center, p.normal, p.u)
            }
            curve => {
                let (surface, seam, forward) = match curve {
                    Curve::Line => {
                        let dr = next.radius - previous.radius;
                        let dz = segment.end[1] - segment.start[1];
                        let surface = if dr == 0. {
                            let axis = previous.axis.expect("nonterminal circle axis");
                            let r = previous.radius;
                            file.entity(format_args!("CYLINDRICAL_SURFACE('',{axis},{r:.16E})"))
                        } else {
                            let (axis, r) = if dr > 0. {
                                (
                                    previous.axis.expect("nonterminal circle axis"),
                                    previous.radius,
                                )
                            } else {
                                (
                                    file.axis(next.center, p.normal.map(|x| -x), join.reference),
                                    next.radius,
                                )
                            };
                            let angle = dr.abs().atan2(dz);
                            file.entity(format_args!(
                                "CONICAL_SURFACE('',{axis},{r:.16E},{angle:.16E})"
                            ))
                        };
                        (
                            surface,
                            line_edge(file, previous.vertex, next.vertex)?,
                            true,
                        )
                    }
                    Curve::Arc {
                        center,
                        radius_mm: r,
                        sweep,
                    } => {
                        let center_z = center[1] - zshift;
                        let center_id = file.point(p.point_f64([0.; 2], center_z));
                        let axis = file.axis(center_id, p.normal, join.reference);
                        let surface = if center[0] == 0. {
                            file.entity(format_args!("SPHERICAL_SURFACE('',{axis},{r:.16E})"))
                        } else {
                            let major = center[0];
                            file.entity(format_args!(
                                "TOROIDAL_SURFACE('',{axis},{major:.16E},{r:.16E})"
                            ))
                        };
                        let radial = join.reference;
                        let center3 = std::array::from_fn(|i| {
                            p.origin[i] + center[0] * radial[i] + center_z * p.normal[i]
                        });
                        let ref3 =
                            std::array::from_fn(|i| (previous.vertex.point[i] - center3[i]) / r);
                        let normal = cross(p.normal, radial).map(|x| -sweep.signum() * x);
                        let circle_center = file.point(center3);
                        let circle_axis = file.axis(circle_center, normal, ref3);
                        let circle = file.entity(format_args!("CIRCLE('',{circle_axis},{r:.16E})"));
                        let a = (segment.start[1] - center[1]).atan2(segment.start[0] - center[0]);
                        let seam = seam(
                            file,
                            surface,
                            circle,
                            [previous.vertex.vertex, next.vertex.vertex],
                            a,
                            sweep.signum(),
                        );
                        (surface, seam, segment.role == Role::Crown)
                    }
                };
                let mut edges = vec![
                    file.oriented_edge(previous.rim.expect("nonterminal circle"), true),
                    file.oriented_edge(seam, true),
                ];
                edges.extend(
                    boundary
                        .into_iter()
                        .rev()
                        .map(|(edge, sense)| file.oriented_edge(edge, !sense)),
                );
                edges.push(file.oriented_edge(seam, false));
                let bound = file.bound(&edges, true);
                file.face_on(&[bound], surface, forward)
            }
        };
        walls.push(face);
        previous = next;
    }
    let top = screw.head_top_z();
    let ([mouth, floor], sides) =
        encode_drive_prism(file, p, drive.corners_mm(), [top, top + socket_depth], true)?;
    let outer = circular_bound(file, top_edge, false, true);
    let down = p.normal.map(|x| -x);
    let mut faces = vec![file.face(&[outer, mouth], top_center, down, p.u)];
    faces.extend(walls);
    let floor_center = file.point(p.point_f64([0.; 2], top + socket_depth));
    faces.push(file.face(&[floor], floor_center, down, p.u));
    faces.extend(sides);
    Ok(faces)
}

fn seam(
    file: &mut StepFile,
    surface: Entity,
    circle: Entity,
    [start, end]: [Entity; 2],
    v: f64,
    sign: f64,
) -> Entity {
    let context=file.entity(format_args!("(GEOMETRIC_REPRESENTATION_CONTEXT(2)PARAMETRIC_REPRESENTATION_CONTEXT()REPRESENTATION_CONTEXT('2D SPACE',''))"));
    let curves = [0., TAU].map(|u| {
        let point = file.entity(format_args!("CARTESIAN_POINT('',({u:.16E},{v:.16E}))"));
        let direction = file.entity(format_args!("DIRECTION('',(0.,{sign:.16E}))"));
        let vector = file.entity(format_args!("VECTOR('',{direction},1.)"));
        let line = file.entity(format_args!("LINE('',{point},{vector})"));
        let representation = file.entity(format_args!(
            "DEFINITIONAL_REPRESENTATION('',({line}),{context})"
        ));
        file.entity(format_args!("PCURVE('',{surface},{representation})"))
    });
    let [a, b] = curves;
    let seam = file.entity(format_args!(
        "SEAM_CURVE('',{circle},({a},{b}),.PCURVE_S1.)"
    ));
    file.entity(format_args!("EDGE_CURVE('',{start},{end},{seam},.T.)"))
}
