//! Shared actual outer head geometry, independent of its recess family.
use super::*;
use crate::{head_profile::Profile, phillips_drive::PhillipsScrew};
#[derive(Clone, Copy)]
pub(super) struct Data<'a> {
    pub profile: &'a Profile,
    pub thread: &'a Thread,
    pub plain: f64,
    pub visible: [f64; 2],
}
impl<'a> From<&'a Screw> for Data<'a> {
    fn from(s: &'a Screw) -> Self {
        Self {
            profile: s.profile(),
            thread: s.thread(),
            plain: s.spec().shank_length_mm,
            visible: s.visible_thread_interval_mm(),
        }
    }
}
impl<'a> From<&'a PhillipsScrew> for Data<'a> {
    fn from(s: &'a PhillipsScrew) -> Self {
        Self {
            profile: s.profile(),
            thread: s.thread(),
            plain: s.spec().shank_length_mm,
            visible: s.visible_thread_interval_mm(),
        }
    }
}
impl Data<'_> {
    pub fn top_z(self) -> f64 {
        self.profile.top()[1] - self.plain
    }
}
pub(super) struct Shell {
    pub top_edge: Entity,
    pub top_center: Entity,
    pub walls: Vec<Entity>,
}
pub(super) fn admission(data: Data<'_>, origin: Point, tolerance: f64) -> Result<(), StepError> {
    if !tolerance.is_finite() || tolerance <= 0. || !origin.into_iter().all(f64::is_finite) {
        return Err(StepError::InvalidGeometry);
    }
    let tolerance = tolerance.min(MAX_ERROR_MM);
    let step = data.profile.end()[0] - data.thread.spec().major_radius_mm;
    if step > 0. && step <= 10. * tolerance {
        return Err(StepError::Limit);
    }
    if step > 0.
        && data.profile.segments().count() == 1
        && data.profile.top()[0] == data.profile.end()[0]
    {
        return Err(StepError::Unsupported);
    }
    let mut scale =
        origin.into_iter().map(f64::abs).fold(0., f64::max) + data.thread.spec().length_mm + 1.;
    for segment in data.profile.segments() {
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
            let a = (segment.start[1] - center[1]).portable_atan2(segment.start[0] - center[0]);
            let b = a + sweep;
            let error = (center[0] + r * b.portable_cos() - segment.end[0])
                .portable_hypot(center[1] + r * b.portable_sin() - segment.end[1]);
            if error + data.profile.rounding_error_mm() > tolerance / 8. {
                return Err(StepError::Limit);
            }
        }
    }
    if 128. * f64::EPSILON * scale + 4. * data.profile.rounding_error_mm() > tolerance / 4. {
        return Err(StepError::Limit);
    }
    Ok(())
}
pub(super) fn placement(data: Data<'_>, p: &Placement, tolerance: f64) -> Result<(), StepError> {
    let norm = |v: Point| v[0].portable_hypot(v[1]).portable_hypot(v[2]);
    let n = p.normal.map(|v| v / norm(p.normal));
    let projection: f64 = p.u.into_iter().zip(n).map(|(a, b)| a * b).sum();
    let u: Point = std::array::from_fn(|i| p.u[i] - projection * n[i]);
    let u = u.map(|v| v / norm(u));
    let v = cross(n, u);
    let error = |a: Point, b: Point| norm(std::array::from_fn(|i| a[i] - b[i]));
    let mut radial = data.thread.spec().major_radius_mm;
    let mut axial = data.thread.spec().length_mm;
    for s in data.profile.segments() {
        radial = radial.max(s.start[0]).max(s.end[0]);
        axial = axial
            .max((s.start[1] - data.plain).abs())
            .max((s.end[1] - data.plain).abs());
        if let Curve::Arc {
            center,
            radius_mm: r,
            ..
        } = s.curve
        {
            radial = radial.max(center[0] + r);
            axial = axial.max((center[1] - data.plain).abs() + r);
        }
    }
    let residual = (error(p.u, u) + error(p.v, v)) * radial + error(p.normal, n) * axial;
    let scale = p.origin.into_iter().map(f64::abs).fold(0., f64::max) + radial + axial + 1.;
    if !residual.is_finite()
        || residual > 4. * data.profile.rounding_error_mm()
        || residual + 128. * f64::EPSILON * scale + 4. * data.profile.rounding_error_mm()
            > tolerance.min(MAX_ERROR_MM) / 4.
    {
        return Err(StepError::InvalidPlacement);
    }
    Ok(())
}
pub(super) fn shell(
    file: &mut StepFile,
    p: &Placement,
    data: Data<'_>,
    join: Junction,
) -> Result<Shell, StepError> {
    let segments: Vec<_> = data.profile.segments().collect();
    let zshift = data.plain;
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
    let first = rim(file, data.profile.top())?;
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
                center: file.point(p.point_f64([0.; 2], data.visible[0])),
                axis: None,
                rim: None,
                vertex: join.vertex,
                radius: data.thread.spec().major_radius_mm,
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
                            let angle = dr.abs().portable_atan2(dz);
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
                        let a = (segment.start[1] - center[1])
                            .portable_atan2(segment.start[0] - center[0]);
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
    Ok(Shell {
        top_edge,
        top_center,
        walls,
    })
}
#[cfg(test)]
mod tests;
