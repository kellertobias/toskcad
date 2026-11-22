//! Partial sweeps of rectilinear meridians: two planar section caps, and one
//! planar or cylindrical wall sector per meridian edge between them.
//!
//! Every meridian vertex sweeps an exact circular arc from the start cap to
//! the end cap. Radial edges sweep annular plane sectors and axial edges
//! cylinder sectors. Holes in the meridian become holes in both caps and
//! inner walls, so no further shell is needed.

use super::*;
use crate::revolution::Revolution;

/// Whether a partial sweep has only axis-parallel and radial meridian edges
/// and keeps every vertex off the axis.
pub(super) fn supports(revolution: &Revolution) -> bool {
    !revolution.is_full()
        && revolution.rings().iter().all(|ring| {
            ring.iter().all(|point| point[0] > 0.)
                && ring
                    .iter()
                    .zip(ring.iter().cycle().skip(1))
                    .all(|(a, b)| (a[0] == b[0]) != (a[1] == b[1]))
        })
}

pub(super) fn encode(
    file: &mut StepFile,
    placement: Placement,
    revolution: &Revolution,
) -> Result<Entity, StepError> {
    if !supports(revolution) {
        return Err(StepError::Unsupported);
    }
    let angle = revolution.angle_rad();
    let turn = angle.signum();
    let (sin, cos) = angle.portable_sin_cos();
    let axis = placement.normal;
    let [u, v] = [placement.u, placement.v];
    // The radial direction of the end section, and the sweep direction there.
    let end_radial: Point = std::array::from_fn(|i| cos * u[i] + sin * v[i]);
    let end_sweep: Point = std::array::from_fn(|i| -sin * u[i] + cos * v[i]);
    let at = |[radius, level]: [f64; 2], end: bool| -> Point {
        let radial = if end { end_radial } else { u };
        std::array::from_fn(|i| placement.origin[i] + radius * radial[i] + level * axis[i])
    };
    // The sweep follows the right-hand rule about the axis times its sign.
    let sweep_axis = axis.map(|n| turn * n);
    let mut faces = Vec::new();
    let mut caps: [Vec<Entity>; 2] = [Vec::new(), Vec::new()];
    for (index, ring) in revolution.rings().iter().enumerate() {
        let starts: Vec<Vertex> = ring.iter().map(|&p| file.vertex(at(p, false))).collect();
        let ends: Vec<Vertex> = ring.iter().map(|&p| file.vertex(at(p, true))).collect();
        let arcs = ring
            .iter()
            .enumerate()
            .map(|(k, &[radius, level])| {
                let center = file.point(placement.point_f64([0.; 2], level));
                let frame = file.axis(center, sweep_axis, u);
                let circle = file.entity(format_args!("CIRCLE('',{frame},{radius:.16E})"));
                file.entity(format_args!(
                    "EDGE_CURVE('',{},{},{circle},.T.)",
                    starts[k].vertex, ends[k].vertex
                ))
            })
            .collect::<Vec<_>>();
        let count = ring.len();
        let mut start_lines = Vec::with_capacity(count);
        let mut end_lines = Vec::with_capacity(count);
        for k in 0..count {
            let next = (k + 1) % count;
            start_lines.push(line_edge(file, starts[k], starts[next])?);
            end_lines.push(line_edge(file, ends[k], ends[next])?);
        }
        for k in 0..count {
            let next = (k + 1) % count;
            let [a, b] = [ring[k], ring[next]];
            // Material lies left of every meridian edge, so the wall's
            // outward normal is the right normal (dz, -dr) of the edge.
            let outward = [b[1] - a[1], a[0] - b[0]];
            // Start line, arc at b, end line back, arc at a back: this loop
            // turns clockwise about the outward normal for a positive sweep.
            let mut edges = vec![
                file.oriented_edge(start_lines[k], true),
                file.oriented_edge(arcs[next], true),
                file.oriented_edge(end_lines[k], false),
                file.oriented_edge(arcs[k], false),
            ];
            if turn > 0. {
                edges = vec![
                    file.oriented_edge(arcs[k], true),
                    file.oriented_edge(end_lines[k], true),
                    file.oriented_edge(arcs[next], false),
                    file.oriented_edge(start_lines[k], false),
                ];
            }
            let bound = file.bound(&edges, true);
            if a[0] == b[0] {
                let center = file.point(placement.point_f64([0.; 2], a[1]));
                let frame = file.axis(center, axis, u);
                let radius = a[0];
                let surface = file.entity(format_args!(
                    "CYLINDRICAL_SURFACE('',{frame},{radius:.16E})"
                ));
                faces.push(file.face_on(&[bound], surface, outward[0] > 0.));
            } else {
                let center = file.point(placement.point_f64([0.; 2], a[1]));
                let normal = axis.map(|n| n * outward[1].signum());
                faces.push(file.face(&[bound], center, normal, u));
            }
        }
        // A ring counterclockwise in (radius, axial) runs counterclockwise
        // about u x axis, the start section's inward side.
        let outer = index == 0;
        for (end, lines) in [(false, &start_lines), (true, &end_lines)] {
            let forward = (turn > 0.) != end;
            let mut edges: Vec<_> = lines
                .iter()
                .map(|&line| file.oriented_edge(line, forward))
                .collect();
            if !forward {
                edges.reverse();
            }
            caps[usize::from(end)].push(file.bound(&edges, outer));
        }
    }
    let first = revolution.rings()[0][0];
    for (end, bounds) in caps.iter().enumerate() {
        let end = end == 1;
        // The start cap faces against the sweep, the end cap along it.
        let normal: Point = if end {
            end_sweep.map(|x| turn * x)
        } else {
            v.map(|x| -turn * x)
        };
        let reference = if end { end_radial } else { u };
        let point = file.point(at(first, end));
        faces.push(file.face(bounds, point, normal, reference));
    }
    Ok(file.solid(&faces))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, PI, TAU};

    const ORIGIN: Point = [13., -7., 20.];
    const U: Point = [1., 0., 0.];
    const V: Point = [0., 0.6, 0.8];

    fn read(text: &str) -> crate::brep::Solid {
        let [solid] = crate::step_read::brep::read_solids(text, 1e-7, None)
            .unwrap()
            .try_into()
            .unwrap();
        solid
    }

    #[test]
    fn partial_sweeps_of_holed_and_stepped_meridians_read_back_exactly() {
        let holed = vec![
            vec![[2., 0.], [5., 0.], [5., 6.], [2., 6.]],
            vec![[3., 2.], [3., 4.], [4., 4.], [4., 2.]],
        ];
        let stepped = vec![vec![
            [1., 0.],
            [6., 0.],
            [6., 2.],
            [3., 2.],
            [3., 5.],
            [1., 5.],
        ]];
        for rings in [holed, stepped] {
            for angle in [PI, -PI, FRAC_PI_2, 1.25 * PI, -0.3] {
                let revolution = Revolution::new(rings.clone(), angle).unwrap();
                assert!(supports_cylindrical_revolution(&revolution));
                let text = write_cylindrical_revolution(&revolution, ORIGIN, U, V).unwrap();
                let edges: usize = rings.iter().map(Vec::len).sum();
                assert_eq!(text.matches("=ADVANCED_FACE(").count(), edges + 2);
                assert_eq!(text.matches("=CIRCLE(").count(), edges);
                let solid = read(&text);
                assert_eq!(solid.faces().len(), edges + 2);
                let expected = revolution.volume_mm3();
                assert!(
                    (solid.volume_mm3() - expected).abs() <= 1e-7 + solid.volume_error_mm3(),
                    "{angle} {} {expected}",
                    solid.volume_mm3()
                );
                let placement = Placement::new(ORIGIN, U, V).unwrap();
                for axis in 0..3 {
                    let direction = [placement.u[axis], placement.v[axis], placement.normal[axis]];
                    let [low, high] =
                        [-1., 1.].map(|sign| revolution.support(direction.map(|x| sign * x)));
                    let bounds = solid.bounds_mm();
                    assert!((bounds[0][axis] - (ORIGIN[axis] - low)).abs() < 1e-6);
                    assert!((bounds[1][axis] - (ORIGIN[axis] + high)).abs() < 1e-6);
                }
            }
        }
    }

    #[test]
    fn slopes_axis_vertices_and_full_sweeps_are_not_partial_rectilinear() {
        for (ring, angle) in [
            (vec![[2., 0.], [5., 0.], [4., 6.], [2., 6.]], PI),
            (vec![[0., 0.], [5., 0.], [5., 6.], [0., 6.]], PI),
            (vec![[2., 0.], [5., 0.], [5., 6.], [2., 6.]], TAU),
        ] {
            let revolution = Revolution::new(vec![ring], angle).unwrap();
            assert!(!supports(&revolution));
        }
    }
}
