//! Exact shared-curve exchange of spline cap additions and recesses.
use super::cubic_trim::CubicTrim;
use super::*;
use crate::normal_relief::precision::{add_up, dot_error, mul_up, product_error};
use crate::{
    normal_sweep::Edge as AnalyticEdge,
    spline_cap_relief::{Curve, Relief, Surface},
    spline_sweep::Edge,
};
fn curve_reversed(edge: &Edge) -> bool {
    matches!(edge,Edge::Cubic(c) if c.reversed())
}
pub fn write_spline_cap_relief(
    relief: &Relief,
    origin: Point,
    u: Point,
    v: Point,
) -> Result<String, StepError> {
    write_one(|writer| writer.append_spline_cap_relief(relief, origin, u, v))
}
pub(super) fn encode(
    file: &mut StepFile,
    relief: &Relief,
    placement: Placement,
) -> Result<Entity, StepError> {
    precision(relief, &placement)?;
    let mut candidate = StepFile::new();
    let solid = encode_shell(&mut candidate, relief, &placement)?;
    candidate.finish(vec![solid]);
    crate::step_read::verify_spline_cap_relief(
        &candidate.data,
        relief,
        placement.origin,
        placement.u,
        placement.v,
        relief.tolerance_mm(),
    )
    .map_err(|_| StepError::Unsupported)?;
    crate::step_read::read_spline_cap_relief(&candidate.data, relief.tolerance_mm())
        .map_err(|_| StepError::Unsupported)?;
    encode_shell(file, relief, &placement)
}
fn encode_shell(
    file: &mut StepFile,
    relief: &Relief,
    placement: &Placement,
) -> Result<Entity, StepError> {
    let graph = relief.boundary();
    let placed = |p: Point| placement.point_f64([p[0], p[1]], p[2]);
    let vertices = graph
        .vertices
        .iter()
        .map(|&point| file.vertex(placed(point)))
        .collect::<Vec<_>>();
    let mut supports = std::collections::BTreeMap::new();
    let mut edges = Vec::with_capacity(graph.edges.len());
    for edge in &graph.edges {
        let ends = edge.vertices.map(|i| vertices[i]);
        let entity = match &edge.curve {
            Curve::Line { .. }
            | Curve::Profile {
                edge: Edge::Analytic(AnalyticEdge::Line { .. }),
                ..
            } => line_edge(file, ends[0], ends[1])?,
            Curve::Profile {
                edge:
                    Edge::Analytic(AnalyticEdge::Circle {
                        center,
                        radius_mm,
                        clockwise,
                    }),
                z_mm,
            } => {
                let center = placement.point_f64(*center, *z_mm);
                let seam = std::array::from_fn(|i| center[i] + radius_mm * placement.u[i]);
                circular_seam(
                    center,
                    seam,
                    ends[0].point,
                    edge.vertices[0] == edge.vertices[1],
                    relief.tolerance_mm() / 8.,
                )?;
                let center = file.point(center);
                let axis = file.axis(
                    center,
                    placement.normal.map(|v| if *clockwise { -v } else { v }),
                    placement.u,
                );
                let circle = file.entity(format_args!("CIRCLE('',{axis},{radius_mm:.16E})"));
                let vertex = ends[0].vertex;
                file.entity(format_args!(
                    "EDGE_CURVE('',{vertex},{vertex},{circle},.T.)"
                ))
            }
            Curve::Profile {
                edge: Edge::Analytic(edge @ AnalyticEdge::Arc { .. }),
                z_mm,
            } => encode_arc(file, ends, placement, *edge, *z_mm)?,
            Curve::Profile {
                edge: Edge::Rational(_),
                ..
            } => return Err(StepError::Unsupported),
            Curve::Profile {
                edge: Edge::Cubic(curve),
                z_mm,
            } => {
                let (entity, support) = CubicTrim::decoded(curve, relief.tolerance_mm() / 4.)
                    .encode(file, ends, placement, *z_mm)?;
                supports.insert((edge.vertices[0], edge.vertices[1]), support);
                entity
            }
        };
        let forward =
            !matches!(&edge.curve, Curve::Profile {edge:Edge::Cubic(curve),..} if curve.reversed());
        edges.push((entity, forward));
    }
    let mut faces = Vec::with_capacity(graph.faces.len());
    for face in &graph.faces {
        let bounds = face
            .loops
            .iter()
            .enumerate()
            .map(|(i, loop_)| {
                let edges = loop_
                    .iter()
                    .map(|edge| {
                        file.oriented_edge(edges[edge.edge].0, edge.forward == edges[edge.edge].1)
                    })
                    .collect::<Vec<_>>();
                file.bound(&edges, i == 0)
            })
            .collect::<Vec<_>>();
        let entity = match &face.surface {
            Surface::Plane { z_mm, upward } => {
                let origin = file.point(placed([0., 0., *z_mm]));
                file.face(
                    &bounds,
                    origin,
                    placement.normal.map(|n| if *upward { n } else { -n }),
                    placement.u,
                )
            }
            Surface::Extrusion {
                edge,
                span_mm,
                inward,
            } => match edge {
                Edge::Analytic(AnalyticEdge::Line { a, b }) => {
                    let tangent = wall_tangent(placement, *a, *b, span_mm[0])?;
                    let origin = file.point(placement.point_f64(*a, span_mm[0]));
                    let normal =
                        cross(tangent, placement.normal).map(|n| if *inward { -n } else { n });
                    file.face(&bounds, origin, normal, tangent)
                }
                Edge::Analytic(
                    edge @ (AnalyticEdge::Arc {
                        center, clockwise, ..
                    }
                    | AnalyticEdge::Circle {
                        center, clockwise, ..
                    }),
                ) => {
                    let radius = edge.radius_mm().expect("circular support");
                    let a = edge.start();
                    let origin = file.point(placement.point_f64(*center, span_mm[0]));
                    let axis = file.axis(
                        origin,
                        placement.normal,
                        placement.vector([
                            (a[0] - center[0]) / radius,
                            (a[1] - center[1]) / radius,
                            0.,
                        ]),
                    );
                    let surface =
                        file.entity(format_args!("CYLINDRICAL_SURFACE('',{axis},{radius:.16E})"));
                    file.face_on(&bounds, surface, (!*clockwise) != *inward)
                }
                Edge::Cubic(_) => {
                    let use_=face.loops[0].iter().find(|use_|matches!(&graph.edges[use_.edge].curve,Curve::Profile{z_mm,..} if *z_mm==span_mm[0])).ok_or(StepError::InvalidGeometry)?;
                    let ends = graph.edges[use_.edge].vertices;
                    let curve = supports[&(ends[0], ends[1])];
                    let axis = file.direction(placement.normal);
                    let vector = file.entity(format_args!("VECTOR('',{axis},1.)"));
                    let surface = file.entity(format_args!(
                        "SURFACE_OF_LINEAR_EXTRUSION('',{curve},{vector})"
                    ));
                    file.face_on(&bounds, surface, (!curve_reversed(edge)) != *inward)
                }
                Edge::Rational(_) => return Err(StepError::Unsupported),
            },
        };
        faces.push(entity);
    }
    Ok(file.solid(&faces))
}

fn circular_seam(
    center: Point,
    seam: Point,
    vertex: Point,
    closed: bool,
    tolerance: f64,
) -> Result<(), StepError> {
    if !center
        .into_iter()
        .chain(seam)
        .chain(vertex)
        .all(f64::is_finite)
        || !tolerance.is_finite()
        || tolerance < 0.
        || center == seam
        || !closed
        || length3(subtract3(vertex, seam)) > tolerance
    {
        return Err(StepError::InvalidGeometry);
    }
    Ok(())
}
fn subtract3(a: Point, b: Point) -> Point {
    std::array::from_fn(|i| a[i] - b[i])
}
fn length3(a: Point) -> f64 {
    a[0].portable_hypot(a[1]).portable_hypot(a[2])
}
fn wall_tangent(
    placement: &Placement,
    a: [f64; 2],
    b: [f64; 2],
    z: f64,
) -> Result<Point, StepError> {
    direction(placement.point_f64(a, z), placement.point_f64(b, z))
}
fn encode_arc(
    file: &mut StepFile,
    ends: [Vertex; 2],
    placement: &Placement,
    edge: AnalyticEdge,
    z: f64,
) -> Result<Entity, StepError> {
    let AnalyticEdge::Arc {
        center,
        a,
        clockwise,
        ..
    } = edge
    else {
        return Err(StepError::InvalidGeometry);
    };
    let radius = edge.radius_mm().expect("arc radius");
    if !radius.is_finite() || radius <= 0. || edge.signed_sweep_rad().expect("arc sweep") == 0. {
        return Err(StepError::InvalidGeometry);
    }
    let radial = [(a[0] - center[0]) / radius, (a[1] - center[1]) / radius, 0.];
    circular_arc_edge(
        file,
        ends,
        placement.point_f64(center, z),
        placement.normal.map(|v| if clockwise { -v } else { v }),
        placement.vector(radial),
        radius,
        edge.signed_sweep_rad().expect("arc sweep"),
    )
    .map(|(edge, _)| edge)
}

fn evaluation_error(placement: &Placement, p: Point) -> f64 {
    let axes = [placement.u, placement.v, placement.normal];
    let mut total = 0.;
    for i in 0..3 {
        let mut value = placement.origin[i];
        let mut error = 0.;
        for (axis, coordinate) in axes.iter().zip(p) {
            let term = axis[i] * coordinate;
            error = add_up(error, product_error(axis[i], coordinate, term));
            let sum = value + term;
            let part = sum - value;
            error = add_up(error, ((value - (sum - part)) + (term - part)).abs());
            value = sum;
        }
        total = add_up(total, error);
    }
    total
}
fn precision(relief: &Relief, placement: &Placement) -> Result<(), StepError> {
    let axes = [placement.u, placement.v, placement.normal];
    let mut gram = 0.;
    for i in 0..3 {
        for j in 0..3 {
            let error = dot_error(axes[i], axes[j], if i == j { 1. } else { 0. });
            gram = add_up(gram, mul_up(error, error));
        }
    }
    let bounds = relief.bounds();
    let mut extent = 0.;
    for (low, high) in bounds[0].into_iter().zip(bounds[1]) {
        let e = low.abs().max(high.abs());
        extent = add_up(extent, mul_up(e, e));
    }
    let sqrt = |value: f64| {
        if value == 0. {
            0.
        } else {
            value.sqrt().next_up()
        }
    };
    if mul_up(sqrt(gram), sqrt(extent)) > relief.tolerance_mm() / 8. {
        return Err(StepError::Unsupported);
    }
    let graph = relief.boundary();
    let controls = graph
        .vertices
        .iter()
        .copied()
        .chain(graph.edges.iter().flat_map(|edge| match &edge.curve {
            Curve::Line { .. } => vec![],
            Curve::Profile {
                edge: Edge::Cubic(c),
                z_mm,
            } => c.poles().iter().map(|p| [p[0], p[1], *z_mm]).collect(),
            Curve::Profile {
                edge: Edge::Rational(c),
                z_mm,
            } => c.poles().iter().map(|p| [p[0], p[1], *z_mm]).collect(),
            Curve::Profile {
                edge:
                    Edge::Analytic(
                        AnalyticEdge::Arc { center, .. } | AnalyticEdge::Circle { center, .. },
                    ),
                z_mm,
            } => vec![[center[0], center[1], *z_mm]],
            Curve::Profile {
                edge: Edge::Analytic(AnalyticEdge::Line { .. }),
                ..
            } => vec![],
        }));
    if controls
        .into_iter()
        .any(|p| evaluation_error(placement, p) > relief.tolerance_mm() / 8.)
    {
        return Err(StepError::Unsupported);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        spline_cap_relief::{Kind, tests::fixture},
        spline_sweep::Sweep,
    };
    #[test]
    fn finite_periodic_trim_and_tangent_contracts_reject_mismatched_supports_and_retry() {
        assert_eq!(
            circular_seam([0.; 3], [0.; 3], [0.; 3], true, 1e-7),
            Err(StepError::InvalidGeometry)
        );
        assert_eq!(
            circular_seam([0.; 3], [1., 0., 0.], [1., 0., 0.], false, 1e-7),
            Err(StepError::InvalidGeometry)
        );
        assert_eq!(
            circular_seam([0.; 3], [1., 0., 0.], [2., 0., 0.], true, 1e-7),
            Err(StepError::InvalidGeometry)
        );
        assert_eq!(
            circular_seam([f64::NAN, 0., 0.], [1., 0., 0.], [1., 0., 0.], true, 1e-7),
            Err(StepError::InvalidGeometry)
        );
        assert_eq!(
            circular_seam([0.; 3], [1., 0., 0.], [1., 0., 0.], true, -1.),
            Err(StepError::InvalidGeometry)
        );
        circular_seam([0.; 3], [1., 0., 0.], [1., 0., 0.], true, 1e-7).unwrap();
        let frame = Placement::new([0.; 3], [1., 0., 0.], [0., 1., 0.]).unwrap();
        assert_eq!(
            wall_tangent(&frame, [0., 0.], [0., 0.], 0.),
            Err(StepError::InvalidGeometry)
        );
        assert_eq!(
            wall_tangent(&frame, [0., 0.], [1., 0.], 0.).unwrap(),
            [1., 0., 0.]
        );
        let mut file = StepFile::new();
        let a = file.vertex([1., 0., 0.]);
        let wrong = file.vertex([0., 2., 0.]);
        let arc = AnalyticEdge::Arc {
            a: [1., 0.],
            b: [0., 1.],
            center: [0.; 2],
            clockwise: false,
        };
        assert_eq!(
            encode_arc(&mut file, [a, wrong], &frame, arc, 0.).err(),
            Some(StepError::InvalidGeometry)
        );
        let invalid = AnalyticEdge::Arc {
            a: [0.; 2],
            b: [0., 1.],
            center: [0.; 2],
            clockwise: false,
        };
        assert_eq!(
            encode_arc(&mut file, [a, wrong], &frame, invalid, 0.).err(),
            Some(StepError::InvalidGeometry)
        );
        let b = file.vertex([0., 1., 0.]);
        encode_arc(&mut file, [a, b], &frame, arc, 0.).unwrap();
        assert_eq!(
            encode_arc(
                &mut file,
                [a, b],
                &frame,
                AnalyticEdge::Line {
                    a: [1., 0.],
                    b: [0., 1.]
                },
                0.
            )
            .err(),
            Some(StepError::InvalidGeometry)
        );
    }
    #[test]
    fn frame_extent_and_actual_placement_roundoff_refuse_before_writer_mutation() {
        let source = fixture(Kind::Raised);
        let square = |s: f64| {
            let p = [[0., 0.], [s, 0.], [s, s], [0., s]];
            (0..4)
                .map(|i| {
                    Edge::Analytic(AnalyticEdge::Line {
                        a: p[i],
                        b: p[(i + 1) % 4],
                    })
                })
                .collect()
        };
        let base = Sweep::new(square(2000.), None, 5., 1e-7).unwrap();
        let geometry =
            Relief::new(base, source.motifs().to_vec(), Kind::Raised, 1e-7, 1e-4).unwrap();
        let mut writer = Writer::new();
        assert_eq!(
            writer.append_spline_cap_relief(&geometry, [0.; 3], [1. + 4e-11, 0., 0.], [0., 1., 0.]),
            Err(StepError::Unsupported)
        );
        assert_eq!(
            writer.append_spline_cap_relief(&source, [1e20, 0., 0.], [1., 0., 0.], [0., 1., 0.]),
            Err(StepError::Unsupported)
        );
        writer
            .append_spline_cap_relief(&source, [0.; 3], [1., 0., 0.], [0., 1., 0.])
            .unwrap();
        let text = writer.finish().unwrap();
        assert!(crate::step_read::read_spline_cap_relief(&text, 1e-7).is_ok());
    }
}

#[cfg(test)]
mod placement_tests {
    use super::*;
    #[test]
    fn two_associations_of_a_far_periodic_seam_cannot_silently_emit_inconsistent_incidence() {
        let source =
            crate::spline_cap_relief::tests::fixture(crate::spline_cap_relief::Kind::Raised);
        let original = &source.motifs()[0];
        let motif = crate::spline_sweep::Sweep::new(
            original.outline().to_vec(),
            Some(AnalyticEdge::Circle {
                center: [6.1, 6.1],
                radius_mm: 0.27,
                clockwise: true,
            }),
            1.,
            1e-7,
        )
        .unwrap();
        let geometry = Relief::new(
            source.base().clone(),
            vec![motif],
            crate::spline_cap_relief::Kind::Raised,
            1e-7,
            1e-4,
        )
        .unwrap();
        let far = Placement::new([1e8, 0., 0.], [1., 0., 0.], [0., 1., 0.]).unwrap();
        assert_eq!(
            encode_shell(&mut StepFile::new(), &geometry, &far).err(),
            Some(StepError::InvalidGeometry)
        );
        let exact = Placement::new([0.; 3], [1., 0., 0.], [0., 1., 0.]).unwrap();
        let mut file = StepFile::new();
        let solid = encode_shell(&mut file, &geometry, &exact).unwrap();
        file.finish(vec![solid]);
        assert_eq!(
            crate::step_read::verify_spline_cap_relief(
                &file.data,
                &geometry,
                [0.; 3],
                [1., 0., 0.],
                [0., 1., 0.],
                1e-7
            )
            .unwrap(),
            geometry.volume_mm3()
        );
    }
}
