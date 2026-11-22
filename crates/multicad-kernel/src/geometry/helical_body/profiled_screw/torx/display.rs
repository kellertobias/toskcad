//! Certified polynomial-to-chord display of the true finite cone intersections.
use super::*;
pub(in crate::geometry::helical_body) fn edges(
    body: &HelicalBody,
    screw: &Screw,
    frame: PlaneFrame,
    deflection: f64,
    remaining: &mut usize,
) -> Result<Vec<DisplayEdge>, GeometryError> {
    let Some(source) = screw.recess().hexalobular_geometry() else {
        return Ok(Vec::new());
    };
    if !deflection.is_finite() || deflection <= 0. {
        return Err(input("Invalid curved drive display deflection"));
    }
    let axes = [frame.u, frame.v, frame.normal()];
    // Entrywise L1 bounds the operator norm of this actual stored placement.
    let scale = axes
        .into_iter()
        .flatten()
        .fold(0_f64, |sum, value| (sum + value.abs()).next_up());
    let request = (deflection / (4. * scale).next_up()).next_down();
    let topology = body.topology(frame, body.feature, 0)?;
    let mut edges = Vec::with_capacity(12);
    for i in 0..12 {
        if *remaining < 2 {
            return Err(error(
                OPERATION,
                GeometryErrorKind::Limit,
                "Curved drive edges exceed the display point budget",
            ));
        }
        let origin = body.origin(format!("screw.drive/rim/1/{i}"));
        let ordinal = topology
            .elements
            .iter()
            .find(|e| e.key.kind == TopologyKind::Edge && e.origins.contains(&origin))
            .expect("retained cone intersection edge")
            .key
            .ordinal;
        let binding = match screw.recess() {
            multicad_rust_kernel::screw_recess::Recess::ConicalHexalobular(recess) => {
                recess.finite_binding()
            }
            _ => unreachable!("curved recess source"),
        };
        let local = match binding {
            Some(binding) => binding.intersection_polyline(i, request, (*remaining).min(100_000)),
            None => source.intersection_polyline(i, request, (*remaining).min(100_000)),
        }
        .map_err(|cause| {
            error(
                OPERATION,
                GeometryErrorKind::Limit,
                format!("Unresolved curved drive display source: {cause:?}"),
            )
        })?;
        // Account for every actual coordinate and for cancellation in adding
        // the head offset; the placement reserve is independent of source shape.
        let local_extent = (local
            .iter()
            .flatten()
            .map(|v| v.abs())
            .fold(0_f64, f64::max)
            .next_up()
            + screw.head_top_z().abs())
        .next_up();
        let world_extent = (frame
            .origin
            .into_iter()
            .map(f64::abs)
            .fold(0_f64, f64::max)
            .next_up()
            + (scale * local_extent).next_up())
        .next_up();
        let rounding = (64. * f64::EPSILON * world_extent).next_up();
        if rounding > deflection / 4. {
            return Err(input(
                "Curved drive display placement precision is unresolved",
            ));
        }
        let points = local
            .into_iter()
            .map(|mut p| {
                p[2] += screw.head_top_z();
                position(frame, p)
            })
            .collect::<Vec<_>>();
        *remaining -= points.len();
        edges.push(DisplayEdge {
            index: ordinal,
            points_mm: points,
        });
    }
    Ok(edges)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn distance(p: Point, a: Point, b: Point) -> f64 {
        let d: Point = std::array::from_fn(|i| b[i] - a[i]);
        let v: Point = std::array::from_fn(|i| p[i] - a[i]);
        let t = (dot(d, v) / dot(d, d)).clamp(0., 1.);
        (0..3)
            .map(|i| (v[i] - t * d[i]).powi(2))
            .sum::<f64>()
            .sqrt()
    }
    #[test]
    fn actual_cone_intersection_display_stays_within_deflection_in_proper_poses() {
        for left in [false, true] {
            let screw = super::super::tests::source(left, 0.);
            let body = HelicalBody::new(
                *screw.thread(),
                EntityId::new(),
                HelicalForm::ProfiledScrew(Box::new(screw.clone())),
                None,
                GeometryPolicy::default(),
            )
            .unwrap();
            let source = screw.recess().hexalobular_geometry().unwrap();
            for frame in [
                PlaneFrame::default(),
                PlaneFrame {
                    origin: [100., -30., 50.],
                    u: [0.8, 0.6, 0.],
                    v: [-0.36, 0.48, 0.8],
                },
            ] {
                let mut remaining = MAX_DISPLAY_EDGE_POINTS;
                let deflection = 0.01;
                let actual = edges(&body, &screw, frame, deflection, &mut remaining).unwrap();
                assert_eq!(actual.len(), 12);
                assert!(remaining < MAX_DISPLAY_EDGE_POINTS);
                let all = body.display_edges(frame, deflection).unwrap();
                for (i, edge) in actual.iter().enumerate() {
                    assert!(
                        all.iter()
                            .any(|e| e.index == edge.index && e.points_mm == edge.points_mm)
                    );
                    assert!(edge.points_mm.len() >= 2);
                    let arc = source.arcs()[i];
                    for sample in 0..=200 {
                        let theta = arc.start_angle_rad + arc.sweep_rad * sample as f64 / 200.;
                        let x = arc.center_mm[0] + arc.radius_mm * theta.portable_cos();
                        let y = arc.center_mm[1] + arc.radius_mm * theta.portable_sin();
                        let z = screw.head_top_z() + source.cone().tip_depth_mm
                            - x.portable_hypot(y) * source.cone().depth_per_radius;
                        let p = position(frame, [x, y, z]);
                        let error = edge
                            .points_mm
                            .windows(2)
                            .map(|pair| distance(p, pair[0], pair[1]))
                            .fold(f64::INFINITY, f64::min);
                        assert!(error <= deflection, "arc{i} actual error{error}");
                    }
                }
            }
            let before = screw.clone();
            let mut precision_budget = MAX_DISPLAY_EDGE_POINTS;
            assert!(edges(&body, &screw, PlaneFrame::default(), 0.01, &mut 1).is_err());
            assert!(
                edges(
                    &body,
                    &screw,
                    PlaneFrame::default(),
                    1e-18,
                    &mut precision_budget
                )
                .is_err()
            );
            assert!(
                edges(
                    &body,
                    &screw,
                    PlaneFrame::default(),
                    f64::NAN,
                    &mut precision_budget
                )
                .is_err()
            );
            assert_eq!(screw, before);
        }
    }
}
