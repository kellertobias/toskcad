use super::*;
#[test]
fn actual_phillips_writer_retains_all_true_support_families_and_shared_shell() {
    let spec = crate::phillips_drive::Spec {
        mouth_diameter_mm: 1.9,
        mouth_radius_mm: 0.3,
        b_width_mm: 0.61,
        g_diameter_mm: 0.81,
        tilt_angle_deg: 7.,
        outer_flank_angle_deg: 26.5,
        floor_angle_deg: 28.,
        outer_corner_angle_deg: 92.,
        inner: crate::phillips_drive::Inner::Rounded { radius_mm: 0.25 },
    };
    let cutter = Cutter::new(
        spec,
        crate::phillips_drive::Domain {
            corner_start_depth_mm: -4.,
            corner_end_depth_mm: 4.,
            corner_reach_mm: 32.,
        },
    )
    .unwrap();
    let output = write_phillips_cutter(&cutter, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7).unwrap();
    assert_eq!(output.matches("ADVANCED_FACE(").count(), 29);
    assert_eq!(output.matches("EDGE_CURVE(").count(), 76);
    assert_eq!(output.matches("VERTEX_POINT(").count(), 49);
    assert_eq!(output.matches("INTERSECTION_CURVE(").count(), 32);
    assert!(output.contains("ELLIPSE("));
    assert!(output.contains("TOROIDAL_SURFACE("));
    assert!(output.contains("CONICAL_SURFACE("));
    assert!(output.contains("CYLINDRICAL_SURFACE("));
}
#[test]
fn all_metric_phillips_tools_export_in_proper_and_far_frames() {
    let u = [0.36, 0.48, 0.8];
    let v = [-0.8, 0.6, 0.];
    for metric in 0..=5 {
        let (spec, domain) = crate::phillips_drive::test_fixture(metric);
        let source = Cutter::new(spec, domain).unwrap();
        for (origin, u, v) in [
            ([0.; 3], [1., 0., 0.], [0., 1., 0.]),
            ([13., -7., 4.], u, v),
            ([1e6, -2e6, 3e6], u, v),
        ] {
            let result = write_phillips_cutter(&source, origin, u, v, 1e-7);
            assert!(result.is_ok(), "metric{metric} at{origin:?}: {result:?}");
        }
    }
}
#[test]
fn cancelled_and_failed_cutter_appends_preserve_prior_members_and_retry_bytes() {
    use std::cell::Cell;
    let (s, d) = crate::phillips_drive::test_fixture(1);
    let source = Cutter::new(s, d).unwrap();
    let u = [1., 0., 0.];
    let v = [0., 1., 0.];
    let origin = [7., -2., 3.];
    let calls = Cell::new(0);
    write_phillips_cutter_with_interrupt(
        &source,
        origin,
        u,
        v,
        1e-7,
        Some(&|| {
            calls.set(calls.get() + 1);
            false
        }),
    )
    .unwrap();
    let total = calls.get();
    assert!(total > 100);
    let mut expected = Writer::new();
    expected
        .append_phillips_cutter(&source, [0.; 3], u, v, 1e-7)
        .unwrap();
    expected
        .append_phillips_cutter(&source, origin, u, v, 1e-7)
        .unwrap();
    let expected = expected.finish().unwrap();
    for at in [0, 1, total / 4, total / 2, total * 3 / 4, total - 1] {
        let mut writer = Writer::new();
        writer
            .append_phillips_cutter(&source, [0.; 3], u, v, 1e-7)
            .unwrap();
        let count = Cell::new(0);
        assert_eq!(
            writer.append_phillips_cutter_with_interrupt(
                &source,
                origin,
                u,
                v,
                1e-7,
                Some(&|| {
                    let n = count.get();
                    count.set(n + 1);
                    n >= at
                })
            ),
            Err(Error::Cancelled)
        );
        writer
            .append_phillips_cutter(&source, origin, u, v, 1e-7)
            .unwrap();
        assert_eq!(writer.finish().unwrap(), expected);
    }
    for (tolerance, origin, u, v, error) in [
        (0., origin, u, v, Error::InvalidGeometry),
        (f64::NAN, origin, u, v, Error::InvalidGeometry),
        (1e-30, origin, u, v, Error::Limit),
        (1e-7, [1e13; 3], u, v, Error::Limit),
        (1e-7, origin, [0.; 3], v, Error::InvalidPlacement),
    ] {
        let mut writer = Writer::new();
        writer
            .append_phillips_cutter(&source, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7)
            .unwrap();
        assert_eq!(
            writer.append_phillips_cutter(&source, origin, u, v, tolerance),
            Err(error)
        );
        writer
            .append_phillips_cutter(&source, [7., -2., 3.], [1., 0., 0.], [0., 1., 0.], 1e-7)
            .unwrap();
        assert_eq!(writer.finish().unwrap(), expected);
    }
}
#[test]
fn actual_changed_corner_event_roles_and_unqualified_finite_caps_stay_explicit() {
    let (mut spec, domain) = crate::phillips_drive::test_fixture(1);
    spec.inner = crate::phillips_drive::Inner::ShallowV {
        angle_deg: 138.,
        e_width_mm: 0.55,
    };
    let source = Cutter::new(spec, domain).unwrap();
    let output = write_phillips_cutter(&source, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7).unwrap();
    assert!(
        source
            .boundary_events()
            .unwrap()
            .trims()
            .iter()
            .any(
                |t| t.envelope_role() == crate::phillips_drive::EnvelopeRole::OuterFlank
                    && t.wall_kind() == crate::phillips_drive::WallKind::Inner
            )
    );
    assert_eq!(output.matches("INTERSECTION_CURVE(").count(), 32);
    let (s, mut d) = crate::phillips_drive::test_fixture(1);
    d.corner_reach_mm = 0.001;
    let source = Cutter::new(s, d).unwrap();
    assert_eq!(
        write_phillips_cutter(&source, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7),
        Err(Error::Unsupported)
    );
}
#[test]
fn actual_literal_curve_buffers_refuse_gaps_wrong_supports_and_exhausted_precision() {
    use curves::{Curve, Literal};
    use support::{Frame, Kind, Pose, Support};
    let (s, d) = crate::phillips_drive::test_fixture(0);
    let source = Cutter::new(s, d).unwrap();
    let topology = source.boundary_topology().unwrap();
    let pose = Pose::new([0.; 3], [1., 0., 0.], [0., 1., 0.]).unwrap();
    let points = topology
        .vertices()
        .iter()
        .map(|v| pose.point(v.point_mm_bounds().map(support::b).map(I::midpoint)))
        .collect::<Vec<_>>();
    let [a, b, c, d] = [0, 1, 2, 3].map(|i| topology.events().trim_carrier(i, 1e-9, 4096).unwrap());
    let spans = [a, b, c, d];
    let circle = topology
        .edges()
        .iter()
        .enumerate()
        .find(|(i, _)| {
            matches!(
                topology.curve_analytic_support(*i).unwrap(),
                Some(AnalyticCurve::CircleArc { .. })
            )
        })
        .unwrap()
        .0;
    assert!(matches!(
        Curve::new(&topology, circle, &pose, &points, &spans, 1e-30),
        Err(Error::Limit)
    ));
    let far = Pose::new([1e308; 3], [1., 0., 0.], [0., 1., 0.]).unwrap();
    assert!(matches!(
        Curve::new(&topology, circle, &far, &points, &spans, 1e-7),
        Err(Error::Unsupported)
    ));
    let line = topology
        .edges()
        .iter()
        .enumerate()
        .find(|(i, _)| {
            matches!(
                topology.curve_analytic_support(*i).unwrap(),
                Some(AnalyticCurve::Line(_))
            )
        })
        .unwrap()
        .0;
    let mut collapsed = points.clone();
    let [a, b] = topology.edges()[line].vertices();
    collapsed[b] = collapsed[a];
    assert!(matches!(
        Curve::new(&topology, line, &pose, &collapsed, &spans, 1e-7),
        Err(Error::Limit)
    ));
    let trim = topology
        .edges()
        .iter()
        .enumerate()
        .find(|(_, e)| {
            matches!(
                e.source_identity(),
                CurveKey::Trim {
                    envelope: crate::phillips_drive::EnvelopeRole::EntryRoundover,
                    ..
                }
            )
        })
        .unwrap()
        .0;
    let mut missing = spans.clone();
    missing[0].clear();
    assert!(matches!(
        Curve::new(&topology, trim, &pose, &points, &missing, 1e-7),
        Err(Error::Unsupported)
    ));
    let mut gap = spans.clone();
    assert!(gap[0].len() > 1);
    gap[0][1].poles_mm[0][0] += 1e-12;
    assert!(matches!(
        Curve::new(&topology, trim, &pose, &points, &gap, 1e-7),
        Err(Error::Unsupported)
    ));
    let mut changed = spans.clone();
    changed[0][0].poles_mm[1][0] += 0.01;
    assert!(matches!(
        Curve::new(&topology, trim, &pose, &points, &changed, 1e-7),
        Err(Error::Limit)
    ));
    assert_eq!(
        curves::trim_index(&topology, CurveKey::Mouth { edge: 0 }),
        Err(Error::Unsupported)
    );
    let analytic = topology.curve_analytic_support(circle).unwrap().unwrap();
    let ends = topology.edges()[circle].vertices().map(|i| points[i]);
    let raw = curves::analytic_literal(analytic, &pose, ends).unwrap();
    let bad = Literal::Line {
        origin: [0.; 3],
        direction: [1., 0., 0.],
        forward: true,
    };
    assert_eq!(
        curves::analytic_distance(analytic, &pose, &bad, source.bounding_box_mm()),
        Err(Error::Unsupported)
    );
    assert!(matches!(raw, Literal::Circle { .. }));
    if let Literal::Circle {
        axis,
        radius,
        forward,
    } = raw.clone()
    {
        for bad in [
            Literal::Circle {
                axis,
                radius,
                forward: !forward,
            },
            Literal::Circle {
                axis,
                radius: 0.,
                forward,
            },
            Literal::Ellipse {
                axis,
                semiaxes: [radius; 2],
                forward,
            },
            Literal::Circle {
                axis,
                radius: 1e300,
                forward,
            },
        ] {
            assert_eq!(
                curves::analytic_distance(analytic, &pose, &bad, source.bounding_box_mm()),
                Err(Error::Unsupported)
            );
        }
    }
    let mut backwards = spans[0].clone();
    backwards[0].poles_mm[1][1] = backwards[0].poles_mm[0][1] + 0.01;
    assert_eq!(
        curves::monotone(&backwards, pose.axes),
        Err(Error::Unsupported)
    );
    let cone = topology
        .faces()
        .iter()
        .enumerate()
        .map(|(i, _)| topology.face_surface(i).unwrap())
        .find(|s| matches!(s, crate::phillips_drive::Surface::Cone { .. }))
        .unwrap();
    let cone = Support::new(cone).unwrap();
    let mut bad = cone.literal(&pose);
    if let Kind::Cone(r, _) = bad.kind {
        bad.kind = Kind::Cone(r, 1e-300);
    }
    assert_eq!(
        cone.distance_bound(&pose, &bad, source.bounding_box_mm()),
        Err(Error::Unsupported)
    );
    let line = topology.curve_analytic_support(line).unwrap().unwrap();
    assert_eq!(
        curves::analytic_distance(line, &pose, &raw, source.bounding_box_mm()),
        Err(Error::Unsupported)
    );
    let plane = Support::new(crate::phillips_drive::Surface::Plane {
        normal_bounds: [[0., 0.], [0., 0.], [1., 1.]],
        offset_mm_bounds: [0., 0.],
    })
    .unwrap();
    let mut literal = plane.literal(&pose);
    literal.axis.normal = [0., 0., -1.];
    assert_eq!(
        plane.distance_bound(&pose, &literal, source.bounding_box_mm()),
        Err(Error::Unsupported)
    );
    literal = plane.literal(&pose);
    literal.kind = Kind::Cylinder(1.);
    assert_eq!(
        plane.distance_bound(&pose, &literal, source.bounding_box_mm()),
        Err(Error::Unsupported)
    );
    assert!(matches!(
        Frame::new([I::point(0.); 3], [I::point(0.); 3], [I::point(1.); 3]),
        Err(Error::InvalidGeometry)
    ));
}
#[test]
fn tiny_positive_roundovers_and_frame_strain_cannot_disappear_in_step() {
    let (mut s, d) = crate::phillips_drive::test_fixture(1);
    s.mouth_radius_mm = 1e-12;
    let source = Cutter::new(s, d).unwrap();
    assert_eq!(
        write_phillips_cutter(&source, [1e6; 3], [1., 0., 0.], [0., 1., 0.], 1e-7),
        Err(Error::Limit)
    );
    s.mouth_radius_mm = 1e-8;
    let source = Cutter::new(s, d).unwrap();
    assert_eq!(
        write_phillips_cutter(&source, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7),
        Err(Error::Limit)
    );
    let (s, d) = crate::phillips_drive::test_fixture(1);
    let source = Cutter::new(s, d).unwrap();
    assert_eq!(
        write_phillips_cutter(&source, [0.; 3], [1., 0., 0.], [5e-11, 1., 0.], 1e-9),
        Err(Error::Limit)
    );
    for precision in [1e-11, 5e-12, 3e-12, 1e-12] {
        assert!(
            write_phillips_cutter(&source, [0.; 3], [1., 0., 0.], [0., 1., 0.], precision).is_err()
        );
    }
}
#[test]
fn generic_writer_and_cutter_errors_preserve_the_export_api_contract() {
    assert_eq!(
        Writer::new().finish().map_err(Error::from),
        Err(Error::InvalidGeometry)
    );
    assert_eq!(Error::from(StepError::Unsupported), Error::Unsupported);
    assert_eq!(Error::from(StepError::Limit), Error::Limit);
    let (mut s, d) = crate::phillips_drive::test_fixture(1);
    s.mouth_diameter_mm = f64::NAN;
    assert_eq!(
        Cutter::new(s, d).map_err(Error::from).unwrap_err(),
        Error::InvalidGeometry
    );
    assert_eq!(
        Error::from(crate::phillips_drive::Error::Limit),
        Error::Limit
    );
}
