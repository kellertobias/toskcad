use super::*;
fn i([a, b]: Bounds) -> I {
    I::new(a, b)
}
fn zero(x: I) {
    assert!(x.valid() && x.lo <= 0. && x.hi >= 0., "{x:?}");
}
fn rotate([x, y]: [I; 2], q: u8) -> [I; 2] {
    match q {
        0 => [x, y],
        1 => [-y, x],
        2 => [-x, -y],
        _ => [y, -x],
    }
}
#[test]
fn every_metric_analytic_edge_retains_its_actual_support_and_finite_source_domain() {
    for metric in 0..=5 {
        let (s, d) = super::super::super::tests::fixture(metric);
        let source = Cutter::new(s, d).unwrap();
        let topology = source.boundary_topology().unwrap();
        let mut kinds = [0usize; 4];
        for (index, edge) in topology.edges().iter().enumerate() {
            let support = topology.curve_analytic_support(index).unwrap();
            if support.is_none() {
                assert!(matches!(edge.source_identity(), CurveKey::Trim { .. }));
                kinds[3] += 1;
                continue;
            }
            for t in [0., 0.25, 0.5, 0.75, 1.] {
                let p = topology.curve_point_mm_bounds(index, t).unwrap().map(i);
                match support.unwrap() {
                    AnalyticCurve::Line(line) => {
                        kinds[0] += 1;
                        let d = line.direction_bounds.map(i);
                        let o = line.origin_mm_bounds.map(i);
                        let v = std::array::from_fn::<_, 3, _>(|j| p[j] - o[j]);
                        for j in 0..3 {
                            zero(v[(j + 1) % 3] * d[(j + 2) % 3] - v[(j + 2) % 3] * d[(j + 1) % 3]);
                        }
                    }
                    AnalyticCurve::CircleArc {
                        radius_mm_bounds,
                        depth_mm_bounds,
                        start_angle_rad_bounds,
                        sweep_rad_bounds,
                    } => {
                        kinds[1] += 1;
                        let r = i(radius_mm_bounds);
                        zero(p[0].square() + p[1].square() - r.square());
                        zero(p[2] - i(depth_mm_bounds));
                        let theta = i(start_angle_rad_bounds) + I::point(t) * i(sweep_rad_bounds);
                        zero(p[0] - r * theta.cos());
                        zero(p[1] - r * theta.sin());
                    }
                    AnalyticCurve::EllipseArc {
                        center_mm_bounds,
                        semiaxes_mm_bounds,
                        quarter_turn,
                        start_angle_rad_bounds,
                        sweep_rad_bounds,
                    } => {
                        kinds[2] += 1;
                        let c = center_mm_bounds.map(i);
                        let [a, b] = semiaxes_mm_bounds.map(i);
                        let [x, y] = rotate([p[0] - c[0], p[1] - c[1]], (4 - quarter_turn) % 4);
                        zero((x / a).square() + (y / b).square() - I::point(1.));
                        zero(p[2] - c[2]);
                        let theta = i(start_angle_rad_bounds) + I::point(t) * i(sweep_rad_bounds);
                        zero(x - a * theta.cos());
                        zero(y - b * theta.sin());
                    }
                }
            }
        }
        assert!(kinds[0] > 0 && kinds[1] > 0);
        assert_eq!(kinds[3], 32);
        assert_eq!(kinds[2] > 0, metric == 0);
        assert_eq!(
            topology.curve_analytic_support(usize::MAX),
            Err(Error::Invalid)
        );
    }
}
