use super::*;
use std::cell::Cell;

fn qualify(input: BoundaryInput) -> RadialBoundary {
    BoundaryOwnership::new(input, 1e-7, None)
        .unwrap()
        .bind_surfaces(1e-7, None)
        .unwrap()
        .certify_axial(1e-7, 1e-6, None)
        .unwrap()
        .certify_caps(None)
        .unwrap()
        .certify_ordering(None)
        .unwrap()
        .certify_material(1e-7, None)
        .unwrap()
}
fn point_segment(p: Point, a: Point, b: Point) -> f64 {
    let d = std::array::from_fn::<_, 3, _>(|k| b[k] - a[k]);
    let t = ((0..3).map(|k| (p[k] - a[k]) * d[k]).sum::<f64>()
        / d.into_iter().map(|v| v * v).sum::<f64>())
    .clamp(0., 1.);
    (0..3)
        .map(|k| (p[k] - a[k] - t * d[k]).powi(2))
        .sum::<f64>()
        .sqrt()
}
fn actual_point(input: &BoundaryInput, edge: &Edge, t: f64) -> Point {
    match &input.curves[edge.curve] {
        Curve::Spline { spans, poles } => {
            let t = if edge.forward { t } else { 1. - t };
            let span = ((t * *spans as f64).floor() as usize).min(spans - 1);
            let t = (t * *spans as f64 - span as f64).clamp(0., 1.);
            let mut row = poles[9 * span..=9 * span + 9].to_vec();
            while row.len() > 1 {
                row = row
                    .windows(2)
                    .map(|p| std::array::from_fn(|k| (1. - t) * p[0][k] + t * p[1][k]))
                    .collect();
            }
            row[0]
        }
        Curve::Circle { axis, radius_mm } => {
            let norm = |v: Point| v.map(|x| x / v.into_iter().map(|x| x * x).sum::<f64>().sqrt());
            let n = norm(axis.normal);
            let dot = |a: Point, b: Point| a.into_iter().zip(b).map(|(a, b)| a * b).sum::<f64>();
            let u = norm(std::array::from_fn(|k| {
                axis.reference[k] - dot(axis.reference, n) * n[k]
            }));
            let v = [
                n[1] * u[2] - n[2] * u[1],
                n[2] * u[0] - n[0] * u[2],
                n[0] * u[1] - n[1] * u[0],
            ];
            let phase = |i: usize| {
                let p = std::array::from_fn(|k| input.vertices[i][k] - axis.origin[k]);
                dot(p, v).atan2(dot(p, u))
            };
            let a = phase(edge.vertices[0]);
            let b = phase(edge.vertices[1]);
            let sweep = if edge.forward {
                (b - a).rem_euclid(std::f64::consts::TAU)
            } else {
                -(a - b).rem_euclid(std::f64::consts::TAU)
            };
            let angle = a + t * sweep;
            std::array::from_fn(|k| {
                axis.origin[k] + radius_mm * (u[k] * angle.cos() + v[k] * angle.sin())
            })
        }
    }
}

#[test]
fn actual_spline_and_circle_edges_keep_source_endpoints_and_whole_curve_chord() {
    let original = qualify(crate::radial_tensor::tests::cam::wedge(&[1, 2]));
    for reverse in [false, true] {
        let mut input = original.input().clone();
        if reverse {
            for curve in &mut input.curves {
                match curve {
                    Curve::Spline { poles, .. } => poles.reverse(),
                    Curve::Circle { axis, .. } => axis.normal = axis.normal.map(|v| -v),
                }
            }
            for edge in &mut input.edges {
                edge.forward = !edge.forward;
            }
        }
        let body = qualify(input);
        for (i, edge) in body.input().edges.iter().enumerate() {
            let line = sample(&body, i, 0.01, 10000, None).unwrap();
            assert_eq!(line[0], body.input().vertices[edge.vertices[0]]);
            assert_eq!(
                *line.last().unwrap(),
                body.input().vertices[edge.vertices[1]]
            );
            for j in 0..=100 {
                let p = actual_point(body.input(), edge, j as f64 / 100.);
                let distance = line
                    .windows(2)
                    .map(|s| point_segment(p, s[0], s[1]))
                    .fold(f64::INFINITY, f64::min);
                assert!(distance <= 0.01 + 1e-7, "{reverse}/{i}: {distance}");
            }
        }
    }
    let thread = crate::tapered_thread::Thread::new(
        crate::helical_thread::Spec {
            major_radius_mm: 4.,
            pitch_mm: 1.25,
            length_mm: 1.25,
            depth_mm: 0.5,
            crest_width_mm: 0.25,
            lower_flank_deg: 12.,
            upper_flank_deg: 18.,
            starts: 1,
            left_hand: false,
            root_radius_mm: None,
        },
        0.03125,
        1e-7,
        1e-6,
    )
    .unwrap();
    let body = thread
        .retained_boundary_with_interrupt([[0.; 3], [1., 0., 0.], [0., 1., 0.]], 1e-7, 1e-6, None)
        .unwrap()
        .boundary;
    for i in 0..body.input().edges.len() {
        assert!(sample(&body, i, 0.025, 10000, None).unwrap().len() >= 2);
    }
}

#[test]
fn source_edge_requests_cancel_atomically_and_keep_explicit_resource_limits() {
    let body = qualify(crate::radial_tensor::tests::cam::wedge(&[1, 2]));
    let original = body.input().clone();
    for i in 0..body.input().edges.len() {
        let count = Cell::new(0);
        let observe = || {
            count.set(count.get() + 1);
            false
        };
        sample(&body, i, 0.1, 10000, Some(&observe)).unwrap();
        let total = count.get();
        for at in 1..=total {
            count.set(0);
            let cancel = || {
                count.set(count.get() + 1);
                count.get() >= at
            };
            assert_eq!(
                sample(&body, i, 0.1, 10000, Some(&cancel)).err(),
                Some(Error::Cancelled)
            );
        }
    }
    for deflection in [0., -1., f64::NAN, f64::INFINITY] {
        assert_eq!(
            sample(&body, 0, deflection, 100, None).err(),
            Some(Error::Invalid)
        );
    }
    for limit in [0, 1] {
        assert_eq!(
            sample(&body, 0, 0.01, limit, None).err(),
            Some(Error::Invalid)
        );
    }
    assert_eq!(
        sample(&body, 0, 0.01, 1_000_001, None).err(),
        Some(Error::Limit)
    );
    assert_eq!(
        sample(&body, usize::MAX, 0.01, 100, None).err(),
        Some(Error::Invalid)
    );
    let circle = body
        .input()
        .edges
        .iter()
        .position(|e| matches!(body.input().curves[e.curve], Curve::Circle { .. }))
        .unwrap();
    assert_eq!(
        sample(&body, circle, 0.001, 2, None).err(),
        Some(Error::Limit)
    );
    assert_eq!(
        split_parameter([1., 1_f64.next_up()]),
        Err(Error::Unresolved)
    );
    assert_eq!(body.input(), &original);
    assert!(sample(&body, circle, 0.001, 10000, None).is_ok());
}
