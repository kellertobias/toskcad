use super::*;

fn line(a: Point, b: Point) -> Edge {
    Edge::Line { a, b }
}

fn arc(center: Point, a: Point, b: Point) -> Edge {
    Edge::Arc {
        center,
        a,
        b,
        clockwise: false,
    }
}

/// A 160 x 100 plate with R3 corners. `skew` moves every arc centre away
/// from its corner along the diagonal, as a rounded three-point arc does.
/// Each centre is written with the diagonal toward its corner.
fn plate(skew: f64) -> Vec<Edge> {
    let c = |x: f64, y: f64, sx: f64, sy: f64| [x - sx * skew, y - sy * skew];
    vec![
        line([3., 0.], [157., 0.]),
        arc(c(157., 3., 1., -1.), [157., 0.], [160., 3.]),
        line([160., 3.], [160., 97.]),
        arc(c(157., 97., 1., 1.), [160., 97.], [157., 100.]),
        line([157., 100.], [3., 100.]),
        arc(c(3., 97., -1., 1.), [3., 100.], [0., 97.]),
        line([0., 97.], [0., 3.]),
        arc(c(3., 3., -1., -1.), [0., 3.], [3., 0.]),
    ]
}

fn hole(center: Point, radius_mm: f64) -> Vec<Edge> {
    vec![Edge::Circle {
        center,
        radius_mm,
        clockwise: true,
    }]
}

fn area(polygon: &[Point]) -> f64 {
    (0..polygon.len())
        .map(|i| {
            let [p, q] = [polygon[i], polygon[(i + 1) % polygon.len()]];
            p[0] * q[1] - q[0] * p[1]
        })
        .sum::<f64>()
        / 2.
}

#[test]
fn rounded_plate_splits_into_corner_quarters_and_one_stepped_polygon() {
    let profile = Profile::new(&[plate(0.)], 1e-7).unwrap();
    assert!(profile.rectangles.is_empty());
    let [body] = &profile.polygons[..] else {
        panic!("one polygon")
    };
    // Each corner is replaced by the chords to its centre.
    assert_eq!(body.len(), 12);
    assert!((area(body) - (160. * 100. - 4. * 9.)).abs() < 1e-9);
    assert_eq!(profile.bends.len(), 4);
    let mut centers: Vec<_> = profile.bends.iter().map(|bend| bend.center).collect();
    centers.sort_by(|a, b| a.partial_cmp(b).unwrap());
    assert_eq!(
        centers,
        [[3., 3.], [3., 97.], [157., 3.], [157., 97.]].to_vec()
    );
    for bend in &profile.bends {
        assert_eq!(bend.inner, None);
        assert_eq!(bend.outer, 3.);
        assert_eq!(bend.authored, bend.center);
        let [s, t] = bend.spokes;
        assert_eq!(t, [-s[1], s[0]]);
    }
}

#[test]
fn rounded_midpoint_arcs_snap_only_within_the_linear_tolerance() {
    // The control-desk board: mid (0.8786797, 0.8786797) puts each centre
    // 1.487e-7 mm off its corner on both axes. The authored arc then strays
    // at most (sqrt 2 - 1) * 1.487e-7 = 6.16e-8 mm from the exact quarter.
    let skew = 1.487e-7;
    let snapped = Profile::new(&[plate(skew)], 1e-7).unwrap();
    for bend in &snapped.bends {
        assert!(bend.center.iter().all(|v| [3., 97., 157.].contains(v)));
        assert_eq!(bend.outer, 3.);
    }
    assert!(Profile::new(&[plate(skew)], 6e-8).is_none());
    assert!(Profile::new(&[plate(skew)], 6.2e-8).is_some());
    // A larger construction error is refused at the default tolerance.
    assert!(Profile::new(&[plate(1e-6)], 1e-7).is_none());
}

/// A clockwise rounded-slot hole ring along +x, ordered like the
/// shared-workflow bracket's slot.
fn slot(centers: [Point; 2], radius: f64) -> Vec<Edge> {
    let [c0, c1] = centers;
    let [top0, bottom0] = [[c0[0], c0[1] + radius], [c0[0], c0[1] - radius]];
    let [top1, bottom1] = [[c1[0], c1[1] + radius], [c1[0], c1[1] - radius]];
    vec![
        line(top0, top1),
        Edge::Arc {
            center: c1,
            a: top1,
            b: bottom1,
            clockwise: true,
        },
        line(bottom1, bottom0),
        Edge::Arc {
            center: c0,
            a: bottom0,
            b: top0,
            clockwise: true,
        },
    ]
}

#[test]
fn holes_must_clear_every_edge_of_one_piece() {
    let profile = Profile::new(
        &[
            plate(0.),
            hole([20., 20.], 1.6),
            slot([[60., 30.], [75., 30.]], 5.),
        ],
        1e-7,
    )
    .unwrap();
    assert_eq!(profile.holes.len(), 2);
    assert_eq!(profile.holes[0].centers, [[20., 20.]; 2]);
    assert_eq!(profile.holes[1].centers, [[75., 30.], [60., 30.]]);
    assert_eq!(profile.holes[1].radius, 5.);
    assert!(profile.holes.iter().all(|hole| hole.host == 0));

    // Near a corner the hole must clear the chords to the corner centre.
    assert!(Profile::new(&[plate(0.), hole([5., 5.], 1.6)], 1e-7).is_some());
    assert!(Profile::new(&[plate(0.), hole([4., 4.], 1.6)], 1e-7).is_none());
    // Touching the outline is refused.
    assert!(Profile::new(&[plate(0.), hole([50., 1.6], 1.6)], 1e-7).is_none());
    // A hole ring must be one circle or a rounded slot.
    assert!(Profile::new(&[plate(0.), plate(0.)], 1e-7).is_none());
    // Slot ends that bulge inward are not a stadium.
    let mut pinched = slot([[60., 30.], [75., 30.]], 5.);
    if let Edge::Arc { clockwise, .. } = &mut pinched[1] {
        *clockwise = false;
    }
    assert!(Profile::new(&[plate(0.), pinched], 1e-7).is_none());
}

#[test]
fn refuses_concave_and_non_quarter_rounds() {
    // An 80 degree convex round cannot snap to a cardinal quarter.
    let (sin, cos) = 80_f64.to_radians().sin_cos();
    let end = [10. * cos, 10. * sin];
    let ring = vec![
        line([0., 0.], [10., 0.]),
        arc([0., 0.], [10., 0.], end),
        line(end, [0., 0.]),
    ];
    assert!(Profile::new(&[ring], 1e-7).is_none());
    // A lone concave (clockwise) quarter has no outer partner.
    let notch = vec![
        line([0., 0.], [10., 0.]),
        line([10., 0.], [10., 7.]),
        Edge::Arc {
            center: [10., 10.],
            a: [10., 7.],
            b: [7., 10.],
            clockwise: true,
        },
        line([7., 10.], [0., 10.]),
        line([0., 10.], [0., 0.]),
    ];
    assert!(Profile::new(&[notch], 1e-7).is_none());
}
