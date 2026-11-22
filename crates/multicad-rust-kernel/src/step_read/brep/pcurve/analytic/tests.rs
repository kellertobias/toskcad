use super::*;
fn frame() -> Frame {
    Frame {
        origin: [1., -2., 3.],
        x: [0.6, 0.8, 0.],
        y: [-0.8, 0.6, 0.],
        z: [0., 0., 1.],
    }
}
fn axis() -> Frame {
    Frame {
        origin: [0.; 3],
        x: [1., 0., 0.],
        y: [0., 1., 0.],
        z: [0., 0., 1.],
    }
}
fn line(origin: Point, direction: Point) -> Curve {
    Curve::Line { origin, direction }
}
fn edge(curve: Curve) -> Edge {
    Edge {
        curve,
        range: [-0.2, 1.4],
        vertices: [0, 1],
    }
}
fn periodic(frame: Frame, radius: f64) -> Curve {
    Curve::Circle { frame, radius }
}
fn check(surface: &Surface, uv: &Curve, curve: Curve) -> f64 {
    let mut budget = Budget::new(16_000_000, None);
    match certify(surface, uv, &edge(curve), 1e-7, &mut budget).unwrap() {
        Proof::Certified(error) => {
            assert!((0. ..=1e-7).contains(&error));
            error
        }
        Proof::Disagreement(error) => {
            panic!("valid analytic support refused by lower witness {error}")
        }
        Proof::Unsupported => panic!("valid analytic support unproved"),
    }
}
fn add_p(a: Point, b: Point) -> Point {
    std::array::from_fn(|i| a[i] + b[i])
}
fn scale_p(a: Point, b: f64) -> Point {
    a.map(|x| x * b)
}
fn radial_p(f: Frame, u: f64) -> Point {
    add_p(scale_p(f.x, u.cos()), scale_p(f.y, u.sin()))
}
fn rephase(mut f: Frame, angle: f64, reverse: bool) -> Frame {
    let (s, c) = angle.sin_cos();
    let (x, y) = (f.x, f.y);
    f.x = add_p(scale_p(x, c), scale_p(y, s));
    f.y = add_p(scale_p(x, -s), scale_p(y, c));
    if reverse {
        f.y = scale_p(f.y, -1.);
        f.z = scale_p(f.z, -1.);
    }
    f
}
#[test]
fn planar_conics_cover_positive_phase_and_oriented_swapped_ellipse_axes() {
    let f = frame();
    let uv = periodic(axis(), 2.5);
    check(
        &Surface::Plane(f),
        &uv,
        periodic(rephase(f, 0.73, false), 2.5),
    );
    let uv = Curve::Ellipse {
        frame: axis(),
        semiaxes: [2., 3.],
    };
    let swapped = Frame {
        x: f.y,
        y: scale_p(f.x, -1.),
        z: f.z,
        ..f
    };
    check(
        &Surface::Plane(f),
        &uv,
        Curve::Ellipse {
            frame: swapped,
            semiaxes: [3., 2.],
        },
    );
    // A source LINE may have any nonzero scale, origin, and orientation.
    let uv = line([10., -7., 0.], [-3., 2., 0.]);
    let p = add_p(f.origin, add_p(scale_p(f.x, 10.), scale_p(f.y, -7.)));
    let d = add_p(scale_p(f.x, -3.), scale_p(f.y, 2.));
    check(
        &Surface::Plane(f),
        &uv,
        line(add_p(p, scale_p(d, -31.)), scale_p(d, 7.)),
    );
}
#[test]
fn cylinder_and_cone_isocurves_prove_complete_natural_images() {
    let f = frame();
    let radius = 2.5;
    for cone in [false, true] {
        let angle = 0.4f64;
        let surface = if cone {
            Surface::Cone {
                frame: f,
                radius,
                half_angle: angle,
            }
        } else {
            Surface::Cylinder { frame: f, radius }
        };
        let v = 1.25;
        let rho = radius + if cone { v * angle.tan() } else { 0. };
        let circle_frame = Frame {
            origin: add_p(f.origin, scale_p(f.z, v)),
            ..f
        };
        check(
            &surface,
            &line([0.37, v, 0.], [-7., 0., 0.]),
            periodic(rephase(circle_frame, -0.71, true), rho),
        );
        let u = 0.71;
        let r = radial_p(f, u);
        let generator = add_p(f.z, scale_p(r, if cone { angle.tan() } else { 0. }));
        let p = add_p(add_p(f.origin, scale_p(r, radius)), scale_p(generator, 3.));
        check(
            &surface,
            &line([u, -5., 0.], [0., -4., 0.]),
            line(p, scale_p(generator, -3.)),
        );
    }
}
#[test]
fn sphere_latitudes_and_both_torus_isocurves_cover_full_circles() {
    let f = frame();
    let radius = 2.5;
    let v = 0.47f64;
    let sphere = Surface::Sphere { frame: f, radius };
    let center = add_p(f.origin, scale_p(f.z, radius * v.sin()));
    check(
        &sphere,
        &line([-4., v, 0.], [3., 0., 0.]),
        periodic(
            rephase(
                Frame {
                    origin: center,
                    ..f
                },
                0.23,
                false,
            ),
            radius * v.cos(),
        ),
    );
    let (major, minor) = (4., 0.75);
    let torus = Surface::Torus {
        frame: f,
        major,
        minor,
    };
    let center = add_p(f.origin, scale_p(f.z, minor * v.sin()));
    check(
        &torus,
        &line([20., v, 0.], [-2., 0., 0.]),
        periodic(
            rephase(
                Frame {
                    origin: center,
                    ..f
                },
                -0.3,
                true,
            ),
            major + minor * v.cos(),
        ),
    );
    let u = 0.71;
    let r = radial_p(f, u);
    let meridian = Frame {
        origin: add_p(f.origin, scale_p(r, major)),
        x: r,
        y: f.z,
        z: [r[1], -r[0], 0.],
    };
    check(
        &torus,
        &line([u, -30., 0.], [0., 0.5, 0.]),
        periodic(rephase(meridian, 0.8, false), minor),
    );
}
#[test]
fn real_pointset_lower_witnesses_do_not_assume_equal_parameters() {
    let surface = Surface::Plane(axis());
    let shifted = Frame {
        origin: [1., 2., 0.],
        ..axis()
    };
    let mut budget = Budget::new(16_000_000, None);
    let proof = certify(
        &surface,
        &periodic(shifted, 3.),
        &edge(periodic(axis(), 3.)),
        1e-7,
        &mut budget,
    )
    .unwrap();
    assert!(matches!(proof,Proof::Disagreement(lower) if lower>1e-7));
    let surface = Surface::Cylinder {
        frame: axis(),
        radius: 2.5,
    };
    let mut budget = Budget::new(16_000_000, None);
    assert!(
        matches!(certify(&surface,&line([7.,0.,0.],[2.,0.,0.]),&edge(periodic(axis(),2.4)),1e-7,&mut budget).unwrap(),Proof::Disagreement(lower) if lower>0.09)
    );
    let surface = Surface::Plane(Frame {
        origin: [0., 0., 1.],
        ..axis()
    });
    let mut budget = Budget::new(16_000_000, None);
    assert!(
        matches!(certify(&surface,&line([0.;3],[1.,0.,0.]),&edge(line([0.;3],[1.,0.,0.])),1e-7,&mut budget).unwrap(),Proof::Disagreement(lower) if lower>0.99)
    );
}
#[test]
fn diagonal_and_spherical_half_meridian_domains_are_not_promoted_to_full_supports() {
    let mut budget = Budget::new(16_000_000, None);
    let cylinder = Surface::Cylinder {
        frame: axis(),
        radius: 2.,
    };
    assert!(matches!(
        certify(
            &cylinder,
            &line([0.; 3], [1., 1e-300, 0.]),
            &edge(periodic(axis(), 2.)),
            1e-7,
            &mut budget
        )
        .unwrap(),
        Proof::Unsupported
    ));
    let sphere = Surface::Sphere {
        frame: axis(),
        radius: 2.,
    };
    for uv in [
        line([0.; 3], [0., 1., 0.]),
        line([0., 2., 0.], [1., 0., 0.]),
    ] {
        let mut budget = Budget::new(16_000_000, None);
        assert!(matches!(
            certify(&sphere, &uv, &edge(periodic(axis(), 2.)), 1e-7, &mut budget).unwrap(),
            Proof::Unsupported
        ));
    }
}
#[test]
fn analytic_proof_work_frontier_and_cancel_preserve_retry() {
    let surface = Surface::Cylinder {
        frame: axis(),
        radius: 2.,
    };
    let uv = line([5., 0., 0.], [3., 0., 0.]);
    let edge = edge(periodic(axis(), 2.));
    let initial = 16_000_000;
    let mut budget = Budget::new(initial, None);
    assert!(matches!(
        certify(&surface, &uv, &edge, 1e-7, &mut budget).unwrap(),
        Proof::Certified(_)
    ));
    let used = initial - budget.remaining;
    let mut exact = Budget::new(used, None);
    assert!(matches!(
        certify(&surface, &uv, &edge, 1e-7, &mut exact).unwrap(),
        Proof::Certified(_)
    ));
    assert_eq!(exact.remaining, 0);
    let mut limited = Budget::new(used - 1, None);
    assert!(certify(&surface, &uv, &edge, 1e-7, &mut limited).is_err());
    let polls = std::cell::Cell::new(0);
    let stop = || {
        polls.set(polls.get() + 1);
        polls.get() == 3
    };
    let mut interrupted = Budget::new(initial, Some(&stop));
    assert!(certify(&surface, &uv, &edge, 1e-7, &mut interrupted).is_err());
    let mut retry = Budget::new(initial, None);
    assert!(matches!(
        certify(&surface, &uv, &edge, 1e-7, &mut retry).unwrap(),
        Proof::Certified(_)
    ));
    assert_eq!(retry.remaining, budget.remaining);
}

#[test]
fn same_point_set_with_opposite_sense_never_certifies_or_fakes_distance() {
    let plane = Surface::Plane(axis());
    for (uv, curve) in [
        (
            periodic(axis(), 2.),
            periodic(rephase(axis(), 0.7, true), 2.),
        ),
        (
            line([0.; 3], [1., 0., 0.]),
            line([10., 0., 0.], [-2., 0., 0.]),
        ),
    ] {
        let mut budget = Budget::new(16_000_000, None);
        assert!(matches!(
            certify(&plane, &uv, &edge(curve), 1e-7, &mut budget).unwrap(),
            Proof::Unsupported
        ));
    }
    // Explicit orientation is required even when a loop is smaller than tol.
    let mut budget = Budget::new(16_000_000, None);
    assert!(matches!(
        certify(
            &plane,
            &periodic(axis(), 1e-12),
            &edge(periodic(rephase(axis(), 0., true), 1e-12)),
            1e-7,
            &mut budget
        )
        .unwrap(),
        Proof::Unsupported
    ));
    let cylinder = Surface::Cylinder {
        frame: axis(),
        radius: 2.,
    };
    let mut budget = Budget::new(16_000_000, None);
    assert!(matches!(
        certify(
            &cylinder,
            &line([0.; 3], [-3., 0., 0.]),
            &edge(periodic(axis(), 2.)),
            1e-7,
            &mut budget
        )
        .unwrap(),
        Proof::Unsupported
    ));
}
#[test]
fn finite_input_overflow_cannot_be_hidden_by_maximum_reduction() {
    let surface = Surface::Plane(axis());
    let uv = line([0.; 3], [1., 0., 0.]);
    let huge = Edge {
        curve: line([1e308, 0., 0.], [1e308, 0., 0.]),
        range: [1., 2.],
        vertices: [0, 1],
    };
    let mut budget = Budget::new(16_000_000, None);
    assert!(matches!(
        certify(&surface, &uv, &huge, 1e-7, &mut budget).unwrap(),
        Proof::Unsupported
    ));
}
