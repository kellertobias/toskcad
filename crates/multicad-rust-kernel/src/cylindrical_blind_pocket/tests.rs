use super::*;
use crate::{Circle, Rational};
fn stock() -> AxialCylinder {
    AxialCylinder::new(
        17,
        Circle::new(
            [Rational::new(0, 1).unwrap(); 2],
            Rational::new(10, 1).unwrap(),
        )
        .unwrap(),
        Rational::new(0, 1).unwrap(),
        Rational::new(10, 1).unwrap(),
    )
    .unwrap()
}
fn spec() -> Spec {
    Spec {
        radius_mm: 2.,
        depth_mm: 4.,
        seat: Some(Seat {
            radius_mm: 4.,
            depth_mm: 2.,
        }),
        minimum_wall_mm: 1.,
        minimum_floor_mm: 1.,
    }
}
fn entry(end: bool) -> Frame {
    Frame {
        origin_mm: [0., 0., if end { 10. } else { 0. }],
        axes: if end {
            [[1., 0., 0.], [0., -1., 0.], [0., 0., -1.]]
        } else {
            [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]
        },
    }
}
#[test]
fn actual_both_cap_pockets_keep_positive_annular_shoulders_and_floors() {
    for end in [false, true] {
        let p = Pocket::new(
            stock(),
            Placement::Axial,
            entry(end),
            spec(),
            1e-7,
            1e-5,
            None,
        )
        .unwrap();
        assert_eq!(p.stock().feature(), 17);
        assert_eq!(p.spec(), spec());
        assert_eq!(p.end(), end);
        assert_eq!(p.roles().len(), p.geometry().rings()[0].len());
        assert!(p.roles().contains(&Role::Shoulder));
        assert!(p.source_error_mm() < 1e-10 && p.placement_error_mm() < 1e-10);
        assert_eq!(p.frame().origin_mm, [0.; 3]);
        assert!((p.geometry().volume_mm3() - (1000. - 40.) * std::f64::consts::PI).abs() < 1e-10);
        let z = |d: f64| if end { 10. - d } else { d };
        assert!(!p.geometry().contains([3., 0., z(1.)], 0.));
        assert!(p.geometry().contains([3., 0., z(3.)], 0.));
        assert!(!p.geometry().contains([1., 0., z(3.)], 0.));
        assert!(p.geometry().contains([1., 0., z(5.)], 0.));
    }
}
#[test]
fn plain_pocket_is_distinct_and_failures_do_not_mutate_original_stock() {
    let original = stock();
    let mut plain = spec();
    plain.seat = None;
    let p = Pocket::new(
        original,
        Placement::Axial,
        entry(false),
        plain,
        1e-7,
        1e-5,
        None,
    )
    .unwrap();
    assert!(!p.roles().contains(&Role::Shoulder));
    assert!((p.geometry().volume_mm3() - 984. * std::f64::consts::PI).abs() < 1e-10);
    for s in [
        Spec {
            depth_mm: 9.5,
            ..spec()
        },
        Spec {
            minimum_wall_mm: 7.,
            ..spec()
        },
        Spec {
            seat: Some(Seat {
                radius_mm: 1.,
                depth_mm: 1.,
            }),
            ..spec()
        },
        Spec {
            radius_mm: f64::NAN,
            ..spec()
        },
    ] {
        assert!(
            Pocket::new(
                original,
                Placement::Axial,
                entry(false),
                s,
                1e-7,
                1e-5,
                None
            )
            .is_err()
        );
    }
    assert_eq!(original.height(), stock().height());
    assert!(
        Pocket::new(
            original,
            Placement::Axial,
            entry(false),
            spec(),
            1e-7,
            1e-5,
            Some(&|| true)
        )
        .is_err()
    );
}

#[test]
fn literal_proper_poses_and_clocking_keep_the_stock_chart_without_normalization() {
    let u = [0.8, 0.6, 0.];
    let v = [-0.36, 0.48, 0.8];
    let n = [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ];
    let frame = Frame {
        origin_mm: [17., -8., 13.],
        axes: [u, v, n],
    };
    for end in [false, true] {
        let entry = Frame {
            origin_mm: std::array::from_fn(|i| {
                frame.origin_mm[i] + if end { 10. * n[i] } else { 0. }
            }),
            axes: if end {
                [u, v.map(|x| -x), n.map(|x| -x)]
            } else {
                frame.axes
            },
        };
        let pocket = Pocket::new(
            stock(),
            Placement::Framed(frame),
            entry,
            spec(),
            1e-7,
            1e-5,
            None,
        )
        .unwrap();
        assert_eq!(pocket.frame(), frame);
        assert_eq!(pocket.end(), end);
    }
    let mut roll = entry(false);
    roll.axes[0] = [-1., 0., 0.];
    roll.axes[1] = [0., -1., 0.];
    assert!(Pocket::new(stock(), Placement::Axial, roll, spec(), 1e-7, 1e-5, None).is_ok());
}
#[test]
fn source_and_placement_budgets_refuse_invalid_reflected_skew_offaxis_and_tight_cases() {
    let build = |placement, entry, tolerance| {
        Pocket::new(stock(), placement, entry, spec(), tolerance, 1e-5, None)
    };
    let mut shifted = entry(false);
    shifted.origin_mm[0] = 0.01;
    assert!(build(Placement::Axial, shifted, 1e-7).is_err());
    let mut scaled = entry(false);
    scaled.axes[0][0] = 2.;
    assert!(build(Placement::Axial, scaled, 1e-7).is_err());
    let mut reflected = entry(false);
    reflected.axes[1][1] = -1.;
    assert!(build(Placement::Axial, reflected, 1e-7).is_err());
    let mut nonfinite = entry(false);
    nonfinite.origin_mm[0] = f64::INFINITY;
    assert!(build(Placement::Axial, nonfinite, 1e-7).is_err());
    let mut mismatched = entry(false);
    mismatched.axes[2][0] = 1e-12;
    assert!(build(Placement::Framed(mismatched), entry(false), 1e-7).is_err());
    assert!(build(Placement::Axial, entry(false), 1e-30).is_err());
    assert!(build(Placement::Axial, entry(false), 1e-4).is_err());
    let calls = std::cell::Cell::new(0);
    let cancel = || {
        let n = calls.get() + 1;
        calls.set(n);
        n > 1
    };
    assert_eq!(
        Pocket::new(
            stock(),
            Placement::Axial,
            entry(false),
            spec(),
            1e-7,
            1e-5,
            Some(&cancel)
        )
        .unwrap_err(),
        Error::Cancelled
    );
}

#[test]
fn nonbinary_exact_stock_conversion_and_physical_frame_margins_are_explicit() {
    let decimal = |x: &str| Rational::parse_decimal(x).unwrap();
    let stock = AxialCylinder::new(
        99,
        Circle::new([decimal("1.1"), decimal("-2.2")], decimal("10.1")).unwrap(),
        decimal("17.3"),
        decimal("10.1"),
    )
    .unwrap();
    let entry = Frame {
        origin_mm: [1.1, -2.2, 27.4],
        ..entry(true)
    };
    let result = Pocket::new(stock, Placement::Axial, entry, spec(), 1e-7, 1e-5, None).unwrap();
    assert_eq!(result.stock().circle().radius(), decimal("10.1"));
    assert_eq!(result.frame().origin_mm, [1.1, -2.2, 17.3]);
    assert!(result.source_error_mm() > 0.);
    assert!(Pocket::new(stock, Placement::Axial, entry, spec(), 1e-15, 1e-5, None).is_err());
    let frame = Frame {
        origin_mm: [0.; 3],
        axes: [[0.9999, 0., 0.], [0., 1., 0.], [0., 0., 0.9999]],
    };
    let s = Spec {
        seat: Some(Seat {
            radius_mm: 8.9999,
            depth_mm: 2.,
        }),
        ..spec()
    };
    assert!(
        Pocket::new(
            super::tests::stock(),
            Placement::Framed(frame),
            frame,
            s,
            1e-7,
            1e-5,
            None
        )
        .is_err()
    );
    let wide = Frame {
        origin_mm: [0.; 3],
        axes: [[1., 0., 0.], [0., 1., 1e-5], [0., -1e-5, 1.]],
    };
    assert!(
        Pocket::new(
            super::tests::stock(),
            Placement::Axial,
            wide,
            Spec {
                depth_mm: 0.001,
                seat: Some(Seat {
                    radius_mm: 4.,
                    depth_mm: 0.0001
                }),
                ..spec()
            },
            1e-7,
            1e-5,
            None
        )
        .is_err()
    );
}

#[test]
fn finite_carrier_resource_and_unresolved_tiny_features_refuse_without_clamping() {
    let q = |v| Rational::new(v, 1).unwrap();
    let large = AxialCylinder::new(
        7,
        Circle::new([q(0); 2], q(100_000_000)).unwrap(),
        q(0),
        q(10),
    )
    .unwrap();
    assert_eq!(
        Pocket::new(
            large,
            Placement::Axial,
            entry(false),
            spec(),
            0.01,
            1.,
            None
        )
        .unwrap_err(),
        Error::Limit
    );
    let broad = AxialCylinder::new(
        7,
        Circle::new([q(0); 2], q(1_000_000)).unwrap(),
        q(0),
        q(10),
    )
    .unwrap();
    assert_eq!(
        Pocket::new(
            broad,
            Placement::Axial,
            entry(false),
            Spec {
                seat: Some(Seat {
                    radius_mm: 4.,
                    depth_mm: 2e-9
                }),
                ..spec()
            },
            1.1e-9,
            1.2e-9,
            None
        )
        .unwrap_err(),
        Error::Unresolved
    );
    assert!(
        Pocket::new(
            stock(),
            Placement::Axial,
            entry(false),
            Spec {
                seat: Some(Seat {
                    radius_mm: 4.,
                    depth_mm: 1e-5
                }),
                ..spec()
            },
            1e-7,
            1e-5,
            None
        )
        .is_err()
    );
    let positive = Pocket::new(
        stock(),
        Placement::Axial,
        entry(false),
        Spec {
            seat: Some(Seat {
                radius_mm: 4.,
                depth_mm: 2e-5,
            }),
            ..spec()
        },
        1e-7,
        1e-5,
        None,
    )
    .unwrap();
    assert!(positive.roles().contains(&Role::SeatWall));
}

#[test]
fn independently_admissible_affine_charts_can_still_fail_relative_circle_qualification() {
    let lower = ((4_f64 - 0.00999).sqrt() - 1.).sqrt();
    let upper = ((4_f64 + 0.00999).sqrt() - 1.).sqrt();
    let frame = |a: f64| Frame {
        origin_mm: [0.; 3],
        axes: [[a, 0., 0.], [0., a, 0.], [0., 0., a * a]],
    };
    assert!(metric(frame(lower)).is_ok() && metric(frame(upper)).is_ok());
    assert_eq!(
        Pocket::new(
            stock(),
            Placement::Framed(frame(lower)),
            frame(upper),
            spec(),
            1e-7,
            1e-5,
            None
        )
        .unwrap_err(),
        Error::Unsupported
    );
}
