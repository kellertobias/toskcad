use super::*;
use crate::orthogonal::{self, Box3, Operation};
fn identity() -> Frame {
    Frame {
        origin_mm: [0.; 3],
        axes: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
    }
}
fn req() -> PassageRequirements {
    PassageRequirements {
        minimum_radial_clearance_mm: 1e-7,
        linear_tolerance_mm: 1e-7,
        minimum_feature_mm: 1e-7,
    }
}
fn stock() -> (OrthogonalSolid, OrthogonalThroughBore) {
    let stock = orthogonal::boolean(
        &[Box3 {
            min: [-10., -10., -4.1],
            max: [10., 10., 0.],
        }],
        &[],
        Operation::Union,
        0.,
    )
    .unwrap()
    .unwrap();
    let boundary = orthogonal_boundary::build(stock.quads()).unwrap();
    let bore = orthogonal_pockets::identify_through_bore(
        &stock,
        &boundary,
        [0., 0., -5.],
        [0., 0., 1.],
        1.7,
        1e-7,
        1e-7,
    )
    .unwrap();
    (stock, bore)
}
fn fixture_source(left: bool, nb: bool) -> Hanger {
    Hanger::new(crate::authored_hanger_bolt::tests::fixture(left, nb), None).unwrap()
}
#[test]
fn actual_clamp_and_washer_voids_clear_whole_prefix_and_exact_wood_interface() {
    for (left, nb, posed) in [
        (false, false, false),
        (true, true, false),
        (false, true, true),
        (true, false, true),
    ] {
        let (stock, bore) = stock();
        let source = fixture_source(left, nb);
        let frame = if posed {
            Frame {
                origin_mm: [11., -7., 5.],
                axes: [[0.36, 0.48, 0.8], [-0.8, 0.6, 0.], [-0.48, -0.64, 0.6]],
            }
        } else {
            identity()
        };
        let p = certify_hanger_through_stock(
            ThroughStock {
                material: &stock,
                bore,
                frame,
            },
            &source,
            frame,
            1.,
            req(),
            &|| false,
        )
        .unwrap();
        assert_eq!(p.wood_side, WoodSide::AfterEnd);
        assert_eq!(p.source_cap_annular_land_mm[0], bore.radius_mm);
        let metric = scale(frame).unwrap();
        assert!(
            q(p.source_cap_annular_land_mm[1]).unwrap() - q(bore.radius_mm).unwrap()
                >= q(1.).unwrap() / q(metric.minimum_scale.lo).unwrap()
        );
        assert_eq!(p.wood_envelope_plane_gap_mm, [0., 0.]);
        assert!(p.radial_clearance_lower_mm > 0.19);
        let ring = Annulus::new(1.7, 3.5, 0.2).unwrap();
        let mut ring_frame = frame;
        ring_frame.origin_mm = std::array::from_fn(|i| frame.origin_mm[i] - 4.3 * frame.axes[2][i]);
        let r =
            certify_hanger_annulus(&ring, ring_frame, &source, frame, req(), &|| false).unwrap();
        assert_eq!(
            r.source_cap_annular_land_mm,
            [ring.inner_radius_mm(), ring.outer_radius_mm()]
        );
        assert!(r.wood_envelope_plane_gap_mm[0] > 4.09);
        assert!(r.radial_clearance_lower_mm > 0.19);
    }
}
#[test]
fn sub_tolerance_actual_interface_penetration_and_radial_interference_are_not_waived() {
    let (stock, bore) = stock();
    let source = fixture_source(false, false);
    for offset in [-1e-12, -f64::from_bits(1)] {
        let mut frame = identity();
        frame.origin_mm[2] = offset;
        assert!(matches!(
            certify_hanger_through_stock(
                ThroughStock {
                    material: &stock,
                    bore,
                    frame: identity()
                },
                &source,
                frame,
                1.,
                req(),
                &|| false
            ),
            Err(Error::Unsupported)
        ));
    }
    let mut off = identity();
    off.origin_mm[0] = 0.21;
    assert!(matches!(
        certify_hanger_through_stock(
            ThroughStock {
                material: &stock,
                bore,
                frame: identity()
            },
            &source,
            off,
            1.,
            req(),
            &|| false
        ),
        Err(Error::Unsupported)
    ));
    let ring = Annulus::new(1.4, 3.5, 0.2).unwrap();
    let mut f = identity();
    f.origin_mm[2] = -4.3;
    assert!(matches!(
        certify_hanger_annulus(&ring, f, &source, identity(), req(), &|| false),
        Err(Error::Unsupported)
    ));
}
#[test]
fn validation_wall_and_cancellation_keep_real_source_and_domains_unchanged() {
    use std::cell::Cell;
    let (stock, bore) = stock();
    let source = fixture_source(false, false);
    assert!(
        certify_hanger_through_stock(
            ThroughStock {
                material: &stock,
                bore,
                frame: identity()
            },
            &source,
            identity(),
            9.,
            req(),
            &|| false
        )
        .is_err()
    );
    let mut bad = req();
    bad.minimum_radial_clearance_mm = f64::NAN;
    assert!(matches!(
        certify_hanger_through_stock(
            ThroughStock {
                material: &stock,
                bore,
                frame: identity()
            },
            &source,
            identity(),
            1.,
            bad,
            &|| false
        ),
        Err(Error::Invalid)
    ));
    let count = Cell::new(0);
    certify_hanger_through_stock(
        ThroughStock {
            material: &stock,
            bore,
            frame: identity(),
        },
        &source,
        identity(),
        1.,
        req(),
        &|| {
            count.set(count.get() + 1);
            false
        },
    )
    .unwrap();
    for stop in 1..=count.get() {
        let calls = Cell::new(0);
        assert!(matches!(
            certify_hanger_through_stock(
                ThroughStock {
                    material: &stock,
                    bore,
                    frame: identity()
                },
                &source,
                identity(),
                1.,
                req(),
                &|| {
                    calls.set(calls.get() + 1);
                    calls.get() == stop
                }
            ),
            Err(Error::Cancelled)
        ));
    }
    assert_eq!(stock.cells().len(), 1);
    assert_eq!(source, fixture_source(false, false));
}

#[test]
fn reverse_and_other_cardinal_axes_preserve_real_source_plane_equality() {
    for axis in 0..3 {
        for reverse in [false, true] {
            let mut lower = [-10.; 3];
            let mut upper = [10.; 3];
            lower[axis] = if reverse { 0. } else { -4.1 };
            upper[axis] = if reverse { 4.1 } else { 0. };
            let stock = orthogonal::boolean(
                &[Box3 {
                    min: lower,
                    max: upper,
                }],
                &[],
                Operation::Union,
                0.,
            )
            .unwrap()
            .unwrap();
            let boundary = orthogonal_boundary::build(stock.quads()).unwrap();
            let mut start = [0.; 3];
            let mut end = [0.; 3];
            start[axis] = -5.;
            end[axis] = 5.;
            let bore = orthogonal_pockets::identify_through_bore(
                &stock, &boundary, start, end, 1.7, 1e-7, 1e-7,
            )
            .unwrap();
            let mut axes = [[0.; 3]; 3];
            axes[0][(axis + 1) % 3] = 1.;
            axes[1][(axis + 2) % 3] = if reverse { -1. } else { 1. };
            axes[2][axis] = if reverse { -1. } else { 1. };
            let placed = Frame {
                origin_mm: [0.; 3],
                axes,
            };
            let source = fixture_source(false, false);
            let result = certify_hanger_through_stock(
                ThroughStock {
                    material: &stock,
                    bore,
                    frame: identity(),
                },
                &source,
                placed,
                1.,
                req(),
                &|| false,
            )
            .unwrap();
            assert_eq!(
                result.wood_side,
                if reverse {
                    WoodSide::BeforeStart
                } else {
                    WoodSide::AfterEnd
                }
            );
            assert_eq!(result.wood_envelope_plane_gap_mm, [0., 0.]);
        }
    }
}
#[test]
fn exact_oblique_envelope_tangency_is_proved_without_near_plane_fitting() {
    let source = fixture_source(false, false);
    let ring = Annulus::new(10., 12., 4.1).unwrap();
    let c = 3_f64.sqrt() / 2.;
    let placed = Frame {
        origin_mm: [0., 0., 5.6],
        axes: [[c, 0., -0.5], [0., 1., 0.], [0.5, 0., c]],
    };
    let proof =
        certify_hanger_annulus(&ring, identity(), &source, placed, req(), &|| false).unwrap();
    assert_eq!(proof.wood_envelope_plane_gap_mm, [0., 0.]);
    let mut below = placed;
    below.origin_mm[2] = below.origin_mm[2].next_down();
    assert!(matches!(
        certify_hanger_annulus(&ring, identity(), &source, below, req(), &|| false),
        Err(Error::Unsupported)
    ));
}
#[test]
fn a_valid_bore_does_not_qualify_unchecked_stock_beyond_its_finite_slab() {
    let source = fixture_source(false, false);
    let stock = orthogonal::boolean(
        &[
            Box3 {
                min: [-10., -10., -4.1],
                max: [10., 10., 0.],
            },
            Box3 {
                min: [8., -2., 0.],
                max: [10., 2., 2.],
            },
        ],
        &[],
        Operation::Union,
        0.,
    )
    .unwrap()
    .unwrap();
    let boundary = orthogonal_boundary::build(stock.quads()).unwrap();
    let bore = orthogonal_pockets::identify_through_bore(
        &stock,
        &boundary,
        [0., 0., -5.],
        [0., 0., 3.],
        1.7,
        1e-7,
        1e-7,
    )
    .unwrap();
    assert!(matches!(
        certify_hanger_through_stock(
            ThroughStock {
                material: &stock,
                bore,
                frame: identity()
            },
            &source,
            identity(),
            1.,
            req(),
            &|| false
        ),
        Err(Error::Unsupported)
    ));
}

#[test]
fn exact_radial_threshold_and_real_finite_strained_maps_are_not_nominal_fits() {
    let (stock, bore) = stock();
    let source = fixture_source(false, false);
    let mut exact = req();
    exact.minimum_radial_clearance_mm = 1.7 - 1.5;
    let proof = certify_hanger_through_stock(
        ThroughStock {
            material: &stock,
            bore,
            frame: identity(),
        },
        &source,
        identity(),
        1.,
        exact,
        &|| false,
    )
    .unwrap();
    assert!(proof.radial_clearance_lower_mm >= exact.minimum_radial_clearance_mm);
    exact.minimum_radial_clearance_mm = exact.minimum_radial_clearance_mm.next_up();
    assert!(matches!(
        certify_hanger_through_stock(
            ThroughStock {
                material: &stock,
                bore,
                frame: identity()
            },
            &source,
            identity(),
            1.,
            exact,
            &|| false
        ),
        Err(Error::Unsupported)
    ));
    let mut common = identity();
    common.axes[0][0] = common.axes[0][0].next_up();
    let proof = certify_hanger_through_stock(
        ThroughStock {
            material: &stock,
            bore,
            frame: common,
        },
        &source,
        common,
        1.,
        req(),
        &|| false,
    )
    .unwrap();
    assert!(proof.radial_clearance_lower_mm > 0.19);
    let s = 1e-8;
    let c = (1_f64 - s * s).sqrt();
    let placed = Frame {
        origin_mm: [0.05, -0.02, 1e-6],
        axes: [[c, 0., -s], [0., 1., 0.], [s, 0., c]],
    };
    let proof = certify_hanger_through_stock(
        ThroughStock {
            material: &stock,
            bore,
            frame: identity(),
        },
        &source,
        placed,
        1.,
        req(),
        &|| false,
    )
    .unwrap();
    assert!(proof.radial_clearance_lower_mm > 0.12);
    assert!(proof.wood_envelope_plane_gap_mm[0] > 9e-7);
}
#[test]
fn annular_width_uses_original_exact_dimensions_and_cancellation_is_atomic() {
    use std::cell::Cell;
    let source = fixture_source(false, false);
    let ring = Annulus::new(1.8, 6., 5.).unwrap();
    let mut frame = identity();
    frame.origin_mm[2] = -9.;
    let mut limit = req();
    limit.minimum_feature_mm = 6. - 1.8;
    assert!(q(6.).unwrap() - q(1.8).unwrap() < q(limit.minimum_feature_mm).unwrap());
    assert!(matches!(
        certify_hanger_annulus(&ring, frame, &source, identity(), limit, &|| false),
        Err(Error::Unsupported)
    ));
    let thin = Annulus::new(1.7, 3.5, 1e-8).unwrap();
    assert!(matches!(
        certify_hanger_annulus(&thin, frame, &source, identity(), req(), &|| false),
        Err(Error::Unsupported)
    ));
    let calls = Cell::new(0);
    certify_hanger_annulus(&ring, frame, &source, identity(), req(), &|| {
        calls.set(calls.get() + 1);
        false
    })
    .unwrap();
    for stop in 1..=calls.get() {
        let count = Cell::new(0);
        assert!(matches!(
            certify_hanger_annulus(&ring, frame, &source, identity(), req(), &|| {
                count.set(count.get() + 1);
                count.get() == stop
            }),
            Err(Error::Cancelled)
        ));
    }
}
#[test]
fn invalid_source_domains_disconnected_stock_and_real_capacity_refuse() {
    let (stock, bore) = stock();
    let source = fixture_source(false, false);
    for wall in [0., f64::NAN] {
        assert!(matches!(
            certify_hanger_through_stock(
                ThroughStock {
                    material: &stock,
                    bore,
                    frame: identity()
                },
                &source,
                identity(),
                wall,
                req(),
                &|| false
            ),
            Err(Error::Invalid)
        ));
    }
    let mut broken = bore;
    broken.faces[0] = usize::MAX;
    assert!(matches!(
        certify_hanger_through_stock(
            ThroughStock {
                material: &stock,
                bore: broken,
                frame: identity()
            },
            &source,
            identity(),
            1.,
            req(),
            &|| false
        ),
        Err(Error::Invalid)
    ));
    let disconnected = orthogonal::boolean(
        &[
            Box3 {
                min: [-10., -10., -4.1],
                max: [10., 10., 0.],
            },
            Box3 {
                min: [20., 20., -4.1],
                max: [21., 21., 0.],
            },
        ],
        &[],
        Operation::Union,
        0.,
    )
    .unwrap()
    .unwrap();
    assert!(matches!(
        certify_hanger_through_stock(
            ThroughStock {
                material: &disconnected,
                bore,
                frame: identity()
            },
            &source,
            identity(),
            1.,
            req(),
            &|| false
        ),
        Err(Error::Unsupported)
    ));
    let mut boxes = vec![Box3 {
        min: [-10., -10., -4.1],
        max: [10., -8., 0.],
    }];
    boxes.extend((0..1400).map(|k| {
        let x = -7. + f64::from(k) * 0.01;
        Box3 {
            min: [x, -8., -4.1],
            max: [x + 0.005, 10., 0.],
        }
    }));
    let complex = orthogonal::boolean(&boxes, &[], Operation::Union, 0.)
        .unwrap()
        .unwrap();
    assert!(complex.cells().len() > 4096);
    assert!(matches!(
        certify_hanger_through_stock(
            ThroughStock {
                material: &complex,
                bore,
                frame: identity()
            },
            &source,
            identity(),
            1.,
            req(),
            &|| false
        ),
        Err(Error::Limit)
    ));
}

#[test]
fn one_sided_exact_enclosures_and_earned_intersections_preserve_large_valid_bounds() {
    use num_bigint::BigInt;
    let maximum = q(f64::MAX).unwrap();
    for value in [
        q(1.).unwrap(),
        Q::new(1.into(), 3.into()),
        Q::new((-1).into(), 3.into()),
        &maximum - q(1.).unwrap(),
        -&maximum + q(1.).unwrap(),
        q(f64::from_bits(1)).unwrap() / q(2.).unwrap(),
    ] {
        let interval = directed(&value).unwrap();
        assert!(interval.checked().is_ok());
        assert!(q(interval.lo).unwrap() <= value && q(interval.hi).unwrap() >= value);
    }
    for value in [
        &maximum + q(1.).unwrap(),
        -&maximum - q(1.).unwrap(),
        Q::from_integer(BigInt::from(1_u8) << 20000),
    ] {
        assert!(matches!(directed(&value), Err(Error::Unsupported)));
    }
    let radius = q(1.5).unwrap();
    let square = q(1.).unwrap();
    let center = &maximum - q(2.).unwrap();
    // Independent exact inequality earns the cap used by this numerical
    // handoff; no manufactured stock/bore or geometry receipt is fabricated.
    assert!((&maximum - &center).pow(2) >= (&radius * &radius * &square));
    let interval = prefix_interval(&center, &radius, &square, &maximum).unwrap();
    let actual = center + radius;
    assert!(interval.checked().is_ok());
    assert!(q(interval.lo).unwrap() <= actual && q(interval.hi).unwrap() >= actual);
}
#[test]
fn actual_extreme_stock_plane_support_recovers_finite_bound_but_world_overflow_refuses() {
    let stock = orthogonal::boolean(
        &[Box3 {
            min: [0., 0., f64::MAX.next_down()],
            max: [1e-100, 1e-100, f64::MAX],
        }],
        &[],
        Operation::Union,
        0.,
    )
    .unwrap()
    .unwrap();
    assert!(stock.volume_mm3().is_finite());
    let source = fixture_source(false, false);
    let span = [
        q(stock.bounds().min[2]).unwrap(),
        q(stock.bounds().max[2]).unwrap(),
    ];
    let domain = Domain {
        axis: 2,
        center: [Q::zero(), Q::zero(), Q::zero()],
        span,
        radius: q(1.7).unwrap(),
        frame: identity(),
    };
    let (o, a) = affine(identity(), identity(), &|| false).unwrap();
    let (side, actual) =
        wood_separation(&source, &o, &a, &domain, &scale(identity()).unwrap()).unwrap();
    assert_eq!(side, WoodSide::BeforeStart);
    assert!(actual.checked().is_ok());
    let true_gap = q(stock.bounds().min[2]).unwrap() - q(source.wood().full_length_mm()).unwrap();
    assert!(q(actual.lo).unwrap() <= true_gap && q(actual.hi).unwrap() >= true_gap);
    // This tests only the actual slab-plane support primitive. The tiny XY
    // face is deliberately NOT claimed to admit a through bore or a stack.
    let s = 1.00000001;
    let scaled = Frame {
        origin_mm: [0.; 3],
        axes: [[s, 0., 0.], [0., s, 0.], [0., 0., s * s]],
    };
    let metric = scale(scaled).unwrap();
    assert!(matches!(
        wood_separation(&source, &o, &a, &domain, &metric),
        Err(Error::Unsupported)
    ));
}
