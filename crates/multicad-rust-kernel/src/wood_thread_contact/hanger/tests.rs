use super::*;
fn identity() -> Frame {
    Frame {
        origin_mm: [0.; 3],
        axes: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
    }
}
fn run(
    stock: &OrthogonalSolid,
    pilot: OrthogonalPocket,
    stock_frame: Frame,
    source: &Hanger,
    frame: Frame,
    req: Requirements,
) -> R<HangerContact> {
    certify_hanger_wood(
        Stock {
            material: stock,
            pilot,
            frame: stock_frame,
        },
        source,
        frame,
        req,
        &|| false,
    )
}
#[test]
fn original_hanger_and_actual_pilot_prove_positive_material_bite_for_both_hands_and_proper_pose() {
    for (left, nb, posed) in [
        (false, false, false),
        (true, true, false),
        (false, true, true),
        (true, false, true),
    ] {
        let (material, pilot, _, mut frame, req) = super::super::tests::fixture(left, 1);
        let source =
            Hanger::new(crate::authored_hanger_bolt::tests::fixture(left, nb), None).unwrap();
        let stock_frame = if posed {
            Frame {
                origin_mm: [11., -7., 5.],
                axes: [[0.36, 0.48, 0.8], [-0.8, 0.6, 0.], [-0.48, -0.64, 0.6]],
            }
        } else {
            identity()
        };
        if posed {
            frame = Frame {
                origin_mm: std::array::from_fn(|i| {
                    stock_frame.origin_mm[i]
                        + (0..3)
                            .map(|j| stock_frame.axes[j][i] * frame.origin_mm[j])
                            .sum::<f64>()
                }),
                axes: std::array::from_fn(|j| {
                    std::array::from_fn(|i| {
                        (0..3)
                            .map(|k| stock_frame.axes[k][i] * frame.axes[j][k])
                            .sum()
                    })
                }),
            };
        }
        let proof = run(&material, pilot, stock_frame, &source, frame, req).unwrap();
        assert!(proof.engagement_mm[0] > 11.9);
        assert!(proof.tip_clearance_mm[0] > 3.6);
        assert!(proof.witness_volume_lower_mm3 > 0.);
        assert_eq!(
            source.classify_box(proof.thread_witness_mm, None).unwrap(),
            Classification::Inside
        );
    }
}
#[test]
fn tolerated_intent_offsets_are_measured_and_independently_safe_not_renamed_exact_planes() {
    let (material, pilot, _, frame, req) = super::super::tests::fixture(false, 1);
    let source = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, false),
        None,
    )
    .unwrap();
    for shift in [-5e-8, 5e-8] {
        let mut actual = frame;
        actual.origin_mm[2] += shift;
        actual.origin_mm[0] += 5e-8;
        let p = run(&material, pilot, identity(), &source, actual, req).unwrap();
        assert!(p.source_eccentricity_mm[0] > 0.);
        assert!(p.source_axial_offset_mm[0] > 0. || p.source_axial_offset_mm[1] < 0.);
    }
    let mut off = frame;
    off.origin_mm[0] += 0.01;
    assert!(run(&material, pilot, identity(), &source, off, req).is_err());
    let mut deep = frame;
    deep.origin_mm[2] -= 0.01;
    assert!(run(&material, pilot, identity(), &source, deep, req).is_err());
}
#[test]
fn actual_wall_floor_prefix_alignment_and_pitch_domain_failures_never_publish_contact() {
    let (material, pilot, _, frame, req) = super::super::tests::fixture(false, 1);
    let source = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, false),
        None,
    )
    .unwrap();
    for case in 0..10 {
        let (mut p, mut f, mut r) = (pilot, frame, req);
        match case {
            0 => r.minimum_engagement_mm = 13.,
            1 => r.tip_clearance_mm = 5.,
            2 => r.minimum_floor_mm = 5.,
            3 => r.minimum_wall_mm = 11.,
            4 => p.radius_mm = 1.,
            5 => p.radius_mm = 3.,
            6 => f.axes = identity().axes,
            7 => f.axes[0][1] = 1e-8,
            8 => p.depth_mm = 10.,
            _ => f.origin_mm[2] += 1.,
        }
        assert!(
            run(&material, p, identity(), &source, f, r).is_err(),
            "case{case}"
        );
    }
    for value in [0., -1., f64::NAN, 1e8] {
        let mut r = req;
        r.minimum_wall_mm = value;
        assert!(matches!(
            run(&material, pilot, identity(), &source, frame, r),
            Err(Error::Invalid)
        ));
    }
    let mut r = req;
    r.linear_tolerance_mm = 1e-8;
    assert!(matches!(
        run(&material, pilot, identity(), &source, frame, r),
        Err(Error::Invalid)
    ));
    let mut f = frame;
    f.origin_mm[0] = f64::NAN;
    assert!(matches!(
        run(&material, pilot, identity(), &source, f, req),
        Err(Error::Invalid)
    ));
    let mut f = frame;
    f.axes[0][0] = 0.;
    assert!(matches!(
        run(&material, pilot, identity(), &source, f, req),
        Err(Error::Unsupported)
    ));
    let mut p = pilot;
    p.radius_mm = 3_f64.next_down();
    assert!(matches!(
        run(&material, p, identity(), &source, frame, req),
        Err(Error::Unsupported)
    ));
    assert!(run(&material, pilot, identity(), &source, frame, req).is_ok());
}
#[test]
fn actual_stock_overhang_disconnection_capacity_and_every_checkpoint_refuse_without_mutation() {
    use crate::orthogonal::{self, Box3, Operation};
    let (material, pilot, _, frame, req) = super::super::tests::fixture(false, 1);
    let source = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, false),
        None,
    )
    .unwrap();
    let boxes = [
        Box3 {
            min: [0.; 3],
            max: [20.; 3],
        },
        Box3 {
            min: [0., 0., 20.],
            max: [1., 1., 22.],
        },
    ];
    let overhang = orthogonal::boolean(&boxes, &[], Operation::Union, 1e-7)
        .unwrap()
        .unwrap();
    let b = orthogonal_boundary::build(overhang.quads()).unwrap();
    let p = orthogonal_pockets::identify_pocket(
        &overhang,
        &b,
        [10., 10., 20.001],
        [10., 10., 4.],
        2.,
        1e-7,
        1e-4,
    )
    .unwrap();
    assert!(matches!(
        run(&overhang, p, identity(), &source, frame, req),
        Err(Error::Unsupported)
    ));
    let separate = orthogonal::boolean(
        &[
            boxes[0],
            Box3 {
                min: [30., 0., 0.],
                max: [31., 1., 1.],
            },
        ],
        &[],
        Operation::Union,
        1e-7,
    )
    .unwrap()
    .unwrap();
    assert!(matches!(
        run(&separate, pilot, identity(), &source, frame, req),
        Err(Error::Unsupported)
    ));
    let mut comb = vec![Box3 {
        min: [0.; 3],
        max: [2800., 1., 1.],
    }];
    comb.extend((0..1400).map(|i| Box3 {
        min: [f64::from(2 * i), 1., 0.],
        max: [f64::from(2 * i + 1), 2., 1.],
    }));
    let large = orthogonal::boolean(&comb, &[], Operation::Union, 1e-7)
        .unwrap()
        .unwrap();
    assert!(large.cells().len() > 4096);
    assert!(matches!(
        run(&large, pilot, identity(), &source, frame, req),
        Err(Error::Limit)
    ));
    let count = std::cell::Cell::new(0);
    certify_hanger_wood(
        Stock {
            material: &material,
            pilot,
            frame: identity(),
        },
        &source,
        frame,
        req,
        &|| {
            count.set(count.get() + 1);
            false
        },
    )
    .unwrap();
    let total = count.get();
    for at in [1, 4, total - 1, total] {
        count.set(0);
        assert!(matches!(
            certify_hanger_wood(
                Stock {
                    material: &material,
                    pilot,
                    frame: identity()
                },
                &source,
                frame,
                req,
                &|| {
                    count.set(count.get() + 1);
                    count.get() == at
                }
            ),
            Err(Error::Cancelled)
        ));
    }
    assert!(run(&material, pilot, identity(), &source, frame, req).is_ok());
}
#[test]
fn a_real_hanger_core_box_inside_the_pilot_void_is_not_a_wood_material_witness() {
    let (_, pilot, _, frame, _) = super::super::tests::fixture(false, 1);
    let source = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, false),
        None,
    )
    .unwrap();
    let witness = [[-0.01, 0.01], [-0.01, 0.01], [5., 5.01]];
    assert_eq!(
        source.classify_box(witness, None).unwrap(),
        Classification::Inside
    );
    let origin = frame.origin_mm.map(|x| q(x).unwrap());
    let columns = [[1, 0, 0], [0, -1, 0], [0, 0, -1]];
    assert_eq!(
        material_witness(
            &source,
            pilot,
            witness,
            WitnessPlacement {
                origin: &origin,
                columns: &columns,
                axis: 2,
                sign: -1.,
                determinant: ip(1.)
            },
            &|| false
        )
        .unwrap(),
        None
    );
    let third = num_rational::BigRational::new(1.into(), 3.into());
    let e = exact_interval(&third).unwrap();
    assert!(q(e.lo).unwrap() <= third && q(e.hi).unwrap() >= third);
    let huge = num_rational::BigRational::from_integer(num_bigint::BigInt::from(1) << 2048);
    assert!(exact_interval(&huge).is_err());
    for (source, target) in [
        (crate::authored_hanger_bolt::Error::Invalid, Error::Invalid),
        (crate::authored_hanger_bolt::Error::Limit, Error::Limit),
        (
            crate::authored_hanger_bolt::Error::Cancelled,
            Error::Cancelled,
        ),
        (
            crate::authored_hanger_bolt::Error::Unsupported,
            Error::Unsupported,
        ),
        (
            crate::authored_hanger_bolt::Error::Unresolved,
            Error::Unsupported,
        ),
    ] {
        assert_eq!(source_error(source), target);
    }
}

#[test]
fn exact_origin_expression_has_a_static_full_binary64_resource_bound() {
    let large = f64::MAX;
    let small = f64::from_bits(1);
    let delta = finite_origin_delta(
        [large, -large, small],
        [small, large, -small],
        [
            [large, small, -large],
            [small, large, small],
            [-large, -small, large],
        ],
        [large, small, -large],
    );
    assert!(
        delta
            .iter()
            .all(|q| q.numer().bits() + q.denom().bits() < 12710)
    );
    let square = dot(&delta, &delta).unwrap();
    assert!(square.numer().bits() + square.denom().bits() < 12710);
    assert!(square > Q::zero());
}
