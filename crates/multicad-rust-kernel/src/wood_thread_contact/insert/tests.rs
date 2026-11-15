use super::*;
use crate::{
    authored_flat_thread::{AuthoredThread, SourceOptions},
    helical_thread::Spec,
};
fn fixture(
    left: bool,
    starts: u32,
) -> (
    OrthogonalSolid,
    OrthogonalPocket,
    Insert,
    Frame,
    Requirements,
) {
    let (material, mut pilot, _, frame, mut req) = super::super::tests::fixture(left, starts);
    pilot.radius_mm = 4.;
    pilot.depth_mm = 6.;
    let source = |radius, pitch, depth, crest| {
        AuthoredThread::new(
            Spec {
                major_radius_mm: radius,
                pitch_mm: pitch,
                length_mm: 4.,
                depth_mm: depth,
                crest_width_mm: crest,
                lower_flank_deg: 30.,
                upper_flank_deg: 30.,
                starts,
                left_hand: left,
                root_radius_mm: None,
            },
            SourceOptions {
                position_limit_mm: 1e-8,
                normal_angle_limit_rad: 1e-6,
            },
        )
        .unwrap()
    };
    let insert = Insert::new(
        source(5., 2., 1., 0.4),
        source(1.5, 0.5, 0.25, 0.05),
        0.5,
        None,
    )
    .unwrap();
    req.minimum_engagement_mm = 2.;
    (material, pilot, insert, frame, req)
}
#[test]
fn original_authored_pairs_prove_real_stock_bite_without_headed_alias() {
    for left in [false, true] {
        for starts in [1, 2] {
            let (material, pilot, insert, frame, req) = fixture(left, starts);
            let proof = certify_insert(
                Stock {
                    material: &material,
                    pilot,
                    frame: Frame {
                        origin_mm: [0.; 3],
                        axes: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
                    },
                },
                &insert,
                frame,
                req,
                &|| false,
            )
            .unwrap();
            assert!(proof.engagement_mm[0] <= 4. && proof.engagement_mm[1] >= 4.);
            assert!(proof.tip_clearance_mm[0] > 1.99);
            assert!(proof.witness_volume_lower_mm3 > 0.);
            assert_eq!(
                insert.classify_box(proof.thread_witness_mm, None).unwrap(),
                Classification::Inside
            );
        }
    }
}
fn identity() -> Frame {
    Frame {
        origin_mm: [0.; 3],
        axes: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
    }
}
fn run(
    material: &OrthogonalSolid,
    pilot: OrthogonalPocket,
    insert: &Insert,
    frame: Frame,
    req: Requirements,
) -> R<InsertContact> {
    certify_insert(
        Stock {
            material,
            pilot,
            frame: identity(),
        },
        insert,
        frame,
        req,
        &|| false,
    )
}
#[test]
fn actual_fit_planes_walls_floors_and_source_owned_witnesses_cannot_be_waived() {
    let (material, pilot, insert, frame, req) = fixture(false, 1);
    for case in 0..13 {
        let (mut p, mut f, mut r) = (pilot, frame, req);
        match case {
            0 => p.radius_mm = 4_f64.next_down(),
            1 => p.radius_mm = 5.,
            2 => r.minimum_engagement_mm = 5.,
            3 => r.tip_clearance_mm = 3.,
            4 => r.minimum_wall_mm = 8.,
            5 => r.minimum_floor_mm = 15.,
            6 => f.origin_mm[0] = 10_f64.next_up(),
            7 => f.origin_mm[2] = 20_f64.next_up(),
            8 => f.origin_mm[2] = 20_f64.next_down(),
            9 => f.axes[0][0] = 1. + 1e-12,
            10 => f.axes[2][2] = 1.,
            11 => p.depth_mm = 3.,
            _ => p.face = usize::MAX,
        }
        assert!(run(&material, p, &insert, f, r).is_err(), "case{case}");
    }
    let mut p = pilot;
    p.radius_mm = 5_f64.next_down();
    assert_eq!(
        run(&material, p, &insert, frame, req).unwrap_err(),
        Error::Unsupported
    );
    let mut p = pilot;
    p.radius_mm = 4.25;
    let mut shifted = frame;
    shifted.origin_mm[0] = 10_f64.next_up();
    assert!(run(&material, p, &insert, shifted, req).is_ok());
    let mut p = pilot;
    p.radius_mm = 4.8;
    let mut eccentric = frame;
    eccentric.origin_mm[0] -= 0.14;
    eccentric.origin_mm[1] += 0.04;
    assert_eq!(
        run(&material, p, &insert, eccentric, req).unwrap_err(),
        Error::Unsupported
    );
    for value in [f64::NAN, -1., 0., 1e8] {
        let mut bad = req;
        bad.minimum_wall_mm = value;
        assert_eq!(
            run(&material, pilot, &insert, frame, bad).unwrap_err(),
            Error::Invalid
        );
    }
    let mut bad = req;
    bad.linear_tolerance_mm = 1e-8;
    assert_eq!(
        run(&material, pilot, &insert, frame, bad).unwrap_err(),
        Error::Invalid
    );
    assert!(run(&material, pilot, &insert, frame, req).is_ok());
}
#[test]
fn cancellation_and_actual_material_capacity_are_bounded() {
    use crate::orthogonal::{self, Box3, Operation};
    use std::cell::Cell;
    let (material, pilot, insert, frame, req) = fixture(false, 1);
    let calls = Cell::new(0);
    certify_insert(
        Stock {
            material: &material,
            pilot,
            frame: identity(),
        },
        &insert,
        frame,
        req,
        &|| {
            calls.set(calls.get() + 1);
            false
        },
    )
    .unwrap();
    for target in 1..=calls.get() {
        let n = Cell::new(0);
        assert_eq!(
            certify_insert(
                Stock {
                    material: &material,
                    pilot,
                    frame: identity()
                },
                &insert,
                frame,
                req,
                &|| {
                    n.set(n.get() + 1);
                    n.get() >= target
                }
            )
            .unwrap_err(),
            Error::Cancelled
        );
    }
    let mut cells = vec![Box3 {
        min: [0.; 3],
        max: [2800., 1., 1.],
    }];
    cells.extend((0..1400).map(|i| Box3 {
        min: [f64::from(i * 2), 1., 0.],
        max: [f64::from(i * 2 + 1), 2., 1.],
    }));
    let large = orthogonal::boolean(&cells, &[], Operation::Union, 1e-7)
        .unwrap()
        .unwrap();
    assert!(large.cells().len() > 4096);
    assert_eq!(
        run(&large, pilot, &insert, frame, req).unwrap_err(),
        Error::Limit
    );
    let split = orthogonal::boolean(
        &[
            Box3 {
                min: [0.; 3],
                max: [20.; 3],
            },
            Box3 {
                min: [30.; 3],
                max: [40.; 3],
            },
        ],
        &[],
        Operation::Union,
        1e-7,
    )
    .unwrap()
    .unwrap();
    assert_eq!(
        run(&split, pilot, &insert, frame, req).unwrap_err(),
        Error::Unsupported
    );
}
#[test]
fn actual_common_frames_preserve_metric_lengths_and_original_inner_authority() {
    let (material, mut pilot, insert, mut local, req) = fixture(true, 2);
    pilot.radius_mm = 4.25;
    let stock = Frame {
        origin_mm: [100., -20., 50.],
        axes: [[0.8, 0.6, 0.], [-0.6, 0.8, 0.], [0., 0., 1.]],
    };
    local.origin_mm = [102., -6., 70.];
    local.axes = [
        stock.axes[0],
        stock.axes[1].map(|x| -x),
        stock.axes[2].map(|x| -x),
    ];
    let proof = certify_insert(
        Stock {
            material: &material,
            pilot,
            frame: stock,
        },
        &insert,
        local,
        req,
        &|| false,
    )
    .unwrap();
    assert!(proof.engagement_mm[0] <= 4. && proof.engagement_mm[1] >= 4.);
    assert!(proof.witness_volume_lower_mm3 > 0.);
    let mut changed = insert.inner().spec();
    changed.left_hand = false;
    let inner = AuthoredThread::new(
        changed,
        SourceOptions {
            position_limit_mm: 1e-8,
            normal_angle_limit_rad: 1e-6,
        },
    )
    .unwrap();
    let paired = Insert::new(*insert.outer(), inner, insert.minimum_wall_mm(), None).unwrap();
    assert!(
        certify_insert(
            Stock {
                material: &material,
                pilot,
                frame: stock
            },
            &paired,
            local,
            req,
            &|| false
        )
        .is_ok()
    );
}
#[test]
fn six_real_stock_faces_retain_their_actual_inward_direction_and_floor() {
    let (material, _, insert, _, req) = fixture(false, 1);
    let boundary = orthogonal_boundary::build(material.quads()).unwrap();
    for axis in 0..3 {
        for sign in [-1., 1.] {
            let mut approach = [10.; 3];
            let mut end = [10.; 3];
            approach[axis] = if sign > 0. { -0.001 } else { 20.001 };
            end[axis] = if sign > 0. { 6. } else { 14. };
            let pilot = orthogonal_pockets::identify_pocket(
                &material, &boundary, approach, end, 4., 1e-7, 1e-4,
            )
            .unwrap();
            let mut axes = [[0.; 3]; 3];
            axes[2][axis] = sign;
            axes[0][(axis + 1) % 3] = 1.;
            axes[1][(axis + 2) % 3] = sign;
            let proof = run(
                &material,
                pilot,
                &insert,
                Frame {
                    origin_mm: pilot.center_mm,
                    axes,
                },
                req,
            )
            .unwrap();
            assert!(proof.tip_clearance_mm[0] > 1.99);
            assert!(proof.witness_volume_lower_mm3 > 0.);
        }
    }
}

#[test]
fn paired_source_error_classes_are_preserved() {
    use crate::wood_insert::Error as E;
    for (cause, expected) in [
        (E::Cancelled, Error::Cancelled),
        (E::Limit, Error::Limit),
        (E::Invalid, Error::Invalid),
        (E::Unsupported, Error::Unsupported),
        (E::Unresolved, Error::Unsupported),
    ] {
        assert_eq!(insert_error(cause), expected);
    }
}
#[test]
fn a_real_overhanging_stock_cell_is_not_waived_by_local_pilot_contact() {
    use crate::orthogonal::{self, Box3, Operation};
    let (_, _, insert, frame, req) = fixture(false, 1);
    let material = orthogonal::boolean(
        &[
            Box3 {
                min: [0.; 3],
                max: [20.; 3],
            },
            Box3 {
                min: [0., 0., 20.],
                max: [1., 20., 25.],
            },
            Box3 {
                min: [0., 0., 24.],
                max: [20., 20., 25.],
            },
        ],
        &[],
        Operation::Union,
        1e-7,
    )
    .unwrap()
    .unwrap();
    assert_eq!(material.solid_count(), 1);
    let boundary = orthogonal_boundary::build(material.quads()).unwrap();
    let pilot = orthogonal_pockets::identify_pocket(
        &material,
        &boundary,
        [10., 10., 20.001],
        [10., 10., 14.],
        4.,
        1e-7,
        1e-4,
    )
    .unwrap();
    assert_eq!(
        run(&material, pilot, &insert, frame, req).unwrap_err(),
        Error::Unsupported
    );
}

#[test]
fn managed_real_inset_preserves_original_source_and_consumes_actual_floor_gap() {
    use std::cell::Cell;
    for left in [false, true] {
        let (material, pilot, insert, original, req) = fixture(left, 1);
        let source = || Stock {
            material: &material,
            pilot,
            frame: identity(),
        };
        let mut placed = original;
        placed.origin_mm[2] -= 1e-9;
        let proof = certify_inset_insert(source(), &insert, placed, req, 1e-7, &|| false).unwrap();
        assert!(proof.inset_mm[0] > 0. && proof.inset_mm[1] < 1e-7);
        assert!(proof.contact.tip_clearance_mm[0] > 1.99 && proof.contact.tip_clearance_mm[1] < 2.);
        assert_eq!(insert.full_length_mm(), 4.);
        assert!(certify_insert(source(), &insert, placed, req, &|| false).is_err());
        for max in [-1., f64::NAN, 1e-3, 1e-10] {
            assert!(certify_inset_insert(source(), &insert, placed, req, max, &|| false).is_err());
        }
        let mut outside = original;
        outside.origin_mm[2] += 1e-9;
        assert!(certify_inset_insert(source(), &insert, outside, req, 1e-7, &|| false).is_err());
        let calls = Cell::new(0);
        certify_inset_insert(source(), &insert, placed, req, 1e-7, &|| {
            calls.set(calls.get() + 1);
            false
        })
        .unwrap();
        for stop in 1..=calls.get() {
            let n = Cell::new(0);
            assert!(matches!(
                certify_inset_insert(source(), &insert, placed, req, 1e-7, &|| {
                    n.set(n.get() + 1);
                    n.get() == stop
                }),
                Err(Error::Cancelled)
            ));
        }
    }
}
