use super::*;
use crate::{
    authored_flat_thread::{AuthoredThread, SourceOptions},
    headed_screw::{ScrewHead, ScrewSpec},
    helical_thread::{Spec, Thread},
};
fn identity() -> Frame {
    Frame {
        origin_mm: [0.; 3],
        axes: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
    }
}
pub(in crate::wood_thread_contact::insert) fn fixture(left: bool) -> (Insert, Screw, Frame) {
    let p = 0.5;
    let major = (2.874 + 2.980) / 4.;
    let pitch_mid = (2.580 + 2.655) / 2.;
    let crest = p / 2. - (2. * major - pitch_mid) / 3_f64.sqrt();
    let root = p / (4. * 3_f64.sqrt());
    let spec = Spec {
        major_radius_mm: major,
        pitch_mm: p,
        length_mm: 8.,
        depth_mm: (p - crest) / (2. / 3_f64.sqrt()) - root,
        crest_width_mm: crest,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: left,
        root_radius_mm: Some(root),
    };
    let male = Thread::new(spec).unwrap();
    let make = |s| {
        AuthoredThread::new(
            s,
            SourceOptions {
                position_limit_mm: 1e-8,
                normal_angle_limit_rad: 1e-6,
            },
        )
        .unwrap()
    };
    let inner = make(Spec {
        major_radius_mm: 1.5,
        length_mm: 4.,
        depth_mm: (3. - (2.459 + 2.599) / 2.) / 2.,
        crest_width_mm: p / 2. - (3. - (2.675 + 2.775) / 2.) / 3_f64.sqrt(),
        root_radius_mm: None,
        ..spec
    });
    let wood = make(Spec {
        major_radius_mm: 5.,
        pitch_mm: 2.,
        length_mm: 4.,
        depth_mm: 1.,
        crest_width_mm: 0.4,
        root_radius_mm: None,
        ..spec
    });
    let insert = Insert::new(wood, inner, 0.5, None).unwrap();
    let screw = Screw::new(
        male,
        ScrewSpec {
            head: ScrewHead::SocketCap {
                radius_mm: 2.75,
                height_mm: 3.,
                socket_across_flats_mm: 2.,
                socket_depth_mm: 1.,
            },
            shank_length_mm: 0.,
            shank_radius_mm: 1.5,
        },
    )
    .unwrap();
    let angle = if left {
        72_f64.to_radians()
    } else {
        -72_f64.to_radians()
    };
    let (s, c) = angle.sin_cos();
    (
        insert,
        screw,
        Frame {
            origin_mm: [0., 0., -3.1],
            axes: [[c, s, 0.], [-s, c, 0.], [0., 0., 1.]],
        },
    )
}
fn options() -> ClockedMatingOptions {
    ClockedMatingOptions {
        requirements: MatingRequirements {
            minimum_engagement_mm: 4.,
            minimum_clearance_mm: 4e-7,
        },
        maximum_cells: 4096,
        linear_tolerance_mm: 1e-7,
    }
}
#[test]
fn actual_original_metric_laws_qualify_general_clock_with_complete_finite_overlap() {
    for left in [false, true] {
        let (insert, screw, frame) = fixture(left);
        assert!(
            certify_insert_mating(
                &insert,
                identity(),
                MatingExternal {
                    source: MatingScrew::Regular(&screw),
                    frame
                },
                options().requirements,
                &|| false
            )
            .is_err()
        );
        let v =
            certify_clocked_insert_mating(&insert, identity(), &screw, frame, options(), &|| false)
                .unwrap();
        assert!(v.mating.engagement_mm[0] <= 4. && v.mating.engagement_mm[1] >= 4.);
        assert!(v.mating.minimum_radial_clearance_lower_mm > 0.007);
        assert!(
            v.mating.source_phase_offset_mm[0] <= -3. && v.mating.source_phase_offset_mm[1] >= -3.
        );
        assert!(
            v.radial_map_error_mm < 1e-10
                && v.axial_map_error_mm < 1e-10
                && v.angle_error_rad < 1e-10
        );
        assert!(!v.mating.entire_thread_within_insert);
        assert_eq!(screw.thread().spec().length_mm, 8.);
        assert_eq!(insert.full_length_mm(), 4.);
        assert!(v.cells_used <= 4096);
    }
}
#[test]
fn actual_source_domain_frame_phase_and_resource_failures_remain_typed() {
    let (insert, screw, frame) = fixture(false);
    for n in 0..8 {
        let mut r = options();
        let mut f = frame;
        match n {
            0 => r.maximum_cells = 0,
            1 => r.maximum_cells = 1,
            2 => r.requirements.minimum_engagement_mm = 4_f64.next_up(),
            3 => f.origin_mm[2] = 0.01,
            4 => f.origin_mm[0] = 0.01,
            5 => f.axes[2] = [0., 0., -1.],
            6 => r.linear_tolerance_mm = f64::NAN,
            _ => r.requirements.minimum_clearance_mm = 1.,
        }
        assert!(
            certify_clocked_insert_mating(&insert, identity(), &screw, f, r, &|| false).is_err(),
            "case{n}"
        );
    }
    let wrong = Thread::new(Spec {
        pitch_mm: 0.6,
        root_radius_mm: None,
        ..screw.thread().spec()
    })
    .unwrap();
    let wrong = Screw::new(wrong, screw.spec()).unwrap();
    assert!(matches!(
        certify_clocked_insert_mating(&insert, identity(), &wrong, frame, options(), &|| false),
        Err(Error::Unsupported)
    ));
    let short = Thread::new(Spec {
        length_mm: 4.,
        ..screw.thread().spec()
    })
    .unwrap();
    let short = Screw::new(short, screw.spec()).unwrap();
    assert!(matches!(
        certify_clocked_insert_mating(&insert, identity(), &short, frame, options(), &|| false),
        Err(Error::Unsupported)
    ));
}
#[test]
fn cancellation_at_every_actual_query_checkpoint_is_atomic_and_retry_is_stable() {
    use std::cell::Cell;
    let (insert, screw, frame) = fixture(false);
    let calls = Cell::new(0);
    let first =
        certify_clocked_insert_mating(&insert, identity(), &screw, frame, options(), &|| {
            calls.set(calls.get() + 1);
            false
        })
        .unwrap();
    for stop in 0..calls.get() {
        let count = Cell::new(0);
        assert!(matches!(
            certify_clocked_insert_mating(&insert, identity(), &screw, frame, options(), &|| {
                let n = count.get();
                count.set(n + 1);
                n == stop
            }),
            Err(Error::Cancelled)
        ));
    }
    let again =
        certify_clocked_insert_mating(&insert, identity(), &screw, frame, options(), &|| false)
            .unwrap();
    assert_eq!(
        again.mating.source_phase_offset_mm,
        first.mating.source_phase_offset_mm
    );
}

#[test]
fn actual_whole_screw_wood_clearance_is_independent_of_mating_and_insert_length() {
    use std::cell::Cell;
    let (material, mut pilot, _, _, _) = crate::wood_thread_contact::tests::fixture(false, 1);
    pilot.radius_mm = 4.;
    pilot.depth_mm = 6.;
    let (_, screw, f) = fixture(false);
    let place = Frame {
        origin_mm: [10., 10., 23.1],
        axes: f.axes.map(|a| [a[0], -a[1], -a[2]]),
    };
    let req = ScrewClearanceRequirements {
        minimum_tip_clearance_mm: 1.,
        linear_tolerance_mm: 1e-7,
        minimum_feature_mm: 1e-4,
    };
    let run = |frame, r, interrupt: &dyn Fn() -> bool| {
        certify_clocked_insert_screw_clearance(
            Stock {
                material: &material,
                pilot,
                frame: identity(),
            },
            &screw,
            frame,
            r,
            interrupt,
        )
    };
    let proof = run(place, req, &|| false).unwrap();
    assert!(proof.tip_clearance_mm[0] > 1.09 && proof.tip_clearance_mm[0] < 1.11);
    assert!(proof.minimum_radial_clearance_mm > 1.3);
    let mut bad = req;
    bad.minimum_tip_clearance_mm = 1.11;
    assert!(matches!(
        run(place, bad, &|| false),
        Err(Error::Unsupported)
    ));
    bad.minimum_tip_clearance_mm = f64::NAN;
    assert!(matches!(run(place, bad, &|| false), Err(Error::Invalid)));
    let mut intruding = place;
    intruding.origin_mm[2] = 19.9;
    assert!(matches!(
        run(intruding, req, &|| false),
        Err(Error::Unsupported)
    ));
    let calls = Cell::new(0);
    assert!(
        run(place, req, &|| {
            calls.set(calls.get() + 1);
            false
        })
        .is_ok()
    );
    for stop in 0..calls.get() {
        let count = Cell::new(0);
        assert!(matches!(
            run(place, req, &|| {
                let n = count.get();
                count.set(n + 1);
                n == stop
            }),
            Err(Error::Cancelled)
        ));
    }
}

#[test]
fn original_large_core_and_small_polar_domain_cannot_borrow_a_clearance_envelope() {
    let (insert, original, frame) = fixture(false);
    let big = Thread::new(Spec {
        major_radius_mm: 10.,
        depth_mm: 0.2,
        root_radius_mm: None,
        ..original.thread().spec()
    })
    .unwrap();
    let big = Screw::new(
        big,
        ScrewSpec {
            head: ScrewHead::ExternalHex {
                across_flats_mm: 30.,
                height_mm: 3.,
            },
            shank_length_mm: 0.,
            shank_radius_mm: 10.,
        },
    )
    .unwrap();
    assert!(matches!(
        certify_clocked_insert_mating(&insert, identity(), &big, frame, options(), &|| false),
        Err(Error::Unsupported)
    ));
    let tiny = Thread::new(Spec {
        major_radius_mm: 0.01,
        depth_mm: 0.009999999,
        root_radius_mm: None,
        ..original.thread().spec()
    })
    .unwrap();
    let tiny = Screw::new(
        tiny,
        ScrewSpec {
            head: ScrewHead::ExternalHex {
                across_flats_mm: 0.04,
                height_mm: 0.01,
            },
            shank_length_mm: 0.,
            shank_radius_mm: 0.01,
        },
    )
    .unwrap();
    let mut moved = frame;
    moved.origin_mm[0] = 1e-8;
    assert!(matches!(
        certify_clocked_insert_mating(&insert, identity(), &tiny, moved, options(), &|| false),
        Err(Error::Unsupported)
    ));
}
#[test]
fn actual_numeric_ray_certificates_refuse_uncertainty_without_inventing_a_body() {
    assert!(matches!(
        ray_angle(I::new(-1., 1.), ip(1.), I::new(1., 2.)),
        Err(Error::Unsupported)
    ));
    assert!(matches!(
        ray_angle(ip(1.), ip(1.), I::new(1e-320, 2.)),
        Err(Error::Unsupported)
    ));
    assert!(exact(&Q::from_integer(num_bigint::BigInt::from(1) << 2000)).is_err());
}

#[test]
fn actual_multibody_and_large_retained_stock_are_not_accepted_by_a_small_query() {
    use crate::orthogonal::{self, Box3, Operation};
    let (_, screw, frame) = fixture(false);
    let req = ScrewClearanceRequirements {
        minimum_tip_clearance_mm: 1.,
        linear_tolerance_mm: 1e-7,
        minimum_feature_mm: 1e-4,
    };
    let query = |material: &OrthogonalSolid, center: [f64; 3], end: [f64; 3], r| {
        let boundary = orthogonal_boundary::build(material.quads()).unwrap();
        let pilot =
            orthogonal_pockets::identify_pocket(material, &boundary, center, end, r, 1e-7, 1e-4)
                .unwrap();
        certify_clocked_insert_screw_clearance(
            Stock {
                material,
                pilot,
                frame: identity(),
            },
            &screw,
            frame,
            req,
            &|| false,
        )
    };
    let separate = orthogonal::boolean(
        &[
            Box3 {
                min: [0.; 3],
                max: [20.; 3],
            },
            Box3 {
                min: [30.; 3],
                max: [50.; 3],
            },
        ],
        &[],
        Operation::Union,
        1e-7,
    )
    .unwrap()
    .unwrap();
    assert!(matches!(
        query(&separate, [10., 10., 20.001], [10., 10., 14.], 4.),
        Err(Error::Unsupported)
    ));
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
    assert!(matches!(
        query(&large, [0.5, 0.5, 1.001], [0.5, 0.5, 0.5], 0.1),
        Err(Error::Limit)
    ));
}

#[test]
fn actual_periodic_phillips_material_qualifies_shared_general_clock_receiver() {
    for left in [false, true] {
        for plain in [0., 0.01, 1.] {
            let (insert, owner) = super::super::mating::tests::periodic_fixture(left, 8., plain);
            let (_, _, frame) = fixture(left);
            let receiver = Receiver::insert(&insert, identity()).unwrap();
            let proof = certify_clocked_receiver(
                &receiver,
                MatingExternal {
                    source: MatingScrew::Periodic(&owner),
                    frame,
                },
                options(),
                &|| false,
            )
            .unwrap();
            assert!(proof.bounds.minimum_radial_clearance_lower_mm > 0.005);
            assert_eq!(proof.bounds.external_indexed_phase_error_mm, 0.);
            assert_eq!(proof.bounds.external_evaluator_phase_error_mm, 0.);
            assert!(proof.cells_used > 0 && proof.cells_used <= 4096);
        }
    }
}
