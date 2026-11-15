use super::*;
use crate::{headed_screw::Screw, helical_thread::Thread};
fn identity() -> Frame {
    Frame {
        origin_mm: [0.; 3],
        axes: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
    }
}
#[test]
fn actual_screw_continuation_has_its_own_floor_gap_independent_of_insert_length() {
    let (material, mut pilot, original, frame, _) = super::super::super::tests::fixture(false, 1);
    pilot.radius_mm = 4.;
    pilot.depth_mm = 6.;
    for (length, gap) in [(3., 3.), (4.5, 1.5), (5., 1.), (5.5, 0.5)] {
        let mut spec = original.thread().spec();
        spec.length_mm = length;
        let screw = Screw::new(Thread::new(spec).unwrap(), original.spec()).unwrap();
        let query = |required| {
            certify_insert_screw_clearance(
                Stock {
                    material: &material,
                    pilot,
                    frame: identity(),
                },
                MatingExternal {
                    source: mating::MatingScrew::Regular(&screw),
                    frame,
                },
                ScrewClearanceRequirements {
                    minimum_tip_clearance_mm: required,
                    linear_tolerance_mm: 1e-7,
                    minimum_feature_mm: 1e-4,
                },
                &|| false,
            )
        };
        let proof = query(0.4).unwrap();
        assert!(proof.tip_clearance_mm[0] <= gap && proof.tip_clearance_mm[1] >= gap);
        assert!(proof.minimum_radial_clearance_mm > 0.99);
        if length == 5.5 {
            assert_eq!(query(1.).unwrap_err(), Error::Unsupported);
        } else {
            assert!(query(1.).unwrap().tip_clearance_mm[0] >= 1.);
        }
    }
}
#[test]
fn whole_screw_stock_proof_refuses_foreign_axes_oversize_radius_bad_prefix_and_resources() {
    use crate::orthogonal::{self, Box3, Operation};
    use std::cell::Cell;
    let (material, mut pilot, original, frame, _) = super::super::super::tests::fixture(false, 1);
    pilot.radius_mm = 4.;
    pilot.depth_mm = 6.;
    let mut s = original.thread().spec();
    s.length_mm = 4.5;
    let screw = Screw::new(Thread::new(s).unwrap(), original.spec()).unwrap();
    let req = ScrewClearanceRequirements {
        minimum_tip_clearance_mm: 1.,
        linear_tolerance_mm: 1e-7,
        minimum_feature_mm: 1e-4,
    };
    let run = |mat: &OrthogonalSolid, pilot, frame, req| {
        certify_insert_screw_clearance(
            Stock {
                material: mat,
                pilot,
                frame: identity(),
            },
            MatingExternal {
                source: mating::MatingScrew::Regular(&screw),
                frame,
            },
            req,
            &|| false,
        )
    };
    for case in 0..5 {
        let (mut p, mut f) = (pilot, frame);
        match case {
            0 => p.radius_mm = 2.5,
            1 => f.origin_mm[2] -= 0.01,
            2 => {
                f.axes = identity().axes;
            }
            3 => p.face = usize::MAX,
            _ => f.origin_mm[0] += 2.,
        }
        assert!(run(&material, p, f, req).is_err(), "case{case}");
    }
    for value in [f64::NAN, -1., 1e8] {
        let mut bad = req;
        bad.minimum_tip_clearance_mm = value;
        assert_eq!(
            run(&material, pilot, frame, bad).unwrap_err(),
            Error::Invalid
        );
    }
    let mut bad = req;
    bad.minimum_feature_mm = 0.;
    assert_eq!(
        run(&material, pilot, frame, bad).unwrap_err(),
        Error::Invalid
    );
    let mut bad = req;
    bad.linear_tolerance_mm = 1e-8;
    assert_eq!(
        run(&material, pilot, frame, bad).unwrap_err(),
        Error::Invalid
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
        run(&split, pilot, frame, req).unwrap_err(),
        Error::Unsupported
    );
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
    assert_eq!(run(&large, pilot, frame, req).unwrap_err(), Error::Limit);
    let calls = Cell::new(0);
    certify_insert_screw_clearance(
        Stock {
            material: &material,
            pilot,
            frame: identity(),
        },
        MatingExternal {
            source: mating::MatingScrew::Regular(&screw),
            frame,
        },
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
            certify_insert_screw_clearance(
                Stock {
                    material: &material,
                    pilot,
                    frame: identity()
                },
                MatingExternal {
                    source: mating::MatingScrew::Regular(&screw),
                    frame
                },
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
}

#[test]
fn actual_periodic_thread_continuation_keeps_independent_blind_floor_proof() {
    let (material, mut pilot, _, frame, _) = super::super::super::tests::fixture(false, 1);
    pilot.radius_mm = 4.;
    pilot.depth_mm = 6.;
    for (length, gap) in [(3., 3.), (4.5, 1.5), (5.5, 0.5)] {
        let (_, owner) = super::super::mating::tests::periodic_fixture(false, length, 1.);
        let requirements = ScrewClearanceRequirements {
            minimum_tip_clearance_mm: 0.4,
            linear_tolerance_mm: 1e-7,
            minimum_feature_mm: 1e-4,
        };
        let query = |required| {
            certify_insert_screw_clearance(
                Stock {
                    material: &material,
                    pilot,
                    frame: identity(),
                },
                MatingExternal {
                    source: mating::MatingScrew::Periodic(&owner),
                    frame,
                },
                ScrewClearanceRequirements {
                    minimum_tip_clearance_mm: required,
                    ..requirements
                },
                &|| false,
            )
        };
        let proof = query(0.4).unwrap();
        assert!(proof.tip_clearance_mm[0] <= gap && proof.tip_clearance_mm[1] >= gap);
        assert!(proof.minimum_radial_clearance_mm > 2.5);
        if length == 5.5 {
            assert_eq!(query(1.).unwrap_err(), Error::Unsupported)
        } else {
            assert!(query(1.).is_ok())
        }
    }
}
