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
    for (length, gap) in [(3., 3.), (4.5, 1.5), (5.5, 0.5)] {
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
            assert!(query(1.).is_ok());
        }
    }
}
