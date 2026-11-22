//! Box probes against drilled orthogonal stock: void inside a bore, the
//! plain cells beside it, and no empty claim when a probe reaches the wall.
use multicad_kernel::{BooleanOperation, GeometryPolicy, PlaneFrame, RigidTransform, Solid};

fn policy() -> GeometryPolicy {
    GeometryPolicy::default()
}

fn block(min: [f64; 3], size: [f64; 3]) -> Solid {
    Solid::box_solid(size, policy())
        .unwrap()
        .transformed(RigidTransform {
            translation_mm: min,
            ..Default::default()
        })
        .unwrap()
}

/// A 20 × 20 × 5 plate less its corner [15, 20]², with a bore of radius 2
/// through it at (5, 5).
fn drilled_l_plate() -> Solid {
    let plate = block([0.; 3], [20., 20., 5.])
        .boolean(
            &block([15., 15., -1.], [6., 6., 7.]),
            BooleanOperation::Subtract,
            policy(),
        )
        .unwrap()
        .unwrap();
    let bore = Solid::cylinder(
        PlaneFrame {
            origin: [5., 5., -1.],
            ..Default::default()
        },
        2.,
        7.,
        policy(),
    )
    .unwrap();
    plate
        .boolean(&bore, BooleanOperation::Subtract, policy())
        .unwrap()
        .unwrap()
}

/// A 2 × 2 column from z = -20 to 10 centred on `(x, y)`.
fn probe(x: f64, y: f64) -> Solid {
    block([x - 1., y - 1., -20.], [2., 2., 30.])
}

#[test]
fn a_column_inside_the_bore_is_empty_and_one_beside_it_meets_the_plate() {
    let plate = drilled_l_plate();
    for (a, b) in [(&plate, &probe(5., 5.)), (&probe(5., 5.), &plate)] {
        assert!(
            a.boolean(b, BooleanOperation::Intersect, policy())
                .unwrap()
                .is_none()
        );
    }
    let beside = probe(12., 5.);
    for (a, b) in [(&plate, &beside), (&beside, &plate)] {
        let common = a
            .boolean(b, BooleanOperation::Intersect, policy())
            .unwrap()
            .unwrap();
        assert_eq!(common.volume_mm3().unwrap(), 20.);
        let bounds = common.bounds().unwrap();
        assert_eq!((bounds.min[2], bounds.max[2]), (0., 5.));
    }
    // Reaching the bore wall: never reported empty.
    assert!(!matches!(
        plate.boolean(&probe(6.5, 5.), BooleanOperation::Intersect, policy()),
        Ok(None)
    ));
}
