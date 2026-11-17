use super::*;
use crate::{
    Circle, Rational,
    authored_flat_thread::SourceOptions,
    cylinder::AxialCylinder,
    helical_thread::{Spec, Thread},
    thread_end_relief,
    thread_entry_chamfer::Chamfer,
};
fn spec(left: bool, major: f64) -> Spec {
    Spec {
        major_radius_mm: major,
        pitch_mm: 0.5,
        length_mm: 4.,
        depth_mm: 0.25,
        crest_width_mm: 0.125,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: left,
        root_radius_mm: None,
    }
}
fn sources(left: bool, end: Rational) -> (Chamfer, Material) {
    let a = Chamfer::new(
        Thread::new(spec(left, 1.375)).unwrap(),
        0.4,
        Some(thread_end_relief::Intent {
            width_mm: 0.5,
            radial_depth_mm: 0.05,
        }),
    )
    .unwrap();
    let thread = AuthoredThread::new(
        spec(left, 1.5),
        SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-7,
        },
    )
    .unwrap();
    let body = Cavity::new_with_exact_receiving_extent(
        thread,
        Intent {
            entry_depth_mm: Some(0.4),
            terminal: Some(Terminal {
                width_mm: 0.5,
                radial_depth_mm: 0.05,
            }),
        },
        end,
    )
    .unwrap();
    let zero = Rational::new(0, 1).unwrap();
    let stock = AxialCylinder::new(
        7,
        Circle::new([zero; 2], Rational::new(7, 2).unwrap()).unwrap(),
        zero,
        end,
    )
    .unwrap();
    (a, Material::new(stock, body, 1.).unwrap())
}
fn near(radius: f64, z: f64) -> [[f64; 2]; 3] {
    [
        [radius - 0.000001, radius + 0.000001],
        [-0.000001, 0.000001],
        [z - 0.000001, z + 0.000001],
    ]
}
#[test]
fn complete_actual_treatment_boxes_and_nonbinary_finite_end_have_strict_material_classification() {
    use BoxClassification::*;
    for left in [false, true] {
        for end in [Rational::new(4, 1).unwrap(), Rational::new(41, 10).unwrap()] {
            let (a, b) = sources(left, end);
            let original = (a, b);
            // Independent analytic section values: male cone.985 atz=.01,
            // neck1.075; female cone1.64 atz=.01 and terminal radius1.55.
            for (r, z, class) in [
                (0.984, 0.01, Inside),
                (0.986, 0.01, Outside),
                (1.074, 3.8, Inside),
                (1.076, 3.8, Outside),
                (1., -0.01, Outside),
                (1., 4.01, Outside),
            ] {
                assert_eq!(a.classify_treatment_box(near(r, z)).unwrap(), class);
            }
            for (r, z, class) in [
                (1.641, 0.01, Inside),
                (1.639, 0.01, Outside),
                (1.551, 3.8, Inside),
                (1.549, 3.8, Outside),
                (4., 0.01, Outside),
                (1.6, -0.01, Outside),
                (1.6, 4.11, Outside),
            ] {
                assert_eq!(b.classify_treatment_box(near(r, z)).unwrap(), class);
            }
            if end == Rational::new(41, 10).unwrap() {
                assert_eq!(b.classify_treatment_box(near(1.551, 4.05)).unwrap(), Inside);
                assert_eq!(
                    b.classify_treatment_box(near(1.549, 4.05)).unwrap(),
                    Outside
                );
                assert_eq!(
                    b.classify_treatment_box(near(1.551, 4.1)).unwrap(),
                    UnresolvedBoundary
                );
            }
            assert_eq!((a, b), original);
        }
    }
}
#[test]
fn whole_boxes_straddling_actual_material_or_unproved_helical_regions_are_never_samples() {
    use BoxClassification::*;
    let (a, b) = sources(false, Rational::new(4, 1).unwrap());
    // Corner samples would miss actual material inside each of these boxes.
    assert_eq!(
        a.classify_treatment_box([[-2., 2.], [-2., 2.], [0.00999, 0.01001]])
            .unwrap(),
        UnresolvedBoundary
    );
    assert_eq!(
        b.classify_treatment_box([[-4., 4.], [-4., 4.], [0.00999, 0.01001]])
            .unwrap(),
        UnresolvedBoundary
    );
    for bounds in [
        near(0.985, 0.01),
        near(1.075, 3.8),
        near(1., 0.),
        near(1., 4.),
        near(1., 1.),
        near(1.2, 0.2),
    ] {
        assert_eq!(
            a.classify_treatment_box(bounds).unwrap(),
            UnresolvedBoundary
        );
    }
    for bounds in [
        near(1.64, 0.01),
        near(1.55, 3.8),
        near(2., 0.),
        near(2., 4.),
        near(2., 1.),
        near(1.4, 0.3),
    ] {
        assert_eq!(
            b.classify_treatment_box(bounds).unwrap(),
            UnresolvedBoundary
        );
    }
    let no_tail = Chamfer::new(*a.source_thread(), 0.4, None).unwrap();
    assert_eq!(
        no_tail.classify_treatment_box(near(1., 3.8)).unwrap(),
        UnresolvedBoundary
    );
    for (entry, tail, z) in [
        (Some(0.4), None, 3.8),
        (None, b.cavity().intent().terminal, 0.01),
    ] {
        let body = Cavity::new(
            *b.cavity().source_thread(),
            Intent {
                entry_depth_mm: entry,
                terminal: tail,
            },
        )
        .unwrap();
        let material = Material::new(*b.stock(), body, 1.).unwrap();
        assert_eq!(
            material.classify_treatment_box(near(1.6, z)).unwrap(),
            UnresolvedBoundary
        );
    }
    for invalid in [
        [[f64::NAN, 1.], [0., 1.], [0., 1.]],
        [[1., 0.], [0., 1.], [0., 1.]],
    ] {
        assert!(a.classify_treatment_box(invalid).is_err());
        assert!(b.classify_treatment_box(invalid).is_err());
    }
    let huge = [[1e308, 1e308], [0., 1.], [0.01, 0.02]];
    assert_eq!(a.classify_treatment_box(huge).unwrap(), UnresolvedBoundary);
    assert_eq!(b.classify_treatment_box(huge).unwrap(), UnresolvedBoundary);
}
