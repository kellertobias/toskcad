use super::*;
use crate::rust_kernel::thread_end_relief;
fn external(major: f64, left: bool, policy: GeometryPolicy) -> Solid {
    let mut request = super::super::tests::cut(PlaneFrame::default(), left, 1);
    request.profile.major_radius_mm = major;
    Solid::rust_helical_thread_entry_chamfer(
        request.profile,
        None,
        0.2,
        Some(thread_end_relief::Intent {
            width_mm: 0.5,
            radial_depth_mm: 0.05,
        }),
        EntityId::new(),
        PlaneFrame::default(),
        policy,
    )
    .unwrap()
}
#[test]
fn actual_roles_laws_policy_and_added_geometry_cannot_borrow_pair_authority() {
    let policy = GeometryPolicy::default();
    let frame = PlaneFrame::default();
    let stock = super::super::tests::fixture_stock(Rational::new(4, 1).unwrap(), frame);
    let request = super::super::tests::cut(frame, false, 1);
    let receiving = stock
        .rust_cut_authored_internal_thread_treatments_with_interrupt(
            request,
            EntityId::new(),
            policy,
            None,
        )
        .unwrap();
    let male = external(1.375, false, policy);
    let positive = male
        .rust_chamfer_receiver_separation(&receiving, Default::default(), None)
        .unwrap();
    assert!(positive.minimum_radial_gap_mm() > 0.);
    let original = (
        male.rust_thread_entry_chamfer_data().unwrap(),
        receiving.rust_internal_thread_treatment_data().unwrap(),
    );
    let overlap = external(1.625, false, policy);
    assert_eq!(
        overlap
            .rust_chamfer_receiver_separation(&receiving, Default::default(), None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert!(
        overlap
            .boolean(&receiving, BooleanOperation::Intersect, policy)
            .is_err()
    );
    let wrong_hand = external(1.375, true, policy);
    assert_eq!(
        wrong_hand
            .rust_chamfer_receiver_separation(&receiving, Default::default(), None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Input
    );
    let wrong_policy = external(
        1.375,
        false,
        GeometryPolicy {
            linear_tolerance_mm: policy.linear_tolerance_mm * 2.,
            ..policy
        },
    );
    assert_eq!(
        wrong_policy
            .rust_chamfer_receiver_separation(&receiving, Default::default(), None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Unsupported
    );
    let cutter = stock
        .rust_cut_external_thread_entry_chamfer_with_interrupt(
            request.profile,
            None,
            0.2,
            Some(thread_end_relief::Intent {
                width_mm: 0.5,
                radial_depth_mm: 0.05,
            }),
            EntityId::new(),
            frame,
            policy,
            None,
        )
        .unwrap();
    assert!(cutter.rust_thread_entry_chamfer_data().is_some());
    assert_eq!(
        cutter
            .rust_chamfer_receiver_separation(&receiving, Default::default(), None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert!(
        cutter
            .boolean(&receiving, BooleanOperation::Intersect, policy)
            .is_err()
    );
    assert_eq!(
        stock
            .rust_chamfer_receiver_separation(&receiving, Default::default(), None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(
        male.rust_chamfer_receiver_separation(&stock, Default::default(), None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Unsupported
    );
    let extra = Solid::box_solid([1.; 3], policy).unwrap();
    let compound = Solid::compound_named(&[&male, &extra], EntityId::new(), policy).unwrap();
    assert_eq!(
        compound
            .rust_chamfer_receiver_separation(&receiving, Default::default(), None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert!(
        compound
            .boolean(&receiving, BooleanOperation::Intersect, policy)
            .is_err()
    );
    assert!(
        male.boolean(
            &receiving,
            BooleanOperation::Intersect,
            GeometryPolicy {
                boolean_fuzzy_mm: 0.001,
                ..policy
            }
        )
        .is_err()
    );
    assert_eq!(
        (
            male.rust_thread_entry_chamfer_data().unwrap(),
            receiving.rust_internal_thread_treatment_data().unwrap()
        ),
        original
    );
    assert_eq!(
        male.rust_chamfer_receiver_separation(&receiving, Default::default(), None)
            .unwrap()
            .minimum_radial_gap_mm(),
        positive.minimum_radial_gap_mm()
    );
}

#[test]
fn distinct_exact_receiving_end_cannot_borrow_a_stock_with_the_same_rounded_height() {
    let correct = Rational::new(4, 1).unwrap();
    let wrong = Rational::new(40_000_000_000_000_001, 10_000_000_000_000_000).unwrap();
    assert_eq!(correct.as_f64(), wrong.as_f64());
    assert_ne!(correct, wrong);
    for left in [false, true] {
        let frame = PlaneFrame::default();
        let stock = super::super::tests::fixture_stock(correct, frame);
        let request = super::super::tests::cut(frame, left, 1);
        let before = (stock.bounds().unwrap(), stock.volume_mm3().unwrap());
        assert_eq!(
            stock
                .rust_cut_authored_internal_thread_treatments_with_exact_extent_and_interrupt(
                    request,
                    wrong,
                    EntityId::new(),
                    GeometryPolicy::default(),
                    None,
                )
                .err()
                .unwrap()
                .kind,
            GeometryErrorKind::Unsupported
        );
        assert_eq!(
            (stock.bounds().unwrap(), stock.volume_mm3().unwrap()),
            before
        );
        let retry = stock
            .rust_cut_authored_internal_thread_treatments_with_exact_extent_and_interrupt(
                request,
                correct,
                EntityId::new(),
                GeometryPolicy::default(),
                None,
            )
            .unwrap();
        assert_eq!(
            retry
                .rust_internal_thread_treatment_data()
                .unwrap()
                .material
                .stock()
                .height(),
            correct
        );
    }
}
