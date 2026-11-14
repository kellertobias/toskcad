use super::*;
#[test]
fn every_actual_stack_owner_uses_outward_source_bounds_and_unknown_geometry_stays_unknown() {
    let p = GeometryPolicy::default();
    let (wood, hanger, nut, _) = super::super::tests::fixture(false, false, p);
    let ring =
        Solid::rust_annulus(1.7, 3.5, 0.2, EntityId::new(), PlaneFrame::default(), p).unwrap();
    let box_source = Solid::box_solid([2.; 3], p).unwrap();
    let stock = box_source.orthogonal_stock().unwrap();
    let far = box_source
        .transformed(RigidTransform {
            translation_mm: [100.; 3],
            ..Default::default()
        })
        .unwrap();
    for source in [&wood, &hanger, &nut, &ring, &box_source, &stock] {
        assert_eq!(
            source.rust_hanger_component_separation(&far, None).unwrap(),
            Some(true)
        );
        assert_eq!(
            source
                .rust_hanger_component_separation(source, None)
                .unwrap(),
            Some(false)
        );
    }
    let unknown = Solid::cylinder(PlaneFrame::default(), 1., 2., p).unwrap();
    assert_eq!(
        unknown
            .rust_hanger_component_separation(&far, None)
            .unwrap(),
        None
    );
    assert_eq!(
        hanger
            .rust_hanger_component_separation(&unknown, None)
            .unwrap(),
        None
    );
    let different = Solid::box_solid(
        [2.; 3],
        GeometryPolicy {
            linear_tolerance_mm: 1e-6,
            minimum_feature_mm: 1e-4,
            ..p
        },
    )
    .unwrap();
    assert_eq!(
        hanger
            .rust_hanger_component_separation(&different, None)
            .unwrap_err()
            .kind,
        GeometryErrorKind::Input
    );
}
#[test]
fn known_and_unknown_source_queries_poll_every_checkpoint_without_changing_source_identity() {
    use std::cell::Cell;
    let p = GeometryPolicy::default();
    let (_, hanger, _, _) = super::super::tests::fixture(true, true, p);
    let before = hanger.rust_hanger_data();
    let far = Solid::box_solid([2.; 3], p)
        .unwrap()
        .transformed(RigidTransform {
            translation_mm: [100.; 3],
            ..Default::default()
        })
        .unwrap();
    let unknown = Solid::cylinder(PlaneFrame::default(), 1., 2., p).unwrap();
    for other in [&far, &unknown] {
        let n = Cell::new(0);
        hanger
            .rust_hanger_component_separation(
                other,
                Some(&|| {
                    n.set(n.get() + 1);
                    false
                }),
            )
            .unwrap();
        for stop in 0..n.get() {
            let at = Cell::new(0);
            assert_eq!(
                hanger
                    .rust_hanger_component_separation(
                        other,
                        Some(&|| {
                            let n = at.get();
                            at.set(n + 1);
                            n == stop
                        })
                    )
                    .unwrap_err()
                    .kind,
                GeometryErrorKind::Cancelled
            );
        }
    }
    assert_eq!(hanger.rust_hanger_data(), before);
}
