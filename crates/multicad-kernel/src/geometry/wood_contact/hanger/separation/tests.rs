use super::*;
fn unknown(p: GeometryPolicy) -> Solid {
    let q = |n| multicad_rust_kernel::Rational::new(n, 1).unwrap();
    Solid::from_rust_prism(
        multicad_rust_kernel::prism::Prism::extrude(19, q(10), q(2), q(2)).unwrap(),
        p,
    )
    .unwrap()
}
fn unhandled_thread(p: GeometryPolicy) -> Solid {
    Solid::rust_helical_thread(
        crate::HelicalThread {
            major_radius_mm: 1.5,
            pitch_mm: 0.5,
            length_mm: 2.,
            depth_mm: 0.2,
            crest_width_mm: 0.1,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            radial_taper: 0.,
            starts: 1,
            left_hand: false,
        },
        None,
        EntityId::new(),
        PlaneFrame::default(),
        p,
    )
    .unwrap()
}
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
    let cylinder = Solid::cylinder(
        PlaneFrame {
            origin: [31., -17., 9.],
            u: [0.8, 0.6, 0.],
            v: [-0.36, 0.48, 0.8],
        },
        1.,
        2.,
        p,
    )
    .unwrap();
    let sphere = Solid::rust_sphere_named(2.3, PlaneFrame::default(), EntityId::new(), p).unwrap();
    for source in [
        &wood,
        &hanger,
        &nut,
        &ring,
        &box_source,
        &stock,
        &cylinder,
        &sphere,
    ] {
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
    for unknown in [unknown(p), unhandled_thread(p)] {
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
    }
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
    let unknown = unknown(p);
    let unhandled_thread = unhandled_thread(p);
    for other in [&far, &unknown, &unhandled_thread] {
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
