use super::*;
fn spec() -> CylindricalBlindPocketSpec {
    CylindricalBlindPocketSpec {
        radius_mm: 2.,
        depth_mm: 16.,
        minimum_wall_mm: 1.,
        minimum_floor_mm: 1.,
        seat: Some(multicad_rust_kernel::cylindrical_blind_pocket::Seat {
            radius_mm: 6.,
            depth_mm: 1.,
        }),
    }
}
#[test]
fn rectangular_receiver_both_ends_keep_real_stock_and_closed_manufacturing_geometry() {
    let policy = GeometryPolicy::default();
    let source =
        crate::geometry::construction::box_solid([20.; 3], EntityId::new(), policy).unwrap();
    let original = source.volume_mm3().unwrap();
    for top in [false, true] {
        let frame = PlaneFrame {
            origin: [10., 10., if top { 20. } else { 0. }],
            u: [1., 0., 0.],
            v: [0., if top { -1. } else { 1. }, 0.],
        };
        let cut = source
            .rust_cut_rectangular_blind_counterbore(spec(), frame, EntityId::new(), policy, None)
            .unwrap();
        assert_eq!(source.volume_mm3().unwrap(), original);
        let expected = 8000. - 96. * std::f64::consts::PI;
        assert!((cut.volume_mm3().unwrap() - expected).abs() < 1e-9);
        assert_eq!(cut.bounds().unwrap(), source.bounds().unwrap());
        assert_eq!(cut.faces().unwrap().len(), 22);
        let (bound, actual_frame) = cut.rust_rectangular_blind_counterbore_source().unwrap();
        assert_eq!(actual_frame, frame);
        assert_eq!(bound.spec().seat_radius_mm, 6.);
        assert!(source.rust_rectangular_blind_counterbore_source().is_none());
        let original_topology = source.topology_snapshot(EntityId::new(), 0).unwrap();
        let cut_topology = cut.topology_snapshot(EntityId::new(), 0).unwrap();
        for anchor in original_topology.elements.iter().flat_map(|e| &e.frames) {
            assert!(
                cut_topology
                    .elements
                    .iter()
                    .flat_map(|e| &e.frames)
                    .any(|a| a.origin == anchor.origin && a.frame == anchor.frame)
            );
        }

        let proof = cut
            .manufacturing_mesh_with_receipt(Default::default())
            .unwrap();
        assert!(proof.vertex_links > 0 && proof.edge_pairs > 0);
        let file = tempfile::NamedTempFile::new().unwrap();
        cut.write_step(file.path()).unwrap();
        let imported = Solid::import_rust_step(
            &std::fs::read_to_string(file.path()).unwrap(),
            EntityId::new(),
            policy,
        )
        .unwrap();
        assert!((imported.volume_mm3().unwrap() - expected).abs() < 1e-7);
    }
}
#[test]
fn rectangular_receiver_refuses_unearned_stock_charts_and_preserves_source() {
    let policy = GeometryPolicy::default();
    let id = EntityId::new();
    let source = crate::geometry::construction::box_solid([20.; 3], id, policy).unwrap();
    let frame = PlaneFrame {
        origin: [10., 10., 0.],
        ..Default::default()
    };
    let attempt = |spec, frame, feature, policy, cancel| {
        source.rust_cut_rectangular_blind_counterbore(spec, frame, feature, policy, cancel)
    };
    assert_eq!(
        attempt(spec(), frame, id, policy, None).err().unwrap().kind,
        GeometryErrorKind::Input
    );
    assert_eq!(
        attempt(spec(), frame, EntityId::new(), policy, Some(&|| true))
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Cancelled
    );
    for bad in [
        PlaneFrame {
            origin: [10., 10., 0.5],
            ..frame
        },
        PlaneFrame {
            u: [0.8, 0.6, 0.],
            v: [-0.6, 0.8, 0.],
            ..frame
        },
        PlaneFrame {
            origin: [2f64.powi(-54), 10., 0.],
            ..frame
        },
    ] {
        assert_eq!(
            attempt(spec(), bad, EntityId::new(), policy, None)
                .err()
                .unwrap()
                .kind,
            GeometryErrorKind::Unsupported
        );
    }
    assert_eq!(
        attempt(
            CylindricalBlindPocketSpec {
                seat: None,
                ..spec()
            },
            frame,
            EntityId::new(),
            policy,
            None
        )
        .err()
        .unwrap()
        .kind,
        GeometryErrorKind::Input
    );
    assert!(
        attempt(
            CylindricalBlindPocketSpec {
                minimum_floor_mm: 5.,
                ..spec()
            },
            frame,
            EntityId::new(),
            policy,
            None
        )
        .is_err()
    );
    assert!(
        attempt(
            CylindricalBlindPocketSpec {
                radius_mm: 0.,
                ..spec()
            },
            frame,
            EntityId::new(),
            policy,
            None
        )
        .is_err()
    );
    assert!(attempt(spec(), frame, EntityId::nil(), policy, None).is_err());
    let decimal =
        crate::geometry::construction::box_solid([20., 20., 20.1], EntityId::new(), policy)
            .unwrap();
    assert_eq!(
        decimal
            .rust_cut_rectangular_blind_counterbore(spec(), frame, EntityId::new(), policy, None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(source.volume_mm3().unwrap(), 8000.);
    assert_eq!(source.faces().unwrap().len(), 6);
}
