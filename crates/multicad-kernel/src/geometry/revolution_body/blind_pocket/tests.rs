use super::*;
use multicad_core::references::{GeometryIntent, Resolution, TopologyKind, resolve};
fn spec() -> CylindricalBlindPocketSpec {
    CylindricalBlindPocketSpec {
        radius_mm: 2.,
        depth_mm: 4.,
        seat: Some(multicad_rust_kernel::cylindrical_blind_pocket::Seat {
            radius_mm: 4.,
            depth_mm: 2.,
        }),
        minimum_wall_mm: 1.,
        minimum_floor_mm: 1.,
    }
}
fn stock(frame: PlaneFrame, feature: EntityId, radius: f64) -> Solid {
    crate::geometry::construction::cylinder(frame, radius, 10., feature, GeometryPolicy::default())
        .unwrap()
}
fn entry(frame: PlaneFrame, end: bool) -> PlaneFrame {
    PlaneFrame {
        origin: position(frame, [0., 0., if end { 10. } else { 0. }]),
        u: frame.u,
        v: if end { frame.v.map(|v| -v) } else { frame.v },
    }
}
#[test]
fn actual_both_cap_pockets_manufacture_and_keep_typed_stock_provenance() {
    for end in [false, true] {
        let frame = PlaneFrame {
            origin: [10., 20., 30.],
            ..Default::default()
        };
        let original = stock(frame, EntityId::new(), 10.);
        let cut = original
            .rust_cut_cylindrical_blind_pocket(
                spec(),
                entry(frame, end),
                EntityId::new(),
                GeometryPolicy::default(),
                None,
            )
            .unwrap();
        let (source, pose) = cut.rust_cylindrical_blind_pocket_geometry().unwrap();
        assert_eq!(pose, frame);
        assert_eq!(
            source.stock().height(),
            original.axial_cylinder().unwrap().height()
        );
        assert_eq!(source.spec(), spec());
        assert_eq!(source.end(), end);
        assert_eq!(cut.bounds().unwrap(), original.bounds().unwrap());
        let topology = cut.topology_snapshot(EntityId::new(), 0).unwrap();
        let shoulder = topology
            .elements
            .iter()
            .find(|e| {
                e.key.kind == TopologyKind::Face
                    && e.origins
                        .iter()
                        .any(|o| o.role == "connection.receiver/seat.shoulder")
            })
            .unwrap();
        let face = cut
            .faces()
            .unwrap()
            .into_iter()
            .find(|face| face.index == shoulder.key.ordinal)
            .unwrap();
        assert!((face.area_mm2 - 12. * std::f64::consts::PI).abs() < 1e-10);
        let (point, normal) = face.plane.unwrap();
        assert!((point[2] - (if end { 38. } else { 32. })).abs() < 1e-10);
        assert!(normal[2] * (if end { 1. } else { -1. }) > 0.99);
        assert!((cut.volume_mm3().unwrap() - 960. * std::f64::consts::PI).abs() < 1e-9);
        let proof = cut
            .manufacturing_mesh_with_receipt(Default::default())
            .unwrap();
        assert!(proof.edge_pairs > 0 && proof.vertex_links > 0);
        assert_eq!(proof.removed_exact_zero_triangles, 0);
        let copy = cut.share().unwrap();
        assert!(std::ptr::eq(
            source,
            copy.rust_cylindrical_blind_pocket_geometry().unwrap().0
        ));
        assert!(cut.supports_step_export());
        let file = tempfile::NamedTempFile::new().unwrap();
        cut.write_step(file.path()).unwrap();
        let imported = Solid::import_rust_step(
            &std::fs::read_to_string(file.path()).unwrap(),
            EntityId::new(),
            GeometryPolicy::default(),
        )
        .unwrap();
        assert!((imported.volume_mm3().unwrap() - cut.volume_mm3().unwrap()).abs() < 1e-9);
        assert!(original.rust_cylindrical_blind_pocket_geometry().is_none());
    }
}
#[test]
fn proper_source_axes_are_retained_and_occurrences_do_not_mutate_metadata() {
    let frame = PlaneFrame {
        origin: [17., -8., 13.],
        u: [0.8, 0.6, 0.],
        v: [-0.36, 0.48, 0.8],
    };
    for end in [false, true] {
        let source = stock(frame, EntityId::new(), 10.);
        let cut = source
            .rust_cut_cylindrical_blind_pocket(
                spec(),
                entry(frame, end),
                EntityId::new(),
                GeometryPolicy::default(),
                None,
            )
            .unwrap();
        let (receipt, pose) = cut.rust_cylindrical_blind_pocket_geometry().unwrap();
        assert_eq!(pose, frame);
        assert_eq!(receipt.frame().axes, [frame.u, frame.v, frame.normal()]);
        cut.manufacturing_mesh_with_receipt(Default::default())
            .unwrap();
        let moved = cut
            .transformed(crate::RigidTransform {
                rotation: [[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]],
                translation_mm: [11., 12., 13.],
            })
            .unwrap();
        let (retained, moved_frame) = moved.rust_cylindrical_blind_pocket_geometry().unwrap();
        assert!(std::ptr::eq(receipt, retained));
        assert_ne!(moved_frame, pose);
    }
}
#[test]
fn source_cap_outer_and_rim_ancestry_survives_radius_and_depth_edits() {
    let stock_feature = EntityId::new();
    let feature = EntityId::new();
    let owner = EntityId::new();
    let frame = PlaneFrame::default();
    let original = stock(frame, stock_feature, 10.);
    let topology = original.topology_snapshot(owner, 0).unwrap();
    let captures = topology
        .elements
        .iter()
        .filter(|e| e.key.kind == TopologyKind::Face || e.key.kind == TopologyKind::Edge)
        .map(|e| GeometryIntent::capture(&topology, e.key).unwrap())
        .collect::<Vec<_>>();
    for (radius, depth) in [(10., 4.), (11., 5.)] {
        let cut = stock(frame, stock_feature, radius)
            .rust_cut_cylindrical_blind_pocket(
                CylindricalBlindPocketSpec {
                    depth_mm: depth,
                    ..spec()
                },
                entry(frame, true),
                feature,
                GeometryPolicy::default(),
                None,
            )
            .unwrap();
        let current = cut.topology_snapshot(owner, 1).unwrap();
        for intent in &captures {
            assert!(
                matches!(resolve(intent, &current, 1), Resolution::Resolved { .. }),
                "{:?}",
                resolve(intent, &current, 1)
            );
        }
        assert_eq!(
            cut.rust_cylindrical_blind_pocket_geometry()
                .unwrap()
                .0
                .spec()
                .depth_mm,
            depth
        );
    }
}
#[test]
fn invalid_prior_cut_and_cancelled_sources_cannot_publish_or_reuse_pocket_metadata() {
    let source = stock(PlaneFrame::default(), EntityId::new(), 10.);
    let feature = EntityId::new();
    let make = |spec, frame, feature, policy, cancel| {
        source.rust_cut_cylindrical_blind_pocket(spec, frame, feature, policy, cancel)
    };
    assert_eq!(
        make(
            spec(),
            entry(PlaneFrame::default(), false),
            EntityId::nil(),
            GeometryPolicy::default(),
            None
        )
        .err()
        .unwrap()
        .kind,
        GeometryErrorKind::Input
    );
    let mut offaxis = entry(PlaneFrame::default(), false);
    offaxis.origin[0] = 0.1;
    assert_eq!(
        make(spec(), offaxis, feature, GeometryPolicy::default(), None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(
        make(
            spec(),
            entry(PlaneFrame::default(), false),
            feature,
            GeometryPolicy::default(),
            Some(&|| true)
        )
        .err()
        .unwrap()
        .kind,
        GeometryErrorKind::Cancelled
    );
    let cut = make(
        spec(),
        entry(PlaneFrame::default(), false),
        feature,
        GeometryPolicy::default(),
        None,
    )
    .unwrap();
    assert!(!cut.rust_cylindrical_stock_available());
    assert_eq!(
        cut.rust_cut_cylindrical_blind_pocket(
            spec(),
            entry(PlaneFrame::default(), false),
            EntityId::new(),
            GeometryPolicy::default(),
            None
        )
        .err()
        .unwrap()
        .kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(
        cut.rust_cylindrical_blind_pocket_geometry()
            .unwrap()
            .0
            .spec(),
        spec()
    );
}

#[test]
fn policy_source_identity_and_delayed_cancellation_refusals_keep_the_original_stock() {
    let identity = EntityId::new();
    let original = stock(PlaneFrame::default(), identity, 10.);
    let frame = entry(PlaneFrame::default(), false);
    assert_eq!(
        original
            .rust_cut_cylindrical_blind_pocket(
                spec(),
                frame,
                identity,
                GeometryPolicy::default(),
                None
            )
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Input
    );
    for limit in 1..=3 {
        let calls = std::cell::Cell::new(0);
        let cancel = || {
            let n = calls.get() + 1;
            calls.set(n);
            n > limit
        };
        assert_eq!(
            original
                .rust_cut_cylindrical_blind_pocket(
                    spec(),
                    frame,
                    EntityId::new(),
                    GeometryPolicy::default(),
                    Some(&cancel)
                )
                .err()
                .unwrap()
                .kind,
            GeometryErrorKind::Cancelled
        );
        assert!(original.rust_cylindrical_blind_pocket_geometry().is_none());
    }
    assert_eq!(
        original
            .rust_cut_cylindrical_blind_pocket(
                CylindricalBlindPocketSpec {
                    radius_mm: 0.,
                    ..spec()
                },
                frame,
                EntityId::new(),
                GeometryPolicy::default(),
                None
            )
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Input
    );
    let invalid_policy = GeometryPolicy {
        maximum_coordinate_mm: 1e9,
        ..Default::default()
    };
    assert_eq!(
        original
            .rust_cut_cylindrical_blind_pocket(spec(), frame, EntityId::new(), invalid_policy, None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Input
    );
    use multicad_rust_kernel::tapered_thread::Error;
    for (error, kind) in [
        (Error::Limit, GeometryErrorKind::Limit),
        (Error::Invalid, GeometryErrorKind::Input),
        (Error::Cancelled, GeometryErrorKind::Cancelled),
        (Error::Unsupported, GeometryErrorKind::Unsupported),
        (Error::Unresolved, GeometryErrorKind::Unsupported),
    ] {
        assert_eq!(source_failure(error).kind, kind);
    }
}
#[test]
fn independently_edited_revolution_never_inherits_a_stale_receiving_certificate() {
    let cut = stock(PlaneFrame::default(), EntityId::new(), 10.)
        .rust_cut_cylindrical_blind_pocket(
            spec(),
            entry(PlaneFrame::default(), true),
            EntityId::new(),
            GeometryPolicy::default(),
            None,
        )
        .unwrap();
    let (source, _) = cut.rust_cylindrical_blind_pocket_geometry().unwrap();
    let mut points = source.geometry().rings()[0].clone();
    for p in &mut points {
        if p[1] == 6. {
            p[1] = 5.;
        }
    }
    let feature = EntityId::new();
    let profile = CurvedProfile {
        frame: PlaneFrame {
            v: [0., 0., 1.],
            ..Default::default()
        },
        outer: points
            .iter()
            .zip(points.iter().cycle().skip(1))
            .map(|(&a, &b)| CurvedEdge2d::Line { a, b })
            .collect(),
        holes: vec![],
    };
    let origins = (0..points.len())
        .map(|i| Origin {
            source: feature,
            role: format!("edited.{i}"),
        })
        .collect::<Vec<_>>();
    let edited = Solid::rust_revolve_curved_profile_named(
        &profile,
        RevolutionAxis {
            origin_mm: [0.; 3],
            direction: [0., 0., 1.],
        },
        std::f64::consts::TAU,
        feature,
        &origins,
        GeometryPolicy::default(),
    )
    .unwrap();
    assert!(edited.rust_cylindrical_blind_pocket_geometry().is_none());
    assert_eq!(source.spec().depth_mm, 4.);
    assert!(edited.volume_mm3().unwrap() < cut.volume_mm3().unwrap());
}

#[test]
fn a_loose_stock_policy_cannot_authorize_a_later_strict_cut_policy() {
    let loose = GeometryPolicy {
        linear_tolerance_mm: 1e-5,
        minimum_feature_mm: 0.001,
        ..Default::default()
    };
    let strict = GeometryPolicy::default();
    let scale = 1. + 4e-11;
    let frame = PlaneFrame {
        origin: [17., -8., 13.],
        u: [0.8 * scale, 0.6 * scale, 0.],
        v: [-0.36, 0.48, 0.8],
    };
    assert!(
        frame
            .checked_metric_extent([1000., 1000., 10.], loose)
            .is_ok()
    );
    assert!(
        frame
            .checked_metric_extent([1000., 1000., 10.], strict)
            .is_err()
    );
    let identity = EntityId::new();
    let cylinder = AxialCylinder::new(
        identity.uuid_bits(),
        Circle::new(
            [Rational::new(0, 1).unwrap(); 2],
            Rational::new(1000, 1).unwrap(),
        )
        .unwrap(),
        Rational::new(0, 1).unwrap(),
        Rational::new(10, 1).unwrap(),
    )
    .unwrap();
    let original = Solid::from_axial_cylinder_sweep(
        cylinder,
        frame,
        [0.; 2],
        false,
        Origin {
            source: identity,
            role: "actual.loose.stock".into(),
        },
        loose,
    )
    .unwrap();
    let rejected = original
        .rust_cut_cylindrical_blind_pocket(
            spec(),
            entry(frame, false),
            EntityId::new(),
            strict,
            None,
        )
        .err()
        .unwrap();
    assert_eq!(rejected.kind, GeometryErrorKind::Input);
    assert!(original.rust_cylindrical_blind_pocket_geometry().is_none());
    let admitted = original
        .rust_cut_cylindrical_blind_pocket(
            spec(),
            entry(frame, false),
            EntityId::new(),
            loose,
            None,
        )
        .unwrap();
    assert_eq!(
        admitted.rust_cylindrical_blind_pocket_geometry().unwrap().1,
        frame
    );
}

#[test]
fn literal_blind_pocket_step_reconstructs_the_independent_material_domain() {
    use multicad_rust_kernel::{revolution::Revolution, step, step_read};
    for frame in [
        PlaneFrame::default(),
        PlaneFrame {
            origin: [30., -20., 10.],
            u: [0.6, 0.8, 0.],
            v: [0., 0., 1.],
        },
        PlaneFrame {
            origin: [100000., -100000., 100000.],
            u: [0.6, 0.8, 0.],
            v: [0., 0., 1.],
        },
    ] {
        for end in [false, true] {
            for seated in [false, true] {
                let mut requested = spec();
                if !seated {
                    requested.seat = None;
                }
                let original = stock(frame, EntityId::new(), 10.);
                let cut = original
                    .rust_cut_cylindrical_blind_pocket(
                        requested,
                        entry(frame, end),
                        EntityId::new(),
                        GeometryPolicy::default(),
                        None,
                    )
                    .unwrap();
                let (source, pose) = cut.rust_cylindrical_blind_pocket_geometry().unwrap();
                assert!(step::supports_cylindrical_revolution(source.geometry()));
                let text = step::write_cylindrical_revolution(
                    source.geometry(),
                    pose.origin,
                    pose.u,
                    pose.v,
                )
                .unwrap();
                assert!(step_read::read_import_solid(&text, 1e-7).is_ok());
                assert!(cut.supports_step_export());
                let file = tempfile::NamedTempFile::new().unwrap();
                cut.write_step(file.path()).unwrap();
                let text = std::fs::read_to_string(file.path()).unwrap();
                let (read, imported) = step_read::read_import_solid(&text, 1e-7).unwrap();
                assert_eq!(read.solid_count, 1);
                assert_eq!(read.face_count, if seated { 7 } else { 5 });
                let expected_volume =
                    std::f64::consts::PI * (1000. - if seated { 40. } else { 16. });
                assert!((read.volume_mm3 - expected_volume).abs() < 1e-8);
                let step_read::ImportProfile::CylindricalRevolution { meridian } = imported.profile
                else {
                    panic!("actual pocket must retain its complete revolution boundary")
                };
                let decoded = Revolution::new(vec![meridian], std::f64::consts::TAU).unwrap();
                let normal = [
                    imported.u[1] * imported.v[2] - imported.u[2] * imported.v[1],
                    imported.u[2] * imported.v[0] - imported.u[0] * imported.v[2],
                    imported.u[0] * imported.v[1] - imported.u[1] * imported.v[0],
                ];
                for radius in [0., 1.9, 2.1, 3.9, 4.1, 9.9, 10.1] {
                    for z in [-0.1, 0.1, 1.9, 2.1, 3.9, 4.1, 9.9, 10.1] {
                        let travel = if end { 10. - z } else { z };
                        let cavity = travel > 0.
                            && travel < requested.depth_mm
                            && radius < requested.radius_mm
                            || seated && travel > 0. && travel < 2. && radius < 4.;
                        let material = radius < 10. && z > 0. && z < 10. && !cavity;
                        let world = position(pose, [radius, 0., z]);
                        let delta: [f64; 3] =
                            std::array::from_fn(|i| world[i] - imported.origin_mm[i]);
                        let local = [imported.u, imported.v, normal]
                            .map(|axis| (0..3).map(|i| delta[i] * axis[i]).sum::<f64>());
                        assert_eq!(
                            decoded.contains(local, 1e-7),
                            material,
                            "end={end} seat={seated} r={radius} z={z}"
                        );
                    }
                }
                let reopened =
                    Solid::import_rust_step(&text, EntityId::new(), GeometryPolicy::default())
                        .unwrap();
                assert!((reopened.volume_mm3().unwrap() - expected_volume).abs() < 1e-8);
                assert!(reopened.rust_cylindrical_blind_pocket_geometry().is_none());
                let proof = reopened
                    .manufacturing_mesh_with_receipt(Default::default())
                    .unwrap();
                assert_eq!(proof.edge_pairs as usize * 2, proof.mesh.indices.len());
                assert!(original.rust_cylindrical_blind_pocket_geometry().is_none());
            }
        }
    }
}
