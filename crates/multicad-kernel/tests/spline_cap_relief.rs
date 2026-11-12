use multicad_core::{
    document::EntityId,
    references::{GeometryIntent, Origin, Resolution, TopologyKind, resolve},
};
use multicad_kernel::rust_kernel::Rational;
use multicad_kernel::*;
fn profile(z: f64) -> CurvedProfile {
    let points = [[3., 4.], [4., 5.], [5., 4.], [4., 3.]];
    let controls = [[3., 5.], [5., 5.], [5., 3.], [3., 3.]];
    let hole = (0..4)
        .map(|i| {
            CurvedEdge2d::BSpline(BSpline2d {
                degree: 2,
                poles: vec![points[i], controls[i], points[(i + 1) % 4]],
                knots: vec![0., 1.],
                multiplicities: vec![3, 3],
                weights: None,
            })
        })
        .collect();
    let p = [[2., 2.], [7., 2.], [7., 7.], [2., 7.]];
    CurvedProfile {
        frame: PlaneFrame {
            origin: [0., 0., z],
            ..Default::default()
        },
        outer: (0..4)
            .map(|i| CurvedEdge2d::Line {
                a: p[i],
                b: p[(i + 1) % 4],
            })
            .collect(),
        holes: vec![
            hole,
            vec![CurvedEdge2d::Circle {
                center: [6., 6.],
                radius_mm: 0.25,
            }],
        ],
    }
}
fn tool(profile: &CurvedProfile, height: f64, policy: GeometryPolicy) -> Solid {
    let feature = EntityId::new();
    let origins = profile
        .outer
        .iter()
        .chain(profile.holes.iter().flatten())
        .enumerate()
        .map(|(i, _)| Origin {
            source: feature,
            role: format!("curve.{i}"),
        })
        .collect::<Vec<_>>();
    Solid::rust_extrude_spline_profile_named(profile, height, feature, &origins, policy).unwrap()
}
fn stock(policy: GeometryPolicy) -> Solid {
    Solid::rust_box(
        EntityId::new().uuid_bits(),
        [20, 20, 5].map(|v| Rational::new(v, 1).unwrap()),
        Default::default(),
        policy,
    )
    .unwrap()
}
#[test]
fn all_counter_islands_retain_physical_material_faces_mesh_step_and_source_identity() {
    let policy = GeometryPolicy::default();
    let dir = tempfile::tempdir().unwrap();
    for operation in [BooleanOperation::Union, BooleanOperation::Subtract] {
        let blank = stock(policy);
        let height = if operation == BooleanOperation::Union {
            -1.
        } else {
            1.
        };
        let tool = tool(&profile(0.), height, policy);
        let feature = EntityId::new();
        let result =
            Solid::boolean_spline_cap_relief(&[&blank, &tool], operation, feature, policy, None)
                .unwrap()
                .unwrap();
        assert_eq!(blank.volume_mm3().unwrap(), 2000.);
        let sign = if operation == BooleanOperation::Union {
            1.
        } else {
            -1.
        };
        assert!(
            (result.volume_mm3().unwrap() - (2000. + sign * tool.volume_mm3().unwrap())).abs()
                < 1e-7
        );
        let (data, _) = result.rust_spline_cap_relief_data().unwrap();
        let enclosure = result.rust_spline_cap_relief_source_enclosure().unwrap();
        assert!(enclosure.0 > 0. && enclosure.1 >= 1e-6);
        let tight = MeshOptions {
            linear_deflection_mm: enclosure.0 / 2.,
            ..Default::default()
        };
        assert!(result.manufacturing_mesh_with_receipt(tight).is_err());
        assert!(result.display_edges(enclosure.0 / 2.).is_err());

        assert_eq!(data.motifs()[0].holes().len(), 2);
        let z = if operation == BooleanOperation::Union {
            -0.5
        } else {
            0.5
        };
        for (p, expected) in [
            ([2.5, 2.5, z], operation == BooleanOperation::Union),
            ([4., 4., z], operation == BooleanOperation::Subtract),
            ([6., 6., z], operation == BooleanOperation::Subtract),
            ([10., 10., 0.5], true),
        ] {
            assert_eq!(result.rust_spline_cap_relief_contains(p).unwrap(), expected);
        }
        let snapshot = result.topology_snapshot(feature, 1).unwrap();
        let captures=snapshot.elements.iter().filter(|e|e.key.kind==TopologyKind::Face).map(|f|{let intent=GeometryIntent::capture(&snapshot,f.key).unwrap();assert!(matches!(resolve(&intent,&snapshot,1),Resolution::Resolved{key,..} if key==f.key));intent}).collect::<Vec<_>>();
        assert_eq!(result.faces().unwrap().len(), data.boundary().faces.len());
        let display = result.display_edges(0.01).unwrap();
        assert_eq!(
            display.len(),
            snapshot
                .elements
                .iter()
                .filter(|e| e.key.kind == TopologyKind::Edge)
                .count()
        );
        result
            .manufacturing_mesh_with_receipt(Default::default())
            .unwrap();
        let path = dir.path().join("counter-relief.step");
        result.write_step(&path).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        let proof = rust_kernel::step_read::read_import_compound(&text, policy.linear_tolerance_mm)
            .unwrap();
        assert!(
            result.step_readback_matches(
                &snapshot,
                &proof
                    .iter()
                    .map(|(p, i)| (p.clone(), Some(i.clone())))
                    .collect::<Vec<_>>(),
                policy.linear_tolerance_mm
            )
        );
        let imported = Solid::import_rust_step(&text, EntityId::new(), policy).unwrap();
        assert!((imported.volume_mm3().unwrap() - result.volume_mm3().unwrap()).abs() < 1e-6);
        assert!(imported.rust_spline_cap_relief_data().is_some());
        assert_eq!(
            imported.rust_spline_cap_relief_source_enclosure(),
            Some((0., 0.))
        );
        imported.write_step(&path).unwrap();
        let pose = RigidTransform {
            rotation: [[0.8, -0.36, 0.48], [0.6, 0.48, -0.64], [0., 0.8, 0.6]],
            translation_mm: [80., 60., 20.],
        };
        let moved = result.transformed(pose).unwrap();
        assert_eq!(
            moved.rust_spline_cap_relief_source_enclosure(),
            Some(enclosure)
        );
        let moved_snapshot = moved.topology_snapshot(feature, 2).unwrap();
        for intent in &captures {
            assert!(matches!(
                resolve(intent, &moved_snapshot, 2),
                Resolution::Resolved { .. }
            ));
        }
        moved
            .manufacturing_mesh_with_receipt(Default::default())
            .unwrap();
        moved.write_step(&path).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        Solid::import_rust_step(&text, EntityId::new(), policy).unwrap();
    }
}
#[test]
fn literal_near_offsets_breakthrough_and_cancellation_refuse_without_changing_sources() {
    let policy = GeometryPolicy::default();
    let blank = stock(policy);
    let before = blank.topology_snapshot(EntityId::new(), 0).unwrap();
    for (z, height) in [(1e-11, 1.), (-1e-11, 1.), (0., 5.), (0., 6.)] {
        let tool = tool(&profile(z), height, policy);
        let error = Solid::boolean_spline_cap_relief(
            &[&blank, &tool],
            BooleanOperation::Subtract,
            EntityId::new(),
            policy,
            None,
        )
        .unwrap()
        .err()
        .unwrap();
        assert_eq!(error.kind, GeometryErrorKind::Unsupported);
    }
    let tool = tool(&profile(0.), 1., policy);
    let count = std::cell::Cell::new(0);
    let f = || {
        count.set(count.get() + 1);
        false
    };
    Solid::boolean_spline_cap_relief(
        &[&blank, &tool],
        BooleanOperation::Subtract,
        EntityId::new(),
        policy,
        Some(&f),
    )
    .unwrap()
    .unwrap();
    for stop in [0, 1, count.get() / 2, count.get() - 1] {
        let n = std::cell::Cell::new(0);
        let f = || {
            let at = n.get();
            n.set(at + 1);
            at == stop
        };
        let error = Solid::boolean_spline_cap_relief(
            &[&blank, &tool],
            BooleanOperation::Subtract,
            EntityId::new(),
            policy,
            Some(&f),
        )
        .unwrap()
        .err()
        .unwrap();
        assert_eq!(error.kind, GeometryErrorKind::Cancelled);
    }
    assert_eq!(
        before.elements,
        blank.topology_snapshot(before.owner, 0).unwrap().elements
    );
    Solid::boolean_spline_cap_relief(
        &[&blank, &tool],
        BooleanOperation::Subtract,
        EntityId::new(),
        policy,
        None,
    )
    .unwrap()
    .unwrap();
}

#[test]
fn authored_untouched_stock_cap_keeps_its_physical_attachment_chart() {
    let policy = GeometryPolicy::default();
    let owner = EntityId::new();
    let blank = Solid::box_named([20., 20., 5.], owner, policy)
        .unwrap()
        .with_feature_cap_roles(owner, "plate.start", "plate.end")
        .unwrap();
    let before = blank.topology_snapshot(owner, 0).unwrap();
    let origin = Origin {
        source: owner,
        role: "plate.end".into(),
    };
    let old = before
        .elements
        .iter()
        .find(|e| e.key.kind == TopologyKind::Face && e.origins.contains(&origin))
        .unwrap();
    let attachment = old
        .frames
        .iter()
        .find(|f| f.origin == origin)
        .unwrap()
        .frame;
    let mut intent = GeometryIntent::capture(&before, old.key).unwrap();
    intent.attachment = Some(multicad_core::references::Attachment {
        anchor: origin,
        offset_mm: [0.; 3],
        angle_rad: 0.,
    });
    let tool = tool(&profile(0.), 1., policy);
    let feature = EntityId::new();
    let result = Solid::boolean_spline_cap_relief(
        &[&blank, &tool],
        BooleanOperation::Subtract,
        feature,
        policy,
        None,
    )
    .unwrap()
    .unwrap();
    let after = result.topology_snapshot(feature, 1).unwrap();
    let Resolution::Resolved {
        key,
        frame: Some(frame),
        ..
    } = resolve(&intent, &after, 1)
    else {
        panic!("retained far cap needs no repair");
    };
    for (actual, expected) in frame
        .origin_mm
        .into_iter()
        .chain(frame.u)
        .chain(frame.v)
        .zip(
            attachment
                .origin_mm
                .into_iter()
                .chain(attachment.u)
                .chain(attachment.v),
        )
    {
        assert!((actual - expected).abs() < 1e-10);
    }
    let face = result
        .faces()
        .unwrap()
        .into_iter()
        .find(|f| f.index == key.ordinal)
        .unwrap();
    let (p, n) = face.plane.unwrap();
    assert!(
        (0..3)
            .map(|i| (frame.origin_mm[i] - p[i]) * n[i])
            .sum::<f64>()
            .abs()
            < policy.linear_tolerance_mm
    );
    let too_many = std::iter::once(&blank)
        .chain(std::iter::repeat_n(&tool, 193))
        .collect::<Vec<_>>();
    assert_eq!(
        Solid::boolean_spline_cap_relief(
            &too_many,
            BooleanOperation::Subtract,
            EntityId::new(),
            policy,
            None
        )
        .unwrap()
        .err()
        .unwrap()
        .kind,
        GeometryErrorKind::Limit
    );
}
