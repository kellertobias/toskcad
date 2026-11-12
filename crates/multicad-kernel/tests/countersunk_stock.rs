use multicad_core::{document::EntityId, references::*};
use multicad_kernel::{
    rust_kernel::{Circle, Rational, chamfered_bores::Countersink, polygon::PolygonPrism},
    *,
};
fn source(depth: f64) -> Solid {
    let q = |n| Rational::new(n, 1).unwrap();
    let prism = PolygonPrism::with_circular_holes(
        vec![
            vec![[0, 0], [20, 0], [20, 16], [0, 16]]
                .into_iter()
                .map(|p| p.map(q))
                .collect(),
        ],
        vec![Circle::new([q(10), q(8)], q(1)).unwrap()],
        q(6),
    )
    .unwrap();
    let feature = EntityId::new();
    let origins: Vec<_> = (0..5)
        .map(|i| Origin {
            source: feature,
            role: format!("source.{i}"),
        })
        .collect();
    Solid::from_polygon_prism(
        prism,
        PlaneFrame {
            origin: [40., -8., 12.],
            u: [0.8, 0.6, 0.],
            v: [-0.36, 0.48, 0.8],
        },
        feature,
        &origins,
        &[false; 5],
        GeometryPolicy::default(),
    )
    .unwrap()
    .rust_countersink_polygon(
        &[Countersink {
            hole: 0,
            end: false,
            radius_mm: 3.,
            depth_mm: depth,
        }],
        &[Origin {
            source: feature,
            role: "seat".into(),
        }],
        GeometryPolicy::default(),
    )
    .unwrap()
}
fn encoded(source: &Solid) -> String {
    let dir = tempfile::tempdir().unwrap();
    let p = dir.path().join("cone.step");
    source.write_step(&p).unwrap();
    std::fs::read_to_string(p).unwrap()
}
#[test]
fn actual_conical_mouths_export_blind_import_and_reexport_with_closed_material_and_face_intents() {
    let policy = GeometryPolicy::default();
    for depth in [0.75, 1.25, 2.5] {
        let source = source(depth);
        let expected =
            (320. - std::f64::consts::PI) * 6. - std::f64::consts::PI * depth * (9. + 3. - 2.) / 3.;
        assert!((source.volume_mm3().unwrap() - expected).abs() < 1e-9);
        source
            .manufacturing_mesh_with_receipt(Default::default())
            .unwrap();
        let owner = EntityId::new();
        let text = encoded(&source);
        let imported = Solid::import_rust_step(&text, owner, policy).unwrap();
        assert_eq!(imported.solid_count().unwrap(), 1);
        assert_eq!(imported.faces().unwrap().len(), 8);
        assert!((imported.volume_mm3().unwrap() - expected).abs() < 1e-7);
        imported
            .manufacturing_mesh_with_receipt(Default::default())
            .unwrap();
        let before = imported.topology_snapshot(owner, 1).unwrap();
        let same = Solid::import_rust_step(&text, owner, policy)
            .unwrap()
            .topology_snapshot(owner, 2)
            .unwrap();
        let pose = imported
            .transformed(RigidTransform {
                translation_mm: [20., -10., 5.],
                ..Default::default()
            })
            .unwrap()
            .topology_snapshot(owner, 2)
            .unwrap();
        for face in before
            .elements
            .iter()
            .filter(|e| e.key.kind == TopologyKind::Face)
        {
            let intent = GeometryIntent::capture(&before, face.key).unwrap();
            for target in [&same, &pose] {
                assert!(matches!(
                    resolve(&intent, target, 2),
                    Resolution::Resolved { .. }
                ));
            }
        }
        let reimported =
            Solid::import_rust_step(&encoded(&imported), EntityId::new(), policy).unwrap();
        assert!((reimported.volume_mm3().unwrap() - expected).abs() < 1e-7);
    }
}
#[test]
fn an_unchanged_outer_wall_keeps_its_import_identity_when_cone_slope_crosses_45_degrees() {
    let owner = EntityId::new();
    let policy = GeometryPolicy::default();
    let a = Solid::import_rust_step(&encoded(&source(1.25)), owner, policy).unwrap();
    let snapshot = a.topology_snapshot(owner, 1).unwrap();
    let wall = GeometryIntent::capture(
        &snapshot,
        TopologyKey {
            kind: TopologyKind::Face,
            ordinal: 2,
        },
    )
    .unwrap();
    for depth in [2., 2.5] {
        let b = Solid::import_rust_step(&encoded(&source(depth)), owner, policy)
            .unwrap()
            .topology_snapshot(owner, 2)
            .unwrap();
        assert!(
            matches!(resolve(&wall, &b, 2), Resolution::Resolved { .. }),
            "depth{depth}"
        );
    }
}

#[test]
fn cap_loops_use_actual_conical_mouths_and_mixed_treatments_refuse_atomically() {
    let source = source(1.25);
    let owner = EntityId::new();
    let before = source.topology_snapshot(owner, 0).unwrap();
    let faces = source.faces().unwrap();
    for ordinal in [0, 1] {
        let element = before
            .elements
            .iter()
            .find(|e| {
                e.key
                    == TopologyKey {
                        kind: TopologyKind::Face,
                        ordinal,
                    }
            })
            .unwrap();
        let AnalyticGeometry::Plane { frame } = element.geometry.unwrap() else {
            panic!("cap must be planar")
        };
        let contours = source
            .planar_face_contours(
                ordinal,
                PlaneFrame {
                    origin: frame.origin_mm,
                    u: frame.u,
                    v: frame.v,
                },
            )
            .unwrap();
        assert_eq!(contours.loops.len(), 2);
        contours.validate(1e-7).unwrap();
        assert!(
            (contours.loops.iter().map(|l| l.signed_area()).sum::<f64>()
                - faces[ordinal as usize].area_mm2)
                .abs()
                < 1e-7
        );
    }
    let circle=before.elements.iter().find(|e|matches!(e.geometry,Some(AnalyticGeometry::CircularArc {radius_mm,..}) if radius_mm==3.)).unwrap();
    let error = source
        .treat_edges_named(
            &[TreatmentEdge {
                key: circle.key,
                origin: Origin {
                    source: EntityId::new(),
                    role: "mixed".into(),
                },
            }],
            EdgeTreatment::Chamfer { distance_mm: 0.1 },
            EntityId::new(),
            GeometryPolicy::default(),
        )
        .err()
        .unwrap();
    // Use one consistent treatment identity to exercise the material guard.
    assert_eq!(error.kind, GeometryErrorKind::Input);
    let feature = EntityId::new();
    let error = source
        .treat_edges_named(
            &[TreatmentEdge {
                key: circle.key,
                origin: Origin {
                    source: feature,
                    role: "mixed".into(),
                },
            }],
            EdgeTreatment::Chamfer { distance_mm: 0.1 },
            feature,
            GeometryPolicy::default(),
        )
        .err()
        .unwrap();
    assert_eq!(error.kind, GeometryErrorKind::Unsupported);
    assert_eq!(source.topology_snapshot(owner, 0).unwrap(), before);
}
