use super::*;
use multicad_rust_kernel::{helical_thread::Spec, tapered_thread::Thread};
use std::sync::OnceLock;

fn source() -> RadialBoundary {
    static SOURCE: OnceLock<RadialBoundary> = OnceLock::new();
    SOURCE
        .get_or_init(|| {
            Thread::new(
                Spec {
                    major_radius_mm: 4.,
                    pitch_mm: 1.25,
                    length_mm: 1.25,
                    depth_mm: 0.5,
                    crest_width_mm: 0.25,
                    lower_flank_deg: 12.,
                    upper_flank_deg: 18.,
                    starts: 1,
                    left_hand: false,
                    root_radius_mm: None,
                },
                0.03125,
                1e-7,
                1e-6,
            )
            .unwrap()
            .retained_boundary_with_interrupt(
                [[0.; 3], [1., 0., 0.], [0., 1., 0.]],
                1e-7,
                1e-6,
                None,
            )
            .unwrap()
            .boundary
        })
        .clone()
}
fn fixture() -> Solid {
    Solid::from_rust_radial_boundary(source(), EntityId::new(), GeometryPolicy::default()).unwrap()
}
fn options() -> MeshOptions {
    MeshOptions {
        linear_deflection_mm: 0.05,
        angular_deflection_rad: 0.4,
        maximum_triangles: 200000,
    }
}

fn qualify(input: multicad_rust_kernel::radial_tensor::BoundaryInput) -> RadialBoundary {
    multicad_rust_kernel::radial_tensor::BoundaryOwnership::new(input, 1e-7, None)
        .unwrap()
        .bind_surfaces(1e-7, None)
        .unwrap()
        .certify_axial(1e-7, 1e-6, None)
        .unwrap()
        .certify_caps(None)
        .unwrap()
        .certify_ordering(None)
        .unwrap()
        .certify_material(1e-7, None)
        .unwrap()
}

#[test]
fn retained_finite_identity_survives_table_reordering_and_occurrence_motion() {
    use multicad_core::references::{GeometryIntent, Resolution, TopologyKind, resolve};
    let original = source();
    let mut input = original.input().clone();
    let nv = input.vertices.len();
    let ne = input.edges.len();
    let nc = input.curves.len();
    let ns = input.surfaces.len();
    input.vertices.reverse();
    input.curves.reverse();
    input.surfaces.reverse();
    input.edges.reverse();
    input.faces.reverse();
    for edge in &mut input.edges {
        edge.vertices = edge.vertices.map(|v| nv - 1 - v);
        edge.curve = nc - 1 - edge.curve;
    }
    for face in &mut input.faces {
        face.surface = ns - 1 - face.surface;
        for (edge, _) in &mut face.outer {
            *edge = ne - 1 - *edge;
        }
        face.outer.rotate_left(1);
    }
    let reordered = qualify(input);
    assert_eq!(source_signature(&original), source_signature(&reordered));
    let compound = |boundary: &RadialBoundary| {
        let mut writer = multicad_rust_kernel::step::Writer::new();
        writer.append_radial_boundary(boundary, 1e-7).unwrap();
        writer
            .append_box(
                [2, 3, 4].map(|v| Rational::new(v, 1).unwrap()),
                [20., 0., 0.],
                [1., 0., 0.],
                [0., 1., 0.],
            )
            .unwrap();
        writer.finish().unwrap()
    };
    let feature = EntityId::new();
    let before =
        Solid::import_rust_step(&compound(&original), feature, Default::default()).unwrap();
    let after =
        Solid::import_rust_step(&compound(&reordered), feature, Default::default()).unwrap();
    fn member(solid: &Solid) -> &Solid {
        solid
            .rust_compound()
            .unwrap()
            .members()
            .find(|member| member.radial_boundary_body().is_some())
            .unwrap()
    }
    assert_eq!(
        member(&before).radial_boundary_body().unwrap().0.feature,
        member(&after).radial_boundary_body().unwrap().0.feature
    );
    let owner = EntityId::new();
    let old = before.topology_snapshot(owner, 0).unwrap();
    let next = after.topology_snapshot(owner, 1).unwrap();
    let cap=old.elements.iter().find(|e| e.key.kind==TopologyKind::Face && matches!(e.geometry,Some(AnalyticGeometry::Plane {frame}) if frame.normal()[2]>0.99 && frame.origin_mm[2]>1.)).unwrap();
    let edge = old
        .elements
        .iter()
        .find(|e| {
            e.key.kind == TopologyKind::Edge
                && matches!(e.geometry, Some(AnalyticGeometry::CircularArc { .. }))
        })
        .unwrap();
    let cap_intent = GeometryIntent::capture(&old, cap.key).unwrap();
    let edge_intent = GeometryIntent::capture(&old, edge.key).unwrap();
    let Resolution::Resolved {
        geometry: Some(AnalyticGeometry::Plane { frame }),
        ..
    } = resolve(&cap_intent, &next, 1)
    else {
        panic!("Reordered cap did not keep its finite source identity")
    };
    assert!((frame.origin_mm[2] - 1.25).abs() < 1e-7);
    assert!(matches!(
        resolve(&edge_intent, &next, 1),
        Resolution::Resolved {
            geometry: Some(AnalyticGeometry::CircularArc { .. }),
            ..
        }
    ));
    let transform = RigidTransform {
        rotation: [[0.8, -0.36, 0.48], [0.6, 0.48, -0.64], [0., 0.8, 0.6]],
        translation_mm: [100., -30., 50.],
    };
    let posed = after.transformed(transform).unwrap();
    let placed = posed.topology_snapshot(owner, 2).unwrap();
    let Resolution::Resolved {
        geometry: Some(AnalyticGeometry::Plane { frame }),
        ..
    } = resolve(&cap_intent, &placed, 2)
    else {
        panic!("Occurrence pose renamed its immutable source cap")
    };
    assert!((frame.origin_mm[2] - 50.75).abs() < 1e-7);
    let Resolution::Resolved {
        geometry: Some(AnalyticGeometry::CircularArc { frame, .. }),
        ..
    } = resolve(&edge_intent, &placed, 2)
    else {
        panic!("Occurrence pose renamed its immutable source edge")
    };
    let Some(AnalyticGeometry::CircularArc {
        frame: original_arc,
        ..
    }) = edge.geometry
    else {
        unreachable!()
    };
    let old = original_arc.origin_mm;
    let expected = [
        100. + 0.8 * old[0] - 0.36 * old[1] + 0.48 * old[2],
        -30. + 0.6 * old[0] + 0.48 * old[1] - 0.64 * old[2],
        50. + 0.8 * old[1] + 0.6 * old[2],
    ];
    assert!(
        frame
            .origin_mm
            .into_iter()
            .zip(expected)
            .all(|(a, b)| (a - b).abs() < 1e-7)
    );
    let mut input = original.input().clone();
    let tensor = input
        .surfaces
        .iter_mut()
        .find_map(|surface| {
            if let Surface::Tensor(tensor) = surface {
                Some(tensor)
            } else {
                None
            }
        })
        .unwrap();
    let mut poles = tensor.poles().to_vec();
    poles[tensor.dimensions()[1] + 1][0] += 1e-4;
    *tensor =
        multicad_rust_kernel::tapered_thread::RetainedTensor::new(tensor.spans(), poles).unwrap();
    let changed = qualify(input);
    assert_ne!(source_signature(&original), source_signature(&changed));
}

#[test]
fn retained_facade_keeps_actual_geometry_metrics_topology_queries_and_manufacture() {
    let solid = fixture();
    assert!(solid.uses_rust_geometry());
    let (boundary, frame) = solid.rust_radial_boundary().unwrap();
    assert_eq!(frame, PlaneFrame::default());
    assert_eq!(boundary.input(), source().input());
    let original = boundary.input().clone();
    let share = solid.share().unwrap();
    assert!(std::ptr::eq(
        boundary,
        share.rust_radial_boundary().unwrap().0
    ));
    let properties = solid.rust_radial_boundary_properties(None).unwrap();
    let volume = solid.volume_mm3().unwrap();
    assert!(properties.volume_mm3[0] <= volume && volume <= properties.volume_mm3[1]);
    assert_eq!(solid.faces().unwrap().len(), original.faces.len());
    assert_eq!(solid.solid_count().unwrap(), 1);
    assert_eq!(
        solid
            .rust_radial_boundary_classify([0., 0., 0.625], None)
            .unwrap(),
        Location::Inside
    );
    assert_eq!(
        solid
            .rust_radial_boundary_classify([10., 0., 0.625], None)
            .unwrap(),
        Location::Outside
    );
    assert_eq!(
        solid
            .rust_radial_boundary_classify(original.vertices[0], None)
            .unwrap(),
        Location::Unresolved
    );
    let topology = solid
        .topology_with_display(EntityId::new(), 7, Some(0.05))
        .unwrap();
    assert_eq!(topology.edges.len(), original.edges.len());
    assert_eq!(topology.vertices.len(), original.vertices.len());
    assert_eq!(
        solid.display_edges(0.05).unwrap().len(),
        original.edges.len()
    );
    assert_eq!(
        solid.display_vertices().unwrap().len(),
        original.vertices.len()
    );
    let mesh = solid.mesh(options()).unwrap();
    let manufactured = solid.manufacturing_mesh_with_receipt(options()).unwrap();
    assert_eq!(mesh.vertices, manufactured.mesh.vertices);
    assert_eq!(mesh.indices, manufactured.mesh.indices);
    assert!(manufactured.edge_pairs > 0 && manufactured.vertex_links > 0);
    assert!(
        solid
            .mesh(MeshOptions {
                maximum_triangles: 5_000_000,
                ..options()
            })
            .is_ok()
    );
    assert_eq!(solid.rust_radial_boundary().unwrap().0.input(), &original);
    assert_eq!(
        solid
            .rust_radial_boundary_properties(Some(&|| true))
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Cancelled
    );
    assert_eq!(
        solid.rust_radial_boundary_properties(None).unwrap(),
        properties
    );
}

#[test]
fn retained_occurrences_qualify_rounded_world_source_once_and_keep_distinct_caches() {
    let solid = fixture();
    let transform = RigidTransform {
        rotation: [[0.8, -0.36, 0.48], [0.6, 0.48, -0.64], [0., 0.8, 0.6]],
        translation_mm: [100., -30., 50.],
    };
    let posed = solid.transformed(transform).unwrap();
    let (body, frame) = posed.radial_boundary_body().unwrap();
    let first = body.placement(frame, None).unwrap();
    assert!(Rc::ptr_eq(&first, &body.placement(frame, None).unwrap()));
    assert!(std::ptr::eq(
        solid.rust_radial_boundary().unwrap().0,
        posed.rust_radial_boundary().unwrap().0
    ));
    assert_eq!(
        posed
            .rust_radial_boundary_classify([100.3, -30.4, 50.375], None)
            .unwrap(),
        Location::Inside
    );
    let properties = posed.rust_radial_boundary_properties(None).unwrap();
    let scalar = posed.volume_mm3().unwrap();
    assert!(properties.volume_mm3[0] <= scalar && scalar <= properties.volume_mm3[1]);
    let original = solid.rust_radial_boundary_properties(None).unwrap();
    assert!(
        properties.volume_mm3[0] <= original.volume_mm3[1]
            && properties.volume_mm3[1] >= original.volume_mm3[0]
    );
    assert!(posed.manufacturing_mesh(options()).is_ok());
    let other = solid
        .transformed(RigidTransform {
            translation_mm: [-50., 10., 3.],
            ..RigidTransform::default()
        })
        .unwrap();
    other.bounds().unwrap();
    assert!(Rc::ptr_eq(&first, &body.placement(frame, None).unwrap()));
    assert_eq!(
        body.placement(frame, Some(&|| true)).err().unwrap().kind,
        GeometryErrorKind::Cancelled
    );
    assert_eq!(
        source().input(),
        solid.rust_radial_boundary().unwrap().0.input()
    );
}

#[test]
fn retained_constructor_and_query_failures_do_not_publish_partial_or_stale_geometry() {
    let policy = GeometryPolicy::default();
    assert_eq!(
        Solid::from_rust_radial_boundary(source(), EntityId::nil(), policy)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Input
    );
    assert_eq!(
        Solid::from_rust_radial_boundary(
            source(),
            EntityId::new(),
            GeometryPolicy {
                boolean_fuzzy_mm: 1e-8,
                ..policy
            }
        )
        .err()
        .unwrap()
        .kind,
        GeometryErrorKind::Input
    );
    assert_eq!(
        Solid::from_rust_radial_boundary(
            source(),
            EntityId::new(),
            GeometryPolicy {
                linear_tolerance_mm: 1e-10,
                ..policy
            }
        )
        .err()
        .unwrap()
        .kind,
        GeometryErrorKind::Input
    );
    let coarse = multicad_rust_kernel::radial_tensor::BoundaryOwnership::new(
        source().input().clone(),
        1e-3,
        None,
    )
    .unwrap()
    .bind_surfaces(1e-3, None)
    .unwrap()
    .certify_axial(1e-3, 1e-4, None)
    .unwrap()
    .certify_caps(None)
    .unwrap()
    .certify_ordering(None)
    .unwrap()
    .certify_material(1e-3, None)
    .unwrap();
    assert_eq!(
        Solid::from_rust_radial_boundary(coarse, EntityId::new(), policy)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(
        Solid::from_rust_radial_boundary_with_interrupt(
            source(),
            EntityId::new(),
            policy,
            Some(&|| true)
        )
        .err()
        .unwrap()
        .kind,
        GeometryErrorKind::Cancelled
    );
    let solid = fixture();
    assert!(
        solid
            .transformed(RigidTransform {
                translation_mm: [1e7, 0., 0.],
                ..RigidTransform::default()
            })
            .is_err()
    );
    assert_eq!(
        solid
            .rust_radial_boundary_classify([f64::NAN, 0., 0.], None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Input
    );
    let ordinary = Solid::rust_box(
        EntityId::new().uuid_bits(),
        [1, 2, 3].map(|v| Rational::new(v, 1).unwrap()),
        PlaneFrame::default(),
        policy,
    )
    .unwrap();
    assert!(ordinary.rust_radial_boundary().is_none());
    assert_eq!(
        ordinary
            .rust_radial_boundary_properties(None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(
        ordinary
            .rust_radial_boundary_classify([0.; 3], None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Unsupported
    );
    let (body, _) = solid.radial_boundary_body().unwrap();
    assert!(body.planar_loops(u32::MAX).is_empty());
    assert!(
        body.sample_edge(
            PlaneFrame::default(),
            u32::MAX,
            Solid::rust_mesh_options(options())
        )
        .is_none()
    );
    assert!(solid.manufacturing_mesh(options()).is_ok());
}

#[test]
fn normal_step_import_admits_actual_retained_members_and_mixed_compounds_without_occt() {
    use multicad_rust_kernel::{step, step_read::SourceImport};
    let policy = GeometryPolicy::default();
    let feature = EntityId::new();
    let text = step::write_radial_boundary(&source(), 1e-7).unwrap();
    let imported = Solid::import_rust_step(&text, feature, policy).unwrap();
    assert!(imported.rust_radial_boundary().is_some());
    assert!(imported.manufacturing_mesh(options()).is_ok());
    let mut writer = step::Writer::new();
    writer.append_radial_boundary(&source(), 1e-7).unwrap();
    writer
        .append_box(
            [2, 3, 4].map(|v| Rational::new(v, 1).unwrap()),
            [20., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
        )
        .unwrap();
    let text = writer.finish().unwrap();
    let compound = Solid::import_rust_step(&text, feature, policy).unwrap();
    assert!(compound.uses_rust_geometry());
    assert_eq!(compound.solid_count().unwrap(), 2);
    assert!(compound.manufacturing_mesh(options()).is_ok());
    let path =
        std::env::temp_dir().join(format!("multicad-retained-facade-{}.step", EntityId::new()));
    compound.write_step(&path).unwrap();
    let exported = std::fs::read_to_string(&path).unwrap();
    let verified = compound
        .verify_retained_step_source(&exported, 1e-7)
        .unwrap()
        .unwrap();
    assert_eq!(verified.solid_count, 2);
    assert_eq!(verified.retained_members.len(), 1);
    assert_eq!(verified.decoded_scalar_members.len(), 1);
    assert!(verified.retained_members[0].decoded_volume_bounds_mm3[0] > 0.);
    std::fs::remove_file(path).unwrap();
    let members = multicad_rust_kernel::step_read::read_source_compound(&exported, 1e-7).unwrap();
    assert_eq!(
        members
            .iter()
            .filter(|m| matches!(m, SourceImport::Radial(_)))
            .count(),
        1
    );
    assert_eq!(
        members
            .iter()
            .filter(|m| matches!(m, SourceImport::Analytic(_)))
            .count(),
        1
    );
    let reread = Solid::import_rust_step(&exported, feature, policy).unwrap();
    assert_eq!(reread.solid_count().unwrap(), 2);
}
