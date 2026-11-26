use super::*;
use crate::geometry::pipe_hardware_body::tests::{frame, spec};

fn authored(rounded: bool, external: bool, posed: bool) -> Solid {
    let kind = if external {
        PipeHardwareKind::HexHead {
            across_flats_mm: 16.,
            height_mm: 3.,
        }
    } else {
        PipeHardwareKind::HexNut {
            across_flats_mm: 16.,
        }
    };
    Solid::rust_authored_axial_pipe_hardware_with_interrupt(
        spec(rounded, external, false),
        kind,
        0.2,
        EntityId::new(),
        frame(posed),
        GeometryPolicy::default(),
        None,
    )
    .unwrap()
}
fn import(text: &str) -> Solid {
    Solid::import_rust_step(text, EntityId::new(), GeometryPolicy::default()).unwrap()
}
fn check(rounded: bool, external: bool, posed: bool) {
    let source = authored(rounded, external, posed);
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("hardware.step");
    source.write_step(&path).unwrap();
    let text = std::fs::read_to_string(&path).unwrap();
    let imported = import(&text);
    let (kind, af, span) = imported.rust_imported_pipe_hex_data().unwrap();
    assert_eq!(
        kind,
        if external {
            HexStockKind::Head
        } else {
            HexStockKind::Nut
        }
    );
    assert!((af - 16.).abs() <= 1e-9, "{af}");
    let height = span[1] - span[0];
    let expected = if external {
        3.
    } else {
        spec(rounded, external, false).length_mm
    };
    assert!((height - expected).abs() <= 1e-9, "{height}");
    // The same material: volume, world bounds and face count agree.
    let area: f64 = source.faces().unwrap().iter().map(|f| f.area_mm2).sum();
    let v = source.volume_mm3().unwrap();
    let w = imported.volume_mm3().unwrap();
    assert!((v - w).abs() <= 1e-6 * v + 1e-7 * area, "{v} {w}");
    let a = source.bounds().unwrap();
    let b = imported.bounds().unwrap();
    // Authored bounds may re-box a placed thread box; the imported support
    // bounds lie inside them and agree exactly at identity.
    for k in 0..3 {
        assert!(a.min[k] - 1e-6 <= b.min[k] && b.max[k] <= a.max[k] + 1e-6);
        if !posed {
            assert!((a.min[k] - b.min[k]).abs() < 1e-6 && (a.max[k] - b.max[k]).abs() < 1e-6);
        }
    }
    let faces = imported.faces().unwrap();
    assert_eq!(faces.len(), source.faces().unwrap().len());
    let face_area: f64 = faces.iter().map(|f| f.area_mm2).sum();
    assert!(
        (face_area - area).abs() <= 1e-6 * area,
        "{face_area} {area}"
    );
    let topology = imported.topology_snapshot(EntityId::new(), 1).unwrap();
    topology.validate().unwrap();
    assert!(
        topology
            .elements
            .iter()
            .any(|e| e.origins.iter().any(|o| o.role.ends_with(".wall/0")))
    );
    assert!(!imported.display_edges(0.05).unwrap().is_empty());
    // The display mesh joins the hex to the decoded thread in the stock
    // chart at every pose; its faces are the reported faces.
    let mesh = imported.mesh(MeshOptions::default()).unwrap();
    assert!(!mesh.indices.is_empty());
    let ids: BTreeSet<u32> = mesh.face_indices.iter().copied().collect();
    let known: BTreeSet<u32> = faces.iter().map(|f| f.index).collect();
    assert_eq!(ids, known);
    // Display vertices (cap fan apices included) are within the requested
    // deflection of the exact material.
    let slack = MeshOptions::default().linear_deflection_mm;
    for p in mesh.vertices.as_chunks::<3>().0 {
        for k in 0..3 {
            assert!(
                b.min[k] - slack <= p[k] && p[k] <= b.max[k] + slack,
                "{p:?} {b:?}"
            );
        }
    }
    let manufacturing = imported
        .manufacturing_mesh_with_receipt(MeshOptions::default())
        .unwrap();
    assert_eq!(
        manufacturing.edge_pairs as usize * 2,
        manufacturing.mesh.indices.len()
    );
    // Re-export writes the retained tables beside the earned stock and
    // reads back to the same material; a moved import exports too.
    assert!(imported.supports_step_export());
    let again = dir.path().join("again.step");
    imported.write_step(&again).unwrap();
    let back = import(&std::fs::read_to_string(&again).unwrap());
    assert!((back.volume_mm3().unwrap() - w).abs() <= 1e-9 * w);
    assert_eq!(back.rust_imported_pipe_hex_data().unwrap().0, kind);
    let moved = imported
        .transformed(crate::RigidTransform {
            translation_mm: [5., 0., 0.],
            ..Default::default()
        })
        .unwrap();
    assert!((moved.volume_mm3().unwrap() - w).abs() <= 1e-9 * w);
    let moved_path = dir.path().join("moved.step");
    moved.write_step(&moved_path).unwrap();
    let moved_back = import(&std::fs::read_to_string(&moved_path).unwrap());
    let mb = moved_back.bounds().unwrap();
    assert!((mb.min[0] - b.min[0] - 5.).abs() < 1e-6);
}
#[test]
fn imported_bspt_hex_head_keeps_its_material_topology_mesh_and_re_exports() {
    check(true, true, false);
}
#[test]
fn imported_posed_npt_hex_nut_keeps_its_material_and_re_exports() {
    check(false, false, true);
}
#[test]
fn imported_posed_bspt_hex_head_meshes_in_its_stock_chart_and_re_exports() {
    check(true, true, true);
}
#[test]
fn imported_bspt_hex_nut_keeps_its_material_topology_mesh_and_re_exports() {
    check(true, false, false);
}
#[test]
fn compound_hex_stock_members_import_re_export_and_verify_member_by_member() {
    let policy = GeometryPolicy::default();
    let source = authored(false, false, false);
    let ordinary = Solid::rust_box(
        EntityId::new().uuid_bits(),
        [1, 2, 3].map(|n| Rational::new(n, 1).unwrap()),
        PlaneFrame {
            origin: [60., 0., 0.],
            ..Default::default()
        },
        policy,
    )
    .unwrap();
    let compound =
        Solid::make_rust_compound(&[&ordinary, &source], EntityId::new(), policy, None).unwrap();
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("compound.step");
    compound.write_step(&path).unwrap();
    let text = std::fs::read_to_string(&path).unwrap();
    let imported = import(&text);
    let members = imported
        .rust_compound()
        .unwrap()
        .members()
        .collect::<Vec<_>>();
    assert_eq!(members.len(), 2);
    let hex = members
        .iter()
        .find(|m| m.rust_imported_pipe_hex_data().is_some())
        .expect("one typed pipe hex member");
    let v = source.volume_mm3().unwrap();
    assert!((hex.volume_mm3().unwrap() - v).abs() <= 1e-6 * v);
    assert!((imported.volume_mm3().unwrap() - compound.volume_mm3().unwrap()).abs() <= 1e-6 * v);
    // The imported compound re-exports, and its file verifies member by
    // member: the hex member by exact table and corner readback.
    assert!(imported.supports_step_export());
    let again = dir.path().join("again.step");
    imported.write_step(&again).unwrap();
    let again_text = std::fs::read_to_string(&again).unwrap();
    let proof = imported
        .verify_partitioned_step_source(&again_text, 1e-7, None)
        .unwrap()
        .unwrap();
    assert_eq!(proof.members.solid_count, 2);
    assert_eq!(proof.members.retained_members.len(), 1);
    let back = import(&again_text);
    assert!((back.volume_mm3().unwrap() - imported.volume_mm3().unwrap()).abs() <= 1e-9 * v);
    // A moved member is a different source for the same file.
    let moved = imported
        .transformed(crate::RigidTransform {
            translation_mm: [1e-3, 0., 0.],
            ..Default::default()
        })
        .unwrap();
    assert!(
        moved
            .verify_partitioned_step_source(&again_text, 1e-7, None)
            .unwrap()
            .is_err()
    );
}
