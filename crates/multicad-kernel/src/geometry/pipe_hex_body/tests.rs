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
    for k in 0..3 {
        assert!((a.min[k] - b.min[k]).abs() < 1e-6 && (a.max[k] - b.max[k]).abs() < 1e-6);
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
    // The display mesh joins the hex to the decoded thread in an exact
    // chart; posed imports may refuse it as unsupported, never silently.
    match imported.mesh(MeshOptions::default()) {
        Ok(mesh) => {
            assert!(!mesh.indices.is_empty());
            let ids: BTreeSet<u32> = mesh.face_indices.iter().copied().collect();
            let known: BTreeSet<u32> = faces.iter().map(|f| f.index).collect();
            assert!(ids.is_subset(&known), "{ids:?} {known:?}");
        }
        Err(e) => {
            assert!(posed, "{e}");
            assert_eq!(e.kind, GeometryErrorKind::Unsupported, "{e}");
        }
    }
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
fn imported_bspt_hex_nut_keeps_its_material_topology_mesh_and_re_exports() {
    check(true, false, false);
}
#[test]
fn a_compound_member_hex_stock_import_is_refused_as_unsupported() {
    let source = authored(false, false, false);
    let ordinary = Solid::rust_box(
        EntityId::new().uuid_bits(),
        [1, 2, 3].map(|n| Rational::new(n, 1).unwrap()),
        PlaneFrame {
            origin: [60., 0., 0.],
            ..Default::default()
        },
        GeometryPolicy::default(),
    )
    .unwrap();
    let compound = Solid::make_rust_compound(
        &[&ordinary, &source],
        EntityId::new(),
        GeometryPolicy::default(),
        None,
    )
    .unwrap();
    let file = tempfile::NamedTempFile::new().unwrap();
    compound.write_step(file.path()).unwrap();
    let text = std::fs::read_to_string(file.path()).unwrap();
    let error = Solid::import_rust_step(&text, EntityId::new(), GeometryPolicy::default())
        .err()
        .unwrap();
    assert_eq!(error.kind, GeometryErrorKind::Unsupported, "{error}");
}
