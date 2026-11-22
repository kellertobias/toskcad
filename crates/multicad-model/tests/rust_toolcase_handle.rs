use multicad_core::model;
use multicad_core::{document::EntityId, project_file, references::AnalyticGeometry};
use multicad_kernel::{PlaneFrame, Solid};
use multicad_model::{Compiler, GeometryBackend};

#[path = "support/axis_solids.rs"]
mod axis_solids;
#[path = "support/saved_prism.rs"]
mod saved_prism;

use axis_solids::{Grid, Tool, operand};

#[test]
fn saved_toolcase_handle_joins_tangent_posts_into_one_manufacturable_solid() {
    let document = project_file::read(
        &mut include_bytes!("../../../examples/toolcase.tcd").as_slice(),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap()
    .project
    .document;
    let feature: EntityId =
        serde_json::from_str(r#""1d27a5a5-e8e2-4834-9868-8bfcdd4020ba""#).unwrap();
    let compiler = Compiler {
        omit_display: true,
        ..Default::default()
    };
    let rust = compiler
        .compile_sources_with_backend(&document, 0, &[feature], GeometryBackend::Rust, None)
        .unwrap();
    let solid = rust.solid(feature).unwrap();
    assert!(solid.uses_rust_geometry());
    assert_eq!(solid.solid_count().unwrap(), 1);
    let near = |a: f64, b: f64| assert!((a - b).abs() < 1e-7 * b.abs().max(1.), "{a} != {b}");
    // Independent expectation from the saved join: the box posts evaluated
    // exactly, with the saved grip cylinder joined in closed form.
    let source = model::read_boolean(&document, feature).unwrap();
    assert_eq!(source.operation, model::Operation::Join);
    let (mut posts, mut grips) = (vec![], vec![]);
    for id in std::iter::once(source.target.producer).chain(source.tools.iter().map(|t| t.producer))
    {
        match operand(&document, id) {
            Tool::Box(b) => posts.push(b),
            Tool::Cylinder(c) => grips.push(c),
            Tool::Prism(_) => panic!("unexpected prism"),
        }
    }
    assert_eq!([posts.len(), grips.len()], [2, 1]);
    let grid = Grid::new(&posts, &[], None);
    let (dv, da) = grips[0].join(&grid);
    let expected_volume = grid.volume_mm3() + dv;
    let expected_area = grid.area_mm2() + da;
    let [mut expected_min, mut expected_max] = grid.bounds();
    let grip = grips[0].bounds();
    for i in 0..3 {
        expected_min[i] = expected_min[i].min(grip[0][i]);
        expected_max[i] = expected_max[i].max(grip[1][i]);
    }
    near(solid.volume_mm3().unwrap(), expected_volume);
    let rust_faces = solid.faces().unwrap();
    let area = |faces: &[multicad_kernel::Face]| faces.iter().map(|face| face.area_mm2).sum();
    near(area(&rust_faces), expected_area);
    let bounds = solid.bounds().unwrap();
    for i in 0..3 {
        near(bounds.min[i], expected_min[i]);
        near(bounds.max[i], expected_max[i]);
    }
    assert!(
        solid
            .topology_snapshot(feature, 0)
            .unwrap()
            .elements
            .iter()
            .any(|element| matches!(element.geometry, Some(AnalyticGeometry::CircularArc { .. })))
    );
    let mut transitions = 0;
    for face in rust_faces {
        let Some((origin, normal)) = face.plane else {
            continue;
        };
        if normal[0].abs() != 1. || origin[0] == bounds.min[0] || origin[0] == bounds.max[0] {
            continue;
        }
        let frame = PlaneFrame {
            origin,
            u: [0., 1., 0.],
            v: [0., 0., normal[0]],
        };
        let contours = solid.planar_face_contours(face.index, frame).unwrap();
        assert_eq!(contours.loops.len(), 1);
        near(contours.loops[0].signed_area(), face.area_mm2);
        transitions += 1;
    }
    assert_eq!(transitions, 12);
    assert!(
        !solid
            .manufacturing_mesh_with_receipt(Default::default())
            .unwrap()
            .mesh
            .indices
            .is_empty()
    );
    let temp = tempfile::tempdir().unwrap();
    let path = temp.path().join("carry-handle.step");
    solid.write_step(&path).unwrap();
    let readers = [
        Solid::read_step(&path, compiler.policy).unwrap(),
        Solid::import_rust_step(
            &std::fs::read_to_string(&path).unwrap(),
            feature,
            compiler.policy,
        )
        .unwrap(),
    ];
    assert!(readers[1].uses_rust_geometry());
    for imported in readers {
        assert_eq!(imported.solid_count().unwrap(), 1);
        near(imported.volume_mm3().unwrap(), solid.volume_mm3().unwrap());
        near(area(&imported.faces().unwrap()), expected_area);
        let imported_bounds = imported.bounds().unwrap();
        for i in 0..3 {
            near(imported_bounds.min[i], bounds.min[i]);
            near(imported_bounds.max[i], bounds.max[i]);
        }
    }
}
