use multicad_core::{
    document::{Document, Entity},
    model, project_file,
};
use multicad_kernel::Solid;
use multicad_model::{Compiler, GeometryBackend};

#[path = "support/axis_solids.rs"]
mod axis_solids;
#[path = "support/saved_prism.rs"]
mod saved_prism;

use axis_solids::{Aabb, Grid, Tool, operand};

/// Independent expectation from the saved document: the target's boxes (a
/// box primitive, or a join of box primitives and rectangle extrusions) less
/// the saved cylinder primitives, evaluated exactly.
fn expected(document: &Document, feature: multicad_core::document::EntityId) -> (f64, f64, Aabb) {
    let source = model::read_boolean(document, feature).unwrap();
    assert_eq!(source.operation, model::Operation::Cut);
    let producers = match model::read_boolean(document, source.target.producer) {
        Ok(join) => {
            assert_eq!(join.operation, model::Operation::Join);
            std::iter::once(join.target.producer)
                .chain(join.tools.iter().map(|t| t.producer))
                .collect()
        }
        Err(_) => vec![source.target.producer],
    };
    let material: Vec<Aabb> = producers
        .into_iter()
        .map(|id| match operand(document, id) {
            Tool::Box(b) => b,
            _ => panic!("the bored stock must be boxes"),
        })
        .collect();
    let grid = Grid::new(&material, &[], None);
    let [mut volume, mut area] = [grid.volume_mm3(), grid.area_mm2()];
    for tool in &source.tools {
        let Tool::Cylinder(bore) = operand(document, tool.producer) else {
            panic!("bores must be cylinders");
        };
        let (dv, da) = bore.cut(&grid);
        volume += dv;
        area += da;
    }
    (volume, area, grid.bounds())
}

/// Rebuild saved bore cuts and compare each with its independent expectation.
fn check_saved_cuts(
    mut project: &[u8],
    select: impl Fn(&Document, &Entity) -> bool,
    expected_count: usize,
) {
    let document = project_file::read(
        &mut project,
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap()
    .project
    .document;
    let features = document
        .entities()
        .filter(|entity| select(&document, entity))
        .map(|entity| entity.id)
        .collect::<Vec<_>>();
    assert_eq!(features.len(), expected_count);
    let compiler = Compiler {
        omit_display: true,
        ..Default::default()
    };
    let rust = compiler
        .compile_sources_with_backend(&document, 0, &features, GeometryBackend::Rust, None)
        .unwrap_or_else(|cause| panic!("Rust: {cause}"));
    let area = |solid: &Solid| {
        solid
            .faces()
            .unwrap()
            .iter()
            .map(|face| face.area_mm2)
            .sum::<f64>()
    };
    let near = |a: f64, b: f64| assert!((a - b).abs() < 1e-7 * b.abs().max(1.), "{a} != {b}");
    let temp = tempfile::tempdir().unwrap();
    for feature in features {
        let source = model::read_boolean(&document, feature).unwrap();
        let expected_bores = source.tools.len();
        let a = rust.solid(feature).unwrap();
        let (volume, surface, [min, max]) = expected(&document, feature);
        assert!(a.uses_rust_geometry());
        assert_eq!(a.solid_count().unwrap(), 1);
        near(a.volume_mm3().unwrap(), volume);
        near(area(a), surface);
        assert_eq!(
            a.faces()
                .unwrap()
                .iter()
                .filter(|face| face.plane.is_none())
                .count(),
            expected_bores
        );
        let a_bounds = a.bounds().unwrap();
        for axis in 0..3 {
            near(a_bounds.min[axis], min[axis]);
            near(a_bounds.max[axis], max[axis]);
        }
        let receipt = a
            .manufacturing_mesh_with_receipt(Default::default())
            .unwrap();
        assert!(!receipt.mesh.indices.is_empty());
        let path = temp.path().join(format!("{feature}.step"));
        a.write_step(&path).unwrap();
        for imported in [
            Solid::read_step(&path, compiler.policy).unwrap(),
            Solid::import_rust_step(
                &std::fs::read_to_string(&path).unwrap(),
                feature,
                compiler.policy,
            )
            .unwrap(),
        ] {
            near(imported.volume_mm3().unwrap(), a.volume_mm3().unwrap());
            near(area(&imported), area(a));
        }
    }
}

#[test]
fn front_panel_end_plates_rebuild_with_orthogonal_blind_and_through_bores() {
    check_saved_cuts(
        include_bytes!("../../../examples/front-panel.tcd"),
        |_, entity| {
            matches!(
                entity.name.as_str(),
                "Right end part holes" | "Left end part holes"
            )
        },
        2,
    );
}

#[test]
fn toolcase_slide_spacers_rebuild_with_side_through_bores() {
    check_saved_cuts(
        include_bytes!("../../../examples/toolcase.tcd"),
        |document, entity| {
            model::read_boolean(document, entity.id)
                .ok()
                .and_then(|source| document.get(source.target.producer))
                .is_some_and(|target| target.name.starts_with("Slide spacer "))
        },
        4,
    );
}
