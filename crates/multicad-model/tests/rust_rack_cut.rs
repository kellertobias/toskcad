use multicad_core::{document::EntityId, project_file, references::AnalyticGeometry};
use multicad_kernel::Solid;
use multicad_model::{Compiler, GeometryBackend, export::verify_rust_step_shape};

#[path = "support/axis_solids.rs"]
mod axis_solids;
#[path = "support/saved_prism.rs"]
mod saved_prism;

use axis_solids::{Aabb, Cylinder, Grid, Prism, Tool, operand};
use multicad_core::{
    document::Document,
    model::{self, Operation},
};

fn disjoint(a: &Aabb, b: &Aabb) -> bool {
    (0..3).any(|i| a[1][i] <= b[0][i] || b[1][i] <= a[0][i])
}

/// Independent expectation for the saved rack: the joined cheeks, plates,
/// shelves and web are boxes, the trim recesses are removed boxes, the
/// dovetail grooves are straight convex prisms and the magnet and insert
/// seats are cylinders, all read from the saved document; an optional crop
/// box intersects the result.
fn expected_rack(document: &Document, rack: EntityId, crop: Option<EntityId>) -> (f64, f64, Aabb) {
    let cut = model::read_boolean(document, rack).unwrap();
    assert_eq!(cut.operation, Operation::Cut);
    let join = model::read_boolean(document, cut.target.producer).unwrap();
    assert_eq!(join.operation, Operation::Join);
    let material: Vec<Aabb> = std::iter::once(join.target.producer)
        .chain(join.tools.iter().map(|t| t.producer))
        .map(|id| match operand(document, id) {
            Tool::Box(b) => b,
            _ => panic!("the rack stock must be boxes"),
        })
        .collect();
    let crop = crop.map(|id| match operand(document, id) {
        Tool::Box(b) => b,
        _ => panic!("the coupon crop must be a box"),
    });
    let (mut removed, mut prisms, mut cylinders) = (vec![], vec![], vec![]);
    for tool in &cut.tools {
        match operand(document, tool.producer) {
            Tool::Box(b) => removed.push(b),
            Tool::Prism(p) => prisms.push(p),
            Tool::Cylinder(c) => cylinders.push(c),
        }
    }
    assert_eq!([removed.len(), prisms.len(), cylinders.len()], [4, 16, 16]);
    // The curved and prismatic cuts are separate from each other and from
    // the trim recesses, so their effects add.
    let shaped: Vec<Aabb> = prisms
        .iter()
        .map(Prism::bounds)
        .chain(cylinders.iter().map(Cylinder::bounds))
        .collect();
    for (i, a) in shaped.iter().enumerate() {
        assert!(removed.iter().all(|b| disjoint(a, b)));
        assert!(shaped[i + 1..].iter().all(|b| disjoint(a, b)));
    }
    let grid = Grid::new(&material, &removed, crop);
    let [mut volume, mut area] = [grid.volume_mm3(), grid.area_mm2()];
    for (dv, da) in prisms
        .iter()
        .map(|p| p.cut(&grid))
        .chain(cylinders.iter().map(|c| c.cut(&grid)))
    {
        volume += dv;
        area += da;
    }
    (volume, area, grid.bounds())
}

#[test]
fn saved_toolcase_rack_retains_dovetail_faces_and_all_sixteen_exact_blind_bores() {
    let document = project_file::read(
        &mut include_bytes!("../../../examples/toolcase.tcd").as_slice(),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap()
    .project
    .document;
    let feature = EntityId::parse("c5978ce8-9734-43e1-a9ed-14d1fdbeced5").unwrap();
    let compiler = Compiler {
        omit_display: true,
        ..Compiler::default()
    };
    let rust = compiler
        .compile_sources_with_backend(&document, 0, &[feature], GeometryBackend::Rust, None)
        .unwrap_or_else(|error| panic!("Rust: {error}"));
    let actual = rust.solid(feature).unwrap();
    assert!(actual.uses_rust_geometry());
    let near = |a: f64, b: f64| assert!((a - b).abs() < 1e-7 * b.abs().max(1.), "{a} != {b}");
    let area = |solid: &Solid| {
        solid
            .faces()
            .unwrap()
            .iter()
            .map(|face| face.area_mm2)
            .sum()
    };
    let topology = actual.topology_snapshot(feature, 0).unwrap();
    for radius in [2., 2.6] {
        assert_eq!(
            topology
                .elements
                .iter()
                .filter(|edge| matches!(
                    edge.geometry,
                    Some(AnalyticGeometry::CircularArc { radius_mm, .. }) if radius_mm == radius
                ))
                .count(),
            64
        );
    }
    assert!(
        !actual
            .manufacturing_mesh_with_receipt(Default::default())
            .unwrap()
            .mesh
            .indices
            .is_empty()
    );
    let scratch = tempfile::tempdir().unwrap();
    let path = scratch.path().join("toolcase-rack.step");
    actual.write_step(&path).unwrap();
    let [independent, reimported] = [
        Solid::read_step(&path, compiler.policy).unwrap(),
        Solid::import_rust_step(
            &std::fs::read_to_string(&path).unwrap(),
            feature,
            compiler.policy,
        )
        .unwrap(),
    ];
    let (volume, surface, [min, max]) = expected_rack(&document, feature, None);
    for solid in [actual, &independent, &reimported] {
        assert_eq!(solid.solid_count().unwrap(), 1);
        near(solid.volume_mm3().unwrap(), volume);
        near(area(solid), surface);
        let bounds = solid.bounds().unwrap();
        for i in 0..3 {
            near(bounds.min[i], min[i]);
            near(bounds.max[i], max[i]);
        }
    }
    let receipt = verify_rust_step_shape(
        &path,
        rust.solid(feature).unwrap(),
        &rust.snapshot().shapes[&feature],
        compiler.policy.linear_tolerance_mm,
    )
    .unwrap();
    near(receipt, actual.volume_mm3().unwrap());
}

#[test]
fn saved_rack_coupon_crops_one_blind_bore_into_an_exact_open_notch() {
    let document = project_file::read(
        &mut include_bytes!("../../../examples/toolcase.tcd").as_slice(),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap()
    .project
    .document;
    let feature = EntityId::parse("f7bf396b-be67-4d17-9bb3-f4e5e7d5ca29").unwrap();
    let compiler = Compiler {
        omit_display: true,
        ..Compiler::default()
    };
    let rust = compiler
        .compile_sources_with_backend(&document, 0, &[feature], GeometryBackend::Rust, None)
        .unwrap_or_else(|error| panic!("Rust: {error}"));
    let actual = rust.solid(feature).unwrap();
    assert!(actual.uses_rust_geometry());
    let near = |a: f64, b: f64| assert!((a - b).abs() < 1e-7 * b.abs().max(1.), "{a} != {b}");
    let area = |solid: &Solid| {
        solid
            .faces()
            .unwrap()
            .iter()
            .map(|face| face.area_mm2)
            .sum()
    };
    let topology = actual.topology_snapshot(feature, 0).unwrap();
    assert!(topology.elements.iter().any(|edge| matches!(
        edge.geometry,
        Some(AnalyticGeometry::CircularArc { radius_mm: 2., .. })
    )));
    assert!(!topology.elements.iter().any(|edge| matches!(
        edge.geometry,
        Some(AnalyticGeometry::CircularArc { radius_mm: 2.6, .. })
    )));
    assert!(
        !actual
            .manufacturing_mesh_with_receipt(Default::default())
            .unwrap()
            .mesh
            .indices
            .is_empty()
    );
    let scratch = tempfile::tempdir().unwrap();
    let path = scratch.path().join("rack-coupon.step");
    actual.write_step(&path).unwrap();
    let [independent, reimported] = [
        Solid::read_step(&path, compiler.policy).unwrap(),
        Solid::import_rust_step(
            &std::fs::read_to_string(&path).unwrap(),
            feature,
            compiler.policy,
        )
        .unwrap(),
    ];
    let coupon = model::read_boolean(&document, feature).unwrap();
    assert_eq!(coupon.operation, Operation::Intersect);
    assert_eq!(coupon.tools.len(), 1);
    let (volume, surface, [min, max]) = expected_rack(
        &document,
        coupon.tools[0].producer,
        Some(coupon.target.producer),
    );
    near(max[1], 30.);
    for solid in [actual, &independent, &reimported] {
        assert_eq!(solid.solid_count().unwrap(), 1);
        near(solid.volume_mm3().unwrap(), volume);
        near(area(solid), surface);
        let bounds = solid.bounds().unwrap();
        for i in 0..3 {
            near(bounds.min[i], min[i]);
            near(bounds.max[i], max[i]);
        }
    }
    let receipt = verify_rust_step_shape(
        &path,
        rust.solid(feature).unwrap(),
        &rust.snapshot().shapes[&feature],
        compiler.policy.linear_tolerance_mm,
    )
    .unwrap();
    near(receipt, actual.volume_mm3().unwrap());
}
