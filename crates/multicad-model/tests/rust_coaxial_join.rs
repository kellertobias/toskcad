use multicad_core::{model, project_file};
use multicad_precision::PortableMath as _;

#[path = "support/axis_solids.rs"]
mod axis_solids;
#[path = "support/saved_prism.rs"]
mod saved_prism;

use axis_solids::{Tool, operand};
use multicad_kernel::{MeshOptions, RigidTransform, Solid};
use multicad_model::{Compiler, GeometryBackend};
use std::f64::consts::PI;

#[test]
fn saved_front_panel_button_union_matches_the_analytic_solid() {
    let mut bytes = include_bytes!("../../../examples/front-panel.tcd").as_slice();
    let opened = project_file::read(
        &mut bytes,
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap();
    let document = opened.project.document;
    let feature = document
        .entities()
        .find(|entity| entity.name == "Button union")
        .unwrap()
        .id;
    let compiler = Compiler::default();
    let compile = |backend| {
        compiler
            .compile_sources_with_backend(&document, 0, &[feature], backend, None)
            .unwrap()
    };
    let rust = compile(GeometryBackend::Rust);
    let actual = rust.solid(feature).unwrap();
    assert!(actual.uses_rust_geometry());
    assert_eq!(actual.solid_count().unwrap(), 1);
    // Independent expectation from the saved button, bushing and nut
    // cylinders: a union of coaxial cylinders is, slab by slab along the
    // axis, a disc of the largest covering radius; its surface is every slab
    // wall plus the exposed annulus at every step in radius.
    let source = model::read_boolean(&document, feature).unwrap();
    assert_eq!(source.operation, model::Operation::Join);
    let cylinders: Vec<_> = std::iter::once(source.target.producer)
        .chain(source.tools.iter().map(|t| t.producer))
        .map(|id| match operand(&document, id) {
            Tool::Cylinder(c) => c,
            _ => panic!("the button parts must be cylinders"),
        })
        .collect();
    assert_eq!(cylinders.len(), 3);
    let (axis, center) = (cylinders[0].axis, cylinders[0].center);
    assert!(
        cylinders
            .iter()
            .all(|c| (c.axis, c.center) == (axis, center))
    );
    let mut stations: Vec<f64> = cylinders.iter().flat_map(|c| c.span).collect();
    stations.sort_by(f64::total_cmp);
    stations.dedup();
    let radius = |a: f64, b: f64| {
        cylinders
            .iter()
            .filter(|c| c.span[0] <= a && b <= c.span[1])
            .map(|c| c.radius)
            .fold(0., f64::max)
    };
    let slabs: Vec<f64> = stations.windows(2).map(|w| radius(w[0], w[1])).collect();
    assert!(slabs.iter().all(|r| *r > 0.), "one connected union");
    let mut expected_volume = 0.;
    let mut expected_area = 0.;
    for (i, r) in slabs.iter().enumerate() {
        let length = stations[i + 1] - stations[i];
        expected_volume += PI * r * r * length;
        expected_area += 2. * PI * r * length;
    }
    for i in 0..=slabs.len() {
        let below = if i == 0 { 0. } else { slabs[i - 1] };
        let above = slabs.get(i).copied().unwrap_or(0.);
        expected_area += PI * (above * above - below * below).abs();
    }
    let volume = actual.volume_mm3().unwrap();
    assert!((volume - 1055. * PI).abs() < 1e-8);
    assert!((volume - expected_volume).abs() < 1e-8);
    let actual_bounds = actual.bounds().unwrap();
    let largest = slabs.iter().copied().fold(0., f64::max);
    let others: Vec<usize> = (0..3).filter(|&i| i != axis).collect();
    let mut expected_bounds = multicad_kernel::Bounds {
        min: [0.; 3],
        max: [0.; 3],
    };
    [expected_bounds.min[axis], expected_bounds.max[axis]] =
        [stations[0], *stations.last().unwrap()];
    for (k, i) in others.into_iter().enumerate() {
        [expected_bounds.min[i], expected_bounds.max[i]] =
            [center[k] - largest, center[k] + largest];
    }
    for axis in 0..3 {
        assert!((actual_bounds.min[axis] - expected_bounds.min[axis]).abs() < 1e-6);
        assert!((actual_bounds.max[axis] - expected_bounds.max[axis]).abs() < 1e-6);
    }
    assert_eq!(actual.faces().unwrap().len(), 9);
    let mesh = actual
        .manufacturing_mesh_with_receipt(MeshOptions::default())
        .unwrap();
    assert!(!mesh.mesh.indices.is_empty());
    assert!(actual.supports_step_export());
    let scratch = tempfile::tempdir().unwrap();
    let path = scratch.path().join("button.step");
    actual.write_step(&path).unwrap();
    let independent = Solid::read_step(&path, compiler.policy).unwrap();
    assert_eq!(independent.solid_count().unwrap(), 1);
    assert_eq!(independent.faces().unwrap().len(), 9);
    assert!((independent.volume_mm3().unwrap() - expected_volume).abs() < 1e-7);
    let area = |solid: &Solid| {
        solid
            .faces()
            .unwrap()
            .iter()
            .map(|face| face.area_mm2)
            .sum::<f64>()
    };
    assert!((area(actual) - expected_area).abs() < 1e-7);
    assert!((area(&independent) - expected_area).abs() < 1e-7);
    let independent_bounds = independent.bounds().unwrap();
    for axis in 0..3 {
        assert!((independent_bounds.min[axis] - expected_bounds.min[axis]).abs() < 1e-6);
        assert!((independent_bounds.max[axis] - expected_bounds.max[axis]).abs() < 1e-6);
    }
    let (sin, cos) = 0.57_f64.portable_sin_cos();
    let oblique = actual
        .transformed(RigidTransform {
            rotation: [[cos, -sin, 0.], [sin, cos, 0.], [0., 0., 1.]],
            translation_mm: [12., -3., 5.],
        })
        .unwrap();
    let path = scratch.path().join("placed-button.step");
    oblique.write_step(&path).unwrap();
    let imported = Solid::import_rust_step(
        &std::fs::read_to_string(path).unwrap(),
        feature,
        compiler.policy,
    )
    .unwrap();
    assert_eq!(imported.faces().unwrap().len(), 9);
    assert!((imported.volume_mm3().unwrap() - volume).abs() < 1e-7);
    assert!((area(&imported) - area(&oblique)).abs() < 1e-7);
}
