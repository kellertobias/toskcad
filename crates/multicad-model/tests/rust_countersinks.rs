use multicad_core::{
    constraint_sketch::Dimension,
    document::EntityData,
    model::{self, EdgeTreatmentProfile},
    project_file,
    references::TopologyKind,
};
use multicad_kernel::{GeometryPolicy, MeshOptions, Solid};
use multicad_model::{Compiler, GeometryBackend};
use std::f64::consts::{PI, SQRT_2};

#[path = "support/saved_prism.rs"]
mod saved_prism;

/// Independent expectation from the saved document: a profile cabinetry
/// panel is its saved sketch region times the board thickness (surface: two
/// faces plus the edge band), a numeric transform only translates it, and
/// each 45° chamfer of size c on a saved bore circle of radius r removes the
/// Pappus ring π c² (r + c/3) and replaces the planar annulus and cylinder
/// strip by a cone. Returns volume, area, bounds and the region's circles.
fn expected(
    document: &multicad_core::document::Document,
    feature: multicad_core::document::EntityId,
) -> (f64, f64, [[f64; 3]; 2]) {
    /// Volume, area, bounds and the saved bore circles (curve id, radius).
    type Evaluated = (f64, f64, [[f64; 3]; 2], Vec<(String, f64)>);
    fn walk(
        document: &multicad_core::document::Document,
        feature: multicad_core::document::EntityId,
    ) -> Evaluated {
        let values = multicad_core::parameters::evaluate(document).unwrap();
        let mm = |d: &Dimension| {
            d.value_in(multicad_core::parameters::Unit::Millimeter, Some(&values))
                .unwrap()
        };
        if let Ok((definition, captures)) = model::read_edge_treatment(document, feature) {
            let (volume, area, bounds, circles) = walk(document, definition.input.producer);
            let c = mm(&definition.parameters.size);
            let (mut volume, mut area) = (volume, area);
            for capture in &captures {
                let r = capture
                    .origins
                    .iter()
                    .filter(|origin| origin.role == "profile.edge")
                    .find_map(|origin| {
                        circles
                            .iter()
                            .find(|(id, _)| *id == origin.source.to_string())
                            .map(|(_, r)| *r)
                    })
                    .expect("a chamfered saved bore circle");
                volume -= PI * c * c * (r + c / 3.);
                area += PI * (2. * r + c) * c * SQRT_2
                    - PI * ((r + c).powi(2) - r * r)
                    - 2. * PI * r * c;
            }
            return (volume, area, bounds, circles);
        }
        if let Ok(transform) = model::read_transform(document, feature) {
            let (volume, area, [lo, hi], circles) = walk(document, transform.input.producer);
            assert_eq!(
                transform
                    .parameters
                    .angle
                    .value_in(multicad_core::parameters::Unit::Degree, Some(&values))
                    .unwrap(),
                0.
            );
            let t = transform.parameters.translation.each_ref().map(mm);
            let shift = |p: [f64; 3]| std::array::from_fn(|i| p[i] + t[i]);
            return (volume, area, [shift(lo), shift(hi)], circles);
        }
        let multicad_core::document::EntityData::Feature { definition, .. } =
            &document.get(feature).unwrap().data
        else {
            panic!("{feature} is not a feature");
        };
        assert_eq!(definition.type_id, "cabinetry.panel");
        let data = serde_json::to_value(definition).unwrap()["data"].clone();
        let source = &data["source"];
        assert_eq!(source["kind"], "profile");
        let stock =
            multicad_core::document::EntityId::parse(source["stock"].as_str().unwrap()).unwrap();
        let thickness = model::stock_thickness(&model::read_panel_stock(document, stock).unwrap());
        let region = saved_prism::sketch_region(document, &source["profile"]);
        let circles = region
            .holes
            .iter()
            .filter_map(|hole| match hole {
                saved_prism::Boundary::Circle { id, radius, .. } => Some((id.clone(), *radius)),
                _ => None,
            })
            .collect();
        (
            region.area_mm2() * thickness,
            2. * region.area_mm2() + region.perimeter_mm() * thickness,
            region.swept_bounds(region.n.map(|x| x * thickness)),
            circles,
        )
    }
    let (volume, area, bounds, _) = walk(document, feature);
    (volume, area, bounds)
}

#[test]
fn toolcase_circular_bore_chamfers_rebuild_with_named_faces_and_exact_step() {
    let mut bytes = include_bytes!("../../../examples/toolcase.tcd").as_slice();
    let document = project_file::read(
        &mut bytes,
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap()
    .project
    .document;
    let sources = document
        .entities()
        .filter_map(|entity| {
            let EntityData::Feature {
                definition,
                suppressed: false,
                ..
            } = &entity.data
            else {
                return None;
            };
            (definition.type_id == model::EDGE_TREATMENT_TYPE).then_some(entity.id)
        })
        .collect::<Vec<_>>();
    assert_eq!(sources.len(), 4);
    let compile = |backend| {
        Compiler::default()
            .compile_sources_with_backend(&document, 0, &sources, backend, None)
            .unwrap_or_else(|cause| panic!("{backend:?}: {cause}"))
    };
    let rust = compile(GeometryBackend::Rust);
    let temporary = tempfile::tempdir().unwrap();
    let mut rims = 0;
    for (index, &source) in sources.iter().enumerate() {
        let actual = rust.solid(source).unwrap();
        assert!(actual.uses_rust_geometry());
        assert_eq!(actual.solid_count().unwrap(), 1);
        let (definition, _) = model::read_edge_treatment(&document, source).unwrap();
        assert_eq!(definition.parameters.profile, EdgeTreatmentProfile::Chamfer);
        let area = |solid: &Solid| {
            solid
                .faces()
                .unwrap()
                .iter()
                .map(|f| f.area_mm2)
                .sum::<f64>()
        };
        let (expected_volume, expected_area, [min, max]) = expected(&document, source);
        assert!((actual.volume_mm3().unwrap() - expected_volume).abs() < 1e-5);
        assert!((area(actual) - expected_area).abs() < 1e-5);
        let a = actual.bounds().unwrap();
        for axis in 0..3 {
            assert!((a.min[axis] - min[axis]).abs() < 1e-6);
            assert!((a.max[axis] - max[axis]).abs() < 1e-6);
        }
        let (_, captures) = model::read_edge_treatment(&document, source).unwrap();
        let topology = actual.topology_snapshot(source, 0).unwrap();
        for capture in captures {
            let role = format!("edge/{}", capture.id);
            assert_eq!(
                topology
                    .elements
                    .iter()
                    .filter(|element| {
                        element.key.kind == TopologyKind::Face
                            && element
                                .origins
                                .iter()
                                .any(|origin| origin.source == source && origin.role == role)
                    })
                    .count(),
                1,
            );
            rims += 1;
        }
        let receipt = actual
            .manufacturing_mesh_with_receipt(MeshOptions::default())
            .unwrap();
        assert!(receipt.edge_pairs > 0 && receipt.vertex_links > 0);
        let path = temporary.path().join(format!("chamfer-{index}.step"));
        actual.write_step(&path).unwrap();
        let policy = GeometryPolicy::default();
        for restored in [
            Solid::read_step(&path, policy).unwrap(),
            Solid::import_rust_step(&std::fs::read_to_string(&path).unwrap(), source, policy)
                .unwrap(),
        ] {
            assert_eq!(restored.solid_count().unwrap(), 1);
            assert!((restored.volume_mm3().unwrap() - expected_volume).abs() < 1e-5);
        }
    }
    assert_eq!(rims, 20);
}
