//! Independent checks of a saved toolcase drawer Boolean against its own
//! operands, each compiled alone by the Rust kernel.
//!
//! These relations hold for any exact result, whatever provider built it:
//! a join spans exactly the union of its operand boxes, a cut and an
//! intersection stay inside their target or common box, and the volume lies
//! between the measure-theoretic bounds of its operation. A separately
//! tessellated manufacturing mesh must enclose the same volume to within its
//! linear deflection over the whole surface.
use multicad_core::{document::Document, document::EntityId, model};
use multicad_kernel::{MeshOptions, Solid};
use multicad_model::{Compiler, GeometryBackend};

pub fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-7 * b.abs().max(1.), "{a} != {b}");
}

pub fn area(solid: &Solid) -> f64 {
    solid
        .faces()
        .unwrap()
        .iter()
        .map(|face| face.area_mm2)
        .sum()
}

/// Volume enclosed by the manufacturing mesh, by the divergence theorem.
pub fn mesh_volume(solid: &Solid, options: MeshOptions) -> f64 {
    let mesh = solid.manufacturing_mesh_with_receipt(options).unwrap().mesh;
    let vertex = |index: u32| -> [f64; 3] {
        std::array::from_fn(|axis| mesh.vertices[3 * index as usize + axis])
    };
    mesh.indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|triangle| {
            let [a, b, c] = triangle.map(vertex);
            (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0]))
                / 6.
        })
        .sum()
}

/// Check `actual`, the Rust result of Boolean `feature`, against its
/// operands and an independent tessellation.
pub fn check_against_operands(
    document: &Document,
    revision: u64,
    feature: EntityId,
    compiler: &Compiler,
    actual: &Solid,
) {
    let definition = model::read_boolean(document, feature).unwrap();
    let operands: Vec<EntityId> = std::iter::once(definition.target.producer)
        .chain(definition.tools.iter().map(|tool| tool.producer))
        .collect();
    let compiled = compiler
        .compile_sources_with_backend(document, revision, &operands, GeometryBackend::Rust, None)
        .unwrap();
    let solids: Vec<&Solid> = operands
        .iter()
        .map(|&id| compiled.solid(id).unwrap())
        .collect();
    let volumes: Vec<f64> = solids.iter().map(|s| s.volume_mm3().unwrap()).collect();
    let boxes: Vec<_> = solids.iter().map(|s| s.bounds().unwrap()).collect();
    let bounds = actual.bounds().unwrap();
    let volume = actual.volume_mm3().unwrap();
    let slack = 1e-7 * volume.max(1.);
    let [target, tools] = [volumes[0], volumes[1..].iter().sum::<f64>()];
    let (expected_min, expected_max): ([f64; 3], [f64; 3]) = match definition.operation {
        model::Operation::Join => {
            // Every operand lies inside the union, so the volume is at
            // least the largest one and at most their sum.
            let largest = volumes.iter().copied().fold(0., f64::max);
            assert!(volume >= largest - slack && volume <= target + tools + slack);
            (
                std::array::from_fn(|i| boxes.iter().map(|b| b.min[i]).fold(f64::MAX, f64::min)),
                std::array::from_fn(|i| boxes.iter().map(|b| b.max[i]).fold(f64::MIN, f64::max)),
            )
        }
        model::Operation::Cut => {
            assert!(volume >= target - tools - slack && volume <= target + slack);
            (boxes[0].min, boxes[0].max)
        }
        model::Operation::Intersect => {
            let smallest = volumes.iter().copied().fold(f64::MAX, f64::min);
            assert!(volume > 0. && volume <= smallest + slack);
            (
                std::array::from_fn(|i| boxes.iter().map(|b| b.min[i]).fold(f64::MIN, f64::max)),
                std::array::from_fn(|i| boxes.iter().map(|b| b.max[i]).fold(f64::MAX, f64::min)),
            )
        }
        other => panic!("{other:?} is not a drawer Boolean"),
    };
    // A join's box is exactly the operand union. These saved cuts and the
    // coupon crop also keep the extreme material of their target or common
    // box, so their boxes match it exactly as well.
    for i in 0..3 {
        near(bounds.min[i], expected_min[i]);
        near(bounds.max[i], expected_max[i]);
    }
    let options = MeshOptions::default();
    let chordal = options.linear_deflection_mm * area(actual);
    let tessellated = mesh_volume(actual, options);
    assert!(
        (tessellated - volume).abs() <= chordal,
        "mesh {tessellated} exact {volume} bound {chordal}"
    );
}
