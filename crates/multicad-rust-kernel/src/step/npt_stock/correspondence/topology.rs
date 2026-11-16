//! Actual owner bijection, earned from source coedges and parsed literal graph.
use super::*;
use crate::{axial_pipe_thread::authored::StockClipStation, step_read::npt_stock::Support};
#[derive(Clone, Copy)]
pub(super) struct Topology {
    pub edges: [usize; 20],
    pub vertices: [usize; 12],
}
pub(super) fn qualify(
    plan: &SurfacePlan<'_>,
    owned: &OwnedMaterial,
    matching: [usize; 8],
    work: &mut Work<'_>,
) -> Result<Option<Topology>, Error> {
    let mut edges = [usize::MAX; 20];
    let mut edge_owners = [usize::MAX; 20];
    for (patch_index, tensor_index) in matching.into_iter().enumerate() {
        work.charge(32)?;
        let face = owned
            .faces()
            .iter()
            .find(|face| matches!(face.support(),Support::Tensor(index) if index==tensor_index))
            .expect("retained tensor owns one face");
        for (side, source) in plan.patches()[patch_index]
            .source_patch()
            .coedges()
            .into_iter()
            .enumerate()
        {
            work.charge(128)?;
            let index = plan
                .curves()
                .curves()
                .iter()
                .position(|curve| curve.source_curve() == source.curve)
                .expect("original source coedge");
            let target = face.coedges()[side];
            if source.forward != target.forward() {
                return Ok(None);
            }
            let target = target.edge_index();
            if (edges[index] != usize::MAX && edges[index] != target)
                || (edge_owners[target] != usize::MAX && edge_owners[target] != index)
            {
                return Ok(None);
            }
            edges[index] = target;
            edge_owners[target] = index;
        }
    }
    let mut vertices = [usize::MAX; 12];
    let mut vertex_owners = [usize::MAX; 12];
    for (index, target) in edges.into_iter().enumerate() {
        work.charge(256)?;
        assert!(
            target != usize::MAX,
            "all original atlas edges are incident"
        );
        let source = plan.curves().curves()[index].source_curve();
        for (owner, target) in source
            .endpoints()
            .into_iter()
            .zip(owned.edges()[target].vertex_indices())
        {
            let index = plan
                .curves()
                .vertices()
                .iter()
                .position(|v| v.owner() == owner)
                .expect("original source vertex");
            // Oriented face/edge bijection glues every connected original
            // vertex link. All actual vertices are edge endpoints; coverage
            // gives a surjection between two12-element registries, hence a
            // bijection. A conflict is an internal graph invariant break.
            assert!(
                (vertices[index] == usize::MAX || vertices[index] == target)
                    && (vertex_owners[target] == usize::MAX || vertex_owners[target] == index),
                "earned vertex-link bijection"
            );
            vertices[index] = target;
            vertex_owners[target] = index;
        }
        let expected = match source {
            StockClipCurve::Ring {
                station: StockClipStation::Start,
                ..
            } => Some(false),
            StockClipCurve::Ring {
                station: StockClipStation::End,
                ..
            } => Some(true),
            _ => None,
        };
        let mut actual = None;
        for face in owned.faces() {
            if matches!(face.support(), Support::Plane { .. })
                && face
                    .coedges()
                    .iter()
                    .any(|edge| edge.edge_index() == target)
            {
                actual = Some(face.forward());
            }
        }
        // Exact UV edge binding and matching orientation preserve the three
        // original station ranks. Actual z_v>0 makes vertical edges strictly
        // increasing; the global cap-height range forces low vertices to rank0
        // and high vertices to rank2. Material already earned their face senses.
        assert_eq!(actual, expected, "earned cap station and sense");
    }
    Ok(Some(Topology { edges, vertices }))
}
