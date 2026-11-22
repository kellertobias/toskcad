//! Convert only a proved finite source graph; meshes never enter this handoff.
use super::*;
use multicad_core::references::{TopologyElement, TopologyKey, TopologyKind};
use std::collections::{BTreeMap, BTreeSet};

pub(super) fn assemble(
    boundary: &source::Boundary,
    mut faces: Vec<TopologyElement>,
    edge_origins: &BTreeMap<u32, Origin>,
    vertex_origins: &BTreeMap<u32, Origin>,
    owner: EntityId,
    revision: u64,
) -> Result<TopologySnapshot, GeometryError> {
    let (vertices, edges) = boundary.topology_counts();
    // Source options cap each table at8192 and face count below40.
    let count = vertices + edges + boundary.face_count() as usize;
    if count > 4096
        || faces.len() != boundary.face_count() as usize
        || edge_origins.len() != edges
        || vertex_origins.len() != vertices
    {
        return Err(unavailable(
            "Internal source topology exceeds its handoff budget",
        ));
    }
    let mut adjacent = BTreeMap::<u32, BTreeSet<Origin>>::new();
    let mut corners = BTreeMap::<u32, BTreeSet<Origin>>::new();
    for (index, face) in faces.iter_mut().enumerate() {
        let loops = boundary.face_loops(index as u32, None).map_err(failure)?;
        for coedge in loops.into_iter().flatten() {
            let id = coedge.curve as u32; // The sealed source table has <=8192 entries.
            let origin = edge_origins
                .get(&id)
                .ok_or_else(|| unavailable("Missing actual source curve identity"))?;
            face.boundary.insert(origin.clone());
            adjacent
                .entry(id)
                .or_default()
                .extend(face.origins.iter().cloned());
        }
    }
    let mut elements = Vec::with_capacity(count); // Resource-qualified <=4096 descriptors.
    elements.append(&mut faces);
    let mut faces = elements;
    for (ordinal, (id, curve)) in boundary.curves().enumerate() {
        // Every physical curve has two proved face uses, already visited above.
        let origin = &edge_origins[&id];
        for endpoint in curve.endpoints() {
            if !vertex_origins.contains_key(&endpoint) {
                return Err(unavailable(
                    "Actual source curve lacks its physical endpoint",
                ));
            }
            corners.entry(endpoint).or_default().insert(origin.clone());
        }
        faces.push(TopologyElement {
            key: TopologyKey {
                kind: TopologyKind::Edge,
                ordinal: ordinal as u32,
            },
            origins: BTreeSet::from([origin.clone()]),
            boundary: adjacent
                .remove(&id)
                .expect("immutable earlier face incidence pass"),
            frames: vec![],
            geometry: None,
        });
    }
    for (ordinal, (id, _)) in boundary.vertices().enumerate() {
        faces.push(TopologyElement {
            key: TopologyKey {
                kind: TopologyKind::Vertex,
                ordinal: ordinal as u32,
            },
            // Boundary vertices are exactly the endpoints of physical curves;
            // the preceding immutable edge pass checked their complete name map.
            origins: BTreeSet::from([vertex_origins[&id].clone()]),
            boundary: corners
                .remove(&id)
                .expect("immutable earlier physical endpoint pass"),
            frames: vec![],
            geometry: None,
        });
    }
    let result = TopologySnapshot {
        owner,
        revision,
        elements: faces,
    };
    result.validate().map_err(|e| unavailable(&e))?;
    Ok(result)
}
