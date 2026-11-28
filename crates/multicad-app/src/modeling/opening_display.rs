//! Door and window display elements for the 3D scene (CAD-429). They come
//! from each displayed wall's canonical source through
//! `ResolvedWall::opening_elements`, the same resolved openings the plan
//! symbols use, and are attached to the wall Part's instance. They are never
//! native solids and never manufacturing geometry.
use multicad_core::{
    document::{Document, EntityData, EntityId},
    model::{self, OpeningElement, OpeningPieceRole},
    references::{Origin, TopologyElement, TopologyKey, TopologyKind},
};
use multicad_model::ModelSnapshot;
use multicad_ui::viewport;
use std::{
    collections::{BTreeMap, BTreeSet},
    sync::Arc,
};

/// Frame and leaf: light timber. Glazing: pale blue.
pub(super) const FRAME_RGBA: [u8; 4] = [178, 136, 92, 255];
pub(super) const GLAZING_RGBA: [u8; 4] = [168, 208, 232, 255];

/// One mesh per door or window: a closed box per piece, one face per piece
/// with its role as provenance, glazing faces colored apart.
pub(super) fn element_geometry(element: &OpeningElement) -> Result<viewport::Geometry, String> {
    let mut vertices = vec![];
    let mut indices = vec![];
    let mut faces = vec![];
    let mut topology = vec![];
    let mut colors = BTreeMap::new();
    for (ordinal, piece) in element.pieces.iter().enumerate() {
        let ordinal = u32::try_from(ordinal).map_err(|_| "Too many element pieces")?;
        let c = piece.corners;
        let center: [f64; 3] = std::array::from_fn(|i| c.iter().map(|p| p[i]).sum::<f64>() / 8.);
        // Quads of the box, each listed once; winding is made outward.
        for quad in [
            [0, 2, 6, 4],
            [1, 3, 7, 5],
            [0, 1, 5, 4],
            [2, 3, 7, 6],
            [0, 1, 3, 2],
            [4, 5, 7, 6],
        ] {
            let p = quad.map(|k| c[k]);
            let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
            let (u, v) = (sub(p[1], p[0]), sub(p[2], p[0]));
            let n = [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ];
            let mid: [f64; 3] = std::array::from_fn(|i| p.iter().map(|q| q[i]).sum::<f64>() / 4.);
            let out = sub(mid, center);
            let outward = n[0] * out[0] + n[1] * out[1] + n[2] * out[2] >= 0.;
            let order = if outward { [0, 1, 2, 3] } else { [0, 3, 2, 1] };
            let first = u32::try_from(vertices.len()).map_err(|_| "Element mesh too large")?;
            vertices.extend(order.map(|k| p[k]));
            indices.extend([first, first + 1, first + 2, first, first + 2, first + 3]);
            faces.extend([ordinal, ordinal]);
        }
        topology.push(TopologyElement {
            key: TopologyKey {
                kind: TopologyKind::Face,
                ordinal,
            },
            origins: BTreeSet::from([Origin {
                source: element.opening,
                role: format!("display.{}", piece.role.label()),
            }]),
            boundary: BTreeSet::new(),
            frames: vec![],
            geometry: None,
        });
        if piece.role == OpeningPieceRole::Glazing {
            colors.insert(ordinal, GLAZING_RGBA);
        }
    }
    viewport::Geometry::new(element.opening, vertices, indices, faces, topology.into())?
        .with_face_colors(colors)
}

/// The display elements of every displayed wall Part in `scene`.
pub(super) fn elements(
    document: &Document,
    model: &ModelSnapshot,
    scene: &viewport::Scene,
) -> Result<Vec<viewport::DisplayElement>, String> {
    let hosts: BTreeSet<EntityId> = scene.reference_instances().iter().map(|i| i.id).collect();
    let walls: Vec<_> = model
        .bodies
        .keys()
        .filter(|body| hosts.contains(body))
        .filter_map(|body| match document.get(*body).map(|e| &e.data) {
            Some(EntityData::Part { features, .. }) => features
                .iter()
                .find(|f| model::read_wall(document, **f).is_ok())
                .map(|wall| (*body, *wall)),
            _ => None,
        })
        .collect();
    if walls.is_empty() {
        return Ok(vec![]);
    }
    let parameters = multicad_core::parameters::evaluate(document).map_err(|e| e.to_string())?;
    let mut out = vec![];
    for (body, feature) in walls {
        // A wall that does not resolve has no solid either; its own
        // evaluation diagnostic explains why.
        let Ok(resolved) = model::read_wall(document, feature)
            .and_then(|source| model::resolve_wall(document, &parameters, &source))
        else {
            continue;
        };
        let Ok(plan) = resolved.plan() else {
            continue;
        };
        for element in resolved.opening_elements(&plan) {
            out.push(viewport::DisplayElement {
                id: element.opening,
                host: body,
                geometry: Arc::new(element_geometry(&element)?),
                color: FRAME_RGBA,
            });
        }
    }
    Ok(out)
}
