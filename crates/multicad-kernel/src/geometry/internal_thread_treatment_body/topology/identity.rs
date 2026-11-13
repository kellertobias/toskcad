//! Reference roles describe actual source constraints, never construction ordinals.
use super::*;

fn constraint(value: source::BoundaryConstraint) -> String {
    use source::BoundaryConstraint;
    match value {
        BoundaryConstraint::Rail(period) => format!("r{period}"),
        BoundaryConstraint::Row(row) => format!("a{row}"),
        BoundaryConstraint::Station(station, owner) => format!("z{station}.{owner}"),
    }
}
fn point_key(key: source::BoundaryPointKey) -> String {
    let mut constraints = key.0.map(constraint);
    constraints.sort();
    format!("{}-{}", constraints[0], constraints[1])
}
fn point_tag(point: &source::BoundaryPoint) -> String {
    use source::BoundaryPointRole;
    match point.role() {
        BoundaryPointRole::Source(key) => format!("p{}", point_key(key)),
        BoundaryPointRole::Tail { original, end } => format!("t{end}.{}", point_key(original)),
        BoundaryPointRole::StockRim { end } => format!("s{end}"),
    }
}
fn face_tag(key: source::FaceKey) -> String {
    use source::{FaceKey, ThreadBandRole};
    match key {
        FaceKey::StartCap => "a".into(),
        FaceKey::EndCap => "b".into(),
        FaceKey::EntryCone => "e".into(),
        FaceKey::TerminalShoulder => "t".into(),
        FaceKey::TerminalCylinder => "c".into(),
        FaceKey::MaterialOuterWall => "o".into(),
        FaceKey::Thread { start, band } => format!(
            "s{start}{}",
            match band {
                ThreadBandRole::Lower => 'l',
                ThreadBandRole::Crest => 'c',
                ThreadBandRole::Upper => 'u',
                ThreadBandRole::Root => 'r',
            }
        ),
    }
}

pub(super) struct Identities {
    pub edges: BTreeMap<u32, Origin>,
    pub vertices: BTreeMap<u32, Origin>,
}
pub(super) fn identities(
    boundary: &source::Boundary,
    feature: EntityId,
) -> Result<Identities, GeometryError> {
    let material = boundary
        .material()
        .ok_or_else(|| unavailable("Internal material topology requires its original stock"))?;
    let mut vertices = BTreeMap::new();
    let mut tags = BTreeMap::new();
    for (id, point) in boundary.vertices() {
        let tag = point_tag(point);
        vertices.insert(id, origin(feature, format!("thread/vertex/{tag}")));
        tags.insert(id, tag);
    }
    unique(&vertices)?;
    let mut groups = BTreeMap::<String, Vec<u32>>::new();
    for (id, curve) in boundary.curves() {
        // The source marks exactly every physical-curve endpoint as a vertex,
        // and qualifies both incident face IDs before publishing the boundary.
        let mut endpoints = curve.endpoints().map(|id| tags[&id].as_str());
        endpoints.sort();
        let mut faces = curve.faces().map(|id| {
            face_tag(
                material
                    .face_key(id)
                    .expect("qualified immutable source face incidence"),
            )
        });
        faces.sort();
        let key = format!(
            "thread/edge/{}.{}/{}|{}",
            faces[0], faces[1], endpoints[0], endpoints[1]
        );
        groups.entry(key).or_default().push(id);
    }
    let edges = edge_roles(groups, feature)?;
    Ok(Identities { edges, vertices })
}

fn edge_roles(
    groups: BTreeMap<String, Vec<u32>>,
    feature: EntityId,
) -> Result<BTreeMap<u32, Origin>, GeometryError> {
    let mut edges = BTreeMap::new();
    for (role, ids) in groups {
        if ids.len() != 1 {
            // Distinct finite fragments cannot inherit one semantic role merely
            // because their supports and endpoint constraints coincide.
            // Refuse ambiguity; no dimensional hash is mistaken for edit continuity.
            return Err(unavailable(
                "Finite source fragments need distinct semantic reference roles",
            ));
        }
        let id = ids[0];
        edges.insert(id, origin(feature, role));
    }
    unique(&edges)?;
    Ok(edges)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn distinct_finite_fragments_cannot_share_a_projected_reference_role() {
        let feature = EntityId::from_uuid_bits(8);
        let roles = BTreeMap::from([("thread/edge/projected".into(), vec![17, 31])]);
        assert!(
            edge_roles(roles, feature)
                .unwrap_err()
                .to_string()
                .contains("distinct semantic")
        );
    }
}
