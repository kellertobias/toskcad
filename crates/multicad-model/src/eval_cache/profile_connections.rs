//! Only original joining-plate outputs have this occurrence-aware cache audit.
use super::*;
#[derive(Clone)]
pub(super) struct PoseInputs {
    pub fingerprint: ContentHash,
    pub geometry: BTreeSet<EntityId>,
    pub omitted: BTreeSet<EntityId>,
}
pub(super) fn source(document: &Document, id: EntityId) -> bool {
    matches!(document.get(id).map(|entity|&entity.data),Some(EntityData::Extension{definition,..})if definition.type_id==model::PROFILE_CONNECTION_TYPE)
}
pub(super) fn output(document: &Document, id: EntityId) -> bool {
    model::read_profile_connection_output(document, id)
        .ok()
        .is_some_and(|output| {
            output.access_bores.is_none()
                && matches!(output.geometry, model::ProfileOutputKind::Part { .. })
                && model::read_profile_connection_request(document, output.connection)
                    .ok()
                    .is_some_and(|request| {
                        request.kind == model::ProfileConnectionKind::JoiningPlate
                    })
        })
}
pub(super) fn project(document: &Document, id: EntityId) -> Option<PoseInputs> {
    let request = model::read_profile_connection_request(document, id).ok()?;
    if request.kind != model::ProfileConnectionKind::JoiningPlate {
        return None;
    }
    let mut geometry = BTreeSet::new();
    let mut omitted = BTreeSet::new();
    let mut bytes = b"multicad.profile-connection-absolute-pose/1".to_vec();
    for site in [&request.first, &request.second] {
        site.member.resolve_producer(document).ok()?;
        model::validate_transfer_path(document, &site.occurrence, site.member.body).ok()?;
        // Each COMPLETE absolute occurrence is projected separately. Cancelling
        // a common parent would incorrectly reuse old world-space hardware.
        let projection =
            model::occurrence_pose_projection(document, &site.occurrence, None).ok()?;
        bytes.extend(projection.fingerprint.0);
        geometry.extend(projection.geometry_dependencies);
        for id in
            std::iter::once(site.occurrence.root).chain(site.occurrence.members.iter().copied())
        {
            if matches!(
                document.get(id).map(|entity| &entity.data),
                Some(EntityData::Instance { .. } | EntityData::Assembly { .. })
            ) {
                omitted.insert(id);
            }
        }
    }
    Some(PoseInputs {
        fingerprint: ContentHash::of(&bytes),
        geometry,
        omitted,
    })
}
