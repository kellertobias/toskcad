//! CAD-426: delete a model entry through the same canonical closure as the
//! model tree. A delete that would also remove dependents is refused with
//! their list unless the caller opts into the cascade explicitly.
use super::{Request, ToolResult, commit, editable, id, revision_guard};
use crate::{
    error::ToolError,
    host::{Host, with},
};
use multicad_core::hierarchy::{self, Action};
use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::json;

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct DeleteEntityArgs {
    pub expected_revision: u64,
    /// ID of the part, feature, sketch, plane, instance, component or group.
    pub id: String,
    /// Also delete every entry that depends on this one (later features of
    /// its part history, Booleans that consume it, instances of it). Without
    /// it, a delete with dependents is refused and lists them.
    #[serde(default)]
    pub cascade: bool,
}

pub(super) fn delete_entity(host: &dyn Host, call: &Request, args: DeleteEntityArgs) -> ToolResult {
    let target = id(&args.id)?;
    let call = call.clone();
    with(host, move |session_target| {
        revision_guard(session_target, args.expected_revision)?;
        editable(session_target)?;
        let session = session_target.session();
        let plan = hierarchy::deletion_plan(session, target)
            .map_err(|message| ToolError::invalid(message).involving(target))?;
        let dependents: Vec<_> = plan
            .dependents
            .iter()
            .map(|d| json!({ "id": d.id, "name": d.name }))
            .collect();
        if !plan.dependents.is_empty() && !args.cascade {
            let mut error = ToolError::new(
                "has_dependents",
                format!(
                    "{} Pass cascade: true to delete them together.",
                    plan.refusal()
                ),
            )
            .involving(target);
            error.involved.extend(plan.dependents.iter().map(|d| d.id));
            return Err(error);
        }
        let action = if args.cascade {
            Action::DeleteWithDependents { id: target }
        } else {
            Action::Delete { id: target }
        };
        let prepared =
            hierarchy::prepare(session, args.expected_revision, call.actor.clone(), action)
                .map_err(|error| ToolError::from(error).involving(target))?;
        let removed: Vec<_> = plan.removed.iter().copied().collect();
        let revision = commit(session_target, &call, prepared.prepared)?;
        Ok(json!({
            "revision": revision,
            "deleted": target,
            "dependents": dependents,
            "removed": removed,
        }))
    })?
}
