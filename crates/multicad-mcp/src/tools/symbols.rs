//! Thin controls over canonical pinned definitions and owned instance overrides.
use super::{
    Quantity, Request, ToolResult, commit, compile, id, imports, revision_guard, snapshot,
};
use crate::{
    error::ToolError,
    host::{Host, with},
};
use multicad_core::{
    commands::Prepared,
    document::{EntityId, Placement},
    project_file, symbols,
};
use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::{Value, json};
use std::collections::BTreeMap;

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ImportSymbolArgs {
    pub expected_revision: u64,
    pub file_name: String,
    /// Required when the file contains more than one reusable definition.
    pub source_definition: Option<String>,
}
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SymbolValueArgs {
    pub value: f64,
    pub unit: String,
}
#[derive(Default, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SymbolPlacementArgs {
    /// Three typed lengths; omission places the underhead/root datum at zero.
    pub translation: Option<[Quantity; 3]>,
    /// Unit quaternion in x/y/z/w order; omission is identity.
    pub rotation_xyzw: Option<[f64; 4]>,
}
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PlaceSymbolArgs {
    pub expected_revision: u64,
    pub definition: String,
    pub name: String,
    /// Complete map keyed by the exposed IDs in the pinned template.
    pub overrides: BTreeMap<String, SymbolValueArgs>,
    #[serde(default)]
    pub placement: SymbolPlacementArgs,
}
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PlaceLibrarySymbolArgs {
    pub expected_revision: u64,
    /// Bare .tcd file name in the authorized exchange folder.
    pub file_name: String,
    /// Required when the source file has multiple reusable definitions.
    pub source_definition: Option<String>,
    pub name: String,
    /// Complete exposed-template override map; empty inherits pinned defaults.
    pub overrides: BTreeMap<String, SymbolValueArgs>,
    #[serde(default)]
    pub placement: SymbolPlacementArgs,
    /// Canonical receiving container; omission explicitly chooses project root.
    pub parent: Option<String>,
}
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EditSymbolInstanceArgs {
    pub expected_revision: u64,
    /// Placed instance ID or its symbol ownership record ID.
    pub instance: String,
    /// Replacement map, not a patch. Empty restores all pinned defaults.
    pub overrides: BTreeMap<String, SymbolValueArgs>,
}
fn overrides(
    values: BTreeMap<String, SymbolValueArgs>,
) -> Result<BTreeMap<multicad_core::document::EntityId, symbols::SymbolValue>, ToolError> {
    values
        .into_iter()
        .map(|(key, value)| {
            if !value.value.is_finite() {
                return Err(ToolError::invalid("Symbol override must be finite"));
            }
            Ok((
                id(&key)?,
                symbols::SymbolValue {
                    value: value.value,
                    unit: value.unit,
                },
            ))
        })
        .collect()
}
fn finish(
    host: &dyn Host,
    call: &Request,
    prepared: Prepared,
    revision: u64,
    mut result: Value,
) -> ToolResult {
    compile(prepared.document(), revision, call)?;
    let owned = call.clone();
    result["revision"] = json!(with(host, move |target| {
        // Shared commit validates the prepared owner session and revision
        // before either trusted commit or review queueing.
        commit(target, &owned, prepared)
    })??);
    result["native_model_preflight"] = json!(true);
    Ok(result)
}
fn library_definition(
    host: &dyn Host,
    file_name: &str,
    source_definition: Option<&str>,
) -> Result<(EntityId, symbols::SymbolDefinition), ToolError> {
    let path = imports::scoped(host, file_name, &[".tcd"])?
        .canonicalize()
        .map_err(|e| ToolError::invalid(e.to_string()))?;
    let folder = host
        .scope()
        .export_directory
        .ok_or_else(|| ToolError::new("out_of_scope", "No exchange folder"))?
        .canonicalize()
        .map_err(|e| ToolError::invalid(e.to_string()))?;
    if path.parent() != Some(folder.as_path()) {
        return Err(ToolError::new(
            "out_of_scope",
            "Symbol library must remain inside the exchange folder",
        ));
    }
    let project = project_file::read(
        &mut std::fs::File::open(path).map_err(|e| ToolError::invalid(e.to_string()))?,
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .map_err(|e| ToolError::invalid(e.to_string()))?
    .project;
    let selected = source_definition.map(id).transpose()?;
    let mut definitions: Vec<_> = symbols::definitions(&project.document)?
        .into_iter()
        .filter(|(key, _)| selected.is_none_or(|id| id == *key))
        .collect();
    if definitions.len() != 1 {
        return Err(ToolError::invalid(
            "Choose exactly one source definition from the library file",
        ));
    }
    Ok(definitions.remove(0))
}

pub(super) fn import_symbol(host: &dyn Host, call: &Request, args: ImportSymbolArgs) -> ToolResult {
    let (source_id, definition) =
        library_definition(host, &args.file_name, args.source_definition.as_deref())?;
    let pin = definition.pin()?;
    let symbol = definition.symbol;
    let exposed = definition.exposed.clone();
    let revision = args.expected_revision;
    let owned = call.clone();
    let prepared = with(host, move |target| {
        revision_guard(target, revision)?;
        Ok::<_, ToolError>(symbols::prepare_import(
            target.session(),
            revision,
            owned.actor.clone(),
            &definition,
        )?)
    })??;
    let key = prepared.definition;
    finish(
        host,
        call,
        prepared.prepared,
        revision,
        json!({"definition":key,"source_definition":source_id,"symbol":symbol,"pin":pin,"exposed":exposed,"native_geometry_checked":false,"geometry_scope":"Pinned source imported; place an instance to evaluate its geometry"}),
    )
}
pub(super) fn place_symbol(host: &dyn Host, call: &Request, args: PlaceSymbolArgs) -> ToolResult {
    let revision = args.expected_revision;
    let definition = id(&args.definition)?;
    let values = overrides(args.overrides)?;
    let owned = call.clone();
    let (prepared, record, instance, member) = with(host, move |target| {
        revision_guard(target, revision)?;
        let session = target.session();
        let mut placement = Placement::default();
        if let Some(translation) = args.placement.translation {
            placement.translation_mm = [
                translation[0]
                    .dimension(session.document(), session.parameters())?
                    .1,
                translation[1]
                    .dimension(session.document(), session.parameters())?
                    .1,
                translation[2]
                    .dimension(session.document(), session.parameters())?
                    .1,
            ];
        }
        if let Some(rotation) = args.placement.rotation_xyzw {
            placement.rotation_xyzw = rotation;
        }
        let request = symbols::SymbolPlaceRequest::new(
            session.document(),
            definition,
            args.name,
            values,
            placement,
        )?;
        let prepared = symbols::prepare_place(session, revision, owned.actor.clone(), &request)?;
        Ok::<_, ToolError>((
            prepared.prepared,
            request.record,
            request.instance,
            request.member,
        ))
    })??;
    finish(
        host,
        call,
        prepared,
        revision,
        json!({"definition":definition,"record":record,"instance":instance,"member":member,"native_geometry_checked":true}),
    )
}
pub(super) fn place_library_symbol(
    host: &dyn Host,
    call: &Request,
    args: PlaceLibrarySymbolArgs,
) -> ToolResult {
    call.check()?;
    let (source_id, content) =
        library_definition(host, &args.file_name, args.source_definition.as_deref())?;
    call.check()?;
    let pin = content.pin()?;
    let symbol = content.symbol;
    let parent = args.parent.as_deref().map(id).transpose()?;
    let values = overrides(args.overrides)?;
    let revision = args.expected_revision;
    let owned = call.clone();
    let (prepared, definition, record, instance, member) = with(host, move |target| {
        revision_guard(target, revision)?;
        let session = target.session();
        let definition = symbols::library_placement_definition(session.document(), &content)?
            .unwrap_or_else(EntityId::new);
        let mut placement = Placement::default();
        if let Some(translation) = args.placement.translation {
            placement.translation_mm = [
                translation[0]
                    .dimension(session.document(), session.parameters())?
                    .1,
                translation[1]
                    .dimension(session.document(), session.parameters())?
                    .1,
                translation[2]
                    .dimension(session.document(), session.parameters())?
                    .1,
            ];
        }
        if let Some(rotation) = args.placement.rotation_xyzw {
            placement.rotation_xyzw = rotation;
        }
        let request = symbols::SymbolPlaceRequest::from_definition(
            definition, &content, args.name, values, placement,
        )?;
        let prepared = symbols::prepare_library_place(
            session,
            revision,
            owned.actor.clone(),
            &content,
            &request,
            parent,
        )?;
        Ok::<_, ToolError>((
            prepared.prepared,
            definition,
            request.record,
            request.instance,
            request.member,
        ))
    })??;
    finish(
        host,
        call,
        prepared,
        revision,
        json!({"definition":definition,"source_file":args.file_name,"source_definition":source_id,"symbol":symbol,"pin":pin,"record":record,"instance":instance,"member":member,"parent":parent,"native_geometry_checked":true}),
    )
}
pub(super) fn edit_symbol_instance(
    host: &dyn Host,
    call: &Request,
    args: EditSymbolInstanceArgs,
) -> ToolResult {
    let revision = args.expected_revision;
    let record = id(&args.instance)?;
    let request = symbols::SymbolEditRequest {
        record,
        overrides: overrides(args.overrides)?,
    };
    let owned = call.clone();
    let prepared = with(host, move |target| {
        revision_guard(target, revision)?;
        Ok::<_, ToolError>(symbols::prepare_edit(
            target.session(),
            revision,
            owned.actor.clone(),
            &request,
        )?)
    })??;
    let instances = prepared.instances.clone();
    finish(
        host,
        call,
        prepared.prepared,
        revision,
        json!({"instances":instances,"overrides_semantics":"Replacement map; omitted entries inherit pinned defaults","native_geometry_checked":true}),
    )
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SymbolPlanArgs {
    /// Embedded definition ID, placed instance ID or symbol ownership record ID.
    pub symbol: String,
}
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SymbolPlanChoiceArgs {
    /// A sketch captured in the definition template (see symbol_plan).
    pub sketch: String,
    /// Shift along the symbol's local X and Y in millimeters; omission is zero.
    pub offset_mm: Option<[f64; 2]>,
    /// Counter-clockwise rotation about local Z, seen from above, applied
    /// before the offset; omission is zero.
    pub rotation_deg: Option<f64>,
}
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EditSymbolPlanArgs {
    pub expected_revision: u64,
    /// Embedded definition ID, placed instance ID or symbol ownership record ID.
    pub symbol: String,
    /// New explicit revision label; omission picks the next free one.
    pub revision: Option<String>,
    /// The fixed plan symbol; null removes it so plans project the 3D model.
    pub plan: Option<SymbolPlanChoiceArgs>,
}
/// The embedded definition an ID names directly or through a placed instance.
fn definition_of(
    doc: &multicad_core::document::Document,
    value: EntityId,
) -> Result<EntityId, ToolError> {
    if symbols::read_definition(doc, value).is_ok() {
        return Ok(value);
    }
    if let Ok(instance) = symbols::read_instance(doc, value) {
        return Ok(instance.definition);
    }
    for entity in doc.entities() {
        if let Ok(instance) = symbols::read_instance(doc, entity.id)
            && instance.instance == value
        {
            return Ok(instance.definition);
        }
    }
    Err(ToolError::invalid(
        "Name an embedded symbol definition, placed symbol instance or ownership record",
    ))
}
pub(super) fn symbol_plan(host: &dyn Host, _: &Request, args: SymbolPlanArgs) -> ToolResult {
    let (doc, revision, _) = snapshot(host)?;
    let definition = definition_of(&doc, id(&args.symbol)?)?;
    let value = symbols::read_definition(&doc, definition)?;
    let sketches: Vec<_> = value
        .template
        .entities()
        .filter(|e| matches!(e.data, multicad_core::document::EntityData::Sketch { .. }))
        .map(|e| json!({"sketch":e.id,"name":e.name}))
        .collect();
    let instances: Vec<_> = symbols::instances(&doc, definition)?
        .into_iter()
        .map(|(record, instance)| {
            let plan = value
                .plan
                .as_ref()
                .and_then(|plan| instance.ids.get(&plan.sketch));
            json!({"record":record,"instance":instance.instance,"plan_sketch":plan})
        })
        .collect();
    Ok(json!({
        "revision": revision,
        "definition": definition,
        "symbol": value.symbol,
        "definition_revision": value.revision,
        "plan": value.plan.as_ref().map(|plan| json!({
            "sketch": plan.sketch,
            "offset_mm": plan.alignment.offset_mm,
            "rotation_deg": plan.alignment.rotation_deg,
        })),
        "template_sketches": sketches,
        "instances": instances,
        "semantics": "Plan projections draw the plan sketch, rotated about local Z then offset, at each instance's pose instead of the detailed 3D content; 3D views and picking keep the model.",
    }))
}
pub(super) fn edit_symbol_plan(
    host: &dyn Host,
    call: &Request,
    args: EditSymbolPlanArgs,
) -> ToolResult {
    let revision = args.expected_revision;
    let symbol = id(&args.symbol)?;
    let plan = args
        .plan
        .map(|plan| -> Result<_, ToolError> {
            Ok(symbols::SymbolPlan::aligned(
                id(&plan.sketch)?,
                symbols::PlanAlignment {
                    offset_mm: plan.offset_mm.unwrap_or_default(),
                    rotation_deg: plan.rotation_deg.unwrap_or_default(),
                },
            ))
        })
        .transpose()?;
    let label = args.revision;
    let owned = call.clone();
    let (prepared, definition, report) = with(host, move |target| {
        revision_guard(target, revision)?;
        let doc = target.session().document();
        let definition = definition_of(doc, symbol)?;
        let current = symbols::read_definition(doc, definition)?;
        let label = match label {
            Some(label) => label,
            None => symbols::next_revision(doc, current.symbol, &current.revision)?,
        };
        // The same typed request the Inspector's plan editor applies.
        let request = symbols::SymbolUpdateRequest::plan_edit(doc, definition, label, plan)?;
        let report = symbols::compare_update_request(doc, &request)?;
        let prepared =
            symbols::prepare_update(target.session(), revision, owned.actor.clone(), &request)?;
        Ok::<_, ToolError>((prepared, definition, report))
    })??;
    let instances = prepared.instances.clone();
    finish(
        host,
        call,
        prepared.prepared,
        revision,
        json!({"definition":definition,"instances":instances,"old_pin":report.old_pin.to_string(),"new_pin":report.new_pin.to_string(),"native_geometry_checked":true}),
    )
}
