//! Source-bound first phase of an extrusion cut. The second phase prepares the
//! final composite transaction on the original Session, never on these tools.
use crate::viewport::{Scene, SceneStamp};
use crate::{controls, egui};
use multicad_core::{
    commands::{Actor, Prepared, Session},
    document::{Document, EntityData, EntityId},
    model::{
        self, AutomaticOperation, BodyInput, CutTargetProbe, CutTargetQuery, CutTargetSelection,
        CutoutFace, EvaluationOptions, ExtrudeBooleanRequest, ExtrudeRequest, Operation,
    },
    project_file::ContentHash,
};
use std::collections::{BTreeMap, BTreeSet};

fn sheet_support(document: &Document, sketch: EntityId) -> bool {
    let Some((_, link)) = multicad_core::sketch_regions::supporting_face(document, sketch) else {
        return false;
    };
    document.entities().any(|entity| {
        matches!(&entity.data, EntityData::Part { features, .. } if features.contains(&link.target))
            && model::sheet_feature(document, entity.id).is_some()
    })
}

pub(super) struct SheetCutTarget {
    pub input: BodyInput,
    pub face: CutoutFace,
}

pub(super) struct OperationControls {
    pub operation: Operation,
    pub keep_tools: bool,
    pub component: Option<EntityId>,
    pub include: BTreeSet<EntityId>,
    pub exclude: BTreeSet<EntityId>,
    pub targets: BTreeSet<EntityId>,
    choices: Vec<(EntityId, String, bool)>,
    selection: Option<CutTargetSelection>,
    fresh: bool,
    native: BTreeMap<EntityId, BodyInput>,
    context_error: Option<String>,
    sheet_sketch: Option<EntityId>,
    developed: bool,
    pub sheet_target: Result<SheetCutTarget, String>,
    /// The operation the automatic default last set (CAD-404). The author's
    /// own choice of any other operation ends automatic mode for this form;
    /// distance and direction edits never change an explicit operation.
    automatic: Option<Operation>,
    decision: Option<AutomaticOperation>,
    names: BTreeMap<EntityId, String>,
}
impl OperationControls {
    pub fn new(document: &Document, sketch: Option<EntityId>) -> Self {
        let index = multicad_core::hierarchy::Index::new(document).ok();
        let mut choices: Vec<_> = document
            .entities()
            .filter_map(|entity| match entity.data {
                EntityData::Part { .. } => {
                    let path = index
                        .as_ref()
                        .map(|index| {
                            crate::model_tree::container_path(
                                document,
                                index,
                                index.parent(entity.id).ok().flatten(),
                            )
                        })
                        .unwrap_or_else(|| "Unknown owner".into());
                    let hidden = index.as_ref().is_some_and(|index| {
                        index.effective_visibility(entity.id).ok() == Some(false)
                    });
                    Some((
                        entity.id,
                        format!(
                            "{} · {path}{}",
                            entity.name,
                            if hidden { " · hidden" } else { "" }
                        ),
                        true,
                    ))
                }
                EntityData::Assembly { .. } => Some((
                    entity.id,
                    format!(
                        "{} (component)",
                        index
                            .as_ref()
                            .map(|index| crate::model_tree::container_path(
                                document,
                                index,
                                Some(entity.id)
                            ))
                            .unwrap_or_else(|| format!("{} · {}", entity.name, entity.id))
                    ),
                    false,
                )),
                _ => None,
            })
            .collect();
        choices.sort_by(|a, b| a.1.cmp(&b.1).then(a.0.cmp(&b.0)));
        let developed = sketch.is_some_and(|id| {
            model::developed_sheet_sketch(document, id).is_ok_and(|chart| chart.is_some())
        });
        let sheet_sketch = sketch.filter(|id| developed || sheet_support(document, *id));
        let names = document
            .entities()
            .filter(|entity| matches!(entity.data, EntityData::Part { .. }))
            .map(|entity| (entity.id, entity.name.clone()))
            .collect();
        Self {
            operation: if developed {
                Operation::Cut
            } else {
                Operation::NewBody
            },
            keep_tools: false,
            component: document.hierarchy().active_component,
            include: BTreeSet::new(),
            exclude: BTreeSet::new(),
            targets: BTreeSet::new(),
            choices,
            selection: None,
            fresh: false,
            native: BTreeMap::new(),
            context_error: Some(
                "Wait for a source-current native model before choosing explicit targets".into(),
            ),
            sheet_sketch,
            developed,
            sheet_target: Err("Wait for the current sheet to resolve its sketch support".into()),
            // Sheet sketches keep their dedicated through-sheet behaviour.
            automatic: (!developed && sheet_sketch.is_none()).then_some(Operation::NewBody),
            decision: None,
            names,
        }
    }
    pub fn set_context(
        &mut self,
        document: &Document,
        revision: u64,
        scene: &Scene,
    ) -> Result<(), String> {
        self.native.clear();
        self.sheet_target = Err("Wait for the current sheet to resolve its sketch support".into());
        let result: Result<(), String> = (|| {
            let info = scene
                .model_evaluation()
                .ok_or("Operation targets require a source-bound native evaluation")?;
            if info.source_hash != ContentHash::document(document)? {
                return Err("Operation targets belong to different source content".into());
            }
            if let Some(sketch) = self.sheet_sketch {
                self.sheet_target = (|| {
                    let (body, face) = if let Some(chart) =
                        model::developed_sheet_sketch(document, sketch)
                            .map_err(|error| error.to_string())?
                    {
                        (chart.part, chart.face)
                    } else {
                        crate::sketch_plane::attached_sheet_face(document, sketch, scene)?.ok_or(
                            "Sketch on a current sheet base or final flange leg before cutting",
                        )?
                    };
                    let producer = scene
                        .source_producer(body)
                        .filter(|producer| {
                            model::sheet_profile_head(document, body) == Some(*producer)
                        })
                        .ok_or("The sheet has later features that cannot be developed")?;
                    Ok(SheetCutTarget {
                        input: BodyInput { body, producer },
                        face,
                    })
                })();
            }
            for choice in super::body_choices::choices(document, revision, scene)? {
                if info.consumed_bodies.contains(&choice.input.body) {
                    return Err("A consumed body was returned as an available target".into());
                }
                if let Some(entry) = self
                    .choices
                    .iter_mut()
                    .find(|entry| entry.0 == choice.input.body)
                {
                    entry.1 = choice.label;
                }
                self.native.insert(choice.input.body, choice.input);
            }
            Ok(())
        })();
        self.context_error = result.as_ref().err().cloned();
        result
    }
    pub fn label(&self, id: EntityId) -> String {
        choice_label(&self.choices, id)
    }
    /// True until the author picks an operation; the operation then follows
    /// the extrusion's contact with existing bodies on every preview.
    pub fn automatic(&self) -> bool {
        self.automatic == Some(self.operation)
    }
    /// The automatic default must ask the native model which bodies the
    /// extrusion touches. A document with no bodies is New body at once.
    pub fn needs_discovery(&self) -> bool {
        self.automatic() && !self.names.is_empty()
    }
    /// Pick an operation explicitly. Later automatic decisions never replace it.
    pub fn choose(&mut self, operation: Operation) {
        self.operation = operation;
        self.automatic = None;
        self.decision = None;
        self.invalidate();
    }
    /// Return to the automatic default; the next preview decides again.
    pub fn choose_automatic(&mut self) {
        if self.developed || self.sheet_sketch.is_some() {
            return;
        }
        self.operation = self
            .decision
            .as_ref()
            .map_or(Operation::NewBody, AutomaticOperation::operation);
        self.automatic = Some(self.operation);
        self.invalidate();
    }
    /// Record an automatic decision. Returns false, changing nothing, when the
    /// author has chosen an operation in the meantime.
    pub fn decide(&mut self, decision: AutomaticOperation, contact: CutTargetSelection) -> bool {
        if !self.automatic() {
            return false;
        }
        self.operation = decision.operation();
        self.automatic = Some(self.operation);
        self.decision = Some(decision);
        self.receive(contact);
        true
    }
    #[cfg(test)]
    pub fn decision(&self) -> Option<&AutomaticOperation> {
        self.decision.as_ref().filter(|_| self.automatic())
    }
    fn name(&self, body: EntityId) -> String {
        self.names
            .get(&body)
            .cloned()
            .unwrap_or_else(|| choice_label(&self.choices, body))
    }
    /// Plain inspector feedback for the automatic choice, e.g. "Join · touches Body 1".
    pub fn automatic_summary(&self) -> Option<String> {
        if !self.automatic() {
            return None;
        }
        Some(match &self.decision {
            None => "New body · Preview checks the bodies this extrusion touches".into(),
            Some(AutomaticOperation::NewBody) => "New body · touches no body".into(),
            Some(AutomaticOperation::Join(target)) => {
                format!("Join · touches {}", self.name(target.body))
            }
            Some(AutomaticOperation::Cut(target)) => {
                format!("Cut · enters {}", self.name(target.body))
            }
            Some(AutomaticOperation::Several(several)) => format!(
                "New body · several bodies touch: {}; choose an operation and targets",
                several
                    .iter()
                    .map(|target| self.name(target.body))
                    .collect::<Vec<_>>()
                    .join(", ")
            ),
        })
    }
    /// Join uses automatic discovery until the author picks explicit targets.
    pub fn automatic_join(&self) -> bool {
        self.operation == Operation::Join && self.targets.is_empty()
    }
    pub fn through_sheet(&self) -> bool {
        self.operation == Operation::Cut && self.sheet_sketch.is_some()
    }
    pub fn invalidate(&mut self) {
        self.fresh = false;
    }
    pub fn receive(&mut self, selection: CutTargetSelection) {
        self.selection = Some(selection);
        self.fresh = true;
    }
    pub fn explicit_targets(&self, document: &Document) -> Result<Vec<BodyInput>, String> {
        if let Some(error) = &self.context_error {
            return Err(error.clone());
        }
        self.targets
            .iter()
            .map(|body| {
                let native = self
                    .native
                    .get(body)
                    .ok_or("The selected target is not a current unconsumed native output")?;
                let Some(EntityData::Part { features, .. }) =
                    document.get(*body).map(|entity| &entity.data)
                else {
                    return Err("An extrusion target Part is missing".into());
                };
                let producer = features
                    .iter()
                    .rev()
                    .copied()
                    .find(|id| {
                        matches!(
                            document.get(*id).map(|entity| &entity.data),
                            Some(EntityData::Feature {
                                suppressed: false,
                                ..
                            })
                        )
                    })
                    .ok_or("An extrusion target has no current producer")?;
                if native.producer != producer {
                    return Err("The selected target producer changed; reopen the extrusion".into());
                }
                Ok(native.clone())
            })
            .collect()
    }
    pub fn show(&mut self, ui: &mut egui::Ui) -> bool {
        controls::inspector_layout(ui);
        let mut changed = false;
        ui.label("Operation");
        if self.developed {
            ui.strong("Cut");
            ui.small("The selected developed profile cuts through its sheet.");
            return false;
        }
        let automatic = self.automatic();
        let can_automate = self.sheet_sketch.is_none();
        let mut picked = None;
        controls::inspector_combo(ui, "extrude-body-operation")
            .selected_text(if automatic {
                format!("Automatic · {}", operation_label(self.operation))
            } else {
                operation_label(self.operation).to_owned()
            })
            .show_ui(ui, |ui| {
                if can_automate && ui.selectable_label(automatic, "Automatic").clicked() {
                    picked = Some(None);
                }
                for operation in [
                    Operation::NewBody,
                    Operation::Join,
                    Operation::Cut,
                    Operation::Intersect,
                ] {
                    if ui
                        .selectable_label(
                            !automatic && self.operation == operation,
                            operation_label(operation),
                        )
                        .clicked()
                    {
                        picked = Some(Some(operation));
                    }
                }
            });
        match picked {
            Some(Some(operation)) if automatic || operation != self.operation => {
                self.choose(operation);
                changed = true;
            }
            Some(None) if !automatic => {
                self.choose_automatic();
                changed = true;
            }
            _ => {}
        }
        if let Some(summary) = self.automatic_summary() {
            ui.strong(summary);
            if !self.fresh && self.decision.is_some() {
                ui.small("Previous automatic choice; Preview to update");
            }
            if self.operation == Operation::NewBody {
                return changed;
            }
            changed |= ui
                .checkbox(&mut self.keep_tools, "Keep extrusion tool bodies")
                .on_hover_text(
                    "Keep the generated extrusion Parts alongside the resulting target Parts",
                )
                .changed();
            if changed {
                self.invalidate();
            }
            return changed;
        }
        if self.operation == Operation::NewBody {
            if changed {
                self.invalidate();
            }
            return changed;
        }
        if self.through_sheet() {
            if changed {
                self.invalidate();
            }
            return changed;
        }
        changed |= ui
            .checkbox(&mut self.keep_tools, "Keep extrusion tool bodies")
            .on_hover_text(
                "Keep the generated extrusion Parts alongside the resulting target Parts",
            )
            .changed();
        if self.operation == Operation::Cut {
            ui.small("Cut finds visible intersecting bodies in the chosen component scope. Preview computes exact solid intersections.");
            ui.label("Cut scope").on_hover_text("Starts with the active component when this extrusion opens; explicit includes can reach other bodies or components");
            controls::inspector_combo(ui, "extrude-cut-component")
                .selected_text(
                    self.component
                        .and_then(|id| {
                            self.choices
                                .iter()
                                .find(|entry| entry.0 == id)
                                .map(|entry| entry.1.as_str())
                        })
                        .unwrap_or("Document root"),
                )
                .show_ui(ui, |ui| {
                    changed |= ui
                        .selectable_value(&mut self.component, None, "Document root")
                        .changed();
                    for (id, label, is_body) in &self.choices {
                        if !is_body {
                            changed |= ui
                                .selectable_value(&mut self.component, Some(*id), label)
                                .changed();
                        }
                    }
                });
            if let Some(selection) = &self.selection {
                ui.strong(if self.fresh {
                    "Exact affected bodies"
                } else {
                    "Previous exact set; Preview to update"
                });
                for target in &selection.targets {
                    ui.label(choice_label(&self.choices, target.body));
                }
                if selection.targets.is_empty() {
                    ui.small("No intersecting targets");
                }
                controls::inspector_collapsing(
                    ui,
                    "extrude-cut-excludes",
                    "Exclude affected bodies",
                    |ui| {
                        let bodies: BTreeSet<_> = selection
                            .defaults
                            .iter()
                            .chain(&selection.targets)
                            .map(|target| target.body)
                            .chain(self.exclude.iter().copied())
                            .collect();
                        for body in bodies {
                            let mut excluded = self.exclude.contains(&body);
                            if ui
                                .checkbox(&mut excluded, choice_label(&self.choices, body))
                                .changed()
                            {
                                if excluded {
                                    self.exclude.insert(body);
                                } else {
                                    self.exclude.remove(&body);
                                }
                                changed = true;
                            }
                        }
                    },
                );
                for id in &selection.not_intersecting {
                    ui.small(format!(
                        "No solid overlap: {}",
                        choice_label(&self.choices, *id)
                    ));
                }
            }
            controls::inspector_collapsing(
                ui,
                "extrude-cut-includes",
                "Include other bodies or components",
                |ui| {
                    for (id, label, is_body) in &self.choices {
                        let mut included = self.include.contains(id);
                        if ui
                            .add_enabled(
                                !is_body || self.native.contains_key(id),
                                egui::Checkbox::new(&mut included, label),
                            )
                            .on_disabled_hover_text(
                                "This Part has no current unconsumed native output",
                            )
                            .changed()
                        {
                            if included {
                                self.include.insert(*id);
                            } else {
                                self.include.remove(id);
                            }
                            changed = true;
                        }
                    }
                },
            );
        } else {
            if self.operation == Operation::Join {
                ui.small("With no explicit targets, Join finds the one body the extrusion touches face-to-face or overlaps, including the body the sketch face belongs to.");
                if let Some(selection) = self.selection.as_ref().filter(|_| self.targets.is_empty())
                {
                    ui.strong(match (self.fresh, selection.targets.len()) {
                        (false, _) => "Previous automatic set; Preview to update",
                        (true, 1) => "Automatic target",
                        (true, 0) => "No touching or overlapping bodies",
                        (true, _) => "Several bodies touch the extrusion; choose explicit targets",
                    });
                    for target in &selection.targets {
                        ui.label(choice_label(&self.choices, target.body));
                    }
                }
            }
            ui.strong("Explicit target bodies");
            if let Some(error) = &self.context_error {
                ui.small(error);
            }
            for (id, label, is_body) in &self.choices {
                if !is_body {
                    continue;
                }
                let mut selected = self.targets.contains(id);
                if ui
                    .add_enabled(
                        self.native.contains_key(id),
                        egui::Checkbox::new(&mut selected, label),
                    )
                    .on_disabled_hover_text("This Part has no current unconsumed native output")
                    .changed()
                {
                    if selected {
                        self.targets.insert(*id);
                    } else {
                        self.targets.remove(id);
                    }
                    changed = true;
                }
            }
        }
        if changed {
            self.invalidate();
        }
        changed
    }
}
fn choice_label(choices: &[(EntityId, String, bool)], id: EntityId) -> String {
    choices
        .iter()
        .find(|entry| entry.0 == id)
        .map(|entry| entry.1.clone())
        .unwrap_or_else(|| id.to_string())
}
fn operation_label(operation: Operation) -> &'static str {
    match operation {
        Operation::NewBody => "New body",
        Operation::Join => "Join / merge",
        Operation::Cut => "Cut / remove material",
        Operation::Intersect => "Intersect",
        Operation::Split => "Split",
    }
}

/// What one discovery asks the native model.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum DiscoveryMode {
    /// Bodies the tools touch face-to-face or overlap.
    Join,
    /// Bodies the tools enter.
    Cut,
    /// Both, in that probe order, to choose the default operation (CAD-404).
    Automatic,
}
impl DiscoveryMode {
    fn contacts(self) -> &'static [bool] {
        match self {
            Self::Join => &[true],
            Self::Cut => &[false],
            Self::Automatic => &[true, false],
        }
    }
}
pub(super) struct CutDiscovery {
    pub prepared: Prepared,
    pub options: EvaluationOptions,
    pub content: ContentHash,
    extrusion: ExtrudeRequest,
    mode: DiscoveryMode,
}
impl CutDiscovery {
    pub fn prepare(
        session: &Session,
        actor: Actor,
        extrusion: ExtrudeRequest,
        component: Option<EntityId>,
        include: Vec<EntityId>,
        exclude: Vec<EntityId>,
        mode: DiscoveryMode,
    ) -> Result<Self, String> {
        let tools = model::prepare_extrude(session, session.revision(), actor, &extrusion)
            .map_err(|error| error.to_string())?;
        let content = ContentHash::document(tools.prepared.document())?;
        let tool_ids: Vec<_> = extrusion
            .outputs
            .iter()
            .map(|output| output.feature)
            .collect();
        let options = EvaluationOptions {
            cut_targets: mode
                .contacts()
                .iter()
                .map(|&contact| CutTargetProbe {
                    id: EntityId::new(),
                    query: CutTargetQuery {
                        tools: tool_ids.clone(),
                        component,
                        include: include.clone(),
                        exclude: exclude.clone(),
                        contact,
                    },
                })
                .collect(),
            ..Default::default()
        };
        options.validate().map_err(|error| error.to_string())?;
        Ok(Self {
            prepared: tools.prepared,
            options,
            content,
            extrusion,
            mode,
        })
    }
    /// One selection per probe, in [`DiscoveryMode`] probe order.
    pub fn selections<'a>(
        &self,
        session: &Session,
        scene: &'a Scene,
        expected: SceneStamp,
    ) -> Result<Vec<&'a CutTargetSelection>, String> {
        session
            .validate_preview(&self.prepared)
            .map_err(|error| error.to_string())?;
        if scene.stamp() != expected
            || expected.document != session.document().id()
            || expected.revision != session.revision()
        {
            return Err("Cut target results belong to an obsolete modeling request".into());
        }
        let info = scene
            .model_evaluation()
            .ok_or("The native provider did not return exact cut targets")?;
        if info.source_hash != self.content {
            return Err("Cut target results describe different extrusion source".into());
        }
        if info.cut_targets.len() != self.options.cut_targets.len()
            || info
                .cut_targets
                .iter()
                .zip(&self.options.cut_targets)
                .any(|(result, probe)| result.id != probe.id)
        {
            return Err("Cut target results do not match the requested probe".into());
        }
        for output in &self.extrusion.outputs {
            if !scene
                .reference_instances()
                .iter()
                .any(|instance| instance.id == output.body && instance.geometry == output.feature)
                || !scene
                    .geometry_snapshot(output.feature)
                    .is_some_and(|snapshot| !snapshot.elements.is_empty())
            {
                return Err(
                    "Native target discovery did not retain the exact extrusion tools".into(),
                );
            }
        }
        let selections = info
            .cut_targets
            .iter()
            .map(|probe| probe.result.as_ref().map_err(|error| error.to_string()))
            .collect::<Result<Vec<_>, _>>()?;
        for target in selections
            .iter()
            .flat_map(|selection| selection.defaults.iter().chain(&selection.targets))
        {
            target
                .resolve_producer(session.document())
                .map_err(|error| error.to_string())?;
            let head = match &session
                .document()
                .get(target.body)
                .ok_or("Cut target body is missing")?
                .data
            {
                EntityData::Part { features, .. } => features.iter().rev().copied().find(|id| {
                    matches!(
                        session.document().get(*id).map(|entity| &entity.data),
                        Some(EntityData::Feature {
                            suppressed: false,
                            ..
                        })
                    )
                }),
                _ => None,
            };
            if head != Some(target.producer) {
                return Err(
                    "Cut target results reference a historical rather than current body producer"
                        .into(),
                );
            }
            if !scene
                .reference_instances()
                .iter()
                .any(|instance| instance.id == target.body && instance.geometry == target.producer)
            {
                return Err("Cut target results do not match current native body producers".into());
            }
        }
        Ok(selections)
    }
    pub fn extrusion(&self) -> &ExtrudeRequest {
        &self.extrusion
    }
    pub fn mode(&self) -> DiscoveryMode {
        self.mode
    }
}

/// Stable final producers across parameter edits, target overrides and repeated
/// previews in one form. Removing then restoring a target reuses its identity.
#[derive(Default)]
pub(super) struct BooleanIdentities {
    features: BTreeMap<EntityId, EntityId>,
    group: Option<model::ExtrudeGroupDraft>,
}
impl BooleanIdentities {
    pub fn new(document: &Document) -> Self {
        let name = super::names::GeneratedNames::new(document).numbered("Extrusion operation", "");
        Self {
            features: BTreeMap::new(),
            group: Some(model::ExtrudeGroupDraft::new(name, None)),
        }
    }
    pub fn group_draft(&mut self) -> model::ExtrudeGroupDraft {
        self.group
            .get_or_insert_with(|| model::ExtrudeGroupDraft::new("Extrusion operation", None))
            .clone()
    }
    pub fn compose(
        &mut self,
        extrusion: ExtrudeRequest,
        operation: Operation,
        targets: Vec<BodyInput>,
        keep_tools: bool,
    ) -> Result<ExtrudeBooleanRequest, String> {
        if !matches!(
            operation,
            Operation::Join | Operation::Cut | Operation::Intersect
        ) {
            return Err("Choose Join, Cut or Intersect for an extrusion body operation".into());
        }
        if targets.is_empty() {
            return Err("No affected bodies; change the extent or choose explicit targets".into());
        }
        let mut request = ExtrudeBooleanRequest::new(extrusion, operation, targets, keep_tools);
        request.group = Some(self.group_draft());
        for boolean in &mut request.booleans {
            boolean.feature = *self.features.entry(boolean.target.body).or_default();
        }
        Ok(request)
    }
}
