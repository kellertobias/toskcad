use super::{
    FormStamp, PreviewState,
    body_choices::{self, BodyChoice},
    form::{self, FormAction, ModelEvent, ModelRequest},
    names::GeneratedNames,
};
use crate::{controls, egui, viewport::Scene};
use multicad_core::{
    document::{Document, EntityId},
    model::{BooleanRequest, MAX_BOOLEAN_TARGETS, MAX_BOOLEAN_TOOLS, Operation},
};
use std::collections::{BTreeMap, BTreeSet};

pub(super) struct BooleanPanel {
    stamp: FormStamp,
    choices: Vec<BodyChoice>,
    features: BTreeMap<EntityId, (EntityId, String)>,
    targets: BTreeSet<EntityId>,
    tools: BTreeSet<EntityId>,
    operation: Operation,
    keep_tools: bool,
    filter: String,
    error: Option<String>,
    /// Bumped on each list edit so the viewer selection can follow it.
    generation: u64,
    /// Join / Intersect: chosen bodies in selection order. The first keeps
    /// its name and identity; the user can move bodies up (CAD-408).
    order: Vec<EntityId>,
}
impl BooleanPanel {
    pub fn available(document: &Document, revision: u64, scene: &Scene) -> Result<usize, String> {
        body_choices::available(document, revision, scene)
    }
    pub fn new(
        document: &Document,
        revision: u64,
        scene: &Scene,
        selected: &BTreeSet<EntityId>,
    ) -> Result<Self, String> {
        let choices = body_choices::choices(document, revision, scene)?;
        if choices.len() < 2 {
            return Err("Boolean modeling needs at least two available native bodies".into());
        }
        let mut names = GeneratedNames::new(document);
        let features = choices
            .iter()
            .map(|choice| {
                (
                    choice.input.body,
                    (EntityId::new(), names.numbered("Boolean feature", "")),
                )
            })
            .collect();
        let order: Vec<_> = choices
            .iter()
            .map(|choice| choice.input.body)
            .filter(|id| selected.contains(id))
            .collect();
        let targets = order.iter().copied().collect();
        Ok(Self {
            stamp: FormStamp::new(),
            choices,
            features,
            targets,
            tools: BTreeSet::new(),
            operation: Operation::Cut,
            keep_tools: false,
            filter: String::new(),
            error: None,
            generation: 0,
            order,
        })
    }
    /// Preset the operation chosen from the ribbon's Boolean menu. Only
    /// implemented Join, Subtract and Intersect operations are accepted.
    pub fn with_operation(mut self, operation: Operation) -> Self {
        if matches!(
            operation,
            Operation::Join | Operation::Cut | Operation::Intersect
        ) {
            self.operation = operation;
        }
        self
    }
    pub fn operation(&self) -> Operation {
        self.operation
    }
    /// Join and Intersect treat every chosen body alike (CAD-408).
    fn symmetric(&self) -> bool {
        self.operation != Operation::Cut
    }
    /// Every chosen body, in either role.
    pub fn selected_bodies(&self) -> BTreeSet<EntityId> {
        self.targets.union(&self.tools).copied().collect()
    }
    pub fn generation(&self) -> u64 {
        self.generation
    }
    /// Chosen bodies in their Join / Intersect order: the remembered order,
    /// then any body chosen another way, in list order.
    pub fn ordered(&self) -> Vec<EntityId> {
        let chosen = self.selected_bodies();
        let mut ordered: Vec<_> = self
            .order
            .iter()
            .copied()
            .filter(|id| chosen.contains(id))
            .collect();
        for choice in &self.choices {
            let id = choice.input.body;
            if chosen.contains(&id) && !ordered.contains(&id) {
                ordered.push(id);
            }
        }
        ordered
    }
    /// Remember newly chosen bodies after the ones already chosen.
    fn refresh_order(&mut self) {
        self.order = self.ordered();
    }
    /// Bodies picked before the form opened keep their pick order.
    pub fn prefer_order(&mut self, picked: &[EntityId]) {
        let chosen = self.selected_bodies();
        let mut order: Vec<_> = picked
            .iter()
            .copied()
            .filter(|id| chosen.contains(id))
            .collect();
        order.dedup();
        self.order = order;
        self.refresh_order();
        self.stamp.changed();
    }
    /// Follow the viewer / browser selection. Join and Intersect take the
    /// bodies as they are; Cut keeps the role of bodies still selected and
    /// gives new ones the role still missing.
    pub fn sync_selection(&mut self, selection: &BTreeSet<EntityId>) -> bool {
        let chosen: BTreeSet<_> = self
            .choices
            .iter()
            .map(|choice| choice.input.body)
            .filter(|id| selection.contains(id))
            .collect();
        if chosen == self.selected_bodies() {
            return false;
        }
        self.targets.retain(|id| chosen.contains(id));
        self.tools.retain(|id| chosen.contains(id));
        for id in &chosen {
            if self.targets.contains(id) || self.tools.contains(id) {
                continue;
            }
            if self.symmetric() || self.targets.is_empty() {
                self.targets.insert(*id);
            } else {
                self.tools.insert(*id);
            }
        }
        self.refresh_order();
        self.stamp.changed();
        self.error = None;
        true
    }
    fn verb(&self) -> &'static str {
        match self.operation {
            Operation::Intersect => "intersect",
            _ => "join",
        }
    }
    /// What to do next while the selection is incomplete; `None` once it
    /// can preview. Internal limits appear only when exceeded.
    pub fn hint(&self) -> Option<String> {
        if self.symmetric() {
            let count = self.selected_bodies().len();
            match count {
                0 => Some(format!("Select two or more bodies to {}", self.verb())),
                1 => Some(format!("Select at least one more body to {}", self.verb())),
                n if n > MAX_BOOLEAN_TOOLS + 1 => Some(format!(
                    "Select at most {} bodies to {}",
                    MAX_BOOLEAN_TOOLS + 1,
                    self.verb()
                )),
                _ => None,
            }
        } else if self.targets.is_empty() {
            Some("Select the body to cut".into())
        } else if self.tools.is_empty() {
            Some("Select one or more bodies to cut away".into())
        } else if self.targets.len() > MAX_BOOLEAN_TARGETS {
            Some(format!(
                "Select at most {MAX_BOOLEAN_TARGETS} bodies to cut"
            ))
        } else if self.tools.len() > MAX_BOOLEAN_TOOLS {
            Some(format!(
                "Select at most {MAX_BOOLEAN_TOOLS} bodies to cut away"
            ))
        } else {
            None
        }
    }
    pub fn stamp(&self) -> FormStamp {
        self.stamp
    }
    pub fn request(
        &self,
        document: &Document,
        enabled: bool,
    ) -> Result<Vec<BooleanRequest>, String> {
        if !enabled {
            return Err("The native modeling preview provider is unavailable".into());
        }
        if !matches!(
            self.operation,
            Operation::Join | Operation::Cut | Operation::Intersect
        ) {
            return Err(
                "Choose Join, Subtract or Intersect; Split outputs are not implemented".into(),
            );
        }
        if !self.symmetric() && !self.targets.is_disjoint(&self.tools) {
            return Err("A body to cut cannot also be cut away".into());
        }
        if let Some(hint) = self.hint() {
            return Err(hint);
        }
        let input = |id| {
            let input = self
                .choices
                .iter()
                .find(|choice| choice.input.body == id)
                .ok_or("Body selection is unavailable; restart Boolean")?
                .input
                .clone();
            input
                .resolve_producer(document)
                .map_err(|error| error.to_string())?;
            Ok::<_, String>(input)
        };
        if self.symmetric() {
            // One result: the first body in list order keeps its identity
            // and every other chosen body is merged into it.
            let ordered = self.ordered();
            let mut bodies = ordered.iter().copied();
            let first = bodies.next().ok_or("Select two or more bodies")?;
            let (feature, name) = self
                .features
                .get(&first)
                .ok_or("Body identity was not allocated; restart Boolean")?;
            return Ok(vec![BooleanRequest {
                feature: *feature,
                name: name.clone(),
                target: input(first)?,
                tools: bodies.map(input).collect::<Result<Vec<_>, _>>()?,
                operation: self.operation,
                keep_tools: false,
            }]);
        }
        let tools = self
            .tools
            .iter()
            .map(|id| input(*id))
            .collect::<Result<Vec<_>, _>>()?;
        self.targets
            .iter()
            .map(|id| {
                let (feature, name) = self
                    .features
                    .get(id)
                    .ok_or("Target identity was not allocated; restart Boolean")?;
                Ok(BooleanRequest {
                    feature: *feature,
                    name: name.clone(),
                    target: input(*id)?,
                    tools: tools.clone(),
                    operation: self.operation,
                    keep_tools: self.keep_tools,
                })
            })
            .collect()
    }
    pub fn show(
        &mut self,
        ui: &mut egui::Ui,
        document: &Document,
        enabled: bool,
        preview: PreviewState<'_>,
    ) -> Option<ModelEvent> {
        controls::inspector_layout(ui);
        ui.strong(operation_label(self.operation));
        ui.small(match self.operation {
            Operation::Intersect => "Keeps only the region shared by the chosen bodies. Choose them in the viewer, the browser or the list below; the first in the order keeps its name.",
            Operation::Cut => "Removes the bodies to cut away from the body to cut. Choose them in the viewer, the browser or the lists below.",
            _ => "Merges the chosen bodies into one. Choose them in the viewer, the browser or the list below; the first in the order keeps its name.",
        });
        ui.label("Operation");
        let mut changed = false;
        let before = self.operation;
        controls::inspector_combo(ui, "boolean-operation")
            .selected_text(operation_label(self.operation))
            .show_ui(ui, |ui| {
                for operation in [Operation::Join, Operation::Cut, Operation::Intersect] {
                    changed |= ui
                        .selectable_value(
                            &mut self.operation,
                            operation,
                            operation_label(operation),
                        )
                        .changed();
                }
                ui.add_enabled(false, egui::Button::new("Split"))
                    .on_disabled_hover_text("Split is not available in this build.");
            });
        if self.operation != before && self.symmetric() {
            // Join and Intersect have one list; keep every chosen body.
            let tools = std::mem::take(&mut self.tools);
            self.targets.extend(tools);
        }
        ui.add(
            egui::TextEdit::singleline(&mut self.filter)
                .desired_width(ui.available_width())
                .hint_text("Filter bodies"),
        )
        .on_hover_text("Filter the lists without changing the chosen bodies");
        ui.separator();
        let mut edited = false;
        if self.symmetric() {
            ui.strong(format!(
                "Bodies to {} · {}",
                self.verb(),
                self.targets.len()
            ));
            let mut none = BTreeSet::new();
            edited |= Self::body_list(
                ui,
                "boolean-targets",
                &self.choices,
                &self.filter,
                &mut self.targets,
                &mut none,
            );
            if edited {
                self.refresh_order();
            }
            changed |= self.show_order(ui);
        } else {
            ui.strong(format!("Body to cut · {}", self.targets.len()));
            edited |= Self::body_list(
                ui,
                "boolean-targets",
                &self.choices,
                &self.filter,
                &mut self.targets,
                &mut self.tools,
            );
            ui.strong(format!("Bodies to cut away · {}", self.tools.len()));
            edited |= Self::body_list(
                ui,
                "boolean-tools",
                &self.choices,
                &self.filter,
                &mut self.tools,
                &mut self.targets,
            );
            changed |= ui
                .checkbox(&mut self.keep_tools, "Keep the bodies cut away")
                .on_hover_text("Keep each body cut away as its own body. Otherwise it is consumed by the cut; its history remains in the document.")
                .changed();
        }
        if edited {
            self.generation += 1;
        }
        changed |= edited;
        let validation = self.request(document, enabled);
        // One message at a time: the next step while the selection is
        // incomplete, otherwise a real error or the preview status.
        if let Some(hint) = self.hint() {
            ui.label(hint);
        } else if let Err(error) = &validation {
            ui.colored_label(egui::Color32::LIGHT_RED, error);
        } else if let Some(error) = &self.error {
            ui.colored_label(egui::Color32::LIGHT_RED, error);
        } else {
            form::preview_status(ui, self.stamp, preview);
        }
        let action = form::actions(ui, validation.is_ok(), changed, self.stamp, preview);
        self.finish_frame(changed, action, document, enabled, preview)
    }
    /// The Join / Intersect order: the first body keeps its name and
    /// identity, and any other can be moved up.
    fn show_order(&mut self, ui: &mut egui::Ui) -> bool {
        let ordered = self.ordered();
        if ordered.len() < 2 {
            return false;
        }
        ui.small("Order");
        let mut raise = None;
        for (index, id) in ordered.iter().enumerate() {
            let label = self
                .choices
                .iter()
                .find(|choice| choice.input.body == *id)
                .map_or("Missing body", |choice| choice.label.as_str());
            ui.push_id(("boolean-order", *id), |ui| {
                ui.horizontal_wrapped(|ui| {
                    ui.label(format!("{}. {label}", index + 1));
                    if index == 0 {
                        ui.small("keeps its name");
                    } else if ui
                        .small_button("Move up")
                        .on_hover_text("Move this body up; the first body keeps its name")
                        .clicked()
                    {
                        raise = Some(index);
                    }
                });
            });
        }
        let Some(index) = raise else {
            return false;
        };
        let mut order = ordered;
        order.swap(index - 1, index);
        self.order = order;
        true
    }
    fn body_list(
        ui: &mut egui::Ui,
        salt: &str,
        choices: &[BodyChoice],
        filter: &str,
        selection: &mut BTreeSet<EntityId>,
        opposite: &mut BTreeSet<EntityId>,
    ) -> bool {
        let mut changed = false;
        let filter = filter.to_lowercase();
        egui::ScrollArea::vertical()
            .id_salt(salt)
            .max_height(150.)
            .show(ui, |ui| {
                for choice in choices.iter().filter(|choice| {
                    filter.is_empty() || choice.label.to_lowercase().contains(&filter)
                }) {
                    let mut checked = selection.contains(&choice.input.body);
                    if ui
                        .push_id(choice.input.body, |ui| {
                            ui.checkbox(&mut checked, &choice.label)
                        })
                        .inner
                        .on_hover_text(format!("{}\nProduced by {}", choice.path, choice.producer))
                        .changed()
                    {
                        if checked {
                            selection.insert(choice.input.body);
                            opposite.remove(&choice.input.body);
                        } else {
                            selection.remove(&choice.input.body);
                        }
                        changed = true;
                    }
                }
            });
        changed
    }
    fn finish_frame(
        &mut self,
        changed: bool,
        action: Option<FormAction>,
        document: &Document,
        enabled: bool,
        preview: PreviewState<'_>,
    ) -> Option<ModelEvent> {
        if action == Some(FormAction::Cancel) {
            return Some(ModelEvent::Cancel);
        }
        if changed {
            self.stamp.changed();
            self.error = None;
            if action != Some(FormAction::Preview) {
                return Some(ModelEvent::Changed);
            }
        }
        match action {
            Some(FormAction::Preview) => match self.request(document, enabled) {
                Ok(requests) => {
                    self.error = None;
                    Some(ModelEvent::Preview {
                        stamp: self.stamp,
                        request: ModelRequest::Booleans(requests),
                    })
                }
                Err(error) => {
                    self.error = Some(error);
                    None
                }
            },
            Some(FormAction::Apply)
                if preview.ready(self.stamp) && self.request(document, enabled).is_ok() =>
            {
                Some(ModelEvent::Apply)
            }
            _ => None,
        }
    }
}
fn operation_label(operation: Operation) -> &'static str {
    match operation {
        Operation::Join => "Join",
        Operation::Cut => "Subtract",
        Operation::Intersect => "Intersect",
        Operation::Split => "Split",
        Operation::NewBody => "New body",
    }
}

#[cfg(test)]
mod tests;
