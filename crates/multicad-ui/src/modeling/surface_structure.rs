//! Editable surface relief settings. File selection only changes this draft.
mod linked;
mod visual;
#[cfg(test)]
mod visual_tests;

use super::{
    FormStamp, PreviewState,
    diagram::{self, Choice, DiagramLink, Key},
    dimension::DimensionField,
    form::{self, FormAction, ModelEvent, ModelRequest},
};
use crate::{controls, egui, viewport::Scene};
use multicad_core::{
    document::{Document, EntityId},
    model::{
        BodyInput, MAX_SURFACE_IMAGE_REPEAT, MAX_SURFACE_STRUCTURE_FACES,
        SurfaceStructureDirection as Direction, SurfaceStructureParameters as Parameters,
        SurfaceStructurePattern as Pattern, SurfaceStructureRequest,
    },
    references::{self, AnalyticGeometry, GeometryIntent, Resolution, TopologyKind},
};
use std::collections::BTreeSet;

pub(super) struct SurfaceStructurePanel {
    stamp: FormStamp,
    feature: EntityId,
    input: Option<BodyInput>,
    faces: Vec<GeometryIntent>,
    pattern: Pattern,
    direction: Direction,
    amplitude: DimensionField,
    pitch: DimensionField,
    size: [DimensionField; 2],
    offset: [DimensionField; 2],
    rotation: DimensionField,
    border: DimensionField,
    picking: bool,
    picking_available: bool,
    image_requested: bool,
    /// Callout ↔ field linking of the relief schematic.
    link: DiagramLink,
}

impl SurfaceStructurePanel {
    pub fn new() -> Result<Self, String> {
        Self::from_source(EntityId::new(), None, vec![], Parameters::default())
    }
    pub fn existing(request: SurfaceStructureRequest) -> Result<Self, String> {
        let mut panel = Self::from_source(
            request.feature,
            Some(request.input),
            request.faces,
            request.parameters,
        )?;
        panel.picking = false;
        Ok(panel)
    }
    fn from_source(
        feature: EntityId,
        input: Option<BodyInput>,
        faces: Vec<GeometryIntent>,
        parameters: Parameters,
    ) -> Result<Self, String> {
        let mm = |value| DimensionField::from_dimension(value, "mm");
        Ok(Self {
            stamp: FormStamp::new(),
            feature,
            input,
            faces,
            pattern: parameters.pattern,
            direction: parameters.direction,
            amplitude: mm(&parameters.amplitude)?,
            pitch: mm(&parameters.pitch)?,
            size: [mm(&parameters.size[0])?, mm(&parameters.size[1])?],
            offset: [mm(&parameters.offset[0])?, mm(&parameters.offset[1])?],
            rotation: DimensionField::from_dimension(&parameters.rotation, "deg")?,
            border: mm(&parameters.border)?,
            picking: true,
            picking_available: false,
            image_requested: false,
            link: DiagramLink::default(),
        })
    }
    pub fn stamp(&self) -> FormStamp {
        self.stamp
    }
    pub fn input(&self) -> Option<&BodyInput> {
        self.input.as_ref()
    }
    pub fn faces(&self) -> &[GeometryIntent] {
        &self.faces
    }
    pub fn picking(&self) -> bool {
        self.picking && self.picking_available
    }
    pub fn set_picking_available(&mut self, available: bool) {
        self.picking_available = available;
    }
    pub fn stop_picking(&mut self) {
        self.picking = false;
    }
    pub fn take_image_request(&mut self) -> Option<FormStamp> {
        std::mem::take(&mut self.image_requested).then_some(self.stamp)
    }
    pub fn set_image(&mut self, pattern: Pattern) -> Result<(), String> {
        if !matches!(pattern, Pattern::HeightMap { .. }) {
            return Err("The selected file did not supply a grayscale height map".into());
        }
        self.pattern = pattern;
        self.stamp.changed();
        Ok(())
    }
    pub fn set_faces(
        &mut self,
        scene: &Scene,
        input: BodyInput,
        faces: Vec<GeometryIntent>,
    ) -> Result<bool, String> {
        if self.input.as_ref().is_some_and(|current| *current != input) {
            return Err("Choose faces on the same source Part".into());
        }
        if faces.len() > MAX_SURFACE_STRUCTURE_FACES {
            return Err(format!(
                "Choose up to {MAX_SURFACE_STRUCTURE_FACES} planar faces on one Part"
            ));
        }
        let snapshot = scene
            .geometry_snapshot(input.producer)
            .ok_or("The retained source Part is unavailable")?;
        let mut keys = BTreeSet::new();
        for face in &faces {
            match references::resolve(face, &snapshot, scene.stamp().revision) {
                Resolution::Resolved { key, .. }
                    if key.kind == TopologyKind::Face
                        && keys.insert(key)
                        && snapshot.elements.iter().any(|item| {
                            item.key == key
                                && item.frames.len() == 1
                                && matches!(item.geometry, Some(AnalyticGeometry::Plane { .. }))
                        }) => {}
                _ => return Err("Choose distinct named planar faces on the retained Part".into()),
            }
        }
        if self.input.as_ref() == Some(&input) && self.faces == faces {
            return Ok(false);
        }
        self.input = Some(input);
        self.faces = faces;
        self.stamp.changed();
        Ok(true)
    }

    pub fn request(
        &self,
        document: &Document,
        enabled: bool,
    ) -> Result<SurfaceStructureRequest, String> {
        if !enabled {
            return Err("Native modeling Preview is unavailable".into());
        }
        let input = self.input.clone().ok_or("Pick a planar face on a Part")?;
        if self.faces.is_empty() {
            return Err("Pick at least one planar face to generate the preview".into());
        }
        let pitch = match self.pitch.resolve(document) {
            Ok(pitch) => pitch,
            Err(_) if matches!(self.pattern, Pattern::HeightMap { .. }) => {
                Parameters::default().pitch
            }
            Err(error) => return Err(error),
        };
        Ok(SurfaceStructureRequest {
            feature: self.feature,
            input,
            faces: self.faces.clone(),
            parameters: Parameters {
                pattern: self.pattern.clone(),
                direction: self.direction,
                amplitude: self.amplitude.resolve(document)?,
                pitch,
                size: [
                    self.size[0].resolve(document)?,
                    self.size[1].resolve(document)?,
                ],
                offset: [
                    self.offset[0].resolve(document)?,
                    self.offset[1].resolve(document)?,
                ],
                rotation: self.rotation.resolve(document)?,
                border: self.border.resolve(document)?,
            },
        })
    }

    pub fn show(
        &mut self,
        ui: &mut egui::Ui,
        document: &Document,
        enabled: bool,
        preview: PreviewState<'_>,
    ) -> Option<ModelEvent> {
        controls::inspector_layout(ui);
        ui.strong("Surface Structure");
        form::preview_status(ui, self.stamp, preview);
        if let Some(input) = &self.input {
            ui.label(
                document
                    .get(input.body)
                    .map_or("Missing Part", |part| part.name.as_str()),
            );
        }
        ui.small(format!(
            "{} of {MAX_SURFACE_STRUCTURE_FACES} planar faces selected",
            self.faces.len()
        ));
        let pick = if self.picking() {
            "Finish picking faces"
        } else {
            "Pick or repair faces"
        };
        if ui
            .add_enabled(self.picking_available, egui::Button::new(pick))
            .clicked()
        {
            self.picking = !self.picking;
        }
        if self.picking() {
            ui.small("Click a planar face. Hold Shift to add or remove a face.");
        }
        ui.separator();
        let mut changed = false;
        let mut link = std::mem::take(&mut self.link);
        link.invalid = diagram::first_invalid(document, self.linked_fields());
        let schematic = visual::diagram(self.kind(), &self.diagram_values(document), self.image());
        let values: Vec<(Key, String)> = schematic
            .keys()
            .into_iter()
            .map(|k| (k, self.callout_value(k)))
            .collect();
        let value = |key: Key| {
            values
                .iter()
                .find(|(k, _)| *k == key)
                .map_or_else(String::new, |(_, v)| v.clone())
        };
        link.diagram(
            ui,
            ("surface-structure", self.feature),
            &schematic,
            240.,
            &value,
        );
        ui.small("One face with its relief window and border · click a callout to edit it");
        ui.label("Pattern");
        let mut kind = self.kind();
        let choices: Vec<Choice<'_, visual::Kind>> = [
            (
                visual::Kind::Ribs,
                "Ribs",
                "Parallel raised or recessed bars",
            ),
            (visual::Kind::Waves, "Waves", "Smooth periodic waves"),
            (visual::Kind::Stipple, "Stipple", "A grid of rounded dots"),
        ]
        .into_iter()
        .map(|(value, label, hint)| Choice {
            value,
            label,
            hint,
            art: visual::tile_art(value),
        })
        .collect();
        if diagram::choice_grid(ui, "surface-structure-pattern", &mut kind, &choices, true) {
            self.pattern = match kind {
                visual::Kind::Ribs => Pattern::Ribs,
                visual::Kind::Waves => Pattern::Waves,
                _ => Pattern::Stipple,
            };
            changed = true;
        }
        if matches!(self.pattern, Pattern::HeightMap { .. }) {
            ui.small("Height map imported");
        }
        if ui.button("Import grayscale PNG…").clicked() {
            self.image_requested = true;
        }
        if let Pattern::HeightMap {
            width,
            height,
            repeat,
            ..
        } = &mut self.pattern
        {
            ui.small(format!("Embedded {width} × {height} pixels · black is zero relief"))
                .on_hover_text("The top image row maps to +V. Repeated tiles fit the relief window and each raised or recessed cell leaves a 10% gap.");
            for (label, value) in ["Repeat U", "Repeat V"].into_iter().zip(repeat) {
                changed |= ui
                    .horizontal(|ui| {
                        ui.label(label);
                        ui.add(egui::DragValue::new(value).range(1..=MAX_SURFACE_IMAGE_REPEAT))
                            .changed()
                    })
                    .inner;
            }
        } else {
            changed |= link.field(
                ui,
                visual::PITCH,
                &mut self.pitch,
                "surface-pitch",
                "Pattern pitch",
                "Physical spacing between repeats.",
            );
        }
        ui.horizontal(|ui| {
            changed |= ui
                .selectable_value(&mut self.direction, Direction::Raised, "Raised")
                .changed();
            changed |= ui
                .selectable_value(&mut self.direction, Direction::Recessed, "Recessed")
                .changed();
        });
        changed |= link.field(
            ui,
            visual::AMPLITUDE,
            &mut self.amplitude,
            "surface-amplitude",
            "Relief height or depth",
            "Maximum physical relief above or below each face.",
        );
        changed |= link.field(
            ui,
            visual::SIZE[0],
            &mut self.size[0],
            "surface-width",
            "Width U",
            "Width of the rotated relief window.",
        );
        changed |= link.field(
            ui,
            visual::SIZE[1],
            &mut self.size[1],
            "surface-height",
            "Height V",
            "Height of the rotated relief window.",
        );
        let title = format!(
            "Window placement · offset {} · {} · border {}",
            self.callout_value(visual::OFFSET[0]),
            diagram::field_value(&self.rotation),
            diagram::field_value(&self.border)
        );
        let reveal = link.reveals(&visual::PLACEMENT);
        diagram::section(ui, "surface-placement-section", &title, reveal, |ui| {
            changed |= link.field(
                ui,
                visual::OFFSET[0],
                &mut self.offset[0],
                "surface-offset-u",
                "Offset U",
                "Shift the window from the face extent midpoint.",
            );
            changed |= link.field(
                ui,
                visual::OFFSET[1],
                &mut self.offset[1],
                "surface-offset-v",
                "Offset V",
                "Shift the window from the face extent midpoint.",
            );
            changed |= link.field(
                ui,
                visual::ROTATION,
                &mut self.rotation,
                "surface-rotation",
                "Rotation",
                "Rotate the window and pattern together about the outward normal.",
            );
            changed |= link.field(
                ui,
                visual::BORDER,
                &mut self.border,
                "surface-border",
                "Face border",
                "Keep this much of each face flat at its boundary.",
            );
        });
        self.link = link;
        ui.small(
            "Each face has its own mapping and flat border. Patterns stop before shared edges.",
        );
        if changed {
            self.stamp.changed();
        }
        let request = self.request(document, enabled);
        if let Err(error) = &request {
            ui.small(error);
        }
        let pending = matches!(preview, PreviewState::Pending { stamp } if stamp == self.stamp);
        match form::actions(
            ui,
            request.is_ok() && !pending,
            changed,
            self.stamp,
            preview,
        ) {
            Some(FormAction::Preview) => request.ok().map(|request| ModelEvent::Preview {
                stamp: self.stamp,
                request: ModelRequest::SurfaceStructure(Box::new(request)),
            }),
            Some(FormAction::Apply) => Some(ModelEvent::Apply),
            Some(FormAction::Cancel) => Some(ModelEvent::Cancel),
            None => changed.then_some(ModelEvent::Changed),
        }
    }
}
