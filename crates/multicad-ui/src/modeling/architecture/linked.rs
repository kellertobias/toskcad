//! Linked slab, pedestal and stair Inspector sections (CAD-395): each shows
//! its schematic first, whose callouts focus the same staged fields listed
//! below it; secondary inputs sit in named sections that validation opens.
use super::{ArchitectureKind, ArchitecturePanel, RoomBoundary, rooms, visual};
use crate::{
    egui,
    modeling::diagram::{self, Choice, Key},
};
use multicad_core::{
    document::{Document, EntityId},
    model, parameters,
};

impl ArchitecturePanel {
    /// The slab, pedestal or stair fields linked to callouts, in display order.
    pub(super) fn linked_fields(&self) -> Vec<(Key, &super::DimensionField)> {
        match self.kind {
            ArchitectureKind::Pedestal => vec![(visual::PEDESTAL_HEIGHT, &self.thickness)],
            ArchitectureKind::Slab if self.override_thickness => {
                vec![(visual::SLAB_THICKNESS, &self.thickness)]
            }
            ArchitectureKind::Stair => vec![
                (visual::STAIR_WIDTH, &self.width),
                (visual::STAIR_GOING, &self.going),
            ],
            _ => vec![],
        }
    }
    fn level_slab_mm(&self, document: &Document) -> Option<f64> {
        let parameters = parameters::evaluate(document).ok()?;
        model::resolve_level(document, &parameters, self.level?)
            .ok()
            .map(|l| l.slab_mm)
    }
    pub(super) fn slab_values(&self, document: &Document) -> visual::SlabValues {
        let evaluated = parameters::evaluate(document).ok();
        visual::SlabValues {
            outline: self.outline.positions(document).ok(),
            pedestal: self.kind == ArchitectureKind::Pedestal,
            override_thickness: self.override_thickness,
            thickness: evaluated
                .as_ref()
                .and_then(|p| diagram::value_of(&self.thickness, document, p)),
            level_thickness: self.level_slab_mm(document),
        }
    }
    pub(super) fn stair_values(&self, document: &Document) -> visual::StairValues {
        let evaluated = parameters::evaluate(document).ok();
        let value = |f| {
            evaluated
                .as_ref()
                .and_then(|p| diagram::value_of(f, document, p))
        };
        let riser = self.stair(document).ok().and_then(|stair| {
            model::resolve_stair(document, evaluated.as_ref()?, &stair)
                .ok()
                .map(|r| r.check.riser_mm)
        });
        visual::StairValues {
            width: value(&self.width),
            going: value(&self.going),
            risers: self.risers,
            riser,
        }
    }
    /// The staged text a callout shows.
    pub(super) fn callout_value(&self, document: &Document, key: Key) -> String {
        if let Some((_, field)) = self.linked_fields().into_iter().find(|(k, _)| *k == key) {
            return diagram::field_value(field);
        }
        match key {
            visual::SLAB_OVERRIDE => self
                .level_slab_mm(document)
                .map_or("from level".into(), |mm| format!("{mm} mm")),
            visual::STAIR_RISERS => self.risers.to_string(),
            visual::STAIR_DIRECTION => format!("{}°", self.direction_deg),
            visual::STAIR_START => format!("{}, {} mm", self.start[0], self.start[1]),
            _ => String::new(),
        }
    }
    fn paint_diagram(&mut self, ui: &mut egui::Ui, document: &Document) {
        let mut link = std::mem::take(&mut self.link);
        link.invalid = diagram::first_invalid(document, self.linked_fields());
        let schematic = if self.kind == ArchitectureKind::Stair {
            visual::stair(self.form, &self.stair_values(document))
        } else {
            visual::slab(&self.slab_values(document))
        };
        let values: Vec<(Key, String)> = schematic
            .keys()
            .into_iter()
            .map(|k| (k, self.callout_value(document, k)))
            .collect();
        let value = |key: Key| {
            values
                .iter()
                .find(|(k, _)| *k == key)
                .map_or_else(String::new, |(_, v)| v.clone())
        };
        link.diagram(
            ui,
            ("architecture", self.kind.label(), self.feature),
            &schematic,
            240.,
            &value,
        );
        ui.small("Schematic of this draft · click a dimension to edit it");
        self.link = link;
    }
    pub(super) fn show_slab(
        &mut self,
        ui: &mut egui::Ui,
        document: &Document,
        levels: &[(EntityId, String)],
    ) -> bool {
        let mut changed = false;
        self.paint_diagram(ui, document);
        ui.label("Level");
        changed |= Self::level_combo(ui, "slab-level", levels, &mut self.level, None);
        let rooms = rooms(document);
        let room = rooms
            .iter()
            .find(|(id, _)| Some(*id) == self.room)
            .map_or("none".to_owned(), |(_, n)| n.clone());
        diagram::section(
            ui,
            "slab-room-section",
            &format!("Room boundary · {room}"),
            false,
            |ui| {
                Self::level_combo(
                    ui,
                    "slab-room",
                    &rooms,
                    &mut self.room,
                    Some("None (draw or enter the outline)"),
                );
                if self.room.is_some() {
                    ui.horizontal_wrapped(|ui| {
                        for b in RoomBoundary::ALL {
                            ui.radio_value(&mut self.room_boundary, b, b.label());
                        }
                    });
                    ui.small(format!(
                        "Only {room}'s boundary (green dashes in the viewport) is used; other rooms are not included. The outline is copied, so later wall edits do not move it."
                    ));
                    if ui.button("Use room boundary").clicked() {
                        match self.use_room(document) {
                            Ok(()) => changed = true,
                            Err(e) => self.error = Some(e),
                        }
                    }
                }
            },
        );
        let plane = format!("the level floor (z = {} mm)", self.plane_z(document));
        changed |= self.outline.show(ui, document, &plane);
        let mut link = std::mem::take(&mut self.link);
        if self.kind == ArchitectureKind::Pedestal {
            changed |= link.field(
                ui,
                visual::PEDESTAL_HEIGHT,
                &mut self.thickness,
                "pedestal-height",
                "Pedestal height",
                "Raised above the level's floor",
            );
        } else {
            let title = if self.override_thickness {
                format!("Thickness · {}", diagram::field_value(&self.thickness))
            } else {
                format!(
                    "Thickness · level {}",
                    self.callout_value(document, visual::SLAB_OVERRIDE)
                )
            };
            let reveal = link.reveals(&[visual::SLAB_OVERRIDE, visual::SLAB_THICKNESS]);
            diagram::section(ui, "slab-thickness-section", &title, reveal, |ui| {
                let response = ui.checkbox(
                    &mut self.override_thickness,
                    "Override the level's slab thickness",
                );
                link.track(visual::SLAB_OVERRIDE, &response);
                changed |= response.changed();
                if self.override_thickness {
                    changed |= link.field(
                        ui,
                        visual::SLAB_THICKNESS,
                        &mut self.thickness,
                        "slab-thickness",
                        "Thickness",
                        "Below the floor",
                    );
                }
            });
        }
        self.link = link;
        changed
    }
    pub(super) fn show_stair(
        &mut self,
        ui: &mut egui::Ui,
        document: &Document,
        levels: &[(EntityId, String)],
    ) -> bool {
        let mut changed = false;
        if ui
            .selectable_label(self.stair_placing, "Place in view")
            .clicked()
        {
            self.stair_placing = !self.stair_placing;
            self.stair_anchor = None;
            changed = true;
        }
        ui.small(if self.stair_placing {
            if self.stair_anchor.is_some() {
                "Click the walking direction from the first riser."
            } else {
                "Click the first riser on the base level, then the walking direction."
            }
        } else {
            "Use Place in view to reposition the stair graphically."
        });
        self.paint_diagram(ui, document);
        ui.label("From level");
        changed |= Self::level_combo(ui, "stair-level", levels, &mut self.level, None);
        ui.label("To level");
        changed |= Self::level_combo(
            ui,
            "stair-to-level",
            levels,
            &mut self.to_level,
            Some("Floor-to-floor height of the base level"),
        );
        ui.label("Form");
        let choices: Vec<Choice<'_, model::StairForm>> = visual::FORMS
            .iter()
            .map(|form| Choice {
                value: *form,
                label: visual::form_label(*form),
                hint: "Plan of the flights and landing; risers split evenly between flights",
                art: visual::form_art(*form),
            })
            .collect();
        changed |= diagram::choice_grid(ui, "stair-form", &mut self.form, &choices, true);
        let mut link = std::mem::take(&mut self.link);
        changed |= link.field(
            ui,
            visual::STAIR_WIDTH,
            &mut self.width,
            "stair-width",
            "Width",
            "Clear stair width",
        );
        let response = crate::numeric_field::NumericField::new("Risers", "count")
            .id(ui.make_persistent_id("stair-risers"))
            .value(ui, &mut self.risers, 2..=64)
            .input;
        link.track(visual::STAIR_RISERS, &response);
        changed |= response.changed();
        changed |= link.field(
            ui,
            visual::STAIR_GOING,
            &mut self.going,
            "stair-going",
            "Going",
            "Horizontal depth of one tread",
        );
        let title = format!(
            "Position · {}, {} mm · {}°",
            self.start[0], self.start[1], self.direction_deg
        );
        let reveal = link.reveals(&visual::STAIR_POSITION);
        diagram::section(ui, "stair-position-section", &title, reveal, |ui| {
            for (axis, (label, v)) in ["First riser start X", "First riser start Y"]
                .into_iter()
                .zip(self.start.iter_mut())
                .enumerate()
            {
                let response = crate::numeric_field::NumericField::new(label, "mm")
                    .id(ui.make_persistent_id(("stair-start", axis)))
                    .value(ui, v, f64::NEG_INFINITY..=f64::INFINITY)
                    .input;
                if axis == 0 {
                    link.track(visual::STAIR_START, &response);
                }
                changed |= response.changed();
            }
            let response = crate::numeric_field::NumericField::new("Direction", "deg")
                .id(ui.make_persistent_id("stair-direction"))
                .value(ui, &mut self.direction_deg, -360.0..=360.0)
                .input;
            link.track(visual::STAIR_DIRECTION, &response);
            changed |= response.changed();
            ui.checkbox(&mut self.stair_snap, "Snap placement to 10 mm");
        });
        self.link = link;
        if let Ok(stair) = self.stair(document)
            && let Ok(parameters) = parameters::evaluate(document)
            && let Ok(resolved) = model::resolve_stair(document, &parameters, &stair)
        {
            let c = resolved.check;
            ui.small(format!(
                "Rise {:.0} mm · riser {:.1} mm · going {:.0} mm · 2R+G {:.0} mm · run {:.0} mm",
                c.rise_mm, c.riser_mm, c.going_mm, c.step_measure_mm, c.run_mm
            ));
            for note in &c.notes {
                ui.small(note);
            }
            ui.small("Geometric measures only; not a building-code or accessibility approval.");
        }
        changed
    }
}
