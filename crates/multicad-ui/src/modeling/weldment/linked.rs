//! Linked weldment member inputs (CAD-395): the stock's actual section and
//! the member side view first; callouts focus the placement fields, path
//! points and end extensions, the latter two in named sections.
use super::{PathKind, WeldmentPanel, visual};
use crate::{
    egui,
    modeling::{
        diagram::{self, DiagramLink, Key},
        dimension::DimensionField,
    },
};
use multicad_core::{
    document::Document,
    materials::{self, LibraryDefinition},
    model::{AxisAnchor, ProfileSection},
    parameters,
};

impl WeldmentPanel {
    /// Linked fields in display order: rotation, then path points (when
    /// typed), then the end extensions.
    pub(super) fn linked_fields(&self) -> Vec<(Key, &DimensionField)> {
        let mut fields = vec![(visual::ROTATION, &self.rotation)];
        if self.typed_points() {
            fields.extend(visual::START.into_iter().zip(self.start.iter()));
            fields.extend(visual::END.into_iter().zip(self.end.iter()));
        }
        fields.extend([
            (visual::EXTEND_START, &self.extend_start),
            (visual::EXTEND_END, &self.extend_end),
        ]);
        fields
    }
    /// The member runs between the typed points (not a kept sketch line).
    fn typed_points(&self) -> bool {
        self.path_kind == PathKind::Points
    }
    /// The selected stock's resolved cross-section.
    pub(super) fn section(&self, document: &Document) -> Option<ProfileSection> {
        let entity = document.get(self.stock?)?;
        match materials::read(entity).ok().flatten()? {
            LibraryDefinition::Stock(stock) => ProfileSection::from_stock(&stock).ok(),
            _ => None,
        }
    }
    pub(super) fn diagram_values(&self, document: &Document) -> visual::Values {
        let evaluated = parameters::evaluate(document).ok();
        let value = |f: &DimensionField| {
            evaluated
                .as_ref()
                .and_then(|p| diagram::value_of(f, document, p))
        };
        let point = |p: &[DimensionField; 3]| -> Option<[f64; 3]> {
            Some([value(&p[0])?, value(&p[1])?, value(&p[2])?])
        };
        let length = if self.typed_points() {
            point(&self.start)
                .zip(point(&self.end))
                .map(|(a, b)| (0..3).map(|i| (b[i] - a[i]).powi(2)).sum::<f64>().sqrt())
                .filter(|l| *l > 0.)
        } else {
            None
        };
        visual::Values {
            section: self.section(document),
            anchor: self.anchor,
            rotation_deg: value(&self.rotation),
            extend: [value(&self.extend_start), value(&self.extend_end)],
            length,
            points: self.typed_points(),
        }
    }
    fn point_text(fields: &[DimensionField; 3]) -> String {
        let text: Vec<_> = fields.iter().map(|f| f.text.trim()).collect();
        format!("{} {}", text.join(", "), fields[0].unit())
    }
    pub(super) fn callout_value(&self, key: Key) -> String {
        match key {
            visual::ANCHOR => format!(
                "{} · {}",
                visual::anchor_label(self.anchor[0]),
                visual::anchor_label(self.anchor[1])
            ),
            k if k == visual::START[0] => Self::point_text(&self.start),
            k if k == visual::END[0] => Self::point_text(&self.end),
            _ => self
                .linked_fields()
                .into_iter()
                .find(|(k, _)| *k == key)
                .map_or_else(String::new, |(_, f)| diagram::field_value(f)),
        }
    }
    pub(super) fn paint_diagram(&mut self, ui: &mut egui::Ui, document: &Document) {
        let mut link = std::mem::take(&mut self.link);
        link.invalid = diagram::first_invalid(document, self.linked_fields());
        let schematic = visual::diagram(&self.diagram_values(document));
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
        link.diagram(ui, ("weldment", self.body), &schematic, 260., &value);
        ui.small("Stock section at the path anchor · click a callout to edit it");
        self.link = link;
    }
    pub(super) fn show_points(&mut self, ui: &mut egui::Ui) -> bool {
        let mut changed = false;
        let mut link = std::mem::take(&mut self.link);
        let keys: Vec<Key> = visual::START.into_iter().chain(visual::END).collect();
        let title = format!(
            "Path points · {} → {}",
            Self::point_text(&self.start),
            Self::point_text(&self.end)
        );
        let reveal = link.reveals(&keys);
        diagram::section(ui, "weld-points-section", &title, reveal, |ui| {
            for (i, axis) in ["X", "Y", "Z"].into_iter().enumerate() {
                changed |= link.field(
                    ui,
                    visual::START[i],
                    &mut self.start[i],
                    visual::START[i],
                    &format!("Start {axis}"),
                    "World path start",
                );
            }
            for (i, axis) in ["X", "Y", "Z"].into_iter().enumerate() {
                changed |= link.field(
                    ui,
                    visual::END[i],
                    &mut self.end[i],
                    visual::END[i],
                    &format!("End {axis}"),
                    "World path end",
                );
            }
        });
        self.link = link;
        changed
    }
    fn anchor_radio(
        ui: &mut egui::Ui,
        link: &mut DiagramLink,
        track: bool,
        label: &str,
        value: &mut AxisAnchor,
    ) -> bool {
        let mut changed = false;
        ui.horizontal(|ui| {
            ui.label(label);
            for a in [AxisAnchor::Min, AxisAnchor::Center, AxisAnchor::Max] {
                let response = ui.radio_value(value, a, visual::anchor_label(a));
                if track && a == AxisAnchor::Min {
                    link.track(visual::ANCHOR, &response);
                }
                changed |= response.changed();
            }
        });
        changed
    }
    pub(super) fn show_section_placement(&mut self, ui: &mut egui::Ui) -> bool {
        let mut changed = false;
        let mut link = std::mem::take(&mut self.link);
        ui.label("Section placement");
        changed |= Self::anchor_radio(ui, &mut link, true, "Width", &mut self.anchor[0]);
        changed |= Self::anchor_radio(ui, &mut link, false, "Height", &mut self.anchor[1]);
        changed |= link.field(
            ui,
            visual::ROTATION,
            &mut self.rotation,
            "weld-rotation",
            "Rotation",
            "About the member axis; height starts along the sketch normal, or up",
        );
        let title = format!(
            "End extensions · {} / {}",
            diagram::field_value(&self.extend_start),
            diagram::field_value(&self.extend_end)
        );
        let reveal = link.reveals(&visual::EXTENSIONS);
        diagram::section(ui, "weld-extensions-section", &title, reveal, |ui| {
            changed |= link.field(
                ui,
                visual::EXTEND_START,
                &mut self.extend_start,
                "weld-extend-start",
                "Extend start",
                "Extra length at the path start; negative trims",
            );
            changed |= link.field(
                ui,
                visual::EXTEND_END,
                &mut self.extend_end,
                "weld-extend-end",
                "Extend end",
                "Extra length at the path end; negative trims",
            );
        });
        self.link = link;
        changed
    }
}
