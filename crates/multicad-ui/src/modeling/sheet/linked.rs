//! Linked base-sheet inputs (CAD-395): the bend section and base plan
//! first; callouts focus the size fields, and the base position and the
//! inside-radius / K-factor overrides sit in named sections that callouts
//! and validation open.
use super::{SheetPanel, plane_label, sheet_stocks, visual};
use crate::{
    controls, egui,
    modeling::{
        diagram::{self, Key},
        dimension::DimensionField,
    },
};
use multicad_core::{document::Document, model::PanelPlane, parameters};

impl SheetPanel {
    /// Linked fields shown in the base section, in display order.
    pub(super) fn linked_fields(&self) -> Vec<(Key, &DimensionField)> {
        let mut fields = vec![];
        if self.kept_base.is_none() {
            fields.extend([(visual::LENGTH, &self.length), (visual::WIDTH, &self.width)]);
            fields.extend(visual::ORIGIN.into_iter().zip(self.origin.iter()));
        }
        if self.radius_override {
            fields.push((visual::RADIUS, &self.radius));
        }
        if self.k_override {
            fields.push((visual::K_FACTOR, &self.k_factor));
        }
        fields
    }
    fn stock_values(&self, document: &Document) -> Option<(f64, f64, f64)> {
        sheet_stocks(document)
            .into_iter()
            .find(|s| Some(s.0) == self.stock)
            .map(|s| (s.2, s.3, s.4))
    }
    pub(super) fn diagram_values(&self, document: &Document) -> visual::Values {
        let evaluated = parameters::evaluate(document).ok();
        let value = |f: &DimensionField| {
            evaluated
                .as_ref()
                .and_then(|p| diagram::value_of(f, document, p))
        };
        visual::Values {
            stock: self.stock_values(document),
            length: value(&self.length),
            width: value(&self.width),
            radius_override: self.radius_override,
            radius: value(&self.radius),
            k_override: self.k_override,
            k_factor: value(&self.k_factor),
            reverse: self.reverse,
            kept_base: self.kept_base.is_some(),
        }
    }
    pub(super) fn callout_value(&self, document: &Document, key: Key) -> String {
        let stock = self.stock_values(document);
        match key {
            visual::RADIUS_OVERRIDE => stock.map_or("from stock".into(), |s| format!("{} mm", s.1)),
            visual::K_OVERRIDE => stock.map_or("from stock".into(), |s| s.2.to_string()),
            _ => self
                .linked_fields()
                .into_iter()
                .find(|(k, _)| *k == key)
                .map_or_else(String::new, |(_, f)| diagram::field_value(f)),
        }
    }
    /// The base inputs below the stock choice. Returns whether they changed.
    pub(super) fn show_base(&mut self, ui: &mut egui::Ui, document: &Document) -> bool {
        let mut changed = false;
        let mut link = std::mem::take(&mut self.link);
        link.invalid = diagram::first_invalid(document, self.linked_fields());
        let schematic = visual::diagram(&self.diagram_values(document));
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
        link.diagram(ui, ("sheet-base", self.feature), &schematic, 220., &value);
        ui.small("Representative bend · neutral axis dashed · click a callout to edit it");
        if self.kept_base.is_some() {
            ui.small("The base comes from a sketch region or a converted solid and is kept.");
        } else {
            ui.label("Plane");
            controls::inspector_combo(ui, "sheet-plane")
                .selected_text(plane_label(self.plane))
                .show_ui(ui, |ui| {
                    for plane in [PanelPlane::Xy, PanelPlane::Xz, PanelPlane::Yz] {
                        changed |= ui
                            .selectable_value(&mut self.plane, plane, plane_label(plane))
                            .changed();
                    }
                });
            changed |= link.field(
                ui,
                visual::LENGTH,
                &mut self.length,
                "sheet-length",
                "Length (u)",
                "Base size along the plane's first axis",
            );
            changed |= link.field(
                ui,
                visual::WIDTH,
                &mut self.width,
                "sheet-width",
                "Width (v)",
                "Base size along the plane's second axis",
            );
        }
        changed |= ui
            .checkbox(&mut self.reverse, "Material on the other side")
            .changed();
        if self.kept_base.is_none() {
            let text: Vec<_> = self.origin.iter().map(|f| f.text.trim()).collect();
            let title = format!(
                "Base position · {} {}",
                text.join(", "),
                self.origin[0].unit()
            );
            let reveal = link.reveals(&visual::ORIGIN);
            diagram::section(ui, "sheet-position-section", &title, reveal, |ui| {
                for (i, axis) in ["X", "Y", "Z"].into_iter().enumerate() {
                    changed |= link.field(
                        ui,
                        visual::ORIGIN[i],
                        &mut self.origin[i],
                        visual::ORIGIN[i],
                        &format!("Origin {axis}"),
                        "World position of the base sheet's minimum corner",
                    );
                }
            });
        }
        let stock = self.stock_values(document);
        let radius = if self.radius_override {
            diagram::field_value(&self.radius)
        } else {
            stock.map_or("stock".into(), |s| format!("stock {} mm", s.1))
        };
        let k = if self.k_override {
            diagram::field_value(&self.k_factor)
        } else {
            stock.map_or("stock".into(), |s| format!("stock {}", s.2))
        };
        let title = format!("Bend overrides · R {radius} · K {k}");
        let reveal = link.reveals(&visual::BEND);
        diagram::section(ui, "sheet-bend-section", &title, reveal, |ui| {
            let response = ui.checkbox(&mut self.radius_override, "Override inside radius");
            link.track(visual::RADIUS_OVERRIDE, &response);
            changed |= response.changed();
            if self.radius_override {
                changed |= link.field(
                    ui,
                    visual::RADIUS,
                    &mut self.radius,
                    "sheet-radius",
                    "Inside radius",
                    "Inside bend radius for this part instead of the stock value",
                );
            }
            let response = ui.checkbox(&mut self.k_override, "Override K-factor");
            link.track(visual::K_OVERRIDE, &response);
            changed |= response.changed();
            if self.k_override {
                changed |= link.field(
                    ui,
                    visual::K_FACTOR,
                    &mut self.k_factor,
                    "sheet-k-factor",
                    "K-factor",
                    "Neutral-axis position as a fraction of thickness, 0–1",
                );
            }
        });
        self.link = link;
        changed
    }
}
