//! Linked FDM aid dimensions (CAD-395): the aid schematic first, its
//! callouts focusing the staged fields below; contact, anchor or footprint
//! coordinates sit in a named section that callouts and validation open.
use super::{FdmAidPanel, visual};
use crate::{
    egui,
    modeling::{
        diagram::{self, Key},
        dimension::DimensionField,
    },
};
use multicad_core::{document::Document, model::AidKind, parameters};

impl FdmAidPanel {
    /// The fields of the current kind with their callout keys, essentials
    /// first and then the named position section, in display order.
    pub(super) fn linked_fields(&self) -> Vec<(Key, &DimensionField)> {
        let o = &self.origin;
        match self.kind {
            AidKind::Tab => vec![
                (visual::LENGTH, &self.length),
                (visual::OVERLAP, &self.overlap),
                (visual::WIDTH, &self.width),
                (visual::THICKNESS, &self.thickness),
                (visual::ORIGIN[0], &o[0]),
                (visual::ORIGIN[1], &o[1]),
                (visual::ORIGIN[2], &o[2]),
                (visual::ANGLE, &self.angle),
            ],
            AidKind::Bridge => vec![
                (visual::WIDTH, &self.width),
                (visual::THICKNESS, &self.thickness),
                (visual::ORIGIN[0], &o[0]),
                (visual::ORIGIN[1], &o[1]),
                (visual::ORIGIN[2], &o[2]),
                (visual::END[0], &self.end[0]),
                (visual::END[1], &self.end[1]),
                (visual::END[2], &self.end[2]),
            ],
            AidKind::Support => vec![
                (visual::LENGTH, &self.length),
                (visual::WIDTH, &self.width),
                (visual::BASE, &self.base),
                (visual::GAP, &self.gap),
                (visual::ORIGIN[0], &o[0]),
                (visual::ORIGIN[1], &o[1]),
                (visual::ANGLE, &self.angle),
            ],
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
        let span = point(&self.origin)
            .zip(point(&self.end))
            .map(|(a, b)| (0..3).map(|i| (b[i] - a[i]).powi(2)).sum::<f64>().sqrt())
            .filter(|s| *s > 0.);
        visual::Values {
            length: value(&self.length),
            width: value(&self.width),
            thickness: value(&self.thickness),
            overlap: value(&self.overlap),
            gap: value(&self.gap),
            span,
            second_part: self.other.is_some(),
        }
    }
    fn coordinates(fields: &[DimensionField]) -> String {
        let text: Vec<_> = fields.iter().map(|f| f.text.trim()).collect();
        format!(
            "{} {}",
            text.join(", "),
            fields.first().map_or("mm", |f| f.unit())
        )
    }
    /// The staged text a callout shows: coordinates read as one point.
    pub(super) fn callout_value(&self, key: Key) -> String {
        if key == visual::ORIGIN[0] {
            let n = if self.kind == AidKind::Support { 2 } else { 3 };
            return Self::coordinates(&self.origin[..n]);
        }
        if key == visual::END[0] {
            return Self::coordinates(&self.end);
        }
        self.linked_fields()
            .into_iter()
            .find(|(k, _)| *k == key)
            .map_or_else(String::new, |(_, f)| diagram::field_value(f))
    }
    pub(super) fn show_linked(&mut self, ui: &mut egui::Ui, document: &Document) -> bool {
        let mut changed = false;
        let mut link = std::mem::take(&mut self.link);
        link.invalid = diagram::first_invalid(document, self.linked_fields());
        let schematic = visual::diagram(self.kind, &self.diagram_values(document));
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
        link.diagram(ui, ("fdm-aid", self.feature), &schematic, 240., &value);
        ui.small("Schematic · finished Part dashed · click a dimension to edit it");
        let mut field = |link: &mut diagram::DiagramLink,
                         ui: &mut egui::Ui,
                         key: Key,
                         field: &mut DimensionField,
                         label: &str,
                         help: &str| {
            changed |= link.field(ui, key, field, key, label, help);
        };
        let kind = self.kind;
        match kind {
            AidKind::Tab => {
                field(
                    &mut link,
                    ui,
                    visual::LENGTH,
                    &mut self.length,
                    "Reach",
                    "Outward from the contact",
                );
                field(
                    &mut link,
                    ui,
                    visual::OVERLAP,
                    &mut self.overlap,
                    "Overlap",
                    "Back under the finished Part",
                );
                field(
                    &mut link,
                    ui,
                    visual::WIDTH,
                    &mut self.width,
                    "Width",
                    "Across the pad",
                );
                field(
                    &mut link,
                    ui,
                    visual::THICKNESS,
                    &mut self.thickness,
                    "Thickness",
                    "Pad height above contact",
                );
            }
            AidKind::Bridge => {
                field(
                    &mut link,
                    ui,
                    visual::WIDTH,
                    &mut self.width,
                    "Width",
                    "Horizontal bar width",
                );
                field(
                    &mut link,
                    ui,
                    visual::THICKNESS,
                    &mut self.thickness,
                    "Thickness",
                    "Bar cross-section thickness",
                );
            }
            AidKind::Support => {
                field(
                    &mut link,
                    ui,
                    visual::LENGTH,
                    &mut self.length,
                    "Length",
                    "Footprint length",
                );
                field(
                    &mut link,
                    ui,
                    visual::WIDTH,
                    &mut self.width,
                    "Width",
                    "Footprint width",
                );
                field(
                    &mut link,
                    ui,
                    visual::BASE,
                    &mut self.base,
                    "Base",
                    "Bottom height, usually the plate",
                );
                field(
                    &mut link,
                    ui,
                    visual::GAP,
                    &mut self.gap,
                    "Top gap",
                    "Distance below the lowest target surface",
                );
            }
        }
        let keys = visual::position_keys(kind);
        let (name, summary) = match kind {
            AidKind::Tab => (
                "Contact and direction",
                format!(
                    "{} · {}",
                    Self::coordinates(&self.origin),
                    diagram::field_value(&self.angle)
                ),
            ),
            AidKind::Bridge => (
                "Anchors",
                format!(
                    "{} → {}",
                    Self::coordinates(&self.origin),
                    Self::coordinates(&self.end)
                ),
            ),
            AidKind::Support => (
                "Footprint position",
                format!(
                    "{} · {}",
                    Self::coordinates(&self.origin[..2]),
                    diagram::field_value(&self.angle)
                ),
            ),
        };
        let reveal = link.reveals(&keys);
        diagram::section(
            ui,
            ("fdm-aid-position", name),
            &format!("{name} · {summary}"),
            reveal,
            |ui| match kind {
                AidKind::Tab => {
                    for (i, axis) in ["X", "Y", "Z"].into_iter().enumerate() {
                        field(
                            &mut link,
                            ui,
                            visual::ORIGIN[i],
                            &mut self.origin[i],
                            &format!("Contact {axis}"),
                            "World contact point at the base edge",
                        );
                    }
                    field(
                        &mut link,
                        ui,
                        visual::ANGLE,
                        &mut self.angle,
                        "Direction",
                        "Outward about +Z",
                    );
                }
                AidKind::Bridge => {
                    for (i, axis) in ["X", "Y", "Z"].into_iter().enumerate() {
                        field(
                            &mut link,
                            ui,
                            visual::ORIGIN[i],
                            &mut self.origin[i],
                            &format!("Start {axis}"),
                            "First surface anchor",
                        );
                    }
                    for (i, axis) in ["X", "Y", "Z"].into_iter().enumerate() {
                        field(
                            &mut link,
                            ui,
                            visual::END[i],
                            &mut self.end[i],
                            &format!("End {axis}"),
                            "Second surface anchor",
                        );
                    }
                }
                AidKind::Support => {
                    for (i, axis) in ["X", "Y"].into_iter().enumerate() {
                        field(
                            &mut link,
                            ui,
                            visual::ORIGIN[i],
                            &mut self.origin[i],
                            &format!("Corner {axis}"),
                            "Support footprint corner",
                        );
                    }
                    field(
                        &mut link,
                        ui,
                        visual::ANGLE,
                        &mut self.angle,
                        "Direction",
                        "Footprint orientation about +Z",
                    );
                }
            },
        );
        if kind == AidKind::Support {
            ui.small("Height is derived from the underside during native Preview.");
        }
        self.link = link;
        changed
    }
}
