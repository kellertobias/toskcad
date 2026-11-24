//! Wood box schematic values and callout text (CAD-395).
use super::{WoodBoxPanel, visual};
use crate::{
    egui,
    modeling::{
        diagram::{self, DiagramLink, Key},
        dimension::DimensionField,
    },
};
use multicad_core::{
    document::Document,
    materials::{self, LibraryDefinition, StockKind},
    parameters,
};

impl WoodBoxPanel {
    /// Linked fields in display order: size, then the placement section.
    pub(super) fn linked_fields(&self) -> Vec<(Key, &DimensionField)> {
        let mut fields: Vec<(Key, &DimensionField)> =
            visual::SIZE.into_iter().zip(self.size.iter()).collect();
        fields.extend(visual::ORIGIN.into_iter().zip(self.origin.iter()));
        fields.push((visual::ROTATION, &self.rotation));
        fields
    }
    fn thickness(&self, document: &Document) -> Option<f64> {
        let entity = document.get(self.stock?)?;
        match materials::read(entity).ok().flatten()? {
            LibraryDefinition::Stock(stock) => match stock.specification {
                StockKind::Panel { thickness_mm, .. } => Some(thickness_mm),
                _ => None,
            },
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
        visual::Values {
            size: [
                value(&self.size[0]),
                value(&self.size[1]),
                value(&self.size[2]),
            ],
            faces: self.source.definition.faces.iter().copied().collect(),
            thickness: self.thickness(document),
        }
    }
    /// The staged text a callout shows; the origin reads as one point.
    pub(super) fn callout_value(&self, key: Key) -> String {
        if key == visual::ORIGIN[0] {
            let text: Vec<_> = self.origin.iter().map(|f| f.text.trim()).collect();
            return format!("{} {}", text.join(", "), self.origin[0].unit());
        }
        self.linked_fields()
            .into_iter()
            .find(|(k, _)| *k == key)
            .map_or_else(String::new, |(_, f)| diagram::field_value(f))
    }
    pub(super) fn paint_diagram(
        &self,
        ui: &mut egui::Ui,
        document: &Document,
        link: &mut DiagramLink,
    ) {
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
        link.diagram(ui, ("wood-box", self.source.id), &schematic, 200., &value);
        ui.small("Outside envelope · omitted faces dashed · click a dimension to edit it");
    }
}
