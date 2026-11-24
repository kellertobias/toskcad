//! Values and callout text for the relief schematic (CAD-395).
use super::{Direction, Pattern, SurfaceStructurePanel, visual};
use crate::modeling::{
    diagram::{self, Key},
    dimension::DimensionField,
};
use multicad_core::{document::Document, parameters};

impl SurfaceStructurePanel {
    pub(super) fn kind(&self) -> visual::Kind {
        match self.pattern {
            Pattern::Ribs => visual::Kind::Ribs,
            Pattern::Waves => visual::Kind::Waves,
            Pattern::Stipple => visual::Kind::Stipple,
            Pattern::HeightMap { .. } => visual::Kind::HeightMap,
        }
    }
    pub(super) fn image(&self) -> Option<[u32; 2]> {
        match &self.pattern {
            Pattern::HeightMap { width, height, .. } => {
                Some([u32::from(*width), u32::from(*height)])
            }
            _ => None,
        }
    }
    /// Linked fields in display order; a height map has no pitch.
    pub(super) fn linked_fields(&self) -> Vec<(Key, &DimensionField)> {
        let mut fields = vec![];
        if self.image().is_none() {
            fields.push((visual::PITCH, &self.pitch));
        }
        fields.extend([
            (visual::AMPLITUDE, &self.amplitude),
            (visual::SIZE[0], &self.size[0]),
            (visual::SIZE[1], &self.size[1]),
            (visual::OFFSET[0], &self.offset[0]),
            (visual::OFFSET[1], &self.offset[1]),
            (visual::ROTATION, &self.rotation),
            (visual::BORDER, &self.border),
        ]);
        fields
    }
    pub(super) fn diagram_values(&self, document: &Document) -> visual::Values {
        let evaluated = parameters::evaluate(document).ok();
        let value = |f: &DimensionField| {
            evaluated
                .as_ref()
                .and_then(|p| diagram::value_of(f, document, p))
        };
        visual::Values {
            pitch: value(&self.pitch),
            amplitude: value(&self.amplitude),
            size: [value(&self.size[0]), value(&self.size[1])],
            offset: [value(&self.offset[0]), value(&self.offset[1])],
            rotation_deg: value(&self.rotation),
            border: value(&self.border),
            recessed: self.direction == Direction::Recessed,
        }
    }
    /// The staged text a callout shows; the offset reads as one U, V pair.
    pub(super) fn callout_value(&self, key: Key) -> String {
        if key == visual::OFFSET[0] {
            return format!(
                "{}, {} {}",
                self.offset[0].text.trim(),
                self.offset[1].text.trim(),
                self.offset[0].unit()
            );
        }
        self.linked_fields()
            .into_iter()
            .find(|(k, _)| *k == key)
            .map_or_else(String::new, |(_, f)| diagram::field_value(f))
    }
}
