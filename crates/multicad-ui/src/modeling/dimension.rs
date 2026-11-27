use crate::egui;
use crate::numeric_field::{Emphasis, NumericField};
use multicad_core::{
    constraint_sketch::Dimension,
    document::Document,
    parameters::{BoundExpression, LengthUnit, Unit, parse_length},
};
/// Up to six decimals without trailing zeros.
fn trimmed(value: f64) -> String {
    let text = format!("{value:.6}");
    let text = text.trim_end_matches('0').trim_end_matches('.');
    if text == "-0" {
        "0".into()
    } else {
        text.into()
    }
}

/// Text stays local, including temporarily incomplete formulas. Stable bindings
/// are retained when editing a previously bound dimension.
#[derive(Clone, Debug)]
pub(crate) struct DimensionField {
    pub text: String,
    unit: String,
    previous: Option<BoundExpression>,
    original: Dimension,
    original_text: String,
    /// The latest resolve error, shown below the field (CAD-416).
    error: LastError,
}
const HINT: &str = "Value or =$parameter";
/// Interior slot so `resolve(&self)` callers need not change; a clone copies
/// the message rather than sharing it.
#[derive(Debug, Default)]
struct LastError(std::sync::Mutex<Option<String>>);
impl LastError {
    fn get(&self) -> Option<String> {
        self.0.lock().map(|error| error.clone()).unwrap_or_default()
    }
    fn set(&self, error: Option<String>) {
        if let Ok(mut slot) = self.0.lock() {
            *slot = error;
        }
    }
}
impl Clone for LastError {
    fn clone(&self) -> Self {
        Self(std::sync::Mutex::new(self.get()))
    }
}
impl DimensionField {
    /// `unit` is the field's canonical unit. Millimetre fields show and read
    /// unqualified input in the Settings length unit (CAD-187); typed values
    /// are stored with that explicit unit, so the preference never changes a
    /// stored dimension.
    pub fn from_dimension(value: &Dimension, unit: &'static str) -> Result<Self, String> {
        let unit = if unit == "mm" {
            crate::preferences::length_unit().symbol()
        } else {
            unit
        };
        let source = Unit::parse(value.unit()).ok_or("Unknown dimension unit")?;
        let target = Unit::parse(unit).ok_or("Unknown field unit")?;
        if source.dimension() != target.dimension() {
            return Err("Incompatible dimension unit".into());
        }
        let text = match value {
            Dimension::Literal { value, .. } if source == target => value.to_string(),
            // A converted literal shows a rounded value; left unchanged it
            // still resolves to the exact authored dimension.
            Dimension::Literal { value, .. } => trimmed(value * source.scale() / target.scale()),
            Dimension::Expression { expression } => expression.text.clone(),
        };
        Ok(Self {
            original_text: text.clone(),
            original: value.clone(),
            text,
            unit: if value.expression().is_some() {
                value.unit()
            } else {
                unit
            }
            .into(),
            previous: value.expression().cloned(),
            error: LastError::default(),
        })
    }
    /// Replace the text with a literal given in `unit`, converted to this
    /// field's own unit (a formula field may keep an authored unit).
    pub fn set_literal(&mut self, value: f64, unit: &str) {
        let scale = match (Unit::parse(unit), Unit::parse(&self.unit)) {
            (Some(from), Some(to)) if from.dimension() == to.dimension() => {
                from.scale() / to.scale()
            }
            _ => 1.,
        };
        self.text = super::value_entry::format_number(value * scale);
    }
    /// Keep an exact authored pick while displaying it in this field's unit.
    /// Unedited rounded display text resolves to the original canonical value,
    /// just as it does when reopening a stored literal.
    pub fn set_exact_literal(&mut self, value: f64, unit: &str) {
        let (Some(from), Some(to)) = (Unit::parse(unit), Unit::parse(&self.unit)) else {
            return;
        };
        if from.dimension() != to.dimension() {
            return;
        }
        self.text = trimmed(value * from.scale() / to.scale());
        self.original_text = self.text.clone();
        self.original = Dimension::Literal {
            value,
            unit: unit.into(),
        };
        self.previous = None;
    }
    pub fn show(&mut self, ui: &mut egui::Ui, label: &str, help: &str) -> bool {
        self.paint(ui, None, label, help, Emphasis::Normal)
            .changed()
    }
    /// Stable identity for forms whose validation or preview messages can insert
    /// widgets before an active edit. Other tools keep the existing auto-ID path.
    pub fn show_named(
        &mut self,
        ui: &mut egui::Ui,
        id_salt: impl egui::AsIdSalt,
        label: &str,
        help: &str,
    ) -> bool {
        let id = ui.make_persistent_id(id_salt);
        self.paint(ui, Some(id), label, help, Emphasis::Normal)
            .changed()
    }
    /// The unit unqualified input is read in.
    pub fn unit(&self) -> &str {
        &self.unit
    }
    /// A named field linked to a diagram callout: `highlighted` marks the
    /// measurement the diagram currently identifies and `invalid` the field
    /// validation refers to. Returns the editor response for focus linking.
    pub fn show_linked(
        &mut self,
        ui: &mut egui::Ui,
        id_salt: impl egui::AsIdSalt,
        label: &str,
        help: &str,
        highlighted: bool,
        invalid: bool,
    ) -> egui::Response {
        let id = ui.make_persistent_id(id_salt);
        let emphasis = if invalid {
            Emphasis::Invalid
        } else if highlighted {
            Emphasis::Highlighted
        } else {
            Emphasis::Normal
        };
        self.paint(ui, Some(id), label, help, emphasis)
    }
    /// Formerly a two-column grid cell; CAD-416 gives every numeric Inspector
    /// field the same full-width layout, so this is [`Self::show_named`].
    pub fn show_compact(
        &mut self,
        ui: &mut egui::Ui,
        id_salt: impl egui::AsIdSalt,
        label: &str,
        help: &str,
    ) -> bool {
        self.show_named(ui, id_salt, label, help)
    }
    /// CAD-416: label above, full-width input, unit inside on the right and
    /// the last resolve error below.
    fn paint(
        &mut self,
        ui: &mut egui::Ui,
        id: Option<egui::Id>,
        label: &str,
        help: &str,
        emphasis: Emphasis,
    ) -> egui::Response {
        let error = self.error.get();
        let unit = self.unit.clone();
        let mut field = NumericField::new(label, &unit);
        if let Some(id) = id {
            field = field.id(id);
        }
        let response = field
            .help(help)
            .hint(HINT)
            .emphasis(emphasis)
            .error(error.as_deref())
            .show(ui, &mut self.text)
            .input;
        note_enter(ui, &response);
        response
    }
    pub fn resolve(&self, document: &Document) -> Result<Dimension, String> {
        let result = self.resolve_text(document);
        self.error.set(result.as_ref().err().cloned());
        result
    }
    fn resolve_text(&self, document: &Document) -> Result<Dimension, String> {
        if self.text.trim_start().starts_with('=') {
            return BoundExpression::bind(
                document,
                self.text.trim(),
                &self.unit,
                self.previous.as_ref(),
            )
            .map(|expression| Dimension::Expression { expression })
            .map_err(|error| error.to_string());
        }
        let (value, unit) = match self.text.trim().parse::<f64>() {
            Ok(value) => (value, self.unit.clone()),
            // Explicitly qualified lengths such as `1 in` keep their own unit.
            Err(_) => Unit::parse(&self.unit)
                .filter(|unit| unit.dimension() == multicad_core::parameters::Dimension::LENGTH)
                .and_then(|_| parse_length(&self.text, LengthUnit::Millimeter).ok())
                .map(|input| (input.value, input.unit.symbol().to_string()))
                .ok_or_else(|| {
                    "Enter a number, a length such as 1 in, or an expression beginning with '='"
                        .to_string()
                })?,
        };
        if !value.is_finite() || value.abs() > 1e9 {
            return Err("Value must be finite and no greater than 1,000,000,000".into());
        }
        if self.text.trim() == self.original_text.trim()
            && matches!(self.original, Dimension::Literal { .. })
        {
            return Ok(self.original.clone());
        }
        Ok(Dimension::Literal { value, unit })
    }
    pub fn remember(&mut self, dimension: &Dimension) {
        self.previous = dimension.expression().cloned();
        self.original = dimension.clone();
        self.original_text = self.text.clone();
    }
}

const FIELD_ENTER: &str = "multicad.modeling.field-enter";
/// CAD-410: Enter that ends editing a numeric field asks the active tool to
/// apply, exactly like Apply, once its value is valid.
fn note_enter(ui: &egui::Ui, response: &egui::Response) {
    if response.lost_focus() && ui.input(|input| input.key_pressed(egui::Key::Enter)) {
        let pass = ui.ctx().cumulative_pass_nr();
        ui.ctx()
            .data_mut(|data| data.insert_temp(egui::Id::new(FIELD_ENTER), pass));
    }
}
/// Whether a numeric field took Enter in this pass; clears the request.
pub(crate) fn take_field_enter(ctx: &egui::Context) -> bool {
    let pass = ctx.cumulative_pass_nr();
    ctx.data_mut(|data| data.remove_temp::<u64>(egui::Id::new(FIELD_ENTER)))
        .is_some_and(|noted| noted == pass)
}

#[cfg(test)]
mod tests {
    use super::*;
    use multicad_core::document::{Entity, EntityData};
    /// CAD-416: every DimensionField unit kind shows inside its input, right
    /// of the value, and never in the label or the edited text.
    #[test]
    fn units_paint_inside_the_input_and_counts_show_none() {
        for (unit, value, symbol) in [
            ("mm", 12., Some("mm")),
            ("deg", 30., Some("°")),
            ("in", 2., Some("in")),
            ("count", 6., None),
        ] {
            let mut field = DimensionField::from_dimension(
                &Dimension::Literal {
                    value,
                    unit: unit.into(),
                },
                unit,
            )
            .unwrap();
            let ctx = egui::Context::default();
            let mut input = egui::Rect::NOTHING;
            let mut out = ctx.run_ui(Default::default(), |ui| {
                ui.set_width(220.);
                input = field
                    .show_linked(ui, ("units", unit), "Size", "", false, false)
                    .rect;
            });
            let texts: Vec<(String, egui::Rect)> = out
                .shapes
                .iter()
                .filter_map(|c| match &c.shape {
                    egui::Shape::Text(t) => Some((
                        t.galley.text().to_owned(),
                        t.galley.rect.translate(t.pos.to_vec2()),
                    )),
                    _ => None,
                })
                .collect();
            out.textures_delta.clear();
            let inside: Vec<&str> = texts
                .iter()
                .filter(|(_, r)| input.expand(0.5).contains_rect(*r))
                .map(|(t, _)| t.as_str())
                .collect();
            assert!(texts.iter().any(|(t, _)| t == "Size"), "{unit}: label");
            assert_eq!(field.text, value.to_string(), "{unit}: value text only");
            match symbol {
                Some(symbol) => assert_eq!(inside, [field.text.as_str(), symbol], "{unit}"),
                None => assert_eq!(inside, [field.text.as_str()], "{unit}"),
            }
        }
    }
    #[test]
    fn unchanged_literal_cm_and_rad_retain_authored_units_until_edited() {
        let doc = Document::new("Unit baseline").unwrap();
        for (original, target) in [
            (
                Dimension::Literal {
                    value: 2.,
                    unit: "cm".into(),
                },
                "mm",
            ),
            (
                Dimension::Literal {
                    value: std::f64::consts::FRAC_PI_2,
                    unit: "rad".into(),
                },
                "deg",
            ),
        ] {
            let mut field = DimensionField::from_dimension(&original, target).unwrap();
            assert_eq!(field.resolve(&doc).unwrap(), original);
            field.remember(&original);
            assert_eq!(field.resolve(&doc).unwrap(), original);
            field.text = "25".into();
            let edited = field.resolve(&doc).unwrap();
            assert_eq!(
                edited,
                Dimension::Literal {
                    value: 25.,
                    unit: target.into()
                }
            );
            field.remember(&edited);
            assert_eq!(field.resolve(&doc).unwrap(), edited);
        }
    }
    #[test]
    fn formula_replacement_remembers_stable_binding_and_still_rejects_deleted_parameters() {
        let mut doc = Document::new("Formula baseline").unwrap();
        let mut parameter = Entity::new(
            "width",
            EntityData::Parameter {
                expression: "2".into(),
                unit: "cm".into(),
                exposed: true,
                inputs: vec![],
            },
        );
        doc.insert(parameter.clone()).unwrap();
        let mut field = DimensionField::from_dimension(
            &Dimension::Literal {
                value: 1.,
                unit: "cm".into(),
            },
            "mm",
        )
        .unwrap();
        field.text = "=$width + 5".into();
        let bound = field.resolve(&doc).unwrap();
        field.remember(&bound);
        parameter.name = "renamed_width".into();
        doc.replace(parameter.clone()).unwrap();
        assert_eq!(field.resolve(&doc).unwrap(), bound);
        doc.remove(parameter.id).unwrap();
        assert!(field.resolve(&doc).is_err());
    }
    #[test]
    fn length_fields_follow_the_settings_unit_without_changing_stored_dimensions() {
        use crate::preferences::{LengthUnit, set_length_unit};
        let doc = Document::new("Inch fields").unwrap();
        let authored = Dimension::Literal {
            value: 25.4,
            unit: "mm".into(),
        };
        set_length_unit(LengthUnit::Inch);
        let mut field = DimensionField::from_dimension(&authored, "mm").unwrap();
        assert_eq!(field.text, "1");
        assert_eq!(field.unit, "in");
        // Untouched, the authored millimetre literal is kept exactly.
        assert_eq!(field.resolve(&doc).unwrap(), authored);
        // Unqualified input is inches; qualified input keeps its unit.
        field.text = "2".into();
        assert_eq!(
            field.resolve(&doc).unwrap(),
            Dimension::Literal {
                value: 2.,
                unit: "in".into()
            }
        );
        field.text = "30 mm".into();
        assert_eq!(
            field.resolve(&doc).unwrap(),
            Dimension::Literal {
                value: 30.,
                unit: "mm".into()
            }
        );
        field.text = "5 deg".into();
        assert!(field.resolve(&doc).is_err());
        // Angles and counts are unaffected by the length preference.
        let angle = DimensionField::from_dimension(
            &Dimension::Literal {
                value: 90.,
                unit: "deg".into(),
            },
            "deg",
        )
        .unwrap();
        assert_eq!(angle.unit, "deg");
        // An inch literal shown in millimetres stays inches until edited.
        set_length_unit(LengthUnit::Millimeter);
        let inch = Dimension::Literal {
            value: 1.,
            unit: "in".into(),
        };
        let field = DimensionField::from_dimension(&inch, "mm").unwrap();
        assert_eq!((field.text.as_str(), field.unit.as_str()), ("25.4", "mm"));
        assert_eq!(field.resolve(&doc).unwrap(), inch);
        assert_eq!(
            inch.value_in(Unit::Millimeter, None).unwrap(),
            authored.value_in(Unit::Millimeter, None).unwrap()
        );
    }
}
