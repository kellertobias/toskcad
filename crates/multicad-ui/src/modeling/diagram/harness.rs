//! Painted-form test harness shared by the linked-diagram forms (CAD-395).
//!
//! Each frame runs a real egui pass at a given Inspector width and records
//! every painted text with its screen rectangle, so tests can click a
//! callout where it is actually drawn and type into the focused field.
use crate::egui;

pub(crate) struct Painted {
    pub texts: Vec<(String, egui::Rect)>,
}
impl Painted {
    /// The first painted text starting with `prefix`.
    pub(crate) fn find(&self, prefix: &str) -> Option<egui::Rect> {
        self.texts
            .iter()
            .find(|(t, _)| t.starts_with(prefix))
            .map(|(_, r)| *r)
    }
    pub(crate) fn has(&self, text: &str) -> bool {
        self.texts.iter().any(|(t, _)| t == text)
    }
    pub(crate) fn has_prefix(&self, prefix: &str) -> bool {
        self.find(prefix).is_some()
    }
}

pub(crate) fn frame(
    ctx: &egui::Context,
    width: f32,
    events: Vec<egui::Event>,
    show: &mut dyn FnMut(&mut egui::Ui),
) -> Painted {
    let output = ctx.run_ui(
        egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(
                egui::Pos2::ZERO,
                egui::vec2(width, 4000.),
            )),
            focused: true,
            events,
            ..Default::default()
        },
        |ui| show(ui),
    );
    let texts = output
        .shapes
        .iter()
        .filter_map(|s| match &s.shape {
            egui::Shape::Text(t) => Some((
                t.galley.text().to_owned(),
                egui::Rect::from_min_size(t.pos, t.galley.size()),
            )),
            _ => None,
        })
        .collect();
    output.drop_without_applying_deltas();
    Painted { texts }
}

/// Enough frames for revealed sections to finish animating open.
pub(crate) fn settle(
    ctx: &egui::Context,
    width: f32,
    show: &mut dyn FnMut(&mut egui::Ui),
) -> Painted {
    for _ in 0..10 {
        frame(ctx, width, vec![], show);
    }
    frame(ctx, width, vec![], show)
}

pub(crate) fn click(pos: egui::Pos2) -> Vec<egui::Event> {
    [true, false]
        .into_iter()
        .flat_map(|pressed| {
            [
                egui::Event::PointerMoved(pos),
                egui::Event::PointerButton {
                    pos,
                    button: egui::PointerButton::Primary,
                    pressed,
                    modifiers: Default::default(),
                },
            ]
        })
        .collect()
}

pub(crate) fn replace_text(text: &str) -> Vec<egui::Event> {
    vec![
        egui::Event::Key {
            key: egui::Key::A,
            physical_key: None,
            pressed: true,
            repeat: false,
            modifiers: egui::Modifiers::COMMAND,
        },
        egui::Event::Text(text.into()),
    ]
}

/// Click the painted callout labelled `label` (asserting it is drawn whole
/// inside the Inspector width), then type `text` into whatever field the
/// click focused. Returns the frame painted afterwards.
pub(crate) fn edit_through_callout(
    ctx: &egui::Context,
    width: f32,
    label: &str,
    text: &str,
    show: &mut dyn FnMut(&mut egui::Ui),
) -> Painted {
    let painted = settle(ctx, width, show);
    let rect = painted.find(label).unwrap_or_else(|| {
        panic!(
            "callout {label:?} is painted at width {width}; texts: {:?}",
            painted.texts.iter().map(|(t, _)| t).collect::<Vec<_>>()
        )
    });
    assert!(
        rect.left() >= 0. && rect.right() <= width,
        "{label:?} is unclipped at width {width}: {rect:?}"
    );
    frame(ctx, width, click(rect.center()), show);
    frame(ctx, width, vec![], show);
    frame(ctx, width, replace_text(text), show);
    settle(ctx, width, show)
}

/// The two Inspector widths every linked form is checked at.
pub(crate) const WIDTHS: [f32; 2] = [200., 420.];
