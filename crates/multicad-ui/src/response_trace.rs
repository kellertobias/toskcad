//! Frame and input records for the opt-in response trace (CAD-347). The
//! writer and stage spans live in `multicad_core::response_trace`.
use eframe::egui;
use multicad_core::response_trace as trace;
use std::sync::{
    Mutex,
    atomic::{AtomicU64, Ordering},
};

static FRAME: AtomicU64 = AtomicU64::new(0);
static SEQUENCE: AtomicU64 = AtomicU64::new(0);
static SHOWN: Mutex<Option<(String, String)>> = Mutex::new(None);

/// Event kind as recorded, or `None` for events that are not operator input.
fn kind(event: &egui::Event) -> Option<(&'static str, String)> {
    Some(match event {
        egui::Event::PointerButton {
            button, pressed, ..
        } => (
            if *pressed {
                "pointer-down"
            } else {
                "pointer-up"
            },
            format!("{button:?}"),
        ),
        egui::Event::Key { key, pressed, .. } => (
            if *pressed { "key-down" } else { "key-up" },
            key.name().to_owned(),
        ),
        egui::Event::Text(text) => ("text", text.chars().count().to_string()),
        egui::Event::Paste(_) => ("paste", String::new()),
        egui::Event::MouseWheel { .. } => ("scroll", String::new()),
        egui::Event::Zoom(_) => ("zoom", String::new()),
        // Assistive technology and accessibility automation act without a
        // pointer event; they are still operator actions.
        egui::Event::AccessKitActionRequest(request) => {
            ("accessibility-action", format!("{:?}", request.action))
        }
        _ => return None,
    })
}

/// Start a frame: number it and record every operator input it received.
/// Pointer motion is counted, not listed, because hover alone is not an action.
pub(crate) fn begin_frame(ctx: &egui::Context) -> u64 {
    let frame = FRAME.fetch_add(1, Ordering::Relaxed) + 1;
    if !trace::enabled() {
        return frame;
    }
    let (events, moves, time) = ctx.input(|input| {
        let moves = input
            .events
            .iter()
            .filter(|event| matches!(event, egui::Event::PointerMoved(_)))
            .count();
        (
            input.events.iter().filter_map(kind).collect::<Vec<_>>(),
            moves,
            input.time,
        )
    });
    trace::mark("frame-begin", |record| {
        record.insert("frame".into(), frame.into());
        record.insert("egui_time".into(), time.into());
        record.insert("pointer_moves".into(), moves.into());
    });
    for (kind, detail) in events {
        let sequence = SEQUENCE.fetch_add(1, Ordering::Relaxed) + 1;
        trace::mark("input", |record| {
            record.insert("sequence".into(), sequence.into());
            record.insert("frame".into(), frame.into());
            record.insert("kind".into(), kind.into());
            if !detail.is_empty() {
                record.insert("detail".into(), detail.into());
            }
        });
    }
    frame
}

/// Finish building a frame. The mode and status line are recorded whenever
/// they change so a trace can be read against what the operator saw.
pub(crate) fn end_frame(ctx: &egui::Context, frame: u64, mode: &str, status: &str) {
    if !trace::enabled() {
        return;
    }
    let shown = (mode.to_owned(), status.to_owned());
    let changed = SHOWN
        .lock()
        .map(|mut last| {
            let changed = last.as_ref() != Some(&shown);
            *last = Some(shown);
            changed
        })
        .unwrap_or(false);
    trace::mark("frame-built", |record| {
        record.insert("frame".into(), frame.into());
        if changed {
            record.insert("mode".into(), mode.into());
            record.insert("status".into(), status.into());
        }
    });
    crate::viewport::paint_frame_marker(ctx, frame);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inputs_are_operator_actions_and_motion_is_not_an_action() {
        let modifiers = egui::Modifiers::NONE;
        let down = egui::Event::PointerButton {
            pos: egui::Pos2::ZERO,
            button: egui::PointerButton::Primary,
            pressed: true,
            modifiers,
        };
        assert_eq!(kind(&down).unwrap(), ("pointer-down", "Primary".into()));
        let key = egui::Event::Key {
            key: egui::Key::Enter,
            physical_key: None,
            pressed: true,
            repeat: false,
            modifiers,
        };
        assert_eq!(kind(&key).unwrap(), ("key-down", "Enter".into()));
        assert_eq!(
            kind(&egui::Event::Text("12".into())).unwrap(),
            ("text", "2".into())
        );
        assert!(kind(&egui::Event::PointerMoved(egui::Pos2::ZERO)).is_none());
        assert!(kind(&egui::Event::WindowFocused(true)).is_none());
    }
}
