//! The Inspector of an active sketch draw tool (CAD-399): the tool's own
//! settings first (Fillet radius, Chamfer distance, Text), then its step
//! guidance and drawing options. Nothing else is shown while a tool is active.
use super::*;

impl SketchEditor {
    pub(super) fn tool_inspector(
        &mut self,
        ui: &mut egui::Ui,
        session: &mut Session,
        tool: DrawTool,
    ) {
        ui.strong(tool.label());
        if tool.corner() {
            self.corner_size_ui(ui, tool);
        }
        if tool == DrawTool::Text {
            self.text_tool_ui(ui, session);
        }
        ui.separator();
        self.tool_guide_ui(ui, tool);
        if !tool.corner() && tool != DrawTool::Text {
            ui.separator();
            ui.label(egui::RichText::new("Options").small().weak());
            if ui
                .checkbox(&mut self.construction, "Construction geometry")
                .on_hover_text(
                    "Construction curves constrain the sketch but are excluded from solid profiles",
                )
                .changed()
            {
                self.message = if self.construction {
                    "Construction drawing on"
                } else {
                    "Construction drawing off"
                }
                .into();
            }
        }
        if tool == DrawTool::Arc {
            ui.checkbox(&mut self.arc_clockwise, "Clockwise arc")
                .on_hover_text("Direction around the sketch's positive normal");
        }
        if tool == DrawTool::Spline {
            ui.small(format!("{} control poles", self.spline_poles.len()));
            if ui
                .add_enabled(
                    !self.spline_poles.is_empty(),
                    egui::Button::new("Remove last pole"),
                )
                .on_hover_text("Remove only the last uncommitted control pole")
                .clicked()
            {
                self.spline_poles.pop();
            }
            use crate::inspector_actions::{self, Button, Owner};
            let clicked = inspector_actions::row(
                ui,
                Owner::Drawing,
                None,
                Button::new(
                    "Finish spline",
                    self.can_finish_spline(),
                    "Commit the pending open spline; needs at least four control poles",
                ),
                Button::new(
                    "Cancel drawing",
                    true,
                    "Discard the uncommitted control poles and end the tool; existing geometry stays",
                ),
            );
            if clicked.apply
                && let Err(error) = self.finish_spline(session, false, Channel::Ui)
            {
                self.message = error;
            }
            if clicked.cancel {
                self.cancel();
            }
        }
        ui.separator();
    }
}
