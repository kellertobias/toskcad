//! CAD-417: leaving Project with unapplied parameter rows asks first.
//!
//! Project keeps explicit Apply and Rebuild and Discard Edit actions, so a
//! pending parameter draft is never applied silently. Leaving Project by a
//! mode switch, another Project tab, another project or Close Project asks to
//! apply, discard or stay. Closing the window is exempt (user decision
//! 2026-11-27) and goes straight to the Quit dialog. Measure is a temporary
//! overlay that keeps the draft, so it does not ask.
use super::*;
#[cfg(test)]
#[path = "parameter_leave_tests.rs"]
mod tests;

/// Where the user was going when the prompt interrupted them.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum ParameterLeave {
    /// A command that leaves the parameter table (mode switch or Close
    /// Project).
    Command(Command),
    /// Another Project tab.
    Tab(usize),
    /// Another open project.
    Project(usize),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum LeaveChoice {
    Apply,
    Discard,
    Stay,
}

pub(super) const PROMPT_TITLE: &str = "Unapplied parameter edits";
pub(super) const APPLY_LABEL: &str = "Apply and continue";
pub(super) const DISCARD_LABEL: &str = "Discard changes";
pub(super) const STAY_LABEL: &str = "Stay";
/// Accent of the pending-edit dot on the Project mode tab.
pub(super) const PENDING_MARKER: Color32 = Color32::from_rgb(255, 196, 64);

impl Application {
    /// The active project shows its parameter table with edits that wait
    /// for Apply and Rebuild.
    pub(super) fn parameter_rows_pending(&self) -> bool {
        !self.welcome && self.mode == Mode::Parameters && self.current().editor.rows_dirty()
    }

    /// Only the parameter draft (and file details, which commit on leave)
    /// is pending, so the prompt alone decides what leaving does.
    pub(super) fn only_parameter_rows_pending(&self) -> bool {
        let project = self.current();
        self.parameter_rows_pending()
            && !project.sketch.has_pending_edits()
            && !project.drawings.dirty()
            && !project.exchange.dirty()
            && !project.plate.dirty()
            && !project.sketch.material_library().dirty()
    }

    /// The Project mode tab carries a pending marker while this holds.
    pub(super) fn project_tab_pending(&self) -> bool {
        !self.welcome && self.current().editor.rows_dirty()
    }

    /// Ask before `target` leaves the parameter draft behind. Returns
    /// whether the prompt now owns the request.
    pub(super) fn guard_parameter_leave(&mut self, target: ParameterLeave) -> bool {
        if !self.only_parameter_rows_pending() {
            return false;
        }
        if self.parameter_leave.is_none() {
            self.parameter_leave = Some(target);
        }
        self.message =
            "Unapplied parameter edits · Apply and continue, Discard changes or Stay".into();
        self.message_error = false;
        true
    }

    /// The UI commands that leave the parameter table.
    pub(super) fn parameter_leave_target(&self, command: Command) -> Option<ParameterLeave> {
        let leaves = match command {
            Command::Mode(mode) => mode != Mode::Parameters,
            Command::CloseProject => true,
            _ => false,
        };
        leaves.then_some(ParameterLeave::Command(command))
    }

    /// Carry out the user's choice; Apply and Discard continue to the
    /// interrupted destination, Stay keeps the table open.
    pub(super) fn resolve_parameter_leave(&mut self, choice: LeaveChoice) -> Result<(), String> {
        let Some(target) = self.parameter_leave.take() else {
            return Ok(());
        };
        match choice {
            LeaveChoice::Stay => {
                if self.mode == Mode::Parameters {
                    self.tab = project::PARAMETERS_TAB;
                }
                self.message = "Parameter edits are still pending".into();
                self.message_error = false;
                return Ok(());
            }
            LeaveChoice::Apply => {
                if let Err(error) = self.invoke(Command::Apply, &[], Channel::Ui) {
                    self.parameter_leave = Some(target);
                    return Err(error);
                }
            }
            LeaveChoice::Discard => self.invoke(Command::Discard, &[], Channel::Ui)?,
        }
        if self.current().editor.rows_dirty() {
            self.parameter_leave = Some(target);
            return Err("The parameter draft is still pending".into());
        }
        match target {
            ParameterLeave::Command(command) => self.invoke(command, &[], Channel::Ui),
            ParameterLeave::Tab(index) => {
                self.select_project_tab(index);
                Ok(())
            }
            ParameterLeave::Project(index) => {
                if index < self.projects.len() {
                    self.activate_project(index);
                }
                Ok(())
            }
        }
    }

    pub(super) fn show_parameter_leave_dialog(&mut self, ctx: &egui::Context) {
        if self.parameter_leave.is_none() {
            return;
        }
        if self.welcome || !self.current().editor.rows_dirty() {
            self.parameter_leave = None;
            return;
        }
        let can_apply = self.enabled(Command::Apply);
        let modal = controls::modal(ctx, "parameter-leave", PROMPT_TITLE, |ui| {
            ui.label("The parameter table has edits that are not applied yet.");
            ui.label("Apply them as one undoable step, discard them, or stay in Parameters.");
            ui.add_space(12.);
            ui.horizontal(|ui| {
                let apply = ui
                    .add_enabled(
                        can_apply,
                        egui::Button::new(APPLY_LABEL).fill(ui.visuals().selection.bg_fill),
                    )
                    .on_hover_text("Apply and Rebuild, then continue")
                    .on_disabled_hover_text("Resolve the draft error before applying")
                    .clicked();
                let discard = ui
                    .button(DISCARD_LABEL)
                    .on_hover_text("Restore the committed parameters, then continue")
                    .clicked();
                let stay = ui
                    .button(STAY_LABEL)
                    .on_hover_text("Keep editing the parameter table · Escape")
                    .clicked();
                if apply {
                    Some(LeaveChoice::Apply)
                } else if discard {
                    Some(LeaveChoice::Discard)
                } else if stay {
                    Some(LeaveChoice::Stay)
                } else {
                    None
                }
            })
            .inner
        });
        let choice = if modal.should_close() {
            Some(LeaveChoice::Stay)
        } else {
            modal.inner
        };
        if let Some(choice) = choice {
            if let Err(error) = self.resolve_parameter_leave(choice) {
                self.message = error;
                self.message_error = true;
            }
            ctx.request_repaint();
        }
    }

    /// Paint the pending-edit dot on the Project mode tab.
    pub(super) fn paint_project_pending_marker(&self, ui: &egui::Ui, tab: &egui::Response) {
        let center = tab.rect.right_top() + egui::vec2(-5., 6.);
        ui.painter().circle_filled(center, 3.5, PENDING_MARKER);
    }
}
