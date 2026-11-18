//! Visual pairing for the joint form: click panel A, then panel B, on their
//! visible faces. The form keeps one waiting operand: a click fills it and
//! moves on to the next empty one, Pick in the Inspector chooses the
//! operand to replace, and once every pair is complete further clicks are
//! explained instead of silently starting another pair. Both chosen panels
//! stay outlined in their role colors. For corner joints the viewport shows
//! the missing material of a refused corner; for loose tenons it draws both
//! receiving panels, the insertion direction, every rounded mortise and the
//! tenons resolved from the current values, or the reason they are refused.
use super::{JointPanel, Kind, guide};
use crate::{
    egui,
    modeling::{TerminationViewport, panel_face},
};
use multicad_core::{
    document::{Document, EntityId},
    model::{self, MachiningFace},
    parameters::EvaluatedParameters,
};

const ROLES: [&str; 2] = ["A", "B"];

impl JointPanel {
    pub fn picking(&self) -> bool {
        self.place
    }
    /// Paint the pairs, track the hovered panel and apply a click. Returns
    /// true when a click changed the pairs.
    pub fn viewport(
        &mut self,
        document: &Document,
        revision: u64,
        view: TerminationViewport<'_>,
    ) -> bool {
        let Ok(parameters) = multicad_core::parameters::evaluate(document) else {
            return false;
        };
        self.paint(document, &parameters, &view);
        if !self.place {
            return false;
        }
        let Some(hit) = panel_face::hover(document, revision, &view) else {
            // Explain an object that is not a rectangular panel.
            if let Some(part) = panel_face::hovered_part(document, revision, &view) {
                let reason = refusal(document, part);
                if let Some(p) = view.pointer {
                    panel_face::label(&view, p + egui::vec2(8., -8.), &reason, panel_face::INVALID);
                }
                if view.clicked {
                    self.pick_note = Some(reason);
                }
            }
            return false;
        };
        let name = document.get(hit.part).map_or("", |e| e.name.as_str());
        let p = panel_face::on_screen(&view, &hit.frame, hit.at);
        let waiting = if self.pairs.is_empty() {
            Some((0, 0))
        } else {
            self.active
        };
        let (text, color) = match waiting {
            _ if self.existing.is_some() => (
                "An applied joint keeps its panel pairs".to_owned(),
                panel_face::INVALID,
            ),
            None => (
                format!("{name} · every pair is chosen; use Pick in the Inspector to replace one"),
                panel_face::INVALID,
            ),
            Some((i, k)) => match self.duplicate(i, k, hit.part) {
                Some(note) => (note, panel_face::INVALID),
                None => (
                    format!(
                        "Set pair {} panel {} · {name} · {}",
                        i + 1,
                        ROLES[k],
                        panel_face::face_label(hit.face)
                    ),
                    panel_face::HOVER,
                ),
            },
        };
        panel_face::outline(&view, &hit.frame, color);
        panel_face::label(&view, p + egui::vec2(8., -8.), &text, color);
        if !view.clicked {
            return false;
        }
        let changed = self.pick(document, hit.part);
        if changed {
            self.stamp.changed();
            self.error = None;
        }
        changed
    }
    /// The first operand without a panel.
    pub(super) fn first_open_slot(&self) -> Option<(usize, usize)> {
        self.next_open_slot(0)
    }
    /// The first operand without a panel, starting at pair `from` and
    /// wrapping around.
    pub(super) fn next_open_slot(&self, from: usize) -> Option<(usize, usize)> {
        let n = self.pairs.len();
        (0..n).map(|d| (from + d) % n).find_map(|i| {
            self.pairs[i]
                .iter()
                .position(Option::is_none)
                .map(|k| (i, k))
        })
    }
    /// Remove pair `i`; the waiting operand stays on its pair, or moves to
    /// the first empty one when its own pair is removed.
    pub(super) fn remove_pair(&mut self, i: usize) {
        self.pairs.remove(i);
        self.active = match self.active {
            Some((j, _)) if j == i => self.first_open_slot(),
            Some((j, k)) if j > i => Some((j - 1, k)),
            other => other,
        };
    }
    /// Why `part` cannot fill operand (i, k), if it cannot.
    fn duplicate(&self, i: usize, k: usize, part: EntityId) -> Option<String> {
        let pair = self.pairs.get(i)?;
        if pair[1 - k] == Some(part) {
            return Some(format!(
                "Already panel {} of pair {}; choose a different panel for {}",
                ROLES[1 - k],
                i + 1,
                ROLES[k]
            ));
        }
        (pair[k] == Some(part)).then(|| format!("Already panel {} of pair {}", ROLES[k], i + 1))
    }
    /// Put a clicked panel into the waiting operand and wait for the next
    /// empty one.
    pub(super) fn pick(&mut self, document: &Document, part: EntityId) -> bool {
        self.pick_note = None;
        let name = document.get(part).map_or(String::new(), |e| e.name.clone());
        if self.existing.is_some() {
            self.pick_note = Some(
                "An applied joint keeps its panel pairs; delete it and create another to change them"
                    .into(),
            );
            return false;
        }
        let (panels, refused) = panel_face::selected_panels(document, [part]);
        if panels.is_empty() {
            self.pick_note = refused.into_iter().next();
            return false;
        }
        if self.pairs.is_empty() {
            self.pairs.push([None, None]);
            self.active = Some((0, 0));
        }
        let Some((i, k)) = self.active.filter(|(i, _)| *i < self.pairs.len()) else {
            self.pick_note = Some(format!(
                "{name} was not used: every pair is chosen. Use Pick beside A or B to replace one, or Add pair for another corner"
            ));
            return false;
        };
        if let Some(note) = self.duplicate(i, k, part) {
            self.pick_note = Some(format!("{name}: {note}"));
            return false;
        }
        self.pairs[i][k] = Some(part);
        self.active = self.next_open_slot(i);
        true
    }
    fn paint(
        &self,
        document: &Document,
        parameters: &EvaluatedParameters,
        view: &TerminationViewport<'_>,
    ) {
        let spec = match (self.kind, self.tenon.spec(document)) {
            (Kind::LooseTenon, Ok(spec)) => Some(spec),
            _ => None,
        };
        let corner = !matches!(self.kind, Kind::LooseTenon | Kind::Dowels);
        for (i, pair) in self.pairs.iter().enumerate() {
            for (k, part) in pair.iter().enumerate() {
                let Some(part) = part else { continue };
                if let Ok(b) = model::machining_panel_box(document, parameters, *part) {
                    let name = document
                        .get(*part)
                        .map_or(String::new(), |e| e.name.clone());
                    let color = [panel_face::TARGET, panel_face::ROLE_B][k];
                    let top = panel_face::box_edges(view, &panel_face::panel_box(&b), 2.5, color);
                    let at = top + egui::vec2(6., -6. - 14. * k as f32);
                    panel_face::label(
                        view,
                        at,
                        &format!("Pair {} · {} · {name}", i + 1, ROLES[k]),
                        color,
                    );
                }
            }
            if corner && let [Some(a), Some(b)] = pair {
                guide::paint(
                    document,
                    parameters,
                    view,
                    &model::PanelJointPair { a: *a, b: *b },
                );
            }
            let (Some(spec), [Some(a), Some(b)]) = (&spec, pair) else {
                continue;
            };
            let pair = model::PanelJointPair { a: *a, b: *b };
            match model::resolve_loose_tenon(document, parameters, spec, &pair) {
                Ok((_, layout)) => {
                    for m in &layout.face_mortises {
                        panel_face::polygon(view, &layout.face_frame, &m.outline, panel_face::CUT);
                    }
                    for m in &layout.edge_mortises {
                        panel_face::polygon(view, &layout.edge_frame, &m.outline, panel_face::CUT);
                    }
                    for t in &layout.tenons {
                        let at = |k: f64| {
                            let p = [0, 1, 2].map(|i| t.center[i] + t.axis[i] * k);
                            let s = view
                                .camera
                                .project([p[0] as f32, p[1] as f32, p[2] as f32], view.rect);
                            egui::pos2(s[0], s[1])
                        };
                        // The insertion direction: face panel into edge panel.
                        view.painter.arrow(
                            at(-t.into_face_mm),
                            at(t.into_edge_mm) - at(-t.into_face_mm),
                            egui::Stroke::new(2., panel_face::HOVER),
                        );
                    }
                }
                Err(error) => {
                    if let Some(frame) =
                        panel_face::frame(document, parameters, *a, MachiningFace::Top)
                    {
                        panel_face::outline(view, &frame, panel_face::INVALID);
                        panel_face::label(
                            view,
                            panel_face::on_screen(view, &frame, [0., 0.]),
                            &error.to_string(),
                            panel_face::INVALID,
                        );
                    }
                }
            }
        }
    }
}

/// Why a clicked Part cannot be a joint panel.
fn refusal(document: &Document, part: EntityId) -> String {
    let (_, refused) = panel_face::selected_panels(document, [part]);
    refused
        .into_iter()
        .next()
        .unwrap_or_else(|| "Choose a rectangular panel".into())
}
