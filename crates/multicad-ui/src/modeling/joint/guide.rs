//! The corner guide for a refused corner joint: the viewport outlines the
//! corner block both panels must fill and the material each panel lacks,
//! and the Inspector draws the corner's cross-section beside the exact
//! change. Nothing here edits a panel; the user resizes or moves it with
//! the ordinary panel controls and the joint previews once they overlap.
use super::{JointPanel, Kind};
use crate::{
    egui,
    modeling::{TerminationViewport, panel_face},
};
use multicad_core::{
    document::Document,
    model::{self, CornerGuide, PairRole, PanelJointPair},
    parameters::EvaluatedParameters,
};

/// The corner block outline.
const BLOCK: egui::Color32 = egui::Color32::from_rgb(255, 230, 90);
/// Material a panel lacks.
const MISSING: egui::Color32 = egui::Color32::from_rgb(80, 220, 140);

fn role_color(role: PairRole) -> egui::Color32 {
    match role {
        PairRole::A => panel_face::TARGET,
        PairRole::B => panel_face::ROLE_B,
    }
}
fn mm(v: f64) -> String {
    let s = format!("{v:.1}");
    s.trim_end_matches('0').trim_end_matches('.').to_owned()
}

impl JointPanel {
    /// The corner diagnosis of every complete pair of a corner joint: None
    /// for a valid or incomplete pair.
    pub(super) fn diagnoses(&self, document: &Document) -> Vec<Option<CornerGuide>> {
        let corner = !matches!(self.kind, Kind::LooseTenon | Kind::Dowels);
        let Ok(parameters) = multicad_core::parameters::evaluate(document) else {
            return vec![];
        };
        self.pairs
            .iter()
            .map(|pair| match pair {
                [Some(a), Some(b)] if corner && a != b => model::diagnose_panel_corner(
                    document,
                    &parameters,
                    &PanelJointPair { a: *a, b: *b },
                )
                .ok()
                .flatten(),
                _ => None,
            })
            .collect()
    }
}

/// Outline the corner block and every panel's missing material, or say why
/// the pair cannot form a corner at all.
pub(super) fn paint(
    document: &Document,
    parameters: &EvaluatedParameters,
    view: &TerminationViewport<'_>,
    pair: &PanelJointPair,
) {
    if pair.a == pair.b {
        return;
    }
    let guide = match model::diagnose_panel_corner(document, parameters, pair) {
        Ok(Some(guide)) => guide,
        Ok(None) => return,
        Err(error) => {
            if let Ok(b) = model::machining_panel_box(document, parameters, pair.a) {
                let top = panel_face::box_edges(
                    view,
                    &panel_face::panel_box(&b),
                    1.,
                    panel_face::INVALID,
                );
                panel_face::label(
                    view,
                    top + egui::vec2(6., 22.),
                    &error.to_string(),
                    panel_face::INVALID,
                );
            }
            return;
        }
    };
    if let Some(block) = &guide.block {
        let top = panel_face::box_edges(view, block, 1.5, BLOCK);
        panel_face::label(
            view,
            top + egui::vec2(6., 18.),
            "Corner block both panels must fill",
            BLOCK,
        );
    }
    for x in &guide.extensions {
        let top = panel_face::box_edges(view, &x.region, 3., MISSING);
        panel_face::label(
            view,
            top + egui::vec2(6., -6.),
            &format!(
                "+{} mm {} of {}",
                mm(x.amount_mm),
                x.dimension(),
                x.role.label()
            ),
            MISSING,
        );
    }
    if guide.extensions.is_empty() {
        for part in [pair.a, pair.b] {
            if let Ok(b) = model::machining_panel_box(document, parameters, part) {
                panel_face::box_edges(view, &panel_face::panel_box(&b), 1., panel_face::INVALID);
            }
        }
    }
}

/// The Inspector's explanation: the corner cross-section with both panels,
/// the corner block and the missing material, then the exact change.
pub(super) fn section(ui: &mut egui::Ui, guide: &CornerGuide) {
    let s = &guide.section;
    // A window around the corner block, wide enough to show which way each
    // panel continues.
    let pad = (s.block[1][0] - s.block[0][0])
        .max(s.block[1][1] - s.block[0][1])
        .max(1.)
        * 2.5;
    let mut lo = [s.block[0][0] - pad, s.block[0][1] - pad];
    let mut hi = [s.block[1][0] + pad, s.block[1][1] + pad];
    for r in s.add.iter().map(|(_, r)| r) {
        for i in 0..2 {
            lo[i] = lo[i].min(r[0][i] - 2.);
            hi[i] = hi[i].max(r[1][i] + 2.);
        }
    }
    let width = ui.available_width().clamp(160., 260.);
    let (rect, response) =
        ui.allocate_exact_size(egui::vec2(width, width * 0.7), egui::Sense::hover());
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 4., ui.visuals().extreme_bg_color);
    let inner = rect.shrink(10.);
    let scale =
        (inner.width() as f64 / (hi[0] - lo[0])).min(inner.height() as f64 / (hi[1] - lo[1]));
    let to = |p: [f64; 2]| {
        let p = [p[0].clamp(lo[0], hi[0]), p[1].clamp(lo[1], hi[1])];
        egui::pos2(
            inner.center().x + ((p[0] - (lo[0] + hi[0]) / 2.) * scale) as f32,
            inner.center().y - ((p[1] - (lo[1] + hi[1]) / 2.) * scale) as f32,
        )
    };
    let rect_of = |r: &[[f64; 2]; 2]| egui::Rect::from_two_pos(to(r[0]), to(r[1]));
    for (r, role) in [(&s.a, PairRole::A), (&s.b, PairRole::B)] {
        let c = role_color(role);
        painter.rect(
            rect_of(r),
            0.,
            c.gamma_multiply(0.35),
            egui::Stroke::new(1.5, c),
            egui::StrokeKind::Inside,
        );
        let r = rect_of(r);
        painter.text(
            r.center(),
            egui::Align2::CENTER_CENTER,
            role.label(),
            egui::FontId::proportional(13.),
            c,
        );
    }
    for (role, r) in &s.add {
        painter.rect(
            rect_of(r),
            0.,
            MISSING.gamma_multiply(0.45),
            egui::Stroke::new(1.5, MISSING),
            egui::StrokeKind::Inside,
        );
        painter.text(
            rect_of(r).center(),
            egui::Align2::CENTER_CENTER,
            format!("+{}", role.label()),
            egui::FontId::proportional(11.),
            egui::Color32::BLACK,
        );
    }
    if guide.block.is_some() {
        painter.rect_stroke(
            rect_of(&s.block),
            0.,
            egui::Stroke::new(1.5, BLOCK),
            egui::StrokeKind::Middle,
        );
    }
    response.on_hover_text(
        "Corner cross-section: A and B as they are, the yellow corner block both must fill, green material each one lacks",
    );
    ui.colored_label(egui::Color32::LIGHT_RED, &guide.message);
}
