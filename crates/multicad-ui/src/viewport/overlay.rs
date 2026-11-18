//! The annotation layer over the 3D viewport. The scene is a GPU paint
//! callback added to the viewport's own layer late in the frame, so 2D
//! labels and leaders painted on that layer before it end up underneath the
//! shaded model. Placement tools paint on this sublayer instead: egui draws it
//! directly after the viewport layer, whatever order the shapes were added in.
//! Only the 2D annotations move up; their anchors stay the projected 3D points.
use crate::egui;

/// A painter clipped to `rect` whose shapes are drawn above the viewport
/// layer of `ui`, including the scene callback added later this frame.
pub fn overlay_painter(ui: &egui::Ui, rect: egui::Rect) -> egui::Painter {
    let parent = ui.layer_id();
    let layer = egui::LayerId::new(parent.order, ui.id().with("viewport-annotations"));
    // egui forgets sublayers at the end of every frame.
    ui.ctx().set_sublayer(parent, layer);
    ui.painter()
        .clone()
        .with_layer_id(layer)
        .with_clip_rect(rect)
}

#[cfg(test)]
mod tests {
    use super::*;
    use egui::epaint::{Shape, StrokeKind};

    fn color(shape: &Shape) -> Option<egui::Color32> {
        match shape {
            Shape::Rect(r) => Some(r.fill),
            _ => None,
        }
    }

    #[test]
    fn annotations_paint_above_a_scene_added_later_on_the_viewport_layer() {
        let ctx = egui::Context::default();
        let label = egui::Color32::from_rgb(255, 230, 90);
        let scene = egui::Color32::from_rgb(40, 40, 40);
        let mut order = vec![];
        // Two frames: the layer order must hold every frame, not only the first.
        for _ in 0..2 {
            let output = ctx.run_ui(Default::default(), |ui| {
                let rect = egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(200., 200.));
                let overlay = overlay_painter(ui, rect);
                overlay.rect(
                    rect.shrink(80.),
                    0.,
                    label,
                    egui::Stroke::NONE,
                    StrokeKind::Inside,
                );
                // The shaded model arrives afterwards on the viewport layer.
                ui.painter()
                    .rect(rect, 0., scene, egui::Stroke::NONE, StrokeKind::Inside);
            });
            order = output
                .shapes
                .iter()
                .filter_map(|c| color(&c.shape))
                .filter(|c| [label, scene].contains(c))
                .collect::<Vec<_>>();
            output.drop_without_applying_deltas();
        }
        assert_eq!(order, vec![scene, label]);
    }

    #[test]
    fn annotations_keep_the_viewport_clip() {
        let ctx = egui::Context::default();
        let rect = egui::Rect::from_min_size(egui::pos2(10., 20.), egui::vec2(50., 40.));
        ctx.run_ui(Default::default(), |ui| {
            let overlay = overlay_painter(ui, rect);
            assert_eq!(overlay.clip_rect(), rect.intersect(ui.clip_rect()));
            assert_ne!(overlay.layer_id(), ui.layer_id());
            assert_eq!(overlay.layer_id().order, ui.layer_id().order);
        })
        .drop_without_applying_deltas();
    }
}
