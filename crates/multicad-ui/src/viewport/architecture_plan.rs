//! Detached native plan graphics, derived on the geometry worker, never CAD source.
use super::Scene;
use crate::{camera::Camera, egui};
use multicad_core::{
    document::EntityId,
    drawings::{DrawingVector, VectorRole},
    project_file::ContentHash,
};
use std::{collections::BTreeSet, sync::Arc};

#[derive(Clone)]
pub struct ArchitecturePlan {
    pub source_hash: ContentHash,
    pub sources: BTreeSet<EntityId>,
    pub vectors: Vec<DrawingVector>,
}
impl ArchitecturePlan {
    pub(super) fn validate(
        &self,
        source: ContentHash,
        roots: &BTreeSet<EntityId>,
    ) -> Result<(), String> {
        if self.source_hash != source || !self.sources.is_subset(roots) {
            return Err("Native plan belongs to another source or occurrence set".into());
        }
        if self.vectors.iter().any(|vector| {
            !self.sources.contains(&vector.source)
                || !vector.line_width_mm.is_finite()
                || vector.line_width_mm <= 0.
                || vector.preview_mm.len() < 2
                || vector
                    .preview_mm
                    .iter()
                    .flatten()
                    .any(|value| !value.is_finite())
        }) {
            return Err("Native plan contains invalid or unattributed graphics".into());
        }
        Ok(())
    }
    /// The original scene remains the hit/reference scene. These sources get
    /// their native section/plan graphics instead of an occluding wall top.
    pub fn surface_scene(&self, scene: &Arc<Scene>) -> Result<Arc<Scene>, String> {
        let visible: BTreeSet<_> = scene
            .instances()
            .iter()
            .filter(|instance| !self.sources.contains(&instance.id))
            .map(|instance| instance.id)
            .collect();
        if visible.len() == scene.instances().len() {
            Ok(scene.clone())
        } else {
            scene.visible_subset(&visible)
        }
    }
    /// Emphasize selected wall-hosted openings in the plan: the footprint
    /// each cuts through its wall, from the canonical wall source. Openings
    /// of walls not drawn in this plan are skipped.
    pub fn paint_openings(
        &self,
        painter: &egui::Painter,
        session: &multicad_core::commands::Session,
        camera: Camera,
        rect: egui::Rect,
        openings: &BTreeSet<EntityId>,
    ) {
        let document = session.document();
        for id in openings {
            let Some(host) = multicad_core::model::opening_host(document, *id) else {
                continue;
            };
            if !self.sources.contains(&host.part) {
                continue;
            }
            let Some(corners) =
                multicad_core::model::opening_footprint(document, session.parameters(), *id)
            else {
                continue;
            };
            let mut points: Vec<_> = corners
                .iter()
                .map(|p| {
                    let pixel = camera.project([p[0] as f32, p[1] as f32, camera.target[2]], rect);
                    egui::pos2(pixel[0], pixel[1])
                })
                .collect();
            points.push(points[0]);
            let color = egui::Color32::from_rgb(90, 170, 255);
            painter.add(egui::Shape::line(points, egui::Stroke::new(2.5, color)));
        }
    }
    pub fn paint(
        &self,
        painter: &egui::Painter,
        scene: &Scene,
        camera: Camera,
        rect: egui::Rect,
        selected: &BTreeSet<EntityId>,
        hovered: Option<EntityId>,
    ) {
        let visible: BTreeSet<_> = scene
            .instances()
            .iter()
            .map(|instance| instance.id)
            .collect();
        for vector in &self.vectors {
            if !visible.contains(&vector.source) {
                continue;
            }
            let points: Vec<_> = vector
                .preview_mm
                .iter()
                .map(|point| {
                    let pixel =
                        camera.project([point[0] as f32, point[1] as f32, camera.target[2]], rect);
                    egui::pos2(pixel[0], pixel[1])
                })
                .collect();
            let highlighted = selected.contains(&vector.source) || hovered == Some(vector.source);
            let color = if highlighted {
                egui::Color32::from_rgb(90, 170, 255)
            } else if let Some([r, g, b]) = vector.color_rgb {
                egui::Color32::from_rgb(r, g, b)
            } else {
                match vector.role {
                    VectorRole::Section => egui::Color32::from_gray(235),
                    VectorRole::Hatch => egui::Color32::from_gray(155),
                    VectorRole::Grain => egui::Color32::from_gray(120),
                    VectorRole::Visible => egui::Color32::from_gray(225),
                    VectorRole::Hidden => egui::Color32::from_gray(150),
                }
            };
            let stroke =
                egui::Stroke::new((vector.line_width_mm as f32 * 4.).clamp(0.5, 2.), color);
            if matches!(
                vector.curve,
                multicad_core::drawings::VectorCurve::HatchTriangle { .. }
            ) {
                painter.add(egui::Shape::convex_polygon(
                    points,
                    color,
                    egui::Stroke::NONE,
                ));
            } else if vector.role == VectorRole::Hidden {
                painter.extend(egui::Shape::dashed_line(&points, stroke, 5., 3.));
            } else {
                painter.add(egui::Shape::line(points, stroke));
            }
        }
    }
}
