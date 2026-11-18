//! Point-to-point Move: a source snap point on the moving targets and a target
//! snap point on other geometry define one staged translation. Points are
//! display-derived picks; only the resulting numeric translation is authored.
use super::super::gizmo::math::{Point, Values, sub, translated};
use crate::{camera::Camera, egui, viewport::Scene};
use multicad_core::{
    document::{Document, EntityData, EntityId},
    model::NumericTransform,
    parameters::EvaluatedParameters,
};
use std::collections::BTreeSet;

/// Screen distance within which a pointer acquires a snap point.
pub(crate) const SNAP_PIXELS: f32 = 10.;

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum PointPick {
    /// Waiting for a point on the moving targets.
    Source,
    /// Waiting for a point on other geometry; `source` is pre-placement world.
    Target { source: Point, label: &'static str },
}

/// The last applied pair. The source stays in pre-placement world
/// coordinates so its marker follows later staged edits of the placement.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct PointAlignment {
    pub source: Point,
    pub source_label: &'static str,
    pub target: Point,
    pub target_label: &'static str,
}

/// One supported snap point: a topology vertex or straight-edge midpoint.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct SnapPoint {
    /// Instance owning the point; `None` for the world origin.
    pub owner: Option<EntityId>,
    /// World point before the staged placement.
    pub canonical: Point,
    /// Where the point is displayed (after the staged placement for sources).
    pub shown: Point,
    pub label: &'static str,
}

/// Instances moved by the targets: the targets themselves and linked
/// Instances of a targeted Part definition, which inherit its movement.
pub(crate) fn moving_instances(
    document: &Document,
    scene: &Scene,
    selected: &BTreeSet<EntityId>,
) -> BTreeSet<EntityId> {
    scene
        .instances()
        .iter()
        .map(|instance| instance.id)
        .filter(|id| {
            selected.contains(id)
                || matches!(
                    document.get(*id).map(|e| &e.data),
                    Some(EntityData::Instance { component, .. }) if selected.contains(component)
                )
        })
        .collect()
}

/// Supported snap points of the scene instances accepted by `owner`. `shown`
/// maps a pre-placement world point to its displayed position.
pub(crate) fn snap_points(
    scene: &Scene,
    owner: impl Fn(EntityId) -> bool,
    shown: impl Fn(Point) -> Point,
) -> Vec<SnapPoint> {
    let mut points = vec![];
    for instance in scene.instances().iter().filter(|i| owner(i.id)) {
        let Some(geometry) = scene.producer_geometry(instance.geometry) else {
            continue;
        };
        let mut push = |local: Point, label| {
            let canonical = instance.transform.point(local);
            points.push(SnapPoint {
                owner: Some(instance.id),
                canonical,
                shown: shown(canonical),
                label,
            });
        };
        for vertex in geometry.points() {
            push(vertex.position_mm, "Vertex");
        }
        for edge in geometry.edges() {
            if let [a, b] = edge.points_mm.as_slice() {
                push(
                    std::array::from_fn(|i| a[i] + (b[i] - a[i]) * 0.5),
                    "Edge midpoint",
                );
            }
        }
    }
    points
}

/// The nearest snap point to the pointer on screen, preferring the nearer
/// one in depth among points that project to the same place.
pub(crate) fn nearest(
    camera: Camera,
    rect: egui::Rect,
    pointer: egui::Pos2,
    points: &[SnapPoint],
) -> Option<(SnapPoint, egui::Pos2)> {
    points
        .iter()
        .filter_map(|point| {
            let [x, y, depth] = camera.project(point.shown.map(|v| v as f32), rect);
            let screen = egui::pos2(x, y);
            let distance = screen.distance(pointer);
            (screen.is_finite()
                && (!camera.perspective || depth > 0.1)
                && rect.contains(screen)
                && distance <= SNAP_PIXELS)
                .then_some((*point, screen, distance, depth))
        })
        .min_by(|a, b| {
            // Within one pixel, the point closest to the viewer wins.
            if (a.2 - b.2).abs() < 1. {
                a.3.total_cmp(&b.3)
            } else {
                a.2.total_cmp(&b.2)
            }
        })
        .map(|(point, screen, ..)| (point, screen))
}

/// Translate the staged placement so the placed `source` coincides with
/// `target`. Rotation, axis and pivot are kept, so dimensions and the
/// orientation chosen so far are unchanged.
pub(crate) fn aligned(
    parameters: &NumericTransform,
    evaluated: &EvaluatedParameters,
    source: Point,
    target: Point,
) -> Result<NumericTransform, String> {
    if source.iter().chain(&target).any(|v| !v.is_finite()) {
        return Err("The chosen point is not finite".into());
    }
    let values = Values::read(parameters, evaluated)?;
    let pose = parameters
        .resolve(evaluated)
        .map_err(|error| error.to_string())?;
    let delta = sub(target, pose.point(source));
    let result = translated(parameters, values, delta)?;
    let placed = result
        .resolve(evaluated)
        .map_err(|error| error.to_string())?
        .point(source);
    if sub(placed, target).iter().any(|v| v.abs() > 1e-6) {
        return Err("The point-to-point translation could not be represented exactly".into());
    }
    Ok(result)
}

/// Paint the source (square) and target (ring) markers with the pending or
/// applied translation between them.
pub(crate) fn paint_markers(
    painter: &egui::Painter,
    camera: Camera,
    rect: egui::Rect,
    source: Option<(Point, &str)>,
    target: Option<(Point, &str)>,
) {
    let screen = |p: Point| {
        let [x, y, depth] = camera.project(p.map(|v| v as f32), rect);
        let screen = egui::pos2(x, y);
        (screen.is_finite() && (!camera.perspective || depth > 0.1)).then_some(screen)
    };
    let source_color = egui::Color32::from_rgb(255, 170, 51);
    let target_color = egui::Color32::from_rgb(92, 214, 120);
    let source = source.and_then(|(p, label)| Some((screen(p)?, label)));
    let target = target.and_then(|(p, label)| Some((screen(p)?, label)));
    if let (Some((a, _)), Some((b, _))) = (source, target)
        && a.distance(b) > 2.
    {
        let stroke = egui::Stroke::new(1.5, egui::Color32::from_gray(230));
        painter.extend(egui::Shape::dashed_line(&[a, b], stroke, 6., 4.));
    }
    if let Some((p, label)) = source {
        painter.rect_filled(
            egui::Rect::from_center_size(p, egui::vec2(9., 9.)),
            1.,
            source_color,
        );
        painter.text(
            p + egui::vec2(-10., -9.),
            egui::Align2::RIGHT_BOTTOM,
            format!("Source · {label}"),
            egui::FontId::proportional(12.),
            source_color,
        );
    }
    if let Some((p, label)) = target {
        painter.circle_stroke(p, 8., egui::Stroke::new(2., target_color));
        painter.circle_filled(p, 2.5, target_color);
        painter.text(
            p + egui::vec2(11., 9.),
            egui::Align2::LEFT_TOP,
            format!("Target · {label}"),
            egui::FontId::proportional(12.),
            target_color,
        );
    }
}
