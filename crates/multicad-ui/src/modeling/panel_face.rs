//! Panel face placement in the viewport: locate the rectangular panel face
//! under the pointer in its own (u, v) coordinates, snap the point, and
//! draw face outlines, local axes, edge offsets, holes and outlines over the
//! scene. Machining and loose-tenon forms use it to place features directly
//! on visible faces; the Inspector keeps the same values.
use super::TerminationViewport;
use crate::{
    egui,
    viewport::{PickFilter, PickOptions},
};
use multicad_core::{
    document::{Document, EntityData, EntityId},
    model::{self, MachiningFace, MachiningFaceFrame, PanelEdgeName},
    parameters::EvaluatedParameters,
};
use multicad_precision::PortableMath as _;

/// The faces of a rectangular panel in hit-test order.
pub(super) const FACES: [MachiningFace; 6] = [
    MachiningFace::Top,
    MachiningFace::Bottom,
    MachiningFace::Edge {
        edge: PanelEdgeName::UStart,
    },
    MachiningFace::Edge {
        edge: PanelEdgeName::UEnd,
    },
    MachiningFace::Edge {
        edge: PanelEdgeName::VStart,
    },
    MachiningFace::Edge {
        edge: PanelEdgeName::VEnd,
    },
];
/// How far from a face plane a picked point may lie.
const ON_FACE_MM: f64 = 0.5;
pub(super) const TARGET: egui::Color32 = egui::Color32::from_rgb(255, 176, 64);
pub(super) const HOVER: egui::Color32 = egui::Color32::from_rgb(255, 230, 90);
pub(super) const CUT: egui::Color32 = egui::Color32::from_rgb(80, 220, 140);
pub(super) const INVALID: egui::Color32 = egui::Color32::from_rgb(255, 90, 90);
/// Panel B of a joint pair; panel A uses `TARGET`.
pub(super) const ROLE_B: egui::Color32 = egui::Color32::from_rgb(90, 190, 255);
pub(super) const AXIS_U: egui::Color32 = egui::Color32::from_rgb(235, 80, 80);
pub(super) const AXIS_V: egui::Color32 = egui::Color32::from_rgb(90, 200, 90);
/// Cosine between view ray and face normal below which a face is edge-on.
const EDGE_ON: f64 = 1e-6;
const LABEL_PLATE: egui::Color32 = egui::Color32::from_rgba_premultiplied(16, 18, 22, 210);

/// A point on a panel face, in that face's coordinates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct FaceHit {
    pub part: EntityId,
    pub face: MachiningFace,
    pub frame: MachiningFaceFrame,
    pub at: [f64; 2],
}

pub(super) fn face_label(face: MachiningFace) -> &'static str {
    match face {
        MachiningFace::Top => "top face",
        MachiningFace::Bottom => "bottom face",
        MachiningFace::Inside => "inside face",
        MachiningFace::Edge { edge } => match edge {
            PanelEdgeName::UStart => "length-start edge",
            PanelEdgeName::UEnd => "length-end edge",
            PanelEdgeName::VStart => "width-start edge",
            PanelEdgeName::VEnd => "width-end edge",
        },
    }
}

/// The Part an entity selection stands for: the Part itself or the Part
/// holding a selected feature.
pub(crate) fn selected_part(document: &Document, id: EntityId) -> Option<EntityId> {
    match document.get(id).map(|e| &e.data) {
        Some(EntityData::Part { .. }) => Some(id),
        Some(EntityData::Feature { .. }) => document.entities().find_map(|e| match &e.data {
            EntityData::Part { features, .. } if features.contains(&id) => Some(e.id),
            _ => None,
        }),
        _ => None,
    }
}
/// The rectangular panel Parts named by a selection, in selection order,
/// and a reason for every selected item that is not one.
pub(crate) fn selected_panels(
    document: &Document,
    selection: impl IntoIterator<Item = EntityId>,
) -> (Vec<EntityId>, Vec<String>) {
    let Ok(parameters) = multicad_core::parameters::evaluate(document) else {
        return (vec![], vec![]);
    };
    let mut panels = vec![];
    let mut refused = vec![];
    for id in selection {
        let name = document.get(id).map_or(String::new(), |e| e.name.clone());
        match selected_part(document, id) {
            Some(part) if model::machining_panel_box(document, &parameters, part).is_ok() => {
                if !panels.contains(&part) {
                    panels.push(part);
                }
            }
            Some(part) if model::panel_feature(document, part).is_some() => refused.push(format!(
                "{name} is a sketch-region panel; choose a rectangular panel"
            )),
            _ => refused.push(format!("{name} is not a cabinet panel")),
        }
    }
    (panels, refused)
}
/// The status line of a tool opened with targets taken from the selection:
/// it names the targets it kept and why any selected item was not used.
pub(crate) fn seeded_message(
    tool: &str,
    document: &Document,
    targets: &[EntityId],
    refused: &[String],
    otherwise: &str,
) -> String {
    let names: Vec<_> = targets
        .iter()
        .map(|t| document.get(*t).map_or(String::new(), |e| e.name.clone()))
        .collect();
    let mut parts = vec![tool.to_owned()];
    if names.is_empty() {
        parts.extend(refused.iter().cloned());
        parts.push(otherwise.into());
    } else {
        parts.push(format!("targets from the selection: {}", names.join(", ")));
        if !refused.is_empty() {
            parts.push(format!("not used: {}", refused.join("; ")));
        }
        parts.push("check them in the viewport and Inspector".into());
    }
    parts.join(" · ")
}
/// Both sides of every carcass whose left or right side is in `panels`, so
/// shelf-pin rows on opposite sides stay aligned.
pub(crate) fn with_opposite_sides(document: &Document, panels: Vec<EntityId>) -> Vec<EntityId> {
    let mut out = panels.clone();
    for part in panels {
        let Some(model::PanelSource::Carcass { carcass, role }) =
            model::panel_feature(document, part)
                .and_then(|f| model::read_panel(document, f).ok())
                .map(|p| p.source)
        else {
            continue;
        };
        if !matches!(
            role,
            model::PanelRole::LeftSide | model::PanelRole::RightSide
        ) {
            continue;
        }
        for (other, _, body) in model::carcass_panels(document, carcass).unwrap_or_default() {
            if matches!(
                other,
                model::PanelRole::LeftSide | model::PanelRole::RightSide
            ) && !out.contains(&body)
            {
                out.push(body);
            }
        }
    }
    out
}

/// The face frame of a panel Part's concrete face.
pub(super) fn frame(
    document: &Document,
    parameters: &EvaluatedParameters,
    part: EntityId,
    face: MachiningFace,
) -> Option<MachiningFaceFrame> {
    let b = model::machining_panel_box(document, parameters, part).ok()?;
    Some(model::panel_face_frame(&b, face))
}

/// The panel face of `part` that contains a world point.
pub(super) fn locate(
    document: &Document,
    parameters: &EvaluatedParameters,
    part: EntityId,
    point: [f64; 3],
) -> Option<FaceHit> {
    let b = model::machining_panel_box(document, parameters, part).ok()?;
    FACES
        .iter()
        .filter_map(|&face| {
            let frame = model::panel_face_frame(&b, face);
            let (at, depth) = frame.local(point);
            let inside = (-ON_FACE_MM..=frame.size[0] + ON_FACE_MM).contains(&at[0])
                && (-ON_FACE_MM..=frame.size[1] + ON_FACE_MM).contains(&at[1]);
            (inside && depth.abs() <= ON_FACE_MM).then_some((
                depth.abs(),
                FaceHit {
                    part,
                    face,
                    frame,
                    at: [
                        at[0].clamp(0., frame.size[0]),
                        at[1].clamp(0., frame.size[1]),
                    ],
                },
            ))
        })
        .min_by(|a, b| a.0.total_cmp(&b.0))
        .map(|(_, hit)| hit)
}

/// Snap face coordinates to `grid` millimeters, and to the face center
/// lines and edges within `pull` millimeters.
pub(super) fn snap(at: [f64; 2], size: [f64; 2], grid: f64, pull: f64) -> [f64; 2] {
    std::array::from_fn(|i| {
        let v = at[i];
        for target in [0., size[i] / 2., size[i]] {
            if (v - target).abs() <= pull {
                return target;
            }
        }
        if grid > 0. {
            ((v / grid).round() * grid).clamp(0., size[i])
        } else {
            v
        }
    })
}

/// The panel face under the pointer, if the scene is current.
pub(super) fn hover(
    document: &Document,
    revision: u64,
    view: &TerminationViewport<'_>,
) -> Option<FaceHit> {
    let scene = view.scene?;
    if scene.stamp().document != document.id() || scene.stamp().revision != revision {
        return None;
    }
    let pointer = view.pointer.filter(|p| view.rect.contains(*p))?;
    let hit = scene.pick_screen_with(
        view.camera,
        view.rect,
        pointer,
        PickOptions {
            filter: PickFilter::Faces,
            include_hidden: false,
            section: view.section,
        },
    )?;
    let parameters = multicad_core::parameters::evaluate(document).ok()?;
    locate(document, &parameters, hit.instance, hit.point_mm)
}

/// The Part under the pointer, panel or not, if the scene is current. Tools
/// that only take panels use it to explain why another object is refused.
pub(super) fn hovered_part(
    document: &Document,
    revision: u64,
    view: &TerminationViewport<'_>,
) -> Option<EntityId> {
    let scene = view.scene?;
    if scene.stamp().document != document.id() || scene.stamp().revision != revision {
        return None;
    }
    let pointer = view.pointer.filter(|p| view.rect.contains(*p))?;
    scene
        .pick_screen_with(
            view.camera,
            view.rect,
            pointer,
            PickOptions {
                filter: PickFilter::Faces,
                include_hidden: false,
                section: view.section,
            },
        )
        .map(|hit| hit.instance)
}

/// A panel's whole box as a guide box.
pub(super) fn panel_box(b: &model::PanelBox) -> model::GuideBox {
    let n = b.normal();
    model::GuideBox {
        origin: b.origin,
        edges: [
            b.u.map(|c| c * b.length_mm),
            b.v.map(|c| c * b.width_mm),
            n.map(|c| c * b.thickness_mm),
        ],
    }
}
/// Draw the twelve edges of a box; returns its highest on-screen corner for
/// a label.
pub(super) fn box_edges(
    view: &TerminationViewport<'_>,
    b: &model::GuideBox,
    width: f32,
    color: egui::Color32,
) -> egui::Pos2 {
    let p = b.corners().map(|c| screen(view, c));
    for (i, j) in model::GuideBox::EDGES {
        view.painter
            .line_segment([p[i], p[j]], egui::Stroke::new(width, color));
    }
    p.into_iter()
        .min_by(|a, b| a.y.total_cmp(&b.y))
        .unwrap_or_default()
}

fn screen(view: &TerminationViewport<'_>, p: [f64; 3]) -> egui::Pos2 {
    let s = view
        .camera
        .project([p[0] as f32, p[1] as f32, p[2] as f32], view.rect);
    egui::pos2(s[0], s[1])
}
/// The screen position of a face point.
pub(super) fn on_screen(
    view: &TerminationViewport<'_>,
    frame: &MachiningFaceFrame,
    at: [f64; 2],
) -> egui::Pos2 {
    screen(view, frame.point(at))
}

pub(super) fn outline(
    view: &TerminationViewport<'_>,
    frame: &MachiningFaceFrame,
    color: egui::Color32,
) {
    let [l, w] = frame.size;
    polygon(view, frame, &[[0., 0.], [l, 0.], [l, w], [0., w]], color);
}
pub(super) fn polygon(
    view: &TerminationViewport<'_>,
    frame: &MachiningFaceFrame,
    points: &[[f64; 2]],
    color: egui::Color32,
) {
    let stroke = egui::Stroke::new(2., color);
    let p: Vec<_> = points.iter().map(|q| on_screen(view, frame, *q)).collect();
    for i in 0..p.len() {
        view.painter
            .line_segment([p[i], p[(i + 1) % p.len()]], stroke);
    }
}
pub(super) fn circle(
    view: &TerminationViewport<'_>,
    frame: &MachiningFaceFrame,
    at: [f64; 2],
    diameter: f64,
    color: egui::Color32,
) {
    let r = diameter / 2.;
    let points: Vec<[f64; 2]> = (0..24)
        .map(|k| {
            let t = k as f64 / 24. * std::f64::consts::TAU;
            [at[0] + r * t.portable_cos(), at[1] + r * t.portable_sin()]
        })
        .collect();
    polygon(view, frame, &points, color);
}
/// Text on a dark plate, so the bright annotation colors read the same over
/// light and dark materials and over the grid.
pub(super) fn label(
    view: &TerminationViewport<'_>,
    at: egui::Pos2,
    text: &str,
    color: egui::Color32,
) {
    plate(view, at, egui::Align2::LEFT_BOTTOM, text, color);
}
fn plate(
    view: &TerminationViewport<'_>,
    at: egui::Pos2,
    align: egui::Align2,
    text: &str,
    color: egui::Color32,
) {
    let galley =
        view.painter
            .layout_no_wrap(text.to_owned(), egui::FontId::proportional(12.), color);
    let rect = align.anchor_size(at, galley.size());
    view.painter.rect_filled(
        rect.expand2(egui::vec2(3., 1.)),
        3.,
        LABEL_PLATE.gamma_multiply(color.a() as f32 / 255.),
    );
    view.painter.galley(rect.min, galley, color);
}
/// Whether the outside of a face points towards the viewer. Annotations sit
/// on top of the shaded model, so faces seen from behind are drawn faded:
/// the overlay must not present a hidden face as the visible one.
pub(super) fn faces_viewer(view: &TerminationViewport<'_>, frame: &MachiningFaceFrame) -> bool {
    let center = on_screen(view, frame, [frame.size[0] / 2., frame.size[1] / 2.]);
    view.camera.ray(center, view.rect).is_none_or(|ray| {
        let d = ray.direction();
        // An edge-on face shows no area; f32 camera angles leave noise there.
        (0..3).map(|i| frame.inward[i] * d[i]).sum::<f64>() > EDGE_ON
    })
}
/// The annotation color of a face seen from behind.
pub(super) fn far_side(color: egui::Color32) -> egui::Color32 {
    color.gamma_multiply(0.4)
}
/// Arrows for the face's u (first) and v (second) axes at its origin.
pub(super) fn axes(view: &TerminationViewport<'_>, frame: &MachiningFaceFrame) {
    let o = on_screen(view, frame, [0., 0.]);
    let len = frame.size[0].min(frame.size[1]).clamp(20., 120.);
    for (to, color, name) in [([len, 0.], AXIS_U, "u"), ([0., len], AXIS_V, "v")] {
        let p = on_screen(view, frame, to);
        view.painter.arrow(o, p - o, egui::Stroke::new(2., color));
        label(view, p, name, color);
    }
}
/// Measured offsets from the face's start edges to a point.
pub(super) fn offsets(view: &TerminationViewport<'_>, frame: &MachiningFaceFrame, at: [f64; 2]) {
    let stroke = egui::Stroke::new(1., HOVER);
    for (from, color, value) in [([0., at[1]], AXIS_U, at[0]), ([at[0], 0.], AXIS_V, at[1])] {
        let (a, b) = (on_screen(view, frame, from), on_screen(view, frame, at));
        view.painter.line_segment([a, b], stroke);
        // Centered on its leader like a dimension value, clear of the
        // placement label beside the pointer.
        plate(
            view,
            a + (b - a) / 2.,
            egui::Align2::CENTER_CENTER,
            &format!("{value:.1}"),
            color,
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Which faces of the 600 × 400 × 18 panel face a camera.
    fn facing(camera: crate::camera::Camera) -> Vec<MachiningFace> {
        let (s, side) = crate::modeling::panel_fixture::flat_panel();
        let parameters = multicad_core::parameters::evaluate(s.document()).unwrap();
        let rect = egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(600., 600.));
        let ctx = egui::Context::default();
        let mut out = vec![];
        ctx.run_ui(Default::default(), |ui| {
            let view = TerminationViewport {
                scene: None,
                report: None,
                camera,
                rect,
                pointer: None,
                clicked: false,
                painter: ui.painter(),
                section: None,
            };
            out = FACES
                .iter()
                .copied()
                .filter(|&f| {
                    let frame = frame(s.document(), &parameters, side, f).unwrap();
                    faces_viewer(&view, &frame)
                })
                .collect();
        })
        .drop_without_applying_deltas();
        out
    }

    #[test]
    fn only_faces_turned_towards_the_camera_read_as_visible() {
        let mut camera = crate::camera::Camera {
            yaw: 0.,
            pitch: std::f32::consts::FRAC_PI_2,
            scale: 0.5,
            target: [300., 200., 0.],
            perspective: false,
        };
        // Straight down: only the top face.
        assert_eq!(facing(camera), vec![MachiningFace::Top]);
        // From below the far side is the bottom face.
        camera.pitch = -std::f32::consts::FRAC_PI_2;
        assert_eq!(facing(camera), vec![MachiningFace::Bottom]);
        // An oblique perspective view from above and in front shows the top
        // and the near edges, never the bottom or both opposite edges.
        camera.pitch = 0.5;
        camera.yaw = 0.6;
        camera.perspective = true;
        let seen = facing(camera);
        assert!(seen.contains(&MachiningFace::Top), "{seen:?}");
        assert!(!seen.contains(&MachiningFace::Bottom), "{seen:?}");
        for pair in [
            [PanelEdgeName::UStart, PanelEdgeName::UEnd],
            [PanelEdgeName::VStart, PanelEdgeName::VEnd],
        ] {
            let n = pair
                .iter()
                .filter(|&&edge| seen.contains(&MachiningFace::Edge { edge }))
                .count();
            assert_eq!(n, 1, "{seen:?}");
        }
    }

    #[test]
    fn snapping_prefers_edges_and_center_lines_then_the_grid() {
        let size = [600., 400.];
        assert_eq!(snap([3., 197.], size, 1., 5.), [0., 200.]);
        assert_eq!(snap([123.4, 57.6], size, 1., 5.), [123., 58.]);
        assert_eq!(snap([598., 401.], size, 1., 5.), [600., 400.]);
        assert_eq!(snap([123.4, 57.6], size, 0., 0.), [123.4, 57.6]);
    }
}
