//! Display handles produce typed rigid-placement controls. They never edit a
//! mesh or commit a document; the controller requires a complete native receipt.
pub(crate) mod math;
use crate::{camera::Camera, egui, viewport::Scene};
use egui::{Color32, Pos2, Rect, Stroke};
use math::*;
use multicad_core::{model::NumericTransform, parameters::EvaluatedParameters};

/// Right-handed orthonormal world directions for World, Object or a datum.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GizmoBasis {
    pub axes: [[f64; 3]; 3],
}
impl Default for GizmoBasis {
    fn default() -> Self {
        Self {
            axes: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        }
    }
}
impl GizmoBasis {
    pub(super) fn validate(self) -> Result<(), String> {
        if self.axes.iter().flatten().any(|v| !v.is_finite())
            || self.axes.iter().any(|v| (dot(*v, *v) - 1.).abs() > 1e-6)
            || dot(cross(self.axes[0], self.axes[1]), self.axes[2]) < 1. - 1e-6
            || dot(self.axes[0], self.axes[1]).abs() > 1e-6
        {
            return Err("Choose an orthonormal right-handed transform coordinate system".into());
        }
        Ok(())
    }
}
/// Default Move/Rotate increments. Shift held during a gesture releases them
/// temporarily; the gesture keeps its handle, frame and pivot either way.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransformSnap {
    /// Active grid spacing, applied to the drag displacement along the
    /// dragged axis or plane directions of the transform coordinate frame.
    pub grid_mm: Option<f64>,
    /// Angular increment about the dragged rotation axis.
    pub angle_deg: Option<f64>,
}
impl TransformSnap {
    pub const DEFAULT_ANGLE_DEG: f64 = 45.;
    pub const OFF: Self = Self {
        grid_mm: None,
        angle_deg: None,
    };
    /// Grid steps at the active spacing and 45-degree rotation steps.
    pub fn grid(spacing_mm: f64) -> Self {
        Self {
            grid_mm: (spacing_mm.is_finite() && spacing_mm > 0.).then_some(spacing_mm),
            angle_deg: Some(Self::DEFAULT_ANGLE_DEG),
        }
    }
}
impl Default for TransformSnap {
    fn default() -> Self {
        Self::grid(10.)
    }
}
pub struct PlacementViewport<'a> {
    /// Complete canonical source scene, never the noninteractive candidate.
    pub scene: Option<&'a Scene>,
    pub camera: Camera,
    pub rect: Rect,
    pub basis: GizmoBasis,
    /// The canvas owns the pointer rather than a popup, orientation widget or
    /// camera control. This gates new gestures; captured gestures still finish.
    pub start_allowed: bool,
    /// Increments for new gestures; a running gesture keeps its own.
    pub snap: TransformSnap,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Handle {
    Axis(usize),
    Plane(usize),
    Rotate(usize),
}
impl Handle {
    fn label(self) -> &'static str {
        match self {
            Self::Axis(0) => "Move X",
            Self::Axis(1) => "Move Y",
            Self::Axis(_) => "Move Z",
            Self::Plane(0) => "Move in YZ plane",
            Self::Plane(1) => "Move in ZX plane",
            Self::Plane(_) => "Move in XY plane",
            Self::Rotate(0) => "Rotate X",
            Self::Rotate(1) => "Rotate Y",
            Self::Rotate(_) => "Rotate Z",
        }
    }
}
struct HandleShape {
    handle: Handle,
    points: Vec<Pos2>,
    color: Color32,
    closed: bool,
}
impl HandleShape {
    fn distance(&self, point: Pos2) -> f32 {
        if self.closed
            && egui::epaint::PathShape::convex_polygon(
                self.points.clone(),
                Color32::WHITE,
                Stroke::NONE,
            )
            .visual_bounding_rect()
            .contains(point)
        {
            let mut signs = self
                .points
                .iter()
                .zip(self.points.iter().cycle().skip(1))
                .take(self.points.len())
                .map(|(a, b)| (*b - *a).rot90().dot(point - *a));
            let first = signs.next().unwrap_or_default();
            if signs.all(|v| v * first >= 0.) {
                return 0.;
            }
        }
        self.points
            .windows(2)
            .map(|p| segment_distance(point, p[0], p[1]))
            .fold(f32::INFINITY, f32::min)
    }
}
fn segment_distance(p: Pos2, a: Pos2, b: Pos2) -> f32 {
    let d = b - a;
    let t = ((p - a).dot(d) / d.length_sq().max(1e-8)).clamp(0., 1.);
    p.distance(a + d * t)
}
fn projected(camera: Camera, rect: Rect, p: Point) -> Option<Pos2> {
    let p = camera.project(p.map(|v| v as f32), rect);
    (p.iter().all(|v| v.is_finite()) && (!camera.perspective || p[2] > 0.1))
        .then(|| Pos2::new(p[0], p[1]))
}
fn shapes(camera: Camera, rect: Rect, origin: Point, basis: GizmoBasis) -> Vec<HandleShape> {
    let Some(center) = projected(camera, rect, origin) else {
        return vec![];
    };
    if !rect.expand(100.).contains(center) {
        return vec![];
    }
    let depth = camera.project(origin.map(|v| v as f32), rect)[2];
    let distance = rect.height().max(1.) / camera.scale * 2.;
    let factor = if camera.perspective {
        distance / depth
    } else {
        1.
    };
    let length = 82. / (camera.scale * factor) as f64;
    let colors = [
        Color32::from_rgb(194, 55, 61),
        Color32::from_rgb(39, 146, 88),
        Color32::from_rgb(52, 112, 211),
    ];
    let mut result = vec![];
    for (i, color) in colors.into_iter().enumerate() {
        let point = |a: usize, s: f64| add(origin, mul(basis.axes[a], s * length));
        if let Some(end) = projected(camera, rect, point(i, 1.))
            && end.distance(center) > 16.
        {
            result.push(HandleShape {
                handle: Handle::Axis(i),
                points: vec![center, end],
                color,
                closed: false,
            });
        }
        let a = (i + 1) % 3;
        let b = (i + 2) % 3;
        let quad = [[0.17, 0.17], [0.34, 0.17], [0.34, 0.34], [0.17, 0.34]]
            .into_iter()
            .map(|[u, v]| {
                projected(
                    camera,
                    rect,
                    add(
                        add(origin, mul(basis.axes[a], u * length)),
                        mul(basis.axes[b], v * length),
                    ),
                )
            })
            .collect::<Option<Vec<_>>>();
        if let Some(points) = quad
            && (points[1] - points[0])
                .rot90()
                .dot(points[2] - points[0])
                .abs()
                > 8.
        {
            result.push(HandleShape {
                handle: Handle::Plane(i),
                points,
                color,
                closed: true,
            });
        }
        let ring = (0..=64)
            .map(|j| {
                let t = j as f64 / 64. * std::f64::consts::TAU;
                projected(
                    camera,
                    rect,
                    add(
                        add(origin, mul(basis.axes[a], t.cos() * length * 0.74)),
                        mul(basis.axes[b], t.sin() * length * 0.74),
                    ),
                )
            })
            .collect::<Option<Vec<_>>>();
        if let Some(points) = ring {
            result.push(HandleShape {
                handle: Handle::Rotate(i),
                points,
                color,
                closed: false,
            });
        }
    }
    result
}
struct Drag {
    handle: Handle,
    source: NumericTransform,
    values: Values,
    camera: Camera,
    rect: Rect,
    basis: GizmoBasis,
    start: Point,
    last_angle: f64,
    /// Unsnapped accumulated rotation; the staged angle is derived from it.
    angle: f64,
    last: NumericTransform,
    changed: bool,
    snap: TransformSnap,
    /// Shift held: increments are released until Shift is let go.
    free: bool,
    readout: Option<String>,
}
impl Drag {
    fn new(
        handle: Handle,
        source: &NumericTransform,
        parameters: &EvaluatedParameters,
        input: &PlacementViewport<'_>,
        pointer: Pos2,
    ) -> Result<Self, String> {
        let values = Values::read(source, parameters)?;
        let ray = input
            .camera
            .ray(pointer, input.rect)
            .ok_or("No finite pointer ray")?;
        let origin = values.origin();
        let start = match handle {
            Handle::Axis(i) => [
                axis_parameter(ray, origin, input.basis.axes[i])
                    .ok_or("This axis is edge-on; orbit the view before dragging")?,
                0.,
                0.,
            ],
            Handle::Plane(i) | Handle::Rotate(i) => {
                plane_point(ray, origin, input.basis.axes[i])
                    .ok_or("This plane is edge-on; orbit the view before dragging")?
            }
        };
        if matches!(handle, Handle::Rotate(_))
            && dot(sub(start, origin), sub(start, origin)) < 1e-12
        {
            return Err("Choose a point on the rotation ring".into());
        }
        Ok(Self {
            handle,
            source: source.clone(),
            last: source.clone(),
            values,
            camera: input.camera,
            rect: input.rect,
            basis: input.basis,
            start,
            last_angle: 0.,
            angle: 0.,
            changed: false,
            snap: input.snap,
            free: false,
            readout: None,
        })
    }
    fn update(
        &mut self,
        pointer: Pos2,
        parameters: &EvaluatedParameters,
    ) -> Result<Option<NumericTransform>, String> {
        let ray = self
            .camera
            .ray(pointer, self.rect)
            .ok_or("No finite pointer ray")?;
        let origin = self.values.origin();
        let next = match self.handle {
            Handle::Axis(i) => {
                let end =
                    axis_parameter(ray, origin, self.basis.axes[i]).ok_or("Axis became edge-on")?;
                let distance = snap_step(end - self.start[0], self.grid());
                self.readout = Some(format!(
                    "{} {} · {}",
                    self.handle.label(),
                    crate::preferences::format_mm(distance),
                    self.grid_label()
                ));
                translated(&self.source, self.values, mul(self.basis.axes[i], distance))?
            }
            Handle::Plane(i) => {
                let end =
                    plane_point(ray, origin, self.basis.axes[i]).ok_or("Plane became edge-on")?;
                let (a, b) = ((i + 1) % 3, (i + 2) % 3);
                let delta = snap_planar(
                    sub(end, self.start),
                    self.basis.axes[a],
                    self.basis.axes[b],
                    self.grid(),
                );
                self.readout = Some(format!(
                    "{} {}, {} {} · {}",
                    ["X", "Y", "Z"][a],
                    crate::preferences::format_mm(dot(delta, self.basis.axes[a])),
                    ["X", "Y", "Z"][b],
                    crate::preferences::format_mm(dot(delta, self.basis.axes[b])),
                    self.grid_label()
                ));
                translated(&self.source, self.values, delta)?
            }
            Handle::Rotate(i) => {
                let end = plane_point(ray, origin, self.basis.axes[i])
                    .ok_or("Rotation plane became edge-on")?;
                let a = sub(self.start, origin);
                let b = sub(end, origin);
                if dot(b, b) < 1e-12 {
                    return Ok(None);
                }
                let angle = dot(self.basis.axes[i], cross(a, b)).atan2(dot(a, b));
                let mut delta = angle - self.last_angle;
                while delta > std::f64::consts::PI {
                    delta -= std::f64::consts::TAU;
                }
                while delta < -std::f64::consts::PI {
                    delta += std::f64::consts::TAU;
                }
                self.angle += delta;
                self.last_angle = angle;
                let step = (!self.free).then_some(self.snap.angle_deg).flatten();
                let degrees = snap_step(self.angle.to_degrees(), step);
                self.readout = Some(format!(
                    "{} {:.1}° · {}",
                    self.handle.label(),
                    degrees,
                    match step {
                        Some(step) => format!("{step:.0}° steps · Shift for free"),
                        None if self.snap.angle_deg.is_some() => "free (Shift)".into(),
                        None => "free".into(),
                    }
                ));
                rotated(
                    &self.source,
                    self.values,
                    self.basis.axes[i],
                    degrees,
                    parameters,
                )?
            }
        };
        if next == self.last {
            return Ok(None);
        }
        self.changed = true;
        self.last = next.clone();
        Ok(Some(next))
    }
}
impl Drag {
    fn grid(&self) -> Option<f64> {
        (!self.free).then_some(self.snap.grid_mm).flatten()
    }
    fn grid_label(&self) -> String {
        match self.grid() {
            Some(step) => format!(
                "grid {} · Shift for free",
                crate::preferences::format_mm(step)
            ),
            None if self.snap.grid_mm.is_some() => "free (Shift)".into(),
            None => "free".into(),
        }
    }
}
#[derive(Default)]
pub(crate) struct Interaction {
    drag: Option<Drag>,
}
pub(crate) struct Gesture {
    pub input_owned: bool,
    pub started: bool,
    pub changed: Option<NumericTransform>,
    pub released: bool,
    pub moved: bool,
    pub copy_modifier: bool,
}
impl Interaction {
    pub fn active(&self) -> bool {
        self.drag.is_some()
    }
    pub fn show(
        &mut self,
        ui: &mut egui::Ui,
        input: &PlacementViewport<'_>,
        source: &NumericTransform,
        parameters: &EvaluatedParameters,
    ) -> Result<Gesture, String> {
        self.show_with_help(
            ui,
            input,
            source,
            parameters,
            "drag to stage, Enter to apply, Escape to cancel. Moves step by the grid and rotations by 45°; hold Shift while dragging for free placement. Option / Alt copies.",
        )
    }
    pub fn show_with_help(
        &mut self,
        ui: &mut egui::Ui,
        input: &PlacementViewport<'_>,
        source: &NumericTransform,
        parameters: &EvaluatedParameters,
        help: &str,
    ) -> Result<Gesture, String> {
        input.basis.validate()?;
        if let Some(drag) = &self.drag
            && (drag.camera != input.camera || drag.rect != input.rect || drag.basis != input.basis)
        {
            return Err(
                "The view or transform coordinates changed; the placement gesture was cancelled"
                    .into(),
            );
        }
        let (pointer, press, pressed, released, down, alt, shift, motion) = ui.input(|i| {
            // Mouse-up may precede an Option/Alt key-up in the same frame.
            // The copy intent belongs to the gesture's pointer event.
            let alt = i
                .events
                .iter()
                .rev()
                .find_map(|event| match event {
                    egui::Event::PointerButton {
                        button: egui::PointerButton::Primary,
                        modifiers,
                        ..
                    } => Some(modifiers.alt),
                    _ => None,
                })
                .unwrap_or(i.modifiers.alt);
            let press = i
                .events
                .iter()
                .enumerate()
                .rev()
                .find_map(|(index, event)| match event {
                    egui::Event::PointerButton {
                        pos,
                        button: egui::PointerButton::Primary,
                        pressed: true,
                        ..
                    } => Some((index, *pos)),
                    _ => None,
                });
            // Input can contain an entire drag between rendered frames. Keep
            // its original down position and ordered motion/rotation samples.
            let motion = i
                .events
                .iter()
                .skip(press.map_or(0, |(index, _)| index + 1))
                .filter_map(|event| match event {
                    egui::Event::PointerMoved(pos)
                    | egui::Event::PointerButton {
                        pos,
                        button: egui::PointerButton::Primary,
                        pressed: false,
                        ..
                    } => Some(*pos),
                    _ => None,
                })
                .collect::<Vec<_>>();
            (
                i.pointer.interact_pos(),
                press.map(|(_, pos)| pos),
                i.pointer.primary_pressed(),
                i.pointer.primary_released(),
                i.pointer.primary_down(),
                alt,
                i.modifiers.shift,
                motion,
            )
        });
        let values = Values::read(source, parameters)?;
        if self.drag.is_some() && !down && !pressed && !released {
            return Err("The pointer gesture was interrupted; the placement was cancelled".into());
        }
        // Paint at the current numeric pivot while solving every pointer delta
        // against the frozen start frame, so translation cannot accumulate.
        let draw_origin = values.origin();
        let draw_basis = self.drag.as_ref().map_or(input.basis, |d| d.basis);
        let geometry = shapes(input.camera, input.rect, draw_origin, draw_basis);
        let hit_pointer = if self.drag.is_none() && pressed {
            press.or(pointer)
        } else {
            pointer
        };
        let hovered = hit_pointer
            .filter(|p| input.rect.contains(*p) && (input.start_allowed || self.drag.is_some()))
            .and_then(|p| {
                geometry
                    .iter()
                    .map(|s| (s, s.distance(p)))
                    .filter(|(_, d)| *d < 7.)
                    .min_by(|a, b| a.1.total_cmp(&b.1))
                    .map(|(s, _)| s.handle)
            });
        let mut started = false;
        if self.drag.is_none()
            && input.start_allowed
            && pressed
            && let (Some(handle), Some(pointer)) = (hovered, press.or(pointer))
        {
            self.drag = Some(Drag::new(handle, source, parameters, input, pointer)?);
            started = true;
        }
        let active = self.drag.as_ref().map(|d| d.handle);
        let painter = ui.painter_at(input.rect);
        for shape in geometry {
            let highlighted = active == Some(shape.handle) || hovered == Some(shape.handle);
            let color = if highlighted {
                Color32::WHITE
            } else {
                shape.color
            };
            if shape.closed {
                painter.add(egui::Shape::convex_polygon(
                    shape.points,
                    color.gamma_multiply(0.18),
                    Stroke::new(1., color),
                ));
            } else if matches!(shape.handle, Handle::Axis(_)) {
                let end = *shape.points.last().unwrap();
                let start = shape.points[0];
                painter.arrow(
                    start,
                    end - start,
                    Stroke::new(if highlighted { 3. } else { 2. }, color),
                );
                if let Handle::Axis(i) = shape.handle {
                    painter.text(
                        end + egui::vec2(5., -5.),
                        egui::Align2::LEFT_BOTTOM,
                        ["X", "Y", "Z"][i],
                        egui::FontId::proportional(11.),
                        color,
                    );
                }
            } else {
                painter.add(egui::Shape::line(
                    shape.points,
                    Stroke::new(if highlighted { 2.5 } else { 1.5 }, color),
                ));
            }
        }
        if let Some(handle) = hovered {
            ui.ctx().set_cursor_icon(egui::CursorIcon::Grab);
            if let Some(pointer) = pointer {
                ui.interact(
                    Rect::from_center_size(pointer, egui::vec2(12., 12.)),
                    ui.id().with("placement-handle-tooltip"),
                    egui::Sense::hover(),
                )
                .on_hover_text(format!("{} · {help}", handle.label()));
            }
        }
        let mut changed = None;
        if let Some(drag) = &mut self.drag
            && (down || released)
        {
            // Shift only releases increments inside this captured gesture.
            // Toggling it re-solves the current pointer without ending it.
            drag.free = shift;
            let points = if motion.is_empty() {
                pointer.into_iter().collect()
            } else {
                motion
            };
            for point in points {
                if let Some(next) = drag.update(point, parameters)? {
                    changed = Some(next);
                }
            }
        }
        if let Some(drag) = &self.drag
            && drag.changed
            && let Some(text) = &drag.readout
            && let Some(anchor) = projected(
                input.camera,
                input.rect,
                values_origin(changed.as_ref().unwrap_or(source), parameters),
            )
        {
            paint_readout(&painter, input.rect, anchor + egui::vec2(18., 18.), text);
        }
        let moved = self.drag.as_ref().is_some_and(|d| d.changed);
        let owned = self.drag.is_some();
        let released = owned && released;
        if released {
            self.drag = None;
        }
        Ok(Gesture {
            input_owned: owned,
            started,
            changed,
            released,
            moved,
            copy_modifier: alt,
        })
    }
}

fn values_origin(source: &NumericTransform, parameters: &EvaluatedParameters) -> Point {
    Values::read(source, parameters).map_or([0.; 3], |values| values.origin())
}
/// Staged-value label beside the gizmo. It reports the value the preview uses.
/// It stays inside the viewport so a narrow canvas never clips the value.
fn paint_readout(painter: &egui::Painter, viewport: Rect, at: Pos2, text: &str) {
    let galley = painter.layout_no_wrap(
        text.to_owned(),
        egui::FontId::proportional(12.),
        Color32::from_gray(235),
    );
    let inner = viewport.shrink(6.);
    let at = Pos2::new(
        at.x.min(inner.right() - galley.size().x).max(inner.left()),
        at.y.min(inner.bottom() - galley.size().y).max(inner.top()),
    );
    let rect = Rect::from_min_size(at, galley.size()).expand(4.);
    painter.rect_filled(rect, 3., Color32::from_black_alpha(190));
    painter.galley(at, galley, Color32::from_gray(235));
}

#[cfg(test)]
mod tests;
