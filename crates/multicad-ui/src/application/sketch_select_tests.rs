//! CAD-412: the Select tool picks sketch geometry in the real canvas.
//! Every click goes through the full application frame and input handling:
//! the ribbon, the canvas, the element list and the Inspector.
use super::ui_harness::Harness;
use super::*;
use crate::sketch_editor::{DrawTool, SelectionSnapshot};
use multicad_core::constraint_sketch::{CurveKind, source::Source};
use multicad_precision::PortableMath as _;

/// Stroke colour of selected sketch curves.
const HIGHLIGHT: egui::Color32 = egui::Color32::from_rgb(255, 195, 80);

pub(super) struct SketchHarness {
    pub(super) h: Harness,
}
impl SketchHarness {
    /// An empty project with a new XY sketch open in Model.
    pub(super) fn new() -> Self {
        let mut h = Harness::new(Document::new("Sketch select").unwrap());
        h.app
            .invoke(Command::NewSketch, &["xy"], Channel::Ui)
            .unwrap();
        let mut s = Self { h };
        s.settle();
        s
    }
    /// Render until the open sketch's solution is current and painted.
    pub(super) fn settle(&mut self) {
        self.h.frame(vec![]);
        self.h.until("the sketch solution", |app| {
            !app.current().sketch.sketch_evaluation_pending()
        });
        self.h.frame(vec![]);
    }
    pub(super) fn screen(&self, p: [f64; 2]) -> egui::Pos2 {
        let project = self.h.app.current();
        project.sketch.sketch_screen(&project.session, p)
    }
    /// Where a user clicks `p`: the sketch must be solved and painted first,
    /// since a tool change can start a new solve on the worker.
    pub(super) fn click_point(&mut self, p: [f64; 2]) -> egui::Pos2 {
        if self.h.app.current().sketch.sketch_evaluation_pending() {
            self.settle();
        }
        self.screen(p)
    }
    pub(super) fn button(
        &mut self,
        at: egui::Pos2,
        button: egui::PointerButton,
        modifiers: egui::Modifiers,
    ) {
        self.h.frame(vec![egui::Event::PointerMoved(at)]);
        for pressed in [true, false] {
            self.h.frame(vec![
                egui::Event::PointerMoved(at),
                egui::Event::PointerButton {
                    pos: at,
                    button,
                    pressed,
                    modifiers,
                },
            ]);
        }
        self.settle();
    }
    pub(super) fn click_with(
        &mut self,
        p: [f64; 2],
        offset: egui::Vec2,
        modifiers: egui::Modifiers,
    ) {
        let at = self.click_point(p) + offset;
        self.button(at, egui::PointerButton::Primary, modifiers);
    }
    pub(super) fn click(&mut self, p: [f64; 2], offset: egui::Vec2) {
        self.click_with(p, offset, egui::Modifiers::NONE);
    }
    /// Choose a drawing tool through its command, as its ribbon button does.
    pub(super) fn tool(&mut self, tool: DrawTool) {
        self.h
            .app
            .invoke(Command::Draw(tool), &[], Channel::Ui)
            .unwrap();
        self.h.frame(vec![]);
    }
    pub(super) fn draw(&mut self, tool: DrawTool, args: &[&str]) {
        self.h
            .app
            .invoke(Command::Draw(tool), args, Channel::Ui)
            .unwrap();
        self.settle();
    }
    pub(super) fn source(&self) -> Source {
        let project = self.h.app.current();
        let id = project.sketch.selected.unwrap();
        project
            .sketch
            .canonical_source(project.session.document(), id)
    }
    pub(super) fn selection(&self) -> SelectionSnapshot {
        let project = self.h.app.current();
        project
            .sketch
            .selection_snapshot(project.session.document())
    }
    pub(super) fn message(&self) -> String {
        self.h.app.current().sketch.status_message().to_owned()
    }
    /// Painted selection highlight: curve segments within 3 px of `p`, or a
    /// selected point marker at `p`.
    pub(super) fn highlighted_at(&mut self, p: [f64; 2]) -> bool {
        let at = self.click_point(p);
        let out = self.h.render(vec![]);
        let hit = out.shapes.iter().any(|clipped| match &clipped.shape {
            egui::Shape::LineSegment {
                points: [a, b],
                stroke,
            } if stroke.color == HIGHLIGHT => {
                let v = *b - *a;
                let t = ((at - *a).dot(v) / v.length_sq().max(f32::EPSILON)).clamp(0., 1.);
                at.distance(*a + v * t) <= 3.
            }
            egui::Shape::Circle(circle) => {
                circle.fill == egui::Color32::YELLOW && circle.center.distance(at) <= 1.
            }
            _ => false,
        });
        out.drop_without_applying_deltas();
        hit
    }
    /// Click the drag-value field painted beside or under `caption`, replace its
    /// value and press Enter, as a user types a value.
    pub(super) fn type_beside(&mut self, caption: &str, text: &str) {
        let out = self.h.render(vec![]);
        let label = out
            .shapes
            .iter()
            .find_map(|s| match &s.shape {
                egui::Shape::Text(t) if t.galley.text() == caption => {
                    Some(egui::Rect::from_min_size(t.pos, t.galley.size()))
                }
                _ => None,
            })
            .unwrap_or_else(|| panic!("Missing {caption}"));
        let field = out
            .platform_output
            .accesskit_update
            .as_ref()
            .unwrap()
            .nodes
            .iter()
            .filter(|(_, n)| {
                matches!(
                    n.role(),
                    egui::accesskit::Role::SpinButton | egui::accesskit::Role::TextInput
                )
            })
            .filter_map(|(_, n)| n.bounds())
            .filter(|b| {
                // Beside the caption on its row, or (Inspector fields,
                // CAD-416) right under it with the label above the input.
                let beside = b.x0 as f32 >= label.max.x - 1.
                    && (b.y0 as f32) <= label.center().y
                    && (b.y1 as f32) >= label.center().y;
                let below = b.y0 as f32 >= label.max.y - 1.
                    && (b.x0 as f32) <= label.max.x
                    && (b.x1 as f32) >= label.min.x;
                beside || below
            })
            .min_by(|a, b| a.y0.total_cmp(&b.y0))
            .unwrap_or_else(|| panic!("a field beside or under {caption}"));
        out.drop_without_applying_deltas();
        let at = egui::pos2(
            ((field.x0 + field.x1) * 0.5) as f32,
            ((field.y0 + field.y1) * 0.5) as f32,
        );
        self.h.click_at(at);
        assert!(self.h.ctx.text_edit_focused(), "{caption} field focused");
        let key = |key, pressed, modifiers| egui::Event::Key {
            key,
            physical_key: None,
            pressed,
            repeat: false,
            modifiers,
        };
        for events in [
            vec![key(egui::Key::A, true, egui::Modifiers::COMMAND)],
            vec![key(egui::Key::A, false, egui::Modifiers::COMMAND)],
            vec![
                egui::Event::Text(text.into()),
                key(egui::Key::Enter, true, egui::Modifiers::NONE),
            ],
            vec![key(egui::Key::Enter, false, egui::Modifiers::NONE)],
        ] {
            self.h.frame(events);
        }
        self.settle();
    }
    pub(super) fn timed(&mut self, events: Vec<egui::Event>, time: f64) {
        let app = &mut self.h.app;
        self.h
            .ctx
            .run_ui(
                egui::RawInput {
                    screen_rect: Some(egui::Rect::from_min_size(
                        egui::Pos2::ZERO,
                        super::ui_harness::SIZE,
                    )),
                    events,
                    time: Some(time),
                    focused: true,
                    ..Default::default()
                },
                |ui| app.show(ui),
            )
            .drop_without_applying_deltas();
    }
}

fn edge_selected(s: &SketchHarness, curve: EntityId, what: &str) {
    assert_eq!(
        s.selection().entities,
        vec![curve],
        "{what}: {}",
        s.message()
    );
}
fn point_at(source: &Source, p: [f64; 2]) -> EntityId {
    source
        .sketch
        .points
        .iter()
        .find(|q| (q.position_mm[0] - p[0]).abs() < 1e-6 && (q.position_mm[1] - p[1]).abs() < 1e-6)
        .unwrap_or_else(|| panic!("no point at {p:?}"))
        .id
}
/// Midpoint of an arc curve, from its center and end points.
fn arc_midpoint(source: &Source, arc: &multicad_core::constraint_sketch::Curve) -> [f64; 2] {
    let CurveKind::Arc {
        center,
        start,
        end,
        clockwise,
    } = arc.kind
    else {
        panic!("not an arc")
    };
    let at = |id| {
        source
            .sketch
            .points
            .iter()
            .find(|p| p.id == id)
            .unwrap()
            .position_mm
    };
    let (c, a, b) = (at(center), at(start), at(end));
    let angle = |p: [f64; 2]| (p[1] - c[1]).portable_atan2(p[0] - c[0]);
    let (from, mut to) = (angle(a), angle(b));
    if clockwise {
        while to > from {
            to -= std::f64::consts::TAU;
        }
    } else {
        while to < from {
            to += std::f64::consts::TAU;
        }
    }
    let r = (a[0] - c[0]).portable_hypot(a[1] - c[1]);
    let mid = (from + to) / 2.;
    [c[0] + r * mid.portable_cos(), c[1] + r * mid.portable_sin()]
}
const SHIFT: egui::Modifiers = egui::Modifiers::SHIFT;
const COMMAND: egui::Modifiers = egui::Modifiers::COMMAND;
const ON_LINE: egui::Vec2 = egui::vec2(0., 1.);

/// Reproduction of the reported state: after drawing, the drawing tool stays
/// active, so canvas clicks keep drawing. The ribbon's Select tool ends it.
#[test]
fn the_select_tool_ends_the_drawing_tool_and_a_canvas_click_selects_the_edge() {
    let mut s = SketchHarness::new();
    s.tool(DrawTool::Rectangle);
    s.click([0., 0.], egui::Vec2::ZERO);
    s.click([60., 40.], egui::Vec2::ZERO);
    let rectangle = s.source().sketch.shapes[0].clone();
    assert!(
        s.h.texts().iter().any(|t| t == "Select"),
        "a Select tool in the sketch ribbon"
    );
    s.h.click_label("Select");
    s.settle();
    assert_eq!(s.h.app.current().sketch.tool, None);
    assert!(s.h.app.current().sketch.select_active());
    let revision = s.h.revision();
    let document = s.h.app.current().session.document().clone();
    s.click([30., 0.], ON_LINE);
    let selection = s.selection();
    assert_eq!(
        selection.entities,
        vec![rectangle.curves[0]],
        "{}",
        s.message()
    );
    assert_eq!(selection.element, Some(rectangle.id));
    assert!(s.highlighted_at([30., 0.]), "visible highlight");
    assert_eq!(s.h.revision(), revision);
    assert_eq!(s.h.app.current().session.document(), &document);
}

/// Line, arc, point, rectangle edge and filleted corner each select with a
/// visible highlight; Shift/Cmd toggle; empty space clears; nothing commits.
#[test]
fn select_picks_lines_arcs_points_edges_and_fillets_toggles_and_clears() {
    let mut s = SketchHarness::new();
    s.draw(DrawTool::Rectangle, &["0", "0", "60", "40"]);
    s.draw(DrawTool::Line, &["0", "60", "60", "60"]);
    s.draw(DrawTool::Arc, &["100", "0", "120", "0", "100", "20"]);
    s.tool(DrawTool::Fillet);
    s.click([60., 40.], egui::Vec2::ZERO);
    let source = s.source();
    let [rectangle, line, arc] = [0, 1, 2].map(|i| source.sketch.shapes[i].clone());
    // The fillet is a shape of its own (CAD-414); its one curve is the arc.
    let fillet_shape = source
        .sketch
        .shapes
        .iter()
        .find(|s| s.kind == multicad_core::constraint_sketch::ShapeKind::Fillet)
        .expect("the fillet shape")
        .clone();
    let fillet_arc = source
        .sketch
        .curves
        .iter()
        .find(|c| c.id == fillet_shape.curves[0])
        .expect("the fillet arc");
    let fillet_mid = arc_midpoint(&source, fillet_arc);
    let fillet = fillet_arc.id;
    let corner = point_at(&source, [0., 0.]);
    s.h.click_label("Select");
    s.settle();
    let revision = s.h.revision();
    let document = s.h.app.current().session.document().clone();
    let diagonal = std::f64::consts::FRAC_1_SQRT_2;
    let arc_mid = [100. + 20. * diagonal, 20. * diagonal];
    for (what, p, offset, entity, element) in [
        ("line", [30., 60.], ON_LINE, line.curves[0], Some(line.id)),
        (
            "arc",
            arc_mid,
            egui::Vec2::ZERO,
            arc.curves[0],
            Some(arc.id),
        ),
        ("point", [0., 0.], egui::Vec2::ZERO, corner, None),
        (
            "rectangle edge",
            [30., 0.],
            ON_LINE,
            rectangle.curves[0],
            Some(rectangle.id),
        ),
    ] {
        s.click(p, offset);
        let selection = s.selection();
        assert_eq!(selection.entities, vec![entity], "{what}: {}", s.message());
        assert_eq!(selection.element, element, "{what} in the element list");
        assert!(s.highlighted_at(p), "{what} visibly highlighted");
    }
    // The 2 mm default fillet spans a few pixels at this zoom; zoomed in as
    // a user would, its arc is picked rather than the virtual sharp beside it.
    let view = s.h.app.current().sketch.camera;
    {
        let camera = &mut s.h.app.current_mut().sketch.camera;
        camera.target = [fillet_mid[0] as f32, fillet_mid[1] as f32, 0.];
        camera.scale *= 8.;
    }
    s.settle();
    s.click(fillet_mid, egui::Vec2::ZERO);
    let selection = s.selection();
    assert_eq!(selection.entities, vec![fillet], "fillet: {}", s.message());
    assert_eq!(
        selection.element,
        Some(fillet_shape.id),
        "fillet in the element list"
    );
    assert!(s.highlighted_at(fillet_mid), "fillet visibly highlighted");
    s.h.app.current_mut().sketch.camera = view;
    s.settle();
    // Shift adds, Cmd removes again, a plain click replaces.
    s.click([30., 0.], ON_LINE);
    s.click_with([30., 60.], ON_LINE, SHIFT);
    assert_eq!(
        s.selection().entities,
        vec![rectangle.curves[0], line.curves[0]]
    );
    assert!(s.highlighted_at([30., 0.]) && s.highlighted_at([30., 60.]));
    s.click_with([30., 0.], ON_LINE, COMMAND);
    assert_eq!(s.selection().entities, vec![line.curves[0]]);
    assert!(!s.highlighted_at([30., 0.]));
    s.click_with(arc_mid, egui::Vec2::ZERO, COMMAND);
    assert_eq!(s.selection().entities, vec![line.curves[0], arc.curves[0]]);
    // Empty space clears everything.
    s.click([30., 50.], egui::Vec2::ZERO);
    let cleared = s.selection();
    assert!(cleared.entities.is_empty(), "{cleared:?}");
    assert_eq!(cleared.element, None);
    assert!(cleared.highlight.is_empty());
    assert!(!s.highlighted_at([30., 60.]));
    assert_eq!(s.h.revision(), revision, "selection never commits");
    assert_eq!(s.h.app.current().session.document(), &document);
}

/// Canvas → list → Inspector and list → canvas stay in sync.
#[test]
fn canvas_selection_the_element_list_and_the_inspector_stay_in_sync() {
    let mut s = SketchHarness::new();
    s.draw(DrawTool::Rectangle, &["0", "0", "60", "40"]);
    s.draw(DrawTool::Line, &["0", "60", "60", "60"]);
    let source = s.source();
    let (rectangle, line) = (
        source.sketch.shapes[0].clone(),
        source.sketch.shapes[1].clone(),
    );
    let revision = s.h.revision();
    // Canvas → list and Inspector.
    s.click([60., 20.], egui::vec2(-1., 0.));
    let selection = s.selection();
    assert_eq!(selection.element, Some(rectangle.id));
    assert_eq!(
        selection.highlight,
        rectangle.curves.iter().copied().collect()
    );
    let texts = s.h.texts();
    assert!(
        texts.iter().any(|t| t == "Selected: Rectangle 1"),
        "{texts:?}"
    );
    for label in ["Width", "Height"] {
        assert!(texts.iter().any(|t| t == label), "{label}: {texts:?}");
    }
    // List → canvas and Inspector.
    s.h.click_label("Line 1");
    s.settle();
    let selection = s.selection();
    assert_eq!(selection.element, Some(line.id));
    assert_eq!(selection.highlight, line.curves.iter().copied().collect());
    assert!(
        s.highlighted_at([30., 60.]),
        "list choice highlights in the canvas"
    );
    assert!(!s.highlighted_at([30., 0.]));
    let texts = s.h.texts();
    assert!(texts.iter().any(|t| t == "Selected: Line 1"), "{texts:?}");
    assert!(texts.iter().any(|t| t == "Length"), "{texts:?}");
    // Canvas again replaces the list choice.
    s.click([30., 0.], ON_LINE);
    assert_eq!(s.selection().element, Some(rectangle.id));
    assert_eq!(s.h.revision(), revision);
}

/// Esc is the other way back to Select; the sketch's origin does not matter:
/// Add Sketch on a datum pick, on a body face, or a reopened sketch.
#[test]
fn select_works_after_esc_in_picked_face_and_reopened_sketches() {
    // Datum pick.
    let mut h = Harness::new(Document::new("Datum pick").unwrap());
    h.app.invoke(Command::NewSketch, &[], Channel::Ui).unwrap();
    h.frame(vec![]);
    let sketch = &h.app.current().sketch;
    let p = sketch
        .camera
        .project([10., 10., 0.], sketch.viewport_rect.unwrap());
    h.click_at(egui::pos2(p[0], p[1]));
    assert!(h.app.current().sketch.editing, "opened on the datum");
    let mut s = SketchHarness { h };
    s.settle();
    s.tool(DrawTool::Rectangle);
    s.click([0., 0.], egui::Vec2::ZERO);
    s.click([60., 40.], egui::Vec2::ZERO);
    let rectangle = s.source().sketch.shapes[0].clone();
    s.h.key(egui::Key::Escape);
    assert!(s.h.app.current().sketch.select_active());
    s.click([30., 0.], ON_LINE);
    edge_selected(&s, rectangle.curves[0], "datum pick");
    // Reopened.
    s.h.app
        .invoke(Command::FinishSketch, &[], Channel::Ui)
        .unwrap();
    for _ in 0..3 {
        s.h.frame(vec![]);
    }
    // CAD-425 (approved by the user) removed the "Open selected sketch"
    // button; reopen through Edit sketch in the tree's context menu.
    let id =
        s.h.app
            .current()
            .sketch
            .selected
            .expect("the finished sketch");
    let name =
        s.h.app
            .current()
            .session
            .document()
            .get(id)
            .unwrap()
            .name
            .clone();
    s.h.right_click_label(&name);
    s.h.frame(vec![]);
    s.h.click_label("Edit sketch");
    s.settle();
    s.click([60., 20.], egui::vec2(-1., 0.));
    edge_selected(&s, rectangle.curves[1], "reopened");

    // Body face.
    let mut session = Session::new(Document::new("Face sketch").unwrap()).unwrap();
    let request = multicad_core::model::PrimitiveRequest::new(
        "Stock",
        multicad_core::model::PrimitiveParameters::Box {
            size: [80., 60., 10.].map(multicad_core::model::length),
            frame: Default::default(),
        },
    );
    session
        .commit(
            multicad_core::model::prepare_primitive(
                &session,
                0,
                Actor::local(Channel::Ui),
                &request,
            )
            .unwrap()
            .prepared,
        )
        .unwrap();
    let mut h = Harness::new(session.document().clone());
    h.app.invoke(Command::NewSketch, &[], Channel::Ui).unwrap();
    h.frame(vec![]);
    let sketch = &h.app.current().sketch;
    let p = sketch
        .camera
        .project([40., 30., 10.], sketch.viewport_rect.unwrap());
    h.click_at(egui::pos2(p[0], p[1]));
    assert!(h.app.current().sketch.editing, "opened on the face");
    let mut s = SketchHarness { h };
    s.settle();
    s.tool(DrawTool::Rectangle);
    s.click([10., 10.], egui::Vec2::ZERO);
    s.click([50., 40.], egui::Vec2::ZERO);
    let rectangle = s.source().sketch.shapes.last().unwrap().clone();
    s.h.key(egui::Key::Escape);
    s.click([30., 10.], ON_LINE);
    edge_selected(&s, rectangle.curves[0], "face sketch");
}

/// Real-time input: the pointer travels across frames and the button is held
/// for 100 ms, as a physical click does.
#[test]
fn a_click_with_real_pointer_travel_and_timing_selects() {
    let mut s = SketchHarness::new();
    s.draw(DrawTool::Rectangle, &["0", "0", "60", "40"]);
    let rectangle = s.source().sketch.shapes[0].clone();
    let target = s.click_point([30., 0.]) + ON_LINE;
    let start = target + egui::vec2(-200., -150.);
    let mut t = 10.;
    for i in 0..=30 {
        s.timed(
            vec![egui::Event::PointerMoved(
                start + (target - start) * (i as f32 / 30.),
            )],
            t,
        );
        t += 1. / 60.;
    }
    let button = |pressed| egui::Event::PointerButton {
        pos: target,
        button: egui::PointerButton::Primary,
        pressed,
        modifiers: egui::Modifiers::NONE,
    };
    s.timed(vec![button(true)], t);
    for _ in 0..6 {
        t += 1. / 60.;
        s.timed(vec![], t);
    }
    s.timed(vec![button(false)], t + 1. / 60.);
    s.settle();
    edge_selected(&s, rectangle.curves[0], "timed click");
}

/// Toggling keeps the element list on what is still selected.
#[test]
fn toggling_keeps_the_element_list_on_the_remaining_selection() {
    let mut s = SketchHarness::new();
    s.draw(DrawTool::Rectangle, &["0", "0", "60", "40"]);
    s.draw(DrawTool::Line, &["0", "60", "60", "60"]);
    let source = s.source();
    let (rectangle, line) = (
        source.sketch.shapes[0].clone(),
        source.sketch.shapes[1].clone(),
    );
    let corner = point_at(&source, [60., 60.]);
    s.click([30., 0.], ON_LINE);
    s.click_with([30., 60.], ON_LINE, SHIFT);
    assert_eq!(s.selection().element, Some(line.id));
    // Removing the line returns the list to the rectangle edge still selected.
    s.click_with([30., 60.], ON_LINE, COMMAND);
    let selection = s.selection();
    assert_eq!(selection.entities, vec![rectangle.curves[0]]);
    assert_eq!(selection.element, Some(rectangle.id));
    // Adding a point keeps it there; removing the edge empties the list.
    s.click_with([60., 60.], egui::Vec2::ZERO, SHIFT);
    let selection = s.selection();
    assert_eq!(selection.entities, vec![rectangle.curves[0], corner]);
    assert_eq!(selection.element, Some(rectangle.id));
    s.click_with([30., 0.], ON_LINE, COMMAND);
    let selection = s.selection();
    assert_eq!(selection.entities, vec![corner]);
    assert_eq!(selection.element, None);
    assert!(selection.highlight.is_empty());
}

/// The product owner's sequence (CAD-412 answers): a sketch opens in Select;
/// draw a rectangle, then a line, and click an edge without pressing Esc. The
/// drawing tool has already ended, so the edge selects. With a tool that keeps
/// going active, the status line names it, so a click that draws is explained.
#[test]
fn after_drawing_without_esc_a_click_on_an_edge_selects_it() {
    let mut s = SketchHarness::new();
    assert!(
        s.h.app.current().sketch.select_active(),
        "a new sketch opens in Select"
    );
    s.tool(DrawTool::Rectangle);
    s.click([0., 0.], egui::Vec2::ZERO);
    s.click([60., 40.], egui::Vec2::ZERO);
    let rectangle = s.source().sketch.shapes[0].clone();
    s.click([30., 0.], ON_LINE);
    edge_selected(&s, rectangle.curves[0], "rectangle edge after drawing");
    assert!(s.highlighted_at([30., 0.]));
    s.tool(DrawTool::Line);
    s.click([0., 60.], egui::Vec2::ZERO);
    s.click([60., 60.], egui::Vec2::ZERO);
    let line = s.source().sketch.shapes[1].clone();
    s.click([30., 60.], ON_LINE);
    edge_selected(&s, line.curves[0], "line after drawing");
    // Polyline keeps going: the status line names it.
    s.tool(DrawTool::Polyline);
    assert!(s.message().starts_with("Polyline · "), "{}", s.message());
}
