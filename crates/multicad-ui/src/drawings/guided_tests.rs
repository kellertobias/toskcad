//! CAD-436: the Drawings empty state, one New drawing action per context and
//! guided Add view, through the real editor frame with headless input. The
//! test provider projects an axis-aligned box onto each view's frame, so
//! orientation, extents and scale follow real projection rules.
use super::move_tests::{click_at, context, controls, frame, run, texts};
use super::*;
use multicad_core::{
    document::{Document, Entity, EntityData},
    drawings::{
        DrawingEvaluation, DrawingStamp, DrawingVector, VectorCurve, VectorRole, ViewEvaluation,
    },
};
use std::sync::Arc;

/// Box corners per source part: min and max in model mm.
pub(crate) type Boxes = Arc<std::collections::BTreeMap<EntityId, ([f64; 3], [f64; 3])>>;
struct BoxPublisher {
    boxes: Boxes,
    reply: Option<DrawingReply>,
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
impl DrawingBackend for BoxPublisher {
    fn request(&mut self, document: &Document, stamp: DrawingStamp) -> Result<(), String> {
        let mut views = std::collections::BTreeMap::new();
        for id in drawings::books(document) {
            for (page, view) in drawings::read(document, id)
                .unwrap()
                .pages
                .iter()
                .flat_map(|p| p.views.iter().map(move |v| (p, v)))
            {
                let mut vectors = vec![];
                let filter = page.part_filter(view.id);
                for source in view.source.entities() {
                    // Like the native provider, hidden parts are not drawn.
                    if filter.is_some_and(|f| !f.shows(document, source)) {
                        continue;
                    }
                    let Some((lo, hi)) = self.boxes.get(&source) else {
                        continue;
                    };
                    let o = view.frame.origin_mm;
                    let (u, v) = (view.frame.u, view.frame.v);
                    let corners: Vec<[f64; 2]> = (0..8)
                        .map(|i| {
                            let p = [
                                if i & 1 == 0 { lo[0] } else { hi[0] },
                                if i & 2 == 0 { lo[1] } else { hi[1] },
                                if i & 4 == 0 { lo[2] } else { hi[2] },
                            ];
                            let d = [p[0] - o[0], p[1] - o[1], p[2] - o[2]];
                            [dot(d, u), dot(d, v)]
                        })
                        .collect();
                    // The silhouette rectangle of the projected corners.
                    let (mut min, mut max) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
                    for c in &corners {
                        for i in 0..2 {
                            min[i] = min[i].min(c[i]);
                            max[i] = max[i].max(c[i]);
                        }
                    }
                    let rect = [min, [max[0], min[1]], max, [min[0], max[1]], min];
                    for pair in rect.windows(2) {
                        vectors.push(DrawingVector {
                            color_rgb: None,
                            source,
                            material_part: source,
                            role: VectorRole::Visible,
                            curve: VectorCurve::Line {
                                a: pair[0],
                                b: pair[1],
                            },
                            preview_mm: pair.to_vec(),
                            line_width_mm: 0.35,
                        });
                    }
                }
                let evaluation = if vectors.is_empty() {
                    ViewEvaluation::Unresolved {
                        message: "No test geometry for this source".into(),
                    }
                } else {
                    ViewEvaluation::Ready {
                        labels: vec![],
                        vectors,
                        warnings: vec![],
                    }
                };
                views.insert(view.id, evaluation);
            }
        }
        self.reply = Some(DrawingReply {
            stamp,
            result: Ok(Arc::new(DrawingEvaluation {
                stamp,
                views,
                annotations: Default::default(),
            })),
        });
        Ok(())
    }
    fn poll(&mut self) -> Option<DrawingReply> {
        self.reply.take()
    }
    fn cancel(&mut self) {
        self.reply = None;
    }
}
pub(crate) fn box_factory(boxes: Boxes) -> DrawingBackendFactory {
    Arc::new(move || {
        Box::new(BoxPublisher {
            boxes: boxes.clone(),
            reply: None,
        })
    })
}
/// A document with parts of these box sizes (x, y, z), side by side on x.
pub(crate) fn boxed_document(sizes: &[[f64; 3]]) -> (Document, Vec<EntityId>, Boxes) {
    tagged_document(&sizes.iter().map(|s| (*s, &[][..])).collect::<Vec<_>>())
}
/// Parts of these box sizes with these tags, side by side on x.
pub(crate) fn tagged_document(parts: &[([f64; 3], &[&str])]) -> (Document, Vec<EntityId>, Boxes) {
    let mut doc = Document::new("Guided views").unwrap();
    let mut boxes = std::collections::BTreeMap::new();
    let mut ids = vec![];
    let mut x = 0.;
    for (i, (size, tags)) in parts.iter().enumerate() {
        let part = Entity::new(
            format!("Block {}", i + 1),
            EntityData::Part {
                features: vec![],
                material: None,
                tags: tags.iter().map(|t| t.to_string()).collect(),
            },
        );
        boxes.insert(part.id, ([x, 0., 0.], [x + size[0], size[1], size[2]]));
        x += size[0] + 50.;
        ids.push(part.id);
        doc.insert(part).unwrap();
    }
    (doc, ids, Arc::new(boxes))
}
pub(super) fn editor(sizes: &[[f64; 3]]) -> (DrawingEditor, Session, Vec<EntityId>) {
    let (doc, ids, boxes) = boxed_document(sizes);
    let mut editor = DrawingEditor::default();
    editor.set_backend_factory(box_factory(boxes));
    (editor, Session::new(doc).unwrap(), ids)
}
/// Run frames until `done`, with a frame budget instead of a fixed count.
pub(super) fn until(
    ctx: &egui::Context,
    editor: &mut DrawingEditor,
    session: &mut Session,
    what: &str,
    done: impl Fn(&DrawingEditor, &Session) -> bool,
) {
    for _ in 0..400 {
        if done(editor, session) {
            return;
        }
        frame(ctx, editor, session, vec![]);
        editor.feed.tick(1e6);
    }
    panic!("Timed out waiting for {what}");
}
pub(super) fn saved(session: &Session, editor: &DrawingEditor) -> drawings::Page {
    drawings::read(session.document(), editor.book.unwrap())
        .unwrap()
        .pages
        .into_iter()
        .find(|p| Some(p.id) == editor.selected)
        .unwrap()
}
pub(super) fn screen(editor: &DrawingEditor, mm: [f64; 2]) -> egui::Pos2 {
    let (paper, scale) = editor.last_paper.expect("paper drawn");
    paper.min + egui::vec2(mm[0] as f32, mm[1] as f32) * scale
}
/// Painted line segments, in screen points.
fn segments(output: &egui::FullOutput) -> Vec<[egui::Pos2; 2]> {
    output
        .shapes
        .iter()
        .filter_map(|s| match &s.shape {
            egui::Shape::Path(path) if path.points.len() == 2 => {
                Some([path.points[0], path.points[1]])
            }
            egui::Shape::LineSegment { points, .. } => Some(*points),
            _ => None,
        })
        .collect()
}
fn new_drawing_controls(output: &egui::FullOutput) -> Vec<(bool, egui::Pos2)> {
    controls(output, "New drawing")
}
fn bounds_of(output: &egui::FullOutput, label: &str) -> egui::Rect {
    output
        .platform_output
        .accesskit_update
        .as_ref()
        .unwrap()
        .nodes
        .iter()
        .find_map(|(_, n)| {
            (n.label() == Some(label))
                .then(|| n.bounds())
                .flatten()
                .map(|b| {
                    egui::Rect::from_min_max(
                        egui::pos2(b.x0 as f32, b.y0 as f32),
                        egui::pos2(b.x1 as f32, b.y1 as f32),
                    )
                })
        })
        .unwrap_or_else(|| panic!("Missing {label}"))
}

#[test]
fn empty_state_is_centred_with_one_large_new_drawing_action() {
    let (mut editor, mut session, _) = editor(&[[100., 40., 60.]]);
    let ctx = context();
    for _ in 0..2 {
        frame(&ctx, &mut editor, &mut session, vec![]);
    }
    let output = run(&ctx, &mut editor, &mut session, vec![]);
    let shown = texts(&output);
    assert!(shown.iter().any(|t| t == "No drawing pages"));
    // Exactly one New drawing control: the empty state's, large and centred
    // in the sheet area (right of the Drawings list, left of the Inspector).
    assert_eq!(
        new_drawing_controls(&output).len(),
        1,
        "One New drawing action"
    );
    let button = bounds_of(&output, "New drawing");
    assert!(
        button.height() >= 44. && button.width() >= 200.,
        "{button:?}"
    );
    let panel = |id: &str| {
        egui::PanelState::load(&ctx, egui::Id::new(id))
            .unwrap_or_else(|| panic!("{id} panel"))
            .outer_rect
    };
    let area_center = (panel("drawing-sidebar").right() + panel("drawing-inspector").left()) * 0.5;
    assert!(
        (button.center().x - area_center).abs() < 40.,
        "Button centred horizontally: {button:?} vs {area_center}"
    );
    assert!(
        (button.center().y - 425.).abs() < 80.,
        "Centred vertically: {button:?}"
    );
    // Without a drawing, Add page is not a second way to start one.
    assert!(!editor.enabled(DrawingTool::AddPage, &session));
    let at = button.center();
    output.drop_without_applying_deltas();
    click_at(&ctx, &mut editor, &mut session, at);
    let book = editor.book.expect("New drawing created a drawing");
    assert_eq!(
        drawings::read(session.document(), book)
            .unwrap()
            .pages
            .len(),
        1
    );
    // With a drawing, the Drawings list keeps the one New drawing action.
    let output = run(&ctx, &mut editor, &mut session, vec![]);
    assert!(!texts(&output).iter().any(|t| t == "No drawing pages"));
    assert_eq!(new_drawing_controls(&output).len(), 1);
    assert!(
        texts(&output)
            .iter()
            .any(|t| t.starts_with("Next: Add view places a front view")),
        "The sheet hints the next step"
    );
    output.drop_without_applying_deltas();
}

#[test]
fn the_drawings_ribbon_offers_no_second_new_drawing_and_guided_view_tools() {
    let registry = crate::registry::Registry::default();
    let workspace = registry.workspace(crate::registry::Mode::Drawings);
    let mut shown = vec![];
    for tab in &workspace.tabs {
        for group in crate::ribbon_ui::groups(crate::registry::Mode::Drawings, tab) {
            for entry in group.entries {
                if let crate::ribbon_ui::Entry::Tool(id) = entry {
                    shown.push(id);
                }
            }
        }
    }
    for hidden in [
        "drawing-new",
        "page-apply",
        "page-discard",
        "drawing-projection",
    ] {
        assert!(!shown.contains(&hidden), "{hidden} is not a ribbon tile");
    }
    for id in ["drawing-add-view", "drawing-standard-views", "sheet-page"] {
        assert!(shown.contains(&id), "{id} is a ribbon tile");
    }
}

#[test]
fn add_view_defaults_to_a_front_view_of_the_whole_model_that_follows_the_pointer() {
    // 200 × 50 × 80 mm: the front view is 200 wide and 80 high.
    let (mut editor, mut session, parts) = editor(&[[200., 50., 80.], [40., 30., 20.]]);
    editor
        .invoke(&mut session, DrawingTool::NewDocument, Channel::Ui)
        .unwrap();
    let ctx = context();
    let revision = session.revision();
    editor
        .invoke(&mut session, DrawingTool::AddView, Channel::Ui)
        .unwrap();
    until(
        &ctx,
        &mut editor,
        &mut session,
        "the fitted preview",
        |e, _| {
            e.next_step_hint()
                .is_some_and(|h| h.starts_with("Click the sheet to place the front view"))
        },
    );
    assert!(
        editor
            .next_step_hint()
            .unwrap()
            .contains("the whole model (2 parts)")
    );
    // Whole model: both blocks, 290 + 40 = 330 mm wide (with the gap) and
    // 80 mm high. A4 landscape fits 1:2 at most.
    let page = editor.draft.as_ref().unwrap();
    let view = page.views.last().unwrap().clone();
    assert_eq!(view.frame.u, [1., 0., 0.]);
    assert_eq!(view.frame.v, [0., 0., 1.], "Front looks along +Y with Z up");
    let layout = drawings::layout_page(page, 1, 1).unwrap();
    let expected = drawings::fit_view(
        &view,
        &layout,
        [0., 0.],
        [290., 80.],
        drawings::FitMode::Standard,
    )
    .unwrap()
    .scale;
    assert_eq!(view.scale, expected);
    assert_eq!(drawings::scale_text(view.scale), "1:2");
    let mut shown = view.source.entities();
    shown.sort();
    let mut all = parts.clone();
    all.sort();
    assert_eq!(shown, all, "The whole model");
    assert_eq!(
        session.revision(),
        revision,
        "Nothing applied while placing"
    );
    // The preview follows the pointer: its outline is drawn around it.
    let pointer = [120., 80.];
    let at = screen(&editor, pointer);
    frame(
        &ctx,
        &mut editor,
        &mut session,
        vec![egui::Event::PointerMoved(at)],
    );
    let output = run(
        &ctx,
        &mut editor,
        &mut session,
        vec![egui::Event::PointerMoved(at)],
    );
    let half = egui::vec2(145. * 0.5, 40. * 0.5) * editor.last_paper.unwrap().1;
    let corner = at + egui::vec2(-half.x, half.y);
    assert!(
        segments(&output)
            .iter()
            .any(|[a, b]| a.distance(corner) < 1.5 || b.distance(corner) < 1.5),
        "The view outline follows the pointer"
    );
    output.drop_without_applying_deltas();
    // One click places it; the saved page then holds the view there.
    click_at(&ctx, &mut editor, &mut session, at);
    until(
        &ctx,
        &mut editor,
        &mut session,
        "the applied view",
        |e, _| !e.view_placement_active(),
    );
    assert_eq!(session.revision(), revision + 1, "One undoable edit");
    let saved = saved(&session, &editor);
    let placed = saved.views.last().unwrap();
    assert!((placed.position_mm[0] - pointer[0]).abs() < 0.5);
    assert!((placed.position_mm[1] - pointer[1]).abs() < 0.5);
    assert_eq!(placed.scale, expected);
    assert!(!editor.dirty());
    let shown = frame(&ctx, &mut editor, &mut session, vec![]);
    assert!(
        shown
            .iter()
            .any(|t| t.starts_with("Front view placed at 1:2"))
    );
}

#[test]
fn add_view_shows_the_selection_and_offers_one_click_orientation_and_scale() {
    let (mut editor, mut session, parts) = editor(&[[200., 50., 80.], [40., 30., 20.]]);
    editor
        .invoke(&mut session, DrawingTool::NewDocument, Channel::Ui)
        .unwrap();
    editor.set_model_selection(vec![parts[1]]);
    let ctx = context();
    editor
        .invoke(&mut session, DrawingTool::AddView, Channel::Ui)
        .unwrap();
    until(
        &ctx,
        &mut editor,
        &mut session,
        "the fitted preview",
        |e, _| {
            e.next_step_hint()
                .is_some_and(|h| h.starts_with("Click the sheet"))
        },
    );
    assert!(
        editor
            .next_step_hint()
            .unwrap()
            .contains("the selection (1 part)")
    );
    // One click on Top re-orients and refits the pending view.
    let output = run(&ctx, &mut editor, &mut session, vec![]);
    let top = controls(&output, "Top");
    assert_eq!(top.len(), 1);
    for label in ["Front", "Right", "Iso"] {
        assert_eq!(controls(&output, label).len(), 1, "{label} is one click");
    }
    output.drop_without_applying_deltas();
    click_at(&ctx, &mut editor, &mut session, top[0].1);
    until(
        &ctx,
        &mut editor,
        &mut session,
        "the refitted preview",
        |e, _| {
            e.next_step_hint()
                .is_some_and(|h| h.starts_with("Click the sheet to place the top view"))
        },
    );
    let view = editor.draft.as_ref().unwrap().views.last().unwrap().clone();
    assert_eq!((view.frame.u, view.frame.v), StandardView::Top.frame());
    assert_eq!(
        view.source,
        drawings::ViewSource::Model {
            entities: vec![parts[1]]
        }
    );
    // A small 40 × 30 part is enlarged to the largest standard scale that
    // fits this sheet.
    let page = editor.draft.as_ref().unwrap();
    let layout = drawings::layout_page(page, 1, 1).unwrap();
    let fit = drawings::fit_view(
        &view,
        &layout,
        [0., 0.],
        [40., 30.],
        drawings::FitMode::Standard,
    )
    .unwrap();
    assert_eq!(view.scale, fit.scale);
    assert!(drawings::scale_factor(view.scale) > 1.);
    // The user can override the scale before placing.
    editor.set_placement_scale(drawings::Scale {
        numerator: 1,
        denominator: 1,
    });
    until(
        &ctx,
        &mut editor,
        &mut session,
        "the overridden scale",
        |e, _| {
            e.next_step_hint()
                .is_some_and(|h| h.starts_with("Click the sheet"))
        },
    );
    editor.place_view([100., 90.]).unwrap();
    until(
        &ctx,
        &mut editor,
        &mut session,
        "the applied view",
        |e, _| !e.view_placement_active(),
    );
    let placed = saved(&session, &editor).views.last().unwrap().clone();
    assert_eq!(drawings::scale_text(placed.scale), "1:1");
    assert_eq!((placed.frame.u, placed.frame.v), StandardView::Top.frame());
}

#[test]
fn escape_cancels_add_view_and_keeps_the_page() {
    let (mut editor, mut session, _) = editor(&[[200., 50., 80.]]);
    editor
        .invoke(&mut session, DrawingTool::NewDocument, Channel::Ui)
        .unwrap();
    let before = session.document().clone();
    let ctx = context();
    editor
        .invoke(&mut session, DrawingTool::AddView, Channel::Ui)
        .unwrap();
    assert!(!editor.enabled(DrawingTool::AddView, &session));
    frame(&ctx, &mut editor, &mut session, vec![]);
    frame(
        &ctx,
        &mut editor,
        &mut session,
        vec![egui::Event::Key {
            key: egui::Key::Escape,
            physical_key: None,
            pressed: true,
            repeat: false,
            modifiers: Default::default(),
        }],
    );
    assert!(!editor.view_placement_active());
    assert!(!editor.dirty());
    assert!(editor.draft.as_ref().unwrap().views.is_empty());
    assert_eq!(session.document(), &before);
    // An empty model cannot add a view.
    let mut empty = DrawingEditor::default();
    let mut s = Session::new(Document::new("Empty").unwrap()).unwrap();
    empty.set_backend_factory(box_factory(Default::default()));
    empty
        .invoke(&mut s, DrawingTool::NewDocument, Channel::Ui)
        .unwrap();
    assert!(
        empty
            .invoke(&mut s, DrawingTool::AddView, Channel::Ui)
            .unwrap_err()
            .contains("Add a part")
    );
}

#[test]
fn standard_views_place_front_top_right_and_iso_aligned_in_one_click() {
    // 300 × 100 × 120: front 300 × 120, top 300 × 100, right 100 × 120.
    let (mut editor, mut session, _) = editor(&[[300., 100., 120.]]);
    editor
        .invoke(&mut session, DrawingTool::NewDocument, Channel::Ui)
        .unwrap();
    let revision = session.revision();
    let ctx = context();
    editor
        .invoke(&mut session, DrawingTool::StandardViews, Channel::Ui)
        .unwrap();
    until(
        &ctx,
        &mut editor,
        &mut session,
        "the arranged views",
        |e, _| !e.view_placement_active(),
    );
    assert_eq!(session.revision(), revision + 1, "One undoable edit");
    let page = saved(&session, &editor);
    let names: Vec<_> = page.views.iter().map(|v| v.name.as_str()).collect();
    assert_eq!(names, ["Front view", "Top view", "Right view", "Iso view"]);
    let [front, top, right, iso] = [0, 1, 2, 3].map(|i| page.views[i].clone());
    assert_eq!(page.projections.len(), 2);
    assert!(page.projections.iter().all(|l| l.base == front.id));
    // All four at one scale; top in the front's column, right in its row,
    // with a shared model center so model points align across them.
    assert!(page.views.iter().all(|v| v.scale == front.scale));
    assert_eq!(top.position_mm[0], front.position_mm[0]);
    assert_eq!(right.position_mm[1], front.position_mm[1]);
    assert_eq!(top.frame.origin_mm, front.frame.origin_mm);
    assert_eq!(right.frame.origin_mm, front.frame.origin_mm);
    assert_eq!(front.frame.origin_mm, [150., 50., 60.]);
    // First angle (the default): top below the front, right view on its left.
    assert!(top.position_mm[1] > front.position_mm[1]);
    assert!(right.position_mm[0] < front.position_mm[0]);
    assert_eq!((iso.frame.u, iso.frame.v), StandardView::Iso.frame());
    // The scale is the standard layout's choice for this sheet.
    let layout = drawings::layout_page(&page, 1, 1).unwrap();
    let iso_extent = {
        let e = editor.feed.result.as_ref().unwrap();
        let Some(ViewEvaluation::Ready { vectors, .. }) = e.views.get(&iso.id) else {
            panic!("iso ready")
        };
        let (mut lo, mut hi) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
        for p in vectors.iter().flat_map(|v| &v.preview_mm) {
            for i in 0..2 {
                lo[i] = lo[i].min(p[i]);
                hi[i] = hi[i].max(p[i]);
            }
        }
        [hi[0] - lo[0], hi[1] - lo[1]]
    };
    let arranged = drawings::standard_layout(
        &layout,
        page.settings.projection,
        [[300., 120.], [300., 100.], [100., 120.], iso_extent],
    )
    .unwrap();
    assert_eq!(front.scale, arranged.scale);
    assert_eq!(front.position_mm, arranged.positions[0]);
    let shown = frame(&ctx, &mut editor, &mut session, vec![]);
    assert!(
        shown
            .iter()
            .any(|t| t.starts_with("Front, top, right and iso views placed"))
    );
    // Dragging the front keeps the projected views aligned.
    let page = editor.draft.as_mut().unwrap();
    drawings::move_view(page, front.id, [5., 3.]).unwrap();
    assert_eq!(page.views[1].position_mm[0], page.views[0].position_mm[0]);
    assert_eq!(page.views[2].position_mm[1], page.views[0].position_mm[1]);
}
