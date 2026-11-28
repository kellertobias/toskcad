//! CAD-425: a finished sketch reopens by double-click in the model tree or on
//! its geometry in the 3D view, or by Edit sketch in the tree's context menu.
//! There is no "Open selected sketch" button.
use super::shape_list_tests::Host;
use super::*;

/// A finished 40 × 20 rectangle sketch, selected but not being edited.
fn finished() -> (Host, EntityId, String) {
    let mut h = Host::new();
    h.editor
        .add(
            &mut h.session,
            DrawTool::Rectangle,
            [0., 0.],
            [40., 20.],
            Channel::Ui,
        )
        .unwrap();
    h.settle();
    let id = h.sketch();
    h.editor.finish().unwrap();
    h.settle();
    assert!(!h.editor.editing);
    let name = h.session.document().get(id).unwrap().name.clone();
    (h, id, name)
}
fn label_center(h: &mut Host, label: &str) -> Pos2 {
    let output = h.frame(vec![]);
    let at = output
        .shapes
        .iter()
        .find_map(|shape| match &shape.shape {
            egui::Shape::Text(text) if text.galley.text() == label => {
                Some(text.pos + text.galley.size() * 0.5)
            }
            _ => None,
        })
        .unwrap_or_else(|| panic!("Missing {label}"));
    output.drop_without_applying_deltas();
    at
}
fn button_press(h: &mut Host, at: Pos2, button: egui::PointerButton, pressed: bool) {
    h.frame(vec![
        egui::Event::PointerMoved(at),
        egui::Event::PointerButton {
            pos: at,
            button,
            pressed,
            modifiers: egui::Modifiers::NONE,
        },
    ])
    .drop_without_applying_deltas();
}
fn double_click(h: &mut Host, at: Pos2) {
    h.frame(vec![egui::Event::PointerMoved(at)])
        .drop_without_applying_deltas();
    for pressed in [true, false, true, false] {
        button_press(h, at, egui::PointerButton::Primary, pressed);
    }
    h.settle();
}
fn assert_editing(h: &Host, id: EntityId, document: &Document, revision: u64) {
    assert!(h.editor.editing, "{}", h.editor.message);
    assert_eq!(h.editor.selected, Some(id));
    assert_eq!(h.session.document(), document, "opening changes no source");
    assert_eq!(h.session.revision(), revision, "opening adds no history");
}

#[test]
fn the_open_selected_sketch_button_is_gone() {
    let (mut h, _, _) = finished();
    let labels = h.labels();
    assert!(
        !labels.iter().any(|l| l.contains("Open selected sketch")),
        "{labels:?}"
    );
}

#[test]
fn double_clicking_the_sketch_in_the_model_tree_opens_it() {
    let (mut h, id, name) = finished();
    let (document, revision) = (h.session.document().clone(), h.session.revision());
    let at = label_center(&mut h, &name);
    // One click only selects it.
    h.press(at);
    assert!(!h.editor.editing);
    h.frame_at(vec![], h.time + 1.)
        .drop_without_applying_deltas();
    double_click(&mut h, at);
    assert_editing(&h, id, &document, revision);
}

#[test]
fn edit_sketch_in_the_tree_context_menu_opens_it() {
    let (mut h, id, name) = finished();
    let (document, revision) = (h.session.document().clone(), h.session.revision());
    let at = label_center(&mut h, &name);
    button_press(&mut h, at, egui::PointerButton::Secondary, true);
    button_press(&mut h, at, egui::PointerButton::Secondary, false);
    h.frame(vec![]).drop_without_applying_deltas();
    assert!(!h.editor.editing);
    h.click_label("Edit sketch");
    assert_editing(&h, id, &document, revision);
}

#[test]
fn double_clicking_sketch_geometry_in_the_3d_view_opens_it() {
    let (mut h, id, _) = finished();
    let (document, revision) = (h.session.document().clone(), h.session.revision());
    // Empty space: nothing opens.
    let empty = h.screen([-80., -80.]);
    double_click(&mut h, empty);
    assert!(!h.editor.editing, "{}", h.editor.message);
    h.frame_at(vec![], h.time + 1.)
        .drop_without_applying_deltas();
    // The rectangle's bottom edge.
    let edge = h.screen([20., 0.]) + Vec2::new(0., 1.);
    double_click(&mut h, edge);
    assert_editing(&h, id, &document, revision);
}
