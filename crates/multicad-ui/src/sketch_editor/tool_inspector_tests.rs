//! Contextual sketch Inspector (CAD-399): the active tool's settings at the
//! top and nothing else, or the selected element's values, in the real host.
use super::shape_list_tests::{Host, key};
use super::*;

/// Painted texts with their top-left positions.
fn texts(h: &mut Host) -> Vec<(String, Pos2)> {
    let output = h.frame(vec![]);
    let texts = output
        .shapes
        .iter()
        .filter_map(|shape| match &shape.shape {
            egui::Shape::Text(text) => Some((text.galley.text().to_string(), text.pos)),
            _ => None,
        })
        .collect();
    output.drop_without_applying_deltas();
    texts
}
/// Texts painted in the right Inspector panel.
fn inspector(h: &mut Host) -> Vec<(String, Pos2)> {
    let right = h.editor.viewport_rect.unwrap().right();
    texts(h).into_iter().filter(|(_, p)| p.x >= right).collect()
}
fn find<'a>(texts: &'a [(String, Pos2)], label: &str) -> Option<&'a Pos2> {
    texts.iter().find(|(t, _)| t == label).map(|(_, p)| p)
}
fn with_rectangle() -> Host {
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
    h
}

#[test]
fn the_fillet_radius_is_the_first_setting_and_is_editable_while_fillet_is_active() {
    let mut h = with_rectangle();
    h.editor.set_tool(DrawTool::Fillet).unwrap();
    let shown = inspector(&mut h);
    let radius = *find(&shown, "Fillet radius").unwrap_or_else(|| panic!("{shown:?}"));
    // Only the Inspector heading and the tool title precede it.
    let above: Vec<_> = shown
        .iter()
        .filter(|(_, p)| p.y < radius.y - 1.)
        .map(|(t, _)| t.as_str())
        .collect();
    assert!(
        above
            .iter()
            .all(|t| ["Inspector", "Sketch fillet"].contains(t)),
        "{above:?}"
    );
    // Visible without scrolling: well inside the window.
    assert!(radius.y < 200., "{radius:?}");
    // Nothing unrelated to the tool is listed.
    for absent in [
        "Construction geometry",
        "Sketch elements",
        "Delete selected element",
        "Coincident",
    ] {
        assert!(find(&shown, absent).is_none(), "{absent}: {shown:?}");
    }
    // The value is editable through the field itself: click, type, Enter.
    let value = shown
        .iter()
        .find(|(t, p)| t == "2.0" && (p.y - radius.y).abs() < 6.)
        .map(|(_, p)| *p)
        .unwrap_or_else(|| panic!("radius value beside its label: {shown:?}"));
    h.press(value + Vec2::new(6., 6.));
    for event in [
        key(egui::Key::A, egui::Modifiers::COMMAND),
        egui::Event::Text("3.5".into()),
        key(egui::Key::Enter, egui::Modifiers::NONE),
    ] {
        h.frame(vec![event]).drop_without_applying_deltas();
    }
    h.frame(vec![]).drop_without_applying_deltas();
    assert_eq!(h.editor.corner_size_mm, 3.5);
    let shown = inspector(&mut h);
    assert!(find(&shown, "3.5").is_some(), "{shown:?}");
    assert!(
        find(&shown, " mm").is_some(),
        "the unit stays beside the value"
    );
}

#[test]
fn the_chamfer_distance_is_first_while_chamfer_is_active() {
    let mut h = with_rectangle();
    h.editor.set_tool(DrawTool::Chamfer).unwrap();
    let shown = inspector(&mut h);
    let distance = *find(&shown, "Chamfer distance").unwrap_or_else(|| panic!("{shown:?}"));
    assert!(
        shown
            .iter()
            .filter(|(_, p)| p.y < distance.y - 1.)
            .all(|(t, _)| ["Inspector", "Sketch chamfer"].contains(&t.as_str())),
        "{shown:?}"
    );
    assert!(find(&shown, "Fillet radius").is_none());
}

#[test]
fn drawing_tools_show_their_steps_and_options_and_no_element_values() {
    let mut h = with_rectangle();
    h.editor.set_tool(DrawTool::Line).unwrap();
    let shown = inspector(&mut h);
    let title = *find(&shown, "Line").unwrap_or_else(|| panic!("{shown:?}"));
    assert!(find(&shown, "Construction geometry").is_some(), "{shown:?}");
    for absent in ["Width", "Height", "Fillet radius", "Selected: Rectangle 1"] {
        assert!(find(&shown, absent).is_none(), "{absent}: {shown:?}");
    }
    assert!(
        shown
            .iter()
            .all(|(t, p)| t == "Inspector" || p.y >= title.y - 1.)
    );
}

#[test]
fn a_selected_element_shows_its_values_without_tool_settings_or_a_second_list() {
    let mut h = with_rectangle();
    h.click_label("Rectangle 1");
    let shown = inspector(&mut h);
    for present in ["Selected: Rectangle 1", "Width", "Height", "X", "Y"] {
        assert!(find(&shown, present).is_some(), "{present}: {shown:?}");
    }
    // Every size and position field carries its unit.
    assert!(
        shown.iter().filter(|(t, _)| t == "mm").count() >= 4,
        "{shown:?}"
    );
    for absent in [
        "Construction geometry",
        "Fillet radius",
        "Rectangle 1",
        "Line 1 (single curve)",
    ] {
        assert!(find(&shown, absent).is_none(), "{absent}: {shown:?}");
    }
    // The left list still names the shape exactly once.
    let all = texts(&mut h);
    assert_eq!(all.iter().filter(|(t, _)| t == "Rectangle 1").count(), 1);
    // Expressions are accepted as in other sketch dimensions.
    let revision = h.session.revision();
    h.enter(super::shape_list::ShapeField::Size(0), "=40 mm + 2 mm");
    assert_eq!(h.session.revision(), revision + 1, "{}", h.editor.message);
}
