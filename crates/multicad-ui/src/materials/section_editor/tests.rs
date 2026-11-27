//! CAD-385: the section modal draws, dimensions, validates and undoes with
//! the sketch command pipeline, and reports Finish/Cancel.
use super::*;

const SNAP: f64 = 1.;

fn click(editor: &mut SectionEditor, points: &[[f64; 2]]) {
    for p in points {
        editor.click(*p, SNAP);
    }
}
fn tube() -> SectionEditor {
    let mut editor = SectionEditor::new(Sketch::default());
    assert_eq!(
        editor.tool(),
        Tool::Rectangle,
        "an empty sketch starts drawing"
    );
    click(&mut editor, &[[0., 0.], [40., 30.]]);
    editor.set_tool(Tool::Circle);
    click(&mut editor, &[[20., 15.], [28., 15.]]);
    editor
}

#[test]
fn drawn_rectangle_and_circle_form_a_hollow_section_with_undo_and_redo() {
    let mut editor = tube();
    assert!(editor.touched());
    let section = editor.status().unwrap();
    assert_eq!(section.size_mm, [40., 30.]);
    assert_eq!(section.holes.len(), 1);
    assert!((section.area_mm2 - (1200. - std::f64::consts::PI * 64.)).abs() < 1e-9);
    editor.undo();
    assert!(editor.status().unwrap().holes.is_empty());
    editor.redo();
    assert_eq!(editor.status().unwrap().holes.len(), 1);
    // Selecting the circle and deleting it removes the hole and its center.
    editor.set_tool(Tool::Select);
    click(&mut editor, &[[28., 15.]]);
    editor.delete_selected();
    let section = editor.status().unwrap();
    assert!(section.holes.is_empty());
    assert_eq!(editor.sketch().points.len(), 4);
}

#[test]
fn a_line_chain_closes_on_its_first_point_and_an_open_chain_is_explained() {
    let mut editor = SectionEditor::new(Sketch::default());
    editor.set_tool(Tool::Line);
    click(
        &mut editor,
        &[[0., 0.], [30., 0.], [30., 20.], [10., 20.], [10., 10.]],
    );
    let open = editor.status().unwrap_err();
    assert!(open.contains("open"), "{open}");
    // Clicking within the snap distance of the first point closes the loop.
    click(&mut editor, &[[0., 10.], [0.4, 0.3]]);
    let section = editor.status().unwrap();
    assert_eq!(section.size_mm, [30., 20.]);
    assert!((section.area_mm2 - (30. * 10. + 20. * 10.)).abs() < 1e-9);
    // An arc rounds a corner of another outline drawn separately.
    let mut arc = SectionEditor::new(Sketch::default());
    arc.set_tool(Tool::Line);
    click(&mut arc, &[[0., 0.], [20., 0.]]);
    arc.set_tool(Tool::Arc);
    click(&mut arc, &[[20., 10.], [20., 0.], [20., 20.]]);
    arc.set_tool(Tool::Line);
    click(&mut arc, &[[20., 20.], [0., 20.], [0., 0.]]);
    let section = arc.status().unwrap();
    assert_eq!(section.size_mm, [30., 20.]);
}

#[test]
fn dimensions_resize_by_number_and_refuse_expressions() {
    let mut editor = tube();
    editor.set_tool(Tool::Select);
    click(&mut editor, &[[20., 0.]]);
    editor.add_dimension();
    let (id, field) = editor.dimensions[0].clone();
    assert_eq!(field.text, "40");
    editor.set_dimension(id, "60");
    let section = editor.status().unwrap();
    assert!((section.size_mm[0] - 60.).abs() < 1e-6);
    assert!((section.size_mm[1] - 30.).abs() < 1e-6);
    let before = editor.sketch().clone();
    editor.set_dimension(id, "=$width");
    assert!(editor.message().unwrap().contains("must be numbers"));
    assert_eq!(
        editor.sketch(),
        &before,
        "a refused dimension changes nothing"
    );
    editor.set_dimension(id, "-5");
    assert_eq!(editor.sketch(), &before);
    editor.add_dimension();
    assert!(editor.message().unwrap().contains("already dimensioned"));
    editor.remove_dimension(id);
    assert!(editor.dimensions.is_empty());
}

#[test]
fn presets_seed_a_dimensioned_editable_sketch_without_touching_it() {
    let spec = |shape, w, h, wall| StockKind::Profile {
        designation: "x".into(),
        length_mm: 6000.,
        width_mm: w,
        height_mm: h,
        wall_mm: wall,
        shape,
        section: None,
        corner_radius_mm: None,
        web_mm: None,
    };
    let tube = SectionEditor::seeded(&spec(ProfileShape::Rectangular, 30., 20., Some(2.)));
    assert!(!tube.touched());
    let section = tube.status().unwrap();
    assert_eq!(section.size_mm, [30., 20.]);
    assert!((section.area_mm2 - (600. - 26. * 16.)).abs() < 1e-9);
    assert_eq!(tube.dimensions.len(), 2);
    let round = SectionEditor::seeded(&spec(ProfileShape::Round, 30., 30., Some(3.)));
    assert_eq!(round.status().unwrap().holes.len(), 1);
    let angle = SectionEditor::seeded(&spec(ProfileShape::Angle, 40., 30., Some(4.)));
    let section = angle.status().unwrap();
    assert!((section.area_mm2 - (40. * 4. + 26. * 4.)).abs() < 1e-9);
    let empty = SectionEditor::seeded(&spec(ProfileShape::ItemLine8_40x40, 40., 40., None));
    assert!(empty.sketch().curves.is_empty());
}

fn frame(
    ctx: &egui::Context,
    editor: &mut SectionEditor,
    events: Vec<egui::Event>,
) -> (Option<Outcome>, Vec<(String, egui::Rect)>) {
    let mut outcome = None;
    let output = ctx.run_ui(
        egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(
                egui::Pos2::ZERO,
                egui::vec2(1000., 800.),
            )),
            focused: true,
            events,
            ..Default::default()
        },
        |ui| outcome = editor.show(ui.ctx()),
    );
    let texts = output
        .shapes
        .iter()
        .filter_map(|s| match &s.shape {
            egui::Shape::Text(t) => Some((
                t.galley.text().to_owned(),
                egui::Rect::from_min_size(t.pos, t.galley.size()),
            )),
            _ => None,
        })
        .collect();
    output.drop_without_applying_deltas();
    (outcome, texts)
}
fn press(at: egui::Pos2) -> Vec<egui::Event> {
    [true, false]
        .into_iter()
        .flat_map(|pressed| {
            [
                egui::Event::PointerMoved(at),
                egui::Event::PointerButton {
                    pos: at,
                    button: egui::PointerButton::Primary,
                    pressed,
                    modifiers: Default::default(),
                },
            ]
        })
        .collect()
}
fn text(texts: &[(String, egui::Rect)], prefix: &str) -> Option<egui::Rect> {
    texts
        .iter()
        .find(|(t, _)| t.starts_with(prefix))
        .map(|(_, r)| *r)
}

#[test]
fn painted_modal_draws_on_the_canvas_reports_status_and_finishes_or_cancels() {
    let ctx = egui::Context::default();
    let mut editor = SectionEditor::new(Sketch::default());
    frame(&ctx, &mut editor, vec![]);
    let (_, texts) = frame(&ctx, &mut editor, vec![]);
    assert!(text(&texts, "Profile section sketch").is_some());
    assert!(
        text(&texts, "Draw a closed section outline").is_some(),
        "{texts:?}"
    );
    // Two canvas clicks with the Rectangle tool draw a rectangle.
    let finish = text(&texts, "Finish").unwrap();
    let below_tools = text(&texts, Tool::Rectangle.hint()).unwrap();
    let a = egui::pos2(below_tools.left() + 80., below_tools.bottom() + 80.);
    let b = egui::pos2(a.x + 120., a.y + 90.);
    frame(&ctx, &mut editor, press(a));
    frame(&ctx, &mut editor, press(b));
    let (_, texts) = frame(&ctx, &mut editor, vec![]);
    assert_eq!(editor.sketch().curves.len(), 4, "{:?}", editor.message());
    assert!(text(&texts, "Valid section").is_some(), "{texts:?}");
    let (outcome, _) = frame(&ctx, &mut editor, press(finish.center()));
    let Some(Outcome::Finished(solved)) = outcome else {
        panic!("finished: {outcome:?}")
    };
    assert_eq!(&solved.sketch, editor.sketch());
    assert_eq!(solved.section, editor.status().unwrap());
    // Escape cancels when no tool input is pending.
    let (outcome, _) = frame(
        &ctx,
        &mut editor,
        vec![egui::Event::Key {
            key: egui::Key::Escape,
            physical_key: None,
            pressed: true,
            repeat: false,
            modifiers: Default::default(),
        }],
    );
    assert_eq!(outcome, Some(Outcome::Cancelled));
}
