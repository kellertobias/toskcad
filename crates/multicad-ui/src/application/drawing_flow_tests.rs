//! CAD-436: the first correct front view takes New drawing, Add view and one
//! click on the sheet, through the full application and its ribbon.
use super::ui_harness::Harness;
use super::*;
use crate::drawings::guided_tests::{box_factory, boxed_document};

#[test]
fn first_front_view_takes_new_drawing_add_view_and_one_sheet_click() {
    // A 120 × 40 × 60 mm block: its front view is 120 wide and 60 high.
    let (document, parts, boxes) = boxed_document(&[[120., 40., 60.]]);
    let mut h = Harness::new(document);
    h.app.set_drawing_backends(box_factory(boxes));
    h.app
        .invoke(Command::Mode(Mode::Drawings), &[], Channel::Ui)
        .unwrap();
    h.frame(vec![]);
    let start = h.revision();
    // 1. New drawing, the one action on the empty sheet area.
    h.click_label("New drawing");
    assert_eq!(h.revision(), start + 1);
    // 2. Add view in the ribbon.
    h.click_label("Add view");
    h.until("the fitted front view", |app| {
        app.current()
            .drawings
            .next_step_hint()
            .is_some_and(|hint| hint.starts_with("Click the sheet to place the front view"))
    });
    assert_eq!(h.revision(), start + 1, "Nothing applied before the click");
    // 3. One click on the sheet.
    let at = h.app.current().drawings.sheet_point([110., 85.]).unwrap();
    h.click_at(at);
    h.until("the applied front view", |app| {
        !app.current().drawings.view_placement_active()
    });
    assert_eq!(
        h.revision(),
        start + 2,
        "One edit for the page, one for the view"
    );
    let document = h.app.current().session.document();
    let book = multicad_core::drawings::books(document).next().unwrap();
    let page = &multicad_core::drawings::read(document, book).unwrap().pages[0];
    assert_eq!(page.views.len(), 1);
    let view = &page.views[0];
    assert_eq!((view.frame.u, view.frame.v), ([1., 0., 0.], [0., 0., 1.]));
    assert_eq!(
        view.source,
        multicad_core::drawings::ViewSource::Model { entities: parts }
    );
    // Centered on the block: x 60, z 30 are the front view's center.
    assert_eq!(view.frame.origin_mm, [60., 0., 30.]);
    // The largest standard scale that fits A4 beside the title block.
    assert_eq!(multicad_core::drawings::scale_text(view.scale), "2:1");
    // Placed where clicked, nudged so the 240 × 120 mm view stays printable.
    let layout = multicad_core::drawings::layout_page(page, 1, 1).unwrap();
    let [x0, y0, x1, y1] = layout.printable_mm;
    let [x, y] = view.position_mm;
    assert!(x - 120. >= x0 - 1e-6 && x + 120. <= x1 + 1e-6, "{x}");
    assert!(y - 60. >= y0 - 1e-6 && y + 60. <= y1 + 1e-6, "{y}");
    assert!(
        (y - 85.).abs() < 0.5,
        "Kept the clicked height where it fits"
    );
    assert!(!h.app.current().drawings.dirty());
}
