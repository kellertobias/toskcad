//! CAD-416: sketch shape Inspector fields use the shared numeric layout and
//! still commit typed values.
use super::shape_list::{ShapeField, field_id};
use super::shape_list_tests::Host;
use super::*;
use multicad_core::constraint_sketch::{ConstraintKind, Dimension};

#[test]
fn shape_fields_are_full_width_with_label_above_unit_inside_and_commit_on_enter() {
    let mut h = Host::new();
    h.editor.set_tool(DrawTool::Rectangle).unwrap();
    h.frame(vec![]).drop_without_applying_deltas();
    h.click([0., 0.], Vec2::ZERO);
    h.click([100., 60.], Vec2::ZERO);
    h.click_label("Rectangle 1");
    let output = h.frame(vec![]);
    let texts: Vec<(String, Rect)> = output
        .shapes
        .iter()
        .filter_map(|shape| match &shape.shape {
            egui::Shape::Text(text) => Some((
                text.galley.text().to_owned(),
                text.galley.rect.translate(text.pos.to_vec2()),
            )),
            _ => None,
        })
        .collect();
    // Full-width separators mark the Inspector content width.
    let separators: Vec<(f32, f32)> = output
        .shapes
        .iter()
        .filter_map(|shape| match &shape.shape {
            egui::Shape::LineSegment { points, .. } if points[0].y == points[1].y => {
                Some((points[0].x.min(points[1].x), points[0].x.max(points[1].x)))
            }
            _ => None,
        })
        .collect();
    output.drop_without_applying_deltas();
    let fields = [
        (ShapeField::Size(0), "Width"),
        (ShapeField::Size(1), "Height"),
        (ShapeField::X, "X"),
        (ShapeField::Y, "Y"),
    ];
    let mut previous_top = f32::NEG_INFINITY;
    for (field, label) in fields {
        let input = h
            .ctx
            .read_response(field_id(field))
            .unwrap_or_else(|| panic!("{label} field"))
            .rect;
        assert!(
            separators
                .iter()
                .any(|(l, r)| (l - input.left()).abs() < 0.5 && (r - input.right()).abs() < 0.5),
            "{label} input {input:?} spans the Inspector content width: {separators:?}"
        );
        let label_rect = texts
            .iter()
            .filter(|(t, r)| t == label && r.bottom() <= input.top() + 0.5)
            .map(|(_, r)| *r)
            .find(|r| input.top() - r.bottom() < 12.)
            .unwrap_or_else(|| panic!("{label} label directly above its input"));
        assert!(label_rect.left() >= input.left() - 0.5);
        let unit = texts
            .iter()
            .find(|(t, r)| t == "mm" && input.expand(0.5).contains_rect(*r))
            .map(|(_, r)| *r)
            .unwrap_or_else(|| panic!("{label} shows mm inside its input"));
        assert!(input.right() - unit.right() < 12.);
        assert!(input.top() > previous_top, "fields stack in one column");
        previous_top = input.top();
    }
    // The edited value is the number alone; Enter commits it as a dimension.
    let revision = h.session.revision();
    h.enter(ShapeField::Size(0), "120");
    assert_eq!(h.session.revision(), revision + 1, "{}", h.editor.message);
    let input = h
        .ctx
        .read_response(field_id(ShapeField::Size(0)))
        .unwrap()
        .rect;
    let output = h.frame(vec![]);
    let shown: Vec<String> = output
        .shapes
        .iter()
        .filter_map(|shape| match &shape.shape {
            egui::Shape::Text(text) if input.contains(text.pos + text.galley.size() * 0.5) => {
                Some(text.galley.text().to_owned())
            }
            _ => None,
        })
        .collect();
    output.drop_without_applying_deltas();
    assert!(
        shown.contains(&"120".to_string()) && shown.contains(&"mm".to_string()),
        "value and unit are separate texts in the field: {shown:?}"
    );
    let source = shape_list::canonical(h.session.document(), h.editor.selected.unwrap()).unwrap();
    assert!(source.sketch.constraints.iter().any(|c| matches!(
        &c.kind,
        ConstraintKind::Length { dimension: Dimension::Literal { value, unit }, .. }
            if *value == 120. && unit == "mm"
    )));
}
