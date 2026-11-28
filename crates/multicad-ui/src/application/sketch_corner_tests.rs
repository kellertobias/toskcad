//! CAD-414: sketch fillets and chamfers through the full application. A size
//! typed in the active corner tool resizes that tool session's treatments
//! through their canonical dimension in one Undo step; an existing fillet or
//! chamfer shows and edits its size when selected; both are typed shapes
//! named with their size and linked to the shape they modify.
use super::sketch_select_tests::SketchHarness;
use super::*;
use crate::sketch_editor::DrawTool;
use multicad_core::constraint_sketch::{
    CurveKind, Dimension, Shape, ShapeKind,
    shape_edit::{self},
};
use multicad_precision::PortableMath as _;

const ZERO: egui::Vec2 = egui::Vec2::ZERO;

impl SketchHarness {
    fn corners(&self, kind: ShapeKind) -> Vec<Shape> {
        self.source()
            .sketch
            .shapes
            .into_iter()
            .filter(|s| s.kind == kind)
            .collect()
    }
    /// Canonical size of a fillet or chamfer, in mm.
    fn size_of(&self, shape: &Shape) -> f64 {
        let source = self.source();
        match shape_edit::corner_size(&source.sketch, shape) {
            Some(Dimension::Literal { value, .. }) => value,
            other => panic!("{other:?}"),
        }
    }
    /// Solved radius of a fillet arc.
    fn solved_radius(&self, shape: &Shape) -> f64 {
        let project = self.h.app.current();
        let sketch = project.sketch.selected.unwrap();
        let solved = project.sketch.solved_for_tests(sketch).unwrap();
        let arc = solved
            .curves
            .iter()
            .find(|c| c.id == shape.curves[0])
            .unwrap();
        let CurveKind::Arc { center, start, .. } = arc.kind else {
            panic!("fillet arc")
        };
        let at = |id| {
            solved
                .points
                .iter()
                .find(|p| p.id == id)
                .unwrap()
                .position_mm
        };
        let (c, a) = (at(center), at(start));
        (a[0] - c[0]).portable_hypot(a[1] - c[1])
    }
    fn undo(&mut self) {
        self.h.app.invoke(Command::Undo, &[], Channel::Ui).unwrap();
        self.settle();
    }
}

/// AC1 and AC4: the reported sequence. Fillet a corner, type 10 in the
/// radius field, press Enter: the existing fillet becomes R10 through its
/// own Radius dimension, as one Undo step, and the solved arc follows.
#[test]
fn typing_a_new_radius_in_the_active_fillet_tool_resizes_the_fillet_just_made() {
    let mut s = SketchHarness::new();
    s.draw(DrawTool::Rectangle, &["0", "0", "60", "40"]);
    s.tool(DrawTool::Fillet);
    s.click([60., 40.], ZERO);
    let fillet = s.corners(ShapeKind::Fillet).remove(0);
    assert!((s.size_of(&fillet) - 2.).abs() < 1e-9, "default 2 mm");
    let constraint = shape_edit::corner_constraints(&s.source().sketch, &fillet);
    let constraints = s.source().sketch.constraints.len();
    let revision = s.h.revision();
    s.type_beside("Fillet radius", "10");
    assert_eq!(
        s.h.revision(),
        revision + 1,
        "one transaction: {}",
        s.message()
    );
    assert_eq!(s.h.app.current().sketch.tool, Some(DrawTool::Fillet));
    let source = s.source();
    assert_eq!(
        source.sketch.constraints.len(),
        constraints,
        "edited in place"
    );
    assert_eq!(
        shape_edit::corner_constraints(&source.sketch, &fillet),
        constraint,
        "the same canonical dimension"
    );
    assert!((s.size_of(&fillet) - 10.).abs() < 1e-9);
    assert!((s.solved_radius(&fillet) - 10.).abs() < 1e-6, "solved");
    assert!(s.message().contains("one Undo step"), "{}", s.message());
    // The next corner uses the new size too.
    s.click([0., 40.], ZERO);
    let fillets = s.corners(ShapeKind::Fillet);
    assert_eq!(fillets.len(), 2);
    assert!((s.size_of(&fillets[1]) - 10.).abs() < 1e-9);
    // Undo removes the second fillet; the next Undo reverts the resize alone.
    s.undo();
    s.undo();
    let fillet = s.corners(ShapeKind::Fillet).remove(0);
    assert!(
        (s.size_of(&fillet) - 2.).abs() < 1e-9,
        "one Undo per resize"
    );
}

/// AC1: every fillet of the session follows; an earlier session's fillet
/// does not. Chamfer distances change on both sides.
#[test]
fn a_new_size_resizes_this_sessions_treatments_only() {
    let mut s = SketchHarness::new();
    s.draw(DrawTool::Rectangle, &["0", "0", "60", "40"]);
    s.tool(DrawTool::Fillet);
    s.click([60., 40.], ZERO);
    s.h.key(egui::Key::Enter);
    s.settle();
    s.tool(DrawTool::Fillet);
    s.click([0., 40.], ZERO);
    s.click([60., 0.], ZERO);
    s.type_beside("Fillet radius", "5");
    let fillets = s.corners(ShapeKind::Fillet);
    assert_eq!(fillets.len(), 3);
    let sizes: Vec<_> = fillets.iter().map(|f| s.size_of(f)).collect();
    assert_eq!(sizes, vec![2., 5., 5.], "only this session's fillets");
    s.h.key(egui::Key::Enter);
    s.settle();
    s.tool(DrawTool::Chamfer);
    s.click([0., 0.], ZERO);
    let chamfer = s.corners(ShapeKind::Chamfer).remove(0);
    let revision = s.h.revision();
    s.type_beside("Chamfer distance", "4");
    assert_eq!(s.h.revision(), revision + 1, "{}", s.message());
    let source = s.source();
    let setbacks = shape_edit::corner_constraints(&source.sketch, &chamfer);
    assert_eq!(setbacks.len(), 2);
    for id in setbacks {
        let c = source
            .sketch
            .constraints
            .iter()
            .find(|c| c.id == id)
            .unwrap();
        assert!(matches!(
            c.kind.dimension(),
            Some(Dimension::Literal { value, .. }) if (*value - 4.).abs() < 1e-9
        ));
    }
}

/// AC2 and AC3: once the tool has ended, a selected fillet or chamfer shows
/// its size and edits it; both are typed shapes linked to the rectangle and
/// listed by kind and size, not as arcs.
#[test]
fn selecting_an_existing_fillet_or_chamfer_shows_and_edits_its_size() {
    let mut s = SketchHarness::new();
    s.draw(DrawTool::Rectangle, &["0", "0", "60", "40"]);
    s.tool(DrawTool::Fillet);
    s.type_beside("Fillet radius", "10");
    s.click([60., 40.], ZERO);
    s.h.key(egui::Key::Enter);
    s.tool(DrawTool::Chamfer);
    s.type_beside("Chamfer distance", "3");
    s.click([0., 0.], ZERO);
    s.h.key(egui::Key::Enter);
    s.settle();
    let rectangle = s.source().sketch.shapes[0].clone();
    let fillet = s.corners(ShapeKind::Fillet).remove(0);
    let chamfer = s.corners(ShapeKind::Chamfer).remove(0);
    assert_eq!(fillet.modifies, Some(rectangle.id));
    assert_eq!(chamfer.modifies, Some(rectangle.id));
    let texts = s.h.texts();
    assert!(texts.iter().any(|t| t == "Fillet R10"), "{texts:?}");
    assert!(texts.iter().any(|t| t == "Chamfer 3"), "{texts:?}");
    assert!(
        !texts.iter().any(|t| t.starts_with("Arc ")),
        "no fillet listed as an arc: {texts:?}"
    );
    // The fillet in the list: its radius is shown and edited in place.
    s.h.click_label("Fillet R10");
    s.settle();
    let texts = s.h.texts();
    assert!(
        texts.iter().any(|t| t == "Selected: Fillet R10"),
        "{texts:?}"
    );
    assert!(texts.iter().any(|t| t == "Radius"), "{texts:?}");
    assert!(texts.iter().any(|t| t == "10"), "{texts:?}");
    let revision = s.h.revision();
    s.type_beside("Radius", "6");
    assert_eq!(s.h.revision(), revision + 1, "{}", s.message());
    assert!((s.size_of(&fillet) - 6.).abs() < 1e-9);
    assert!((s.solved_radius(&fillet) - 6.).abs() < 1e-6);
    assert!(s.h.texts().iter().any(|t| t == "Fillet R6"));
    // The chamfer: both setbacks follow its Distance.
    s.h.click_label("Chamfer 3");
    s.settle();
    let texts = s.h.texts();
    assert!(texts.iter().any(|t| t == "Distance"), "{texts:?}");
    assert!(texts.iter().any(|t| t == "3"), "{texts:?}");
    s.type_beside("Distance", "2.5");
    assert!((s.size_of(&chamfer) - 2.5).abs() < 1e-9);
    let source = s.source();
    for id in shape_edit::corner_constraints(&source.sketch, &chamfer) {
        let c = source
            .sketch
            .constraints
            .iter()
            .find(|c| c.id == id)
            .unwrap();
        assert!(matches!(
            c.kind.dimension(),
            Some(Dimension::Literal { value, .. }) if (*value - 2.5).abs() < 1e-9
        ));
    }
    assert!(s.h.texts().iter().any(|t| t == "Chamfer 2.5"));
}
