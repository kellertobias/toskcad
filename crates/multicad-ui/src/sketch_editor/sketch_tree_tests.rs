//! Sketch-mode sidebar (CAD-415), driven through the real editor host: the
//! Model / Sketch tabs, the element tree's icons and nesting, and selection
//! sync between tree, canvas and Inspector.
use super::shape_list::{ElementKey, ListedElement};
use super::shape_list_tests::Host;
use super::sketch_tree::{ElementType, SidebarTab, TreeEntry, element_rows, tab_id};
use super::*;
use multicad_core::constraint_sketch::{
    Curve, CurveKind, Point, SplineBasis,
    text::{TextOperation, TextStyle},
};
use multicad_core::sketch_draw::{DrawInput, DrawRequest, PointInput};
use std::collections::BTreeSet;

fn click_tab(h: &mut Host, tab: SidebarTab) {
    h.frame(vec![]).drop_without_applying_deltas();
    let rect = h
        .ctx
        .read_response(tab_id(tab))
        .unwrap_or_else(|| panic!("Missing {tab:?} tab"))
        .rect;
    h.press(rect.center());
    h.frame(vec![]).drop_without_applying_deltas();
}
fn tab_shown(h: &mut Host, tab: SidebarTab) -> bool {
    h.frame(vec![]).drop_without_applying_deltas();
    h.frame(vec![]).drop_without_applying_deltas();
    h.ctx.read_response(tab_id(tab)).is_some()
}
fn draw(h: &mut Host, input: DrawInput) {
    h.editor
        .draw(
            &mut h.session,
            &DrawRequest {
                input,
                construction: false,
            },
            Channel::Ui,
        )
        .unwrap_or_else(|e| panic!("{e}"));
    h.settle();
}
/// Commit an edited canonical source as one transaction.
fn commit(h: &mut Host, source: Source) {
    let edit = multicad_solver::prepare_canonical_edit(
        &h.session,
        h.session.revision(),
        Actor::local(Channel::Ui),
        h.sketch(),
        source,
        None,
    )
    .unwrap();
    h.session.commit(edit.preview.unwrap()).unwrap();
    h.settle();
}
fn element_rows_painted(h: &Host) -> Vec<(ElementType, String, usize)> {
    h.editor
        .sketch_sidebar
        .rows
        .iter()
        .filter_map(|row| match row {
            TreeEntry::Element {
                kind, label, depth, ..
            } => Some((*kind, label.clone(), *depth)),
            _ => None,
        })
        .collect()
}
/// Whether `label` is painted over the tree's selection fill.
fn painted_selected(h: &mut Host, label: &str) -> bool {
    let output = h.frame(vec![]);
    let fill = h.ctx.global_style().visuals.selection.bg_fill;
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
    let selected = output.shapes.iter().any(|shape| {
        matches!(&shape.shape, egui::Shape::Rect(r) if r.fill == fill && r.rect.contains(at))
    });
    output.drop_without_applying_deltas();
    selected
}

#[test]
fn sketch_tab_opens_with_the_sketch_switches_to_the_model_tree_and_leaves_with_it() {
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
    // Entering the sketch shows both tabs with Sketch selected.
    assert!(tab_shown(&mut h, SidebarTab::Model));
    assert!(tab_shown(&mut h, SidebarTab::Sketch));
    assert_eq!(h.editor.sketch_sidebar.tab, SidebarTab::Sketch);
    let labels = h.labels();
    assert!(labels.iter().any(|l| l == "Sketch elements"), "{labels:?}");
    assert!(labels.iter().any(|l| l == "Rectangle 1"), "{labels:?}");
    assert!(
        !labels.iter().any(|l| l == "Show sources"),
        "the model tree is hidden behind its tab: {labels:?}"
    );
    // Model tab: the normal model tree, no sketch elements.
    click_tab(&mut h, SidebarTab::Model);
    assert_eq!(h.editor.sketch_sidebar.tab, SidebarTab::Model);
    let labels = h.labels();
    assert!(labels.iter().any(|l| l == "Show sources"), "{labels:?}");
    assert!(!labels.iter().any(|l| l == "Sketch elements"), "{labels:?}");
    assert!(!labels.iter().any(|l| l == "Rectangle 1"), "{labels:?}");
    // And back.
    click_tab(&mut h, SidebarTab::Sketch);
    assert!(h.labels().iter().any(|l| l == "Rectangle 1"));
    // Leaving the sketch: no tabs, only the model tree.
    click_tab(&mut h, SidebarTab::Model);
    h.editor.finish().unwrap();
    assert!(!tab_shown(&mut h, SidebarTab::Model));
    assert!(!tab_shown(&mut h, SidebarTab::Sketch));
    let labels = h.labels();
    assert!(labels.iter().any(|l| l == "Show sources"), "{labels:?}");
    assert!(!labels.iter().any(|l| l == "Sketch elements"), "{labels:?}");
    // Reopening selects the Sketch tab again, whatever was chosen before.
    let sketch = h.sketch();
    h.editor
        .open_at(h.session.document(), sketch, h.session.revision())
        .unwrap();
    assert!(tab_shown(&mut h, SidebarTab::Sketch));
    assert_eq!(h.editor.sketch_sidebar.tab, SidebarTab::Sketch);
    assert!(h.labels().iter().any(|l| l == "Rectangle 1"));
}

/// Rectangle, circle, arc, line, polyline, spline, an unrecorded curve, text
/// and a standalone point, each drawn the way the editor draws it.
fn mixed_sketch() -> Host {
    let mut h = Host::new();
    let p = PointInput::Position;
    draw(
        &mut h,
        DrawInput::Rectangle {
            a: p([0., 0.]),
            b: p([40., 20.]),
        },
    );
    draw(
        &mut h,
        DrawInput::Circle {
            center: p([80., 10.]),
            radius_mm: 8.,
        },
    );
    draw(
        &mut h,
        DrawInput::Arc {
            center: p([120., 10.]),
            start: p([130., 10.]),
            end: p([110., 10.]),
            clockwise: false,
        },
    );
    draw(
        &mut h,
        DrawInput::Line {
            start: p([0., 50.]),
            end: p([40., 60.]),
        },
    );
    draw(
        &mut h,
        DrawInput::Polyline {
            vertices: vec![p([60., 50.]), p([70., 60.]), p([80., 50.])],
            closed: false,
        },
    );
    draw(
        &mut h,
        DrawInput::BSpline {
            poles: vec![
                p([100., 50.]),
                p([105., 60.]),
                p([115., 45.]),
                p([125., 55.]),
            ],
            basis: SplineBasis::open_uniform(3, 4).unwrap(),
            closed: false,
        },
    );
    // A curve with no shape record (older sources, operation results).
    let mut source = h.source();
    let (a, b) = (EntityId::new(), EntityId::new());
    source.sketch.points.extend([
        Point {
            id: a,
            position_mm: [0., -40.],
        },
        Point {
            id: b,
            position_mm: [30., -40.],
        },
    ]);
    source.sketch.curves.push(Curve {
        id: EntityId::new(),
        construction: false,
        kind: CurveKind::Line { start: a, end: b },
    });
    // A standalone point.
    source.sketch.points.push(Point {
        id: EntityId::new(),
        position_mm: [60., -40.],
    });
    commit(&mut h, source);
    let (source, _) = TextOperation::Add {
        style: TextStyle::new("Hi", 6., [90., -40.]),
    }
    .apply(&h.source())
    .unwrap();
    commit(&mut h, source);
    h
}

#[test]
fn a_mixed_sketch_lists_every_element_with_its_own_icon_and_label() {
    let mut h = mixed_sketch();
    let labels = h.labels();
    let rows = element_rows_painted(&h);
    let expected = [
        (ElementType::Rectangle, "Rectangle 1"),
        (ElementType::Circle, "Circle 1"),
        (ElementType::Arc, "Arc 1"),
        (ElementType::Line, "Line 1"),
        (ElementType::Polyline, "Polyline 1"),
        (ElementType::Spline, "Spline 1"),
        (ElementType::Line, "Line 1 (single curve)"),
        (ElementType::Text, "Text 1 \u{201c}Hi\u{201d}"),
        (ElementType::Point, "Point 1"),
    ];
    assert_eq!(
        rows,
        expected
            .iter()
            .map(|(kind, label)| (*kind, label.to_string(), 0))
            .collect::<Vec<_>>()
    );
    for (_, label) in expected {
        assert!(labels.iter().any(|l| l == label), "{label}: {labels:?}");
    }
}

#[test]
fn every_element_type_has_a_distinct_icon_and_label() {
    let icons: BTreeSet<_> = ElementType::ALL
        .iter()
        .map(|kind| format!("{:?}", std::mem::discriminant(&kind.icon())))
        .collect();
    assert_eq!(icons.len(), ElementType::ALL.len());
    let labels: BTreeSet<_> = ElementType::ALL.iter().map(|k| k.label()).collect();
    assert_eq!(labels.len(), ElementType::ALL.len());
    // Each icon paints strokes inside its own tree-sized button.
    let ctx = egui::Context::default();
    for kind in ElementType::ALL {
        let mut rect = None;
        let output = ctx.run_ui(egui::RawInput::default(), |ui| {
            rect = Some(crate::icons::model_tree_icon(ui, kind.icon(), kind.label(), false).rect);
        });
        let rect = rect.unwrap().expand(1.);
        let painted = output
            .shapes
            .iter()
            .filter(|s| !matches!(s.shape, egui::Shape::Noop))
            .count();
        assert!(painted > 0, "{kind:?} paints nothing");
        // The icons added for the tree stay inside their button; older
        // tool icons (Spline's control poles) are unchanged.
        let added = matches!(
            kind,
            ElementType::Slot
                | ElementType::Text
                | ElementType::Point
                | ElementType::Projected
                | ElementType::Fillet
                | ElementType::Chamfer
        );
        for shape in output.shapes.iter().filter(|_| added) {
            let bounds = shape.shape.visual_bounding_rect();
            assert!(
                !bounds.is_positive() || rect.contains_rect(bounds),
                "{kind:?} paints outside its icon: {bounds:?} in {rect:?}"
            );
        }
        output.drop_without_applying_deltas();
    }
}

fn listed(
    key: ElementKey,
    kind: ElementType,
    label: &str,
    parent: Option<ElementKey>,
) -> ListedElement {
    ListedElement {
        key,
        label: label.into(),
        curves: vec![],
        kind,
        parent,
    }
}

#[test]
fn fillets_and_chamfers_are_children_of_the_shape_they_modify() {
    let (rect, fillet, chamfer, circle, orphan) = (
        ElementKey::Shape(EntityId::new()),
        ElementKey::Shape(EntityId::new()),
        ElementKey::Shape(EntityId::new()),
        ElementKey::Shape(EntityId::new()),
        ElementKey::Shape(EntityId::new()),
    );
    let elements = vec![
        listed(rect, ElementType::Rectangle, "Rectangle 1", None),
        listed(circle, ElementType::Circle, "Circle 1", None),
        listed(fillet, ElementType::Fillet, "Fillet 1", Some(rect)),
        listed(chamfer, ElementType::Chamfer, "Chamfer 1", Some(rect)),
        // Its shape was deleted: it stays visible at the top level.
        listed(
            orphan,
            ElementType::Fillet,
            "Fillet 2",
            Some(ElementKey::Shape(EntityId::new())),
        ),
    ];
    let shown = |collapsed: &BTreeSet<ElementKey>| -> Vec<(String, usize, bool)> {
        element_rows(&elements, collapsed)
            .into_iter()
            .filter_map(|row| match row {
                TreeEntry::Element {
                    label,
                    depth,
                    children,
                    ..
                } => Some((label, depth, children)),
                _ => None,
            })
            .collect()
    };
    let row = |l: &str, d, c| (l.to_string(), d, c);
    assert_eq!(
        shown(&BTreeSet::new()),
        vec![
            row("Rectangle 1", 0, true),
            row("Fillet 1", 1, false),
            row("Chamfer 1", 1, false),
            row("Circle 1", 0, false),
            row("Fillet 2", 0, false),
        ]
    );
    assert_eq!(
        shown(&BTreeSet::from([rect])),
        vec![
            row("Rectangle 1", 0, true),
            row("Circle 1", 0, false),
            row("Fillet 2", 0, false),
        ],
        "collapsing the rectangle hides its corners"
    );
}

#[test]
fn constraints_and_dimensions_are_one_collapsible_group() {
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
    let count = h.source().sketch.constraints.len();
    assert!(
        count > 0,
        "a rectangle carries its horizontal/vertical constraints"
    );
    let header = format!("Constraints and dimensions ({count})");
    let labels = h.labels();
    assert!(labels.contains(&header), "{labels:?}");
    assert!(
        !labels.iter().any(|l| l == "Horizontal 1"),
        "collapsed by default: {labels:?}"
    );
    h.click_label(&header);
    let labels = h.labels();
    assert!(labels.iter().any(|l| l == "Horizontal 1"), "{labels:?}");
    assert!(
        h.editor
            .sketch_sidebar
            .rows
            .iter()
            .filter(|r| matches!(r, TreeEntry::Constraint { .. }))
            .count()
            == count
    );
    h.click_label(&header);
    assert!(!h.labels().iter().any(|l| l == "Horizontal 1"));
}

#[test]
fn tree_canvas_and_inspector_selection_stay_in_sync() {
    let mut h = mixed_sketch();
    let source = h.source();
    let (rect, circle) = (
        source.sketch.shapes[0].clone(),
        source.sketch.shapes[1].clone(),
    );
    // Tree → canvas and Inspector.
    h.click_label("Circle 1");
    assert_eq!(h.editor.shape, Some(ElementKey::Shape(circle.id)));
    assert_eq!(
        h.editor.element_highlight(h.session.document()),
        circle.curves.iter().copied().collect()
    );
    let labels = h.labels();
    assert!(
        labels.iter().any(|l| l == "Selected: Circle 1"),
        "{labels:?}"
    );
    assert!(labels.iter().any(|l| l == "Diameter"), "{labels:?}");
    assert!(painted_selected(&mut h, "Circle 1"));
    assert!(!painted_selected(&mut h, "Rectangle 1"));
    // Canvas → tree: the rectangle's bottom edge selects and highlights its
    // row, and the Inspector follows.
    h.click([20., 0.], Vec2::new(0., 1.));
    assert_eq!(
        h.editor.shape,
        Some(ElementKey::Shape(rect.id)),
        "{}",
        h.editor.message
    );
    assert!(painted_selected(&mut h, "Rectangle 1"));
    assert!(!painted_selected(&mut h, "Circle 1"));
    assert!(h.labels().iter().any(|l| l == "Selected: Rectangle 1"));
    // A standalone point: canvas pick selects its row; tree pick selects the
    // point for constraints.
    h.click([60., -40.], Vec2::ZERO);
    let point = match h.editor.shape {
        Some(ElementKey::Point(point)) => point,
        other => panic!("{other:?}: {}", h.editor.message),
    };
    assert!(painted_selected(&mut h, "Point 1"));
    h.click([-80., -80.], Vec2::ZERO);
    assert_eq!(h.editor.shape, None, "empty canvas clears the tree");
    assert!(!painted_selected(&mut h, "Point 1"));
    h.click_label("Point 1");
    assert_eq!(h.editor.shape, Some(ElementKey::Point(point)));
    assert_eq!(h.editor.constraints.selection(), &[point]);
    assert!(h.labels().iter().any(|l| l == "Selected: Point 1"));
    // Clicking the selected row again clears it.
    h.click_label("Point 1");
    assert_eq!(h.editor.shape, None);
}

#[test]
fn a_projected_model_edge_is_listed_as_a_projected_reference() {
    let fixture = crate::external_references::tests::projected_fixture();
    let mut editor = SketchEditor::default();
    editor
        .open_at(
            fixture.session.document(),
            fixture.sketch,
            fixture.session.revision(),
        )
        .unwrap();
    let mut h = Host {
        ctx: egui::Context::default(),
        editor,
        session: fixture.session,
        feed: Default::default(),
        time: 0.,
    };
    h.settle();
    assert!(tab_shown(&mut h, SidebarTab::Sketch));
    let rows = element_rows_painted(&h);
    assert_eq!(
        rows,
        vec![(ElementType::Projected, "Projected reference 1".into(), 0)],
        "{}",
        h.editor.message
    );
    assert!(h.labels().iter().any(|l| l == "Projected reference 1"));
    h.click_label("Projected reference 1");
    assert!(matches!(h.editor.shape, Some(ElementKey::Curve(_))));
    assert!(
        h.labels()
            .iter()
            .any(|l| l == "Selected: Projected reference 1")
    );
}
