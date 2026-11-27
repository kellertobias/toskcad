//! CAD-440: the Frame tool driven through the full application by real
//! input: typed segments with units and `=$width`, clicks on the drawing
//! plane, keys, joint badges and the joint editor. Waits are on readiness
//! (native preview ready, model current), never on frame counts.
use super::ui_harness::Harness;
use super::*;
use crate::modeling::ModelingTool;
use multicad_core::{
    constraint_sketch::Dimension,
    document::{Entity, EntityData, Mutation},
    materials::{LibraryDefinition, ProfileShape, Stock, StockKind, prepare_library},
    model::{self, ConnectionKind, NodePosition},
};

const LENGTH: &str = "Segment length (mm)";
const TURN: &str = "Turn from the previous member (deg, counterclockwise)";

fn document() -> Document {
    let mut session = Session::new(Document::new("Frame tool").unwrap()).unwrap();
    let steel = Entity::new(
        "S235 steel",
        EntityData::Material {
            density_kg_m3: Some(7850.),
            color_rgba: [120, 124, 130, 255],
            section_hatch: None,
        },
    );
    let shs = LibraryDefinition::Stock(Stock {
        material: steel.id,
        specification: StockKind::Profile {
            designation: "SHS 40×40×2".into(),
            length_mm: 6000.,
            width_mm: 40.,
            height_mm: 40.,
            wall_mm: Some(2.),
            shape: ProfileShape::Rectangular,
            section: None,
            corner_radius_mm: None,
            web_mm: None,
        },
    })
    .entity(EntityId::new(), "SHS 40×40×2")
    .unwrap();
    let actor = Actor::local(Channel::Ui);
    session
        .commit(prepare_library(&session, 0, actor.clone(), &[steel, shs]).unwrap())
        .unwrap();
    let width = Entity::new(
        "width",
        EntityData::Parameter {
            expression: "1000".into(),
            unit: "mm".into(),
            exposed: true,
            inputs: vec![],
        },
    );
    let p = session
        .prepare(multicad_core::commands::Request {
            expected_revision: session.revision(),
            actor,
            label: "Add width".into(),
            edits: vec![Mutation::Insert(width)],
        })
        .unwrap();
    session.commit(p).unwrap();
    session.document().clone()
}
fn open_tool(h: &mut Harness) {
    h.app
        .invoke(Command::Modeling(ModelingTool::Frame), &[], Channel::Ui)
        .unwrap();
    h.until("the Frame form", |app| app.current().sketch.model_active());
    h.frame(vec![]);
}
fn apply(h: &mut Harness) {
    h.app.invoke(Command::Preview, &[], Channel::Ui).unwrap();
    h.until("the native frame preview", |app| {
        let p = app.current();
        p.sketch.model_can_apply(&p.session, &p.resources)
    });
    h.app.invoke(Command::Apply, &[], Channel::Ui).unwrap();
    h.until("the applied model", |app| {
        !app.current().sketch.model_active()
    });
}
fn the_frame(document: &Document) -> EntityId {
    model::frames(document).next().expect("an applied frame").id
}
fn cut_list(document: &Document) -> Vec<multicad_model::MemberRecord> {
    multicad_model::Compiler::default()
        .compile(document, 0)
        .unwrap()
        .member_cut_list(document)
        .unwrap()
}
/// Draw the closed 1000 × 600 rectangle by typed segments.
fn typed_rectangle(h: &mut Harness) {
    open_tool(h);
    h.type_into_field_and_enter(LENGTH, "=$width");
    h.type_into_field_and_enter(TURN, "90");
    h.type_into_field_and_enter(LENGTH, "60 cm");
    h.type_into_field_and_enter(LENGTH, "1 m");
    h.click_label("Close loop");
}

#[test]
fn a_typed_rectangle_applies_one_transaction_with_exact_mitred_members() {
    let mut h = Harness::new(document());
    let before = h.app.current().session.document().clone();
    let revision = h.revision();
    typed_rectangle(&mut h);
    // Inner placement is the default; the badges show four mitres.
    let texts = h.texts();
    assert_eq!(texts.iter().filter(|t| *t == "M").count(), 4, "{texts:?}");
    apply(&mut h);
    assert_eq!(h.revision(), revision + 1, "one transaction");
    let document = h.app.current().session.document().clone();
    let frame = model::read_frame(&document, the_frame(&document)).unwrap();
    assert_eq!(frame.members.len(), 4);
    // The typed =$width is kept as an expression.
    let first = frame
        .nodes
        .iter()
        .find_map(|n| match &n.position {
            NodePosition::Polar { length, .. } => Some(length.clone()),
            _ => None,
        })
        .unwrap();
    assert!(
        matches!(&first, Dimension::Expression { expression } if expression.text == "=$width"),
        "{first:?}"
    );
    let list = cut_list(&document);
    let mut lengths: Vec<f64> = list.iter().map(|r| r.length_mm.round()).collect();
    lengths.sort_by(f64::total_cmp);
    assert_eq!(lengths, [600., 600., 1000., 1000.]);
    for r in &list {
        assert!(r.issues.is_empty(), "{:?}", r.issues);
        for end in &r.ends {
            assert_eq!(end.kind, Some(ConnectionKind::Miter));
            assert!((end.angle_deg - 45.).abs() < 1e-6);
        }
    }
    // Undo and Redo restore exactly.
    let after = document;
    h.app.invoke(Command::Undo, &[], Channel::Ui).unwrap();
    assert_eq!(h.app.current().session.document(), &before);
    h.app.invoke(Command::Redo, &[], Channel::Ui).unwrap();
    assert_eq!(h.app.current().session.document(), &after);
}

#[test]
fn clicks_draw_and_close_a_chain_backspace_and_esc_edit_it_and_discard_keeps_the_source() {
    let mut h = Harness::new(document());
    let source = h.app.current().session.document().clone();
    open_tool(&mut h);
    h.click_label("Draw in viewport");
    let viewport = h.app.current().sketch.viewport_rect.expect("viewport");
    let at = [
        viewport.center() + egui::vec2(-120., 60.),
        viewport.center() + egui::vec2(120., 60.),
        viewport.center() + egui::vec2(0., -90.),
        viewport.center() + egui::vec2(-60., -150.),
    ];
    for p in &at {
        h.click_at(*p);
    }
    let members = |h: &mut Harness| {
        h.texts()
            .into_iter()
            .find(|t| t.ends_with(" nodes"))
            .unwrap_or_default()
    };
    assert!(
        members(&mut h).starts_with("3 members"),
        "{}",
        members(&mut h)
    );
    // Backspace removes the last node and the member that reached it.
    h.key(egui::Key::Backspace);
    assert!(
        members(&mut h).starts_with("2 members"),
        "{}",
        members(&mut h)
    );
    // Clicking the first node closes the loop.
    h.click_at(at[0]);
    assert!(
        members(&mut h).starts_with("3 members"),
        "{}",
        members(&mut h)
    );
    assert!(
        h.texts().iter().any(|t| t == "Draw in viewport"),
        "closing ends the chain"
    );
    // Esc while drawing cancels the pending segment, keeping the form.
    h.click_label("Draw in viewport");
    h.key(egui::Key::Escape);
    assert!(h.modeling(), "Esc cancels the segment, not the form");
    assert!(h.texts().iter().any(|t| t == "Draw in viewport"));
    // Discard restores the source.
    h.app.invoke(Command::Discard, &[], Channel::Ui).unwrap();
    h.until("the discarded form", |app| {
        !app.current().sketch.model_active()
    });
    assert_eq!(h.app.current().session.document(), &source);
}

#[test]
fn a_reopened_frame_shows_joint_badges_and_the_joint_editor_overrides_a_corner() {
    let mut h = Harness::new(document());
    typed_rectangle(&mut h);
    apply(&mut h);
    let document = h.app.current().session.document().clone();
    let frame_id = the_frame(&document);
    let ids = model::frame_member_features(&document, frame_id);
    // Reopen through one of its member features.
    {
        let p = h.app.current_mut();
        p.sketch
            .begin_feature_edit(&p.session, &p.resources, &p.scene, ids[0].1, Channel::Ui)
            .unwrap();
    }
    h.until("the reopened frame", |app| {
        app.current().sketch.model_active()
    });
    assert_eq!(h.texts().iter().filter(|t| *t == "M").count(), 4);
    // Click the first corner on the plane: its joint editor opens.
    let frame = model::read_frame(&document, frame_id).unwrap();
    let corner = frame.members[0].start;
    let world = frame
        .evaluate(h.app.current().session.parameters())
        .unwrap()
        .nodes[&corner];
    let viewport = h.app.current().sketch.viewport_rect.unwrap();
    let s = h
        .app
        .current()
        .sketch
        .camera
        .project(world.map(|v| v as f32), viewport);
    h.click_at(egui::pos2(s[0], s[1]));
    assert!(h.texts().iter().any(|t| t == "Joint"), "joint editor opens");
    h.click_label("Butt");
    assert!(h.texts().iter().any(|t| t == "B→*"), "{:?}", h.texts());
    apply(&mut h);
    let after = h.app.current().session.document().clone();
    assert_eq!(
        model::frame_member_features(&after, frame_id),
        ids,
        "member feature ids are kept"
    );
    let stored = model::read_frame(&after, frame_id).unwrap();
    assert_eq!(
        stored.joint(corner).and_then(|j| j.kind),
        Some(ConnectionKind::Butt)
    );
    let butts = cut_list(&after)
        .iter()
        .flat_map(|r| r.ends.clone())
        .filter(|e| e.kind == Some(ConnectionKind::Butt))
        .count();
    assert_eq!(butts, 2, "both ends at the corner are butted");
}

fn outside(document: &Document) -> ([f64; 2], [f64; 2]) {
    let model = multicad_model::Compiler::default()
        .compile(document, 0)
        .unwrap();
    let (mut lo, mut hi) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
    for r in model.member_cut_list(document).unwrap() {
        let b = model
            .solid(model.snapshot().bodies[&r.part])
            .unwrap()
            .bounds()
            .unwrap();
        for k in 0..2 {
            lo[k] = lo[k].min(b.min[k]);
            hi[k] = hi[k].max(b.max[k]);
        }
    }
    (lo, hi)
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-6, "{a} != {b}");
}
fn reopen(h: &mut Harness) {
    let document = h.app.current().session.document().clone();
    let frame = the_frame(&document);
    let feature = model::frame_member_features(&document, frame)[0].1;
    let p = h.app.current_mut();
    p.sketch
        .begin_feature_edit(&p.session, &p.resources, &p.scene, feature, Channel::Ui)
        .unwrap();
    h.until("the reopened frame", |app| {
        app.current().sketch.model_active()
    });
}

#[test]
fn inner_and_outer_put_the_sections_inside_or_outside_the_drawn_loop() {
    let mut h = Harness::new(document());
    typed_rectangle(&mut h);
    apply(&mut h);
    let (lo, hi) = outside(h.app.current().session.document());
    near(lo[0], 0.);
    near(hi[0], 1000.);
    near(hi[1], 600.);
    reopen(&mut h);
    h.click_label("Inner edge on the line");
    h.click_label("Outer edge on the line");
    apply(&mut h);
    let (lo, hi) = outside(h.app.current().session.document());
    near(lo[0], -40.);
    near(lo[1], -40.);
    near(hi[0], 1040.);
    near(hi[1], 640.);
}

#[test]
fn joints_a_section_cannot_take_and_clicks_off_the_plane_are_refused_with_the_reason() {
    let mut h = Harness::new(document());
    open_tool(&mut h);
    h.type_into_field_and_enter(LENGTH, "500");
    h.type_into_field_and_enter(TURN, "170");
    h.type_into_field_and_enter(LENGTH, "500");
    let texts = h.texts();
    assert!(texts.iter().any(|t| t == "!"), "a red badge: {texts:?}");
    // The joint editor at that node (500, 0) names the reason.
    h.app.current_mut().sketch.camera.target = [500., 0., 0.];
    h.frame(vec![]);
    let viewport = h.app.current().sketch.viewport_rect.unwrap();
    let s = h
        .app
        .current()
        .sketch
        .camera
        .project([500., 0., 0.], viewport);
    h.click_at(egui::pos2(s[0], s[1]));
    assert!(
        h.texts().iter().any(|t| t.contains("at least 15°")),
        "{:?}",
        h.texts()
    );
    h.app.invoke(Command::Discard, &[], Channel::Ui).unwrap();
    h.until("the discarded form", |app| {
        !app.current().sketch.model_active()
    });

    // A straight top view looks along the Front (XZ) plane: clicks cannot
    // hit it and are refused with the reason.
    open_tool(&mut h);
    h.click_label("Top (XY)");
    h.click_label("Front (XZ)");
    {
        let camera = &mut h.app.current_mut().sketch.camera;
        camera.yaw = crate::camera::Camera::PLAN_YAW;
        camera.pitch = std::f32::consts::FRAC_PI_2;
    }
    h.click_label("Draw in viewport");
    let viewport = h.app.current().sketch.viewport_rect.unwrap();
    h.click_at(viewport.center());
    assert!(
        h.texts()
            .iter()
            .any(|t| t.starts_with("That click does not hit the drawing plane")),
        "{:?}",
        h.texts()
    );
}

#[test]
fn new_profile_stock_is_created_in_the_project_without_leaving_the_tool() {
    let mut h = Harness::new(document());
    open_tool(&mut h);
    let stocks = |h: &Harness| {
        h.app
            .current()
            .session
            .document()
            .entities()
            .filter(|e| {
                matches!(
                    multicad_core::materials::read(e).ok().flatten(),
                    Some(LibraryDefinition::Stock(_))
                )
            })
            .count()
    };
    assert_eq!(stocks(&h), 1);
    // The picker lists the Project stock with its plain name.
    h.click_label("Square tube 40×40×2 mm");
    assert!(
        h.texts().iter().any(|t| t == "SHS 40×40×2"),
        "the catalogue name stays"
    );
    h.click_label("New profile stock…");
    h.click_label("Apply");
    h.until("the new stock", |app| {
        app.current()
            .session
            .document()
            .entities()
            .filter(|e| {
                matches!(
                    multicad_core::materials::read(e).ok().flatten(),
                    Some(LibraryDefinition::Stock(_))
                )
            })
            .count()
            == 2
    });
    assert!(h.modeling(), "the Frame tool stays open");
    assert!(h.texts().iter().any(|t| t == "Draw in viewport"));
}

#[test]
fn moving_a_node_of_an_applied_frame_keeps_member_and_part_ids() {
    let mut h = Harness::new(document());
    typed_rectangle(&mut h);
    apply(&mut h);
    let document = h.app.current().session.document().clone();
    let frame_id = the_frame(&document);
    let ids = model::frame_member_features(&document, frame_id);
    let parts: Vec<EntityId> = ids
        .iter()
        .map(|(_, f)| model::read_member(&document, *f).unwrap().output.body)
        .collect();
    reopen(&mut h);
    let frame = model::read_frame(&document, frame_id).unwrap();
    let corner = frame.members[0].start;
    let project = |h: &Harness, p: [f64; 3]| {
        let viewport = h.app.current().sketch.viewport_rect.unwrap();
        let s = h
            .app
            .current()
            .sketch
            .camera
            .project(p.map(|v| v as f32), viewport);
        egui::pos2(s[0], s[1])
    };
    let world = frame
        .evaluate(h.app.current().session.parameters())
        .unwrap()
        .nodes[&corner];
    h.click_at(project(&h, world));
    h.click_label("Move node");
    // The new place: 100 mm along −u, −v, snapped to whole millimetres.
    h.click_at(project(&h, [-100., -100., 0.]));
    apply(&mut h);
    let after = h.app.current().session.document().clone();
    assert_eq!(model::frame_member_features(&after, frame_id), ids);
    let parts_after: Vec<EntityId> = ids
        .iter()
        .map(|(_, f)| model::read_member(&after, *f).unwrap().output.body)
        .collect();
    assert_eq!(parts_after, parts);
    let moved = model::read_frame(&after, frame_id).unwrap();
    let p = moved
        .evaluate(h.app.current().session.parameters())
        .unwrap()
        .nodes[&corner];
    assert!(p[0] < -50. && p[1] < -50., "{p:?}");
}
