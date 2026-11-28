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
    settle_view(&mut h);
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
    settle_view(h);
}
/// Run frames until the camera and viewport stop changing (reopening may
/// fit the view), so clicks land where the test projects them.
fn settle_view(h: &mut Harness) {
    // The first evaluated bounds of a document frame the view once (CAD-175);
    // wait for the scene of the current revision so that fit has happened.
    h.until("the evaluated model of the current revision", |app| {
        let p = app.current();
        p.scene
            .is_current(p.session.document().id(), p.session.revision())
    });
    let deadline = std::time::Instant::now() + std::time::Duration::from_secs(20);
    let mut previous = None;
    let mut still = 0;
    while still < 3 {
        assert!(
            std::time::Instant::now() < deadline,
            "timed out waiting for the view to settle"
        );
        h.frame(vec![]);
        let sketch = &h.app.current().sketch;
        let now = (sketch.camera, sketch.viewport_rect);
        if previous == Some(now) {
            still += 1;
        } else {
            still = 0;
        }
        previous = Some(now);
    }
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

#[test]
fn inches_a_gap_set_in_the_joint_editor_and_reversing_a_member_apply_as_typed() {
    // An open chain typed in inches.
    let mut h = Harness::new(document());
    open_tool(&mut h);
    h.type_into_field_and_enter(LENGTH, "10 in");
    apply(&mut h);
    let list = cut_list(h.app.current().session.document());
    assert_eq!(list.len(), 1);
    near(list[0].length_mm, 254.);

    // A rectangle with a 2 mm gap typed at one corner.
    let mut h = Harness::new(document());
    typed_rectangle(&mut h);
    apply(&mut h);
    reopen(&mut h);
    let document = h.app.current().session.document().clone();
    let frame_id = the_frame(&document);
    let frame = model::read_frame(&document, frame_id).unwrap();
    // The corner at the origin, in view.
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
    eprintln!("DEBUG {:?}", h.texts());
    h.type_into_field_and_enter("Weld gap (mm, empty = frame default 0 mm)", "2");
    assert!(
        h.texts().iter().any(|t| t == "M*"),
        "the badge shows the override"
    );
    apply(&mut h);
    let after = h.app.current().session.document().clone();
    let gapped: Vec<f64> = cut_list(&after)
        .iter()
        .map(|r| r.length_mm)
        .filter(|l| (l.round() - l).abs() > 1e-3)
        .collect();
    assert_eq!(gapped.len(), 2, "both members at the corner are shortened");
    for l in gapped {
        near(l.fract(), (1000f64 - std::f64::consts::SQRT_2).fract());
    }

    // Reversing a member keeps its section on the same side.
    reopen(&mut h);
    let front = model::read_frame(&after, frame_id).unwrap().members[0].clone();
    let nodes = model::read_frame(&after, frame_id)
        .unwrap()
        .evaluate(h.app.current().session.parameters())
        .unwrap()
        .nodes;
    // A point on the front rail near the origin selects it.
    let a = nodes[&front.start];
    let mid: [f64; 3] = [a[0] + 60., a[1], a[2]];
    let viewport = h.app.current().sketch.viewport_rect.unwrap();
    let s = h
        .app
        .current()
        .sketch
        .camera
        .project(mid.map(|v| v as f32), viewport);
    h.click_at(egui::pos2(s[0], s[1]));
    h.click_label("Reverse");
    apply(&mut h);
    let reversed = h.app.current().session.document().clone();
    let stored = model::read_frame(&reversed, frame_id).unwrap();
    assert_eq!(
        stored.members[0].start, front.end,
        "the member was reversed"
    );
    let (lo, hi) = outside(&reversed);
    near(lo[0], 0.);
    near(lo[1], 0.);
    near(hi[0], 1000.);
    near(hi[1], 600.);
}

#[test]
fn dragging_a_node_updates_members_live_and_apply_keeps_member_ids() {
    let mut h = Harness::new(document());
    open_tool(&mut h);
    h.type_into_field_and_enter(LENGTH, "200");
    h.type_into_field_and_enter(TURN, "90");
    h.type_into_field_and_enter(LENGTH, "100");
    apply(&mut h);
    let document = h.app.current().session.document().clone();
    let frame_id = the_frame(&document);
    let ids = model::frame_member_features(&document, frame_id);
    reopen(&mut h);
    let frame = model::read_frame(&document, frame_id).unwrap();
    let tip = frame.members[1].end;
    let screen = |h: &Harness, p: [f64; 3]| {
        let viewport = h.app.current().sketch.viewport_rect.unwrap();
        let s = h
            .app
            .current()
            .sketch
            .camera
            .project(p.map(|v| v as f32), viewport);
        egui::pos2(s[0], s[1])
    };
    let from = screen(&h, [200., 100., 0.]);
    let to = screen(&h, [250., 220., 0.]);
    let button = |pos, pressed| egui::Event::PointerButton {
        pos,
        button: egui::PointerButton::Primary,
        pressed,
        modifiers: egui::Modifiers::NONE,
    };
    h.frame(vec![egui::Event::PointerMoved(from)]);
    h.frame(vec![button(from, true)]);
    for k in 1..=10 {
        let p = from + (to - from) * (k as f32 / 10.);
        h.frame(vec![egui::Event::PointerMoved(p)]);
    }
    // While dragging, the second member's live length follows the node:
    // from (200, 0) to (250, 220).
    assert!(h.texts().iter().any(|t| t == "226 mm"), "{:?}", h.texts());
    h.frame(vec![button(to, false)]);
    apply(&mut h);
    let after = h.app.current().session.document().clone();
    assert_eq!(
        model::frame_member_features(&after, frame_id),
        ids,
        "ids kept"
    );
    let p = model::read_frame(&after, frame_id)
        .unwrap()
        .evaluate(h.app.current().session.parameters())
        .unwrap()
        .nodes[&tip];
    near(p[0], 250.);
    near(p[1], 220.);
    assert!(
        cut_list(&after).iter().all(|r| r.issues.is_empty()),
        "both members rebuild with their new joint"
    );
}

#[test]
fn a_frame_drawn_on_a_level_sits_at_its_elevation() {
    let mut session = Session::new(document()).unwrap();
    let (_, prepared) = model::prepare_level(
        &session,
        session.revision(),
        Actor::local(Channel::Ui),
        "First floor",
        &model::LevelDefinition {
            elevation: model::length(3000.),
            height: model::length(2800.),
            slab_thickness: model::length(200.),
        },
    )
    .unwrap();
    session.commit(prepared).unwrap();
    let choices = crate::modeling::frame_plane_choices(session.document());
    assert!(
        choices
            .iter()
            .any(|(label, plane)| label == "Level First floor"
                && plane.origin[2] == model::length(3000.))
    );
    let mut h = Harness::new(session.document().clone());
    open_tool(&mut h);
    h.click_label("Top (XY)");
    h.click_label("Level First floor");
    h.type_into_field_and_enter(LENGTH, "500");
    apply(&mut h);
    let document = h.app.current().session.document().clone();
    let model = multicad_model::Compiler::default()
        .compile(&document, 0)
        .unwrap();
    let r = &model.member_cut_list(&document).unwrap()[0];
    let b = model
        .solid(model.snapshot().bodies[&r.part])
        .unwrap()
        .bounds()
        .unwrap();
    near((b.min[2] + b.max[2]) / 2., 3000.);
}

#[test]
fn a_picked_planar_face_becomes_the_drawing_plane() {
    let mut session = Session::new(document()).unwrap();
    let request = model::PrimitiveRequest::new(
        "Block",
        model::PrimitiveParameters::Box {
            size: [model::length(30.), model::length(20.), model::length(10.)],
            frame: Default::default(),
        },
    );
    let prepared = model::prepare_primitive(
        &session,
        session.revision(),
        Actor::local(Channel::Ui),
        &request,
    )
    .unwrap()
    .prepared;
    session.commit(prepared).unwrap();
    let mut h = Harness::new(session.document().clone());
    // Where the box's top face is, from its own solid.
    let top = {
        let document = h.app.current().session.document().clone();
        let model = multicad_model::Compiler::default()
            .compile(&document, 0)
            .unwrap();
        let b = model
            .solid(model.snapshot().bodies[&request.body])
            .unwrap()
            .bounds()
            .unwrap();
        [
            (b.min[0] + b.max[0]) / 2.,
            (b.min[1] + b.max[1]) / 2.,
            b.max[2],
        ]
    };
    open_tool(&mut h);
    h.click_label("Pick a planar face");
    let viewport = h.app.current().sketch.viewport_rect.unwrap();
    let s = h
        .app
        .current()
        .sketch
        .camera
        .project(top.map(|v| v as f32), viewport);
    h.click_at(egui::pos2(s[0], s[1]));
    assert!(
        h.texts().iter().any(|t| t == "Face of Block"),
        "{:?}",
        h.texts()
    );
    h.type_into_field_and_enter(LENGTH, "20");
    apply(&mut h);
    let document = h.app.current().session.document().clone();
    let frame = model::read_frame(&document, the_frame(&document)).unwrap();
    // The plane lies in the top face: its normal is ±Z at the top height.
    assert!(
        (frame.plane.normal[2].abs() - 1.).abs() < 1e-9,
        "{:?}",
        frame.plane
    );
    let model = multicad_model::Compiler::default()
        .compile(&document, 0)
        .unwrap();
    let r = model
        .member_cut_list(&document)
        .unwrap()
        .into_iter()
        .find(|r| matches!(r.path, model::MemberPath::FrameMember { .. }))
        .unwrap();
    let b = model
        .solid(model.snapshot().bodies[&r.part])
        .unwrap()
        .bounds()
        .unwrap();
    near((b.min[2] + b.max[2]) / 2., top[2]);
}

/// The frame tool document with a table frame: the inner 1000 × 600
/// rectangle of SHS 40×40×2 and four 700 mm legs down from its corners,
/// built through the same typed transaction the tool commits (drawing legs
/// in 3D is CAD-441). Returns the document and the rectangle's corners.
fn table_document() -> (Document, Vec<EntityId>) {
    use model::{ChainSegment, FrameDefinition, FrameMemberRequest, FramePlane, FrameRequest};
    let mut session = Session::new(document()).unwrap();
    let stock = session
        .document()
        .entities()
        .find(|e| e.name == "SHS 40×40×2")
        .unwrap()
        .id;
    let mut frame = FrameDefinition::new(FramePlane::xy());
    frame.side = model::FrameSide::Inner;
    let start = frame.add_node(NodePosition::Absolute {
        u: model::length(0.),
        v: model::length(0.),
        w: None,
    });
    let seg = |l: f64, a: f64, rise: Option<f64>| ChainSegment {
        length: model::length(l),
        angle: model::angle(a),
        rise: rise.map(model::length),
    };
    let parameters = session.parameters().clone();
    let rails = frame
        .add_chain(
            &parameters,
            start,
            &[
                seg(1000., 0., None),
                seg(600., 90., None),
                seg(1000., 180., None),
            ],
            true,
        )
        .unwrap();
    let corners: Vec<EntityId> = rails
        .iter()
        .map(|m| frame.member(*m).unwrap().start)
        .collect();
    let mut members: Vec<FrameMemberRequest> = rails
        .iter()
        .zip(["Front", "Right", "Back", "Left"])
        .map(|(m, n)| FrameMemberRequest::new(*m, n, stock))
        .collect();
    for (i, corner) in corners.iter().enumerate() {
        let leg = frame
            .add_chain(&parameters, *corner, &[seg(0., 0., Some(-700.))], false)
            .unwrap()[0];
        members.push(FrameMemberRequest::new(
            leg,
            format!("Leg {}", i + 1),
            stock,
        ));
    }
    let prepared = model::prepare_frame(
        &session,
        session.revision(),
        Actor::local(Channel::Ui),
        &FrameRequest {
            frame: EntityId::new(),
            name: "Table".into(),
            definition: frame,
            members,
        },
    )
    .unwrap()
    .prepared;
    session.commit(prepared).unwrap();
    (session.document().clone(), corners)
}
/// The Part of the leg standing on `corner`.
fn leg_part(document: &Document, corner: EntityId) -> EntityId {
    let frame = the_frame(document);
    let definition = model::read_frame(document, frame).unwrap();
    let leg = definition
        .members
        .iter()
        .filter(|m| m.start == corner)
        .find(|m| {
            definition
                .node(m.end)
                .is_some_and(|n| matches!(n.position, NodePosition::Polar { .. }))
                && !definition.members.iter().any(|o| o.start == m.end)
        })
        .unwrap()
        .id;
    let feature = model::frame_member_features(document, frame)
        .into_iter()
        .find(|(m, _)| *m == leg)
        .unwrap()
        .1;
    model::read_member(document, feature).unwrap().output.body
}
/// Wait for the native preview of the current form and return the bounds
/// of `part` in the scene the viewport shows.
fn previewed_bounds(h: &mut Harness, part: EntityId) -> [[f64; 3]; 2] {
    h.until("the native corner preview", |app| {
        let p = app.current();
        p.sketch.model_can_apply(&p.session, &p.resources)
    });
    let p = h.app.current();
    let scene = p
        .sketch
        .model_for_tests()
        .display_scene(&p.session, &p.resources)
        .expect("the preview scene is shown");
    let [lo, hi] = scene
        .native_instance_bounds(part)
        .expect("the previewed leg has native bounds");
    [lo, hi]
}
fn click_node(h: &mut Harness, document: &Document, node: EntityId) {
    let frame = model::read_frame(document, the_frame(document)).unwrap();
    let world = frame
        .evaluate(h.app.current().session.parameters())
        .unwrap()
        .nodes[&node];
    let viewport = h.app.current().sketch.viewport_rect.unwrap();
    let s = h
        .app
        .current()
        .sketch
        .camera
        .project(world.map(|v| v as f32), viewport);
    h.click_at(egui::pos2(s[0], s[1]));
}
const LEG_THROUGH: &str = "Leg through, rails butt onto it";
const MITRE_ON_LEG: &str = "Rails mitred, leg below";
const MITRE_AROUND_LEG: &str = "Rails mitred, leg outside";

#[test]
fn a_leg_corner_dropdown_cycles_previews_and_applies_its_construction_in_one_step() {
    let (document, corners) = table_document();
    let mut h = Harness::new(document.clone());
    let leg = leg_part(&document, corners[0]);
    let before = h.app.current().session.document().clone();
    let revision = h.revision();
    reopen(&mut h);
    // A plain rectangle corner of the table is a leg corner: clicking it
    // opens the dropdown there with every construction, leg through first.
    click_node(&mut h, &document, corners[0]);
    let texts = h.texts();
    for label in ["Corner", LEG_THROUGH, MITRE_ON_LEG, MITRE_AROUND_LEG] {
        assert!(texts.iter().any(|t| t == label), "{label}: {texts:?}");
    }
    assert!(
        texts.iter().any(|t| t == &format!("Corner: {LEG_THROUGH}")),
        "{texts:?}"
    );
    // The highlighted entry previews in 3D: the leg runs to the rail tops.
    let [lo, hi] = previewed_bounds(&mut h, leg);
    assert!((hi[2] - 20.).abs() < 1e-3, "{hi:?}");
    assert!(lo[0].abs() < 1e-3 && lo[1].abs() < 1e-3, "{lo:?}");

    // Down: rails mitred, leg below; the leg stops under the rails.
    h.key(egui::Key::ArrowDown);
    assert!(
        h.texts()
            .iter()
            .any(|t| t == &format!("Corner: {MITRE_ON_LEG}"))
    );
    let [_, hi] = previewed_bounds(&mut h, leg);
    assert!((hi[2] + 20.).abs() < 1e-3, "{hi:?}");
    // Down: rails mitred, leg outside; the leg is centred on the node.
    h.key(egui::Key::ArrowDown);
    let [lo, hi] = previewed_bounds(&mut h, leg);
    assert!(
        (lo[0] + 20.).abs() < 1e-3 && (lo[1] + 20.).abs() < 1e-3,
        "{lo:?}"
    );
    assert!((hi[2] - 20.).abs() < 1e-3, "{hi:?}");
    // Down wraps to the first entry; Up comes back to the last.
    h.key(egui::Key::ArrowDown);
    let [_, hi] = previewed_bounds(&mut h, leg);
    assert!((hi[2] - 20.).abs() < 1e-3, "{hi:?}");
    h.key(egui::Key::ArrowUp);
    h.key(egui::Key::ArrowUp);
    assert!(
        h.texts()
            .iter()
            .any(|t| t == &format!("Corner: {MITRE_ON_LEG}"))
    );
    previewed_bounds(&mut h, leg);

    // Enter applies the highlighted construction as the node's override.
    h.key(egui::Key::Enter);
    h.until("the applied corner", |app| {
        app.current().session.revision() > revision
    });
    assert_eq!(h.revision(), revision + 1, "one transaction");
    let after = h.app.current().session.document().clone();
    let stored = model::read_frame(&after, the_frame(&after)).unwrap();
    assert_eq!(
        stored.joint(corners[0]).and_then(|j| j.corner),
        Some(model::FrameCorner::MitreOnLeg)
    );
    assert!(stored.joints.iter().all(|j| j.node == corners[0]));
    let list = cut_list(&after);
    let legs: Vec<f64> = list
        .iter()
        .filter(|r| r.name.starts_with("Leg"))
        .map(|r| r.length_mm.round())
        .collect();
    assert_eq!(legs.iter().filter(|l| **l == 680.).count(), 1, "{legs:?}");
    assert_eq!(legs.iter().filter(|l| **l == 720.).count(), 3, "{legs:?}");
    // One Undo step restores the table exactly.
    h.app.invoke(Command::Undo, &[], Channel::Ui).unwrap();
    assert_eq!(h.app.current().session.document(), &before);
}

#[test]
fn esc_closes_the_corner_dropdown_and_restores_the_joint() {
    let (document, corners) = table_document();
    let mut h = Harness::new(document.clone());
    let leg = leg_part(&document, corners[0]);
    let revision = h.revision();
    reopen(&mut h);
    click_node(&mut h, &document, corners[0]);
    h.key(egui::Key::ArrowDown);
    let [_, hi] = previewed_bounds(&mut h, leg);
    assert!((hi[2] + 20.).abs() < 1e-3, "{hi:?}");
    h.key(egui::Key::Escape);
    let texts = h.texts();
    assert!(!texts.iter().any(|t| t == MITRE_ON_LEG), "{texts:?}");
    assert!(h.app.current().sketch.model_active(), "the tool stays open");
    // The preview is back to the automatic corner and nothing was written.
    let [_, hi] = previewed_bounds(&mut h, leg);
    assert!((hi[2] - 20.).abs() < 1e-3, "{hi:?}");
    assert_eq!(h.revision(), revision);
    // Choosing a joint kind in the editor replaces a corner construction.
    click_node(&mut h, &document, corners[0]);
    h.key(egui::Key::ArrowDown);
    h.click_label("Automatic");
    assert!(!h.texts().iter().any(|t| t == MITRE_ON_LEG));
    apply(&mut h);
    let after = h.app.current().session.document().clone();
    let stored = model::read_frame(&after, the_frame(&after)).unwrap();
    assert!(
        stored.joint(corners[0]).is_none_or(|j| j.corner.is_none()),
        "{:?}",
        stored.joints
    );
}
