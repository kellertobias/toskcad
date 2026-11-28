//! CAD-428: doors, windows and cutouts are listed under the wall that hosts
//! them in the model tree, from the canonical wall source. Selecting one
//! selects it in the 3D view (its reveal faces), the 2D plan (its footprint)
//! and the Inspector; adding, removing and moving it to another wall update
//! the tree, and so does reopening the saved file.
use super::ui_harness::Harness;
use super::*;
use crate::{icons::Icon, model_tree::NodeKind};
use multicad_core::model::{
    self, BodyOutput, Hinge, Justification, NewType, NewTypeKind, Opening, OpeningKind,
    OpeningTypeDefinition, Swing, WallChange, WallDefinition, WallOpeningChange, WallRequest,
    WallTypeDefinition, length,
};
use std::collections::BTreeSet;

struct Room {
    document: Document,
    north: EntityId,
    north_part: EntityId,
    south: EntityId,
    south_part: EntityId,
    door: EntityId,
    window: EntityId,
    window_type: EntityId,
}
/// Two parallel 6 m walls 3 m apart; the north wall hosts a door and a
/// window.
fn room() -> Room {
    let mut session = Session::new(Document::new("Hosted openings").unwrap()).unwrap();
    let actor = || Actor::local(Channel::Ui);
    let wall_type = NewType::new(
        "Block 200",
        NewTypeKind::Wall(WallTypeDefinition {
            thickness: length(200.),
            height: length(2500.),
            justification: Justification::Center,
            material: None,
        }),
    );
    let door_type = NewType::new(
        "Door 900",
        NewTypeKind::Opening(OpeningTypeDefinition {
            kind: OpeningKind::Door {
                hinge: Hinge::Start,
                swing: Swing::Left,
            },
            width: length(900.),
            height: length(2100.),
            sill: length(0.),
        }),
    );
    let window_type = NewType::new(
        "Window 1200",
        NewTypeKind::Opening(OpeningTypeDefinition {
            kind: OpeningKind::Window,
            width: length(1200.),
            height: length(1200.),
            sill: length(900.),
        }),
    );
    let opening = |ty: EntityId, offset: f64| Opening {
        id: EntityId::new(),
        opening_type: ty,
        segment: 0,
        offset: length(offset),
        width: None,
        height: None,
        sill: None,
    };
    let door = opening(door_type.id, 1000.);
    let window = opening(window_type.id, 3500.);
    let mut created = vec![];
    let mut types = vec![wall_type.clone(), door_type, window_type.clone()];
    for (name, y, openings) in [
        ("North wall", 3000., vec![door.clone(), window.clone()]),
        ("South wall", 0., vec![]),
    ] {
        let request = WallRequest::new(
            name,
            WallDefinition {
                layer_placement: None,
                output: BodyOutput {
                    body: EntityId::new(),
                },
                wall_type: wall_type.id,
                path: vec![[length(0.), length(y)], [length(6000.), length(y)]],
                closed: false,
                level: None,
                base_z: length(0.),
                thickness: None,
                height: None,
                justification: None,
                openings,
                cutouts: vec![],
            },
        );
        let prepared = model::prepare_wall_change(
            &session,
            session.revision(),
            actor(),
            &types,
            &WallChange::Create(&request),
        )
        .unwrap();
        session.commit(prepared.prepared).unwrap();
        types.clear();
        created.push((request.feature, request.body));
    }
    Room {
        document: session.document().clone(),
        north: created[0].0,
        north_part: created[0].1,
        south: created[1].0,
        south_part: created[1].1,
        door: door.id,
        window: window.id,
        window_type: window_type.id,
    }
}

/// (row, parent row) for every row of the model tree.
fn tree(h: &mut Harness) -> Vec<(crate::model_tree::TreeRow, Option<EntityId>)> {
    let p = h.app.current_mut();
    let rows = p.sketch.tree_rows_for_tests(&p.session);
    let mut out = vec![];
    let mut parents: Vec<EntityId> = vec![];
    for row in rows {
        parents.truncate(row.depth);
        let parent = row
            .depth
            .checked_sub(1)
            .and_then(|d| parents.get(d))
            .copied();
        parents.push(row.id);
        out.push((row, parent));
    }
    out
}
/// The hosted openings under `part`: (identity, type label, icon, name).
fn children(h: &mut Harness, part: EntityId) -> Vec<(EntityId, &'static str, Icon, String)> {
    tree(h)
        .into_iter()
        .filter(|(row, parent)| *parent == Some(part) && row.kind == NodeKind::Opening)
        .map(|(row, _)| (row.id, row.type_label(), row.icon(), row.name.clone()))
        .collect()
}
fn ready(h: &mut Harness, what: &str) {
    h.until(what, |app| {
        let p = app.current();
        p.scene
            .is_current(p.session.document().id(), p.session.revision())
    });
}
fn apply(h: &mut Harness) {
    h.until("an applicable opening preview", |app| {
        let p = app.current();
        p.sketch.model_can_apply(&p.session, &p.resources)
    });
    let revision = h.app.current().session.revision();
    h.click_label("Apply");
    h.until("the applied opening change", |app| {
        app.current().session.revision() == revision + 1
    });
    ready(h, "the rebuilt walls");
}
/// Every painted line, flattened.
fn lines(shapes: &[egui::epaint::ClippedShape]) -> Vec<egui::epaint::PathShape> {
    fn walk(shape: &egui::Shape, out: &mut Vec<egui::epaint::PathShape>) {
        match shape {
            egui::Shape::Vec(shapes) => shapes.iter().for_each(|s| walk(s, out)),
            egui::Shape::Path(path) => out.push(path.clone()),
            _ => {}
        }
    }
    let mut out = vec![];
    for clipped in shapes {
        walk(&clipped.shape, &mut out);
    }
    out
}

#[test]
fn hosted_doors_and_windows_list_under_their_wall_select_everywhere_and_follow_edits() {
    let r = room();
    let mut h = Harness::new(r.document.clone());
    let source = h.app.current().session.document().clone();

    // Each wall lists the openings its source hosts, with type and glyph.
    assert_eq!(
        children(&mut h, r.north_part),
        vec![
            (r.door, "Door", Icon::Door, "Door 900".into()),
            (r.window, "Window", Icon::Window, "Window 1200".into()),
        ]
    );
    assert!(children(&mut h, r.south_part).is_empty());

    // A click on the door row selects the opening: its reveal faces on the
    // wall in 3D and the Inspector, without changing the document.
    h.click_label("Door 900");
    h.until("the door selection", |app| {
        app.current().sketch.tree_selection_for_tests() == &BTreeSet::from([r.door])
    });
    h.frame(vec![]);
    let solid = h.app.current().sketch.solid_for_tests();
    assert_eq!(solid.selected_openings(), &BTreeSet::from([r.door]));
    assert!(solid.selected_ids().is_empty(), "no Part is selected");
    let faces = solid.selected_member_faces().to_vec();
    assert!(
        faces.len() >= 3 && faces.iter().all(|(instance, _)| *instance == r.north_part),
        "the door's jambs and head on the north wall: {faces:?}"
    );
    let texts = h.texts();
    assert!(texts.iter().any(|t| t == "Door · Door 900"), "{texts:?}");
    assert!(
        texts
            .iter()
            .any(|t| t == "Hosted by North wall · segment 1"),
        "{texts:?}"
    );
    assert_eq!(h.app.current().session.document(), &source);

    // Adding an opening to the south wall lists it there.
    let added = Opening {
        id: EntityId::new(),
        opening_type: r.window_type,
        segment: 0,
        offset: length(500.),
        width: None,
        height: None,
        sill: None,
    };
    {
        let session = &mut h.app.current_mut().session;
        let prepared = model::prepare_wall_opening(
            session,
            session.revision(),
            Actor::local(Channel::Ui),
            r.south,
            &[],
            &WallOpeningChange::AddOpening(added.clone()),
        )
        .unwrap();
        session.commit(prepared.prepared).unwrap();
    }
    ready(&mut h, "the added window");
    assert_eq!(
        children(&mut h, r.south_part)
            .iter()
            .map(|c| c.0)
            .collect::<Vec<_>>(),
        vec![added.id]
    );

    // Edit door → Move to wall: one Apply moves it under the south wall.
    h.app
        .current_mut()
        .sketch
        .set_tree_selection_for_tests(BTreeSet::from([r.door]));
    h.frame(vec![]);
    h.click_label("Edit door");
    h.until("the door form", |app| {
        app.current().sketch.model_interaction_active()
    });
    {
        let p = h.app.current_mut();
        let document = p.session.document().clone();
        p.sketch
            .model_for_tests_mut()
            .rehost_opening_for_tests(&document, Some(r.south))
            .unwrap();
    }
    let texts = h.texts();
    assert!(
        texts.iter().any(|t| t == "Move to wall"),
        "the door form offers Move to wall: {texts:?}"
    );
    apply(&mut h);
    let session = &h.app.current().session;
    assert_eq!(
        session.activity().last().map(|a| a.label.as_str()),
        Some("Move door to South wall")
    );
    assert_eq!(
        model::opening_host(session.document(), r.door)
            .unwrap()
            .wall,
        r.south
    );
    assert_eq!(
        children(&mut h, r.north_part)
            .iter()
            .map(|c| c.0)
            .collect::<Vec<_>>(),
        vec![r.window]
    );
    assert_eq!(
        children(&mut h, r.south_part)
            .iter()
            .map(|c| c.0)
            .collect::<Vec<_>>(),
        vec![added.id, r.door]
    );

    // Remove the window through its form: it leaves the tree.
    h.app
        .current_mut()
        .sketch
        .set_tree_selection_for_tests(BTreeSet::from([r.window]));
    h.frame(vec![]);
    h.click_label("Edit window");
    h.until("the window form", |app| {
        app.current().sketch.model_interaction_active()
    });
    h.click_label("Remove this window on Apply");
    apply(&mut h);
    assert!(model::opening_host(h.app.current().session.document(), r.window).is_none());
    assert!(children(&mut h, r.north_part).is_empty());
    assert!(
        tree(&mut h).iter().all(|(row, _)| row.id != r.window),
        "a removed opening leaves no row"
    );

    // Saving and reopening keeps the hierarchy: it is the wall source.
    let file =
        multicad_core::project_file::ProjectFile::new(h.app.current().session.document().clone())
            .unwrap();
    let mut bytes = vec![];
    multicad_core::project_file::write(&file, &mut bytes, Default::default()).unwrap();
    let reopened = multicad_core::project_file::read(
        &mut bytes.as_slice(),
        Default::default(),
        &multicad_core::project_file::supported_capabilities(),
    )
    .unwrap()
    .project
    .document;
    let before_south = children(&mut h, r.south_part);
    let mut again = Harness::new(reopened);
    assert_eq!(children(&mut again, r.south_part), before_south);
    assert!(children(&mut again, r.north_part).is_empty());
    let _ = r.north;
}

/// The 2D plan emphasizes a selected opening's footprint on its drawn wall,
/// from the canonical wall source; openings of walls not in the plan and
/// unselected openings are not drawn.
#[test]
fn architecture_plan_outlines_the_selected_opening_footprint() {
    let r = room();
    let session = Session::new(r.document.clone()).unwrap();
    let ctx = egui::Context::default();
    let rect = egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(800., 600.));
    let camera = crate::camera::Camera {
        yaw: 0.,
        pitch: std::f32::consts::FRAC_PI_2,
        perspective: false,
        target: [3000., 1500., 0.],
        ..Default::default()
    };
    let paint = |sources: BTreeSet<EntityId>, selected: BTreeSet<EntityId>| {
        let plan = crate::viewport::ArchitecturePlan {
            source_hash: multicad_core::project_file::ContentHash::of(b"plan"),
            sources,
            vectors: vec![],
        };
        let out = ctx.run_ui(
            egui::RawInput {
                screen_rect: Some(rect),
                ..Default::default()
            },
            |ui| {
                plan.paint_openings(ui.painter(), &session, camera, rect, &selected);
            },
        );
        let shapes = lines(&out.shapes);
        out.drop_without_applying_deltas();
        shapes
    };
    let drawn = paint(BTreeSet::from([r.north_part]), BTreeSet::from([r.door]));
    assert_eq!(drawn.len(), 1);
    let footprint = model::opening_footprint(&r.document, session.parameters(), r.door).unwrap();
    // The door spans x 1000–1900 through the 200 mm north wall at y = 3000.
    assert_eq!(
        footprint,
        [
            [1000., 2900.],
            [1900., 2900.],
            [1900., 3100.],
            [1000., 3100.]
        ]
    );
    let expected: Vec<_> = footprint
        .iter()
        .chain(&footprint[..1])
        .map(|p| {
            let q = camera.project([p[0] as f32, p[1] as f32, camera.target[2]], rect);
            egui::pos2(q[0], q[1])
        })
        .collect();
    assert_eq!(drawn[0].points, expected);
    assert!(matches!(drawn[0].stroke.color,
        egui::epaint::ColorMode::Solid(c) if c == egui::Color32::from_rgb(90, 170, 255)));
    assert!(paint(BTreeSet::from([r.south_part]), BTreeSet::from([r.door])).is_empty());
    assert!(paint(BTreeSet::from([r.north_part]), BTreeSet::new()).is_empty());
    let _ = (r.south, r.window);
}
