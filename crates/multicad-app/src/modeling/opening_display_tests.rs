//! CAD-429: a wall's doors and windows are drawn in the 3D scene inside
//! their openings, not only as an empty cutout, from the same canonical
//! definition as the 2D plan; editing them updates both views.
use super::*;
use multicad_core::{
    commands::{Actor, Channel, Session},
    drawings::VectorCurve,
    model::{
        BodyOutput, Hinge, Justification, NewType, NewTypeKind, Opening, OpeningKind,
        OpeningTypeDefinition, Swing, WallChange, WallDefinition, WallOpeningChange, WallRequest,
        WallTypeDefinition, length, prepare_wall_change, prepare_wall_opening,
    },
};
use std::time::{Duration, Instant};

struct Room {
    session: Session,
    wall: EntityId,
    body: EntityId,
    door: Opening,
    window: EntityId,
}
/// A 4 m wall along +X, 200 mm thick and centered on y = 0, with a 900 ×
/// 2100 door at x = 500 (hinged at its start jamb, opening to the +Y face)
/// and a 1200 × 1200 window at x = 2200 on a 900 mm sill.
fn room() -> Room {
    let mut session = Session::new(Document::new("Door in 3D").unwrap()).unwrap();
    let wall_type = NewType::new(
        "Wall 200",
        NewTypeKind::Wall(WallTypeDefinition {
            thickness: length(200.),
            height: length(2500.),
            justification: Justification::Center,
            material: None,
        }),
    );
    let door_type = NewType::new(
        "Door",
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
        "Window",
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
    let door = opening(door_type.id, 500.);
    let window = opening(window_type.id, 2200.);
    let request = WallRequest::new(
        "Wall with door",
        WallDefinition {
            output: BodyOutput {
                body: EntityId::new(),
            },
            wall_type: wall_type.id,
            layer_placement: None,
            path: vec![[length(0.), length(0.)], [length(4000.), length(0.)]],
            closed: false,
            level: None,
            base_z: length(0.),
            thickness: None,
            height: None,
            justification: None,
            cutouts: vec![],
            openings: vec![door.clone(), window.clone()],
        },
    );
    let prepared = prepare_wall_change(
        &session,
        session.revision(),
        Actor::local(Channel::Ui),
        &[wall_type, door_type, window_type],
        &WallChange::Create(&request),
    )
    .unwrap();
    session.commit(prepared.prepared).unwrap();
    Room {
        session,
        wall: request.feature,
        body: request.body,
        door,
        window: window.id,
    }
}
fn publish(session: &Session) -> Arc<viewport::Scene> {
    let mut backend = Backend::new();
    backend
        .request(
            session.document(),
            SceneStamp {
                document: session.document().id(),
                revision: session.revision(),
                generation: 1,
            },
        )
        .unwrap();
    let deadline = Instant::now() + Duration::from_secs(20);
    loop {
        if let Some(reply) = backend.poll() {
            return reply.result.unwrap();
        }
        assert!(Instant::now() < deadline, "Native backend did not publish");
        std::thread::sleep(Duration::from_millis(2));
    }
}
/// The centroid of every drawn triangle, with the identity that draws it:
/// scene instances and display elements alike.
fn drawn(scene: &viewport::Scene) -> Vec<(EntityId, [f64; 3])> {
    let mut out = vec![];
    let mut add = |id: EntityId, geometry: &viewport::Geometry, transform: viewport::Transform| {
        for t in geometry.indices().chunks(3) {
            let p = t.iter().map(|i| geometry.vertices()[*i as usize]);
            let c = p.fold([0.; 3], |a, v| [a[0] + v[0], a[1] + v[1], a[2] + v[2]]);
            out.push((id, transform.point(c.map(|v| v / 3.))));
        }
    };
    for instance in scene.instances() {
        add(
            instance.id,
            &scene.geometries()[&instance.geometry],
            instance.transform,
        );
    }
    for element in scene.display_elements() {
        add(element.id, &element.geometry, Default::default());
    }
    out
}
fn inside(p: [f64; 3], lo: [f64; 3], hi: [f64; 3]) -> bool {
    (0..3).all(|i| p[i] > lo[i] && p[i] < hi[i])
}
/// The extent drawn by one identity.
fn extent(points: &[(EntityId, [f64; 3])], id: EntityId) -> ([f64; 3], [f64; 3]) {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for (_, p) in points.iter().filter(|(owner, _)| *owner == id) {
        for i in 0..3 {
            lo[i] = lo[i].min(p[i]);
            hi[i] = hi[i].max(p[i]);
        }
    }
    (lo, hi)
}

#[test]
fn a_wall_door_and_window_are_drawn_in_3d_inside_their_openings() {
    let r = room();
    let scene = publish(&r.session);
    let points = drawn(&scene);
    // The door opening x 500–1400, through the wall, floor to head; the
    // window opening x 2200–3400, sill 900 to head 2100.
    let door_void = ([501., -99., 1.], [1399., 99., 2099.]);
    let window_void = ([2201., -99., 901.], [3399., 99., 2099.]);
    let in_door: Vec<_> = points
        .iter()
        .filter(|(_, p)| inside(*p, door_void.0, door_void.1))
        .collect();
    assert!(
        !in_door.is_empty(),
        "the door opening is an empty cutout in 3D: nothing is drawn inside it"
    );
    assert!(in_door.iter().all(|(id, _)| *id == r.door.id));
    let in_window: Vec<_> = points
        .iter()
        .filter(|(_, p)| inside(*p, window_void.0, window_void.1))
        .collect();
    assert!(
        !in_window.is_empty(),
        "the window opening is an empty cutout in 3D: nothing is drawn inside it"
    );
    assert!(in_window.iter().all(|(id, _)| *id == r.window));

    // The door leaf stands open on the +Y face at the hinge jamb, a full
    // door width long, as in the plan.
    let leaf = points
        .iter()
        .filter(|(id, p)| *id == r.door.id && p[1] > 100.)
        .map(|(_, p)| *p)
        .collect::<Vec<_>>();
    assert!(!leaf.is_empty(), "the door leaf is drawn open");
    let (lo, hi) = extent(&points, r.door.id);
    assert!(
        (lo[0] - 500.).abs() < 1e-6 && (hi[0] - 1400.).abs() < 1e-6,
        "{lo:?} {hi:?}"
    );
    assert!(
        (lo[1] + 100.).abs() < 1e-6 && (hi[1] - 1000.).abs() < 1e-6,
        "{lo:?} {hi:?}"
    );
    assert!(
        lo[2].abs() < 1e-6 && (hi[2] - 2100.).abs() < 1e-6,
        "{lo:?} {hi:?}"
    );
    let (lo, hi) = extent(&points, r.window);
    assert!(
        (lo[0] - 2200.).abs() < 1e-6 && (hi[0] - 3400.).abs() < 1e-6,
        "{lo:?} {hi:?}"
    );
    assert!(
        (lo[2] - 900.).abs() < 1e-6 && (hi[2] - 2100.).abs() < 1e-6,
        "{lo:?} {hi:?}"
    );
    // The plan leaf line is the 3D leaf's hinge edge: (500, 100) to (500, 1000).
    let plan = scene.architecture_plan().unwrap().as_ref().unwrap();
    let leaf_line = |plan: &viewport::ArchitecturePlan, end: f64| {
        plan.vectors.iter().any(|v| {
            matches!(v.curve, VectorCurve::Line { a, b }
                if (a[0] - 500.).abs() < 1e-6 && (b[0] - 500.).abs() < 1e-6
                    && (a[1] - 100.).abs() < 1e-6 && (b[1] - end).abs() < 1e-6)
        })
    };
    assert!(leaf_line(plan, 1000.));

    // Display elements share their wall's visibility.
    let hidden = scene
        .visible_subset(&std::collections::BTreeSet::new())
        .unwrap();
    assert_eq!(hidden.display_elements().count(), 0);
    assert_eq!(scene.display_elements().count(), 2);
    assert!(scene.display_elements().all(|e| e.host == r.body));

    // The opening cut is unchanged: the wall solid keeps both voids.
    let native = Compiler::default()
        .compile(r.session.document(), r.session.revision())
        .unwrap();
    let volume = native.snapshot().shapes[&r.wall].volume_mm3;
    let expected = 4000. * 200. * 2500. - 900. * 200. * 2100. - 1200. * 200. * 1200.;
    assert!((volume - expected).abs() < 1e-3, "{volume}");
}

#[test]
fn editing_a_door_updates_its_2d_symbol_and_3d_element_together() {
    let mut r = room();
    let mut wider = r.door.clone();
    wider.width = Some(length(1000.));
    let prepared = prepare_wall_opening(
        &r.session,
        r.session.revision(),
        Actor::local(Channel::Ui),
        r.wall,
        &[],
        &WallOpeningChange::ReplaceOpening(wider),
    )
    .unwrap();
    r.session.commit(prepared.prepared).unwrap();
    let scene = publish(&r.session);
    let points = drawn(&scene);
    let (lo, hi) = extent(&points, r.door.id);
    assert!(
        (hi[0] - 1500.).abs() < 1e-6,
        "3D door is 1000 mm wide: {lo:?} {hi:?}"
    );
    assert!(
        (hi[1] - 1100.).abs() < 1e-6,
        "3D leaf is 1000 mm long: {lo:?} {hi:?}"
    );
    let plan = scene.architecture_plan().unwrap().as_ref().unwrap();
    assert!(
        plan.vectors
            .iter()
            .any(|v| matches!(v.curve, VectorCurve::Line { a, b }
            if (a[0] - 500.).abs() < 1e-6 && (b[1] - 1100.).abs() < 1e-6)),
        "2D leaf line is 1000 mm long"
    );
}
