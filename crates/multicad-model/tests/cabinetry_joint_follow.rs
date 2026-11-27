//! CAD-431: a finger joint follows its panels. Two perpendicular panels meet
//! in a T (the upright stands on the centre of the base); after either panel
//! is moved or rotated, both joint cuts match a joint created fresh at the new
//! location. A move that separates the panels is refused with the joint named,
//! never left with stale cuts. Undo, Redo and save/reopen stay consistent.
use multicad_core::{
    commands::{Actor, Channel, Prepared, Session},
    document::{Document, Entity, EntityData, EntityId},
    materials::*,
    model::{
        self, BodyInput, BodyOutput, GrainAxis, NumericTransform, PanelDefinition, PanelEdges,
        PanelJointDefinition, PanelJointKind, PanelJointPair, PanelPlane, PanelRequest,
        PanelSource, TransformRequest, angle, length,
    },
    project_file::{self, ProjectFile},
};
use multicad_model::{Compiler, EvaluationCache, Model};

fn actor() -> Actor {
    Actor::local(Channel::Ui)
}
fn compile(document: &Document) -> Model {
    Compiler::default().compile(document, 0).unwrap()
}
fn checked(s: &mut Session, prepared: Prepared) {
    Compiler::default()
        .compile(prepared.document(), s.revision() + 1)
        .unwrap();
    s.commit(prepared).unwrap();
}
fn panel(
    stock: EntityId,
    name: &str,
    plane: PanelPlane,
    origin: [f64; 3],
    l: f64,
    w: f64,
) -> PanelRequest {
    PanelRequest::new(
        name,
        PanelDefinition {
            output: BodyOutput {
                body: EntityId::new(),
            },
            source: PanelSource::Rectangle {
                stock,
                plane,
                origin: origin.map(length),
                length: length(l),
                width: length(w),
            },
            grain: GrainAxis::U,
            edges: PanelEdges::default(),
            oversize: None,
        },
    )
}

struct Fixture {
    s: Session,
    base: EntityId,
    upright: EntityId,
    joint: EntityId,
}
fn fingers() -> PanelJointKind {
    PanelJointKind::Finger {
        count: 5,
        fit: length(0.),
    }
}
/// A 600 × 400 base on XY and an 18 mm upright on YZ standing on its centre
/// line, with a finger joint between them. Panel origins are the arguments.
fn fixture(base_origin: [f64; 3], upright_origin: [f64; 3]) -> Fixture {
    let mut s = Session::new(Document::new("Joint follow").unwrap()).unwrap();
    let birch = Entity::new(
        "Birch plywood",
        EntityData::Material {
            density_kg_m3: None,
            color_rgba: [214, 186, 140, 255],
            section_hatch: None,
        },
    );
    let stock = LibraryDefinition::Stock(Stock {
        material: birch.id,
        specification: StockKind::Panel {
            thickness_mm: 18.,
            length_mm: 2500.,
            width_mm: 1250.,
            grain: Grain::AlongLength,
        },
    })
    .entity(EntityId::new(), "Birch ply 18")
    .unwrap();
    let id = stock.id;
    s.commit(prepare_library(&s, 0, actor(), &[birch, stock]).unwrap())
        .unwrap();
    let base = panel(id, "Panel 1", PanelPlane::Xy, base_origin, 600., 400.);
    let upright = panel(id, "Panel 2", PanelPlane::Yz, upright_origin, 400., 300.);
    for p in [&base, &upright] {
        let prepared = model::prepare_panel(&s, s.revision(), actor(), p).unwrap();
        checked(&mut s, prepared.prepared);
    }
    let def = PanelJointDefinition {
        generated_by: None,
        kind: fingers(),
        pairs: vec![PanelJointPair {
            a: base.body,
            b: upright.body,
        }],
    };
    let p = model::prepare_panel_joint(&s, s.revision(), actor(), "Finger joint", &def).unwrap();
    checked(&mut s, p.model.prepared);
    Fixture {
        s,
        base: base.body,
        upright: upright.body,
        joint: p.joint,
    }
}
fn last_feature(document: &Document, part: EntityId) -> EntityId {
    let Some(EntityData::Part { features, .. }) = document.get(part).map(|e| &e.data) else {
        panic!("not a Part")
    };
    *features.last().unwrap()
}
fn transform(s: &Session, part: EntityId, parameters: NumericTransform) -> Prepared {
    let request = TransformRequest::new(
        BodyInput {
            body: part,
            producer: last_feature(s.document(), part),
        },
        parameters,
    );
    model::prepare_transform(s, s.revision(), actor(), &request)
        .unwrap()
        .prepared
}
fn shift(by: [f64; 3]) -> NumericTransform {
    NumericTransform {
        translation: by.map(length),
        ..Default::default()
    }
}

/// The faces of a Part's solid, rounded and sorted, so two solids built by
/// different feature chains compare by their geometry.
fn faces(m: &Model, part: EntityId) -> Vec<[i64; 4]> {
    let solid = m.solid(m.snapshot().bodies[&part]).unwrap();
    let mut faces: Vec<_> = solid
        .faces()
        .unwrap()
        .into_iter()
        .map(|f| {
            let r = |x: f64| (x * 1000.).round() as i64;
            [
                r(f.center_mm[0]),
                r(f.center_mm[1]),
                r(f.center_mm[2]),
                r(f.area_mm2),
            ]
        })
        .collect();
    faces.sort();
    faces
}
fn volume(m: &Model, part: EntityId) -> f64 {
    m.solid(m.snapshot().bodies[&part])
        .unwrap()
        .volume_mm3()
        .unwrap()
}
fn same_geometry(moved: &Fixture, fresh: &Fixture) {
    let (a, b) = (compile(moved.s.document()), compile(fresh.s.document()));
    for (m, f) in [(moved.base, fresh.base), (moved.upright, fresh.upright)] {
        assert!((volume(&a, m) - volume(&b, f)).abs() < 1e-3);
        assert_eq!(faces(&a, m), faces(&b, f));
    }
}

const BASE: [f64; 3] = [0., 0., 0.];
const UPRIGHT: [f64; 3] = [291., 0., 0.];

#[test]
fn moving_the_base_re_evaluates_the_finger_joint_at_the_new_location() {
    let mut f = fixture(BASE, UPRIGHT);
    let before = compile(f.s.document());
    // Slide the base 100 mm along X: the upright still stands on it, now
    // 191 mm from the base's start instead of 291 mm.
    let p = transform(&f.s, f.base, shift([100., 0., 0.]));
    checked(&mut f.s, p);
    same_geometry(&f, &fixture([100., 0., 0.], UPRIGHT));
    // The upright's slots stay where the base's fingers now are: nothing of
    // the upright changes for a move along the base's length.
    assert_eq!(
        faces(&compile(f.s.document()), f.upright),
        faces(&before, f.upright)
    );
}

#[test]
fn moving_the_upright_re_evaluates_the_finger_joint_at_the_new_location() {
    let mut f = fixture(BASE, UPRIGHT);
    let p = transform(&f.s, f.upright, shift([-150., 0., 0.]));
    checked(&mut f.s, p);
    same_geometry(&f, &fixture(BASE, [141., 0., 0.]));
    // A second move of the same panel composes with the first.
    let p = transform(&f.s, f.upright, shift([50., 0., 0.]));
    checked(&mut f.s, p);
    same_geometry(&f, &fixture(BASE, [191., 0., 0.]));
}

#[test]
fn rotating_both_panels_together_keeps_the_joint_on_them() {
    let mut f = fixture(BASE, UPRIGHT);
    // A quarter turn of both panels about Z through the world origin.
    let turn = NumericTransform {
        axis: [0., 0., 1.],
        angle: angle(90.),
        ..Default::default()
    };
    // Both in one transaction: turning one first would separate them.
    let requests: Vec<_> = [f.base, f.upright]
        .map(|part| {
            TransformRequest::new(
                BodyInput {
                    body: part,
                    producer: last_feature(f.s.document(), part),
                },
                turn.clone(),
            )
        })
        .into();
    let p = model::prepare_transforms(&f.s, f.s.revision(), actor(), &requests).unwrap();
    checked(&mut f.s, p.prepared);
    let rotated = compile(f.s.document());
    // Same solids, turned: each face centre is the original one turned 90°.
    let turned = |m: &Model, part| {
        let mut v: Vec<_> = faces(m, part)
            .into_iter()
            .map(|[x, y, z, a]| [-y, x, z, a])
            .collect();
        v.sort();
        v
    };
    let g = fixture(BASE, UPRIGHT);
    let gm = compile(g.s.document());
    assert_eq!(faces(&rotated, f.base), turned(&gm, g.base));
    assert_eq!(faces(&rotated, f.upright), turned(&gm, g.upright));
}

#[test]
fn rotating_one_panel_out_of_square_or_apart_is_refused_with_the_joint_named() {
    let f = fixture(BASE, UPRIGHT);
    // Lift the upright clear of the base: the panels no longer touch.
    let lifted = transform(&f.s, f.upright, shift([0., 0., 100.]));
    let Err(e) = Compiler::default().compile(lifted.document(), 0) else {
        panic!("a joint between separated panels must not compile")
    };
    let text = e.to_string();
    assert!(text.contains("Finger joint"), "{text}");
    assert!(text.contains("overlap") || text.contains("touch"), "{text}");
    // The fix names the upright's own origin field, not where it now is:
    // its panel origin Z 0 plus the 100 mm move must reach the base at 0.
    assert!(text.contains("Z 0 → -100 mm"), "{text}");
    // Slide the upright off the end of the base.
    let off = transform(&f.s, f.upright, shift([1000., 0., 0.]));
    assert!(Compiler::default().compile(off.document(), 0).is_err());
    // A 30° turn leaves the panels out of square.
    let skew = transform(
        &f.s,
        f.upright,
        NumericTransform {
            axis: [0., 0., 1.],
            angle: angle(30.),
            pivot: [300., 200., 0.].map(length),
            ..Default::default()
        },
    );
    let Err(e) = Compiler::default().compile(skew.document(), 0) else {
        panic!("an out-of-square joint must not compile")
    };
    assert!(e.to_string().contains("Finger joint"), "{e}");
}

#[test]
fn undo_redo_and_save_reopen_keep_the_moved_joint() {
    let mut f = fixture(BASE, UPRIGHT);
    let original = f.s.document().clone();
    let original_model = compile(&original);
    let p = transform(&f.s, f.base, shift([100., 0., 0.]));
    checked(&mut f.s, p);
    let moved = f.s.document().clone();
    let moved_faces = faces(&compile(&moved), f.base);
    // One Undo reverts the move and the joint cuts with it.
    f.s.undo(f.s.revision(), actor()).unwrap();
    assert_eq!(f.s.document(), &original);
    let undone = compile(f.s.document());
    for part in [f.base, f.upright] {
        assert_eq!(faces(&undone, part), faces(&original_model, part));
    }
    f.s.redo(f.s.revision(), actor()).unwrap();
    assert_eq!(f.s.document(), &moved);
    assert_eq!(faces(&compile(f.s.document()), f.base), moved_faces);
    // The joint still names the same panels; save and reopen rebuild the
    // same moved geometry.
    let joint = model::read_panel_joint(f.s.document(), f.joint).unwrap();
    assert_eq!(joint.pairs[0].a, f.base);
    let mut bytes = vec![];
    project_file::write(
        &ProjectFile::new(f.s.document().clone()).unwrap(),
        &mut bytes,
        Default::default(),
    )
    .unwrap();
    let opened = project_file::read(
        &mut bytes.as_slice(),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap()
    .project;
    assert_eq!(&opened.document, f.s.document());
    let reopened = compile(&opened.document);
    for part in [f.base, f.upright] {
        assert_eq!(faces(&reopened, part), faces(&compile(&moved), part));
    }
}

#[test]
fn a_cached_compile_rebuilds_the_partner_cut_after_a_move() {
    // The evaluation cache must not reuse the upright's cut when only the
    // base's transform changed.
    let mut f = fixture(BASE, UPRIGHT);
    let compiler = Compiler::default();
    let mut cache = EvaluationCache::default();
    compiler
        .compile_cached(&mut cache, f.s.document(), f.s.revision(), None)
        .unwrap();
    let p = transform(&f.s, f.base, shift([-50., 0., 0.]));
    checked(&mut f.s, p);
    let cached = compiler
        .compile_cached(&mut cache, f.s.document(), f.s.revision(), None)
        .unwrap();
    // The upright's own panel feature is reused; its joint cut is rebuilt.
    let stats = cache.last_compile();
    let upright = model::panel_feature(f.s.document(), f.upright).unwrap();
    let cut = last_feature(f.s.document(), f.upright);
    assert!(stats.reused.contains(&upright), "{stats:?}");
    assert!(stats.evaluated.contains(&cut), "{stats:?}");
    let fresh = fixture([-50., 0., 0.], UPRIGHT);
    let fm = compile(fresh.s.document());
    assert_eq!(faces(&cached, f.upright), faces(&fm, fresh.upright));
    assert_eq!(faces(&cached, f.base), faces(&fm, fresh.base));
}
