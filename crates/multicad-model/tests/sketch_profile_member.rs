//! CAD-385: profile stock whose section is an owned constraint sketch
//! produces real native member solids, follows edits of that sketch through
//! the normal transaction/rebuild flow, and survives Undo/Redo and
//! save/reopen. Invalid sections never reach the document.
use multicad_core::{
    commands::{Actor, Channel, Prepared, Session},
    constraint_sketch::{Curve, CurveKind, Point, Sketch},
    document::{Document, Entity, EntityData, EntityId},
    materials::*,
    model::{
        self, AxisAnchor, BodyOutput, ConnectionDefinition, ConnectionKind, MemberDefinition,
        MemberPath, MemberRequest, length,
    },
    project_file::{self, ProjectFile},
};
use multicad_model::{Compiler, MemberRecord, Model};
use std::f64::consts::PI;

fn actor() -> Actor {
    Actor::local(Channel::Ui)
}
fn rel(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-7 * b.abs().max(1.), "{a} != {b}");
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-6, "{a} != {b}");
}
fn checked(s: &mut Session, prepared: Prepared) {
    Compiler::default()
        .compile(prepared.document(), s.revision() + 1)
        .unwrap();
    s.commit(prepared).unwrap();
}
fn compile(s: &Session) -> Model {
    Compiler::default()
        .compile(s.document(), s.revision())
        .unwrap()
}
fn record<'a>(list: &'a [MemberRecord], name: &str) -> &'a MemberRecord {
    list.iter().find(|r| r.name == name).unwrap()
}

#[derive(Default)]
struct Draft(Sketch);
impl Draft {
    fn point(&mut self, p: [f64; 2]) -> EntityId {
        if let Some(q) = self.0.points.iter().find(|q| q.position_mm == p) {
            return q.id;
        }
        let id = EntityId::new();
        self.0.points.push(Point { id, position_mm: p });
        id
    }
    fn curve(&mut self, kind: CurveKind) {
        self.0.curves.push(Curve {
            id: EntityId::new(),
            construction: false,
            kind,
        });
    }
    fn line(&mut self, a: [f64; 2], b: [f64; 2]) {
        let (start, end) = (self.point(a), self.point(b));
        self.curve(CurveKind::Line { start, end });
    }
    fn arc(&mut self, c: [f64; 2], a: [f64; 2], b: [f64; 2]) {
        let (center, start, end) = (self.point(c), self.point(a), self.point(b));
        self.curve(CurveKind::Arc {
            center,
            start,
            end,
            clockwise: false,
        });
    }
    fn circle(&mut self, c: [f64; 2], r: f64) {
        let center = self.point(c);
        self.curve(CurveKind::Circle {
            center,
            radius_mm: r,
        });
    }
}

/// A 40 × 20 mm slot (20 mm straight sides, R10 ends) with a bore of
/// radius `bore`, drawn away from the sketch origin.
fn slot(bore: f64) -> Sketch {
    let mut d = Draft::default();
    let (x, y) = (100., 50.);
    d.line([x, y], [x + 20., y]);
    d.arc([x + 20., y + 10.], [x + 20., y], [x + 20., y + 20.]);
    d.line([x + 20., y + 20.], [x, y + 20.]);
    d.arc([x, y + 10.], [x, y + 20.], [x, y]);
    d.circle([x + 10., y + 10.], bore);
    d.0
}
fn slot_area(bore: f64) -> f64 {
    20. * 20. + PI * 100. - PI * bore * bore
}
/// An open T section: a 40 × 5 flange and a 5 × 25 web.
fn tee() -> Sketch {
    let mut d = Draft::default();
    let p = [
        [0., 0.],
        [40., 0.],
        [40., 5.],
        [22.5, 5.],
        [22.5, 30.],
        [17.5, 30.],
        [17.5, 5.],
        [0., 5.],
    ];
    for i in 0..p.len() {
        d.line(p[i], p[(i + 1) % p.len()]);
    }
    d.0
}
fn sketch_stock(
    id: EntityId,
    material: EntityId,
    name: &str,
    sketch: Sketch,
) -> Result<Entity, String> {
    let size = model::resolve_sketch_section(&sketch)
        .map_err(|e| e.to_string())?
        .size_mm;
    LibraryDefinition::Stock(Stock {
        material,
        specification: StockKind::Profile {
            designation: name.into(),
            length_mm: 6000.,
            width_mm: size[0],
            height_mm: size[1],
            wall_mm: None,
            shape: ProfileShape::Sketch,
            section: Some(Box::new(ProfileSketch { sketch })),
        },
    })
    .entity(id, name)
    .map_err(|e| e.to_string())
}
fn member(name: &str, stock: EntityId, a: [f64; 3], b: [f64; 3]) -> MemberRequest {
    MemberRequest::new(
        name,
        MemberDefinition {
            output: BodyOutput {
                body: EntityId::new(),
            },
            stock,
            path: MemberPath::Points {
                start: a.map(length),
                end: b.map(length),
            },
            anchor: [AxisAnchor::Center; 2],
            rotation: None,
            extend_start: None,
            extend_end: None,
        },
    )
}

struct Fixture {
    s: Session,
    steel: EntityId,
    slot: EntityId,
    tee: EntityId,
    requests: Vec<MemberRequest>,
}
fn fixture() -> Fixture {
    let mut s = Session::new(Document::new("Sketch profile").unwrap()).unwrap();
    let steel = Entity::new(
        "S235 steel",
        EntityData::Material {
            density_kg_m3: Some(7850.),
            color_rgba: [120, 124, 130, 255],
            section_hatch: None,
        },
    );
    let slot = sketch_stock(EntityId::new(), steel.id, "Slot tube", slot(5.)).unwrap();
    let tee = sketch_stock(EntityId::new(), steel.id, "Custom T", tee()).unwrap();
    let (slot_id, tee_id) = (slot.id, tee.id);
    let p = prepare_library(&s, 0, actor(), &[steel.clone(), slot, tee]).unwrap();
    s.commit(p).unwrap();
    let requests = vec![
        member("Rail", slot_id, [0., 0., 0.], [500., 0., 0.]),
        member("Post", tee_id, [0., 300., 0.], [0., 300., 400.]),
    ];
    let p = model::prepare_members(&s, s.revision(), actor(), &requests, None).unwrap();
    checked(&mut s, p.prepared);
    Fixture {
        s,
        steel: steel.id,
        slot: slot_id,
        tee: tee_id,
        requests,
    }
}

#[test]
fn hollow_and_open_sketch_sections_produce_real_member_solids() {
    let f = fixture();
    let m = compile(&f.s);
    let list = m.member_cut_list(f.s.document()).unwrap();
    let rail = record(&list, "Rail");
    assert_eq!(rail.shape, ProfileShape::Sketch);
    rel(rail.volume_mm3, slot_area(5.) * 500.);
    near(rail.length_mm, 500.);
    assert_eq!(rail.section_mm, [40., 20.]);
    near(rail.measured_section_mm[0], 40.);
    near(rail.measured_section_mm[1], 20.);
    assert!(rail.issues.is_empty(), "{:?}", rail.issues);
    let post = record(&list, "Post");
    rel(post.volume_mm3, (40. * 5. + 5. * 25.) * 400.);
    assert_eq!(post.section_mm, [40., 30.]);
    assert!(post.issues.is_empty(), "{:?}", post.issues);
    // The bore is real material removal: a solid slot weighs more.
    let solid = m.solid(m.snapshot().bodies[&f.requests[0].body]).unwrap();
    rel(solid.volume_mm3().unwrap(), slot_area(5.) * 500.);
}

#[test]
fn editing_the_section_sketch_rebuilds_members_with_undo_redo_and_reopen() {
    let mut f = fixture();
    let before = f.s.document().clone();
    let edited = sketch_stock(f.slot, f.steel, "Slot tube", slot(7.)).unwrap();
    let p = prepare_library(&f.s, f.s.revision(), actor(), &[edited]).unwrap();
    checked(&mut f.s, p);
    let list = compile(&f.s).member_cut_list(f.s.document()).unwrap();
    let rail = record(&list, "Rail");
    assert_eq!(
        rail.part, f.requests[0].body,
        "the member keeps its identity"
    );
    rel(rail.volume_mm3, slot_area(7.) * 500.);
    let after = f.s.document().clone();
    f.s.undo(f.s.revision(), actor()).unwrap();
    assert_eq!(f.s.document(), &before);
    let list = compile(&f.s).member_cut_list(f.s.document()).unwrap();
    rel(record(&list, "Rail").volume_mm3, slot_area(5.) * 500.);
    f.s.redo(f.s.revision(), actor()).unwrap();
    assert_eq!(f.s.document(), &after);
    // Save and reopen: the owned sketch, its capability and the geometry.
    let project = ProjectFile::new(after.clone()).unwrap();
    assert!(project.capabilities().contains(SKETCH_PROFILE_CAPABILITY));
    let mut bytes = vec![];
    project_file::write(&project, &mut bytes, Default::default()).unwrap();
    let reopened = project_file::read(
        &mut bytes.as_slice(),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap()
    .project
    .document;
    assert_eq!(reopened, after);
    let s = Session::new(reopened).unwrap();
    let list = compile(&s).member_cut_list(s.document()).unwrap();
    rel(record(&list, "Rail").volume_mm3, slot_area(7.) * 500.);
    // A reader without the capability refuses the file instead of losing it.
    let mut older = project_file::supported_capabilities();
    older.remove(SKETCH_PROFILE_CAPABILITY);
    assert!(project_file::read(&mut bytes.as_slice(), Default::default(), &older).is_err());
}

#[test]
fn invalid_section_edits_never_reach_the_document_or_geometry() {
    let f = fixture();
    let before = f.s.document().clone();
    let revision = f.s.revision();
    let mut open = slot(5.);
    open.curves.remove(1);
    let error = sketch_stock(f.slot, f.steel, "Slot tube", open.clone()).unwrap_err();
    assert!(error.contains("open"), "{error}");
    // Even a hand-built definition with a stale size is refused at the
    // library boundary.
    let stale = LibraryDefinition::Stock(Stock {
        material: f.steel,
        specification: StockKind::Profile {
            designation: "Slot tube".into(),
            length_mm: 6000.,
            width_mm: 40.,
            height_mm: 20.,
            wall_mm: None,
            shape: ProfileShape::Sketch,
            section: Some(Box::new(ProfileSketch { sketch: open })),
        },
    });
    assert!(stale.entity(f.slot, "Slot tube").is_err());
    let mut bowtie = Draft::default();
    bowtie.line([0., 0.], [10., 10.]);
    bowtie.line([10., 10.], [10., 0.]);
    bowtie.line([10., 0.], [0., 10.]);
    bowtie.line([0., 10.], [0., 0.]);
    let error = sketch_stock(f.slot, f.steel, "Slot tube", bowtie.0).unwrap_err();
    assert!(error.contains("crosses"), "{error}");
    assert_eq!(f.s.document(), &before);
    assert_eq!(f.s.revision(), revision);
}

fn connect(f: &mut Fixture, name: &str, kind: ConnectionKind, members: [EntityId; 2]) {
    let (_, p) = model::prepare_connection(
        &f.s,
        f.s.revision(),
        actor(),
        name,
        &ConnectionDefinition {
            kind,
            members,
            gap: None,
        },
    )
    .unwrap();
    checked(&mut f.s, p);
}

#[test]
fn custom_sections_miter_and_butt_but_refuse_copes() {
    let mut f = fixture();
    let shs = LibraryDefinition::Stock(Stock {
        material: f.steel,
        specification: StockKind::Profile {
            designation: "SHS 40×40×2".into(),
            length_mm: 6000.,
            width_mm: 40.,
            height_mm: 40.,
            wall_mm: Some(2.),
            shape: ProfileShape::Rectangular,
            section: None,
        },
    })
    .entity(EntityId::new(), "SHS 40×40×2")
    .unwrap();
    let shs_id = shs.id;
    let p = prepare_library(&f.s, f.s.revision(), actor(), &[shs]).unwrap();
    f.s.commit(p).unwrap();
    let requests = [
        member("Beam", shs_id, [0., 1000., 0.], [500., 1000., 0.]),
        member("Branch", f.slot, [250., 1000., 0.], [250., 1000., 300.]),
        member("Corner A", f.tee, [0., 2000., 0.], [500., 2000., 0.]),
        member("Corner B", f.tee, [500., 2000., 0.], [500., 2500., 0.]),
        member("Slot A", f.slot, [0., 3000., 0.], [500., 3000., 0.]),
        member("Slot B", f.slot, [500., 3000., 0.], [500., 3500., 0.]),
        member("Stub", f.slot, [250., 0., 0.], [250., 0., 300.]),
    ];
    let p = model::prepare_members(&f.s, f.s.revision(), actor(), &requests, None).unwrap();
    checked(&mut f.s, p.prepared);
    // A cope follows the receiving outside; a custom receiver is refused.
    let (_, cope) = model::prepare_connection(
        &f.s,
        f.s.revision(),
        actor(),
        "Cope",
        &ConnectionDefinition {
            kind: ConnectionKind::Cope,
            members: [requests[6].body, f.requests[0].body],
            gap: None,
        },
    )
    .unwrap();
    let error = Compiler::default()
        .compile(cope.document(), f.s.revision() + 1)
        .err()
        .expect("cope refused")
        .message;
    assert!(error.contains("custom sketch section"), "{error}");
    // The custom branch butts onto the tube's flat 40 mm face.
    connect(
        &mut f,
        "Butt",
        ConnectionKind::Butt,
        [requests[1].body, requests[0].body],
    );
    // Two identical symmetric straight-sided custom members miter at a corner.
    connect(
        &mut f,
        "Miter",
        ConnectionKind::Miter,
        [requests[2].body, requests[3].body],
    );
    // Arc-sided miters are a native evaluator limit, refused explicitly.
    let (_, arc_miter) = model::prepare_connection(
        &f.s,
        f.s.revision(),
        actor(),
        "Arc miter",
        &ConnectionDefinition {
            kind: ConnectionKind::Miter,
            members: [requests[4].body, requests[5].body],
            gap: None,
        },
    )
    .unwrap();
    let error = Compiler::default()
        .compile(arc_miter.document(), f.s.revision() + 1)
        .err()
        .expect("arc-sided miter refused")
        .message;
    assert!(error.contains("straight-sided"), "{error}");
    let list = compile(&f.s).member_cut_list(f.s.document()).unwrap();
    let branch = record(&list, "Branch");
    near(branch.span_mm[0], 20.);
    rel(branch.volume_mm3, slot_area(5.) * 280.);
    assert_eq!(branch.ends[0].kind, Some(ConnectionKind::Butt));
    for name in ["Corner A", "Corner B"] {
        let corner = record(&list, name);
        rel(corner.volume_mm3, (40. * 5. + 5. * 25.) * 500.);
        assert!(
            corner
                .ends
                .iter()
                .any(|e| e.kind == Some(ConnectionKind::Miter))
        );
    }
}
