//! CAD-434 (Q11): parametric channel, T and I/H members. The native solid is
//! the exact section swept along the path, mitres between matching members
//! and square butts as the terminating member are built, and the joints the
//! frame slice does not handle yet are refused with the reason.
use multicad_core::{
    commands::{Actor, Channel, Prepared, Session},
    document::{Document, Entity, EntityData, EntityId},
    materials::*,
    model::{
        self, AxisAnchor, BodyOutput, ConnectionDefinition, ConnectionKind, MemberDefinition,
        MemberPath, MemberRequest, length,
    },
    project_file::{self, ProjectFile},
};
use multicad_model::{Compiler, Model};

fn actor() -> Actor {
    Actor::local(Channel::Ui)
}
fn rel(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-7 * b.abs().max(1.), "{a} != {b}");
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-6, "{a} != {b}");
}
fn compile(s: &Session) -> Model {
    Compiler::default()
        .compile(s.document(), s.revision())
        .unwrap()
}
fn checked(s: &mut Session, prepared: Prepared) {
    Compiler::default()
        .compile(prepared.document(), s.revision() + 1)
        .unwrap();
    s.commit(prepared).unwrap();
}
struct Size {
    shape: ProfileShape,
    w: f64,
    h: f64,
    flange: Option<f64>,
    web: Option<f64>,
}
const CHANNEL: Size = Size {
    shape: ProfileShape::Channel,
    w: 55.,
    h: 100.,
    flange: Some(7.5),
    web: Some(4.5),
};
const TEE: Size = Size {
    shape: ProfileShape::Tee,
    w: 50.,
    h: 50.,
    flange: Some(6.),
    web: Some(6.),
};
const BEAM: Size = Size {
    shape: ProfileShape::Beam,
    w: 100.,
    h: 200.,
    flange: Some(8.5),
    web: Some(5.6),
};
const SHS: Size = Size {
    shape: ProfileShape::Rectangular,
    w: 40.,
    h: 40.,
    flange: Some(2.),
    web: None,
};
impl Size {
    fn area(&self) -> f64 {
        let (w, h) = (self.w, self.h);
        let (tf, tw) = (self.flange.unwrap(), self.web.unwrap_or(0.));
        match self.shape {
            ProfileShape::Channel => w * h - (w - tw) * (h - 2. * tf),
            ProfileShape::Tee => w * tf + tw * (h - tf),
            ProfileShape::Beam => 2. * w * tf + tw * (h - 2. * tf),
            _ => w * h - (w - 2. * tf) * (h - 2. * tf),
        }
    }
    fn stock(&self, material: EntityId, name: &str) -> Entity {
        LibraryDefinition::Stock(Stock {
            material,
            specification: StockKind::Profile {
                designation: name.into(),
                length_mm: 6000.,
                width_mm: self.w,
                height_mm: self.h,
                wall_mm: self.flange,
                shape: self.shape,
                section: None,
                corner_radius_mm: None,
                web_mm: self.web,
            },
        })
        .entity(EntityId::new(), name)
        .unwrap()
    }
}
/// A session with steel stock of every size, by name.
fn session() -> (Session, [EntityId; 4]) {
    let mut s = Session::new(Document::new("Open sections").unwrap()).unwrap();
    let steel = Entity::new(
        "S235 steel",
        EntityData::Material {
            density_kg_m3: Some(7850.),
            color_rgba: [120, 124, 130, 255],
            section_hatch: None,
        },
    );
    let stocks = [
        CHANNEL.stock(steel.id, "UPE 100"),
        TEE.stock(steel.id, "T 50"),
        BEAM.stock(steel.id, "IPE 200"),
        SHS.stock(steel.id, "SHS 40×40×2"),
    ];
    let ids = stocks.each_ref().map(|e| e.id);
    let mut entities = vec![steel];
    entities.extend(stocks);
    s.commit(prepare_library(&s, 0, actor(), &entities).unwrap())
        .unwrap();
    (s, ids)
}
fn member(stock: EntityId, a: [f64; 3], b: [f64; 3]) -> MemberDefinition {
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
    }
}
fn create(s: &mut Session, requests: &[MemberRequest]) {
    let p = model::prepare_members(s, s.revision(), actor(), requests, None).unwrap();
    checked(s, p.prepared);
}
fn volume(m: &Model, part: EntityId) -> f64 {
    m.solid(m.snapshot().bodies[&part])
        .unwrap()
        .volume_mm3()
        .unwrap()
}
fn connect(s: &Session, kind: ConnectionKind, members: [EntityId; 2]) -> Result<Prepared, String> {
    let (_, prepared) = model::prepare_connection(
        s,
        s.revision(),
        actor(),
        "Joint",
        &ConnectionDefinition {
            kind,
            members,
            gap: None,
        },
    )
    .map_err(|e| e.to_string())?;
    Ok(prepared)
}
/// The reason a joint is refused, whether at preparation or at evaluation.
fn refusal(s: &Session, kind: ConnectionKind, members: [EntityId; 2]) -> String {
    match connect(s, kind, members) {
        Err(message) => message,
        Ok(prepared) => {
            Compiler::default()
                .compile(prepared.document(), s.revision() + 1)
                .err()
                .expect("the joint is refused")
                .message
        }
    }
}

#[test]
fn channel_tee_and_beam_members_are_their_exact_sections() {
    let (mut s, [channel, tee, beam, _]) = session();
    let requests = [
        MemberRequest::new("Channel", member(channel, [0.; 3], [700., 0., 0.])),
        MemberRequest::new("Tee", member(tee, [0., 300., 0.], [700., 300., 0.])),
        MemberRequest::new("Beam", member(beam, [0., 600., 0.], [700., 600., 0.])),
    ];
    let parts = requests.each_ref().map(|r| r.body);
    create(&mut s, &requests);
    let m = compile(&s);
    let list = m.member_cut_list(s.document()).unwrap();
    for (part, size, name) in [
        (parts[0], &CHANNEL, "Channel"),
        (parts[1], &TEE, "Tee"),
        (parts[2], &BEAM, "Beam"),
    ] {
        rel(volume(&m, part), size.area() * 700.);
        let bounds = m
            .solid(m.snapshot().bodies[&part])
            .unwrap()
            .bounds()
            .unwrap();
        near(bounds.max[0] - bounds.min[0], 700.);
        let mut across = [bounds.max[1] - bounds.min[1], bounds.max[2] - bounds.min[2]];
        across.sort_by(f64::total_cmp);
        let mut expected = [size.w, size.h];
        expected.sort_by(f64::total_cmp);
        near(across[0], expected[0]);
        near(across[1], expected[1]);
        let record = list.iter().find(|r| r.name == name).unwrap();
        assert!(record.issues.is_empty(), "{:?}", record.issues);
        near(record.length_mm, 700.);
        rel(record.volume_mm3, size.area() * 700.);
    }

    // Save and reopen; an application without the capability refuses the
    // file rather than guessing the section.
    let project = ProjectFile::new(s.document().clone()).unwrap();
    assert!(
        project
            .capabilities()
            .contains(OPEN_SECTION_PROFILE_CAPABILITY)
    );
    let mut bytes = vec![];
    project_file::write(&project, &mut bytes, Default::default()).unwrap();
    let opened = project_file::read(
        &mut bytes.as_slice(),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap()
    .project;
    let reopened = Compiler::default().compile(&opened.document, 0).unwrap();
    rel(volume(&reopened, parts[2]), BEAM.area() * 700.);
    let mut older = project_file::supported_capabilities();
    older.remove(OPEN_SECTION_PROFILE_CAPABILITY);
    assert!(project_file::read(&mut bytes.as_slice(), Default::default(), &older).is_err());
}

#[test]
fn beams_mitre_and_tees_butt_as_the_terminating_member() {
    let (mut s, [_, tee, beam, shs]) = session();
    // A square corner of two I/H beams: the doubly symmetric section keeps
    // its centroid on the mitre pivot, so each member keeps area × length.
    let a = MemberRequest::new("A", member(beam, [0.; 3], [1000., 0., 0.]));
    let b = MemberRequest::new("B", member(beam, [1000., 0., 0.], [1000., 800., 0.]));
    let (pa, pb) = (a.body, b.body);
    create(&mut s, &[a, b]);
    let prepared = connect(&s, ConnectionKind::Miter, [pa, pb]).unwrap();
    checked(&mut s, prepared);
    let m = compile(&s);
    rel(volume(&m, pa), BEAM.area() * 1000.);
    rel(volume(&m, pb), BEAM.area() * 800.);
    let list = m.member_cut_list(s.document()).unwrap();
    for name in ["A", "B"] {
        let record = list.iter().find(|r| r.name == name).unwrap();
        assert!(record.issues.is_empty(), "{:?}", record.issues);
        assert!(
            record
                .ends
                .iter()
                .any(|e| e.kind == Some(ConnectionKind::Miter))
        );
    }

    // A T section post butted square onto a tube rail, 20 mm above its axis.
    let rail = MemberRequest::new("Rail", member(shs, [0., 2000., 0.], [800., 2000., 0.]));
    let post = MemberRequest::new("Post", member(tee, [400., 2000., 0.], [400., 2500., 0.]));
    let (rail_part, post_part) = (rail.body, post.body);
    create(&mut s, &[rail, post]);
    let prepared = connect(&s, ConnectionKind::Butt, [post_part, rail_part]).unwrap();
    checked(&mut s, prepared);
    let m = compile(&s);
    rel(volume(&m, post_part), TEE.area() * 480.);
    let list = m.member_cut_list(s.document()).unwrap();
    let post = list.iter().find(|r| r.name == "Post").unwrap();
    assert!(post.issues.is_empty(), "{:?}", post.issues);
    near(post.length_mm, 480.);
}

#[test]
fn joints_frames_cannot_build_on_open_sections_yet_are_refused_with_the_reason() {
    let (mut s, [channel, tee, beam, shs]) = session();
    let rail = MemberRequest::new("Beam rail", member(beam, [0.; 3], [800., 0., 0.]));
    let post = MemberRequest::new("Post", member(shs, [400., 0., 0.], [400., 600., 0.]));
    let channel_rail = MemberRequest::new(
        "Channel rail",
        member(channel, [0., 2000., 0.], [800., 2000., 0.]),
    );
    let tube_post = MemberRequest::new(
        "Tube post",
        member(shs, [400., 2000., 0.], [400., 2600., 0.]),
    );
    let tube_rail =
        MemberRequest::new("Tube rail", member(shs, [0., 4000., 0.], [800., 4000., 0.]));
    let tee_post = MemberRequest::new(
        "Tee post",
        member(tee, [400., 4000., 0.], [400., 4600., 0.]),
    );
    let ids = [
        &rail,
        &post,
        &channel_rail,
        &tube_post,
        &tube_rail,
        &tee_post,
    ]
    .map(|r| r.body);
    create(
        &mut s,
        &[rail, post, channel_rail, tube_post, tube_rail, tee_post],
    );
    let before = s.document().clone();
    for (kind, members, needle) in [
        (
            ConnectionKind::Butt,
            [ids[1], ids[0]],
            "Butt joints with an I/H beam running through are not supported yet",
        ),
        (
            ConnectionKind::Cope,
            [ids[3], ids[2]],
            "Copes against a channel are not supported yet; use a butt joint onto a flat face",
        ),
        (
            ConnectionKind::Cope,
            [ids[5], ids[4]],
            "Coping a T section is not supported yet",
        ),
    ] {
        let message = refusal(&s, kind, members);
        assert!(message.contains(needle), "{kind:?}: {message}");
    }
    assert_eq!(s.document(), &before, "refused joints change nothing");
}
