//! CAD-434: rectangular tubes with an outside corner radius. The native
//! member solid uses the stock radius exactly (inside radius = radius less
//! wall), miters keep their geometry, and files keep the radius.
use multicad_core::{
    commands::{Actor, Channel, Prepared, Request, Session},
    document::{Document, Entity, EntityData, EntityId},
    materials::*,
    model::{
        self, AxisAnchor, BodyOutput, ConnectionDefinition, ConnectionKind, MemberDefinition,
        MemberPath, MemberRequest, length,
    },
    project_file::{self, ProjectFile},
};
use multicad_model::{Compiler, Model};

const PI: f64 = std::f64::consts::PI;

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
fn shs(material: EntityId, name: &str, radius: Option<f64>) -> Entity {
    LibraryDefinition::Stock(Stock {
        material,
        specification: StockKind::Profile {
            designation: name.into(),
            length_mm: 6000.,
            width_mm: 40.,
            height_mm: 40.,
            wall_mm: Some(2.),
            shape: ProfileShape::Rectangular,
            section: None,
            corner_radius_mm: radius,
        },
    })
    .entity(EntityId::new(), name)
    .unwrap()
}
/// SHS 40×40×2 with R4 outside corners: inside R2.
fn rounded_area() -> f64 {
    (1600. - (4. - PI) * 16.) - (36. * 36. - (4. - PI) * 4.)
}
fn sharp_area() -> f64 {
    1600. - 36. * 36.
}
fn session() -> (Session, EntityId, EntityId) {
    let mut s = Session::new(Document::new("Rounded tubes").unwrap()).unwrap();
    let steel = Entity::new(
        "S235 steel",
        EntityData::Material {
            density_kg_m3: Some(7850.),
            color_rgba: [120, 124, 130, 255],
            section_hatch: None,
        },
    );
    let rounded = shs(steel.id, "SHS 40×40×2 R4", Some(4.));
    let sharp = shs(steel.id, "SHS 40×40×2", None);
    let ids = (rounded.id, sharp.id);
    s.commit(prepare_library(&s, 0, actor(), &[steel, rounded, sharp]).unwrap())
        .unwrap();
    (s, ids.0, ids.1)
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

#[test]
fn a_rounded_tube_member_has_exactly_the_rounded_section() {
    let (mut s, rounded, sharp) = session();
    let a = MemberRequest::new("Rounded", member(rounded, [0.; 3], [500., 0., 0.]));
    let b = MemberRequest::new("Sharp", member(sharp, [0., 200., 0.], [500., 200., 0.]));
    let (ra, rb) = (a.body, b.body);
    create(&mut s, &[a, b]);
    let m = compile(&s);
    rel(volume(&m, ra), rounded_area() * 500.);
    rel(volume(&m, rb), sharp_area() * 500.);
    // The radius changes the solid: four corners of (4 − π)(R² − r²) each.
    rel(
        volume(&m, rb) - volume(&m, ra),
        (4. - PI) * (16. - 4.) * 500.,
    );
    // The bounding section is still 40 × 40 and the cut list agrees.
    let bounds = m.solid(m.snapshot().bodies[&ra]).unwrap().bounds().unwrap();
    near(bounds.max[1] - bounds.min[1], 40.);
    near(bounds.max[2] - bounds.min[2], 40.);
    let list = m.member_cut_list(s.document()).unwrap();
    let record = list.iter().find(|r| r.name == "Rounded").unwrap();
    assert!(record.issues.is_empty(), "{:?}", record.issues);
    near(record.length_mm, 500.);
    rel(record.volume_mm3, rounded_area() * 500.);

    // Save and reopen: the radius is kept and its capability is declared.
    let mut bytes = vec![];
    let project = ProjectFile::new(s.document().clone()).unwrap();
    assert!(
        project
            .capabilities()
            .contains(CORNER_RADIUS_PROFILE_CAPABILITY)
    );
    project_file::write(&project, &mut bytes, Default::default()).unwrap();
    let opened = project_file::read(
        &mut bytes.as_slice(),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap()
    .project;
    let reopened = Compiler::default().compile(&opened.document, 0).unwrap();
    rel(volume(&reopened, ra), rounded_area() * 500.);
    // An application without the capability refuses the file instead of
    // silently building sharp corners.
    let mut older = project_file::supported_capabilities();
    older.remove(CORNER_RADIUS_PROFILE_CAPABILITY);
    assert!(project_file::read(&mut bytes.as_slice(), Default::default(), &older).is_err());
}

#[test]
fn rounded_tubes_butt_onto_a_flat_face_and_refuse_rust_miters_at_the_connection() {
    let (mut s, rounded, _) = session();
    let steel = s
        .document()
        .entities()
        .find(|e| matches!(e.data, EntityData::Material { .. }))
        .unwrap()
        .id;
    let mut small = shs(steel, "SHS 30×30×2 R3", Some(3.));
    if let EntityData::Extension { definition, .. } = &mut small.data {
        let mut stock: LibraryDefinition = serde_json::from_value(definition.data.clone()).unwrap();
        if let LibraryDefinition::Stock(Stock {
            specification:
                StockKind::Profile {
                    width_mm,
                    height_mm,
                    ..
                },
            ..
        }) = &mut stock
        {
            (*width_mm, *height_mm) = (30., 30.);
        }
        small = stock.entity(small.id, "SHS 30×30×2 R3").unwrap();
    }
    let small_id = small.id;
    s.commit(prepare_library(&s, s.revision(), actor(), &[small]).unwrap())
        .unwrap();
    // A T: the branch ends on the through rail and is butted square onto
    // its flat face, 20 mm above the rail axis.
    let rail = MemberRequest::new("Rail", member(rounded, [0.; 3], [800., 0., 0.]));
    let post = MemberRequest::new("Post", member(small_id, [400., 0., 0.], [400., 500., 0.]));
    let (rail_part, post_part) = (rail.body, post.body);
    create(&mut s, &[rail, post]);
    let (_, prepared) = model::prepare_connection(
        &s,
        s.revision(),
        actor(),
        "Tee",
        &ConnectionDefinition {
            kind: ConnectionKind::Butt,
            members: [post_part, rail_part],
            gap: None,
        },
    )
    .unwrap();
    checked(&mut s, prepared);
    let m = compile(&s);
    let small_area = (900. - (4. - PI) * 9.) - (26. * 26. - (4. - PI) * 1.);
    rel(volume(&m, post_part), small_area * 480.);
    rel(volume(&m, rail_part), rounded_area() * 800.);
    let list = m.member_cut_list(s.document()).unwrap();
    let post = list.iter().find(|r| r.name == "Post").unwrap();
    assert!(post.issues.is_empty(), "{:?}", post.issues);
    near(post.length_mm, 480.);
    assert!(
        post.ends
            .iter()
            .any(|e| e.kind == Some(ConnectionKind::Butt))
    );

    // Mitres of arc-sided sections are not built by the Rust evaluator yet;
    // the refusal names the connection instead of producing wrong geometry.
    let a = MemberRequest::new("A", member(rounded, [0., 1000., 0.], [800., 1000., 0.]));
    let b = MemberRequest::new("B", member(rounded, [800., 1000., 0.], [800., 1500., 0.]));
    let (pa, pb) = (a.body, b.body);
    create(&mut s, &[a, b]);
    let (miter, prepared) = model::prepare_connection(
        &s,
        s.revision(),
        actor(),
        "Corner",
        &ConnectionDefinition {
            kind: ConnectionKind::Miter,
            members: [pa, pb],
            gap: None,
        },
    )
    .unwrap();
    let error = Compiler::default()
        .compile(prepared.document(), s.revision() + 1)
        .err()
        .expect("rounded miters are refused");
    assert_eq!(error.entity, Some(miter));
    assert!(
        error.message.contains("straight-sided"),
        "{}",
        error.message
    );
}
