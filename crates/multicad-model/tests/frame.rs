//! CAD-439: frames drawn as lines. A `weldment.frame` holds nodes, members
//! and joint overrides; each member is a `weldment.member` feature with its
//! own Part; joints are derived at every node and built natively. These
//! tests measure the native solids and the cut list.
use multicad_core::{
    commands::{Actor, Channel, Prepared, Request, Session},
    constraint_sketch::Dimension,
    document::{Document, Entity, EntityData, EntityId, Mutation},
    materials::*,
    model::{
        self, ChainSegment, ConnectionKind, FrameDefinition, FrameMemberRequest, FramePlane,
        FrameRequest, FrameSide, JointOverride, NodePosition, angle, length,
    },
    parameters::BoundExpression,
    project_file::{self, ProjectFile},
};
use multicad_model::{Compiler, MemberRecord, Model};

fn actor() -> Actor {
    Actor::local(Channel::Ui)
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-6, "{a} != {b}");
}
fn rel(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-7 * b.abs().max(1.), "{a} != {b}");
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
/// The evaluator's refusal of a prepared change.
fn refused(s: &Session, prepared: &Prepared) -> String {
    Compiler::default()
        .compile(prepared.document(), s.revision() + 1)
        .err()
        .expect("the change is refused")
        .message
}
fn record<'a>(list: &'a [MemberRecord], name: &str) -> &'a MemberRecord {
    list.iter()
        .find(|r| r.name == name)
        .unwrap_or_else(|| panic!("no member {name}"))
}

struct Stocks {
    shs: EntityId,
    small: EntityId,
    chs: EntityId,
    big_chs: EntityId,
    channel: EntityId,
}
fn stock(material: EntityId, name: &str, shape: ProfileShape, w: f64, h: f64, wall: f64) -> Entity {
    LibraryDefinition::Stock(Stock {
        material,
        specification: StockKind::Profile {
            designation: name.into(),
            length_mm: 6000.,
            width_mm: w,
            height_mm: h,
            wall_mm: Some(wall),
            shape,
            section: None,
            corner_radius_mm: None,
            web_mm: (shape == ProfileShape::Channel).then_some(4.5),
        },
    })
    .entity(EntityId::new(), name)
    .unwrap()
}
fn session() -> (Session, Stocks) {
    let mut s = Session::new(Document::new("Frames").unwrap()).unwrap();
    let steel = Entity::new(
        "S235 steel",
        EntityData::Material {
            density_kg_m3: Some(7850.),
            color_rgba: [120, 124, 130, 255],
            section_hatch: None,
        },
    );
    let entities = vec![
        stock(
            steel.id,
            "SHS 40×40×2",
            ProfileShape::Rectangular,
            40.,
            40.,
            2.,
        ),
        stock(
            steel.id,
            "SHS 30×30×2",
            ProfileShape::Rectangular,
            30.,
            30.,
            2.,
        ),
        stock(
            steel.id,
            "CHS 33.7×2.6",
            ProfileShape::Round,
            33.7,
            33.7,
            2.6,
        ),
        stock(steel.id, "CHS 60.3×3", ProfileShape::Round, 60.3, 60.3, 3.),
        stock(steel.id, "UPE 100", ProfileShape::Channel, 55., 100., 7.5),
    ];
    let ids: Vec<EntityId> = entities.iter().map(|e| e.id).collect();
    let mut all = vec![steel];
    all.extend(entities);
    s.commit(prepare_library(&s, 0, actor(), &all).unwrap())
        .unwrap();
    (
        s,
        Stocks {
            shs: ids[0],
            small: ids[1],
            chs: ids[2],
            big_chs: ids[3],
            channel: ids[4],
        },
    )
}
fn seg(l: f64, a: f64) -> ChainSegment {
    ChainSegment {
        length: length(l),
        angle: angle(a),
        rise: None,
    }
}
fn origin(frame: &mut FrameDefinition, u: f64, v: f64) -> EntityId {
    frame.add_node(NodePosition::Absolute {
        u: length(u),
        v: length(v),
        w: None,
    })
}
/// Commit a frame (new or edited), giving every new member `stock` and a
/// name from `names` in order.
fn commit_frame(
    s: &mut Session,
    id: EntityId,
    definition: &FrameDefinition,
    new: &[(EntityId, &str, EntityId)],
) {
    let prepared = prepare(s, id, definition, new);
    checked(s, prepared);
}
fn prepare(
    s: &Session,
    id: EntityId,
    definition: &FrameDefinition,
    new: &[(EntityId, &str, EntityId)],
) -> Prepared {
    model::prepare_frame(
        s,
        s.revision(),
        actor(),
        &FrameRequest {
            frame: id,
            name: "Frame".into(),
            definition: definition.clone(),
            members: new
                .iter()
                .map(|(m, name, stock)| FrameMemberRequest::new(*m, *name, *stock))
                .collect(),
        },
    )
    .unwrap()
    .prepared
}
/// A closed 1000 × 600 rectangle of `stock` drawn from the origin by
/// lengths and angles, side = inner.
fn rectangle(s: &mut Session, stock: EntityId, width: Dimension) -> (EntityId, FrameDefinition) {
    let mut frame = FrameDefinition::new(FramePlane::xy());
    frame.side = FrameSide::Inner;
    let start = origin(&mut frame, 0., 0.);
    let members = frame
        .add_chain(
            s.parameters(),
            start,
            &[
                ChainSegment {
                    length: width.clone(),
                    angle: angle(0.),
                    rise: None,
                },
                seg(600., 90.),
                ChainSegment {
                    length: width,
                    angle: angle(180.),
                    rise: None,
                },
            ],
            true,
        )
        .unwrap();
    assert_eq!(members.len(), 4);
    let names = ["Front", "Right", "Back", "Left"];
    let new: Vec<_> = members
        .iter()
        .zip(names)
        .map(|(m, n)| (*m, n, stock))
        .collect();
    let id = EntityId::new();
    commit_frame(s, id, &frame, &new);
    (id, frame)
}
fn shs_area() -> f64 {
    1600. - 36. * 36.
}

#[test]
fn a_closed_inner_rectangle_is_exactly_1000_by_600_outside_with_four_mitres() {
    let (mut s, st) = session();
    let (frame, _) = rectangle(&mut s, st.shs, length(1000.));
    let m = compile(&s);
    let list = m.member_cut_list(s.document()).unwrap();
    assert_eq!(list.len(), 4, "four member Parts");
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for r in &list {
        assert!(r.issues.is_empty(), "{}: {:?}", r.name, r.issues);
        let b = m
            .solid(m.snapshot().bodies[&r.part])
            .unwrap()
            .bounds()
            .unwrap();
        for k in 0..3 {
            lo[k] = lo[k].min(b.min[k]);
            hi[k] = hi[k].max(b.max[k]);
        }
        let nominal = if ["Front", "Back"].contains(&r.name.as_str()) {
            1000.
        } else {
            600.
        };
        near(r.length_mm, nominal);
        for end in &r.ends {
            assert_eq!(end.kind, Some(ConnectionKind::Miter), "{}", r.name);
            assert_eq!(end.connection, Some(frame));
            near(end.angle_deg, 45.);
            // Long points on the outer face: the section width runs inward.
            assert_eq!(end.long_point.unwrap().face, "-x", "{}", r.name);
        }
        // A 45° mitre at each end: area × length at the section centroid.
        rel(r.volume_mm3, shs_area() * (nominal - 40.));
        rel(r.mass_kg.unwrap(), r.volume_mm3 * 1e-9 * 7850.);
    }
    near(lo[0], 0.);
    near(lo[1], 0.);
    near(hi[0], 1000.);
    near(hi[1], 600.);
    near(hi[2] - lo[2], 40.);
}

#[test]
fn a_parameter_drives_a_polar_length_keeping_member_ids_with_undo_and_redo() {
    let (mut s, st) = session();
    let parameter = Entity::new(
        "width",
        EntityData::Parameter {
            expression: "1000".into(),
            unit: "mm".into(),
            exposed: true,
            inputs: vec![],
        },
    );
    let parameter_id = parameter.id;
    let p = s
        .prepare(Request {
            expected_revision: s.revision(),
            actor: actor(),
            label: "Add width".into(),
            edits: vec![Mutation::Insert(parameter)],
        })
        .unwrap();
    s.commit(p).unwrap();
    let width = Dimension::Expression {
        expression: BoundExpression::bind(s.document(), "=$width", "mm", None).unwrap(),
    };
    let (frame, _) = rectangle(&mut s, st.shs, width);
    let ids = |s: &Session| {
        let mut v = model::frame_member_features(s.document(), frame);
        v.sort();
        v.into_iter()
            .map(|(m, f)| {
                (
                    m,
                    f,
                    model::read_member(s.document(), f).unwrap().output.body,
                )
            })
            .collect::<Vec<_>>()
    };
    let before_ids = ids(&s);
    let before = s.document().clone();
    let mut p = s.document().get(parameter_id).unwrap().clone();
    let EntityData::Parameter { expression, .. } = &mut p.data else {
        unreachable!()
    };
    *expression = "1200".into();
    let prepared = s
        .prepare(Request {
            expected_revision: s.revision(),
            actor: actor(),
            label: "Edit width".into(),
            edits: vec![Mutation::Replace(p)],
        })
        .unwrap();
    checked(&mut s, prepared);
    let after = s.document().clone();
    let list = compile(&s).member_cut_list(s.document()).unwrap();
    near(record(&list, "Front").length_mm, 1200.);
    near(record(&list, "Back").length_mm, 1200.);
    near(record(&list, "Right").length_mm, 600.);
    assert_eq!(ids(&s), before_ids, "member feature and Part ids are kept");
    s.undo(s.revision(), actor()).unwrap();
    assert_eq!(s.document(), &before);
    near(
        record(&compile(&s).member_cut_list(s.document()).unwrap(), "Front").length_mm,
        1000.,
    );
    s.redo(s.revision(), actor()).unwrap();
    assert_eq!(s.document(), &after);
}

/// An L corner of two members drawn from one node.
fn corner(
    s: &mut Session,
    a: EntityId,
    b: EntityId,
) -> (EntityId, FrameDefinition, EntityId, [EntityId; 2]) {
    let mut frame = FrameDefinition::new(FramePlane::xy());
    let start = origin(&mut frame, 0., 0.);
    let first = frame
        .add_chain(
            s.parameters(),
            start,
            &[seg(800., 0.), seg(500., 90.)],
            false,
        )
        .unwrap();
    let node = frame.member(first[0]).unwrap().end;
    let id = EntityId::new();
    commit_frame(s, id, &frame, &[(first[0], "A", a), (first[1], "B", b)]);
    (id, frame, node, [first[0], first[1]])
}

#[test]
fn defaults_cope_round_branches_butt_smaller_ends_and_overrides_persist() {
    // A round branch onto a larger round tube is coped, the larger tube
    // running through. (The Rust kernel copes round tubes only onto a
    // receiving tube that runs past the joint, so this is a T.)
    let (mut s, st) = session();
    let mut frame = FrameDefinition::new(FramePlane::xy());
    let start = origin(&mut frame, 0., 0.);
    let rail = frame
        .add_chain(s.parameters(), start, &[seg(800., 0.)], false)
        .unwrap()[0];
    let station = frame.add_station(rail, length(400.), false).unwrap();
    let branch = frame
        .add_chain(s.parameters(), station, &[seg(500., 90.)], false)
        .unwrap()[0];
    commit_frame(
        &mut s,
        EntityId::new(),
        &frame,
        &[(rail, "A", st.big_chs), (branch, "B", st.chs)],
    );
    let list = compile(&s).member_cut_list(s.document()).unwrap();
    let (a, b) = (record(&list, "A"), record(&list, "B"));
    assert_eq!(b.ends[0].kind, Some(ConnectionKind::Cope));
    assert_eq!(b.ends[0].partner, Some(a.part));
    assert!(
        a.ends.iter().all(|e| e.kind.is_none()),
        "the larger tube runs through"
    );
    assert!(b.issues.is_empty() && a.issues.is_empty());

    // A 30×30 end on a 40×40 is butted with the 40×40 through.
    let (mut s, st) = session();
    let (frame, mut definition, node, [ma, mb]) = corner(&mut s, st.shs, st.small);
    let list = compile(&s).member_cut_list(s.document()).unwrap();
    let (a, b) = (record(&list, "A"), record(&list, "B"));
    assert_eq!(b.ends[0].kind, Some(ConnectionKind::Butt));
    assert!(b.ends[0].terminates);
    assert_eq!(a.ends[1].kind, Some(ConnectionKind::Butt));
    assert!(!a.ends[1].terminates, "the 40×40 runs through");
    near(b.length_mm, 500. - 20.);
    near(a.length_mm, 800. + 15.);

    // Override: the 30×30 runs through. It is stored and survives an edit.
    definition.set_joint(JointOverride {
        node,
        through: Some(mb),
        ..Default::default()
    });
    commit_frame(&mut s, frame, &definition, &[]);
    let tail = definition.member(mb).unwrap().end;
    definition
        .set_node(
            tail,
            NodePosition::Polar {
                from: node,
                length: length(650.),
                angle: angle(90.),
                rise: None,
            },
        )
        .unwrap();
    commit_frame(&mut s, frame, &definition, &[]);
    let stored = model::read_frame(s.document(), frame).unwrap();
    assert_eq!(stored.joint(node).unwrap().through, Some(mb));
    let list = compile(&s).member_cut_list(s.document()).unwrap();
    let (a, b) = (record(&list, "A"), record(&list, "B"));
    assert!(a.ends[1].terminates && !b.ends[0].terminates);
    near(b.length_mm, 650. + 20.);
    near(a.length_mm, 800. - 15.);

    // Deleting the through member reports the stale override, never
    // replacing it with an automatic joint.
    definition.remove_member(mb).unwrap();
    let prepared = prepare(&s, frame, &definition, &[]);
    let message = refused(&s, &prepared);
    assert!(message.contains("joint override"), "{message}");
    let _ = ma;
}

#[test]
fn a_gap_on_one_mitre_shortens_both_members_by_their_share() {
    let (mut s, st) = session();
    let (frame, mut definition) = rectangle(&mut s, st.shs, length(1000.));
    // The corner where Front meets Right.
    let front = definition.members[0].clone();
    definition.set_joint(JointOverride {
        node: front.end,
        gap: Some(length(2.)),
        ..Default::default()
    });
    commit_frame(&mut s, frame, &definition, &[]);
    let list = compile(&s).member_cut_list(s.document()).unwrap();
    // Each cut moves back 1 mm square to the mitre plane: √2 mm along a
    // member at 45°, measured at the long point.
    let share = std::f64::consts::SQRT_2;
    near(record(&list, "Front").length_mm, 1000. - share);
    near(record(&list, "Right").length_mm, 600. - share);
    near(record(&list, "Back").length_mm, 1000.);
    near(record(&list, "Left").length_mm, 600.);
}

#[test]
fn stations_make_t_and_x_joints_and_a_straight_chain_stays_one_member() {
    let (mut s, st) = session();
    let mut frame = FrameDefinition::new(FramePlane::xy());
    let start = origin(&mut frame, 0., 0.);
    // Two collinear segments: one member.
    let rail = frame
        .add_chain(
            s.parameters(),
            start,
            &[seg(600., 0.), seg(400., 0.)],
            false,
        )
        .unwrap();
    assert_eq!(rail.len(), 1, "a straight chain stays one member");
    let rail = rail[0];
    // A T at 300 mm and an X at 700 mm.
    let t = frame.add_station(rail, length(300.), false).unwrap();
    let x = frame.add_station(rail, length(700.), false).unwrap();
    let branch = frame
        .add_chain(s.parameters(), t, &[seg(500., 90.)], false)
        .unwrap()[0];
    let up = frame
        .add_chain(s.parameters(), x, &[seg(500., 90.)], false)
        .unwrap()[0];
    let down = frame
        .add_chain(s.parameters(), x, &[seg(500., 270.)], false)
        .unwrap()[0];
    let id = EntityId::new();
    commit_frame(
        &mut s,
        id,
        &frame,
        &[
            (rail, "Rail", st.shs),
            (branch, "Branch", st.small),
            (up, "Up", st.small),
            (down, "Down", st.small),
        ],
    );
    let m = compile(&s);
    let list = m.member_cut_list(s.document()).unwrap();
    assert_eq!(list.len(), 4);
    let rail_r = record(&list, "Rail");
    near(rail_r.length_mm, 1000.);
    assert!(
        rail_r.ends.iter().all(|e| e.kind.is_none()),
        "the rail runs through"
    );
    rel(rail_r.volume_mm3, shs_area() * 1000.);
    for name in ["Branch", "Up", "Down"] {
        let r = record(&list, name);
        assert_eq!(r.ends[0].kind, Some(ConnectionKind::Butt), "{name}");
        assert_eq!(r.ends[0].partner, Some(rail_r.part));
        near(r.length_mm, 500. - 20.);
        rel(r.volume_mm3, (900. - 26. * 26.) * 480.);
    }
}

#[test]
fn a_table_frame_puts_legs_through_its_rectangle_corners() {
    let (mut s, st) = session();
    let (frame, mut definition) = rectangle(&mut s, st.shs, length(1000.));
    let corners: Vec<EntityId> = definition.members.iter().map(|m| m.start).collect();
    let mut new = vec![];
    for (i, corner) in corners.iter().enumerate() {
        let leg = definition
            .add_chain(
                s.parameters(),
                *corner,
                &[ChainSegment {
                    length: length(0.),
                    angle: angle(0.),
                    rise: Some(length(-700.)),
                }],
                false,
            )
            .unwrap()[0];
        new.push((leg, ["Leg 1", "Leg 2", "Leg 3", "Leg 4"][i], st.shs));
    }
    commit_frame(&mut s, frame, &definition, &new);
    let m = compile(&s);
    let list = m.member_cut_list(s.document()).unwrap();
    assert_eq!(list.len(), 8, "eight separate Parts");
    let mut lo = [f64::INFINITY; 2];
    let mut hi = [f64::NEG_INFINITY; 2];
    for r in &list {
        assert!(r.issues.is_empty(), "{}: {:?}", r.name, r.issues);
        let b = m
            .solid(m.snapshot().bodies[&r.part])
            .unwrap()
            .bounds()
            .unwrap();
        for k in 0..2 {
            lo[k] = lo[k].min(b.min[k]);
            hi[k] = hi[k].max(b.max[k]);
        }
    }
    // The outside stays exactly 1000 × 600.
    near(lo[0], 0.);
    near(lo[1], 0.);
    near(hi[0], 1000.);
    near(hi[1], 600.);
    for i in 1..=4 {
        let leg = record(&list, &format!("Leg {i}"));
        // The leg runs through to the top of the rails: 700 mm plus the
        // 20 mm of rail above the frame plane, both ends square.
        near(leg.length_mm, 720.);
        assert!(leg.ends.iter().all(|e| e.angle_deg.abs() < 1e-9));
        assert!(!leg.ends[0].terminates);
        rel(leg.volume_mm3, shs_area() * 720.);
    }
    for (name, span) in [
        ("Front", 1000.),
        ("Right", 600.),
        ("Back", 1000.),
        ("Left", 600.),
    ] {
        let r = record(&list, name);
        // Each rail butts square onto the two legs at its ends.
        for end in &r.ends {
            assert_eq!(end.kind, Some(ConnectionKind::Butt), "{name}");
            assert!(end.terminates);
        }
        near(r.length_mm, span - 80.);
        rel(r.volume_mm3, shs_area() * (span - 80.));
    }
}

#[test]
fn joints_a_section_or_angle_cannot_take_are_refused_at_the_node() {
    // Members meeting at 10°: refused with the reason.
    let (mut s, st) = session();
    let mut frame = FrameDefinition::new(FramePlane::xy());
    let start = origin(&mut frame, 0., 0.);
    let members = frame
        .add_chain(
            s.parameters(),
            start,
            &[seg(500., 0.), seg(500., 170.)],
            false,
        )
        .unwrap();
    let new: Vec<_> = members.iter().map(|m| (*m, "Bar", st.shs)).collect();
    let new = [(new[0].0, "Bar 1", st.shs), (new[1].0, "Bar 2", st.shs)];
    let prepared = prepare(&s, EntityId::new(), &frame, &new);
    let message = refused(&s, &prepared);
    assert!(message.contains("at least 15°"), "{message}");
    assert!(message.contains("frame node"), "{message}");

    // A tube branch on a channel would be coped against it: refused.
    let mut frame = FrameDefinition::new(FramePlane::xy());
    let start = origin(&mut frame, 0., 0.);
    let rail = frame
        .add_chain(s.parameters(), start, &[seg(1000., 0.)], false)
        .unwrap()[0];
    let station = frame.add_station(rail, length(500.), false).unwrap();
    let branch = frame
        .add_chain(s.parameters(), station, &[seg(400., 90.)], false)
        .unwrap()[0];
    let prepared = prepare(
        &s,
        EntityId::new(),
        &frame,
        &[(rail, "Channel", st.channel), (branch, "Post", st.shs)],
    );
    let message = refused(&s, &prepared);
    assert!(message.contains("Copes against a channel"), "{message}");

    // Two separate nodes at one place: their members do not meet.
    let mut frame = FrameDefinition::new(FramePlane::xy());
    let a = origin(&mut frame, 0., 0.);
    let first = frame
        .add_chain(s.parameters(), a, &[seg(500., 0.)], false)
        .unwrap();
    let b = origin(&mut frame, 500., 0.);
    let error = frame
        .add_chain(s.parameters(), b, &[seg(300., 90.)], false)
        .unwrap_err()
        .to_string();
    assert!(error.contains("do not meet"), "{error}");
    let _ = (first, &mut s);
}

#[test]
fn frames_round_trip_in_files_requiring_their_capability() {
    let (mut s, st) = session();
    let project = ProjectFile::new(s.document().clone()).unwrap();
    assert!(!project.capabilities().contains(model::FRAME_CAPABILITY));
    let (frame, definition) = rectangle(&mut s, st.shs, length(1000.));
    let project = ProjectFile::new(s.document().clone()).unwrap();
    assert!(project.capabilities().contains(model::FRAME_CAPABILITY));
    let mut bytes = vec![];
    project_file::write(&project, &mut bytes, Default::default()).unwrap();
    let opened = project_file::read(
        &mut bytes.as_slice(),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap()
    .project;
    assert_eq!(
        model::read_frame(&opened.document, frame).unwrap(),
        definition
    );
    let reopened = Compiler::default().compile(&opened.document, 0).unwrap();
    let list = reopened.member_cut_list(&opened.document).unwrap();
    near(record(&list, "Front").length_mm, 1000.);
    let mut older = project_file::supported_capabilities();
    older.remove(model::FRAME_CAPABILITY);
    assert!(project_file::read(&mut bytes.as_slice(), Default::default(), &older).is_err());
}

#[test]
fn centroid_placement_offset_and_roll_place_the_section_exactly() {
    let (mut s, _) = session();
    let steel = s
        .document()
        .entities()
        .find(|e| matches!(e.data, EntityData::Material { .. }))
        .unwrap()
        .id;
    let angle_stock = stock(steel, "L 40×40×4", ProfileShape::Angle, 40., 40., 4.);
    let angle_id = angle_stock.id;
    s.commit(prepare_library(&s, s.revision(), actor(), &[angle_stock]).unwrap())
        .unwrap();
    let build = |s: &mut Session, side: FrameSide, offset: Option<[f64; 2]>, roll: Option<f64>| {
        let mut frame = FrameDefinition::new(FramePlane::xy());
        let a = origin(&mut frame, 0., 0.);
        let m = frame
            .add_chain(s.parameters(), a, &[seg(500., 0.)], false)
            .unwrap()[0];
        let member = frame.members.iter_mut().find(|x| x.id == m).unwrap();
        member.side = Some(side);
        member.offset = offset.map(|o| o.map(length));
        member.roll = roll.map(angle);
        let id = EntityId::new();
        commit_frame(s, id, &frame, &[(m, "Angle", angle_id)]);
        let model = compile(s);
        let r = model
            .member_cut_list(s.document())
            .unwrap()
            .into_iter()
            .find(|r| {
                r.path
                    == model::MemberPath::FrameMember {
                        frame: id,
                        member: m,
                    }
            })
            .unwrap();
        model
            .solid(model.snapshot().bodies[&r.part])
            .unwrap()
            .bounds()
            .unwrap()
    };
    // An equal angle 40×40×4: its centroid from the heel on both legs.
    let area = 40. * 4. + 36. * 4.;
    let centroid = (40. * 4. * 20. + 36. * 4. * 2.) / area;
    // Section x runs along +y (world), section y along +z for a +x member.
    let b = build(&mut s, FrameSide::Centroid, None, None);
    near(-b.min[1], centroid);
    near(-b.min[2], centroid);
    // An offset of 10 mm along the section width moves it 10 mm along +y.
    let b = build(&mut s, FrameSide::Centroid, Some([10., 0.]), None);
    near(-b.min[1], centroid - 10.);
    near(-b.min[2], centroid);
    // A 90° roll about the axis through the centroid: the section width
    // turns to +z and its height to −y.
    let b = build(&mut s, FrameSide::Centroid, None, Some(90.));
    near(-b.min[1], 40. - centroid);
    near(-b.min[2], centroid);
}
