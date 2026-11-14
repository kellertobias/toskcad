use super::*;
use crate::step::{
    Writer,
    insert::{CurveOptions, DomainOptions, DomainPlan, SurfacePlan, TopologyPlan},
};
pub(super) fn compound() -> String {
    let body = super::super::tests::graph(false, 4.);
    let d = DomainPlan::prepare(&body, DomainOptions::default(), None).unwrap();
    let s = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
    let t = TopologyPlan::prepare(&s, 16_000_000, 64 * 1024 * 1024, None).unwrap();
    let mut writer = Writer::new();
    writer
        .append_sphere(
            crate::sphere::Sphere::new(1.).unwrap(),
            [20., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
        )
        .unwrap();
    for _ in 0..2 {
        writer
            .append_owned_insert_local(&t, 16_000_000, 64 * 1024 * 1024, None)
            .unwrap();
    }
    writer.finish().unwrap()
}
#[test]
fn actual_compound_members_share_one_arena_and_earn_only_their_own_carriers() {
    let text = compound();
    let file = File::parse(&text, Options::default(), None).unwrap();
    assert_eq!(file.literal_input(), text);
    assert_eq!(file.represented_solids().len(), 3);
    let shell = Shell::new(&file.entities, 1e-7, false);
    let ids = file
        .represented_solids()
        .iter()
        .copied()
        .filter(|id| shell.member_faces(*id).unwrap().len() > 2)
        .collect::<Vec<_>>();
    assert_eq!(ids.len(), 2);
    let a = file.member(ids[0], Options::default(), None).unwrap();
    let b = file.member(ids[1], Options::default(), None).unwrap();
    assert!(file.shares_arena(&a) && file.shares_arena(&b));
    assert!(Arc::ptr_eq(&a.entities, &b.entities));
    assert_eq!(a.member_entity(), ids[0]);
    assert_eq!(b.member_entity(), ids[1]);
    for p in [&a, &b] {
        assert_eq!(p.planes().len(), 2);
        let bound = p.support_bindings(Options::default(), None).unwrap();
        let material = bound.material_bounds(Options::default(), None).unwrap();
        assert_eq!(material.bounds_mm()[2], [0., 4.]);
        eprintln!(
            "MEMBER work{} bytes{}",
            material.work_used(),
            material.retained_budget_bytes()
        );
    }
    assert!(matches!(
        file.member(u32::MAX, Options::default(), None),
        Err(Error::Unsupported)
    ));
    let sphere = file
        .represented_solids()
        .iter()
        .copied()
        .find(|id| !ids.contains(id))
        .unwrap();
    assert!(matches!(
        file.member(sphere, Options::default(), None),
        Err(Error::Unsupported)
    ));
    // The strict standalone API keeps its earned one-solid contract.
    assert!(matches!(
        parse(&text, Options::default(), None),
        Err(Error::Unsupported)
    ));
    drop(file);
    assert_eq!(a.literal_input(), text);
    assert_eq!(b.member_entity(), ids[1]);
}
#[test]
fn actual_file_and_selected_query_limits_cancellation_and_dangling_references_are_atomic() {
    let text = super::super::tests::file(false);
    let file = File::parse(&text, Options::default(), None).unwrap();
    let id = file.represented_solids()[0];
    let another = File::parse(&text, Options::default(), None).unwrap();
    let p = file.member(id, Options::default(), None).unwrap();
    assert!(!another.shares_arena(&p));
    assert!(matches!(
        File::parse(
            &text,
            Options {
                retained_budget_bytes: file.retained_budget_bytes() - 1,
                ..Options::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    for options in [
        Options {
            maximum_work: 0,
            ..Options::default()
        },
        Options {
            retained_budget_bytes: 0,
            ..Options::default()
        },
    ] {
        assert!(matches!(
            file.member(id, options, None),
            Err(Error::Invalid)
        ));
    }
    assert!(matches!(
        file.member(
            id,
            Options {
                maximum_work: file.work_used(),
                ..Options::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        file.member(
            id,
            Options {
                retained_budget_bytes: file.retained_budget_bytes(),
                ..Options::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        file.member(
            id,
            Options {
                maximum_work: file.work_used() + 100,
                ..Options::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    let polls = std::cell::Cell::new(0);
    let count = || {
        polls.set(polls.get() + 1);
        false
    };
    file.member(id, Options::default(), Some(&count)).unwrap();
    let total = polls.get();
    assert!(total > 100);
    for at in [1, total - 2] {
        let polls = std::cell::Cell::new(0);
        let stop = || {
            polls.set(polls.get() + 1);
            polls.get() == at
        };
        assert!(matches!(
            file.member(id, Options::default(), Some(&stop)),
            Err(Error::Cancelled)
        ));
    }
    let plane = &p.planes()[0];
    let broken = super::super::tests::argument_change(&text, plane.axis_entity(), 1, "#4294967295");
    let f = File::parse(&broken, Options::default(), None).unwrap();
    assert!(matches!(
        f.member(id, Options::default(), None),
        Err(Error::Invalid)
    ));
    let typed = super::super::tests::argument_change(
        &text,
        plane.entity(),
        0,
        "UNSUPPORTED_LABEL(#4294967295)",
    );
    let f = File::parse(&typed, Options::default(), None).unwrap();
    assert!(matches!(
        f.member(id, Options::default(), None),
        Err(Error::Invalid)
    ));
    // A failed query has not altered the shared actual source arena.
    let again = file.member(id, Options::default(), None).unwrap();
    assert!(file.shares_arena(&again));
    assert_eq!(again.member_entity(), id);
}
pub(super) fn mapped_file() -> (String, u32, u32, u32) {
    let text = super::super::tests::file(false);
    let file = File::parse(&text, Options::default(), None).unwrap();
    let solid = file.represented_solids()[0];
    let (&rep, records) = file
        .entities
        .iter()
        .find(|(_, r)| matches!(r.as_slice(),[e]if e.name=="ADVANCED_BREP_SHAPE_REPRESENTATION"))
        .unwrap();
    let e = &records[0];
    let context = reference(&e.args[2]).unwrap();
    let origin = refs(&e.args[1])
        .unwrap()
        .into_iter()
        .find(|&id| id != solid)
        .unwrap();
    let records = format!(
        "#1000000=CARTESIAN_POINT('',(10.,20.,30.));\n#1000001=DIRECTION('',(1.,0.,0.));\n#1000002=DIRECTION('',(0.,1.,0.));\n#1000003=DIRECTION('',(0.,0.,1.));\n#1000004=(CARTESIAN_TRANSFORMATION_OPERATOR(#1000001,#1000002,#1000000,1.)CARTESIAN_TRANSFORMATION_OPERATOR_3D(#1000003)FUNCTIONALLY_DEFINED_TRANSFORMATION('','')GEOMETRIC_REPRESENTATION_ITEM()REPRESENTATION_ITEM(''));\n#1000005=SHAPE_REPRESENTATION('',(#{solid},#{origin}),#{context});\n#1000006=REPRESENTATION_MAP(#{origin},#1000005);\n#1000007=MAPPED_ITEM('',#1000006,#1000004);\nENDSEC;\nEND-ISO-10303-21;"
    );
    let text = text.replacen("ENDSEC;\nEND-ISO-10303-21;", &records, 1);
    let text = super::super::tests::argument_change(&text, rep, 1, "(#1000007)");
    (text, solid, origin, rep)
}
#[test]
fn actual_mapped_member_preserves_local_body_and_declared_target_ownership() {
    let (text, solid, origin, _) = mapped_file();
    let file = File::parse(&text, Options::default(), None).unwrap();
    assert!(file.represented_solids().is_empty());
    assert_eq!(file.represented_members().len(), 1);
    let member = file.represented_members()[0];
    assert_eq!(member.item_entity(), 1000007);
    assert_eq!(member.local_solid_entity(), solid);
    let link = member.mapping_link().unwrap();
    assert_eq!(link.map_entity(), 1000006);
    assert_eq!(link.source_representation_entity(), 1000005);
    assert_eq!(link.source_origin_entity(), origin);
    assert_eq!(link.operator_entity(), 1000004);
    let p = file
        .member(member.item_entity(), Options::default(), None)
        .unwrap();
    assert_eq!(p.member_entity(), member.item_entity());
    assert_eq!(p.local_solid_entity(), solid);
    assert_eq!(p.mapping_link(), Some(link));
    let b = p.support_bindings(Options::default(), None).unwrap();
    let material = b.material_bounds(Options::default(), None).unwrap();
    assert_eq!(material.bounds_mm()[2], [0., 4.]);
    let op = p
        .physical_mapping(link.operator_entity(), Options::default(), None)
        .unwrap();
    assert_eq!(op.handedness().unwrap(), 1);
    let bounds = op.map_box_mm(material.bounds_mm()).unwrap();
    assert!(bounds[2][0] <= 30. && bounds[2][1] >= 34.);
    let volume = material
        .occupied_volume(FluxOptions::default(), None)
        .unwrap();
    assert!(volume.volume_mm3()[0] > 0.);
    assert!(op.volume_scale_bounds().unwrap()[0] <= 1.);
    // Strict standalone parsing must never drop an actual represented transform.
    assert!(matches!(
        parse(&text, Options::default(), None),
        Err(Error::Unsupported)
    ));
}
#[test]
fn actual_nonidentity_source_origin_and_nested_map_are_honest_refusals() {
    let (text, _, origin, _) = mapped_file();
    let file = File::parse(&text, Options::default(), None).unwrap();
    let axis = kind(&file.entities, origin, "AXIS2_PLACEMENT_3D").unwrap();
    let point = reference(&axis.args[1]).unwrap();
    let shifted = super::super::tests::argument_change(&text, point, 1, "(1.E-1000,0.,0.)");
    assert!(matches!(
        File::parse(&shifted, Options::default(), None),
        Err(Error::Unsupported)
    ));
    let nested =
        super::super::tests::argument_change(&text, 1000005, 1, &format!("(#1000007,#{origin})"));
    assert!(matches!(
        File::parse(&nested, Options::default(), None),
        Err(Error::Unsupported)
    ));
}
