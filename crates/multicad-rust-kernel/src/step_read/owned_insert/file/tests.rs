use super::*;
use crate::step::{
    Writer,
    insert::{CurveOptions, DomainOptions, DomainPlan, SurfacePlan, TopologyPlan},
};
fn compound() -> String {
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
