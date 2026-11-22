use super::*;
use crate::step::{
    Writer,
    insert::{CurveOptions, DomainOptions, DomainPlan, SurfacePlan, TopologyPlan},
};
use multicad_precision::PortableMath as _;
pub(super) fn mixed(left: bool) -> String {
    let body = crate::step_read::owned_insert::tests::graph(left, 4.);
    let d = DomainPlan::prepare(&body, DomainOptions::default(), None).unwrap();
    let s = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
    let t = TopologyPlan::prepare(&s, 16_000_000, 64 * 1024 * 1024, None).unwrap();
    let gear = gear();
    let mut w = Writer::new();
    w.append_herringbone(&gear, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7)
        .unwrap();
    w.append_owned_insert_mapped_candidate(
        &t,
        [31., -17., 9.],
        [0.8, 0.6, 0.],
        [-0.36, 0.48, 0.8],
        16_000_000,
        64 * 1024 * 1024,
        None,
    )
    .unwrap();
    w.finish().unwrap()
}
#[test]
fn actual_twenty_four_tooth_compound_has_a_separate_global_and_selected_budget() {
    for left in [false, true] {
        let text = mixed(left);
        assert!(matches!(
            File::parse(&text, Options::default(), None),
            Err(Error::Limit)
        ));
        let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
        eprintln!(
            "DOCUMENT input={} work={} bytes={}",
            text.len(),
            document.work_used(),
            document.retained_budget_bytes()
        );
        assert!(document.retained_budget_bytes() > 64 * 1024 * 1024);
        assert!(document.work_used() < 16_000_000); // Proves old rejection was retained bytes, not work.
        assert_eq!(document.literal_input(), text);
        let member = document.represented_members()[1];
        assert!(member.mapping_link().is_some());
        let packet = document
            .member_packet(member.item_entity(), Options::default(), None)
            .unwrap();
        assert!(std::ptr::eq(packet.document(), &document));
        assert_eq!(packet.parsed().mapping_link(), member.mapping_link());
        assert_eq!(packet.parsed().literal_input(), text);
        assert!(packet.query_budget_bytes() < 64 * 1024 * 1024);
        let parsed = packet.parsed();
        let bindings = parsed.support_bindings(Options::default(), None).unwrap();
        let material = bindings.material_bounds(Options::default(), None).unwrap();
        let volume = material
            .occupied_volume(FluxOptions::default(), None)
            .unwrap();
        assert!(volume.volume_mm3()[0] > 0.);
        let link = member.mapping_link().unwrap();
        let mapping = parsed
            .physical_mapping(link.operator_entity(), Options::default(), None)
            .unwrap();
        assert_eq!(mapping.handedness().unwrap(), 1);
        let view = document
            .local_member_document(
                document.represented_members()[0].item_entity(),
                Options::default(),
                None,
            )
            .unwrap();
        assert!(view.query_budget_bytes() < 64 * 1024 * 1024);
        assert!(view.work_used() < 16_000_000);
        let ordinary = crate::step_read::verify_herringbone(
            view.as_str(),
            &gear(),
            [0.; 3],
            [1., 0., 0.],
            [0., 1., 0.],
            1e-7,
        );
        assert!(ordinary.is_ok(), "{ordinary:?}");
    }
}
fn gear() -> crate::spur_gear::herringbone::Herringbone {
    let beta = 15_f64.to_radians();
    let pressure = (20_f64.to_radians().portable_tan() / beta.portable_cos()).portable_atan();
    let module = 2. / beta.portable_cos();
    let pitch = module * 24. / 2.;
    crate::spur_gear::herringbone::Herringbone::new_helical(
        crate::spur_gear::Spec {
            teeth: 24,
            pitch_radius_mm: pitch,
            base_radius_mm: pitch * pressure.portable_cos(),
            tip_radius_mm: pitch + 2.,
            root_radius_mm: pitch - 2.5,
            root_fillet_radius_mm: 0.4,
            half_base_angle_rad: (std::f64::consts::PI * module / 2. - 0.05) / (2. * pitch)
                + pressure.portable_tan()
                - pressure,
            face_width_mm: 4.,
            bore_radius_mm: 4.,
        },
        beta.portable_tan() / pitch,
        1e-5,
        1e7,
    )
    .unwrap()
}
#[test]
fn actual_global_and_selected_query_frontiers_and_cancellation_stay_separate() {
    let text = super::super::tests::compound();
    let d = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let member = d.represented_members()[1].item_entity();
    let packet = d.member_packet(member, Options::default(), None).unwrap();
    // An actual selected-AST threshold, after closure scratch fits and before
    // allocating copied records, must refuse independently of final arrays.
    let count = d.source.entities.len();
    let set_bytes = (1 + count / 5) * (20 * size_of::<usize>() + 11 * size_of::<u32>());
    let fixed = size_of::<File<'_>>()
        + size_of::<MemberPacket<'_, '_>>()
        + size_of::<BTreeMap<u32, Vec<Entity<'_>>>>()
        + 2 * size_of::<usize>();
    let selected = &packet.parsed.entities;
    let mut w = Ledger {
        used: 0,
        remaining: 16_000_000,
        interrupt: None,
    };
    let mut preflight = fixed
        + (1 + selected.len() / 5)
            * (20 * size_of::<usize>() + 11 * (size_of::<u32>() + size_of::<Vec<Entity<'_>>>()));
    for id in selected.keys() {
        preflight += record_bytes(&d.source.entities[id], &mut w).unwrap();
    }
    assert!(preflight > fixed + count * size_of::<u32>());
    assert!(matches!(
        d.member_packet(
            member,
            Options {
                retained_budget_bytes: preflight + set_bytes - 1,
                ..Default::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    let measured = budget::arena("", selected, &mut w).unwrap() + fixed;
    assert!(measured <= preflight);
    for n in [1, 17, count] {
        let mut v = Vec::<u32>::new();
        v.try_reserve_exact(n).unwrap();
        assert_eq!(v.capacity(), n);
    }

    let view = d
        .local_member_document(
            d.represented_members()[0].item_entity(),
            Options::default(),
            None,
        )
        .unwrap();
    for options in [
        DocumentOptions {
            maximum_work: 0,
            ..Default::default()
        },
        DocumentOptions {
            maximum_work: 512_000_001,
            ..Default::default()
        },
        DocumentOptions {
            retained_budget_bytes: 0,
            ..Default::default()
        },
        DocumentOptions {
            retained_budget_bytes: 512 * 1024 * 1024 + 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            CompoundDocument::parse(&text, options, None),
            Err(Error::Invalid)
        ));
    }
    for options in [
        DocumentOptions {
            maximum_work: d.work_used() - 1,
            ..Default::default()
        },
        DocumentOptions {
            retained_budget_bytes: d.retained_budget_bytes() - 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            CompoundDocument::parse(&text, options, None),
            Err(Error::Limit)
        ));
    }
    CompoundDocument::parse(
        &text,
        DocumentOptions {
            maximum_work: d.work_used(),
            retained_budget_bytes: d.retained_budget_bytes(),
        },
        None,
    )
    .unwrap();
    for options in [
        Options {
            maximum_work: packet.work_used() - 1,
            ..Default::default()
        },
        Options {
            retained_budget_bytes: packet.query_budget_bytes() - 1,
            ..Default::default()
        },
        Options {
            retained_budget_bytes: 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            d.member_packet(member, options, None),
            Err(Error::Limit)
        ));
    }
    d.member_packet(
        member,
        Options {
            maximum_work: packet.work_used(),
            retained_budget_bytes: packet.query_budget_bytes(),
        },
        None,
    )
    .unwrap();
    for options in [
        Options {
            maximum_work: view.work_used() - 1,
            ..Default::default()
        },
        Options {
            retained_budget_bytes: view.query_budget_bytes() - 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            d.local_member_document(view.member().item_entity(), options, None),
            Err(Error::Limit)
        ));
    }
    d.local_member_document(
        view.member().item_entity(),
        Options {
            maximum_work: view.work_used(),
            retained_budget_bytes: view.query_budget_bytes(),
        },
        None,
    )
    .unwrap();
    assert!(matches!(
        d.member_packet(u32::MAX, Options::default(), None),
        Err(Error::Unsupported)
    ));
    assert!(matches!(
        d.member_packet(
            member,
            Options {
                maximum_work: 0,
                ..Default::default()
            },
            None
        ),
        Err(Error::Invalid)
    ));
    let calls = std::cell::Cell::new(0usize);
    let poll = || {
        calls.set(calls.get() + 1);
        false
    };
    CompoundDocument::parse(&text, DocumentOptions::default(), Some(&poll)).unwrap();
    let parse_calls = calls.get();
    calls.set(0);
    d.member_packet(member, Options::default(), Some(&poll))
        .unwrap();
    let packet_calls = calls.get();
    for stop in [0, parse_calls / 2, parse_calls - 1] {
        calls.set(0);
        assert!(matches!(
            CompoundDocument::parse(
                &text,
                DocumentOptions::default(),
                Some(&|| {
                    let n = calls.get();
                    calls.set(n + 1);
                    n >= stop
                })
            ),
            Err(Error::Cancelled)
        ));
    }
    for stop in [0, packet_calls / 2, packet_calls - 1] {
        calls.set(0);
        assert!(matches!(
            d.member_packet(
                member,
                Options::default(),
                Some(&|| {
                    let n = calls.get();
                    calls.set(n + 1);
                    n >= stop
                })
            ),
            Err(Error::Cancelled)
        ));
    }
    d.member_packet(member, Options::default(), None).unwrap();
}
#[test]
fn original_dependency_records_are_preserved_and_malformed_ownership_refuses() {
    let text = crate::step_read::owned_insert::tests::file(false);
    let d = CompoundDocument::parse(&text, Default::default(), None).unwrap();
    let item = d.represented_members()[0].item_entity();
    let p = d.member_packet(item, Options::default(), None).unwrap();
    assert!(!Arc::ptr_eq(&d.source.entities, &p.parsed.entities));
    for (&id, records) in p.parsed.entities.iter() {
        let original = &d.source.entities[&id];
        assert_eq!(records, original);
        for (copy, raw) in records.iter().zip(original) {
            assert_eq!(copy.source, raw.source);
        }
    }
    let vertex = p.parsed.vertices()[0].entity;
    let raw = d.source.entities[&vertex][0].source.unwrap();
    let changed = text.replacen(raw, "VERTEX_POINT('',#4294967295)", 1);
    let bad = CompoundDocument::parse(&changed, Default::default(), None).unwrap();
    assert!(matches!(
        bad.member_packet(item, Options::default(), None),
        Err(Error::Invalid)
    ));
    let context = text
        .lines()
        .find(|l| l.contains("ADVANCED_BREP_SHAPE_REPRESENTATION("))
        .unwrap();
    let changed = text.replacen(
        context,
        &context.replace(
            "ADVANCED_BREP_SHAPE_REPRESENTATION",
            "UNRELATED_REPRESENTATION",
        ),
        1,
    );
    assert!(matches!(
        CompoundDocument::parse(&changed, Default::default(), None),
        Err(Error::Unsupported)
    ));
    let mut oversized = String::new();
    oversized.try_reserve_exact(64 * 1024 * 1024 + 1).unwrap();
    oversized.extend(std::iter::repeat_n(' ', 64 * 1024 * 1024 + 1));
    assert!(matches!(
        CompoundDocument::parse(&oversized, Default::default(), None),
        Err(Error::Limit)
    ));
}
#[test]
fn original_metadata_copy_keeps_every_token_kind_and_interrupts() {
    let text = "ISO-10303-21;HEADER;FILE_DESCRIPTION(() , '2;1');FILE_NAME('','',(),(),'','','');FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));ENDSEC;DATA;#1=EXAMPLE((ITEM(*,$,'x',.T.,#1,1.0)));ENDSEC;END-ISO-10303-21;";
    let raw = part21::parse(text).unwrap();
    let mut w = Ledger {
        used: 0,
        remaining: 10000,
        interrupt: None,
    };
    let bytes = record_bytes(&raw[&1], &mut w).unwrap();
    assert!(bytes > 0);
    let copied = copy_records(&raw[&1], &mut w).unwrap();
    assert_eq!(copied, raw[&1]);
    assert_eq!(copied[0].source, raw[&1][0].source);
    let mut w = Ledger {
        used: 0,
        remaining: 0,
        interrupt: None,
    };
    assert!(matches!(copy_records(&raw[&1], &mut w), Err(Error::Limit)));
}
