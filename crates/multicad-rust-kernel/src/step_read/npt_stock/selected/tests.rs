use super::*;
use crate::step_read::owned_insert::DocumentOptions;
fn compound(text: &str) -> (String, u32, u32) {
    let document = CompoundDocument::parse(text, DocumentOptions::default(), None).unwrap();
    let original = document.represented_members()[0].item_entity();
    let record = text
        .lines()
        .find(|line| line.starts_with(&format!("#{original}=")))
        .unwrap();
    let new_solid = 900_001u32;
    let unrelated = 900_002u32;
    let copied = record.replacen(&format!("#{original}="), &format!("#{new_solid}="), 1);
    let representation = text
        .lines()
        .find(|line| line.contains("=ADVANCED_BREP_SHAPE_REPRESENTATION("))
        .unwrap();
    let changed = representation.replacen(
        &format!("#{original}"),
        &format!("#{original},#{new_solid}"),
        1,
    );
    let text=text.replacen(representation,&changed,1).replacen("ENDSEC;\nEND-ISO-10303-21;",&format!("{copied}\n#{unrelated}=CARTESIAN_POINT('',(9999999999.,0.,0.));\nENDSEC;\nEND-ISO-10303-21;"),1);
    (text, new_solid, unrelated)
}
#[test]
fn selected_actual_members_ignore_unrelated_numeric_leaves_without_resetting_budgets() {
    let original = super::super::tests::fixture(false);
    let (text, second, unrelated) = compound(&original);
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    assert_eq!(document.represented_members().len(), 2);
    assert!(matches!(
        Parsed::parse(
            &document,
            document.represented_members()[0].item_entity(),
            Default::default(),
            None,
        ),
        Err(Error::Unsupported)
    ));
    let mut counts = Vec::new();
    for member in document.represented_members() {
        let selected =
            Selected::parse(&document, member.item_entity(), Default::default(), None).unwrap();
        assert_eq!(selected.member(), *member);
        assert!(selected.entity_ids().contains(&member.item_entity()));
        assert!(!selected.entity_ids().contains(&unrelated));
        assert_eq!(
            selected.parsed().solid_entity(),
            member.local_solid_entity()
        );
        assert!(
            selected
                .parsed()
                .points()
                .iter()
                .all(|p| p.entity() != unrelated)
        );
        counts.push((
            selected.parsed().points().len(),
            selected.parsed().work_used(),
        ));
        let occupied = selected
            .parsed()
            .projection(Default::default(), None)
            .unwrap()
            .occupied(Default::default(), None)
            .unwrap();
        assert!(occupied.volume_mm3()[0] > 470. && occupied.volume_mm3()[1] < 471.);
    }
    assert_eq!(counts[0], counts[1]);
    assert!(matches!(
        Selected::parse(&document, unrelated, Default::default(), None),
        Err(Error::Unsupported)
    ));
    assert!(matches!(
        Selected::parse(
            &document,
            second,
            Options {
                maximum_work: 1,
                ..Default::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        Selected::parse(&document, second, Default::default(), Some(&|| true)),
        Err(Error::Cancelled)
    ));
}
#[test]
fn selected_dependency_walk_checks_real_missing_typed_and_resource_boundaries() {
    let text = super::super::tests::fixture(false);
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let item = document.represented_members()[0].item_entity();
    for bytes in [
        super::super::SCRATCH + 327_680 * size_of::<u32>(),
        super::super::SCRATCH + 327_680 * size_of::<u32>() + 204,
    ] {
        assert!(matches!(
            Selected::parse(
                &document,
                item,
                Options {
                    maximum_bytes: bytes,
                    ..Default::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
    }
    let parsed = Parsed::parse(&document, item, Default::default(), None).unwrap();
    let point = parsed.points()
        [parsed.tensors()[0].point_indices()[11 * parsed.tensors()[0].spans()[1] + 2]];
    let old = text
        .lines()
        .find(|line| line.starts_with(&format!("#{}=", point.entity())))
        .unwrap();
    for coordinates in ["#999999,0.,0.", "LENGTH_MEASURE(#1),0.,0."] {
        let changed = format!("#{}=CARTESIAN_POINT('',({coordinates}));", point.entity());
        let wrong = text.replacen(old, &changed, 1);
        let document = CompoundDocument::parse(&wrong, DocumentOptions::default(), None).unwrap();
        assert!(matches!(
            Selected::parse(&document, item, Default::default(), None),
            Err(Error::Invalid)
        ));
    }
    let baseline = Selected::parse(&document, item, Default::default(), None).unwrap();
    let calls = std::cell::Cell::new(0usize);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    Selected::parse(&document, item, Default::default(), Some(&count)).unwrap();
    for checkpoint in [2, calls.get() / 4, calls.get() / 2, calls.get() - 1] {
        let calls = std::cell::Cell::new(0usize);
        let stop = || {
            calls.set(calls.get() + 1);
            calls.get() >= checkpoint
        };
        assert!(matches!(
            Selected::parse(&document, item, Default::default(), Some(&stop)),
            Err(Error::Cancelled)
        ));
    }
    assert!(matches!(
        Selected::parse(
            &document,
            item,
            Options {
                maximum_work: baseline.parsed().work_used() - 1,
                ..Default::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert_eq!(document.literal_input(), text);
}
pub(super) fn mapped(text: &str, y: &str, scale: &str) -> String {
    let document = CompoundDocument::parse(text, DocumentOptions::default(), None).unwrap();
    let solid = document.represented_members()[0].item_entity();
    let representation = document
        .records()
        .values()
        .flat_map(|r| r.iter())
        .find(|e| e.name == "ADVANCED_BREP_SHAPE_REPRESENTATION")
        .unwrap();
    let Value::Reference(context) = representation.args[2] else {
        unreachable!()
    };
    let old = text
        .lines()
        .find(|line| line.contains("=ADVANCED_BREP_SHAPE_REPRESENTATION("))
        .unwrap();
    let changed = old.replacen(&format!("#{solid}"), "#900020", 1);
    let records = format!(
        "#900010=CARTESIAN_POINT('',(0.,0.,0.));\n#900011=AXIS2_PLACEMENT_3D('',#900010,$,$);\n#900012=SHAPE_REPRESENTATION('',(#{solid},#900011),#{context});\n#900013=REPRESENTATION_MAP(#900011,#900012);\n#900014=CARTESIAN_POINT('',(11.,-5.,3.));\n#900015=DIRECTION('',(0.,2.,0.));\n#900016=DIRECTION('',{y});\n#900017=DIRECTION('',(0.,0.,5.));\n#900018=(CARTESIAN_TRANSFORMATION_OPERATOR(#900015,#900016,#900014,{scale})CARTESIAN_TRANSFORMATION_OPERATOR_3D(#900017)FUNCTIONALLY_DEFINED_TRANSFORMATION('','')GEOMETRIC_REPRESENTATION_ITEM()REPRESENTATION_ITEM(''));\n#900020=MAPPED_ITEM('',#900013,#900018);\n"
    );
    text.replacen(old, &changed, 1).replacen(
        "ENDSEC;\nEND-ISO-10303-21;",
        &format!("{records}ENDSEC;\nEND-ISO-10303-21;"),
        1,
    )
}
pub(super) fn contains_point(bounds: [[f64; 2]; 3], point: [f64; 3]) {
    for (bound, point) in bounds.into_iter().zip(point) {
        assert!(
            bound[0] <= point && point <= bound[1] && bound[1] - bound[0] < 1e-10,
            "{bound:?} {point}"
        );
    }
}
#[test]
fn actual_selected_step_mapping_retains_express_semantics_inverse_and_provenance() {
    let local = super::super::tests::fixture(false);
    for mapped_source in [false, true] {
        let text = if mapped_source {
            mapped(&local, "(-3.,0.,0.)", "1.")
        } else {
            local.clone()
        };
        let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
        let selected = Selected::parse(
            &document,
            document.represented_members()[0].item_entity(),
            Default::default(),
            None,
        )
        .unwrap();
        let placement = selected.placement(Default::default(), None).unwrap();
        assert!(std::ptr::eq(placement.selected(), &selected));
        assert_eq!(placement.mapping_link().is_some(), mapped_source);
        assert!(placement.work_used() > selected.parsed().work_used());
        assert!(placement.query_budget_bytes() > selected.parsed().query_budget_bytes());
        assert!(placement.determinant_bounds()[0] <= 1. && placement.determinant_bounds()[1] >= 1.);
        assert!(
            placement
                .linear_columns_bounds()
                .iter()
                .flatten()
                .all(|i| i[0] <= i[1])
        );
        assert!(
            placement
                .inverse_columns_bounds()
                .iter()
                .flatten()
                .all(|i| i[0] <= i[1])
        );
        let point = [1., 2., 3.];
        let world = if mapped_source { [9., -4., 6.] } else { point };
        contains_point(placement.map_box_mm(point.map(|x| [x, x])).unwrap(), world);
        contains_point(
            placement.inverse_box_mm(world.map(|x| [x, x])).unwrap(),
            point,
        );
        contains_point(
            placement.origin_mm(),
            if mapped_source {
                [11., -5., 3.]
            } else {
                [0.; 3]
            },
        );
        contains_point(
            placement.inverse_origin_mm(),
            if mapped_source {
                [5., 11., -3.]
            } else {
                [0.; 3]
            },
        );
        assert!(matches!(
            placement.map_box_mm([[f64::NAN, 0.]; 3]),
            Err(Error::Invalid)
        ));
        assert!(matches!(
            placement.inverse_box_mm([[0., f64::INFINITY]; 3]),
            Err(Error::Invalid)
        ));
    }
    for (y, scale) in [("(3.,0.,0.)", "1."), ("(-3.,0.,0.)", "0.")] {
        let text = mapped(&local, y, scale);
        let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
        let selected = Selected::parse(
            &document,
            document.represented_members()[0].item_entity(),
            Default::default(),
            None,
        )
        .unwrap();
        assert!(matches!(
            selected.placement(Default::default(), None),
            Err(Error::Unsupported)
        ));
    }
}
#[test]
fn actual_selected_mapping_refuses_real_schema_numeric_and_same_ledger_frontiers() {
    let local = super::super::tests::fixture(false);
    let text = mapped(&local, "(-3.,0.,0.)", "1.");
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let selected = Selected::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap();
    for options in [
        Options {
            maximum_work: 0,
            ..Default::default()
        },
        Options {
            maximum_bytes: 0,
            ..Default::default()
        },
        Options {
            maximum_work: selected.parsed.work,
            ..Default::default()
        },
        Options {
            maximum_work: selected.parsed.work + 1,
            ..Default::default()
        },
        Options {
            maximum_bytes: 1,
            ..Default::default()
        },
        Options {
            maximum_bytes: selected.parsed.bytes,
            ..Default::default()
        },
        Options {
            maximum_bytes: selected.parsed.bytes + size_of::<Placement<'_, '_, '_>>(),
            ..Default::default()
        },
        Options {
            maximum_bytes: selected.parsed.bytes + size_of::<Placement<'_, '_, '_>>() + 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            selected.placement(options, None),
            Err(Error::Invalid | Error::Limit)
        ));
    }
    let baseline = selected.placement(Default::default(), None).unwrap();
    assert!(matches!(
        selected.placement(
            Options {
                maximum_work: baseline.work_used() - 1,
                ..Default::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        baseline.map_box_mm([[f64::MAX, f64::MAX]; 3]),
        Err(Error::Unsupported)
    ));
    assert!(matches!(
        baseline.inverse_box_mm([[f64::MAX, f64::MAX]; 3]),
        Err(Error::Unsupported)
    ));
    let calls = std::cell::Cell::new(0usize);
    selected
        .placement(
            Default::default(),
            Some(&|| {
                calls.set(calls.get() + 1);
                false
            }),
        )
        .unwrap();
    for checkpoint in [1, 2, calls.get() / 2, calls.get() - 1] {
        let calls = std::cell::Cell::new(0usize);
        let stop = || {
            calls.set(calls.get() + 1);
            calls.get() >= checkpoint
        };
        assert!(matches!(
            selected.placement(Default::default(), Some(&stop)),
            Err(Error::Cancelled)
        ));
    }
    let wrong = text.replace("#900016,#900014,1.)", "#900016,3.,1.)");
    assert_ne!(wrong, text);
    let document = CompoundDocument::parse(&wrong, DocumentOptions::default(), None).unwrap();
    let selected = Selected::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap();
    assert!(matches!(
        selected.placement(Default::default(), None),
        Err(Error::Invalid)
    ));
    for scale in ["1.E308", "1.E-308"] {
        let text = mapped(&local, "(-3.,0.,0.)", scale);
        let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
        let selected = Selected::parse(
            &document,
            document.represented_members()[0].item_entity(),
            Default::default(),
            None,
        )
        .unwrap();
        assert!(matches!(
            selected.placement(Default::default(), None),
            Err(Error::Unsupported)
        ));
    }
}
#[test]
fn singleton_queries_isolate_unrelated_records_and_report_additional_walk_work() {
    for quarter in [false, true] {
        let text = super::super::tests::fixture(quarter);
        let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
        let parsed = Parsed::parse(
            &document,
            document.represented_members()[0].item_entity(),
            Default::default(),
            None,
        )
        .unwrap();
        let baseline = (
            parsed.work_used(),
            parsed.query_budget_bytes(),
            parsed.points().len(),
        );
        eprintln!(
            "NPT SINGLETON SELECTED quarter={quarter} work={} querybytes={} points={} documentbytes={}",
            baseline.0,
            baseline.1,
            baseline.2,
            document.retained_budget_bytes()
        );
        if !quarter {
            use std::fmt::Write;
            let mut extra = String::new();
            for id in 1_000_000..1_004_000 {
                writeln!(extra, "#{id}=CARTESIAN_POINT('',(1.,2.,3.));").unwrap();
            }
            let large = text.replacen(
                "ENDSEC;\nEND-ISO-10303-21;",
                &format!("{extra}ENDSEC;\nEND-ISO-10303-21;"),
                1,
            );
            if let Ok(path) = std::env::var("MULTICAD_NPT_READER_BENCHMARK_OUTPUT") {
                std::fs::write(path, &large).unwrap();
            }
            let document =
                CompoundDocument::parse(&large, DocumentOptions::default(), None).unwrap();
            let parsed = Parsed::parse(
                &document,
                document.represented_members()[0].item_entity(),
                Default::default(),
                None,
            )
            .unwrap();
            eprintln!(
                "NPT SINGLETON PLUS4K work={} querybytes={} points={} documentbytes={}",
                parsed.work_used(),
                parsed.query_budget_bytes(),
                parsed.points().len(),
                document.retained_budget_bytes()
            );
            assert_eq!(
                (
                    parsed.work_used(),
                    parsed.query_budget_bytes(),
                    parsed.points().len()
                ),
                baseline
            );
            let outside = text.replacen(
                "ENDSEC;\nEND-ISO-10303-21;",
                "#1000000=CARTESIAN_POINT('',(99999999999.,0.,0.));\nENDSEC;\nEND-ISO-10303-21;",
                1,
            );
            let document =
                CompoundDocument::parse(&outside, DocumentOptions::default(), None).unwrap();
            assert!(
                Parsed::parse(
                    &document,
                    document.represented_members()[0].item_entity(),
                    Default::default(),
                    None
                )
                .is_ok()
            );
        }
    }
}
