use super::*;
use crate::step_read::owned_insert::DocumentOptions;
fn source(mapped: bool) -> OwnedMaterial {
    source_size(false, mapped)
}
fn source_size(quarter: bool, mapped: bool) -> OwnedMaterial {
    let text = super::super::super::super::tests::fixture(quarter);
    let text = if mapped {
        super::super::super::tests::mapped(&text, "(-3.,0.,0.)", "1.")
    } else {
        text
    };
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    Selected::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap()
    .retain(Default::default(), None)
    .unwrap()
}
fn surface(source: &OwnedMaterial) -> [f64; 3] {
    source
        .tensors()
        .iter()
        .enumerate()
        .map(|(index, _)| {
            source
                .cell(index, [0, 0])
                .unwrap()
                .point_bounds([[0.5, 0.5]; 2])
                .unwrap()
                .map(|b| b[0].midpoint(b[1]))
        })
        .max_by(|a, b| (a[0] * a[0] + a[1] * a[1]).total_cmp(&(b[0] * b[0] + b[1] * b[1])))
        .unwrap()
}
#[test]
fn actual_literal_boxes_use_original_surfaces_and_actual_step_mapping() {
    for (quarter, mapped) in [(false, false), (false, true), (true, false), (true, true)] {
        let source = source_size(quarter, mapped);
        let p = surface(&source);
        for (point, expected) in [
            ([0., 0., p[2]], Location::Inside),
            ([p[0] * 0.98, p[1] * 0.98, p[2]], Location::Inside),
            ([p[0] * 1.02, p[1] * 1.02, p[2]], Location::Outside),
            ([100., 0., p[2]], Location::Outside),
            ([0., 0., -1.], Location::Outside),
            ([0., 0., 0.], Location::Unresolved),
        ] {
            let bounds = point.map(|v| [v, v]);
            let local = source
                .classify_local_box(bounds, Default::default(), None)
                .unwrap();
            assert_eq!(local.location(), expected, "{point:?}");
            assert!(std::ptr::eq(local.source(), &source));
            assert_eq!(local.local_bounds_mm(), bounds);
            assert!(
                local.work_used() > 0
                    && local.query_budget_bytes() > source.retained_budget_bytes()
            );
            let world = source.map_box_mm(bounds).unwrap();
            assert_eq!(
                source
                    .classify_box(world, Default::default(), None)
                    .unwrap()
                    .location(),
                expected
            );
            eprintln!(
                "NPT LITERAL QUERY mapped={mapped} point={point:?} {expected:?} work={} nodes={}",
                local.work_used(),
                local.nodes()
            );
        }
        assert_eq!(
            source
                .classify_local_box(p.map(|v| [v, v]), Default::default(), None)
                .unwrap()
                .location(),
            Location::Unresolved
        );
    }
}
#[test]
fn actual_query_frontiers_preserve_unresolved_and_resource_errors() {
    let source = source(false);
    let p = surface(&source).map(|v| [v, v]);
    let standard = QueryOptions::default();
    for options in [
        QueryOptions {
            maximum_work: 0,
            ..standard
        },
        QueryOptions {
            maximum_work: 16_000_001,
            ..standard
        },
        QueryOptions {
            maximum_bytes: 0,
            ..standard
        },
        QueryOptions {
            maximum_bytes: 64 * 1024 * 1024 + 1,
            ..standard
        },
        QueryOptions {
            maximum_nodes: 65_537,
            ..standard
        },
        QueryOptions {
            maximum_depth: 65,
            ..standard
        },
    ] {
        assert!(matches!(
            source.classify_local_box(p, options, None),
            Err(Error::Invalid)
        ));
    }
    for options in [
        QueryOptions {
            maximum_work: 1,
            ..standard
        },
        QueryOptions {
            maximum_bytes: 1,
            ..standard
        },
        QueryOptions {
            maximum_nodes: 0,
            ..standard
        },
    ] {
        assert!(matches!(
            source.classify_local_box(p, options, None),
            Err(Error::Limit)
        ));
    }
    assert_eq!(
        source
            .classify_local_box(
                p,
                QueryOptions {
                    maximum_depth: 0,
                    ..standard
                },
                None
            )
            .unwrap()
            .location(),
        Location::Unresolved
    );
    assert!(matches!(
        source.classify_local_box([[f64::NAN, 0.]; 3], standard, None),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        source.classify_box(p, standard, Some(&|| true)),
        Err(Error::Cancelled)
    ));
    let wide = [[-1e308, 1e308], [-1e308, 1e308], p[2]];
    assert_eq!(
        source
            .classify_local_box(wide, standard, None)
            .unwrap()
            .location(),
        Location::Unresolved
    );
}
#[test]
fn changed_actual_decimal_xy_geometry_overrides_any_nominal_thread_recipe() {
    use std::collections::{BTreeMap, BTreeSet};
    let text = super::super::super::super::tests::fixture(false);
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let mut cap_points = BTreeSet::new();
    for values in document.records().values() {
        if let [entity] = values.as_slice()
            && entity.name == "B_SPLINE_CURVE_WITH_KNOTS"
            && matches!(entity.args.get(1),Some(Value::Number(n)) if *n==11.)
            && let Value::List(poles) = &entity.args[2]
        {
            for pole in poles {
                if let Value::Reference(id) = pole {
                    cap_points.insert(*id);
                }
            }
        }
    }
    let twice = |text: &str| {
        let (mantissa, exponent) = text.split_once('E').unwrap_or((text, "0"));
        let fraction = mantissa.split_once('.').map_or(0, |(_, s)| s.len());
        let digits = mantissa
            .replace('.', "")
            .parse::<num_bigint::BigInt>()
            .unwrap()
            * 2;
        format!(
            "{digits}E{}",
            exponent.parse::<i32>().unwrap() - fraction as i32
        )
    };
    let mut replacements = BTreeMap::new();
    for (id, values) in document.records() {
        if let [entity] = values.as_slice()
            && entity.name == "CARTESIAN_POINT"
            && let Some(Value::List(coords)) = entity.args.get(1)
            && (coords.len() == 3 || cap_points.contains(id))
        {
            let numbers = part21::numeric_spellings(entity).unwrap();
            let mut result = vec![twice(numbers[0]), twice(numbers[1])];
            if coords.len() == 3 {
                result.push(numbers[2].to_string());
            }
            replacements.insert(
                *id,
                format!("#{id}=CARTESIAN_POINT('',({}));", result.join(",")),
            );
        }
    }
    let scaled = text
        .lines()
        .map(|line| {
            line.strip_prefix('#')
                .and_then(|s| s.split_once('='))
                .and_then(|(id, _)| id.parse::<u32>().ok())
                .and_then(|id| replacements.get(&id))
                .map_or(line, String::as_str)
        })
        .collect::<Vec<_>>()
        .join("\n");
    let original = Selected::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap()
    .retain(Default::default(), None)
    .unwrap();
    let document = CompoundDocument::parse(&scaled, DocumentOptions::default(), None).unwrap();
    let changed = Selected::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap()
    .retain(Default::default(), None)
    .unwrap();
    let height = original.local_enclosure_mm()[2];
    let point = [5.8, 0., height[0].midpoint(height[1])].map(|v| [v, v]);
    assert_eq!(
        original
            .classify_local_box(point, Default::default(), None)
            .unwrap()
            .location(),
        Location::Outside
    );
    assert_eq!(
        changed
            .classify_local_box(point, Default::default(), None)
            .unwrap()
            .location(),
        Location::Inside
    );
}

#[test]
fn whole_query_boxes_and_exact_resource_frontiers_remain_conservative() {
    use std::cell::Cell;
    let source = source(false);
    let surface = surface(&source);
    for (factor, expected) in [
        (0.98, Location::Inside),
        (1.02, Location::Outside),
        (1.00001, Location::Outside),
    ] {
        let point = [factor * surface[0], factor * surface[1], surface[2]];
        let input = point.map(|v| [v - 1e-5, v + 1e-5]);
        let calls = Cell::new(0);
        let count = || {
            calls.set(calls.get() + 1);
            false
        };
        let receipt = source
            .classify_local_box(input, Default::default(), Some(&count))
            .unwrap();
        assert_eq!(receipt.location(), expected);
        let work = receipt.work_used();
        let bytes = receipt.query_budget_bytes();
        let nodes = receipt.nodes();
        assert_eq!(
            source
                .classify_local_box(
                    input,
                    QueryOptions {
                        maximum_work: work,
                        maximum_bytes: bytes,
                        maximum_nodes: nodes,
                        ..Default::default()
                    },
                    None
                )
                .unwrap()
                .location(),
            expected
        );
        for options in [
            QueryOptions {
                maximum_work: work - 1,
                ..Default::default()
            },
            QueryOptions {
                maximum_bytes: bytes - 1,
                ..Default::default()
            },
            QueryOptions {
                maximum_nodes: nodes - 1,
                ..Default::default()
            },
        ] {
            assert!(matches!(
                source.classify_local_box(input, options, None),
                Err(Error::Limit)
            ));
        }
        for stop in [2, calls.get() / 2, calls.get()] {
            let cursor = Cell::new(0);
            let cancel = || {
                cursor.set(cursor.get() + 1);
                cursor.get() == stop
            };
            assert!(matches!(
                source.classify_local_box(input, Default::default(), Some(&cancel)),
                Err(Error::Cancelled)
            ));
        }
    }
    let crossed = [surface[0], surface[1], surface[2]].map(|v| [v - 1e-3, v + 1e-3]);
    assert_eq!(
        source
            .classify_local_box(
                crossed,
                QueryOptions {
                    maximum_depth: 6,
                    ..Default::default()
                },
                None
            )
            .unwrap()
            .location(),
        Location::Unresolved
    );
}

#[test]
fn original_subnormal_decimal_geometry_keeps_queries_conservative() {
    let text = super::super::super::super::tests::fixture(false);
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let scaled = super::super::super::super::material::tests::scale_xy(&text, &document, 9);
    let document = CompoundDocument::parse(&scaled, DocumentOptions::default(), None).unwrap();
    let source = Selected::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap()
    .retain(Default::default(), None)
    .unwrap();
    let z = surface(&source)[2];
    for point in [
        [0., 0., z],
        [1e-200, 0., z],
        [4e-161, 0., z],
        [4.8e-161, 0., z],
    ] {
        let result = source.classify_local_box(
            point.map(|v| [v, v]),
            QueryOptions {
                maximum_depth: 6,
                ..Default::default()
            },
            None,
        );
        eprintln!(
            "NPT TINY QUERY {point:?}: {:?}",
            result.as_ref().map(|r| r.location())
        );
        assert!(result.is_ok() || matches!(result, Err(Error::Limit)));
    }
}
#[test]
fn broad_actual_boxes_do_not_turn_ambiguous_rays_into_material_admission() {
    let source = source(false);
    let z = surface(&source)[2];
    for input in [
        [[-5.9, 5.9], [4.8, 4.9], [z, z]],
        [[-5.9, 5.9], [5., 5.], [z, z]],
        [[-5.9, 5.9], [5.2, 5.2], [z, z]],
        [[4., 5.], [4., 5.], [z, z]],
        [[-6., 6.], [-6., 6.], [z, z]],
    ] {
        let result = source.classify_local_box(
            input,
            QueryOptions {
                maximum_depth: 6,
                ..Default::default()
            },
            None,
        );
        eprintln!(
            "NPT BROAD QUERY {input:?}: {:?}",
            result.as_ref().map(|r| r.location())
        );
        assert!(
            matches!(result,Ok(ref r) if r.location()==Location::Unresolved || r.location()==Location::Outside)
                || matches!(result, Err(Error::Limit))
        );
    }
}

#[test]
fn non_dyadic_actual_surface_rays_have_independent_strict_inside_outside_queries() {
    let source = source(false);
    let point = source
        .cell(0, [0, 0])
        .unwrap()
        .point_bounds([[0.371, 0.371], [0.417, 0.417]])
        .unwrap()
        .map(|b| b[0].midpoint(b[1]));
    for (factor, expected) in [(0.9999, Location::Inside), (1.0001, Location::Outside)] {
        let input = [factor * point[0], factor * point[1], point[2]].map(|v| [v, v]);
        let result = source
            .classify_local_box(input, Default::default(), None)
            .unwrap();
        eprintln!(
            "NPT NONDYADIC QUERY {factor}: {:?}, {} work {} nodes",
            result.location(),
            result.work_used(),
            result.nodes()
        );
        assert_eq!(result.location(), expected);
    }
}
