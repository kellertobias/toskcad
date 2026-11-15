//! Defensive readers consume authentic neutral records, never forged receipts.
use super::*;
use crate::step_read::npt_stock::tests::{change_record, fields, fixture, record_id, with_fields};

#[test]
fn original_pcurve_context_schema_and_real_read_frontiers_are_validated() {
    let text = fixture(false);
    let records = part21::parse(&text).unwrap();
    let id = record_id(&text, "=PCURVE(", 0);
    let Value::Reference(rep) = records[&id][0].args[2] else {
        panic!()
    };
    let Value::Reference(context) = records[&rep][0].args[2] else {
        panic!()
    };
    let mut baseline = Budget::new(Default::default(), None).unwrap();
    let expected = pcurve(&records, id, &mut baseline).unwrap();
    let changed = [
        (
            change_record(&text, rep, |line| {
                let mut args = fields(line);
                args[1] = "()";
                with_fields(line, &args)
            }),
            Error::Unsupported,
        ),
        (
            change_record(&text, context, |line| {
                line.replace("PARAMETRIC_REPRESENTATION_CONTEXT()", "")
            }),
            Error::Unsupported,
        ),
        (
            change_record(&text, context, |line| {
                line.replace(
                    "PARAMETRIC_REPRESENTATION_CONTEXT",
                    "FOREIGN_REPRESENTATION_CONTEXT",
                )
            }),
            Error::Unsupported,
        ),
        (
            change_record(&text, context, |line| {
                line.replace(
                    "GEOMETRIC_REPRESENTATION_CONTEXT(2)",
                    "GEOMETRIC_REPRESENTATION_CONTEXT(2,2)",
                )
            }),
            Error::Invalid,
        ),
        (
            change_record(&text, context, |line| {
                line.replace(
                    "GEOMETRIC_REPRESENTATION_CONTEXT(2)",
                    "GEOMETRIC_REPRESENTATION_CONTEXT(3)",
                )
            }),
            Error::Unsupported,
        ),
    ];
    for (changed, expected_error) in changed {
        let actual_records = part21::parse(&changed).unwrap();
        let mut budget = Budget::new(Default::default(), None).unwrap();
        assert_eq!(
            pcurve(&actual_records, id, &mut budget),
            Err(expected_error)
        );
    }
    // The first entity lookup costs32; the actual representation lookup must
    // honor the same caller budget rather than starting an independent one.
    let mut tight = Budget::new(
        Options {
            maximum_work: 63,
            ..Default::default()
        },
        None,
    )
    .unwrap();
    assert_eq!(pcurve(&records, id, &mut tight), Err(Error::Limit));
    let polls = std::cell::Cell::new(0usize);
    let cancel = || {
        polls.set(polls.get() + 1);
        polls.get() == 2
    };
    let mut cancelled = Budget::new(Default::default(), Some(&cancel)).unwrap();
    assert_eq!(pcurve(&records, id, &mut cancelled), Err(Error::Cancelled));
    assert_eq!(polls.get(), 2);
    let mut recovered = Budget::new(Default::default(), None).unwrap();
    assert_eq!(pcurve(&records, id, &mut recovered).unwrap(), expected);
    assert_eq!(recovered.used, baseline.used);
}

#[test]
fn original_cap_curve_schema_decimal_binding_and_height_are_checked() {
    let text = fixture(false);
    let document = CompoundDocument::parse(&text, Default::default(), None).unwrap();
    let parsed = Parsed::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap();
    let face = parsed
        .faces
        .iter()
        .find(|f| matches!(f.support, Support::Plane { .. }))
        .unwrap();
    let Support::Plane { height } = face.support else {
        unreachable!()
    };
    let edge = &parsed.edges[face.coedges[0].edge];
    let curve_id = parsed.curves[edge.curve].id;
    let height_id = parsed.points[height].id;
    let records = document.records();
    let Value::Reference(support) = records[&face.id][0].args[2] else {
        panic!()
    };
    let mut query = Budget::new(Default::default(), None).unwrap();
    let geometry = edge
        .pcurves
        .iter()
        .find_map(|id| {
            let (surface, curve) = pcurve(records, *id, &mut query).unwrap();
            (surface == support).then_some(curve)
        })
        .unwrap();
    let read = |input: &str| {
        let records = part21::parse(input).unwrap();
        let mut budget = Budget::new(Default::default(), None).unwrap();
        let (points, _) = super::super::selection::points(
            &records,
            document.represented_members()[0].item_entity(),
            &mut budget,
        )
        .unwrap();
        let curve = super::super::curve(&records, &points, curve_id, &mut budget).unwrap();
        let height = point_index(&points, height_id, &mut budget).unwrap();
        cap_curve(&records, &points, geometry, &curve, height, &mut budget)
    };
    assert_eq!(read(&text), Ok(()));
    let changed = [
        (
            change_record(&text, geometry, |line| {
                let mut args = fields(line);
                args[0] = "0";
                with_fields(line, &args)
            }),
            Error::Unsupported,
        ),
        (
            change_record(&text, geometry, |line| {
                let mut args = fields(line);
                args[1] = "10";
                with_fields(line, &args)
            }),
            Error::Unsupported,
        ),
        (
            change_record(&text, geometry, |line| {
                let mut args = fields(line);
                args[2] = "()";
                with_fields(line, &args)
            }),
            Error::Invalid,
        ),
    ];
    for (input, expected) in changed {
        assert_eq!(read(&input), Err(expected));
    }
    let Value::List(poles) = &records[&geometry][0].args[2] else {
        panic!()
    };
    let Value::Reference(uv_point) = poles[0] else {
        panic!()
    };
    let changed = change_record(&text, uv_point, |line| {
        let mut args = fields(line);
        args[1] = "(123.,456.)";
        with_fields(line, &args)
    });
    assert_eq!(read(&changed), Err(Error::Unsupported));
    let wrong_dimension = change_record(&text, uv_point, |line| {
        let mut args = fields(line);
        args[1] = "(0.,0.,0.)";
        with_fields(line, &args)
    });
    assert_eq!(read(&wrong_dimension), Err(Error::Unsupported));
    let point_id = parsed.points[parsed.curves[edge.curve].poles[0]].id;
    let changed = change_record(&text, point_id, |line| {
        let mut args = fields(line);
        let Value::List(xyz) = &records[&point_id][0].args[1] else {
            panic!()
        };
        let Value::Number(z) = xyz[2] else { panic!() };
        let list = args[1]
            .trim_matches(['(', ')'])
            .split(',')
            .collect::<Vec<_>>();
        let xyz = format!("({},{},{:.16E})", list[0], list[1], z + 0.001);
        args[1] = &xyz;
        with_fields(line, &args)
    });
    assert_eq!(read(&changed), Err(Error::Unsupported));
    let curve = &parsed.curves[edge.curve];
    let polls = std::cell::Cell::new(0usize);
    let count = || {
        polls.set(polls.get() + 1);
        false
    };
    let mut baseline = Budget::new(Default::default(), Some(&count)).unwrap();
    cap_curve(
        records,
        &parsed.points,
        geometry,
        curve,
        height,
        &mut baseline,
    )
    .unwrap();
    let mut tight = Budget::new(
        Options {
            maximum_work: baseline.used - 1,
            ..Default::default()
        },
        None,
    )
    .unwrap();
    assert_eq!(
        cap_curve(records, &parsed.points, geometry, curve, height, &mut tight),
        Err(Error::Limit)
    );
    let seen = std::cell::Cell::new(0usize);
    let cancel = || {
        seen.set(seen.get() + 1);
        seen.get() == polls.get()
    };
    let mut cancelled = Budget::new(Default::default(), Some(&cancel)).unwrap();
    assert_eq!(
        cap_curve(
            records,
            &parsed.points,
            geometry,
            curve,
            height,
            &mut cancelled
        ),
        Err(Error::Cancelled)
    );
    assert_eq!(seen.get(), polls.get());
    assert_eq!(read(&text), Ok(()));
    assert_eq!(document.literal_input(), text);
}
