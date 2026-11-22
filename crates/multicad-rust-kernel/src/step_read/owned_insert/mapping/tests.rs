use super::super::tests::file;
use super::*;
const OP: u32 = 2_000_000;
fn fixture(axes: [&str; 3], scale: &str) -> String {
    let mut text = file(false);
    let at = text.rfind("ENDSEC;").unwrap();
    text.insert_str(
        at,
        &format!(
            "#2000001=CARTESIAN_POINT('',(11.,-5.,3.));\n\
         #2000002=DIRECTION('',{});\n\
         #2000003=DIRECTION('',{});\n\
         #2000004=DIRECTION('',{});\n\
         #2000000=(CARTESIAN_TRANSFORMATION_OPERATOR(#2000002,#2000003,#2000001,{scale})\
         CARTESIAN_TRANSFORMATION_OPERATOR_3D(#2000004)\
         FUNCTIONALLY_DEFINED_TRANSFORMATION('','')\
         GEOMETRIC_REPRESENTATION_ITEM() REPRESENTATION_ITEM(''));\n",
            axes[0], axes[1], axes[2]
        ),
    );
    text
}
fn contains(actual: [[f64; 2]; 3], expected: [f64; 3]) {
    for (bound, x) in actual.into_iter().zip(expected) {
        assert!(
            bound[0] <= x && x <= bound[1],
            "{bound:?} does not enclose{x}"
        );
        assert!(bound[1] - bound[0] < 1e-10);
    }
}
#[test]
fn actual_projected_directions_preserve_reflections_and_ignore_ratio_magnitude() {
    for (axes, expected) in [
        (["(2.,0.,2.)", "(2.,3.,5.)", "(0.,0.,4.)"], [12., -3., 6.]),
        (["(2.,0.,2.)", "(2.,-3.,5.)", "(0.,0.,4.)"], [12., -7., 6.]),
        (
            ["(1.,0.,0.)", "(0.,4.,-3.)", "(0.,3.,4.)"],
            [12., -1.6, 4.2],
        ),
        (
            [
                "(2.E300,0.,2.E300)",
                "(2.E300,3.E300,5.E300)",
                "(0.,0.,4.E300)",
            ],
            [12., -3., 6.],
        ),
        (
            [
                "(2.E-300,0.,2.E-300)",
                "(2.E-300,3.E-300,5.E-300)",
                "(0.,0.,4.E-300)",
            ],
            [12., -3., 6.],
        ),
    ] {
        let text = fixture(axes, "$");
        let p = parse(&text, Options::default(), None).unwrap();
        let map = p.physical_mapping(OP, Options::default(), None).unwrap();
        assert!(std::ptr::eq(map.parsed(), &p));
        assert_eq!(map.operator_entity(), OP);
        contains(map.origin_mm(), [11., -5., 3.]);
        contains(
            map.map_box_mm([[1., 1.], [2., 2.], [3., 3.]]).unwrap(),
            expected,
        );
        let volume = map.volume_scale_bounds().unwrap();
        assert!(volume[0] <= 1. && volume[1] >= 1.);
        assert_eq!(map.scale_bounds(), [1., 1.]);
        let axes = map.axes();
        assert!(axes[0][0][0] > 0.99999999999);
        assert_eq!(
            map.handedness().unwrap(),
            if expected[1] == -7. { -1 } else { 1 }
        );
        assert_eq!(p.literal_input(), text);
    }
}
#[test]
fn actual_opposite_handed_operators_do_not_have_one_sign_for_their_shared_interval_box() {
    let positive = fixture(["(1.,0.,0.)", "(0.,1.,0.)", "(0.,0.,1.)"], "$");
    let negative = fixture(["(1.,0.,0.)", "(0.,-1.,0.)", "(0.,0.,1.)"], "$");
    let p = parse(&positive, Options::default(), None).unwrap();
    let n = parse(&negative, Options::default(), None).unwrap();
    let p = p.physical_mapping(OP, Options::default(), None).unwrap();
    let n = n.physical_mapping(OP, Options::default(), None).unwrap();
    assert_eq!(p.handedness().unwrap(), 1);
    assert_eq!(n.handedness().unwrap(), -1);
    // A real enclosure of these two distinct decoded operator bases cannot
    // certify one sign. Exercise the numerical sign primitive, without
    // constructing a fictitious Parsed packet or PhysicalMapping receipt.
    let family = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            I::new(
                p.axes[i][j].lo.min(n.axes[i][j].lo),
                p.axes[i][j].hi.max(n.axes[i][j].hi),
            )
        })
    });
    assert!(matches!(basis_handedness(family), Err(Error::Unsupported)));
}
#[test]
fn actual_optional_axes_and_uniform_scale_follow_express_defaults() {
    let original = fixture(["(1.,0.,0.)", "(0.,1.,0.)", "(0.,0.,1.)"], "2.");
    let text = original
        .replace(
            "OPERATOR(#2000002,#2000003,#2000001,",
            "OPERATOR($,$,#2000001,",
        )
        .replace("OPERATOR_3D(#2000004)", "OPERATOR_3D($)");
    let p = parse(&text, Options::default(), None).unwrap();
    let map = p.physical_mapping(OP, Options::default(), None).unwrap();
    let box_bounds = map.map_box_mm([[-2., 3.], [1., 4.], [0., 2.]]).unwrap();
    for (actual, expected) in box_bounds.into_iter().zip([[7., 17.], [-3., 3.], [3., 7.]]) {
        assert!(actual[0] <= expected[0] && actual[1] >= expected[1]);
    }
    let volume = map.volume_scale_bounds().unwrap();
    assert!(volume[0] <= 8. && volume[1] >= 8.);
    // Exactly +/-X for axis3 chooses the alternate default first axis +Y.
    for sign in ["1.", "-1."] {
        let text = fixture(
            ["(1.,0.,0.)", "(0.,0.,1.)", &format!("({sign},0.,0.)")],
            "$",
        )
        .replace("OPERATOR(#2000002,", "OPERATOR($,");
        let p = parse(&text, Options::default(), None).unwrap();
        let map = p.physical_mapping(OP, Options::default(), None).unwrap();
        let expected = if sign == "1." {
            [14., -4., 5.]
        } else {
            [8., -4., 5.]
        };
        contains(
            map.map_box_mm([[1., 1.], [2., 2.], [3., 3.]]).unwrap(),
            expected,
        );
    }
}
#[test]
fn degenerate_uncertain_and_unrepresentable_physical_queries_refuse() {
    for axes in [
        ["(0.,0.,0.)", "(0.,1.,0.)", "(0.,0.,1.)"],
        ["(0.,0.,1.)", "(0.,1.,0.)", "(0.,0.,1.)"],
        ["(1.,0.,0.)", "(1.,0.,0.)", "(0.,0.,1.)"],
        ["(1.,0.,0.)", "(0.,1.,0.)", "(0.,0.,1.E-1000)"],
    ] {
        let text = fixture(axes, "$");
        let p = parse(&text, Options::default(), None).unwrap();
        assert!(matches!(
            p.physical_mapping(OP, Options::default(), None),
            Err(Error::Unsupported)
        ));
    }
    let text = fixture(["(1.,0.,0.)", "(0.,1.,0.)", "(0.,0.,1.)"], "1.E200");
    let p = parse(&text, Options::default(), None).unwrap();
    let map = p.physical_mapping(OP, Options::default(), None).unwrap();
    assert!(matches!(map.volume_scale_bounds(), Err(Error::Unsupported)));
    assert!(matches!(
        map.map_box_mm([[1e200, 1e200]; 3]),
        Err(Error::Unsupported)
    ));
    assert!(matches!(
        map.map_box_mm([[f64::NAN, 1.]; 3]),
        Err(Error::Invalid)
    ));
    assert!(matches!(map.map_box_mm([[2., 1.]; 3]), Err(Error::Invalid)));
    for scale in ["0.", "-1.", "1.E-1000"] {
        let text = fixture(["(1.,0.,0.)", "(0.,1.,0.)", "(0.,0.,1.)"], scale);
        let p = parse(&text, Options::default(), None).unwrap();
        assert!(matches!(
            p.physical_mapping(OP, Options::default(), None),
            Err(Error::Unsupported)
        ));
    }
}
#[test]
fn actual_operator_schema_is_checked_without_trusting_unused_metadata() {
    let text = fixture(["(1.,0.,0.)", "(0.,1.,0.)", "(0.,0.,1.)"], "$");
    for (from, to, expected) in [
        (
            "CARTESIAN_TRANSFORMATION_OPERATOR(",
            "UNSUPPORTED_OPERATOR(",
            Error::Unsupported,
        ),
        (
            "CARTESIAN_TRANSFORMATION_OPERATOR_3D(",
            "UNSUPPORTED_DIMENSION(",
            Error::Unsupported,
        ),
        (
            "FUNCTIONALLY_DEFINED_TRANSFORMATION(",
            "UNSUPPORTED_FUNCTION(",
            Error::Unsupported,
        ),
        (
            "GEOMETRIC_REPRESENTATION_ITEM()",
            "UNSUPPORTED_GEOMETRY()",
            Error::Unsupported,
        ),
        (
            " REPRESENTATION_ITEM('')",
            " UNSUPPORTED_ITEM('')",
            Error::Unsupported,
        ),
        (
            " REPRESENTATION_ITEM('')",
            " REPRESENTATION_ITEM($)",
            Error::Unsupported,
        ),
        (
            "GEOMETRIC_REPRESENTATION_ITEM()",
            "GEOMETRIC_REPRESENTATION_ITEM(1.)",
            Error::Unsupported,
        ),
        (
            "FUNCTIONALLY_DEFINED_TRANSFORMATION('','')",
            "FUNCTIONALLY_DEFINED_TRANSFORMATION('')",
            Error::Unsupported,
        ),
        ("#2000001,$)", "#2000001)", Error::Unsupported),
        ("OPERATOR_3D(#2000004)", "OPERATOR_3D()", Error::Unsupported),
        (
            " REPRESENTATION_ITEM('')",
            " EXTRA() REPRESENTATION_ITEM('')",
            Error::Unsupported,
        ),
        (
            "#2000002,#2000003,#2000001",
            "#2000002,#2000003,#4294967295",
            Error::Invalid,
        ),
        (
            "#2000002=DIRECTION('',(1.,0.,0.))",
            "#2000002=DIRECTION($,(1.,0.,0.))",
            Error::Invalid,
        ),
        (
            "#2000002=DIRECTION('',(1.,0.,0.))",
            "#2000002=DIRECTION('',(1.,0.))",
            Error::Unsupported,
        ),
        (
            "#2000002=DIRECTION('',(1.,0.,0.))",
            "#2000002=DIRECTION('',$)",
            Error::Invalid,
        ),
        (
            "#2000002=DIRECTION('',(1.,0.,0.))",
            "#2000002=DIRECTION('',(1.7976931348623157E308,0.,0.))",
            Error::Unsupported,
        ),
        ("#2000001,$)", "#2000001,#2000002)", Error::Invalid),
    ] {
        assert!(text.contains(from), "missingmutation{from}");
        let changed = text.replacen(from, to, 1);
        let p = parse(&changed, Options::default(), None).unwrap();
        assert_eq!(
            p.physical_mapping(OP, Options::default(), None).err(),
            Some(expected),
            "{from}→{to}"
        );
    }
    let p = parse(&text, Options::default(), None).unwrap();
    assert!(matches!(
        p.physical_mapping(u32::MAX, Options::default(), None),
        Err(Error::Invalid)
    ));
}
#[test]
fn exact_operator_resource_frontiers_and_every_cancel_checkpoint_preserve_input() {
    use std::cell::Cell;
    let text = fixture(["(1.,0.,0.)", "(0.,1.,0.)", "(0.,0.,1.)"], "$");
    let p = parse(&text, Options::default(), None).unwrap();
    let map = p.physical_mapping(OP, Options::default(), None).unwrap();
    let exact = Options {
        maximum_work: map.work_used(),
        retained_budget_bytes: map.query_budget_bytes(),
    };
    assert_eq!(
        p.physical_mapping(OP, exact, None).unwrap().axes(),
        map.axes()
    );
    for options in [
        Options {
            maximum_work: exact.maximum_work - 1,
            ..exact
        },
        Options {
            maximum_work: 1,
            ..exact
        },
        Options {
            retained_budget_bytes: exact.retained_budget_bytes - 1,
            ..exact
        },
    ] {
        assert!(matches!(
            p.physical_mapping(OP, options, None),
            Err(Error::Limit)
        ));
    }
    for options in [
        Options {
            maximum_work: 0,
            ..exact
        },
        Options {
            maximum_work: 16_000_001,
            ..exact
        },
        Options {
            retained_budget_bytes: 0,
            ..exact
        },
        Options {
            retained_budget_bytes: 64 * 1024 * 1024 + 1,
            ..exact
        },
    ] {
        assert!(matches!(
            p.physical_mapping(OP, options, None),
            Err(Error::Invalid)
        ));
    }
    let calls = Cell::new(0);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    p.physical_mapping(OP, exact, Some(&count)).unwrap();
    let last = calls.get();
    for at in 1..=last {
        calls.set(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() == at
        };
        assert!(
            matches!(
                p.physical_mapping(OP, exact, Some(&cancel)),
                Err(Error::Cancelled)
            ),
            "checkpoint{at}"
        );
    }
    assert_eq!(map.parsed().literal_input(), text);
    assert_eq!(
        p.physical_mapping(OP, exact, None).unwrap().axes(),
        map.axes()
    );
}
#[test]
fn neutral_operator_entry_keeps_its_document_independent_resource_contract() {
    let text = fixture(["(0.,2.,0.)", "(-3.,0.,0.)", "(0.,0.,5.)"], "1.");
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let baseline =
        evaluate_neutral_operator(document.records(), OP, Default::default(), None).unwrap();
    assert_eq!(baseline.handedness, 1);
    assert!(baseline.work > 512 && baseline.bytes > SCRATCH);
    for options in [
        Options {
            maximum_work: 0,
            ..Default::default()
        },
        Options {
            maximum_work: 16_000_001,
            ..Default::default()
        },
        Options {
            retained_budget_bytes: 0,
            ..Default::default()
        },
        Options {
            retained_budget_bytes: 64 * 1024 * 1024 + 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            evaluate_neutral_operator(document.records(), OP, options, None),
            Err(Error::Invalid)
        ));
    }
    for options in [
        Options {
            maximum_work: baseline.work - 1,
            ..Default::default()
        },
        Options {
            retained_budget_bytes: baseline.bytes - 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            evaluate_neutral_operator(document.records(), OP, options, None),
            Err(Error::Limit)
        ));
    }
    assert!(matches!(
        evaluate_neutral_operator(document.records(), OP, Default::default(), Some(&|| true)),
        Err(Error::Cancelled)
    ));
}
