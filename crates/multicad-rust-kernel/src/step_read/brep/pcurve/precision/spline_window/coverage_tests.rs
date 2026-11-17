//! Genuine parsed-input contracts, including helper-only refusals. No sealed
//! CurveCells, Source fields, successful receipts or work counters are edited.
use super::*;
fn append(text: &str, record: &str) -> String {
    let at = text.rfind("ENDSEC;").unwrap();
    format!("{}\n{}\n{}", &text[..at], record, &text[at..])
}
#[test]
fn actual_metadata_numeric_nodes_and_typed_children_have_bounded_preflight() {
    let positive = input("0.0625", "1.");
    let typed = append(&positive, "#990001=WORKSPACE_AUDIT(WRAPPED(1.,(2.,3.)));");
    let source = Source::parse(&typed, 1e-7).unwrap();
    numeric_preflight(&source, 990001, &mut Budget::new(512_000_000, None)).unwrap();
    for (items, count) in [("0.", 2049), ("$", 8193)] {
        let extra = format!(
            "#990001=WORKSPACE_AUDIT(({}));",
            vec![items; count].join(",")
        );
        let text = append(&positive, &extra);
        let source = Source::parse(&text, 1e-7).unwrap();
        assert!(matches!(
            numeric_preflight(&source, 990001, &mut Budget::new(512_000_000, None)),
            Err(ReadError::Limit)
        ));
    }
    let source = Source::parse(&typed, 1e-7).unwrap();
    numeric_preflight(&source, 990001, &mut Budget::new(512_000_000, None)).unwrap();
}
fn actual_target(source: &Source<'_>) -> crate::brep::NurbsCurve {
    let Curve::Nurbs { curve, .. } = source.curve(900001, 0).unwrap() else {
        panic!("real parsed spline")
    };
    curve
}
#[test]
fn original_bezier_unclamped_and_c0_join_schemas_refuse_the_stronger_window_contract() {
    let original = input("0.0625", "1.");
    let header = "#900001=B_SPLINE_CURVE_WITH_KNOTS('',2,(#900010,#900011,#900012),.UNSPECIFIED.,.F.,.F.,(3,3),(0.,1.),.UNSPECIFIED.);";
    assert!(original.contains(header));
    let bezier = original.replace(
        header,
        "#900001=BEZIER_CURVE('',2,(#900010,#900011,#900012),.UNSPECIFIED.,.F.,.F.);",
    );
    let unclamped = original.replace("(3,3),(0.,1.)", "(1,1,1,1,1,1),(-2.,-1.,0.,1.,2.,3.)");
    let corner=append(&original.replace(header,"#900001=B_SPLINE_CURVE_WITH_KNOTS('',2,(#900010,#900011,#900012,#900013,#900014),.UNSPECIFIED.,.F.,.F.,(3,2,3),(0.,0.5,1.),.UNSPECIFIED.);"),"#900013=CARTESIAN_POINT('',(1.5,0.,0.));\n#900014=CARTESIAN_POINT('',(2.,1.,0.));");
    for text in [bezier, unclamped, corner] {
        assert_ne!(text, original);
        let source = Source::parse(&text, 1e-7).unwrap();
        let curve = actual_target(&source);
        assert!(matches!(
            knot_domain(
                &source,
                900001,
                &curve,
                true,
                &mut Budget::new(512_000_000, None)
            ),
            Err(ReadError::Unsupported)
        ));
    }
    let source = Source::parse(&original, 1e-7).unwrap();
    let curve = actual_target(&source);
    assert!(
        knot_domain(
            &source,
            900001,
            &curve,
            true,
            &mut Budget::new(512_000_000, None)
        )
        .is_ok()
    );
}
#[test]
fn genuine_nonmonotone_and_conditioned_original_curves_refuse_global_derivative_bounds() {
    let original = input("0.0625", "1.");
    // Both endpoints differ, but no Cartesian derivative has a strict global
    // lower bound: X changes sign and Y has a zero endpoint derivative.
    let nonmonotone = original
        .replace(
            "#900011=CARTESIAN_POINT('',(0.5,0.,0.))",
            "#900011=CARTESIAN_POINT('',(1.,0.,0.))",
        )
        .replace(
            "#900012=CARTESIAN_POINT('',(1.,1.,0.))",
            "#900012=CARTESIAN_POINT('',(0.,1.,0.))",
        );
    // Bounded real XYZ with an admissible extremely small parameter span.
    // Its derivative norm overflows; this is a real conditioning refusal.
    let conditioned = input("0.0625", "1E-300");
    for text in [nonmonotone, conditioned] {
        assert_ne!(text, original);
        let source = Source::parse(&text, 1e-7).unwrap();
        let curve = actual_target(&source);
        let mut b = Budget::new(512_000_000, None);
        let packet = super::super::super::bezier::original_curve(
            &source,
            900001,
            &curve,
            Coordinates::Raw(3),
            &mut b,
        )
        .unwrap();
        assert!(matches!(
            bounds(&packet, &mut b),
            Err(ReadError::Unsupported)
        ));
    }
}
#[test]
fn actual_original_evaluator_rejects_outside_domain_and_malformed_vertex_metadata() {
    let text = input("0.0625", "1.");
    let source = Source::parse(&text, 1e-7).unwrap();
    let (_, _, packet) = setup(&source);
    let mut b = Budget::new(512_000_000, None);
    assert!(matches!(
        evaluate(&packet, I::point(2.), &mut b),
        Err(ReadError::Unsupported)
    ));
    let p = evaluate(&packet, I::point(0.5), &mut b).unwrap();
    assert!(p[0].lo <= 0.5 && p[0].hi >= 0.5 && p[1].lo <= 0.25 && p[1].hi >= 0.25);
    let malformed = append(&text, "#990001=VERTEX_POINT('',#900010,$);");
    let source = Source::parse(&malformed, 1e-7).unwrap();
    assert!(matches!(
        raw_vertex(&source, 990001, &mut b),
        Err(ReadError::Invalid)
    ));
}
#[test]
fn every_observed_window_cancellation_checkpoint_is_atomic_and_retries() {
    let text = input("0.0625", "1.");
    let source = Source::parse(&text, 1e-7).unwrap();
    let (d, e, packet) = setup(&source);
    let calls = std::cell::Cell::new(0);
    let observe = || {
        calls.set(calls.get() + 1);
        false
    };
    prepare(
        &source,
        900001,
        &d,
        &e,
        &packet,
        &mut Budget::new(512_000_000, Some(&observe)),
    )
    .unwrap();
    let total = calls.get();
    assert!(total > 10);
    for checkpoint in 1..=total {
        calls.set(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() == checkpoint
        };
        let mut b = Budget::new(512_000_000, Some(&cancel));
        assert!(
            matches!(
                prepare(&source, 900001, &d, &e, &packet, &mut b),
                Err(ReadError::Limit)
            ),
            "checkpoint{checkpoint}"
        );
        assert!(b.cancelled.get());
    }
    assert!(
        prepare(
            &source,
            900001,
            &d,
            &e,
            &packet,
            &mut Budget::new(512_000_000, None)
        )
        .is_ok()
    );
}

#[test]
fn real_positive_rational_packet_is_not_a_polynomial_derivative_certificate() {
    let original = input("0.0625", "1.");
    let header = "#900001=B_SPLINE_CURVE_WITH_KNOTS('',2,(#900010,#900011,#900012),.UNSPECIFIED.,.F.,.F.,(3,3),(0.,1.),.UNSPECIFIED.);";
    let text=original.replace(header,"#900001=(BOUNDED_CURVE() B_SPLINE_CURVE(2,(#900010,#900011,#900012),.UNSPECIFIED.,.F.,.F.) B_SPLINE_CURVE_WITH_KNOTS((3,3),(0.,1.),.UNSPECIFIED.) CURVE() GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_CURVE((1.,2.,1.)) REPRESENTATION_ITEM(''));");
    assert_ne!(text, original);
    let source = Source::parse(&text, 1e-7).unwrap();
    let curve = actual_target(&source);
    assert!(curve.weights().is_some());
    let mut b = Budget::new(512_000_000, None);
    let packet = super::super::super::bezier::original_curve(
        &source,
        900001,
        &curve,
        Coordinates::Raw(3),
        &mut b,
    )
    .unwrap();
    assert!(matches!(
        bounds(&packet, &mut b),
        Err(ReadError::Unsupported)
    ));
}
#[test]
fn actual_foreign_parser_owner_packet_refuses_and_original_owner_retries() {
    let text = input("0.0625", "1.");
    let a = Source::parse(&text, 1e-7).unwrap();
    let b = Source::parse(&text, 1e-7).unwrap();
    let (da, ea, pa) = setup(&a);
    let (db, eb, _) = setup(&b);
    assert!(matches!(
        prepare(
            &b,
            900001,
            &db,
            &eb,
            &pa,
            &mut Budget::new(512_000_000, None)
        ),
        Err(ReadError::Unsupported)
    ));
    assert!(
        prepare(
            &a,
            900001,
            &da,
            &ea,
            &pa,
            &mut Budget::new(512_000_000, None)
        )
        .is_ok()
    );
}
#[test]
fn actual_original_vertex_cannot_earn_a_collapsed_or_reversed_closest_window() {
    let original = input("0.0625", "1.");
    let text = input("100.", "1.");
    assert_ne!(text, original);
    let source = Source::parse(&text, 1e-7).unwrap();
    let (d, e, packet) = setup(&source);
    assert!(matches!(
        prepare(
            &source,
            900001,
            &d,
            &e,
            &packet,
            &mut Budget::new(512_000_000, None)
        ),
        Err(ReadError::Unsupported)
    ));
}
