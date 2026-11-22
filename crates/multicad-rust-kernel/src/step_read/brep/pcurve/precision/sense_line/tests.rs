use super::*;
const BASE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../multicad-app/tests/fixtures/exchange/corpus/occt-cone-frustum.step"
));
const TOL: f64 = 1e-7;
/// An independently authored literal bilinear boundary. These raw operators
/// are parsed normally; no Source, Definition, spline or certificate is forged.
fn input(finite: bool) -> String {
    let uv = if finite {
        "#900030=B_SPLINE_CURVE_WITH_KNOTS('',1,(#900031,#900032),.UNSPECIFIED.,.F.,.F.,(2,2),(0.,4.),.UNSPECIFIED.);\n#900032=CARTESIAN_POINT('',(4.,0.));"
    } else {
        "#900030=LINE('',#900031,#900033);\n#900033=VECTOR('',#900034,1.);\n#900034=DIRECTION('',(1.,0.));"
    };
    let added = format!(
        r#"
#900001=B_SPLINE_SURFACE_WITH_KNOTS('',1,1,((#900010,#900011),(#900012,#900013)),.UNSPECIFIED.,.F.,.F.,.F.,(2,2),(2,2),(0.,4.),(0.,1.),.UNSPECIFIED.);
#900010=CARTESIAN_POINT('',(0.,0.,0.));
#900011=CARTESIAN_POINT('',(0.,0.,2.));
#900012=CARTESIAN_POINT('',(4.,0.,0.));
#900013=CARTESIAN_POINT('',(4.,0.,2.));
#900014=CARTESIAN_POINT('',(4.,0.,0.));
#900020=LINE('',#900010,#900021);
#900021=VECTOR('',#900022,1.);
#900022=DIRECTION('',(1.,0.,0.));
#900031=CARTESIAN_POINT('',(0.,0.));
{uv}
#900040=PCURVE('',#900001,#900041);
#900041=DEFINITIONAL_REPRESENTATION('',(#900030),#900042);
#900042=(GEOMETRIC_REPRESENTATION_CONTEXT(2) PARAMETRIC_REPRESENTATION_CONTEXT() REPRESENTATION_CONTEXT('',''));
"#
    );
    let at = BASE.rfind("ENDSEC;").unwrap();
    format!("{}{}{}", &BASE[..at], added, &BASE[at..])
}
fn inputs(source: &Source<'_>) -> (super::super::super::uv::Definition, Surface, Edge) {
    let mut budget = Budget::new(512_000_000, None);
    let definition = super::super::super::uv::definition(source, 900040, &mut budget).unwrap();
    let surface = source.surface(900001).unwrap();
    let curve = source.curve(900020, 0).unwrap();
    let mut vertices = [900010, 900014];
    let mut points = [source.point(900010).unwrap(), source.point(900014).unwrap()];
    let Curve::Line { direction, .. } = &curve else {
        panic!()
    };
    // Build the real finite LINE edge in its own positive parameter direction,
    // including the negative-target-sense fixture. This is input preparation;
    // no sampled/cached ordering grants the separate OriginalREAL IP2 proof.
    let forward: f64 = (0..3)
        .map(|i| (points[1][i] - points[0][i]) * direction[i])
        .sum();
    if forward < 0. {
        vertices.swap(0, 1);
        points.swap(0, 1);
    }
    let edge = Edge::new(curve, vertices, points, TOL).unwrap();
    (definition, surface, edge)
}
#[test]
fn actual_original_bilinear_boundaries_keep_positive_sense_for_line_and_finite_uv() {
    for finite in [false, true] {
        let input = input(finite);
        let source = Source::parse(&input, TOL).unwrap();
        let (definition, surface, edge) = inputs(&source);
        certify(
            &source,
            900020,
            &definition,
            &surface,
            &edge,
            &mut Budget::new(512_000_000, None),
        )
        .unwrap();
    }
}
#[test]
fn original_reverse_zero_and_hidden_underflow_coordinates_never_infer_forward_sense() {
    let original = input(false);
    for (from, to) in [
        (
            "#900034=DIRECTION('',(1.,0.));",
            "#900034=DIRECTION('',(-1.,0.));",
        ),
        (
            "#900034=DIRECTION('',(1.,0.));",
            "#900034=DIRECTION('',(1.,1E-1000));",
        ),
        (
            "#900012=CARTESIAN_POINT('',(4.,0.,0.));",
            "#900012=CARTESIAN_POINT('',(0.,0.,0.));",
        ),
        (
            "#900022=DIRECTION('',(1.,0.,0.));",
            "#900022=DIRECTION('',(-1.,0.,0.));",
        ),
    ] {
        let input = original.replace(from, to);
        assert_ne!(input, original);
        let source = Source::parse(&input, TOL).unwrap();
        let (definition, surface, edge) = inputs(&source);
        // The degenerate-row test changes the boundary operator, not the
        // unchanged actual LINE's used endpoint data.
        assert!(
            certify(
                &source,
                900020,
                &definition,
                &surface,
                &edge,
                &mut Budget::new(512_000_000, None)
            )
            .is_err()
        );
    }
}
#[test]
fn literal_positive_rational_boundary_weights_preserve_the_whole_derivative_sign() {
    let base = input(true);
    let old = "#900001=B_SPLINE_SURFACE_WITH_KNOTS('',1,1,((#900010,#900011),(#900012,#900013)),.UNSPECIFIED.,.F.,.F.,.F.,(2,2),(2,2),(0.,4.),(0.,1.),.UNSPECIFIED.);";
    let new = "#900001=(BOUNDED_SURFACE() B_SPLINE_SURFACE(1,1,((#900010,#900011),(#900012,#900013)),.UNSPECIFIED.,.F.,.F.,.F.) B_SPLINE_SURFACE_WITH_KNOTS((2,2),(2,2),(0.,4.),(0.,1.),.UNSPECIFIED.) GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_SURFACE(((2.,3.),(5.,7.))) REPRESENTATION_ITEM('') SURFACE());";
    let input = base.replace(old, new);
    assert_ne!(base, input);
    let source = Source::parse(&input, TOL).unwrap();
    let (definition, surface, edge) = inputs(&source);
    certify(
        &source,
        900020,
        &definition,
        &surface,
        &edge,
        &mut Budget::new(512_000_000, None),
    )
    .unwrap();
}
#[test]
fn actual_sense_work_frontier_late_cancel_and_retry_are_source_stable() {
    let input = input(true);
    let source = Source::parse(&input, TOL).unwrap();
    let (definition, surface, edge) = inputs(&source);
    let run = |b: &mut Budget<'_>| certify(&source, 900020, &definition, &surface, &edge, b);
    let mut full = Budget::new(512_000_000, None);
    run(&mut full).unwrap();
    let used = full.used();
    let mut exact = Budget::new(used, None);
    run(&mut exact).unwrap();
    assert_eq!(exact.remaining, 0);
    assert!(matches!(
        run(&mut Budget::new(used - 1, None)),
        Err(ReadError::Limit)
    ));
    let calls = std::cell::Cell::new(0);
    let stop = || {
        calls.set(calls.get() + 1);
        calls.get() == 20
    };
    let mut cancelled = Budget::new(512_000_000, Some(&stop));
    assert!(run(&mut cancelled).is_err());
    assert!(cancelled.cancelled.get());
    assert_eq!(calls.get(), 20);
    let mut retry = Budget::new(512_000_000, None);
    run(&mut retry).unwrap();
    assert_eq!(retry.used(), used);
}

#[test]
fn actual_noncardinal_and_other_axis_boundaries_use_original_vector_flow() {
    let original = input(false);
    let posed = original
        .replace(
            "#900012=CARTESIAN_POINT('',(4.,0.,0.));",
            "#900012=CARTESIAN_POINT('',(2.4,3.2,0.));",
        )
        .replace(
            "#900013=CARTESIAN_POINT('',(4.,0.,2.));",
            "#900013=CARTESIAN_POINT('',(2.4,3.2,2.));",
        )
        .replace(
            "#900014=CARTESIAN_POINT('',(4.,0.,0.));",
            "#900014=CARTESIAN_POINT('',(2.4,3.2,0.));",
        )
        .replace(
            "#900022=DIRECTION('',(1.,0.,0.));",
            "#900022=DIRECTION('',(0.6,0.8,0.));",
        );
    let other_axis = original
        .replace(
            "#900034=DIRECTION('',(1.,0.));",
            "#900034=DIRECTION('',(0.,1.));",
        )
        .replace(
            "#900022=DIRECTION('',(1.,0.,0.));",
            "#900022=DIRECTION('',(0.,0.,1.));",
        )
        .replace(
            "#900014=CARTESIAN_POINT('',(4.,0.,0.));",
            "#900014=CARTESIAN_POINT('',(0.,0.,2.));",
        );
    for input in [posed, other_axis] {
        let source = Source::parse(&input, TOL).unwrap();
        let (definition, surface, edge) = inputs(&source);
        certify(
            &source,
            900020,
            &definition,
            &surface,
            &edge,
            &mut Budget::new(512_000_000, None),
        )
        .unwrap();
    }
}
