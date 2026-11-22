use super::*;
const TOL: f64 = 1e-7;
#[test]
fn actual_occt_supplementary_trim_matrix() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let mut failures = Vec::new();
    let mut pairs = 0;
    for name in [
        "box-void",
        "cone-frustum",
        "cylinder",
        "placed-box",
        "pocketed-plate",
        "revolved-vase",
        "sphere",
        "spline-prism",
        "torus",
        "two-bodies",
    ] {
        let input = std::fs::read_to_string(root.join(format!(
            "../multicad-app/tests/fixtures/exchange/corpus/occt-{name}.step"
        )))
        .unwrap();
        match validate_pcurves(&input, TOL) {
            Ok(report) => {
                eprintln!(
                    "PCURVE {name}: {} curves, {} work",
                    report.pcurves_checked, report.work_used
                );
                pairs += report.pcurves_checked;
            }
            Err(error) => failures.push((name, error)),
        }
    }
    for name in ["regular-torus", "spline-prism", "annulus"] {
        let input = std::fs::read_to_string(
            root.join(format!("src/step_read/tests/fixtures/occt-{name}.step")),
        )
        .unwrap();
        match validate_pcurves(&input, TOL) {
            Ok(report) => {
                eprintln!(
                    "PCURVE pure/{name}: {} curves, {} work",
                    report.pcurves_checked, report.work_used
                );
                pairs += report.pcurves_checked;
            }
            Err(error) => failures.push((name, error)),
        }
    }
    assert!(
        failures.is_empty(),
        "actual supplementary trim failures: {failures:?}"
    );
    assert!(pairs > 100);
}
#[test]
fn byte_scan_and_cancellation_receipt_are_bound_to_actual_text_and_tolerance() {
    let input = "ISO-10303-21; no supplementary curves";
    let receipt = validate_pcurve_input(input, TOL, None).unwrap();
    assert!(std::ptr::eq(receipt.input(), input));
    assert_eq!(receipt.tolerance(), TOL);
    assert_eq!(receipt.report().work_used, input.len() * 4);
    assert_eq!(receipt.report().pcurves_checked, 0);
    let cancelled = validate_pcurve_input(input, TOL, Some(&|| true)).unwrap_err();
    assert!(cancelled.cancelled);
    for tolerance in [0., -1., f64::NAN, f64::INFINITY] {
        assert_eq!(
            validate_pcurves(input, tolerance).unwrap_err().kind,
            ReadError::Invalid
        );
    }
    assert_eq!(validate_pcurves(input, TOL).unwrap(), *receipt.report());
}
const PLACED: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../multicad-app/tests/fixtures/exchange/corpus/occt-placed-box.step"
));
fn replace_once(input: &str, from: &str, to: &str) -> String {
    assert_eq!(input.matches(from).count(), 1);
    input.replacen(from, to, 1)
}
#[test]
fn selected_and_auxiliary_trims_have_same_geometric_and_sense_requirements() {
    for master in [".PCURVE_S1.", ".PCURVE_S2.", ".CURVE_3D."] {
        let input = replace_once(
            PLACED,
            "#26 = SURFACE_CURVE('',#27,(#31,#43),.PCURVE_S1.);",
            &format!("#26 = SURFACE_CURVE('',#27,(#31,#43),{master});"),
        );
        assert_eq!(validate_pcurves(&input, TOL).unwrap().pcurves_checked, 24);
        let shifted = replace_once(
            &input,
            "#51 = CARTESIAN_POINT('',(16.,0.));",
            "#51 = CARTESIAN_POINT('',(16.1,0.));",
        );
        let error = validate_pcurves(&shifted, TOL).unwrap_err();
        assert_eq!(error.kind, ReadError::Invalid);
        assert!(error.disagreements.iter().any(|d| d.edge_entity == 21
            && d.surface_entity == 44
            && d.pcurve_entity == 43
            && d.deviation_mm > 0.099999));
        let reversed = replace_once(
            &input,
            "#53 = DIRECTION('',(0.,-1.));",
            "#53 = DIRECTION('',(0.,1.));",
        );
        assert_eq!(
            validate_pcurves(&reversed, TOL).unwrap_err().kind,
            ReadError::Unsupported
        );
        // Positive affine reparameterization, including a different line
        // origin along the same support, retains the complete image and sense.
        let moved = replace_once(
            &input,
            "#51 = CARTESIAN_POINT('',(16.,0.));",
            "#51 = CARTESIAN_POINT('',(16.,5.));",
        );
        let scaled = replace_once(
            &moved,
            "#52 = VECTOR('',#53,1.);",
            "#52 = VECTOR('',#53,3.);",
        );
        assert_eq!(validate_pcurves(&scaled, TOL).unwrap().pcurves_checked, 24);
    }
}
#[test]
fn all_proved_auxiliary_disagreements_are_reported_before_refusal() {
    let input = replace_once(
        PLACED,
        "#51 = CARTESIAN_POINT('',(16.,0.));",
        "#51 = CARTESIAN_POINT('',(16.1,0.));",
    );
    let input = replace_once(
        &input,
        "#39 = CARTESIAN_POINT('',(0.,0.));",
        "#39 = CARTESIAN_POINT('',(0.2,0.));",
    );
    let error = validate_pcurves(&input, TOL).unwrap_err();
    assert_eq!(error.kind, ReadError::Invalid);
    assert!(
        error
            .disagreements
            .iter()
            .any(|d| d.edge_entity == 21 && d.pcurve_entity == 31 && d.deviation_mm > 0.199999)
    );
    assert!(
        error
            .disagreements
            .iter()
            .any(|d| d.edge_entity == 21 && d.pcurve_entity == 43 && d.deviation_mm > 0.099999)
    );
    // A new unowned construction PCURVE is not an owned finite edge trim.
    let unused = replace_once(
        PLACED,
        "ENDSEC;\nEND-ISO-10303-21;",
        "#900001 = PCURVE('',#44,#49);\nENDSEC;\nEND-ISO-10303-21;",
    );
    assert_eq!(validate_pcurves(&unused, TOL).unwrap().pcurves_checked, 24);
}
#[test]
fn weighted_finite_uv_isocurve_keeps_its_domain_and_positive_sense() {
    use crate::brep::{NurbsCurve, NurbsSurface};
    let source = Surface::Nurbs {
        surface: NurbsSurface::new(
            [2, 1],
            [3, 2],
            vec![
                [0., 0., 0.],
                [0., 1., 0.],
                [1., 0., 1.],
                [1., 1., 1.],
                [2., 0., 0.],
                [2., 1., 0.],
            ],
            None,
            [vec![0., 1.], vec![0., 1.]],
            [vec![3, 3], vec![2, 2]],
        )
        .unwrap(),
        closed: [false, false],
    };
    let target = Curve::Nurbs {
        curve: NurbsCurve::new(
            2,
            vec![[0., 0., 0.], [1., 0., 1.], [2., 0., 0.]],
            None,
            vec![4., 10.],
            vec![3, 3],
        )
        .unwrap(),
        closed: false,
    };
    let edge = Edge::new(target, [0, 1], [[0., 0., 0.], [2., 0., 0.]], TOL).unwrap();
    let uv = |a: f64, b: f64| Curve::Nurbs {
        curve: NurbsCurve::new(
            1,
            vec![[a, 0., 0.], [b, 0., 0.]],
            Some(vec![1., 3.]),
            vec![-9., 2.],
            vec![2, 2],
        )
        .unwrap(),
        closed: false,
    };
    let run = |curve: &Curve| {
        affine::certify(&source, curve, &edge, TOL, &mut Budget::new(MAX_WORK, None)).unwrap()
    };
    assert!(matches!(run(&uv(0., 1.)), Proof::Certified(_)));
    assert!(matches!(run(&uv(0., 0.5)), Proof::Unsupported));
    assert!(matches!(run(&uv(1., 0.)), Proof::Unsupported));
    assert!(matches!(
        affine::at(&uv(0., 1.), -10., &mut Budget::new(MAX_WORK, None)),
        Err(ReadError::Unsupported)
    ));
    let wrapped = match uv(0., 1.) {
        Curve::Nurbs { curve, .. } => Curve::Nurbs {
            curve,
            closed: true,
        },
        _ => unreachable!(),
    };
    assert!(matches!(run(&wrapped), Proof::Unsupported));
}
#[test]
fn rational_straight_uv_spline_can_cover_a_differently_scaled_line_but_not_reverse_it() {
    let line = Curve::Line {
        origin: [0., 0., 0.],
        direction: [2., 0., 0.],
    };
    let edge = Edge::new(line, [0, 1], [[0., 0., 0.], [4., 0., 0.]], TOL).unwrap();
    let controls = [[0., 0., 0.], [4., 0., 0.]].map(vector);
    let run = |c: &[[I; 3]], forward| {
        affine::nurbs(
            1,
            (&[-7., 3.], &[2, 2]),
            Some(&[1., 9.]),
            c,
            forward,
            &edge,
            TOL,
            &mut Budget::new(MAX_WORK, None),
        )
        .unwrap()
    };
    assert!(matches!(run(&controls, true), Proof::Certified(_)));
    assert!(matches!(run(&controls, false), Proof::Unsupported));
    let shifted = [[0., 0.1, 0.], [4., 0.1, 0.]].map(vector);
    assert!(matches!(run(&shifted,true),Proof::Disagreement(lower) if lower>0.09999));
}
#[test]
fn complete_curve_proof_budget_and_observed_cancel_retry_preserve_atomic_receipt() {
    let source = Source::parse(PLACED, TOL).unwrap();
    let mut budget = Budget::new(MAX_WORK, None);
    let original = validate_source(&source, TOL, &mut budget).unwrap();
    let used = budget.used();
    assert!(used > 1000);
    assert_eq!(
        validate_source(&source, TOL, &mut Budget::new(used - 1, None))
            .unwrap_err()
            .kind,
        ReadError::Limit
    );
    assert_eq!(
        validate_source(&source, TOL, &mut Budget::new(used, None)).unwrap(),
        original
    );
    let calls = std::cell::Cell::new(0_usize);
    let stop = || {
        let next = calls.get() + 1;
        calls.set(next);
        next > 40
    };
    let error = validate_pcurve_input(PLACED, TOL, Some(&stop)).unwrap_err();
    assert!(error.cancelled);
    assert!(calls.get() > 40);
    let retry = validate_pcurve_input(PLACED, TOL, None).unwrap();
    assert_eq!(retry.report().pcurves_checked, 24);
    // Original REAL 2D context dimension is checked, not a rounded cache alias.
    let malformed = replace_once(
        PLACED,
        "#54 = ( GEOMETRIC_REPRESENTATION_CONTEXT(2)",
        "#54 = ( GEOMETRIC_REPRESENTATION_CONTEXT(2.0000000000000001)",
    );
    assert_eq!(
        validate_pcurves(&malformed, TOL).unwrap_err().kind,
        ReadError::Unsupported
    );
}
#[test]
fn complete_nonuniform_nurbs_basis_requires_exact_positive_knot_map() {
    use crate::brep::NurbsCurve;
    let points = vec![[0., 0., 0.], [1., 2., 0.], [2., -1., 0.], [3., 0., 0.]];
    let target = |middle| {
        let curve = Curve::Nurbs {
            curve: NurbsCurve::new(
                2,
                points.clone(),
                None,
                vec![2., middle, 10.],
                vec![3, 1, 3],
            )
            .unwrap(),
            closed: false,
        };
        Edge::new(curve, [0, 1], [points[0], points[3]], TOL).unwrap()
    };
    let controls: Vec<_> = points.iter().copied().map(vector).collect();
    let run = |edge: &Edge| {
        affine::nurbs(
            2,
            (&[0., 0.25, 1.], &[3, 1, 3]),
            None,
            &controls,
            true,
            edge,
            TOL,
            &mut Budget::new(MAX_WORK, None),
        )
        .unwrap()
    };
    assert!(matches!(run(&target(4.)), Proof::Certified(_)));
    // A rounded equality of derived normalized knots is not original
    // coefficient correspondence. Exact dyadic cross-products reject it.
    assert!(matches!(run(&target(4_f64.next_up())), Proof::Unsupported));
    let mut bumped = controls.clone();
    bumped[1][2] = I::point(0.3);
    assert!(!matches!(
        affine::nurbs(
            2,
            (&[0., 0.25, 1.], &[3, 1, 3]),
            None,
            &bumped,
            true,
            &target(4.),
            TOL,
            &mut Budget::new(MAX_WORK, None)
        )
        .unwrap(),
        Proof::Certified(_)
    ));
}
#[test]
fn nominated_master_must_exist_and_be_an_actual_pcurve() {
    let original = "#26 = SURFACE_CURVE('',#27,(#31,#43),.PCURVE_S1.);";
    for replacement in [
        "#26 = SURFACE_CURVE('',#27,(#31),.PCURVE_S2.);",
        "#26 = SURFACE_CURVE('',#27,(#32,#43),.PCURVE_S1.);",
    ] {
        let input = replace_once(PLACED, original, replacement);
        assert_eq!(
            validate_pcurves(&input, TOL).unwrap_err().kind,
            ReadError::Invalid
        );
    }
    // This valid direct-support association contains no supplementary curve
    // at all. The 3D geometry still belongs to the ordinary retained reader.
    let input = replace_once(
        PLACED,
        original,
        "#26 = SURFACE_CURVE('',#27,(#32,#44),.CURVE_3D.);",
    );
    assert_eq!(validate_pcurves(&input, TOL).unwrap().pcurves_checked, 22);
}
fn metadata_file(records: &str) -> String {
    format!(
        "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('test'),'2;1');\nFILE_NAME('','',(''),(''),'','','');\nFILE_SCHEMA(('AP242'));\nENDSEC;\nDATA;\n{records}\nENDSEC;\nEND-ISO-10303-21;"
    )
}
#[test]
fn quoted_or_unowned_curve_names_do_not_impose_a_retained_solid_provider() {
    for record in [
        "#1=GEOMETRIC_SET('PCURVE SURFACE_CURVE',());",
        "/* PCURVE SURFACE_CURVE */ #1=GEOMETRIC_SET('',());",
        // An unused construction definition is not an edge-owned trim.
        "#1=PCURVE('',#2,#3);#2=PLANE('',#4);",
    ] {
        let input = metadata_file(record);
        assert!(Source::parse(&input, TOL).is_err());
        assert_eq!(validate_pcurves(&input, TOL).unwrap().pcurves_checked, 0);
    }
    // Even without a single PCURVE entity or keyword in a label, a selected
    // missing supplementary representation is refused by the exact AST pass.
    let input = metadata_file(
        "#1=EDGE_CURVE('',#8,#9,#2,.T.);#2=SURFACE_CURVE('',#3,(#4),.PCURVE_S1.);#3=LINE('',#5,#6);#4=PLANE('',#7);",
    );
    assert_eq!(
        validate_pcurves(&input, TOL).unwrap_err().kind,
        ReadError::Invalid
    );
}
#[test]
fn recursive_geometry_value_visits_are_billed_even_without_any_references() {
    let numbers = std::iter::repeat_n("1.", 2048)
        .collect::<Vec<_>>()
        .join(",");
    let input = replace_once(
        PLACED,
        "ENDSEC;\nEND-ISO-10303-21;",
        &format!("#900001=UNSUPPORTED_SURFACE('',({numbers}));\nENDSEC;\nEND-ISO-10303-21;"),
    );
    let source = Source::parse(&input, TOL).unwrap();
    let mut budget = Budget::new(MAX_WORK, None);
    geometry_preflight(&source, 900001, &mut budget).unwrap();
    let used = budget.used();
    assert!(used >= 2048 * 16);
    assert_eq!(
        geometry_preflight(&source, 900001, &mut Budget::new(used - 1, None)),
        Err(ReadError::Limit)
    );
    assert_eq!(
        geometry_preflight(&source, 900001, &mut Budget::new(used, None)),
        Ok(())
    );
    let calls = std::cell::Cell::new(0);
    let stop = || {
        calls.set(calls.get() + 1);
        calls.get() > 100
    };
    let mut budget = Budget::new(MAX_WORK, Some(&stop));
    assert_eq!(
        geometry_preflight(&source, 900001, &mut budget),
        Err(ReadError::Limit)
    );
    assert!(budget.cancelled.get());
}
