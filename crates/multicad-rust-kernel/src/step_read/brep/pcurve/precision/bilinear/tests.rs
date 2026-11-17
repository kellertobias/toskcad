use super::*;
const TOL: f64 = 1e-7;
const BASE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../multicad-app/tests/fixtures/exchange/corpus/occt-cone-frustum.step"
));
#[test]
fn bernstein_products_and_homogeneous_derivatives_enclose_independent_polynomials() {
    let mut budget = Budget::new(512_000_000, None);
    let a = Poly::from_controls(&[[I::point(1.); 4], [I::point(2.); 4]], 0).unwrap();
    let b = Poly::from_controls(&[[I::point(3.); 4], [I::point(5.); 4]], 0).unwrap();
    let p = product(a, b, &mut budget).unwrap();
    // (1+t)(3+2t)=3+5t+2t² has quadratic Bernstein coefficients
    // [3, 11/2, 10], independently of the product implementation.
    for (coefficient, exact) in p.c[..=2].iter().zip([3., 5.5, 10.]) {
        assert!(coefficient.lo <= exact && exact <= coefficient.hi);
    }
    // X/W=(1+t)/(3+2t): X'W-XW'=1 identically.
    let z = Poly::zero(1).unwrap();
    let derivative = derivative_numerator([a, z, z, b], &mut budget).unwrap();
    for c in &derivative[0].c[..=derivative[0].degree] {
        assert!(c.lo <= 1. && 1. <= c.hi);
    }
    {
        let [left, right] = p.split(&mut budget).unwrap();
        assert!(left.c[0].lo <= 3. && left.c[0].hi >= 3.);
        assert!(left.c[2].lo <= 6. && left.c[2].hi >= 6.);
        assert!(right.c[0].lo <= 6. && right.c[0].hi >= 6.);
        assert!(right.c[2].lo <= 10. && right.c[2].hi >= 10.);
    }
}
fn input(spline: bool, rational: bool) -> String {
    let uv = if rational {
        "#900030=(B_SPLINE_CURVE(2,(#900031,#900032,#900033),.UNSPECIFIED.,.F.,.F.) B_SPLINE_CURVE_WITH_KNOTS((3,3),(2.,7.),.UNSPECIFIED.) RATIONAL_B_SPLINE_CURVE((1.,2.,1.)) REPRESENTATION_ITEM('')) ;"
    } else {
        "#900030=B_SPLINE_CURVE_WITH_KNOTS('',2,(#900031,#900032,#900033),.UNSPECIFIED.,.F.,.F.,(3,3),(2.,7.),.UNSPECIFIED.);"
    };
    let target = if spline && rational {
        "#900020=(B_SPLINE_CURVE(2,(#900010,#900015,#900013),.UNSPECIFIED.,.F.,.F.) B_SPLINE_CURVE_WITH_KNOTS((3,3),(-1.,4.),.UNSPECIFIED.) RATIONAL_B_SPLINE_CURVE((1.,2.,1.)) REPRESENTATION_ITEM(''));"
    } else if spline {
        "#900020=B_SPLINE_CURVE_WITH_KNOTS('',2,(#900010,#900015,#900013),.UNSPECIFIED.,.F.,.F.,(3,3),(-1.,4.),.UNSPECIFIED.);"
    } else {
        "#900020=LINE('',#900010,#900021);\n#900021=VECTOR('',#900022,1.);\n#900022=DIRECTION('',(1.,1.,0.));"
    };
    let extra = format!(
        r#"
#900001=B_SPLINE_SURFACE_WITH_KNOTS('',1,1,((#900010,#900011),(#900012,#900013)),.UNSPECIFIED.,.F.,.F.,.F.,(2,2),(2,2),(0.,1.),(0.,1.),.UNSPECIFIED.);
#900010=CARTESIAN_POINT('',(0.,0.,0.));
#900011=CARTESIAN_POINT('',(0.,1.,0.));
#900012=CARTESIAN_POINT('',(1.,0.,0.));
#900013=CARTESIAN_POINT('',(1.,1.,0.));
#900015=CARTESIAN_POINT('',(0.25,0.25,0.));
{target}
#900031=CARTESIAN_POINT('',(0.,0.));
#900032=CARTESIAN_POINT('',(0.25,0.25));
#900033=CARTESIAN_POINT('',(1.,1.));
{uv}
#900040=PCURVE('',#900001,#900041);
#900041=DEFINITIONAL_REPRESENTATION('',(#900030),#900042);
#900042=(GEOMETRIC_REPRESENTATION_CONTEXT(2) PARAMETRIC_REPRESENTATION_CONTEXT() REPRESENTATION_CONTEXT('',''));
#900050=VERTEX_POINT('',#900010);
#900051=VERTEX_POINT('',#900013);
"#
    );
    let at = BASE.rfind("ENDSEC;").unwrap();
    format!("{}{}{}", &BASE[..at], extra, &BASE[at..])
}
fn run(source: &Source<'_>, budget: &mut Budget<'_>) -> Result<Option<Proof>> {
    let d = super::super::super::uv::definition(source, 900040, budget)?;
    let surface = source.surface(900001)?;
    let target = source.curve(900020, 0)?;
    let edge = Edge::new(
        target,
        [900050, 900051],
        [source.point(900010)?, source.point(900013)?],
        TOL,
    )?;
    certify(source, 900020, &d, &surface, &edge, TOL, budget)
}
#[test]
fn actual_original_quadratic_traces_compose_bilinear_surfaces_with_positive_maps() {
    for spline in [false, true] {
        for rational in [false, true] {
            let text = input(spline, rational);
            let source = Source::parse(&text, TOL).unwrap();
            assert!(
                matches!(run(&source, &mut Budget::new(512_000_000, None)).unwrap(), Some(Proof::Certified(error)) if error <= TOL)
            );
        }
    }
}
#[test]
fn actual_nonplanar_bilinear_image_uses_full_composed_quadratic_coefficients() {
    let text = input(true, false)
        .replace("#900013=CARTESIAN_POINT('',(1.,1.,0.));", "#900013=CARTESIAN_POINT('',(1.,1.,1.));")
        .replace("#900015=CARTESIAN_POINT('',(0.25,0.25,0.));", "#900015=CARTESIAN_POINT('',(0.5,0.5,0.));")
        .replace("#900030=B_SPLINE_CURVE_WITH_KNOTS('',2,(#900031,#900032,#900033),.UNSPECIFIED.,.F.,.F.,(3,3),(2.,7.),.UNSPECIFIED.);", "#900030=B_SPLINE_CURVE_WITH_KNOTS('',1,(#900031,#900033),.UNSPECIFIED.,.F.,.F.,(2,2),(2.,7.),.UNSPECIFIED.);");
    let source = Source::parse(&text, TOL).unwrap();
    // Actual surface P(u,v)=(u,v,uv), trace (t,t), and target
    // (t,t,t²) have different degrees but the same complete image.
    assert!(
        matches!(run(&source, &mut Budget::new(512_000_000, None)).unwrap(), Some(Proof::Certified(error)) if error <= TOL)
    );
}
#[test]
fn original_bent_curve_and_reversed_target_never_gain_a_positive_certificate() {
    let original = input(false, false);
    for text in [
        original.replace("(0.25,0.25));", "(0.25,0.4));"),
        original.replace("DIRECTION('',(1.,1.,0.))", "DIRECTION('',(-1.,-1.,0.))"),
    ] {
        let source = Source::parse(&text, TOL).unwrap();
        // A descending LINE is an actual source refusal before composition;
        // a bent source remains unqualified by the whole coefficient bound.
        assert!(!matches!(
            run(&source, &mut Budget::new(512_000_000, None)),
            Ok(Some(Proof::Certified(_)))
        ));
    }
}
#[test]
fn actual_resource_frontier_and_cancel_leave_original_owner_retryable() {
    let text = input(true, true);
    let source = Source::parse(&text, TOL).unwrap();
    let mut full = Budget::new(512_000_000, None);
    assert!(matches!(
        run(&source, &mut full).unwrap(),
        Some(Proof::Certified(_))
    ));
    let used = full.used();
    assert!(matches!(
        run(&source, &mut Budget::new(used - 1, None)),
        Err(ReadError::Limit)
    ));
    assert!(matches!(
        run(&source, &mut Budget::new(used, None)).unwrap(),
        Some(Proof::Certified(_))
    ));
    let cancel = || true;
    let mut interrupted = Budget::new(512_000_000, Some(&cancel));
    assert!(matches!(
        run(&source, &mut interrupted),
        Err(ReadError::Limit)
    ));
    assert!(interrupted.cancelled.get());
    assert!(matches!(
        run(&source, &mut Budget::new(512_000_000, None)).unwrap(),
        Some(Proof::Certified(_))
    ));
}
