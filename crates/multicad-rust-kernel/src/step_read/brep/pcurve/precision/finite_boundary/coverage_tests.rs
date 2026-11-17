//! Real parsed schemas and observed callback frontiers; no certificate edits.
use super::*;
fn record(text: &str, id: u32, new: &str) -> String {
    let begin = text.find(&format!("#{id}=")).unwrap();
    let end = begin + text[begin..].find(';').unwrap() + 1;
    format!("{}{}{}", &text[..begin], new, &text[end..])
}
#[test]
fn actual_plane_and_higher_degree_surface_decline_the_narrow_provider() {
    let original = input(false);
    let plane = original.replace("#900100=PCURVE('',#900040,", "#900100=PCURVE('',#31,");
    let thirdrow = "#900040=B_SPLINE_SURFACE_WITH_KNOTS('',2,2,((#900050,#900051,#900052),(#900050,#900051,#900052),(#900053,#900054,#900055)),.UNSPECIFIED.,.F.,.F.,.F.,(3,3),(3,3),(0.,1.),(0.,6.323173064558),.UNSPECIFIED.);";
    let three = record(&original, 900040, thirdrow);
    for text in [plane, three] {
        assert_ne!(text, original);
        let source = Source::parse(&text, 1e-7).unwrap();
        let mut b = Budget::new(512_000_000, None);
        let (d, s, e) = pair(&source, &mut b);
        assert!(
            certify(&source, 900001, &d, &s, &e, 1e-7, &mut b)
                .unwrap()
                .is_none()
        );
    }
}
#[test]
fn real_bezier_metadata_and_interior_isoline_never_bypass_boundary_conditions() {
    let original = input(false);
    let uv = record(
        &original,
        900102,
        "#900102=BEZIER_CURVE('',1,(#900110,#900111),.UNSPECIFIED.,.F.,.F.);",
    );
    let surface = record(
        &original,
        900040,
        "#900040=BEZIER_SURFACE('',1,2,((#900050,#900051,#900052),(#900053,#900054,#900055)),.UNSPECIFIED.,.F.,.F.,.F.);",
    );
    let interior = original
        .replace(
            "#900110=CARTESIAN_POINT('',(0.,0.))",
            "#900110=CARTESIAN_POINT('',(0.5,0.))",
        )
        .replace(
            "#900111=CARTESIAN_POINT('',(0.,3.161586532279))",
            "#900111=CARTESIAN_POINT('',(0.5,3.161586532279))",
        );
    for text in [uv, surface, interior] {
        assert_ne!(text, original);
        let source = Source::parse(&text, 1e-7).unwrap();
        let mut b = Budget::new(512_000_000, None);
        let (d, s, e) = pair(&source, &mut b);
        assert!(
            certify(&source, 900001, &d, &s, &e, 1e-7, &mut b)
                .unwrap()
                .is_none()
        );
    }
    let source = Source::parse(&original, 1e-7).unwrap();
    let mut b = Budget::new(512_000_000, None);
    let (d, s, e) = pair(&source, &mut b);
    assert!(
        certify(&source, 900001, &d, &s, &e, 1e-7, &mut b)
            .unwrap()
            .is_some()
    );
}
#[test]
fn actual_horizontal_uv_declines_the_private_v_parameter_operator() {
    let original = input(false);
    let text = original.replace(
        "#900111=CARTESIAN_POINT('',(0.,3.161586532279))",
        "#900111=CARTESIAN_POINT('',(1.,0.))",
    );
    assert_ne!(text, original);
    let source = Source::parse(&text, 1e-7).unwrap();
    let mut b = Budget::new(512_000_000, None);
    let (d, _, _) = pair(&source, &mut b);
    // Helper-only restriction: the general public route declines this valid
    // U-moving curve earlier. The private V(t)=t operator must also refuse.
    assert!(matches!(
        original_uv(&source, &d, &mut b),
        Err(ReadError::Unsupported)
    ));
}
#[test]
fn every_observed_finite_boundary_checkpoint_cancels_without_publication_and_retries() {
    let text = input(false);
    let source = Source::parse(&text, 1e-7).unwrap();
    let (d, s, e) = pair(&source, &mut Budget::new(512_000_000, None));
    let calls = std::cell::Cell::new(0);
    let observe = || {
        calls.set(calls.get() + 1);
        false
    };
    assert!(
        certify(
            &source,
            900001,
            &d,
            &s,
            &e,
            1e-7,
            &mut Budget::new(512_000_000, Some(&observe))
        )
        .unwrap()
        .is_some()
    );
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
                certify(&source, 900001, &d, &s, &e, 1e-7, &mut b),
                Err(ReadError::Limit)
            ),
            "checkpoint{checkpoint}"
        );
        assert!(b.cancelled.get());
    }
    assert!(
        certify(
            &source,
            900001,
            &d,
            &s,
            &e,
            1e-7,
            &mut Budget::new(512_000_000, None)
        )
        .unwrap()
        .is_some()
    );
}

#[test]
fn independently_parsed_unsupported_coefficient_schemas_refuse_before_any_bound() {
    let original = input(false);
    let three = record(
        &original,
        900040,
        "#900040=B_SPLINE_SURFACE_WITH_KNOTS('',2,2,((#900050,#900051,#900052),(#900050,#900051,#900052),(#900053,#900054,#900055)),.UNSPECIFIED.,.F.,.F.,.F.,(3,3),(3,3),(0.,1.),(0.,6.323173064558),.UNSPECIFIED.);",
    );
    let bezier = record(
        &original,
        900001,
        "#900001=BEZIER_CURVE('',2,(#900010,#900011,#900012),.UNSPECIFIED.,.F.,.F.);",
    );
    let unclamped = record(
        &original,
        900001,
        "#900001=B_SPLINE_CURVE_WITH_KNOTS('',2,(#900010,#900011,#900012),.UNSPECIFIED.,.F.,.F.,(2,1,3),(0.,0.5,6.323173064558),.UNSPECIFIED.);",
    );
    for text in [three, bezier, unclamped] {
        assert_ne!(text, original);
        let source = Source::parse(&text, 1e-7).unwrap();
        let mut b = Budget::new(512_000_000, None);
        let d = super::super::super::super::uv::definition(&source, 900100, &mut b).unwrap();
        let support = source.surface(d.surface).unwrap();
        let Surface::Nurbs { surface, .. } = support else {
            panic!("actual parsed tensor")
        };
        let Curve::Nurbs { curve: target, .. } = source.curve(900001, 0).unwrap() else {
            panic!("actual parsed spline")
        };
        let iso = original_uv(&source, &d, &mut b).unwrap();
        // Helper-only schema contract: the public path declines the 3-row
        // surface or missing/unclamped target knots before this comparison.
        assert!(matches!(
            coefficient_errors(&source, 900040, &surface, 900001, &target, &iso, &mut b),
            Err(ReadError::Unsupported)
        ));
    }
}
#[test]
fn actual_degree_twelve_boundary_is_refused_before_degree_eleven_extraction() {
    let original = input(false);
    let row0 = (0..=12)
        .map(|i| format!("#{}", 930000 + i))
        .collect::<Vec<_>>()
        .join(",");
    let row1 = (0..=12)
        .map(|i| format!("#{}", 940000 + i))
        .collect::<Vec<_>>()
        .join(",");
    let text = record(
        &original,
        900001,
        &format!(
            "#900001=B_SPLINE_CURVE_WITH_KNOTS('',12,({row0}),.UNSPECIFIED.,.F.,.F.,(13,13),(0.,6.323173064558),.UNSPECIFIED.);"
        ),
    );
    let text = record(
        &text,
        900040,
        &format!(
            "#900040=B_SPLINE_SURFACE_WITH_KNOTS('',1,12,(({row0}),({row1})),.UNSPECIFIED.,.F.,.F.,.F.,(2,2),(13,13),(0.,1.),(0.,6.323173064558),.UNSPECIFIED.);"
        ),
    );
    // Binomial p12 controls (i,i²) give C(.5)=(6,39) exactly.
    let text = text.replace(
        "#900021=CARTESIAN_POINT('',(0.5,0.25,0.))",
        "#900021=CARTESIAN_POINT('',(6.,39.,0.))",
    );
    let mut points = String::new();
    for i in 0..=12 {
        points += &format!(
            "#{}=CARTESIAN_POINT('',({i}.,{}.,0.));\n#{}=CARTESIAN_POINT('',({i}.,{}.,1.));\n",
            930000 + i,
            i * i,
            940000 + i,
            i * i
        );
    }
    let at = text.rfind("ENDSEC;").unwrap();
    let text = format!("{}{}{}", &text[..at], points, &text[at..]);
    let source = Source::parse(&text, 1e-7).unwrap();
    let mut b = Budget::new(512_000_000, None);
    let (d, s, e) = pair(&source, &mut b);
    assert!(matches!(&e.curve,Curve::Nurbs{curve,..}if curve.degree()==12));
    assert!(
        certify(&source, 900001, &d, &s, &e, 1e-7, &mut b)
            .unwrap()
            .is_none()
    );
}
