use super::*;
const BASE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../multicad-app/tests/fixtures/exchange/corpus/occt-cone-frustum.step"
));
fn input(end: &str) -> String {
    let extra = format!(
        "\n#910001=B_SPLINE_CURVE_WITH_KNOTS('',2,(#910010,#910011,#910012),.UNSPECIFIED.,.F.,.F.,(3,3),(0.,{end}),.UNSPECIFIED.);\n#910010=CARTESIAN_POINT('',(0.,0.,0.));\n#910011=CARTESIAN_POINT('',(0.5,0.,0.));\n#910012=CARTESIAN_POINT('',(1.,1.,0.));\n#910020=CARTESIAN_POINT('',(0.,0.,0.));\n#910021=CARTESIAN_POINT('',(0.5,0.25,0.));\n#910030=VERTEX_POINT('',#910020);\n#910031=VERTEX_POINT('',#910021);\n#910100=PCURVE('',#31,#910101);\n#910101=DEFINITIONAL_REPRESENTATION('',(#910102),#910103);\n#910103=(GEOMETRIC_REPRESENTATION_CONTEXT(2)PARAMETRIC_REPRESENTATION_CONTEXT()REPRESENTATION_CONTEXT('',''));\n#910102=B_SPLINE_CURVE_WITH_KNOTS('',1,(#910110,#910111),.UNSPECIFIED.,.F.,.F.,(2,2),(0.,3.161586532279),.UNSPECIFIED.);\n#910110=CARTESIAN_POINT('',(0.,0.));\n#910111=CARTESIAN_POINT('',(1.,3.161586532279));\n"
    );
    let at = BASE.rfind("ENDSEC;").unwrap();
    format!("{}{}{}", &BASE[..at], extra, &BASE[at..])
}
#[test]
fn original_nondyadic_natural_domain_is_mapped_not_identified_with_its_cache() {
    for text in ["6.323173064558", "6.3231730645580000000001"] {
        let input = input(text);
        let source = Source::parse(&input, 1e-7).unwrap();
        let mut b = Budget::new(512_000_000, None);
        let d = super::super::super::uv::definition(&source, 910100, &mut b).unwrap();
        let curve = source.curve(910001, 0).unwrap();
        let Curve::Nurbs { curve: ref c, .. } = curve else {
            panic!()
        };
        let data =
            super::super::bezier::original_curve(&source, 910001, c, Coordinates::Raw(3), &mut b)
                .unwrap();
        let e = Edge {
            curve,
            vertices: [910030, 910031],
            range: [0., 3.161586532279],
        };
        assert!(
            matches!(
                prepare(&source, 910001, &d, &e, &data, &mut b),
                Err(ReadError::Unsupported)
            ),
            "legacy exact-domain contract changed"
        );
        let mut b = Budget::new(512_000_000, None);
        let w = prepare_normalized(&source, 910001, &d, &e, &data, &mut b).unwrap();
        assert!(w.matches(&source, 910001, &d, &e));
        assert!(w.gamma() > 0.15);
        assert!(w.closest_delta() < 1e-11);
        // Exact original map of cache-q=.5 is B_raw/2. The decimal shadow may
        // round to the SAME cache; its exact nonidentity still lies in the bound.
        let digits = text.replace('.', "");
        let fraction = text.split_once('.').unwrap().1.len() as u32;
        let raw = BigRational::new(digits.parse().unwrap(), BigInt::from(10).pow(fraction));
        let half = &raw / BigRational::from_integer(2.into());
        let uv_end = BigRational::new(3161586532279_i64.into(), 1_000_000_000_000_i64.into());
        let delta = if half >= uv_end {
            half - uv_end
        } else {
            uv_end - half
        };
        assert!(BigRational::from_float(w.cached_window_delta()).unwrap() >= delta);
        let used = b.used();
        assert!(matches!(
            prepare_normalized(
                &source,
                910001,
                &d,
                &e,
                &data,
                &mut Budget::new(used - 1, None)
            ),
            Err(ReadError::Limit)
        ));
        assert!(
            prepare_normalized(&source, 910001, &d, &e, &data, &mut Budget::new(used, None))
                .is_ok()
        );
        let foreign = Source::parse(&input, 1e-7).unwrap();
        assert!(!w.matches(&foreign, 910001, &d, &e));
    }
}
#[test]
fn normalized_finite_reference_cannot_leave_the_cached_natural_domain() {
    let input = input("6.323173064558");
    let source = Source::parse(&input, 1e-7).unwrap();
    let mut b = Budget::new(512_000_000, None);
    let d = super::super::super::uv::definition(&source, 910100, &mut b).unwrap();
    let curve = source.curve(910001, 0).unwrap();
    let Curve::Nurbs { curve: ref c, .. } = curve else {
        panic!()
    };
    let data =
        super::super::bezier::original_curve(&source, 910001, c, Coordinates::Raw(3), &mut b)
            .unwrap();
    let edge = Edge {
        curve,
        vertices: [910030, 910031],
        range: [-1e-20, 3.161586532279],
    };
    assert!(matches!(
        prepare_normalized(&source, 910001, &d, &edge, &data, &mut b),
        Err(ReadError::Unsupported)
    ));
}
