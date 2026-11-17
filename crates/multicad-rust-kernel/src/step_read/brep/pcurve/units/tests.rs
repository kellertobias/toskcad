use super::*;
const INPUT: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../multicad-app/tests/fixtures/exchange/corpus/occt-cone-frustum.step"
));
const TOL: f64 = 1e-7;
fn fixture(factor: &str, angle: &str, turn: &str) -> String {
    let old = "#115 = ( NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.) );";
    assert_eq!(INPUT.matches(old).count(), 1);
    INPUT.replace(old, &format!(
        "#115 = ( CONVERSION_BASED_UNIT('DEGREE',#119) NAMED_UNIT(*) PLANE_ANGLE_UNIT() );\n#119 = PLANE_ANGLE_MEASURE_WITH_UNIT(PLANE_ANGLE_MEASURE({factor}),#120);\n#120 = ( NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.) );"
    )).replace("8.,0.380506377112)", &format!("8.,{angle})"))
      .replace("(-6.28318530718,-0.)", &format!("({turn},-0.)"))
}

#[test]
fn actual_converted_charts_pass_and_unchanged_angular_seam_names_its_disagreement() {
    for (factor, angle, turn) in [
        ("1.745329251994E-02", "21.8014094863518", "-360."),
        ("1.570796326795E-02", "24.2237883181687", "-400."),
    ] {
        let good = fixture(factor, angle, turn);
        let report = validate_pcurves(&good, TOL).unwrap();
        assert_eq!(report.pcurves_checked, 6);
        assert!(report.certified_maximum_deviation_mm <= TOL);
        let stale = fixture(factor, angle, "-6.28318530718");
        let error = validate_pcurves(&stale, TOL).unwrap_err();
        assert_eq!(error.kind, ReadError::Invalid);
        assert_eq!(error.disagreements.len(), 1);
        let d = &error.disagreements[0];
        assert_eq!(
            (d.surface_entity, d.edge_entity, d.pcurve_entity),
            (43, 59, 74)
        );
        assert!(d.deviation_mm > 0.1);
    }
}

#[test]
fn a_near_degree_unit_is_used_literally_and_cannot_heal_inconsistent_geometry() {
    let input = fixture("1.745329301994E-02", "21.8014094863518", "-360.");
    let source = Source::parse(&input, TOL).unwrap();
    let (p, _) = source.surface(43).unwrap().evaluate([0., -10.]);
    let expected = 8. - 10. * (21.8014094863518_f64 * 1.745329301994E-02).tan();
    assert!((p[0] - expected).abs() < 1e-12);
    assert!((p[0] - 4.).abs() > TOL);
    let error = validate_pcurves(&input, TOL).unwrap_err();
    assert_eq!(error.kind, ReadError::Invalid);
    assert!(error.disagreements.iter().any(|d| d.deviation_mm > TOL));
}

#[test]
fn wrong_unit_bases_underflow_directions_and_reverse_sense_never_earn_a_certificate() {
    let good = fixture("1.745329251994E-02", "21.8014094863518", "-360.");
    for bad in [
        good.replace(
            "#120 = ( NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.) );",
            "#120 = ( NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.METRE.) );",
        ),
        good.replace(
            "#79 = DIRECTION('',(-0.,-1.));",
            "#79 = DIRECTION('',(1E-1000,-1.));",
        ),
        good.replace(
            "#79 = DIRECTION('',(-0.,-1.));",
            "#79 = DIRECTION('',(-0.,1.));",
        ),
    ] {
        assert_ne!(good, bad);
        assert!(validate_pcurves(&bad, TOL).is_err());
    }
    let cancelled = validate_pcurves_with_interrupt(&good, TOL, Some(&|| true)).unwrap_err();
    assert!(cancelled.cancelled);
    assert!(cancelled.disagreements.is_empty());
}
