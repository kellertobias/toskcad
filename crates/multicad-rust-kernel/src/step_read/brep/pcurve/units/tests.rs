use super::*;
use multicad_precision::PortableMath as _;
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
    let expected = 8. - 10. * (21.8014094863518_f64 * 1.745329301994E-02).portable_tan();
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

/// A proper cyclic rotation R(x,y,z)=(y,z,x), followed by (11,-7,5).
/// Every actual 3D point/direction record moves; all original 2D charts stay
/// unchanged. Original direction decimal spellings are merely permuted.
fn posed(input: &str) -> String {
    let source = Source::parse(input, TOL).unwrap();
    let mut edits = Vec::new();
    for parts in source.entities.values() {
        let [entity] = parts.as_slice() else { continue };
        if !matches!(entity.name, "CARTESIAN_POINT" | "DIRECTION") {
            continue;
        }
        let [_, coordinates] = entity.args.as_slice() else {
            panic!()
        };
        let values = list(coordinates).unwrap();
        if values.len() != 3 {
            continue;
        }
        let original = entity.source.unwrap();
        let start = original.rfind('(').unwrap() + 1;
        let end = original[start..].find(')').unwrap() + start;
        let coordinates: Vec<_> = original[start..end].split(',').map(str::trim).collect();
        assert_eq!(coordinates.len(), 3);
        let ordered = [coordinates[1], coordinates[2], coordinates[0]];
        let changed = if entity.name == "DIRECTION" {
            ordered.join(",")
        } else {
            // All 3D points in this actual fixture have exact integer leaves.
            // Thus this decimal rewrite is an exact rigid transform, not an
            // f64 reevaluation of the original support or an author recipe.
            ordered
                .into_iter()
                .zip([11_i64, -7, 5])
                .map(|(value, shift)| {
                    let n = value.parse::<f64>().unwrap();
                    assert!(n.is_finite() && n.fract() == 0. && n.abs() < 1e6);
                    format!("{}.", n as i64 + shift)
                })
                .collect::<Vec<_>>()
                .join(",")
        };
        let replacement = format!("{}{}{}", &original[..start], changed, &original[end..]);
        let offset = original.as_ptr() as usize - input.as_ptr() as usize;
        assert_eq!(&input[offset..offset + original.len()], original);
        edits.push((offset, original.len(), replacement));
    }
    edits.sort_by_key(|e| e.0);
    let mut output = input.to_owned();
    for (offset, length, replacement) in edits.into_iter().rev() {
        output.replace_range(offset..offset + length, &replacement);
    }
    output
}

#[test]
fn actual_proper_pose_preserves_converted_raw_cone_whole_trace_and_seam_rejection() {
    let input = fixture("1.745329251994E-02", "21.8014094863518", "-360.");
    let posed = posed(&input);
    assert_ne!(input, posed);
    let original = validate_pcurves(&input, TOL).unwrap();
    let report = validate_pcurves(&posed, TOL).unwrap();
    assert_eq!(report.pcurves_checked, original.pcurves_checked);
    assert_eq!(report.edges_checked, original.edges_checked);
    assert_eq!(report.pcurves_checked, 6);
    assert!(report.certified_maximum_deviation_mm <= TOL);
    // Independent metrology checks the actual transformed support; the
    // positive authority above remains the complete OriginalREAL proof.
    let source = Source::parse(&posed, TOL).unwrap();
    let (point, _) = source.surface(43).unwrap().evaluate([0., -10.]);
    for (actual, expected) in point.into_iter().zip([11., 3., 9.]) {
        assert!((actual - expected).abs() < 1e-10);
    }
    let stale = posed.replace("(-360.,-0.)", "(-6.28318530718,-0.)");
    assert_ne!(stale, posed);
    let error = validate_pcurves(&stale, TOL).unwrap_err();
    assert_eq!(error.kind, ReadError::Invalid);
    assert!(
        error
            .disagreements
            .iter()
            .any(
                |d| (d.surface_entity, d.edge_entity, d.pcurve_entity) == (43, 59, 74)
                    && d.deviation_mm > 0.1
            )
    );
}

#[test]
fn actual_original_converted_seam_has_exact_work_frontier_late_cancel_and_retry() {
    let input = posed(&fixture("1.745329251994E-02", "21.8014094863518", "-360."));
    let source = Source::parse(&input, TOL).unwrap();
    // Preparation is separate from the actual leaf contract under test.
    // Its inputs are earned from the unchanged public parser/binding path.
    let mut preparation = Budget::new(MAX_WORK, None);
    let binding = bindings(&source, &mut preparation)
        .unwrap()
        .into_iter()
        .find(|b| b.pcurve == 74)
        .unwrap();
    let edge = edge(&source, binding, &mut preparation).unwrap();
    let definition = uv::definition(&source, binding.pcurve, &mut preparation).unwrap();
    let surface = source.surface(definition.surface).unwrap();
    let run = |budget: &mut Budget<'_>| {
        certify(
            &source,
            binding.curve,
            &definition,
            &surface,
            &edge,
            TOL,
            budget,
        )
    };
    let mut full = Budget::new(MAX_WORK, None);
    let Some(Proof::Certified(first)) = run(&mut full).unwrap() else {
        panic!()
    };
    assert!(first <= TOL);
    let used = full.used();
    assert!(used > 16384 && used < MAX_WORK);
    let mut exact = Budget::new(used, None);
    let Some(Proof::Certified(same)) = run(&mut exact).unwrap() else {
        panic!()
    };
    assert_eq!(same, first);
    assert_eq!(exact.remaining, 0);
    assert!(matches!(
        run(&mut Budget::new(used - 1, None)),
        Err(ReadError::Limit)
    ));
    let calls = std::cell::Cell::new(0);
    let interrupt = || {
        calls.set(calls.get() + 1);
        calls.get() == 20
    };
    let mut cancelled = Budget::new(MAX_WORK, Some(&interrupt));
    assert!(matches!(run(&mut cancelled), Err(ReadError::Limit)));
    assert!(cancelled.cancelled.get());
    assert_eq!(calls.get(), 20);
    assert!(cancelled.used() > 0);
    let mut retry = Budget::new(MAX_WORK, None);
    let Some(Proof::Certified(after)) = run(&mut retry).unwrap() else {
        panic!()
    };
    assert_eq!(after, first);
    assert_eq!(retry.used(), used);
    // A failed private attempt never changes the actual input or public result.
    assert_eq!(validate_pcurves(&input, TOL).unwrap().pcurves_checked, 6);
}
