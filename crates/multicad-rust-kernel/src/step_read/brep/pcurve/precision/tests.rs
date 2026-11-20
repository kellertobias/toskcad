use super::*;
use crate::brep::Edge;
use crate::step_read::brep::pcurve::{Proof, affine, sweep, uv};
const TOL: f64 = 1e-7;
const INPUT: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/src/step_read/tests/fixtures/occt-spline-prism.step"
));
fn budget() -> Budget<'static> {
    Budget::new(512_000_000, None)
}
fn leaf<'a, 'i>(source: &'a Source<'i>, id: u32, index: usize) -> &'a Value<'i> {
    let [_, values] = source.args(id, "CARTESIAN_POINT").unwrap() else {
        panic!()
    };
    &list(values).unwrap()[index]
}
fn pair(source: &Source<'_>, pc: u32, edge_id: u32) -> (uv::Definition, Surface, Edge, u32) {
    let definition = uv::definition(source, pc, &mut budget()).unwrap();
    let args = source.args(edge_id, "EDGE_CURVE").unwrap();
    let points = [&args[1], &args[2]].map(|v| {
        let args = source.args(reference(v).unwrap(), "VERTEX_POINT").unwrap();
        source.point(reference(&args[1]).unwrap()).unwrap()
    });
    let id = reference(&args[3]).unwrap();
    let curve = source.curve(id, 0).unwrap();
    (
        uv::Definition {
            surface: definition.surface,
            curve_id: definition.curve_id,
            geometry: definition.geometry,
        },
        source.surface(definition.surface).unwrap(),
        Edge::new(curve, [0, 1], points, TOL).unwrap(),
        id,
    )
}
#[test]
fn original_numeric_binding_distinguishes_equal_cache_aliases_and_foreign_arenas() {
    let original = Source::parse(INPUT, TOL).unwrap();
    let value = scalar(&original, 46, leaf(&original, 46, 0), &mut budget()).unwrap();
    assert_eq!((value.lo, value.hi), (20., 20.));
    let changed = INPUT.replace(
        "#46 = CARTESIAN_POINT('',(20.,0.,0.));",
        "#46 = CARTESIAN_POINT('',(20.0000000000000001,0.,0.));",
    );
    assert_ne!(changed, INPUT);
    let changed = Source::parse(&changed, TOL).unwrap();
    let value = scalar(&changed, 46, leaf(&changed, 46, 0), &mut budget()).unwrap();
    assert!(value.lo < 20. && value.hi > 20.);
    let raw = rational(&changed, 46, leaf(&changed, 46, 0), &mut budget()).unwrap();
    assert!(
        BigRational::from_float(value.lo).unwrap() < raw
            && raw < BigRational::from_float(value.hi).unwrap()
    );
    assert!(matches!(
        scalar(&original, 46, leaf(&changed, 46, 0), &mut budget()),
        Err(ReadError::Invalid)
    ));
    let tiny = INPUT.replace(
        "#46 = CARTESIAN_POINT('',(20.,0.,0.));",
        "#46 = CARTESIAN_POINT('',(20.,1E-1000,0.));",
    );
    let tiny = Source::parse(&tiny, TOL).unwrap();
    assert!(matches!(
        scalar(&tiny, 46, leaf(&tiny, 46, 1), &mut budget()),
        Err(ReadError::Unsupported)
    ));
    assert_eq!(
        scalar(&original, 46, leaf(&original, 46, 1), &mut budget())
            .unwrap()
            .hi,
        0.
    );
}
#[test]
fn actual_frames_earn_cardinal_or_conditioned_real_axes_and_refuse_near_parallel_uncertainty() {
    let original = Source::parse(INPUT, TOL).unwrap();
    let f = frame(&original, 33, 3, &mut budget()).unwrap();
    for v in f.axes.into_iter().flatten() {
        assert_eq!(v.lo, v.hi);
        assert!([-1., 0., 1.].contains(&v.lo));
    }
    let changed = INPUT
        .replace(
            "#35 = DIRECTION('',(0.,1.,0.));",
            "#35 = DIRECTION('',(0.6,0.8,0.));",
        )
        .replace(
            "#36 = DIRECTION('',(1.,0.,0.));",
            "#36 = DIRECTION('',(0.,0.,1.));",
        );
    let actual = Source::parse(&changed, TOL).unwrap();
    let f = frame(&actual, 33, 3, &mut budget()).unwrap();
    let expected = [[0., 0., 1.], [0.8, -0.6, 0.], [0.6, 0.8, 0.]];
    for (a, e) in f.axes.into_iter().zip(expected) {
        for (a, e) in a.into_iter().zip(e) {
            let q = Rational::parse_decimal(&e.to_string()).unwrap();
            let (n, d) = q.ratio();
            let q = BigRational::new(BigInt::from(n), BigInt::from(d));
            assert!(
                BigRational::from_float(a.lo).unwrap() <= q
                    && q <= BigRational::from_float(a.hi).unwrap()
            );
        }
    }
    let near = changed.replace(
        "#36 = DIRECTION('',(0.,0.,1.));",
        "#36 = DIRECTION('',(0.6000000000000001,0.8,0.));",
    );
    let near = Source::parse(&near, TOL).unwrap();
    assert!(matches!(
        frame(&near, 33, 3, &mut budget()),
        Err(ReadError::Unsupported)
    ));
    assert!(frame(&original, 33, 3, &mut budget()).is_ok());
}
#[test]
fn all_actual_extrusion_and_planar_cubic_pairs_keep_the_original_limit() {
    let source = Source::parse(INPUT, TOL).unwrap();
    for (pc, id) in [(43, 21), (102, 87), (147, 140), (171, 164), (154, 140)] {
        let (definition, surface, edge, curve) = pair(&source, pc, id);
        let proof = if matches!(surface, Surface::Extrusion { .. }) {
            sweep::certify(&surface, &definition.geometry, &edge, TOL, &mut budget()).unwrap()
        } else {
            affine::certify(&surface, &definition.geometry, &edge, TOL, &mut budget()).unwrap()
        };
        let Proof::Certified(upper) = proof else {
            panic!("decoded actual pair unproved")
        };
        let raw = reserve(
            &source,
            curve,
            &definition,
            &surface,
            &edge,
            &mut budget(),
            &mut Cache::new(&source),
        )
        .unwrap();
        sense(&source, curve, &definition, &surface, &edge, &mut budget()).unwrap();
        assert!(
            raw.is_finite() && raw >= 0. && (I::point(raw) + I::point(upper)).hi < TOL,
            "pc={pc} numeric={raw} decoded={upper}"
        );
    }
    // Equivalent ORIGINAL open tensor support, actual rows share the retained
    // cubic edge poles. No inferred authored recipe or triangle surface.
    let tensor=INPUT.replace("#44 = SURFACE_OF_LINEAR_EXTRUSION('',#45,#50);","#44 = B_SPLINE_SURFACE_WITH_KNOTS('',3,1,((#143,#46),(#144,#47),(#145,#48),(#146,#49)),.UNSPECIFIED.,.F.,.F.,.F.,(4,4),(2,2),(0.,1.),(-6.,0.),.UNSPECIFIED.);");
    assert_ne!(tensor, INPUT);
    let source = Source::parse(&tensor, TOL).unwrap();
    for (pc, id) in [(147, 140), (171, 164)] {
        let (definition, surface, edge, curve) = pair(&source, pc, id);
        assert!(matches!(
            affine::certify(&surface, &definition.geometry, &edge, TOL, &mut budget()).unwrap(),
            Proof::Certified(_)
        ));
        assert!(
            reserve(
                &source,
                curve,
                &definition,
                &surface,
                &edge,
                &mut budget(),
                &mut Cache::new(&source)
            )
            .unwrap()
                < TOL
        );
        sense(&source, curve, &definition, &surface, &edge, &mut budget()).unwrap();
    }
}
#[test]
fn raw_knot_boundary_aliases_and_noninteger_degree_aliases_do_not_become_original_authority() {
    let changed = INPUT.replacen(
        "(4,4),(0.,1.),.PIECEWISE_BEZIER_KNOTS.",
        "(4,4),(0.,1.0000000000000001),.PIECEWISE_BEZIER_KNOTS.",
        1,
    );
    assert_ne!(changed, INPUT);
    let source = Source::parse(&changed, TOL).unwrap();
    let (definition, surface, edge, curve) = pair(&source, 43, 21);
    assert!(matches!(
        reserve(
            &source,
            curve,
            &definition,
            &surface,
            &edge,
            &mut budget(),
            &mut Cache::new(&source)
        ),
        Err(ReadError::Unsupported)
    ));
    let changed = INPUT.replace(
        "#142 = B_SPLINE_CURVE_WITH_KNOTS('',3,",
        "#142 = B_SPLINE_CURVE_WITH_KNOTS('',3.0000000000000001,",
    );
    assert_ne!(changed, INPUT);
    let source = Source::parse(&changed, TOL).unwrap();
    let decoded = source.curve(142, 0).unwrap();
    assert!(matches!(
        curve_nurbs(&source, 142, &decoded, 3, &mut budget()),
        Err(ReadError::Unsupported)
    ));
}
#[test]
fn original_decimal_affine_cancellation_preserves_large_gauge_and_reports_real_offset() {
    // The physical resulting point stays small; the original plane origin
    // and the original UV offset each have exact long decimal authority.
    let changed = INPUT
        .replace(
            "#34 = CARTESIAN_POINT('',(0.,0.,0.));",
            "#34 = CARTESIAN_POINT('',(1000000000.123456789,0.,-2000000000.125));",
        )
        .replace(
            "#39 = CARTESIAN_POINT('',(0.,0.));",
            "#39 = CARTESIAN_POINT('',(-1000000000.123456789,-2000000000.125));",
        );
    let source = Source::parse(&changed, TOL).unwrap();
    let f = frame(&source, 33, 3, &mut budget()).unwrap();
    let p = translated_point(&source, 34, 39, f.axes, &mut budget()).unwrap();
    for x in p {
        assert_eq!((x.lo, x.hi), (0., 0.));
    }
    let wrong = changed.replace(
        "(-1000000000.123456789,-2000000000.125)",
        "(-1000000000.123456788,-2000000000.125)",
    );
    assert_ne!(wrong, changed);
    let source = Source::parse(&wrong, TOL).unwrap();
    let f = frame(&source, 33, 3, &mut budget()).unwrap();
    let p = translated_point(&source, 34, 39, f.axes, &mut budget()).unwrap();
    let q = BigRational::new(BigInt::from(1), BigInt::from(1_000_000_000));
    assert!(
        p[0].lo > 0.
            && BigRational::from_float(p[0].lo).unwrap() <= q
            && q <= BigRational::from_float(p[0].hi).unwrap()
    );
}
#[test]
fn actual_numeric_proof_work_and_observed_cancellation_frontiers_preserve_retry() {
    use std::cell::Cell;
    let source = Source::parse(INPUT, TOL).unwrap();
    let (definition, surface, edge, curve) = pair(&source, 147, 140);
    let calls = Cell::new(0);
    let callback = || {
        calls.set(calls.get() + 1);
        false
    };
    let mut b = Budget::new(512_000_000, Some(&callback));
    let original = reserve(
        &source,
        curve,
        &definition,
        &surface,
        &edge,
        &mut b,
        &mut Cache::new(&source),
    )
    .unwrap();
    let (used, polls) = (b.used(), calls.get());
    assert!(used > 1024 && polls > 10);
    assert!(matches!(
        reserve(
            &source,
            curve,
            &definition,
            &surface,
            &edge,
            &mut Budget::new(used - 1, None),
            &mut Cache::new(&source)
        ),
        Err(ReadError::Limit)
    ));
    assert_eq!(
        reserve(
            &source,
            curve,
            &definition,
            &surface,
            &edge,
            &mut Budget::new(used, None),
            &mut Cache::new(&source)
        )
        .unwrap()
        .to_bits(),
        original.to_bits()
    );
    for stop in [1, polls / 2, polls] {
        let seen = Cell::new(0);
        let callback = || {
            seen.set(seen.get() + 1);
            seen.get() >= stop
        };
        assert!(matches!(
            reserve(
                &source,
                curve,
                &definition,
                &surface,
                &edge,
                &mut Budget::new(512_000_000, Some(&callback)),
                &mut Cache::new(&source)
            ),
            Err(ReadError::Limit)
        ));
        assert!(seen.get() >= stop);
        assert_eq!(
            reserve(
                &source,
                curve,
                &definition,
                &surface,
                &edge,
                &mut budget(),
                &mut Cache::new(&source)
            )
            .unwrap()
            .to_bits(),
            original.to_bits()
        );
    }
}

#[test]
fn direct_original_plane_cubic_certificate_keeps_large_legal_gauge_and_rejects_real_shift() {
    let gauge = INPUT
        .replace(
            "#77 = CARTESIAN_POINT('',(0.,0.,6.));",
            "#77 = CARTESIAN_POINT('',(1000000000.123456789,0.,6.));",
        )
        .replace(
            "#157 = CARTESIAN_POINT('',(20.,0.));",
            "#157 = CARTESIAN_POINT('',(-999999980.123456789,0.));",
        )
        .replace(
            "#158 = CARTESIAN_POINT('',(20.,12.));",
            "#158 = CARTESIAN_POINT('',(-999999980.123456789,12.));",
        )
        .replace(
            "#159 = CARTESIAN_POINT('',(0.,12.));",
            "#159 = CARTESIAN_POINT('',(-1000000000.123456789,12.));",
        )
        .replace(
            "#160 = CARTESIAN_POINT('',(0.,0.));",
            "#160 = CARTESIAN_POINT('',(-1000000000.123456789,0.));",
        );
    assert_ne!(gauge, INPUT);
    let source = Source::parse(&gauge, TOL).unwrap();
    let (definition, surface, edge, curve) = pair(&source, 154, 140);
    assert!(
        matches!(plane_nurbs_certificate(&source,curve,&definition,&surface,&edge,TOL,&mut budget()).unwrap(),Some(Proof::Certified(upper))if upper<TOL)
    );
    let wrong = gauge
        .replace("-999999980.123456789", "-999999980.023456789")
        .replace("-1000000000.123456789", "-1000000000.023456789");
    assert_ne!(wrong, gauge);
    let source = Source::parse(&wrong, TOL).unwrap();
    let (definition, surface, edge, curve) = pair(&source, 154, 140);
    assert!(
        matches!(plane_nurbs_certificate(&source,curve,&definition,&surface,&edge,TOL,&mut budget()).unwrap(),Some(Proof::Disagreement(lower))if lower>0.05)
    );
    let source = Source::parse(&gauge, TOL).unwrap();
    let (definition, surface, edge, curve) = pair(&source, 154, 140);
    assert!(matches!(
        plane_nurbs_certificate(
            &source,
            curve,
            &definition,
            &surface,
            &edge,
            TOL,
            &mut budget()
        )
        .unwrap(),
        Some(Proof::Certified(_))
    ));
}

#[test]
fn direct_original_plane_line_certificate_uses_real_compensated_origin_and_real_mismatch() {
    let gauge = INPUT
        .replace(
            "#34 = CARTESIAN_POINT('',(0.,0.,0.));",
            "#34 = CARTESIAN_POINT('',(1000000000.123456789,0.,-2000000000.125));",
        )
        .replace(
            "#39 = CARTESIAN_POINT('',(0.,0.));",
            "#39 = CARTESIAN_POINT('',(-1000000000.123456789,-2000000000.125));",
        );
    let source = Source::parse(&gauge, TOL).unwrap();
    let (definition, surface, edge, curve) = pair(&source, 31, 21);
    assert!(
        matches!(plane_analytic_certificate(&source,curve,&definition,&surface,&edge,TOL,&mut budget()).unwrap(),Some(Proof::Certified(upper))if upper<TOL)
    );
    let wrong = gauge.replace(
        "(-1000000000.123456789,-2000000000.125)",
        "(-1000000000.023456789,-2000000000.125)",
    );
    assert_ne!(wrong, gauge);
    let source = Source::parse(&wrong, TOL).unwrap();
    let (definition, surface, edge, curve) = pair(&source, 31, 21);
    assert!(
        matches!(plane_analytic_certificate(&source,curve,&definition,&surface,&edge,TOL,&mut budget()).unwrap(),Some(Proof::Disagreement(lower))if lower>0.05)
    );
}

#[test]
fn original_finite_uv_splines_cover_exact_full_isocurve_and_refuse_partial_or_aliased_endpoints() {
    let tensor=INPUT.replace("#44 = SURFACE_OF_LINEAR_EXTRUSION('',#45,#50);","#44 = B_SPLINE_SURFACE_WITH_KNOTS('',3,1,((#143,#46),(#144,#47),(#145,#48),(#146,#49)),.UNSPECIFIED.,.F.,.F.,.F.,(4,4),(2,2),(0.,1.),(-6.,0.),.UNSPECIFIED.);");
    let segment=tensor.replace("#149 = LINE('',#150,#151);","#149 = B_SPLINE_CURVE_WITH_KNOTS('',1,(#150,#900000),.UNSPECIFIED.,.F.,.F.,(2,2),(0.,1.),.UNSPECIFIED.);");
    let end = segment.rfind("ENDSEC;").unwrap();
    let segment = format!(
        "{}#900000 = CARTESIAN_POINT('',(1.,-6.));\n{}",
        &segment[..end],
        &segment[end..]
    );
    let source = Source::parse(&segment, TOL).unwrap();
    let (definition, surface, edge, curve) = pair(&source, 147, 140);
    assert!(matches!(
        affine::certify(&surface, &definition.geometry, &edge, TOL, &mut budget()).unwrap(),
        Proof::Certified(_)
    ));
    assert!(
        reserve(
            &source,
            curve,
            &definition,
            &surface,
            &edge,
            &mut budget(),
            &mut Cache::new(&source)
        )
        .unwrap()
            < TOL
    );
    sense(&source, curve, &definition, &surface, &edge, &mut budget()).unwrap();
    // Positive rational weights change only the monotone finite segment
    // parameterization; they do not turn it into an unbounded LINE.
    let rational=segment.replace("#149 = B_SPLINE_CURVE_WITH_KNOTS('',1,(#150,#900000),.UNSPECIFIED.,.F.,.F.,(2,2),(0.,1.),.UNSPECIFIED.);","#149 = (BOUNDED_CURVE() B_SPLINE_CURVE(1,(#150,#900000),.UNSPECIFIED.,.F.,.F.) B_SPLINE_CURVE_WITH_KNOTS((2,2),(0.,1.),.UNSPECIFIED.) CURVE() GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_CURVE((2.,5.)) REPRESENTATION_ITEM(''));");
    assert_ne!(rational, segment);
    let source = Source::parse(&rational, TOL).unwrap();
    let (definition, surface, edge, curve) = pair(&source, 147, 140);
    assert!(matches!(
        affine::certify(&surface, &definition.geometry, &edge, TOL, &mut budget()).unwrap(),
        Proof::Certified(_)
    ));
    assert!(
        reserve(
            &source,
            curve,
            &definition,
            &surface,
            &edge,
            &mut budget(),
            &mut Cache::new(&source)
        )
        .unwrap()
            < TOL
    );
    for replacement in [
        "(0.5,-6.)",
        "(1.0000000000000001,-6.)",
        "(1.,-6.0000000000000001)",
    ] {
        let changed = segment.replace("(1.,-6.));", &format!("{replacement});"));
        assert_ne!(changed, segment);
        let source = Source::parse(&changed, TOL).unwrap();
        let (definition, surface, edge, curve) = pair(&source, 147, 140);
        assert!(matches!(
            reserve(
                &source,
                curve,
                &definition,
                &surface,
                &edge,
                &mut budget(),
                &mut Cache::new(&source)
            ),
            Err(ReadError::Unsupported)
        ));
    }
}

#[test]
fn invocation_surface_cache_reuses_only_whole_same_source_proof_and_rechecks_each_finite_trim() {
    let tensor=INPUT.replace("#44 = SURFACE_OF_LINEAR_EXTRUSION('',#45,#50);","#44 = B_SPLINE_SURFACE_WITH_KNOTS('',3,1,((#143,#46),(#144,#47),(#145,#48),(#146,#49)),.UNSPECIFIED.,.F.,.F.,.F.,(4,4),(2,2),(0.,1.),(-6.,0.),.UNSPECIFIED.);");
    let partial=tensor.replace("#173 = LINE('',#174,#175);","#173 = B_SPLINE_CURVE_WITH_KNOTS('',1,(#174,#900001),.UNSPECIFIED.,.F.,.F.,(2,2),(0.,1.),.UNSPECIFIED.);");
    let end = partial.rfind("ENDSEC;").unwrap();
    let partial = format!(
        "{}#900001 = CARTESIAN_POINT('',(0.5,0.));\n{}",
        &partial[..end],
        &partial[end..]
    );
    let source = Source::parse(&partial, TOL).unwrap();
    let (definition, surface, edge, curve) = pair(&source, 147, 140);
    let mut cache = Cache::new(&source);
    let mut cold = budget();
    let first = reserve(
        &source,
        curve,
        &definition,
        &surface,
        &edge,
        &mut cold,
        &mut cache,
    )
    .unwrap();
    let mut warm = budget();
    let second = reserve(
        &source,
        curve,
        &definition,
        &surface,
        &edge,
        &mut warm,
        &mut cache,
    )
    .unwrap();
    assert_eq!(first.to_bits(), second.to_bits());
    assert!(warm.used() < cold.used());
    assert!(cache.retained_bytes() <= cache.peak_bytes() && cache.peak_bytes() <= CACHE_PEAK_BYTES);
    assert_eq!(cache.surfaces.len(), 1);
    // This genuine second trim is incomplete in the SAME immutable source.
    // Its surface's cached uniform bound cannot replace finite UV coverage.
    let (bad_definition, bad_surface, bad_edge, bad_curve) = pair(&source, 171, 164);
    assert!(matches!(
        reserve(
            &source,
            bad_curve,
            &bad_definition,
            &bad_surface,
            &bad_edge,
            &mut budget(),
            &mut cache
        ),
        Err(ReadError::Unsupported)
    ));
    // Identical text in another actual parser owner does not own this cache.
    let foreign = Source::parse(&partial, TOL).unwrap();
    let (d, s, e, c) = pair(&foreign, 147, 140);
    assert!(matches!(
        reserve(&foreign, c, &d, &s, &e, &mut budget(), &mut cache),
        Err(ReadError::Invalid)
    ));
    let mut own = Cache::new(&foreign);
    assert_eq!(
        reserve(&foreign, c, &d, &s, &e, &mut budget(), &mut own)
            .unwrap()
            .to_bits(),
        first.to_bits()
    );
    assert_eq!(
        reserve(
            &source,
            curve,
            &definition,
            &surface,
            &edge,
            &mut budget(),
            &mut cache
        )
        .unwrap()
        .to_bits(),
        first.to_bits()
    );
}
#[test]
fn original_tiny_converter_coordinates_compare_exactly_without_cached_aliases() {
    for (a, b, equal) in [
        ("8.881784197001E-16", "8.881784197001E-16", true),
        ("1.110223024625E-015", "1.1102230246250E-15", true),
        ("2.287302224786E-033", "2.2873022247860E-33", true),
        (
            "8.881784197001E-16",
            "8.881784197001000000000001E-16",
            false,
        ),
    ] {
        let text = INPUT
            .replace(
                "#46 = CARTESIAN_POINT('',(20.,0.,0.));",
                &format!("#46 = CARTESIAN_POINT('',({a},0.,0.));"),
            )
            .replace(
                "#47 = CARTESIAN_POINT('',(20.,12.,0.));",
                &format!("#47 = CARTESIAN_POINT('',({b},12.,0.));"),
            );
        let source = Source::parse(&text, TOL).unwrap();
        assert_eq!(source.point(46).unwrap()[0], source.point(47).unwrap()[0]);
        let mut first = budget();
        assert_eq!(
            exact_same(
                &source,
                (46, leaf(&source, 46, 0)),
                (47, leaf(&source, 47, 0)),
                &mut first
            )
            .unwrap(),
            equal
        );
        let used = first.used();
        assert!(matches!(
            exact_same(
                &source,
                (46, leaf(&source, 46, 0)),
                (47, leaf(&source, 47, 0)),
                &mut Budget::new(used - 1, None)
            ),
            Err(ReadError::Limit)
        ));
        assert_eq!(
            exact_same(
                &source,
                (46, leaf(&source, 46, 0)),
                (47, leaf(&source, 47, 0)),
                &mut Budget::new(used, None)
            )
            .unwrap(),
            equal
        );
        let foreign = Source::parse(&text, TOL).unwrap();
        assert!(matches!(
            exact_same(
                &source,
                (46, leaf(&source, 46, 0)),
                (47, leaf(&foreign, 47, 0)),
                &mut budget()
            ),
            Err(ReadError::Invalid)
        ));
    }
    let text = INPUT.replace(
        "#46 = CARTESIAN_POINT('',(20.,0.,0.));",
        "#46 = CARTESIAN_POINT('',(1E-1000,0.,0.));",
    );
    let source = Source::parse(&text, TOL).unwrap();
    assert!(matches!(
        exact_same(
            &source,
            (46, leaf(&source, 46, 0)),
            (46, leaf(&source, 46, 0)),
            &mut budget()
        ),
        Err(ReadError::Unsupported)
    ));
}
#[test]
fn numeric_owner_scan_debits_even_empty_complex_components() {
    let original = Source::parse(INPUT, TOL).unwrap();
    let old = original.simple(45).unwrap().source.unwrap();
    let body = "(BOUNDED_CURVE() B_SPLINE_CURVE(3,(#46,#47,#48,#49),.UNSPECIFIED.,.F.,.F.) B_SPLINE_CURVE_WITH_KNOTS((4,4),(0.,1.),.PIECEWISE_BEZIER_KNOTS.) CURVE() GEOMETRIC_REPRESENTATION_ITEM() REPRESENTATION_ITEM(''))";
    let plain = INPUT.replacen(old, body, 1);
    let extra = plain.replacen("(BOUNDED_CURVE()", "(AUDIT_COMPONENT() BOUNDED_CURVE()", 1);
    assert_ne!(plain, extra);
    let mut costs = Vec::new();
    for input in [&plain, &extra] {
        // Both owners come from the actual neutral parser; the extra unknown
        // empty component grants no schema or geometry permission.
        let source = Source::parse(input, TOL).unwrap();
        let fields = source
            .component(45, "B_SPLINE_CURVE_WITH_KNOTS")
            .unwrap()
            .unwrap();
        let value = &list(&fields[1]).unwrap()[1];
        let mut work = budget();
        assert_eq!(spelling(&source, 45, value, &mut work).unwrap(), "1.");
        costs.push(work.used());
    }
    assert_eq!(costs[1], costs[0] + 16);
}

fn numeric_label_owner(name: &str, values: &str, numbers: usize) -> String {
    let label = std::iter::repeat_n("1.", numbers)
        .collect::<Vec<_>>()
        .join(",");
    let record = format!("\n#990001={name}(TYPED(({label})),{values});\n");
    let at = INPUT.rfind("ENDSEC;").unwrap();
    format!("{}{}{}", &INPUT[..at], record, &INPUT[at..])
}

#[test]
fn small_numeric_owners_bound_typed_label_tokens_before_allocation() {
    for (name, values) in [
        ("CARTESIAN_POINT", "(1.,2.,3.)"),
        ("DIRECTION", "(0.,0.,1.)"),
        ("VECTOR", "#1,2."),
        ("CIRCLE", "#1,2."),
        ("ELLIPSE", "#1,2.,1."),
        ("CYLINDRICAL_SURFACE", "#1,2."),
        ("CONICAL_SURFACE", "#1,2.,0.25"),
        ("SPHERICAL_SURFACE", "#1,2."),
        ("TOROIDAL_SURFACE", "#1,3.,1."),
    ] {
        let text = numeric_label_owner(name, values, 2048);
        let source = Source::parse(&text, TOL).unwrap();
        let args = source.args(990001, name).unwrap();
        let value = if matches!(name, "CARTESIAN_POINT" | "DIRECTION") {
            &list(&args[1]).unwrap()[0]
        } else {
            &args[args.len() - 1]
        };
        let mut work = budget();
        assert_eq!(
            spelling(&source, 990001, value, &mut work),
            Err(ReadError::Limit),
            "{name}"
        );
        // The complete charged literal traversal happens before token allocation.
        let entity = source.simple(990001).unwrap();
        let (mut count, mut found) = (0, None);
        let mut prefix = budget();
        prefix.charge(16).unwrap();
        locate(&entity.args, value, &mut count, &mut found, &mut prefix).unwrap();
        prefix.charge(entity.source.unwrap().len() + count).unwrap();
        assert_eq!(work.used(), prefix.used(), "{name}");
        assert!(found.is_some());
    }
}

#[test]
fn bounded_point_and_direction_labels_keep_real_reads_resource_and_cancel_contracts() {
    use std::cell::Cell;
    for name in ["CARTESIAN_POINT", "DIRECTION"] {
        let text = numeric_label_owner(name, "(1.,2.,3.)", 2045);
        let source = Source::parse(&text, TOL).unwrap();
        let args = source.args(990001, name).unwrap();
        let value = &list(&args[1]).unwrap()[0];
        let run = |b: &mut Budget<'_>| rational(&source, 990001, value, b);
        let mut full = budget();
        assert_eq!(run(&mut full).unwrap(), BigRational::from_integer(1.into()));
        let used = full.used();
        assert_eq!(run(&mut Budget::new(used - 1, None)), Err(ReadError::Limit));
        assert_eq!(
            run(&mut Budget::new(used, None)).unwrap(),
            BigRational::from_integer(1.into())
        );
        let calls = Cell::new(0);
        let count = || {
            calls.set(calls.get() + 1);
            false
        };
        run(&mut Budget::new(used, Some(&count))).unwrap();
        let last = calls.get();
        calls.set(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() == last
        };
        let mut stopped = Budget::new(used, Some(&cancel));
        assert_eq!(run(&mut stopped), Err(ReadError::Limit));
        assert!(stopped.cancelled.get());
        assert_eq!(
            run(&mut Budget::new(used, None)).unwrap(),
            BigRational::from_integer(1.into())
        );
        let ordinary = numeric_label_owner(name, "(1.,2.,3.)", 0);
        let ordinary = ordinary.replace("TYPED(())", "'ordinary label'");
        let source = Source::parse(&ordinary, TOL).unwrap();
        let actual = coordinates(&source, 990001, name, 3, &mut budget()).unwrap();
        for (a, expected) in actual.into_iter().zip([1., 2., 3.]) {
            assert!(a.lo <= expected && expected <= a.hi);
        }
    }
}

#[test]
fn cached_spellings_are_the_same_tokens_cheaper_and_restart_cold() {
    let source = Source::parse(INPUT, TOL).unwrap();
    // A record with a list of at least two numbers, and its second number.
    let (id, value) = source
        .entities
        .iter()
        .find_map(|(&id, parts)| {
            parts[0].args.iter().find_map(|arg| match arg {
                Value::List(items)
                    if items.len() >= 2 && items.iter().all(|v| matches!(v, Value::Number(_))) =>
                {
                    Some((id, &items[1]))
                }
                _ => None,
            })
        })
        .unwrap();
    let cost = |source: &Source<'_>| {
        let mut b = Budget::new(1_000_000, None);
        let token = spelling(source, id, value, &mut b).unwrap().to_owned();
        (token, b.used())
    };
    // Off by default: every lookup costs the same.
    let (token, first) = cost(&source);
    assert_eq!(cost(&source), (token.clone(), first));
    // On: the first lookup builds the table, later ones are cheaper and
    // return the same original spelling.
    source.cache_spellings();
    let (cached, building) = cost(&source);
    assert_eq!(cached, token);
    assert!(building >= first);
    let (again, hit) = cost(&source);
    assert_eq!(again, token);
    assert!(hit < first, "{hit} {first}");
    // A leaf is never answered for a record that does not own it.
    let other = *source.entities.keys().find(|&&k| k != id).unwrap();
    let mut b = Budget::new(1_000_000, None);
    assert!(spelling(&source, other, value, &mut b).is_err());
    // Restarting the cache repeats the first run's work exactly.
    source.cache_spellings();
    assert_eq!(cost(&source).1, building);
}
