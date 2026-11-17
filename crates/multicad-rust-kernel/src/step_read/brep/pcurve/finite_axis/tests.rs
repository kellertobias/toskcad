use super::*;
const TOL: f64 = 1e-7;
const INPUT: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/src/step_read/tests/fixtures/occt-annulus.step"
));
fn fixture(end: f64, weighted: bool, alias: bool) -> String {
    let fixed = if alias {
        "-1.2500000000000000001"
    } else {
        "-1.25"
    };
    let points = format!(
        "#900000 = CARTESIAN_POINT('',(0.,-1.25));\n#900001 = CARTESIAN_POINT('',({end:.17},{fixed}));\n"
    );
    let curve = if weighted {
        "#60 = ( BOUNDED_CURVE() B_SPLINE_CURVE(1,(#900000,#900001),.UNSPECIFIED.,.F.,.F.) B_SPLINE_CURVE_WITH_KNOTS((2,2),(0.,1.),.UNSPECIFIED.) CURVE() GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_CURVE((2.,3.)) REPRESENTATION_ITEM('') );"
    } else {
        "#60 = B_SPLINE_CURVE_WITH_KNOTS('',1,(#900000,#900001),.UNSPECIFIED.,.F.,.F.,(2,2),(0.,1.),.UNSPECIFIED.);"
    };
    let text = INPUT.replace("#60 = LINE('',#61,#62);", curve);
    assert_ne!(text, INPUT);
    text.replacen(
        "ENDSEC;\nEND-ISO-10303-21;",
        &format!("{points}ENDSEC;\nEND-ISO-10303-21;"),
        1,
    )
}
fn pair(source: &Source<'_>) -> (uv::Definition, Surface, Edge) {
    let d = uv::definition(source, 58, &mut Budget::new(512_000_000, None)).unwrap();
    let curve = source.curve(53, 0).unwrap();
    let p = source.point(25).unwrap();
    let e = Edge::new(curve, [24, 24], [p, p], TOL).unwrap();
    let s = source.surface(d.surface).unwrap();
    (d, s, e)
}
fn run(source: &Source<'_>, budget: &mut Budget<'_>) -> Result<Option<Proof>> {
    let (d, s, e) = pair(source);
    certify(source, 53, &d, &s, &e, TOL, budget)
}
#[test]
fn original_finite_rational_cylinder_span_covers_actual_master_with_resource_retry() {
    for weighted in [false, true] {
        let text = fixture(-std::f64::consts::TAU, weighted, false);
        let source = Source::parse(&text, TOL).unwrap();
        let mut work = Budget::new(512_000_000, None);
        assert!(matches!(run(&source,&mut work).unwrap(),Some(Proof::Certified(e))if e<=TOL));
        let used = work.used();
        assert!(used > 0);
        assert!(matches!(
            run(&source, &mut Budget::new(used - 1, None)),
            Err(ReadError::Limit)
        ));
        assert!(matches!(
            run(&source, &mut Budget::new(used, None)).unwrap(),
            Some(Proof::Certified(_))
        ));
        let calls = std::cell::Cell::new(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() > 20
        };
        assert!(matches!(
            run(&source, &mut Budget::new(512_000_000, Some(&cancel))),
            Err(ReadError::Limit)
        ));
        assert!(calls.get() > 20);
        let all = std::cell::Cell::new(0);
        let observe = || {
            all.set(all.get() + 1);
            false
        };
        assert!(matches!(
            run(&source, &mut Budget::new(512_000_000, Some(&observe))).unwrap(),
            Some(Proof::Certified(_))
        ));
        let last = all.get();
        let visits = std::cell::Cell::new(0);
        let late = || {
            visits.set(visits.get() + 1);
            visits.get() == last
        };
        assert!(matches!(
            run(&source, &mut Budget::new(512_000_000, Some(&late))),
            Err(ReadError::Limit)
        ));
        assert_eq!(visits.get(), last);
        assert!(matches!(
            run(&source, &mut Budget::new(512_000_000, None)).unwrap(),
            Some(Proof::Certified(_))
        ));
    }
}
#[test]
fn finite_extent_direction_and_original_fixed_axis_cannot_be_erased() {
    for (end, alias) in [
        (-1., false),
        (std::f64::consts::TAU, false),
        (-std::f64::consts::TAU, true),
    ] {
        let text = fixture(end, true, alias);
        let source = Source::parse(&text, TOL).unwrap();
        assert!(matches!(
            run(&source, &mut Budget::new(512_000_000, None)),
            Ok(Some(Proof::Unsupported)) | Err(ReadError::Unsupported)
        ));
    }
    let text = fixture(-std::f64::consts::TAU, true, false);
    let source = Source::parse(&text, TOL).unwrap();
    assert!(matches!(
        run(&source, &mut Budget::new(512_000_000, None)).unwrap(),
        Some(Proof::Certified(_))
    ));
}
#[test]
fn actual_occt_analytic_families_keep_finite_uv_domains_after_literal_spline_encoding() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let mut transformed = 0;
    for name in ["cylinder", "cone-frustum", "sphere", "torus"] {
        let input = std::fs::read_to_string(root.join(format!(
            "../multicad-app/tests/fixtures/exchange/corpus/occt-{name}.step"
        )))
        .unwrap();
        let source = Source::parse(&input, TOL).unwrap();
        let mut output = input.clone();
        let mut records = String::new();
        let mut done = std::collections::BTreeSet::new();
        for (&id, parts) in &source.entities {
            if !matches!(parts.as_slice(),[e]if e.name=="PCURVE") {
                continue;
            }
            let d = uv::definition(&source, id, &mut Budget::new(512_000_000, None)).unwrap();
            let surface = source.surface(d.surface).unwrap();
            if !matches!(
                surface,
                Surface::Cylinder { .. }
                    | Surface::Cone { .. }
                    | Surface::Sphere { .. }
                    | Surface::Torus { .. }
            ) {
                continue;
            }
            let Curve::Line { origin, direction } = d.geometry else {
                continue;
            };
            if !done.insert(d.curve_id) {
                continue;
            }
            let first = 900000 + 2 * done.len();
            let last = first + 1;
            for (point_id, t) in [(first, -256.), (last, 256.)] {
                let p: [f64; 2] = std::array::from_fn(|i| origin[i] + t * direction[i]);
                records.push_str(&format!(
                    "#{point_id} = CARTESIAN_POINT('',({:.17},{:.17}));\n",
                    p[0], p[1]
                ));
            }
            let old = format!(
                "#{} = {};",
                d.curve_id,
                source.simple(d.curve_id).unwrap().source.unwrap()
            );
            let new = format!(
                "#{} = B_SPLINE_CURVE_WITH_KNOTS('',1,(#{first},#{last}),.UNSPECIFIED.,.F.,.F.,(2,2),(0.,1.),.UNSPECIFIED.);",
                d.curve_id
            );
            assert!(output.contains(&old));
            output = output.replacen(&old, &new, 1);
        }
        if name == "sphere" {
            // This tracked fixture contains no supplementary curves; it is
            // not evidence for the finite SphereV branch.
            assert!(done.is_empty());
            assert!(
                source
                    .entities
                    .values()
                    .flatten()
                    .all(|e| e.name != "PCURVE")
            );
        } else {
            assert!(
                !done.is_empty(),
                "{name} must exercise actual finite traces"
            );
        }
        transformed += done.len();
        output = output.replacen(
            "ENDSEC;\nEND-ISO-10303-21;",
            &format!("{records}ENDSEC;\nEND-ISO-10303-21;"),
            1,
        );
        let report = crate::step_read::brep::validate_pcurves(&output, TOL)
            .unwrap_or_else(|e| panic!("{name}: {e:?}"));
        if name == "sphere" {
            assert_eq!(output, input);
            assert_eq!(report.pcurves_checked, 0);
        } else {
            assert!(report.pcurves_checked >= done.len());
        }
        assert!(report.certified_maximum_deviation_mm <= TOL);
    }
    assert!(transformed >= 12);
}
fn cardinal_target(vertex: &str) -> String {
    // This is a numeric finite-endpoint contract on actual parsed operators,
    // not an assertion that the independently modified shell is a new solid.
    INPUT
        .replace("#53 = CIRCLE('',#54,6.125);", "#53 = CIRCLE('',#54,5.);")
        .replace(
            "#55 = CARTESIAN_POINT('',(5.452019289956,-10.165409182,21.));",
            "#55 = CARTESIAN_POINT('',(0.,0.,0.));",
        )
        .replace(
            "#56 = DIRECTION('',(0.361615431965,-0.932327345606,0.));",
            "#56 = DIRECTION('',(0.,0.,1.));",
        )
        .replace(
            "#57 = DIRECTION('',(0.932327345606,0.361615431965,0.));",
            "#57 = DIRECTION('',(1.,0.,0.));",
        )
        .replace(
            "#25 = CARTESIAN_POINT('',(11.162524281793,-7.950514661222,21.));",
            &format!("#25 = CARTESIAN_POINT('',({vertex}));"),
        )
}
#[test]
fn original_vertex_real_shadow_is_bounded_or_refused_not_just_cached_inversion() {
    let original = cardinal_target("3.,4.,0.");
    let changed = cardinal_target("3.0000000000000001,4.,0.");
    assert_ne!(original, changed);
    let a = Source::parse(&original, TOL).unwrap();
    let b = Source::parse(&changed, TOL).unwrap();
    assert_eq!(a.point(25).unwrap(), b.point(25).unwrap());
    for source in [&a, &b] {
        let (_, _, edge) = pair(source);
        let reserve =
            endpoint_reserve(source, 53, &edge, &mut Budget::new(512_000_000, None)).unwrap();
        assert!(reserve > 1e-17 && reserve < TOL);
    }
    // A nonzero original REAL hidden behind cached zero has no qualified
    // adjacent finite enclosure. It must not reuse the ordinary endpoint.
    let zero = cardinal_target("5.,0.,0.");
    let shadow = cardinal_target("5.,1E-1000,0.");
    let zero = Source::parse(&zero, TOL).unwrap();
    let shadow = Source::parse(&shadow, TOL).unwrap();
    assert_eq!(zero.point(25).unwrap(), shadow.point(25).unwrap());
    let (_, _, ordinary) = pair(&zero);
    let (_, _, hidden) = pair(&shadow);
    assert!(endpoint_reserve(&zero, 53, &ordinary, &mut Budget::new(512_000_000, None)).is_ok());
    assert!(matches!(
        endpoint_reserve(&shadow, 53, &hidden, &mut Budget::new(512_000_000, None)),
        Err(ReadError::Unsupported)
    ));
}
#[test]
fn true_period_enclosure_has_independent_exact_machin_bounds() {
    use num_bigint::BigInt;
    use num_rational::BigRational as Q;
    fn sum(d: u32, n: usize) -> Q {
        (0..n).fold(Q::from_integer(0.into()), |s, k| {
            let t = Q::new(
                1.into(),
                BigInt::from(2 * k + 1) * BigInt::from(d).pow((2 * k + 1) as u32),
            );
            if k % 2 == 0 { s + t } else { s - t }
        })
    }
    let lower = Q::from_integer(2.into())
        * (Q::from_integer(16.into()) * sum(5, 20) - Q::from_integer(4.into()) * sum(239, 21));
    let upper = Q::from_integer(2.into())
        * (Q::from_integer(16.into()) * sum(5, 21) - Q::from_integer(4.into()) * sum(239, 20));
    assert!(Q::from_float(std::f64::consts::TAU.next_down()).unwrap() < lower);
    assert!(lower < upper);
    assert!(upper < Q::from_float(std::f64::consts::TAU.next_up()).unwrap());
}
