use super::*;
const BASE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../multicad-app/tests/fixtures/exchange/corpus/occt-cone-frustum.step"
));
const XY: [(&str, &str); 8] = [
    ("10.", "20."),
    ("10.5", "21."),
    ("11.", "23."),
    ("11.5", "20.25"),
    ("12.", "22.5"),
    ("12.5", "23."),
    ("13.", "21."),
    ("14.", "24."),
];
const WEIGHTS: [i64; 8] = [1, 2, 3, 4, 4, 3, 2, 1];
fn input(rational_weights: bool) -> String {
    let points = XY
        .iter()
        .enumerate()
        .map(|(i, (x, y))| format!("#{}=CARTESIAN_POINT('',({x},{y}));", 900010 + i))
        .collect::<Vec<_>>()
        .join("\n");
    let curve = if rational_weights {
        "#900001=(BOUNDED_CURVE() B_SPLINE_CURVE(4,(#900010,#900011,#900012,#900013,#900014,#900015,#900016,#900017),.UNSPECIFIED.,.F.,.F.) B_SPLINE_CURVE_WITH_KNOTS((5,3,5),(0.,0.5,1.),.UNSPECIFIED.) CURVE() GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_CURVE((1.,2.,3.,4.,4.,3.,2.,1.)) REPRESENTATION_ITEM(''));"
    } else {
        "#900001=B_SPLINE_CURVE_WITH_KNOTS('',4,(#900010,#900011,#900012,#900013,#900014,#900015,#900016,#900017),.UNSPECIFIED.,.F.,.F.,(5,3,5),(0.,0.5,1.),.UNSPECIFIED.);"
    };
    let added = format!(
        "\n{curve}\n{points}\n#900100=PCURVE('',#31,#900101);\n#900101=DEFINITIONAL_REPRESENTATION('',(#900001),#900102);\n#900102=(GEOMETRIC_REPRESENTATION_CONTEXT(2) PARAMETRIC_REPRESENTATION_CONTEXT() REPRESENTATION_CONTEXT('',''));\n"
    );
    let at = BASE.rfind("ENDSEC;").unwrap();
    format!("{}{}{}", &BASE[..at], added, &BASE[at..])
}
fn definition(source: &Source<'_>) -> super::super::super::uv::Definition {
    super::super::super::uv::definition(source, 900100, &mut Budget::new(512_000_000, None))
        .unwrap()
}
fn decimal(text: &str) -> BigRational {
    let (whole, fraction) = text
        .trim_end_matches('.')
        .split_once('.')
        .unwrap_or((text.trim_end_matches('.'), ""));
    let denominator = BigInt::from(10).pow(fraction.len() as u32);
    let numerator = whole.parse::<BigInt>().unwrap() * &denominator
        + if fraction.is_empty() {
            BigInt::from(0)
        } else {
            fraction.parse::<BigInt>().unwrap()
        };
    BigRational::new(numerator, denominator)
}
/// Independent exact Cox-de Boor recursion, unrelated to knot insertion.
fn basis(i: usize, p: usize, t: &BigRational, u: &[BigRational]) -> BigRational {
    if p == 0 {
        return BigRational::from_integer(if u[i] <= *t && *t < u[i + 1] {
            1.into()
        } else {
            0.into()
        });
    }
    let zero = BigRational::from_integer(0.into());
    let left = if u[i + p] > u[i] {
        (t - &u[i]) / (&u[i + p] - &u[i]) * basis(i, p - 1, t, u)
    } else {
        zero.clone()
    };
    let right = if u[i + p + 1] > u[i + 1] {
        (&u[i + p + 1] - t) / (&u[i + p + 1] - &u[i + 1]) * basis(i + 1, p - 1, t, u)
    } else {
        zero
    };
    left + right
}
fn original(t: &BigRational, rational_weights: bool, normalized: bool) -> [BigRational; 2] {
    let u = [
        vec![decimal("0"); 5],
        vec![decimal("0.5"); 3],
        vec![decimal("1"); 5],
    ]
    .concat();
    let mut h: [BigRational; 3] = std::array::from_fn(|_| decimal("0"));
    for i in 0..8 {
        let n = if *t == decimal("1") {
            BigRational::from_integer(if i == 7 { 1.into() } else { 0.into() })
        } else {
            basis(i, 4, t, &u)
        };
        let w = BigRational::from_integer(if rational_weights { WEIGHTS[i] } else { 1 }.into());
        let mut x = decimal(XY[i].0);
        let mut y = decimal(XY[i].1);
        if normalized {
            x = (x - decimal("10")) / decimal("4");
            y = (y - decimal("20")) / decimal("4");
        }
        h[0] += &n * &w * x;
        h[1] += &n * &w * y;
        h[2] += n * w;
    }
    [&h[0] / &h[2], &h[1] / &h[2]]
}
fn evaluate(cell: &Cell, s: f64) -> V {
    let mut h = cell.controls().to_vec();
    for n in (1..h.len()).rev() {
        for i in 0..n {
            h[i] = std::array::from_fn(|j| {
                (I::point(1.) - I::point(s)) * h[i][j] + I::point(s) * h[i + 1][j]
            });
        }
    }
    [h[0][0] / h[0][3], h[0][1] / h[0][3], h[0][2] / h[0][3]]
}
fn assert_exact_contains(i: I, q: &BigRational) {
    assert!(i.valid());
    assert!(BigRational::from_float(i.lo).unwrap() <= *q);
    assert!(*q <= BigRational::from_float(i.hi).unwrap());
}
#[test]
fn polynomial_and_rational_quartics_with_internal_multiplicity_three_match_exact_cox_de_boor() {
    for rational_weights in [false, true] {
        let text = input(rational_weights);
        let source = Source::parse(&text, 1e-7).unwrap();
        let definition = definition(&source);
        let Curve::Nurbs { curve, .. } = &definition.geometry else {
            panic!()
        };
        let min = [decimal("10"), decimal("20")];
        let max = [decimal("14"), decimal("24")];
        let chart = if rational_weights {
            Coordinates::UnitBox2 {
                min: &min,
                max: &max,
            }
        } else {
            Coordinates::Raw(2)
        };
        let result = original_curve(
            &source,
            900001,
            curve,
            chart,
            &mut Budget::new(512_000_000, None),
        )
        .unwrap();
        assert!(result.matches(&source, 900001));
        assert!(!result.matches(&Source::parse(&text, 1e-7).unwrap(), 900001));
        assert!(!result.matches(&source, 900002));
        assert_eq!(result.degree(), 4);
        assert_eq!(result.cells().len(), 2);
        assert!(result.retained_bytes() <= result.peak_bytes());
        assert!(result.peak_bytes() <= MAX_BYTES);
        assert_eq!(
            matches!(result.coordinates(), Coordinates::UnitBox2 { .. }),
            rational_weights
        );
        for cell in result.cells() {
            let [a, b] = cell.knot_indices();
            let knots = [decimal("0"), decimal("0.5"), decimal("1")];
            assert_exact_contains(cell.domain()[0], &knots[a]);
            assert_exact_contains(cell.domain()[1], &knots[b]);
            for s in [0., 0.125, 0.5, 0.875, 1.] {
                let q = BigRational::from_float(s).unwrap();
                let t = &knots[a] + (&knots[b] - &knots[a]) * q;
                let expected = original(&t, rational_weights, rational_weights);
                let actual = evaluate(cell, s);
                for i in 0..2 {
                    assert_exact_contains(actual[i], &expected[i]);
                }
                assert_exact_contains(actual[2], &decimal("0"));
            }
        }
    }
}
#[test]
fn original_decimal_integer_alias_weight_and_knot_failures_do_not_trust_decoded_structure() {
    let text = input(true);
    let source = Source::parse(&text, 1e-7).unwrap();
    let definition = definition(&source);
    let Curve::Nurbs { curve, .. } = &definition.geometry else {
        panic!()
    };
    for (a, b) in [
        (
            "B_SPLINE_CURVE(4,",
            "B_SPLINE_CURVE(4.0000000000000000000001,",
        ),
        (
            ".UNSPECIFIED.,.F.,.F.) B_SPLINE_CURVE_WITH_KNOTS",
            ".UNSPECIFIED.,.T.,.F.) B_SPLINE_CURVE_WITH_KNOTS",
        ),
        (
            "RATIONAL_B_SPLINE_CURVE((1.,2.",
            "RATIONAL_B_SPLINE_CURVE((1E-1000,2.",
        ),
        ("(0.,0.5,1.),", "(0.,0.5,1E-1000),"),
    ] {
        let bad = text.replace(a, b);
        assert_ne!(bad, text);
        let actual = Source::parse(&bad, 1e-7).unwrap();
        // A real decoded object is only a structural proposal. The new actual
        // source must independently validate EVERY original numeric premise.
        assert!(
            original_curve(
                &actual,
                900001,
                curve,
                Coordinates::Raw(2),
                &mut Budget::new(512_000_000, None)
            )
            .is_err()
        );
    }
    let min = [decimal("10"), decimal("20")];
    let max = [decimal("13"), decimal("24")];
    assert!(
        original_curve(
            &source,
            900001,
            curve,
            Coordinates::UnitBox2 {
                min: &min,
                max: &max
            },
            &mut Budget::new(512_000_000, None)
        )
        .is_err()
    );
}
#[test]
fn extraction_exact_work_frontier_late_cancel_and_retry_keep_actual_source() {
    let text = input(true);
    let source = Source::parse(&text, 1e-7).unwrap();
    let definition = definition(&source);
    let Curve::Nurbs { curve, .. } = &definition.geometry else {
        panic!()
    };
    let run = |b: &mut Budget<'_>| original_curve(&source, 900001, curve, Coordinates::Raw(2), b);
    let mut full = Budget::new(512_000_000, None);
    let first = run(&mut full).unwrap();
    let used = full.used();
    let mut exact = Budget::new(used, None);
    let second = run(&mut exact).unwrap();
    assert_eq!(exact.remaining, 0);
    assert_eq!(first.peak_bytes(), second.peak_bytes());
    assert!(matches!(
        run(&mut Budget::new(used - 1, None)),
        Err(ReadError::Limit)
    ));
    let calls = std::cell::Cell::new(0);
    let stop = || {
        calls.set(calls.get() + 1);
        calls.get() == 30
    };
    let mut cancelled = Budget::new(512_000_000, Some(&stop));
    assert!(run(&mut cancelled).is_err());
    assert!(cancelled.cancelled.get());
    let mut retry = Budget::new(512_000_000, None);
    let after = run(&mut retry).unwrap();
    assert_eq!(retry.used(), used);
    assert_eq!(first.cells().len(), after.cells().len());
    for (a, b) in first.cells().iter().zip(after.cells()) {
        for (a, b) in a.controls().iter().zip(b.controls()) {
            for i in 0..4 {
                assert_eq!((a[i].lo, a[i].hi), (b[i].lo, b[i].hi));
            }
        }
    }
}

#[test]
fn actual_quartic_101_pole_source_extracts_every_original_span_with_shared_end_coefficients() {
    let references = (0..101)
        .map(|i| format!("#{}", 910000 + i))
        .collect::<Vec<_>>()
        .join(",");
    let points = (0..101)
        .map(|i| format!("#{}=CARTESIAN_POINT('',({i}.,{}.));", 910000 + i, i % 7))
        .collect::<Vec<_>>()
        .join("\n");
    let weights = (0..101)
        .map(|i| format!("{}.", 1 + i % 5))
        .collect::<Vec<_>>()
        .join(",");
    let mult = (0..34)
        .map(|i| if i == 0 || i == 33 { "5" } else { "3" })
        .collect::<Vec<_>>()
        .join(",");
    let knots = (0..34)
        .map(|i| format!("{i}."))
        .collect::<Vec<_>>()
        .join(",");
    let records = format!(
        "\n#900001=(BOUNDED_CURVE() B_SPLINE_CURVE(4,({references}),.UNSPECIFIED.,.F.,.F.) B_SPLINE_CURVE_WITH_KNOTS(({mult}),({knots}),.UNSPECIFIED.) CURVE() GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_CURVE(({weights})) REPRESENTATION_ITEM(''));\n{points}\n#900100=PCURVE('',#31,#900101);\n#900101=DEFINITIONAL_REPRESENTATION('',(#900001),#900102);\n#900102=(GEOMETRIC_REPRESENTATION_CONTEXT(2) PARAMETRIC_REPRESENTATION_CONTEXT() REPRESENTATION_CONTEXT('',''));\n"
    );
    let at = BASE.rfind("ENDSEC;").unwrap();
    let text = format!("{}{}{}", &BASE[..at], records, &BASE[at..]);
    let source = Source::parse(&text, 1e-7).unwrap();
    let definition = definition(&source);
    let Curve::Nurbs { curve, .. } = &definition.geometry else {
        panic!()
    };
    let min = [decimal("0"), decimal("0")];
    let max = [decimal("100"), decimal("6")];
    let mut work = Budget::new(512_000_000, None);
    let cells = original_curve(
        &source,
        900001,
        curve,
        Coordinates::UnitBox2 {
            min: &min,
            max: &max,
        },
        &mut work,
    )
    .unwrap();
    assert_eq!(cells.cells().len(), 33);
    assert!(work.used() < 512_000_000);
    assert!(cells.peak_bytes() < MAX_BYTES);
    for (i, cell) in cells.cells().iter().enumerate() {
        assert_eq!(cell.knot_indices(), [i, i + 1]);
    }
    for pair in cells.cells().windows(2) {
        for i in 0..4 {
            let (a, b) = (
                pair[0].controls().last().unwrap()[i],
                pair[1].controls()[0][i],
            );
            assert_eq!((a.lo, a.hi), (b.lo, b.hi));
        }
    }
    let mut original_knots = vec![decimal("0"); 5];
    for i in 1..33 {
        original_knots.extend(std::iter::repeat_n(BigRational::from_integer(i.into()), 3));
    }
    original_knots.extend(std::iter::repeat_n(decimal("33"), 5));
    for index in [0, 16, 32] {
        for s in [0., 0.25, 0.5, 0.75, 1.] {
            let t = BigRational::from_integer(index.into()) + BigRational::from_float(s).unwrap();
            let mut h: [BigRational; 3] = std::array::from_fn(|_| decimal("0"));
            for i in 0..101 {
                let n = if t == decimal("33") {
                    BigRational::from_integer(if i == 100 { 1.into() } else { 0.into() })
                } else {
                    basis(i, 4, &t, &original_knots)
                };
                let w = BigRational::from_integer((1 + i % 5).into());
                h[0] += &n * &w * BigRational::new(i.into(), 100.into());
                h[1] += &n * &w * BigRational::new((i % 7).into(), 6.into());
                h[2] += n * w;
            }
            let actual = evaluate(&cells.cells()[index as usize], s);
            assert_exact_contains(actual[0], &(&h[0] / &h[2]));
            assert_exact_contains(actual[1], &(&h[1] / &h[2]));
        }
    }
}

fn basis_table(p: usize, t: &BigRational, knots: &[BigRational]) -> Vec<BigRational> {
    let zero = decimal("0");
    let mut current = knots
        .windows(2)
        .map(|w| {
            BigRational::from_integer(if w[0] <= *t && *t < w[1] {
                1.into()
            } else {
                0.into()
            })
        })
        .collect::<Vec<_>>();
    for level in 1..=p {
        let mut next = Vec::with_capacity(current.len() - 1);
        for i in 0..current.len() - 1 {
            let a = if current[i].is_zero() || knots[i + level] == knots[i] {
                zero.clone()
            } else {
                (t - &knots[i]) / (&knots[i + level] - &knots[i]) * &current[i]
            };
            let b = if current[i + 1].is_zero() || knots[i + level + 1] == knots[i + 1] {
                zero.clone()
            } else {
                (&knots[i + level + 1] - t) / (&knots[i + level + 1] - &knots[i + 1])
                    * &current[i + 1]
            };
            next.push(a + b);
        }
        current = next;
    }
    current
}
#[test]
fn actual_dense_cubic_and_degree_eleven_sources_keep_raw_three_space_and_capacity_bounds() {
    for (degree, span_count, multiplicity, unit_weights) in [
        (11_usize, 28_usize, 11_usize, false),
        (11, 49, 11, false),
        (11, 49, 10, false),
        (11, 49, 10, true),
        (3, 256, 3, true),
        (11, 256, 1, true),
        (3, 257, 3, true),
    ] {
        let n = degree + 1 + (span_count - 1) * multiplicity;
        let references = (0..n)
            .map(|i| format!("#{}", 910000 + i))
            .collect::<Vec<_>>()
            .join(",");
        let points = (0..n)
            .map(|i| {
                format!(
                    "#{}=CARTESIAN_POINT('',({i}.,{}.,-{}.{}));",
                    910000 + i,
                    i % 7,
                    i / 2,
                    if i % 2 == 0 { 0 } else { 5 }
                )
            })
            .collect::<Vec<_>>()
            .join("\n");
        let weights = (0..n)
            .map(|i| format!("{}.", if unit_weights { 1 } else { 1 + i % 5 }))
            .collect::<Vec<_>>()
            .join(",");
        let mult = (0..=span_count)
            .map(|i| {
                if i == 0 || i == span_count {
                    degree + 1
                } else {
                    multiplicity
                }
                .to_string()
            })
            .collect::<Vec<_>>()
            .join(",");
        let knots = (0..=span_count)
            .map(|i| format!("{i}."))
            .collect::<Vec<_>>()
            .join(",");
        let records = format!(
            "\n#900001=(BOUNDED_CURVE() B_SPLINE_CURVE({degree},({references}),.UNSPECIFIED.,.F.,.F.) B_SPLINE_CURVE_WITH_KNOTS(({mult}),({knots}),.UNSPECIFIED.) CURVE() GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_CURVE(({weights})) REPRESENTATION_ITEM(''));\n{points}\n"
        );
        let at = BASE.rfind("ENDSEC;").unwrap();
        let text = format!("{}{}{}", &BASE[..at], records, &BASE[at..]);
        let source = Source::parse(&text, 1e-7).unwrap();
        let Curve::Nurbs { curve, .. } = source.curve(900001, 0).unwrap() else {
            panic!()
        };
        assert_eq!(curve.points().len(), n);
        let mut budget = Budget::new(512_000_000, None);
        let extracted = original_curve(&source, 900001, &curve, Coordinates::Raw(3), &mut budget);
        if span_count == 257 {
            assert!(matches!(extracted, Err(ReadError::Unsupported)));
            continue;
        }
        let result = extracted.unwrap();
        assert_eq!(result.degree(), degree);
        assert_eq!(result.cells().len(), span_count);
        if unit_weights {
            assert!(result.cells().iter().all(|cell| {
                cell.controls()
                    .iter()
                    .all(|p| p[3].lo == 1. && p[3].hi == 1.)
            }));
            let almost = text.replacen(
                "RATIONAL_B_SPLINE_CURVE((1.,",
                "RATIONAL_B_SPLINE_CURVE((1.0000000000000000000001,",
                1,
            );
            assert_ne!(almost, text);
            let almost_source = Source::parse(&almost, 1e-7).unwrap();
            let Curve::Nurbs {
                curve: almost_curve,
                ..
            } = almost_source.curve(900001, 0).unwrap()
            else {
                panic!()
            };
            let actual = original_curve(
                &almost_source,
                900001,
                &almost_curve,
                Coordinates::Raw(3),
                &mut Budget::new(512_000_000, None),
            )
            .unwrap();
            assert!(actual.cells().iter().any(|cell| {
                cell.controls()
                    .iter()
                    .any(|p| p[3].lo != 1. || p[3].hi != 1.)
            }));
        }

        assert!(result.peak_bytes() < MAX_BYTES);
        assert!(budget.used() < 512_000_000);
        if span_count == 256 {
            assert!(result.peak_bytes() < 1024 * 1024);
            assert_eq!(result.cells().first().unwrap().knot_indices(), [0, 1]);
            assert_eq!(result.cells().last().unwrap().knot_indices(), [255, 256]);
            let used = budget.used();
            let run = |b: &mut Budget<'_>| {
                original_curve(&source, 900001, &curve, Coordinates::Raw(3), b)
            };
            assert!(matches!(
                run(&mut Budget::new(used - 1, None)),
                Err(ReadError::Limit)
            ));
            assert!(run(&mut Budget::new(used, None)).is_ok());
            let calls = std::cell::Cell::new(0usize);
            let observe = || {
                calls.set(calls.get() + 1);
                false
            };
            assert!(run(&mut Budget::new(used, Some(&observe))).is_ok());
            let total = calls.get();
            calls.set(0);
            let cancel = || {
                calls.set(calls.get() + 1);
                calls.get() == total
            };
            let mut interrupted = Budget::new(used, Some(&cancel));
            assert!(matches!(run(&mut interrupted), Err(ReadError::Limit)));
            assert!(interrupted.cancelled.get());
            assert!(run(&mut Budget::new(used, None)).is_ok());
        }
        let mut u = vec![decimal("0"); degree + 1];
        for i in 1..span_count {
            u.extend(std::iter::repeat_n(
                BigRational::from_integer(i.into()),
                multiplicity,
            ));
        }
        u.extend(std::iter::repeat_n(
            BigRational::from_integer(span_count.into()),
            degree + 1,
        ));
        for index in [0, span_count / 2, span_count - 1] {
            let t = BigRational::from_integer(index.into()) + decimal("0.5");
            let basis = basis_table(degree, &t, &u);
            assert_eq!(basis.len(), n);
            let mut h: [BigRational; 4] = std::array::from_fn(|_| decimal("0"));
            for (i, coefficient) in basis.iter().enumerate() {
                let w =
                    BigRational::from_integer((if unit_weights { 1 } else { 1 + i % 5 }).into());
                let nw = coefficient * w;
                h[0] += &nw * BigRational::from_integer(i.into());
                h[1] += &nw * BigRational::from_integer((i % 7).into());
                h[2] += &nw * BigRational::new(-(BigInt::from(i)), 2.into());
                h[3] += nw;
            }
            let actual = evaluate(&result.cells()[index], 0.5);
            for j in 0..3 {
                assert_exact_contains(actual[j], &(&h[j] / &h[3]));
            }
        }
    }
}

#[test]
fn raw_polynomial_and_rational_quartics_still_match_exact_cox_de_boor() {
    for rational_weights in [false, true] {
        let text = input(rational_weights);
        let source = Source::parse(&text, 1e-7).unwrap();
        let definition = definition(&source);
        let Curve::Nurbs { curve, .. } = &definition.geometry else {
            panic!()
        };
        let chart = Coordinates::Raw(2);
        let result = original_curve(
            &source,
            900001,
            curve,
            chart,
            &mut Budget::new(512_000_000, None),
        )
        .unwrap();
        assert!(result.matches(&source, 900001));
        assert!(!result.matches(&Source::parse(&text, 1e-7).unwrap(), 900001));
        assert!(!result.matches(&source, 900002));
        assert_eq!(result.degree(), 4);
        assert_eq!(result.cells().len(), 2);
        assert!(result.retained_bytes() <= result.peak_bytes());
        assert!(result.peak_bytes() <= MAX_BYTES);
        assert!(matches!(result.coordinates(), Coordinates::Raw(2)));
        for cell in result.cells() {
            let [a, b] = cell.knot_indices();
            let knots = [decimal("0"), decimal("0.5"), decimal("1")];
            assert_exact_contains(cell.domain()[0], &knots[a]);
            assert_exact_contains(cell.domain()[1], &knots[b]);
            for s in [0., 0.125, 0.5, 0.875, 1.] {
                let q = BigRational::from_float(s).unwrap();
                let t = &knots[a] + (&knots[b] - &knots[a]) * q;
                let expected = original(&t, rational_weights, false);
                let actual = evaluate(cell, s);
                for i in 0..2 {
                    assert_exact_contains(actual[i], &expected[i]);
                }
                assert_exact_contains(actual[2], &decimal("0"));
            }
        }
    }
}

#[test]
fn selected_complex_record_numeric_scratch_includes_ignored_and_typed_values() {
    let ordinary = input(true);
    for count in [128, 2049] {
        let numbers = std::iter::repeat_n("1.", count)
            .collect::<Vec<_>>()
            .join(",");
        // These are actual parsed tokens, not forged Source/Value fields.
        // Ignored component arguments still belong to spelling's traversal.
        for text in [
            ordinary.replace(
                "BOUNDED_CURVE()",
                &format!("BOUNDED_CURVE(NESTED(({numbers})))"),
            ),
            ordinary.replacen(
                "),.UNSPECIFIED.,.F.,.F.) B_SPLINE_CURVE_WITH_KNOTS",
                &format!("),NESTED(({numbers})),.F.,.F.) B_SPLINE_CURVE_WITH_KNOTS"),
                1,
            ),
        ] {
            assert_ne!(text, ordinary);
            let source = Source::parse(&text, 1e-7).unwrap();
            let definition = definition(&source);
            let Curve::Nurbs { curve, .. } = &definition.geometry else {
                panic!()
            };
            let run = |budget: &mut Budget<'_>| {
                original_curve(&source, 900001, curve, Coordinates::Raw(2), budget)
            };
            let mut full = Budget::new(512_000_000, None);
            if count == 2049 {
                assert!(matches!(run(&mut full), Err(ReadError::Limit)));
            } else {
                assert!(run(&mut full).is_ok());
                let used = full.used();
                assert!(matches!(
                    run(&mut Budget::new(used - 1, None)),
                    Err(ReadError::Limit)
                ));
                assert!(run(&mut Budget::new(used, None)).is_ok());
            }
            let cancel = || true;
            let mut interrupted = Budget::new(512_000_000, Some(&cancel));
            assert!(matches!(run(&mut interrupted), Err(ReadError::Limit)));
            assert!(interrupted.cancelled.get());
        }
    }
}

#[test]
fn raw_quartic_101_pole_source_retains_every_original_span() {
    let references = (0..101)
        .map(|i| format!("#{}", 910000 + i))
        .collect::<Vec<_>>()
        .join(",");
    let points = (0..101)
        .map(|i| format!("#{}=CARTESIAN_POINT('',({i}.,{}.));", 910000 + i, i % 7))
        .collect::<Vec<_>>()
        .join("\n");
    let weights = (0..101)
        .map(|i| format!("{}.", 1 + i % 5))
        .collect::<Vec<_>>()
        .join(",");
    let mult = (0..34)
        .map(|i| if i == 0 || i == 33 { "5" } else { "3" })
        .collect::<Vec<_>>()
        .join(",");
    let knots = (0..34)
        .map(|i| format!("{i}."))
        .collect::<Vec<_>>()
        .join(",");
    let records = format!(
        "\n#900001=(BOUNDED_CURVE() B_SPLINE_CURVE(4,({references}),.UNSPECIFIED.,.F.,.F.) B_SPLINE_CURVE_WITH_KNOTS(({mult}),({knots}),.UNSPECIFIED.) CURVE() GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_CURVE(({weights})) REPRESENTATION_ITEM(''));\n{points}\n#900100=PCURVE('',#31,#900101);\n#900101=DEFINITIONAL_REPRESENTATION('',(#900001),#900102);\n#900102=(GEOMETRIC_REPRESENTATION_CONTEXT(2) PARAMETRIC_REPRESENTATION_CONTEXT() REPRESENTATION_CONTEXT('',''));\n"
    );
    let at = BASE.rfind("ENDSEC;").unwrap();
    let text = format!("{}{}{}", &BASE[..at], records, &BASE[at..]);
    let source = Source::parse(&text, 1e-7).unwrap();
    let definition = definition(&source);
    let Curve::Nurbs { curve, .. } = &definition.geometry else {
        panic!()
    };
    let mut work = Budget::new(512_000_000, None);
    let cells = original_curve(&source, 900001, curve, Coordinates::Raw(2), &mut work).unwrap();
    assert_eq!(cells.cells().len(), 33);
    assert!(work.used() < 512_000_000);
    assert!(cells.peak_bytes() < MAX_BYTES);
    for (i, cell) in cells.cells().iter().enumerate() {
        assert_eq!(cell.knot_indices(), [i, i + 1]);
    }
    for pair in cells.cells().windows(2) {
        for i in 0..4 {
            let (a, b) = (
                pair[0].controls().last().unwrap()[i],
                pair[1].controls()[0][i],
            );
            assert_eq!((a.lo, a.hi), (b.lo, b.hi));
        }
    }
    let mut original_knots = vec![decimal("0"); 5];
    for i in 1..33 {
        original_knots.extend(std::iter::repeat_n(BigRational::from_integer(i.into()), 3));
    }
    original_knots.extend(std::iter::repeat_n(decimal("33"), 5));
    for index in [0, 16, 32] {
        for s in [0., 0.25, 0.5, 0.75, 1.] {
            let t = BigRational::from_integer(index.into()) + BigRational::from_float(s).unwrap();
            let mut h: [BigRational; 3] = std::array::from_fn(|_| decimal("0"));
            for i in 0..101 {
                let n = if t == decimal("33") {
                    BigRational::from_integer(if i == 100 { 1.into() } else { 0.into() })
                } else {
                    basis(i, 4, &t, &original_knots)
                };
                let w = BigRational::from_integer((1 + i % 5).into());
                h[0] += &n * &w * BigRational::from_integer(i.into());
                h[1] += &n * &w * BigRational::from_integer((i % 7).into());
                h[2] += n * w;
            }
            let actual = evaluate(&cells.cells()[index as usize], s);
            assert_exact_contains(actual[0], &(&h[0] / &h[2]));
            assert_exact_contains(actual[1], &(&h[1] / &h[2]));
        }
    }
}

#[test]
fn original_point_typed_label_is_bounded_for_raw_and_unitbox_extraction() {
    let ordinary = input(false);
    let oversized = std::iter::repeat_n("1.", 2049)
        .collect::<Vec<_>>()
        .join(",");
    let bad = ordinary.replace(
        "#900011=CARTESIAN_POINT('',",
        &format!("#900011=CARTESIAN_POINT(TYPED(({oversized})),"),
    );
    assert_ne!(bad, ordinary);
    let min = [decimal("10."), decimal("20.")];
    let max = [decimal("14."), decimal("24.")];
    for text in [&bad, &ordinary] {
        let source = Source::parse(text, 1e-7).unwrap();
        let definition = definition(&source);
        let Curve::Nurbs { curve, .. } = &definition.geometry else {
            panic!()
        };
        for chart in [
            Coordinates::Raw(2),
            Coordinates::UnitBox2 {
                min: &min,
                max: &max,
            },
        ] {
            let result = original_curve(
                &source,
                900001,
                curve,
                chart,
                &mut Budget::new(512_000_000, None),
            );
            if std::ptr::eq(text, &bad) {
                assert!(matches!(result, Err(ReadError::Limit)));
            } else {
                assert!(result.is_ok());
            }
        }
    }
}
