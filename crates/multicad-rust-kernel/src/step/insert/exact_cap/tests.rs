use super::*;
#[test]
fn exact_file_integer_and_scale_matching_never_uses_floating_aliases() {
    let mut w = Work::new(16_000_000, None).unwrap();
    let mut memory = Memory::new(0, SCRATCH).unwrap();
    for text in ["9007199254740993.", "+9007199254740993.E0"] {
        let value = BigInt::from(9_007_199_254_740_993_u64);
        assert_eq!(text.parse::<f64>().unwrap(), 9_007_199_254_740_992.);
        assert!(equals_dyadic(text, 1, &value.to_bytes_le().1, 0, &mut memory, &mut w).unwrap());
        assert!(
            !equals_dyadic(
                text,
                1,
                &BigInt::from(9_007_199_254_740_992_u64).to_bytes_le().1,
                0,
                &mut memory,
                &mut w
            )
            .unwrap()
        );
    }
    for (text, matched) in [
        ("0.125", true),
        ("1.25E-1", true),
        ("0.125000000000000001", false),
    ] {
        assert_eq!(text.parse::<f64>().unwrap(), 0.125);
        assert_eq!(
            equals_dyadic(text, 1, &[1], 3, &mut memory, &mut w).unwrap(),
            matched
        );
    }
    assert!(equals_dyadic("-0.", 0, &[], 0, &mut memory, &mut w).unwrap());
    assert!(equals_dyadic("-2.", -1, &[2], 0, &mut memory, &mut w).unwrap());
    assert_eq!(memory.used, 0);
}
#[test]
fn literal_sub_ulp_ruled_gaps_and_owned_corner_factors_keep_their_true_jacobian() {
    for factor in [false, true] {
        let values: [[[String; 3]; 12]; 2] = std::array::from_fn(|r| {
            std::array::from_fn(|i| {
                [
                    "1.".into(),
                    i.to_string(),
                    if r == 0 {
                        "4.".into()
                    } else if factor {
                        format!("4.{i:016}")
                    } else {
                        "4.0000000000000001".into()
                    },
                ]
            })
        });
        assert_eq!(values[1][1][2].parse::<f64>().unwrap(), 4.);
        let tokens = std::array::from_fn(|r| {
            std::array::from_fn(|i| LiteralPoint {
                index: r * 12 + i,
                coordinates: std::array::from_fn(|axis| values[r][i][axis].as_str()),
            })
        });
        let mut memory = Memory::new(0, SCRATCH + 1024 * 1024).unwrap();
        let mut w = Work::new(16_000_000, None).unwrap();
        let mut memo = DecimalMemo::new(24, &mut memory, &mut w).unwrap();
        let [d, da, dd] =
            tensor_differentials(tokens, factor, &mut memo, &mut memory, &mut w).unwrap();
        let a = std::array::from_fn(|axis| {
            projection::polynomial::P::controls(&std::array::from_fn::<_, 12, _>(|i| {
                I::point(if axis == 0 {
                    1.
                } else if axis == 1 {
                    i as f64
                } else {
                    4.
                })
            }))
        });
        let j = projection::jacobian_differentials(a, d, da, dd, factor, &mut w).unwrap();
        // y=11t, with positive z gap h (or h*t before the owned
        // factor). The independently known signed determinant is -11h.
        let expected = if factor { 121e-16 } else { 11e-16 };
        for p in j {
            let h = p.hull();
            assert!(h.lo <= -expected && h.hi >= -expected);
            assert!(h.hi < 0.);
        }
        memo.release(&mut memory);
        assert_eq!(memory.used, 0);
    }
}
fn line(scale: &str) -> [[String; 2]; 12] {
    std::array::from_fn(|i| [format!("1E{scale}"), format!("{i}E{scale}")])
}
fn run(
    values: &[[String; 2]; 12],
    reverse: bool,
    work: usize,
    bytes: usize,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<f64, Error> {
    let tokens = std::array::from_fn(|i| std::array::from_fn(|a| values[i][a].as_str()));
    let mut m = Memory::new(0, bytes)?;
    let mut w = Work::new(work, interrupt)?;
    angular_lower(tokens, reverse, &mut m, &mut w)
}
#[test]
fn actual_decimal_polynomial_direction_and_numeric_limits_are_not_float_cache_claims() {
    let values = line("-15");
    let lower = run(&values, false, 16_000_000, SCRATCH, None).unwrap();
    assert!(lower > 0.);
    assert!(matches!(
        run(&values, true, 16_000_000, SCRATCH, None),
        Err(Error::Unresolved)
    ));
    let mut reversed = values.clone();
    reversed.reverse();
    assert!(run(&reversed, true, 16_000_000, SCRATCH, None).unwrap() > 0.);
    let constant = std::array::from_fn(|_| ["1.".into(), "0.".into()]);
    assert!(matches!(
        run(&constant, false, 16_000_000, SCRATCH, None),
        Err(Error::Unresolved)
    ));
    // Strictly positive real derivatives can be below the smallest binary64
    // bound; these valid decimal curves must refuse a claimed positive f64.
    assert!(matches!(
        run(&line("-170"), false, 16_000_000, SCRATCH, None),
        Err(Error::Unresolved)
    ));
    let huge = std::array::from_fn(|i| ["1E308".into(), format!("{i}E307")]);
    assert!(matches!(
        run(&huge, false, 16_000_000, SCRATCH, None),
        Err(Error::Unresolved)
    ));
}
#[test]
fn actual_proof_work_storage_cancellation_and_decimal_size_are_bounded() {
    let values = line("-15");
    assert!(matches!(
        run(&values, false, 16_000_000, SCRATCH - 1, None),
        Err(Error::Limit)
    ));
    assert!(matches!(
        run(&values, false, 1, SCRATCH, None),
        Err(Error::Limit)
    ));
    assert!(matches!(
        run(&values, false, 16_000_000, SCRATCH, Some(&|| true)),
        Err(Error::Cancelled)
    ));
    for token in [
        format!("1.{}", "0".repeat(511)),
        "1E-4097".into(),
        "1E-2147483649".into(),
        "1.0E-2147483648".into(),
        format!("{}E-4096", "1".repeat(506)),
    ] {
        assert!(token.parse::<f64>().unwrap().is_finite());
        let mut changed = values.clone();
        changed[0][0] = token;
        assert!(matches!(
            run(&changed, false, 16_000_000, SCRATCH, None),
            Err(Error::Limit)
        ));
    }
    let mut changed = values.clone();
    changed[0][0] = "not-a-number".into();
    assert!(matches!(
        run(&changed, false, 16_000_000, SCRATCH, None),
        Err(Error::Invalid)
    ));
    assert!(run(&values, false, 16_000_000, SCRATCH, None).unwrap() > 0.);
}
