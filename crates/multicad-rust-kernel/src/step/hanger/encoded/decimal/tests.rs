use super::*;
use std::cell::Cell;

fn exact_text(text: &str) -> BigRational {
    let (mantissa, exponent) = text.split_once("E-").unwrap_or((text, "0"));
    let n = BigInt::parse_bytes(mantissa.strip_suffix('.').unwrap().as_bytes(), 10).unwrap();
    BigRational::new(n, BigInt::from(10_u32).pow(exponent.parse().unwrap()))
}
fn run(value: &N, bytes: usize, limit: usize) -> Result<(FiniteDecimal, usize), Error> {
    let mut memory = Memory::new(19, bytes)?;
    let mut work = Work::new(limit, None)?;
    let text = decimal(value, &mut memory, &mut work)?;
    assert_eq!(memory.used, 19 + text.retained_bytes());
    assert_eq!(exact_text(text.as_str()), value.0);
    Ok((text, work.used))
}
#[test]
fn actual_finite_rationals_round_trip_without_binary64_conversion() {
    for twos in [0, 1, 7, 128, 1074] {
        for fives in [0, 1, 9, 64] {
            for numerator in [-37, 0, 1, 1234567] {
                let denominator = (BigInt::one() << twos) * BigInt::from(5_u32).pow(fives);
                let value = N(BigRational::new(numerator.into(), denominator));
                run(&value, 64 * 1024 * 1024, 16_000_000).unwrap();
            }
        }
    }
    assert_eq!(
        run(&N::integer(0), 1_000_000, 16_000_000)
            .unwrap()
            .0
            .as_str(),
        "0."
    );
    let q = N(BigRational::new(1.into(), 8.into()));
    assert_eq!(
        run(&q, 1_000_000, 16_000_000).unwrap().0.as_str(),
        "125.E-3"
    );
    let mut memory = Memory::new(0, 64 * 1024 * 1024).unwrap();
    let mut work = Work::new(16_000_000, None).unwrap();
    let scalar = EncodedScalar::retain(&q, &mut memory, &mut work).unwrap();
    let text = scalar
        .exact_step_real(DomainOptions::default(), None)
        .unwrap();
    assert_eq!(text.as_str(), "125.E-3");
    assert!(text.retained_bytes() >= size_of::<FiniteDecimal>() + text.as_str().len());
    assert!(matches!(
        scalar.exact_step_real(DomainOptions::default(), Some(&|| true)),
        Err(Error::Cancelled)
    ));
    assert!(matches!(
        scalar.exact_step_real(
            DomainOptions {
                maximum_work: 0,
                ..Default::default()
            },
            None
        ),
        Err(Error::Invalid)
    ));
    let calls = Cell::new(0);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    scalar
        .exact_step_real(DomainOptions::default(), Some(&count))
        .unwrap();
    for stop in 1..=calls.get() {
        let current = Cell::new(0);
        let cancel = || {
            current.set(current.get() + 1);
            current.get() == stop
        };
        assert!(matches!(
            scalar.exact_step_real(DomainOptions::default(), Some(&cancel)),
            Err(Error::Cancelled)
        ));
    }
}
#[test]
fn refusal_resources_and_every_cancellation_preserve_the_caller_ledger() {
    let q = N(BigRational::new(
        7.into(),
        (BigInt::one() << 128) * BigInt::from(5_u32).pow(9),
    ));
    let (_, used) = run(&q, 1_000_000, 16_000_000).unwrap();
    run(&q, 1_000_000, used).unwrap();
    assert!(matches!(run(&q, 1_000_000, used - 1), Err(Error::Limit)));
    let mut low: usize = 19;
    let mut high = 1_000_000;
    while low + 1 < high {
        let middle = low.midpoint(high);
        if run(&q, middle, 16_000_000).is_ok() {
            high = middle;
        } else {
            low = middle;
        }
    }
    run(&q, high, 16_000_000).unwrap();
    assert!(matches!(run(&q, high - 1, 16_000_000), Err(Error::Limit)));
    for invalid in [
        N(BigRational::new(1.into(), 3.into())),
        N(BigRational::new(7.into(), 30.into())),
    ] {
        assert!(matches!(
            run(&invalid, 1_000_000, 16_000_000),
            Err(Error::Unresolved)
        ));
    }
    let large = N(BigRational::from_integer(BigInt::one() << 16384));
    assert!(matches!(
        run(&large, 1_000_000, 16_000_000),
        Err(Error::Limit)
    ));
    let calls = Cell::new(0);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    let mut memory = Memory::new(19, 1_000_000).unwrap();
    decimal(
        &q,
        &mut memory,
        &mut Work::new(16_000_000, Some(&count)).unwrap(),
    )
    .unwrap();
    for stop in 1..=calls.get() {
        let current = Cell::new(0);
        let cancel = || {
            current.set(current.get() + 1);
            current.get() == stop
        };
        let mut memory = Memory::new(19, 1_000_000).unwrap();
        let result = decimal(
            &q,
            &mut memory,
            &mut Work::new(16_000_000, Some(&cancel)).unwrap(),
        );
        assert!(matches!(result, Err(Error::Cancelled)));
        assert_eq!(memory.used, 19);
    }
    run(&q, 1_000_000, 16_000_000).unwrap();
}
