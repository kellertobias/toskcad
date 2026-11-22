use crate::{Error, GridCoordinate, Quantum, Rational, parse_length_mm};

#[test]
fn unary_negation_preserves_exact_normalization_and_cancels() {
    let zero = Rational::new(0, 1).unwrap();
    for (numerator, denominator) in [(0, 7), (1, 3), (-7, 2), (25, 10)] {
        let value = Rational::new(numerator, denominator).unwrap();
        let negated = -value;
        let (n, d) = value.ratio();
        assert_eq!(negated.ratio(), (-n, d));
        assert_eq!(-negated, value);
        assert_eq!(negated.checked_add(value).unwrap(), zero);
    }
    assert_eq!((-zero).ratio(), (0, 1));
}
fn q(n: i64, d: i64) -> Rational {
    Rational::new(n, d).unwrap()
}
#[test]
fn equivalent_fraction_unit_positions_keep_exact_values_and_reject_invalid_divisors() {
    for (a, b, expected) in [
        ("1/2 mm", "1 mm/2", q(1, 2)),
        ("1/2 in", "1 in/2", q(127, 10)),
    ] {
        assert_eq!(parse_length_mm(a).unwrap(), expected);
        assert_eq!(parse_length_mm(b).unwrap(), expected);
    }
    for text in ["1/0 mm", "1 mm/0"] {
        assert_eq!(parse_length_mm(text), Err(Error::DivisionByZero));
    }
    for text in ["x/2 mm", "1/x mm", "x mm/2", "1 mm/x"] {
        assert_eq!(parse_length_mm(text), Err(Error::InvalidLength));
    }
}
#[test]
fn grid_classification_never_rounds_off_grid_lengths_or_accepts_nonpositive_steps() {
    for step in [q(0, 1), q(-1, 2)] {
        assert_eq!(Quantum::new(step), Err(Error::NonPositiveQuantum));
    }
    let grid = Quantum::new(q(1, 2)).unwrap();
    for (length, ticks) in [(q(0, 1), 0), (q(3, 2), 3), (q(-3, 2), -3)] {
        assert_eq!(grid.classify(length).unwrap(), GridCoordinate::Ticks(ticks));
    }
    let authored = q(2, 3);
    assert_eq!(
        grid.classify(authored).unwrap(),
        GridCoordinate::OffGrid(authored)
    );
}
#[test]
fn exact_tick_overflow_refuses_instead_of_saturating_positive_or_negative_lengths() {
    let grid = Quantum::new(q(1, 10)).unwrap();
    for sign in [-1, 1] {
        let last = sign * 922_337_203_685_477_580;
        assert_eq!(
            grid.classify(q(last, 1)).unwrap(),
            GridCoordinate::Ticks(last * 10)
        );
        assert_eq!(grid.classify(q(last + sign, 1)), Err(Error::OutOfRange));
    }
}
