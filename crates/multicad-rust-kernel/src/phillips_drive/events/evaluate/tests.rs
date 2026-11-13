use super::*;
#[test]
fn numeric_values_and_jets_reject_nonfinite_and_unknown_roots() {
    assert!(Value::checked(I::new(1., 0.)).is_err());
    assert!(Value::checked(Series::<5>::interval(I::new(1., 0.))).is_err());
    assert!(root(I::new(0., 1.)).is_err());
    assert!(root(Series::<5>::interval(I::new(0., 1.))).is_err());
    let [sin] = [Value::sine(I::point(0.3))];
    assert!(sin.lo <= 0.3_f64.sin() && sin.hi >= 0.3_f64.sin());
}
