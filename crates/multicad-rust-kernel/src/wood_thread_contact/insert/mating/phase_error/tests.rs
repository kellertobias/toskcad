use super::*;
use std::cell::Cell;
#[test]
fn actual_nonbinary_finite_stations_and_literal_evaluator_operations_are_enclosed() {
    for left in [false, true] {
        let (_, thread) = super::super::tests::fixture(left, true);
        let mut spec = thread.spec();
        spec.pitch_mm = 0.7;
        spec.length_mm = 19.3;
        spec.root_radius_mm = None;
        let actual = Thread::new(spec).unwrap();
        let proof = measure(&actual, &|| false).unwrap();
        assert!(proof.indexed_mm > 0. && proof.evaluator_mm > 0.);
        assert!(proof.total_mm >= proof.indexed_mm + proof.evaluator_mm);
        let p = q(spec.pitch_mm).unwrap();
        for k in [-6_i64, -1, 0, 1, 28, 34] {
            for knot in actual.profile_breaks() {
                let exact = Q::from_integer(k.into()) * &p + q(knot).unwrap();
                for value in [
                    (k as f64).mul_add(spec.pitch_mm, knot),
                    k as f64 * spec.pitch_mm + knot,
                ] {
                    let difference = (q(value).unwrap() - &exact).abs();
                    // Only the sealed finite indexed range is asserted. This
                    // oracle uses exact binary rationals, independent of I math.
                    let margin = 2 * i64::from(spec.starts) + 4;
                    let last = (spec.length_mm / spec.pitch_mm).ceil() as i64 + margin;
                    if (-margin..=last).contains(&k) {
                        assert!(difference <= q(proof.indexed_mm).unwrap());
                    }
                }
            }
        }
        let rate = actual.helical_rate();
        for z in [0., 0.1, 1.3, spec.length_mm.next_down(), spec.length_mm] {
            for theta in [
                -std::f64::consts::TAU.next_down(),
                -0.13,
                0.,
                0.37,
                std::f64::consts::TAU,
                std::f64::consts::TAU.next_up(),
            ] {
                let raw = q(z).unwrap() - q(rate).unwrap() * q(theta).unwrap();
                let exact = &raw - (&raw / &p).floor() * &p;
                let evaluated = (z.rem_euclid(spec.pitch_mm)
                    - rate * theta.rem_euclid(std::f64::consts::TAU))
                .rem_euclid(spec.pitch_mm);
                let distance = (q(evaluated).unwrap() - exact).abs();
                let circular = distance.clone().min((&p - distance).abs());
                assert!(circular <= q(proof.total_mm).unwrap(), "{z} {theta}");
            }
        }
        let calls = Cell::new(0);
        measure(&actual, &|| {
            calls.set(calls.get() + 1);
            false
        })
        .unwrap();
        for target in [1, 2, calls.get()] {
            let n = Cell::new(0);
            assert!(matches!(
                measure(&actual, &|| {
                    n.set(n.get() + 1);
                    n.get() >= target
                }),
                Err(Error::Cancelled)
            ));
        }
    }
    assert!(ulp_bound(ip(f64::MIN_POSITIVE)).unwrap().hi > 0.);
    assert!(ulp_bound(ip(f64::from_bits(1))).unwrap().hi > 0.);
    assert!(ulp_bound(ip(f64::NAN)).is_err());
    assert!(ulp_bound(ip(f64::MAX)).is_err());
}
