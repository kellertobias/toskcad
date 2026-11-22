use super::*;
use crate::authored_hanger_bolt::{BoundaryOptions, Hanger};
#[test]
fn source_and_actual_closed_field_volume_receipts_remain_distinct() {
    for (left, nb) in [(false, false), (true, true)] {
        let s = Hanger::new(crate::authored_hanger_bolt::tests::fixture(left, nb), None).unwrap();
        let b = s.boundary(BoundaryOptions::default(), None).unwrap();
        let f = EncodedField::prepare(&b, CurveOptions::default(), None).unwrap();
        let d = FiniteDomains::prepare(&f, DomainOptions::default(), None).unwrap();
        let c = ClosedField::prepare(&d, DomainOptions::default(), None).unwrap();
        let v = c
            .volume_bounds(1e-5, DomainOptions::default(), None)
            .unwrap();
        assert!(v.encoded_mm3[0] <= v.source_mm3[0] && v.encoded_mm3[1] >= v.source_mm3[1]);
        assert!(v.difference_bound_mm3 > 0. && v.encoded_mm3[1] - v.encoded_mm3[0] < 1e-5);
        eprintln!(
            "fieldvolume left={left} nb={nb} source={:?} encoded={:?} reserve={} work={}",
            v.source_mm3, v.encoded_mm3, v.difference_bound_mm3, v.work_used
        );
    }
}
fn choose(n: u64, k: u64) -> u64 {
    let mut v = 1;
    for j in 0..k {
        v = v * (n - j) / (j + 1);
    }
    v
}
#[test]
fn independent_actual_bernstein_area_and_rational_meridian_integral_fit_the_enclosure() {
    use super::super::super::cap::tests::number;
    use num_rational::BigRational as Q;
    use num_traits::{ToPrimitive, Zero};
    let s = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, false),
        None,
    )
    .unwrap();
    let b = s.boundary(BoundaryOptions::default(), None).unwrap();
    let f = EncodedField::prepare(&b, CurveOptions::default(), None).unwrap();
    let d = FiniteDomains::prepare(&f, DomainOptions::default(), None).unwrap();
    let c = ClosedField::prepare(&d, DomainOptions::default(), None).unwrap();
    let mut area = Q::zero();
    // Exact integral of actual degree10 Bernstein B cross B'. Physical
    // piece widths cancel its derivative scale, so no source π is inserted.
    for piece in f.angular().coefficients() {
        for i in 0..=10 {
            for j in 0..10 {
                let dx = &piece[j + 1][0] - &piece[j][0];
                let dy = &piece[j + 1][1] - &piece[j][1];
                let det = &piece[i][0] * dy - &piece[i][1] * dx;
                area += det
                    * Q::new(
                        (10 * choose(10, i as u64) * choose(9, j as u64)).into(),
                        (40 * choose(19, (i + j) as u64)).into(),
                    );
            }
        }
    }
    let area = area.to_f64().unwrap();
    assert!((area - std::f64::consts::PI).abs() < 1e-7);
    let mut integral = [0_f64; 2];
    for meridian in f.meridians() {
        let q = meridian
            .phase_coefficients()
            .iter()
            .map(|v| number(v).to_f64().unwrap())
            .collect::<Vec<_>>();
        let r = meridian
            .radial_coefficients()
            .iter()
            .map(|v| number(v).to_f64().unwrap())
            .collect::<Vec<_>>();
        let w = meridian
            .weights()
            .iter()
            .map(|v| number(v).to_f64().unwrap())
            .collect::<Vec<_>>();
        let density = |u: f64| {
            let evaluate = |c: &[f64]| {
                if c.len() == 2 {
                    c[0] * (1. - u) + c[1] * u
                } else {
                    c[0] * (1. - u) * (1. - u) + 2. * c[1] * u * (1. - u) + c[2] * u * u
                }
            };
            let derivative = |c: &[f64]| {
                if c.len() == 2 {
                    c[1] - c[0]
                } else {
                    2. * ((c[1] - c[0]) * (1. - u) + (c[2] - c[1]) * u)
                }
            };
            let weight = evaluate(&w);
            let radius = evaluate(&r) / weight;
            radius * radius * (derivative(&q) * weight - evaluate(&q) * derivative(&w))
                / (weight * weight)
        };
        let count: usize = 4096;
        let mut value = density(0.) + density(1.);
        for i in 1..count {
            value += if i.is_multiple_of(2) { 2. } else { 4. } * density(i as f64 / count as f64);
        }
        integral[usize::from(meridian.component() == Component::Wood)] +=
            value / (3. * count as f64);
    }
    for (i, length, pitch) in [
        (0, s.metric().full_length_mm(), s.metric().spec().pitch_mm),
        (1, s.wood().full_length_mm(), s.wood().spec().pitch_mm),
    ] {
        assert!((Q::from_float(length).unwrap() / Q::from_float(pitch).unwrap()).is_integer());
        integral[i] *= length / pitch;
    }
    let plain = s.spec();
    let numerical =
        area * (integral[0] + integral[1] + plain.plain_radius_mm.powi(2) * plain.plain_length_mm);
    let receipt = c
        .volume_bounds(1e-5, DomainOptions::default(), None)
        .unwrap();
    assert!(receipt.encoded_mm3[0] < numerical && numerical < receipt.encoded_mm3[1]);
    eprintln!(
        "actual chosen-scalar field numerical volume={numerical}; bounded {:?}",
        receipt.encoded_mm3
    );
}
#[test]
fn volume_precision_and_all_observed_publication_checkpoints_remain_atomic() {
    let s = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, false),
        None,
    )
    .unwrap();
    let b = s.boundary(BoundaryOptions::default(), None).unwrap();
    let f = EncodedField::prepare(&b, CurveOptions::default(), None).unwrap();
    let d = FiniteDomains::prepare(&f, DomainOptions::default(), None).unwrap();
    let c = ClosedField::prepare(&d, DomainOptions::default(), None).unwrap();
    let o = DomainOptions::default();
    assert!(matches!(c.volume_bounds(0., o, None), Err(Error::Invalid)));
    assert!(matches!(
        c.volume_bounds(f64::from_bits(1), o, None),
        Err(Error::Unresolved)
    ));
    assert!(matches!(
        c.volume_bounds(1e-14, o, None),
        Err(Error::Unresolved)
    ));
    assert!(matches!(
        c.volume_bounds(1e-8, o, None),
        Err(Error::Unresolved)
    ));
    assert!(matches!(
        c.volume_bounds(
            1e-5,
            DomainOptions {
                maximum_work: 1,
                ..o
            },
            None
        ),
        Err(Error::Limit)
    ));
    let n = std::cell::Cell::new(0);
    c.volume_bounds(
        1e-5,
        o,
        Some(&|| {
            n.set(n.get() + 1);
            false
        }),
    )
    .unwrap();
    for prefix in 0..n.get() {
        let seen = std::cell::Cell::new(0);
        assert!(matches!(
            c.volume_bounds(
                1e-5,
                o,
                Some(&|| {
                    seen.set(seen.get() + 1);
                    seen.get() > prefix
                })
            ),
            Err(Error::Cancelled)
        ));
    }
    let base = crate::authored_hanger_bolt::QueryOptions {
        volume_error_mm3: 1e-5,
        maximum_cells: 11,
    };
    assert_eq!(
        source_error(
            s.volume(
                crate::authored_hanger_bolt::QueryOptions {
                    volume_error_mm3: 0.,
                    ..base
                },
                None
            )
            .unwrap_err()
        ),
        Error::Invalid
    );
    assert_eq!(
        source_error(
            s.volume(
                crate::authored_hanger_bolt::QueryOptions {
                    maximum_cells: 1,
                    ..base
                },
                None
            )
            .unwrap_err()
        ),
        Error::Limit
    );
    assert_eq!(
        source_error(s.volume(base, Some(&|| true)).unwrap_err()),
        Error::Cancelled
    );
    let spec = s.spec();
    assert_eq!(
        source_error(
            Hanger::new(
                crate::authored_hanger_bolt::Spec {
                    plain_radius_mm: spec.metric.maximum_radius_mm(),
                    ..spec
                },
                None
            )
            .unwrap_err()
        ),
        Error::Unresolved
    );
}
