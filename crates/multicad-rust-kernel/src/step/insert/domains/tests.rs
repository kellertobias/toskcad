use super::*;
fn graph(left: bool, starts: u32) -> Boundary {
    super::super::tests::graph(left, starts)
}
#[test]
fn actual_paired_cells_share_one_exact_angular_station_set_and_original_curve_rows() {
    for (left, starts) in [(false, 1), (true, 1), (false, 2)] {
        let body = graph(left, starts);
        let plan = DomainPlan::prepare(&body, DomainOptions::default(), None).unwrap();
        assert!(std::ptr::eq(plan.boundary(), &body));
        for side in [Side::Outer, Side::Inner] {
            let mut work = Work::new(16_000_000, None).unwrap();
            for pair in plan.stations(side).windows(2) {
                assert_eq!(
                    pair[0]
                        .number(&mut work)
                        .unwrap()
                        .cmp(&pair[1].number(&mut work).unwrap(), &mut work)
                        .unwrap(),
                    Ordering::Less
                );
            }
            for cell in 0..body.cells(side).len() {
                let slabs = plan
                    .slabs()
                    .iter()
                    .filter(|s| s.side == side && s.cell == cell)
                    .collect::<Vec<_>>();
                assert!(!slabs.is_empty());
                for pair in slabs.windows(2) {
                    assert_eq!(pair[0].stations[1], pair[1].stations[0]);
                }
                assert_eq!(
                    plan.exact_cell(side, cell).unwrap().row_count(),
                    body.cells(side)[cell].points().len()
                );
                for s in slabs {
                    assert_ne!(s.rows[0], s.rows[1]);
                    assert!(
                        s.rows
                            .into_iter()
                            .all(|key| key.owner() == side && body.curve(key).is_some())
                    );
                    assert_ne!(s.collapsed, [true, true]);
                }
            }
            assert!(plan.exact_cell(side, usize::MAX).is_none());
        }
        assert!(plan.work_used() <= 16_000_000);
        assert!(plan.retained_bytes() <= 64 * 1024 * 1024);
        eprintln!(
            "domain plan {left}/{starts}: {}+{} stations, {} slabs, {} work",
            plan.stations(Side::Outer).len(),
            plan.stations(Side::Inner).len(),
            plan.slabs().len(),
            plan.work_used()
        );
    }
}
#[test]
fn domain_resources_and_cancellation_refuse_without_publication() {
    let body = graph(false, 1);
    for o in [
        DomainOptions {
            maximum_stations: 0,
            ..DomainOptions::default()
        },
        DomainOptions {
            maximum_stations: 16385,
            ..DomainOptions::default()
        },
        DomainOptions {
            maximum_slabs: 0,
            ..DomainOptions::default()
        },
        DomainOptions {
            maximum_slabs: 131073,
            ..DomainOptions::default()
        },
        DomainOptions {
            bytes: 0,
            ..DomainOptions::default()
        },
        DomainOptions {
            bytes: 64 * 1024 * 1024 + 1,
            ..DomainOptions::default()
        },
    ] {
        assert!(matches!(
            DomainPlan::prepare(&body, o, None),
            Err(Error::Invalid)
        ));
    }
    for o in [
        DomainOptions {
            maximum_stations: 1,
            ..DomainOptions::default()
        },
        DomainOptions {
            maximum_slabs: 1,
            ..DomainOptions::default()
        },
        DomainOptions {
            work: 1,
            ..DomainOptions::default()
        },
        DomainOptions {
            bytes: body.retained_bytes(),
            ..DomainOptions::default()
        },
    ] {
        assert!(matches!(
            DomainPlan::prepare(&body, o, None),
            Err(Error::Limit)
        ));
    }
    assert!(matches!(
        DomainPlan::prepare(&body, DomainOptions::default(), Some(&|| true)),
        Err(Error::Cancelled)
    ));
    let polls = std::cell::Cell::new(0);
    let stop = || {
        polls.set(polls.get() + 1);
        polls.get() > 100
    };
    assert!(matches!(
        DomainPlan::prepare(&body, DomainOptions::default(), Some(&stop)),
        Err(Error::Cancelled)
    ));
}
fn raw(
    v: crate::authored_internal_treatment::ExactBoundaryRational<'_>,
) -> num_rational::BigRational {
    use num_bigint::{BigInt, Sign};
    let sign = match v.signum() {
        -1 => Sign::Minus,
        0 => Sign::NoSign,
        _ => Sign::Plus,
    };
    num_rational::BigRational::new(
        BigInt::from_bytes_le(sign, v.numerator_magnitude_le()),
        BigInt::from_bytes_le(Sign::Plus, v.denominator_magnitude_le()),
    )
}
#[test]
fn retained_curve_rows_match_actual_source_band_planes_over_the_entire_slab() {
    use num_rational::BigRational as Q;
    for left in [false, true] {
        let body = graph(left, 1);
        let plan = DomainPlan::prepare(&body, DomainOptions::default(), None).unwrap();
        let mut work = Work::new(16_000_000, None).unwrap();
        for side in [Side::Outer, Side::Inner] {
            for station in plan.stations(side) {
                let sign = match station.signum() {
                    -1 => num_bigint::Sign::Minus,
                    0 => num_bigint::Sign::NoSign,
                    _ => num_bigint::Sign::Plus,
                };
                assert_eq!(
                    station.number(&mut work).unwrap().0,
                    Q::new(
                        num_bigint::BigInt::from_bytes_le(sign, station.numerator_magnitude_le()),
                        num_bigint::BigInt::from_bytes_le(
                            num_bigint::Sign::Plus,
                            station.denominator_magnitude_le()
                        )
                    )
                );
            }
        }
        for slab in plan.slabs() {
            let source = body.source_boundary(slab.side);
            let cell = &source.cells()[slab.cell];
            let spec = source.source().spec();
            let band = source.source().bands()[cell.band()];
            let f = |x| Q::from_float(x).unwrap();
            let pitch = f(spec.pitch_mm);
            let lead = &pitch * Q::from_integer(spec.starts.into());
            let q0 = f(band.phase_mm[0]) + &pitch * Q::from_integer(cell.period().into());
            let slope = (f(band.radius_mm[1]) - f(band.radius_mm[0]))
                / (f(band.phase_mm[1]) - f(band.phase_mm[0]));
            let intercept = f(band.radius_mm[0]) - &slope * q0;
            let a = plan.stations(slab.side)[slab.stations[0] as usize]
                .number(&mut work)
                .unwrap()
                .0;
            let b = plan.stations(slab.side)[slab.stations[1] as usize]
                .number(&mut work)
                .unwrap()
                .0;
            let mut values = Vec::new();
            for row in slab.rows {
                let c = source
                    .exact_curve_coefficients(row.source_index(), 16_000_000, 0, None)
                    .unwrap();
                let base = std::array::from_fn::<_, 3, _>(|k| raw(c.coordinate(0, k).unwrap()));
                let delta = std::array::from_fn::<_, 3, _>(|k| raw(c.coordinate(1, k).unwrap()));
                assert_ne!(delta[2], Q::from_integer(0.into()));
                let mut samples = Vec::new();
                for v in [&a, &((&a + &b) / Q::from_integer(2.into())), &b] {
                    let t = (v - &base[2]) / &delta[2];
                    assert!(t >= Q::from_integer(0.into()) && t <= Q::from_integer(1.into()));
                    let r = &base[0] + &delta[0] * &t;
                    let z = &base[1] + &delta[1] * t;
                    let q = z - &lead * v;
                    assert_eq!(r, &slope * &q + &intercept);
                    samples.push(q);
                }
                values.push(samples);
            }
            assert!(values[0][1] < values[1][1]);
            for i in [0, 2] {
                assert!(values[0][i] <= values[1][i]);
                assert_eq!(values[0][i] == values[1][i], slab.collapsed[i / 2]);
            }
        }
    }
}
#[test]
fn bounded_rational_primitive_and_source_error_mapping_preserve_failures() {
    for (input, expected) in [
        (SourceError::Invalid, Error::Invalid),
        (SourceError::Unresolved, Error::Unresolved),
        (SourceError::Limit, Error::Limit),
        (SourceError::Cancelled, Error::Cancelled),
    ] {
        assert_eq!(failure(input), expected);
    }
    let mut work = Work::new(16_000_000, None).unwrap();
    let one = N::integer(1);
    let zero = N::integer(0);
    assert!(matches!(one.div(&zero, &mut work), Err(Error::Unresolved)));
    let huge = N(num_rational::BigRational::from_integer(
        num_bigint::BigInt::from(1) << 16384,
    ));
    assert!(matches!(huge.add(&one, &mut work), Err(Error::Limit)));
}
