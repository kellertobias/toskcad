use super::*;
#[test]
fn actual_opt_in_plan_retains_authored_source_and_original_tiny_curve_parameters() {
    for (left, length) in [
        (false, 4.),
        (false, 4_f64.next_up()),
        (true, 4_f64.next_up()),
        (false, 4_f64.next_down()),
        (true, 8.),
    ] {
        let source = crate::step_read::owned_insert::tests::graph(left, length);
        let plan = CoupledPlan::prepare(&source, CurveOptions::default(), None).unwrap();
        assert!(std::ptr::eq(plan.source(), &source));
        assert!(std::ptr::eq(plan.domain().boundary(), &source));
        assert!(plan.work_used() < 16_000_000 && plan.retained_bytes() < 64 * 1024 * 1024);
        assert!(plan.rows().len() > 500);
        for row in plan.rows() {
            assert!(row.position_error_mm() <= 1e-7);
            assert!(row.original_curve_derivative_error_mm() <= 1e-6);
            assert_eq!(row.coefficients().len(), 12);
            for point in row.coefficients() {
                for c in point {
                    assert_eq!(c.fractional_bits(), 80);
                    assert!(!c.numerator_magnitude_le().is_empty());
                    let [a, b] = c.numeric_bounds_mm();
                    assert!(a.is_finite() && b.is_finite() && a <= b);
                }
            }
        }
        let first = &plan.rows()[0].coefficients()[0][0];
        let decimal = first.decimal(16_000_000, 64 * 1024 * 1024, None).unwrap();
        assert!(decimal.contains('.'));
        eprintln!(
            "opt-in coupled {left}/{length:?}: {}rows work {} bytes {}",
            plan.rows().len(),
            plan.work_used(),
            plan.retained_bytes()
        );
    }
}
fn rational(parameter: &Parameter) -> BigRational {
    BigRational::new(
        BigInt::from_bytes_le(
            if parameter.signum() < 0 {
                Sign::Minus
            } else if parameter.signum() > 0 {
                Sign::Plus
            } else {
                Sign::NoSign
            },
            parameter.numerator_magnitude_le(),
        ),
        BigInt::from_bytes_le(Sign::Plus, parameter.denominator_magnitude_le()),
    )
}
fn coefficient(c: &Coefficient) -> BigRational {
    BigRational::new(
        BigInt::from_bytes_le(
            if c.signum() < 0 {
                Sign::Minus
            } else if c.signum() > 0 {
                Sign::Plus
            } else {
                Sign::NoSign
            },
            c.numerator_magnitude_le(),
        ),
        BigInt::one() << c.fractional_bits(),
    )
}
fn decimal(text: &str) -> BigRational {
    let negative = text.starts_with('-');
    let text = text.strip_prefix('-').unwrap_or(text);
    let (a, b) = text.split_once('.').unwrap();
    let n = BigInt::parse_bytes(format!("{a}{b}").as_bytes(), 10).unwrap();
    BigRational::new(
        if negative { -n } else { n },
        BigInt::from(10).pow(b.len() as u32),
    )
}
#[test]
fn actual_original_parameters_decimal_coefficients_and_tiny_derivatives_keep_source_authority() {
    let body = crate::step_read::owned_insert::tests::graph(false, 4_f64.next_up());
    let plan = CoupledPlan::prepare(&body, CurveOptions::default(), None).unwrap();
    let mut signs = [false; 3];
    let mut tiny = 0;
    for row in plan.rows() {
        let key = row.key();
        assert!(body.curve(key.curve).is_some());
        let parameters = row.original_parameters();
        let a = rational(parameters[0]);
        let b = rational(parameters[1]);
        assert!(
            a >= BigRational::zero()
                && a <= BigRational::one()
                && b >= BigRational::zero()
                && b <= BigRational::one()
        );
        let gain = b - a;
        assert!(!gain.is_zero());
        for owner in row.endpoint_owners() {
            match owner {
                surfaces::VertexOwner::Original(key) => {
                    assert_eq!(
                        match key {
                            BoundaryVertexKey::Outer(_) => Side::Outer,
                            BoundaryVertexKey::Inner(_) => Side::Inner,
                        },
                        row.key().curve.owner()
                    )
                }
                surfaces::VertexOwner::Station { curve, .. } => assert_eq!(curve, row.key().curve),
            }
        }
        for point in row.coefficients() {
            for c in point {
                let index = (c.signum() + 1) as usize;
                if !signs[index] {
                    let text = c.decimal(16_000_000, 64 * 1024 * 1024, None).unwrap();
                    assert_eq!(decimal(&text), coefficient(c));
                    signs[index] = true;
                }
            }
        }
        if gain.abs() < BigRational::new(1.into(), 1_000_000_000_000_i64.into()) {
            let component = body.source_boundary(key.curve.owner());
            let data = component
                .exact_curve_coefficients(key.curve.source_index(), 16_000_000, 0, None)
                .unwrap();
            let mut w = Work::new(16_000_000, None).unwrap();
            let source = N::view(data.coordinate(1, 1).unwrap(), &mut w).unwrap();
            for pair in row.coefficients().windows(2) {
                let derivative = (coefficient(&pair[1][2]) - coefficient(&pair[0][2]))
                    * BigInt::from(11)
                    / &gain;
                let error = (derivative - &source.0).abs();
                assert!(
                    error
                        <= BigRational::from_float(row.original_curve_derivative_error_mm())
                            .unwrap()
                );
            }
            tiny += 1;
        }
    }
    assert_eq!(signs, [true; 3]);
    assert!(tiny >= 10);
}
#[test]
fn actual_source_plan_resource_cancel_precision_refusals_allow_unchanged_retry() {
    use std::cell::Cell;
    let body = crate::step_read::owned_insert::tests::graph(false, 4.);
    let callbacks = Cell::new(0);
    let count = || {
        callbacks.set(callbacks.get() + 1);
        false
    };
    let baseline = CoupledPlan::prepare(&body, CurveOptions::default(), Some(&count)).unwrap();
    let total = callbacks.get();
    for point in (0..128).map(|i| total * i / 128).chain([total - 1]) {
        let calls = Cell::new(0);
        let stop = || {
            let n = calls.get();
            calls.set(n + 1);
            n >= point
        };
        assert!(matches!(
            CoupledPlan::prepare(&body, CurveOptions::default(), Some(&stop)),
            Err(Error::Cancelled)
        ));
    }
    assert!(matches!(
        CoupledPlan::prepare(
            &body,
            CurveOptions {
                work: baseline.work_used() - 1,
                ..CurveOptions::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    // Measured actual default-source prerequisite frontier: 216726 weighted
    // operations admit the original angular/domain proofs but cannot retain
    // the next shared-row cache request. This is a real caller budget refusal.
    assert!(matches!(
        CoupledPlan::prepare(
            &body,
            CurveOptions {
                work: 216726,
                ..Default::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    for bytes in [body.retained_bytes() - 1, baseline.retained_bytes()] {
        assert!(matches!(
            CoupledPlan::prepare(
                &body,
                CurveOptions {
                    bytes,
                    ..CurveOptions::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
    }
    assert!(matches!(
        CoupledPlan::prepare(
            &body,
            CurveOptions {
                maximum_spans_per_curve: 1,
                position_error_mm: 1e-10,
                derivative_error_mm: 1e-10,
                ..CurveOptions::default()
            },
            None
        ),
        Err(Error::Unresolved)
    ));
    assert!(matches!(
        CoupledPlan::prepare(
            &body,
            CurveOptions {
                work: 0,
                ..CurveOptions::default()
            },
            None
        ),
        Err(Error::Invalid)
    ));
    let c = &baseline.rows()[0].coefficients()[0][0];
    assert!(matches!(
        c.decimal(1, 64 * 1024 * 1024, None),
        Err(Error::Limit)
    ));
    assert!(matches!(c.decimal(16_000_000, 1, None), Err(Error::Limit)));
    assert!(matches!(
        c.decimal(16_000_000, 0, None),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        c.decimal(16_000_000, 64 * 1024 * 1024 + 1, None),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        c.decimal(16_000_000, 64 * 1024 * 1024, Some(&|| true)),
        Err(Error::Cancelled)
    ));
    let retry = CoupledPlan::prepare(&body, CurveOptions::default(), None).unwrap();
    assert_eq!(baseline.rows().len(), retry.rows().len());
    for (a, b) in baseline.rows().iter().zip(retry.rows()) {
        assert_eq!(a.key(), b.key());
        for (a, b) in a
            .coefficients()
            .iter()
            .flatten()
            .zip(b.coefficients().iter().flatten())
        {
            assert_eq!(a.numerator_magnitude_le(), b.numerator_magnitude_le());
            assert_eq!(a.signum(), b.signum());
        }
    }
}
#[test]
fn actual_nonbinary_multistart_sources_keep_exact_owned_coupling() {
    for (left, starts) in [(false, 2), (true, 4), (false, 3), (true, 3)] {
        let body = super::super::tests::graph(left, starts);
        let plan = CoupledPlan::prepare(&body, CurveOptions::default(), None).unwrap();
        assert!(std::ptr::eq(plan.source(), &body));
        assert!(plan.work_used() < 16_000_000);
        eprintln!(
            "coupled nonbinary {left}/{starts}: {}rows {}work {}bytes",
            plan.rows().len(),
            plan.work_used(),
            plan.retained_bytes()
        );
    }
}
#[test]
fn actual_fractional_and_integer_signed_coordinates_format_as_the_same_exact_dyadic() {
    let body = crate::step_read::owned_insert::tests::graph(true, 4.);
    let plan = CoupledPlan::prepare(&body, CurveOptions::default(), None).unwrap();
    let mut seen = [false; 5];
    for c in plan.rows().iter().flat_map(|r| r.coefficients()).flatten() {
        let value = coefficient(c);
        let index = if value.is_zero() {
            0
        } else if value.is_integer() {
            1
        } else if value.abs() > BigRational::one() {
            2
        } else if value.is_negative() {
            3
        } else {
            4
        };
        if !seen[index] {
            assert_eq!(
                decimal(&c.decimal(16_000_000, 64 * 1024 * 1024, None).unwrap()),
                value
            );
            seen[index] = true;
        }
        if seen.iter().all(|v| *v) {
            break;
        }
    }
    assert_eq!(seen, [true; 5]);
}
#[test]
fn actual_point_inventory_rejects_incomplete_or_wrong_literal_aliases_without_changing_source() {
    let body = crate::step_read::owned_insert::tests::graph(false, 4.);
    let plan = CoupledPlan::prepare(&body, CurveOptions::default(), None).unwrap();
    let omitted = plan
        .rows()
        .iter()
        .flat_map(|r| r.endpoint_owners())
        .find_map(|o| {
            if let surfaces::VertexOwner::Original(k) = o {
                Some(k)
            } else {
                None
            }
        })
        .unwrap();
    for omit in [true, false] {
        let mut memory = Memory::new(plan.retained_bytes(), 64 * 1024 * 1024).unwrap();
        let mut w = Work::new(16_000_000, None).unwrap();
        let mut registry = aliases::Registry::prepare(&body, &mut memory, &mut w).unwrap();
        for (index, row) in plan.rows().iter().enumerate() {
            for (endpoint, owner) in row.endpoint_owners().iter().enumerate() {
                if let surfaces::VertexOwner::Original(vertex) = owner
                    && (!omit || *vertex != omitted)
                {
                    registry
                        .bind(*vertex, row, index, endpoint, plan.rows(), &mut w)
                        .unwrap();
                }
            }
        }
        assert_eq!(
            registry.complete(&mut w),
            if omit { Err(Error::Unresolved) } else { Ok(()) }
        );
        // A real alias can be registered repeatedly by the same source row;
        // a different occupied row cannot rebind it to unrelated coefficients.
        let (row_index, endpoint, vertex) = plan
            .rows()
            .iter()
            .enumerate()
            .find_map(|(i, r)| {
                r.endpoint_owners().iter().enumerate().find_map(|(e, o)| {
                    if let surfaces::VertexOwner::Original(k) = o {
                        Some((i, e, *k))
                    } else {
                        None
                    }
                })
            })
            .unwrap();
        let mut separate = aliases::Registry::prepare(&body, &mut memory, &mut w).unwrap();
        let actual_row = &plan.rows()[row_index];
        separate
            .bind(vertex, actual_row, row_index, endpoint, plan.rows(), &mut w)
            .unwrap();
        separate
            .bind(vertex, actual_row, row_index, endpoint, plan.rows(), &mut w)
            .unwrap();
        let wrong = plan
            .rows()
            .iter()
            .enumerate()
            .find(|(_, r)| {
                r.coefficients()[0]
                    .iter()
                    .zip(&actual_row.coefficients()[endpoint * 11])
                    .any(|(a, b)| {
                        a.signum() != b.signum()
                            || a.numerator_magnitude_le() != b.numerator_magnitude_le()
                    })
            })
            .unwrap();
        assert_eq!(
            separate.bind(vertex, wrong.1, wrong.0, 0, plan.rows(), &mut w),
            Err(Error::Unresolved)
        );
        separate.release(&mut memory);
        let row = &plan.rows()[0];
        assert!(matches!(
            registry.bind(
                BoundaryVertexKey::Outer(u32::MAX),
                row,
                0,
                0,
                plan.rows(),
                &mut w
            ),
            Err(Error::Invalid)
        ));
        registry.release(&mut memory);
    }
    let first = &plan.rows()[0].coefficients()[0];
    let different = plan
        .rows()
        .iter()
        .flat_map(|r| r.coefficients())
        .find(|p| {
            p.iter().zip(first).any(|(a, b)| {
                a.signum() != b.signum() || a.numerator_magnitude_le() != b.numerator_magnitude_le()
            })
        })
        .unwrap();
    let mut w = Work::new(16_000_000, None).unwrap();
    assert_eq!(same_point(first, different, &mut w), Err(Error::Unresolved));
    assert_eq!(same_point(first, first, &mut w), Ok(()));
}
#[test]
fn real_domain_capabilities_reject_wrong_source_and_missing_angular_refinement() {
    let first = crate::step_read::owned_insert::tests::graph(false, 4.);
    let second = crate::step_read::owned_insert::tests::graph(false, 4_f64.next_up());
    let d1 = DomainPlan::prepare(&first, DomainOptions::default(), None).unwrap();
    let d2 = DomainPlan::prepare(&second, DomainOptions::default(), None).unwrap();
    let mut memory =
        Memory::new(d1.retained_bytes() + d2.retained_bytes(), 64 * 1024 * 1024).unwrap();
    memory.add(SCRATCH).unwrap();
    let mut w = Work::new(16_000_000, None).unwrap();
    let r1 = admitted::Registry::prepare(&d1, &mut memory, &mut w).unwrap();
    let r2 = admitted::Registry::prepare(&d2, &mut memory, &mut w).unwrap();
    let a = r1.iter().next().unwrap();
    let b = r2.iter().next().unwrap();
    let curve = admitted::Curve::load(&a, &mut memory, &mut w, None).unwrap();
    assert!(matches!(b.original(&curve, &mut w), Err(Error::Invalid)));
    let unit = super::super::reduced_unit::fit_reduced_unit(
        false,
        CurveOptions {
            position_error_mm: 1e-9,
            derivative_error_mm: 1e-8,
            ..CurveOptions::default()
        },
        &mut memory,
        &mut w,
    )
    .unwrap();
    // Genuine capabilities from two immutable source arenas cannot be mixed.
    // A valid same-source row also refuses a tighter real coefficient budget.
    let other_curve = admitted::Curve::load(&b, &mut memory, &mut w, None).unwrap();
    let other_original = b.original(&other_curve, &mut w).unwrap();
    let refined =
        DomainPlan::prepare_with_knots(&first, unit.knots(), DomainOptions::default(), None)
            .unwrap();
    memory
        .add(refined.retained_bytes() - first.retained_bytes())
        .unwrap();
    let refined_registry = admitted::Registry::prepare(&refined, &mut memory, &mut w).unwrap();
    let mut checked = false;
    for request in refined_registry.iter() {
        if let Ok(span) = request.angular(&unit, &mut w) {
            let angular = rows::angular_slice(span, &mut w).unwrap();
            assert!(matches!(
                rows::row(&angular, &other_original, CurveOptions::default(), &mut w),
                Err(Error::Invalid)
            ));
            let own_curve = admitted::Curve::load(&request, &mut memory, &mut w, None).unwrap();
            let original = request.original(&own_curve, &mut w).unwrap();
            assert!(matches!(
                rows::row(
                    &angular,
                    &original,
                    CurveOptions {
                        position_error_mm: 1e-30,
                        derivative_error_mm: 1e-30,
                        ..Default::default()
                    },
                    &mut w
                ),
                Err(Error::Unresolved)
            ));
            checked = true;
            break;
        }
    }
    assert!(checked);
    let mut refused = false;
    for row in r1.iter() {
        if matches!(row.angular(&unit, &mut w), Err(Error::Unresolved)) {
            refused = true;
            break;
        }
    }
    assert!(refused);
}

#[test]
fn exact_grid_numeric_domain_and_real_domain_station_limits_are_bounded() {
    let mut w = Work::new(16_000_000, None).unwrap();
    let huge = N(BigRational::new(
        BigInt::from(1) << 16304_usize,
        BigInt::from(1),
    ));
    assert!(matches!(
        number::Grid::rational(&huge, BITS, &mut w),
        Err(Error::Limit)
    ));
    let ordinary = N(BigRational::new(BigInt::from(3), BigInt::from(2)));
    assert_eq!(
        number::Grid::rational(&ordinary, 1, &mut w).unwrap().0,
        BigInt::from(3)
    );
    let body = crate::step_read::owned_insert::tests::graph(false, 4.);
    for knots in [
        vec![0., 0.],
        vec![0., f64::NAN],
        vec![-1., 1.],
        vec![0.; 1026],
    ] {
        assert!(matches!(
            DomainPlan::prepare_with_knots(&body, &knots, DomainOptions::default(), None),
            Err(Error::Invalid)
        ));
    }
    // These are additional actual requested subdivision stations, not source
    // replacements. The same sealed graph must refuse insufficient storage.
    let knots: Vec<_> = (0..=1024).map(|i| i as f64 / 1024.).collect();
    assert!(matches!(
        DomainPlan::prepare_with_knots(
            &body,
            &knots,
            DomainOptions {
                maximum_stations: 32,
                ..Default::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    let plan = CoupledPlan::prepare(&body, CurveOptions::default(), None).unwrap();
    for bytes in [0, 64 * 1024 * 1024 + 1] {
        assert!(matches!(
            topology::CoupledTopology::prepare(&plan, 16_000_000, bytes, None),
            Err(Error::Invalid)
        ));
    }
}
