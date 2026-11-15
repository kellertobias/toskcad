use super::*;
fn with_head(f: &mut dyn FnMut(&ClosedRepairedHead<'_, '_, '_, '_, '_>)) {
    super::super::super::degree11::tests::with_head(f);
}
#[test]
fn actual_m2_original_degree11_rim_contacts_keep_source_ownership_and_costs() {
    with_head(&mut |head| {
        let carrier = head
            .degree11_thread_carrier(Default::default(), None)
            .unwrap();
        let contacts = head
            .degree11_thread_contacts(&carrier, Default::default(), None)
            .unwrap();
        assert!(std::ptr::eq(contacts.carrier(), &carrier));
        assert_eq!(
            contacts.head_circle_radius_mm(),
            head.virtual_bottom_station()[0]
        );
        assert!(!contacts.sign_domains().is_empty());
        for c in contacts.contacts() {
            assert!(
                carrier
                    .faces()
                    .iter()
                    .any(|f| f.transfer().face_entity == c.face_entity
                        && f.transfer().surface_entity == c.surface_entity)
            );
            assert!(
                c.isolating_parameter[0]
                    .cmp(c.isolating_parameter[1])
                    .is_le()
            );
        }
        for d in contacts.sign_domains() {
            assert!(
                carrier
                    .faces()
                    .iter()
                    .any(|f| f.transfer().face_entity == d.face_entity
                        && f.transfer().surface_entity == d.surface_entity)
            );
        }
        assert!(contacts.new_work_used() <= WORK);
        assert_eq!(
            contacts.total_work_used(),
            carrier.total_work_used() + contacts.new_work_used()
        );
        assert!(contacts.retained_budget_bytes() < BYTES);
        eprintln!(
            "D11 CONTACT roots{} domains{} newwork{} total{} bytes{}",
            contacts.contacts().len(),
            contacts.sign_domains().len(),
            contacts.new_work_used(),
            contacts.total_work_used(),
            contacts.retained_budget_bytes()
        );
    });
}
#[test]
fn original_degree11_contacts_honor_new_and_aggregate_limits_and_cancel() {
    with_head(&mut |head| {
        let carrier = head
            .degree11_thread_carrier(Default::default(), None)
            .unwrap();
        let options = Degree11ContactOptions::default();
        let baseline = head
            .degree11_thread_contacts(&carrier, options, None)
            .unwrap();
        let used = baseline.new_work_used();
        let total = baseline.total_work_used();
        for options in [
            Degree11ContactOptions {
                maximum_work: used - 1,
                ..options
            },
            Degree11ContactOptions {
                maximum_total_work: total - 1,
                ..options
            },
            Degree11ContactOptions {
                maximum_total_work: carrier.total_work_used() - 1,
                ..options
            },
            Degree11ContactOptions {
                maximum_bytes: 0,
                ..options
            },
        ] {
            assert!(matches!(
                head.degree11_thread_contacts(&carrier, options, None),
                Err(PhillipsReadError::Limit)
            ));
        }
        let exact = head
            .degree11_thread_contacts(
                &carrier,
                Degree11ContactOptions {
                    maximum_work: used,
                    maximum_total_work: total,
                    maximum_contacts: baseline.contacts().len(),
                    ..options
                },
                None,
            )
            .unwrap();
        assert_eq!(exact.new_work_used(), used);
        assert_eq!(exact.total_work_used(), total);
        for options in [
            Degree11ContactOptions {
                maximum_work: WORK + 1,
                ..options
            },
            Degree11ContactOptions {
                maximum_total_work: TOTAL_WORK + 1,
                ..options
            },
            Degree11ContactOptions {
                maximum_bytes: BYTES + 1,
                ..options
            },
            Degree11ContactOptions {
                maximum_contacts: ROOTS + 1,
                ..options
            },
        ] {
            assert!(matches!(
                head.degree11_thread_contacts(&carrier, options, None),
                Err(PhillipsReadError::Invalid)
            ));
        }
        let n = Cell::new(0usize);
        let poll = || {
            n.set(n.get() + 1);
            false
        };
        head.degree11_thread_contacts(&carrier, options, Some(&poll))
            .unwrap();
        for at in [1, 2, 10, n.get() / 2, n.get()] {
            let calls = Cell::new(0usize);
            let cancel = || {
                calls.set(calls.get() + 1);
                calls.get() == at
            };
            assert!(matches!(
                head.degree11_thread_contacts(&carrier, options, Some(&cancel)),
                Err(PhillipsReadError::Cancelled)
            ));
        }
    });
}
#[test]
fn independent_degree22_power_and_factorial_oracle_retains_nondyadic_root() {
    let mut work = Work {
        remaining: WORK,
        interrupt: None,
        cancelled: Cell::new(false),
    };
    let mut a = Arithmetic { work: &mut work };
    // x=u^11, y=0, R=1/2 gives EXACT primitive4u^22-1.
    let poles = std::array::from_fn(|i| [if i == 11 { 1. } else { 0. }, 0., 0.]);
    let p = a.radius_polynomial11(poles, 0.5).unwrap();
    assert_eq!(p.degree(), 22);
    let coefficients = a.bernstein22(&p).unwrap();
    let factorial = (1_u64..=22).map(Z::from).product::<Z>();
    assert!(coefficients[..22].iter().all(|x| *x == -&factorial));
    assert_eq!(coefficients[22], 3 * &factorial);
    let certificate = roots::certify22(p, &mut a, &mut roots::Budget::new()).unwrap();
    assert_eq!(certificate.roots.len(), 1);
    assert_eq!(
        certificate
            .domains
            .iter()
            .map(|d| d.sign)
            .collect::<Vec<_>>(),
        vec![-1, 1]
    );
    // Independent rational sign oracle; no evaluation by Arithmetic or source
    // samples can stand in for the exact algebraic isolator.
    let evaluate = |t: DyadicParameter| {
        let denominator = Z::one() << t.denominator_exponent();
        Z::from(4) * Z::from(t.numerator()).pow(22) - denominator.pow(22)
    };
    let bracket = certificate.roots[0].bracket;
    assert!(!evaluate(bracket[0]).is_positive());
    assert!(!evaluate(bracket[1]).is_negative());
}
#[test]
fn degree22_constant_coincidence_and_real_double_contacts_are_not_sign_welded() {
    let mut work = Work {
        remaining: WORK,
        interrupt: None,
        cancelled: Cell::new(false),
    };
    let mut a = Arithmetic { work: &mut work };
    for offset in [-11., -3.] {
        // Exact integer Bernstein controls give x=22u+offset: one dyadic
        // and one nondyadic true double root at -offset/22.
        let poles = std::array::from_fn(|i| [2. * i as f64 + offset, 0.5, 0.]);
        let p = a.radius_polynomial11(poles, 0.5).unwrap();
        let c = roots::certify22(p, &mut a, &mut roots::Budget::new()).unwrap();
        assert_eq!(c.roots.len(), 1);
        assert!(c.domains.iter().all(|d| d.sign > 0));
        let bracket = c.roots[0].bracket;
        // Exact root rational: (-offset)/22, not a sampled trim.
        for (side, t) in bracket.into_iter().enumerate() {
            let value = Z::from(22) * Z::from(t.numerator())
                + Z::from(offset as i64) * (Z::one() << t.denominator_exponent());
            assert!(if side == 0 {
                !value.is_positive()
            } else {
                !value.is_negative()
            });
        }
    }
    let constant = a.radius_polynomial11([[0.; 3]; 12], 1.).unwrap();
    let c = roots::certify22(constant, &mut a, &mut roots::Budget::new()).unwrap();
    assert!(c.roots.is_empty());
    assert_eq!(c.domains[0].sign, -1);
    let coincident = a.radius_polynomial11([[1., 0., 0.]; 12], 1.).unwrap();
    assert!(matches!(
        roots::certify22(coincident, &mut a, &mut roots::Budget::new()),
        Err(ReadError::Unsupported)
    ));
}

#[test]
fn actual_nominal_short_h0_degree11_event_inventory_uses_original_body() {
    for (left, posed) in [(false, false), (true, true)] {
        let text = super::super::tests::actual_h0_text(left, posed);
        super::super::super::super::common_mouth::tests::derived_thread_probe::with_head(
            &text,
            &mut |head| {
                let carrier = head
                    .degree11_thread_carrier(Default::default(), None)
                    .unwrap();
                let contacts = head
                    .degree11_thread_contacts(&carrier, Default::default(), None)
                    .unwrap();
                eprintln!(
                    "H0 D11 left{left} posed{posed} roots{} domains{} radius{} new{}",
                    contacts.contacts().len(),
                    contacts.sign_domains().len(),
                    contacts.head_circle_radius_mm(),
                    contacts.new_work_used()
                );
                assert!(std::ptr::eq(contacts.carrier().body(), head.body()));
                assert!(!contacts.sign_domains().is_empty());
                if !contacts.contacts().is_empty() {
                    assert!(matches!(
                        head.degree11_thread_contacts(
                            &carrier,
                            Degree11ContactOptions {
                                maximum_contacts: 0,
                                ..Default::default()
                            },
                            None
                        ),
                        Err(PhillipsReadError::Limit)
                    ));
                }
                for c in contacts.contacts() {
                    assert!(
                        carrier
                            .faces()
                            .iter()
                            .any(|f| f.transfer().surface_entity == c.surface_entity)
                    );
                }
            },
        );
    }
}

#[test]
fn identical_coordinates_in_another_body_do_not_transfer_contact_ownership() {
    with_head(&mut |head| {
        let carrier = head
            .degree11_thread_carrier(Default::default(), None)
            .unwrap();
        with_head(&mut |other| {
            assert!(!std::ptr::eq(head.body(), other.body()));
            assert!(matches!(
                other.degree11_thread_contacts(&carrier, Default::default(), None),
                Err(PhillipsReadError::Invalid)
            ));
        });
    });
}

#[test]
fn actual_pythagorean_contact_uses_shared_original_vertex_and_both_adjacent_domains() {
    use num_rational::BigRational as Q;
    let text = super::super::tests::offnominal_text(5., 3., 4.);
    super::super::super::super::common_mouth::tests::derived_thread_probe::with_head(
        &text,
        &mut |head| {
            let carrier = head
                .degree11_thread_carrier(Default::default(), None)
                .unwrap();
            let contacts = head
                .degree11_thread_contacts(&carrier, Default::default(), None)
                .unwrap();
            let knots = contacts
                .contacts()
                .iter()
                .enumerate()
                .filter(|(_, c)| matches!(c.source_identity, Degree11RimContactId::Knot(_)))
                .collect::<Vec<_>>();
            assert!(!knots.is_empty());
            for (index, c) in knots {
                let f = carrier
                    .faces()
                    .iter()
                    .find(|f| f.transfer().surface_entity == c.surface_entity)
                    .unwrap();
                let [_, height] = f.source().dimensions();
                let u = 11
                    * (c.original_u_cell
                        + usize::from(c.isolating_parameter == [DyadicParameter::ONE; 2]));
                let p = f.source().poles()[u * height];
                let [x, y, _] = p.map(|p| Q::from_float(p).unwrap());
                assert_eq!(&x * &x + &y * &y, Q::from_integer(Z::from(25)));
                let references = contacts
                    .sign_domains()
                    .iter()
                    .flat_map(|d| d.ends)
                    .filter(|end| matches!(end,Degree11DomainEnd::Contact(i) if *i==index))
                    .count();
                assert_eq!(references, 2);
            }
            for (rim, f) in carrier.faces().iter().enumerate() {
                for k in 1..f.source().spans()[0] {
                    assert_eq!(
                        knot11(head, &carrier, rim, k),
                        Degree11RimKnot::SurfaceKnot {
                            surface_entity: f.transfer().surface_entity,
                            original_u_knot: k
                        }
                    );
                }
            }
            assert!(matches!(
                head.degree11_thread_contacts(
                    &carrier,
                    Degree11ContactOptions {
                        maximum_contacts: 0,
                        ..Default::default()
                    },
                    None
                ),
                Err(PhillipsReadError::Limit)
            ));
            eprintln!(
                "PYTHAGORAS D11 roots{} domains{} new{}",
                contacts.contacts().len(),
                contacts.sign_domains().len(),
                contacts.new_work_used()
            );
        },
    );
}
