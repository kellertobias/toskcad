use super::exact::*;
use super::roots::DyadicParameter as Dyad;
use super::roots::*;
use super::*;
use num_bigint::BigInt as Z;
fn actual_h0_text(left: bool, posed: bool) -> String {
    use crate::{
        head_profile::Profile, helical_thread::Thread, phillips_drive::PhillipsScrew,
        profiled_screw::SourceSpec,
    };
    let original = super::super::super::super::tests::fixture(0, left, true, 1);
    let mut spec = original.thread().spec();
    spec.length_mm = 0.8;
    spec.depth_mm = 5. * 3_f64.sqrt() * spec.pitch_mm / 16.;
    spec.crest_width_mm = spec.pitch_mm / 8.;
    let head = original.spec().head.unwrap();
    let source = PhillipsScrew::new(
        Thread::new(spec).unwrap(),
        Profile::new_for_thread_union(head, 1., 5.2).unwrap(),
        SourceSpec {
            head: Some(head),
            shank_length_mm: 5.2,
            minimum_wall_mm: 1e-5,
        },
        original.cavity().clone(),
        100_000,
        &|| false,
    )
    .unwrap();
    let (o, u, v) = if posed {
        ([10., 20., 30.], [0.8, 0.6, 0.], [-0.36, 0.48, 0.8])
    } else {
        ([0.; 3], [1., 0., 0.], [0., 1., 0.])
    };
    crate::step::write_phillips_screw(&source, o, u, v, 1e-7).unwrap()
}
#[test]
fn actual_h0_nominal_packet_contact_cost_probe() {
    for (left, posed) in [(false, false), (true, false), (false, true), (true, true)] {
        let text = actual_h0_text(left, posed);
        let start = std::time::Instant::now();
        super::super::super::common_mouth::tests::derived_thread_probe::with_head(
            &text,
            &mut |head| {
                let source = head
                    .contact_thread_carrier(Default::default(), None)
                    .unwrap();
                eprintln!(
                    "CONTACT PASS left{left} pose{posed} roots{} domains{} work{} bytes{} elapsed{:?}",
                    source.contacts().len(),
                    source.sign_domains().len(),
                    source.work_used(),
                    source.retained_budget_bytes(),
                    start.elapsed()
                );
                assert!(std::ptr::eq(source.body(), head.body()));
                assert_eq!(source.thread_carrier().source().faces.len(), 6);
            },
        );
    }
}

#[test]
fn integer_oracle_retains_double_and_endpoint_roots_without_sign_welding() {
    let mut work = Work {
        remaining: 16_000_000,
        interrupt: None,
        cancelled: Cell::new(false),
    };
    let mut arithmetic = Arithmetic { work: &mut work };
    // t(t-1)(4t-1)^2(4t-3): one double contact and two endpoint contacts.
    let mut budget = Budget::new();
    let certificate = certify(
        [0, 3, -31, 108, -144, 64]
            .into_iter()
            .map(Z::from)
            .collect(),
        &mut arithmetic,
        &mut budget,
    )
    .unwrap();
    let independent = certify_sturm(
        [0, 3, -31, 108, -144, 64]
            .into_iter()
            .map(Z::from)
            .collect(),
        &mut arithmetic,
        &mut budget,
    )
    .unwrap();
    assert_eq!(independent.roots.len(), certificate.roots.len());
    assert_eq!(certificate.roots.len(), 4);
    assert_eq!(
        certificate
            .domains
            .iter()
            .map(|d| d.sign)
            .collect::<Vec<_>>(),
        vec![1, 1, -1]
    );
    let expected = [
        Dyad::ZERO,
        Dyad {
            numerator: 1,
            depth: 2,
        },
        Dyad {
            numerator: 3,
            depth: 2,
        },
        Dyad::ONE,
    ];
    for (root, expected) in certificate.roots.iter().zip(expected) {
        assert!(root.bracket[0].cmp(expected).is_le());
        assert!(root.bracket[1].cmp(expected).is_ge());
    }
}

#[test]
fn exact_nondyadic_double_contact_and_distinct_subulp_events_remain_owned() {
    use num_traits::{One, Zero};
    let mut work = Work {
        remaining: 16_000_000,
        interrupt: None,
        cancelled: Cell::new(false),
    };
    let mut arithmetic = Arithmetic { work: &mut work };
    let mut budget = Budget::new();
    let repeated = certify(
        [1, -6, 9].into_iter().map(Z::from).collect(),
        &mut arithmetic,
        &mut budget,
    )
    .unwrap();
    assert_eq!(repeated.roots.len(), 1);
    assert_eq!(
        repeated.domains.iter().map(|d| d.sign).collect::<Vec<_>>(),
        [1, 1]
    );
    let [lo, hi] = repeated.roots[0].bracket;
    assert!(u128::from(lo.numerator) * 3 <= 1u128 << lo.depth);
    assert!(u128::from(hi.numerator) * 3 >= 1u128 << hi.depth);
    let c = Z::one() << 58;
    let one = Z::one();
    let half = &c / 2;
    let p = vec![(&half + &one) * (&half + 2), -&c * (&c + 3), &c * &c];
    let separated = certify(p, &mut arithmetic, &mut budget).unwrap();
    assert_eq!(separated.roots.len(), 2);
    assert!(
        separated.roots[0].bracket[1]
            .cmp(separated.roots[1].bracket[0])
            .is_lt()
    );
    let rounded = |r: &Root| r.bracket[0].numerator as f64 / 2_f64.powi(r.bracket[0].depth as i32);
    assert_eq!(rounded(&separated.roots[0]), rounded(&separated.roots[1]));
    assert_eq!(
        certify(vec![Z::zero()], &mut arithmetic, &mut budget)
            .err()
            .unwrap(),
        ReadError::Unsupported
    );
    for (constant, sign) in [(Z::one(), 1), (-Z::one(), -1)] {
        let cell = certify(vec![constant], &mut arithmetic, &mut budget).unwrap();
        assert!(cell.roots.is_empty());
        assert_eq!(cell.domains.len(), 1);
        assert_eq!(cell.domains[0].sign, sign);
    }
}
#[test]
fn actual_contact_receipt_public_guards_budgets_and_cancellation_are_atomic() {
    let text = actual_h0_text(false, true);
    super::super::super::common_mouth::tests::derived_thread_probe::with_head(&text, &mut |head| {
        let baseline = head
            .contact_thread_carrier(Default::default(), None)
            .unwrap();
        assert_eq!(baseline.contacts().len(), 2);
        assert!(
            baseline
                .sign_domains()
                .iter()
                .any(|d| d.role() == RimDomainRole::ThreadTop)
        );
        assert!(
            baseline
                .sign_domains()
                .iter()
                .any(|d| d.role() == RimDomainRole::HeadBottom)
        );
        for options in [
            DerivedThreadOptions {
                maximum_work: 0,
                ..Default::default()
            },
            DerivedThreadOptions {
                maximum_bytes: EXACT_SCRATCH,
                ..Default::default()
            },
            DerivedThreadOptions {
                maximum_controls: 1,
                ..Default::default()
            },
        ] {
            assert_eq!(
                head.contact_thread_carrier(
                    ThreadContactOptions {
                        source: options,
                        ..Default::default()
                    },
                    None
                )
                .err()
                .unwrap(),
                PhillipsReadError::Limit
            );
        }
        assert_eq!(
            head.contact_thread_carrier(
                ThreadContactOptions {
                    source: DerivedThreadOptions {
                        maximum_work: WORK + 1,
                        ..Default::default()
                    },
                    ..Default::default()
                },
                None
            )
            .err()
            .unwrap(),
            PhillipsReadError::Invalid
        );
        assert_eq!(
            head.contact_thread_carrier(Default::default(), Some(&|| true))
                .err()
                .unwrap(),
            PhillipsReadError::Cancelled
        );
        let polls = Cell::new(0usize);
        let observe = || {
            polls.set(polls.get() + 1);
            false
        };
        let success = head
            .contact_thread_carrier(Default::default(), Some(&observe))
            .unwrap();
        for stop in [1, polls.get() / 2, polls.get() - 2, polls.get() - 1] {
            let count = Cell::new(0usize);
            let cancel = || {
                let n = count.get();
                count.set(n + 1);
                n == stop
            };
            assert_eq!(
                head.contact_thread_carrier(Default::default(), Some(&cancel))
                    .err()
                    .unwrap(),
                PhillipsReadError::Cancelled
            );
        }
        assert_eq!(success.contacts().len(), baseline.contacts().len());
        assert_eq!(head.body().literal_input(), text);
    });
}

#[test]
fn actual_retained_rim_metric_has_independent_exact_whole_polynomial_oracle() {
    use num_rational::BigRational as Q;
    use num_traits::{One, Signed, Zero};
    let text = actual_h0_text(false, true);
    super::super::super::common_mouth::tests::derived_thread_probe::with_head(&text, &mut |head| {
        let receipt = head
            .contact_thread_carrier(Default::default(), None)
            .unwrap();
        let input = receipt.thread_carrier().source();
        assert_eq!(
            receipt.head_circle_radius_mm(),
            head.virtual_bottom_station()[0]
        );
        let radius = Q::from_float(receipt.head_circle_radius_mm()).unwrap();
        let mut work = Work {
            remaining: 16_000_000,
            interrupt: None,
            cancelled: Cell::new(false),
        };
        let mut arithmetic = Arithmetic { work: &mut work };
        for domain in receipt.sign_domains() {
            assert!(domain.source_curve() < input.curves.len());
            for end in domain.ends() {
                if let RimDomainEnd::Contact(i) = end {
                    assert!(i < receipt.contacts().len());
                }
            }
        }
        for contact in receipt.contacts() {
            let id = contact.source_identity();
            let curve = contact.source_curve();
            let span = contact.source_span();
            assert!(matches!(id,RimContactId::Interior{curve:c,span:s,..} if c==curve&&s==span));
            let interval = contact.isolating_parameter();
            assert!(interval[0].cmp(interval[1]).is_le());
            for p in interval {
                assert!(u128::from(p.numerator()) <= 1u128 << p.denominator_exponent());
            }
        }
        for &(edge, _) in &input.faces[0].outer {
            let edge = &input.edges[edge];
            let radial::Curve::Spline { spans, poles } = &input.curves[edge.curve] else {
                panic!("actual emitted rim must retain its spline");
            };
            assert_eq!(poles[0], input.vertices[edge.vertices[0]]);
            assert_eq!(*poles.last().unwrap(), input.vertices[edge.vertices[1]]);
            for span in 0..*spans {
                let controls = std::array::from_fn::<_, 10, _>(|i| poles[span * 9 + i]);
                let polynomial = arithmetic
                    .radius_polynomial(controls, receipt.head_circle_radius_mm())
                    .unwrap();
                let mut scale = None;
                // Both expressions have degree<=18. Agreement at19 distinct
                // EXACT rational nodes proves polynomial identity, not a
                // numerical sampling bound; the independent evaluator uses
                // literal retained controls and rational de Casteljau.
                for j in 0..19 {
                    let t = Q::new(Z::from(j), Z::from(32));
                    let one_minus = Q::one() - &t;
                    let mut xy = controls
                        .map(|p| [Q::from_float(p[0]).unwrap(), Q::from_float(p[1]).unwrap()]);
                    for level in (1..10).rev() {
                        for i in 0..level {
                            xy[i] = std::array::from_fn(|axis| {
                                &xy[i][axis] * &one_minus + &xy[i + 1][axis] * &t
                            });
                        }
                    }
                    let direct = &xy[0][0] * &xy[0][0] + &xy[0][1] * &xy[0][1] - &radius * &radius;
                    let power = polynomial.iter().rev().fold(Q::zero(), |value, c| {
                        value * &t + Q::from_integer(c.clone())
                    });
                    if !direct.is_zero() {
                        let k = &power / &direct;
                        assert!(k.is_positive());
                        if let Some(old) = &scale {
                            assert_eq!(&k, old);
                        } else {
                            scale = Some(k);
                        }
                    } else {
                        assert!(power.is_zero());
                    }
                }
                assert!(scale.is_some());
            }
        }
    });
}
#[test]
fn actual_event_quota_work_frontier_and_final_publication_capacity_are_qualified() {
    let text = actual_h0_text(false, true);
    super::super::super::common_mouth::tests::derived_thread_probe::with_head(&text, &mut |head| {
        for quota in [0, 1] {
            assert_eq!(
                head.contact_thread_carrier(
                    ThreadContactOptions {
                        maximum_contacts: quota,
                        ..Default::default()
                    },
                    None
                )
                .err()
                .unwrap(),
                PhillipsReadError::Limit
            );
        }
        let exact = head
            .contact_thread_carrier(
                ThreadContactOptions {
                    maximum_contacts: 2,
                    ..Default::default()
                },
                None,
            )
            .unwrap();
        assert_eq!(exact.contacts().len(), 2);
        assert_eq!(
            head.contact_thread_carrier(
                ThreadContactOptions {
                    maximum_contacts: ROOTS + 1,
                    ..Default::default()
                },
                None
            )
            .err()
            .unwrap(),
            PhillipsReadError::Invalid
        );
        let needed = exact.work_used();
        assert_eq!(
            head.contact_thread_carrier(
                ThreadContactOptions {
                    source: DerivedThreadOptions {
                        maximum_work: needed - 1,
                        ..Default::default()
                    },
                    maximum_contacts: 2
                },
                None
            )
            .err()
            .unwrap(),
            PhillipsReadError::Limit
        );
        let ready = head
            .contact_thread_carrier(
                ThreadContactOptions {
                    source: DerivedThreadOptions {
                        maximum_work: needed,
                        ..Default::default()
                    },
                    maximum_contacts: 2,
                },
                None,
            )
            .unwrap();
        assert_eq!(ready.work_used(), needed);
        for lower in [true, false] {
            let earned = head
                .contact_thread_carrier(Default::default(), None)
                .unwrap();
            let prefix = earned.work_used();
            let ContactThreadCarrier {
                carrier,
                report,
                head_radius,
                retained,
                ..
            } = earned;
            let limit = retained + EXACT_SCRATCH - usize::from(lower);
            let work = Work {
                remaining: 1,
                interrupt: None,
                cancelled: Cell::new(false),
            };
            let result = publish(carrier, report, head_radius, limit, prefix + 1, work);
            if lower {
                assert_eq!(result.err().unwrap(), PhillipsReadError::Limit);
            } else {
                assert_eq!(result.unwrap().work_used(), prefix + 1);
            }
        }
    });
}

#[test]
fn exact_numeric_stage_refuses_bits_depth_and_nodes_before_publication() {
    use num_traits::One;
    let mut work = Work {
        remaining: 16_000_000,
        interrupt: None,
        cancelled: Cell::new(false),
    };
    let mut arithmetic = Arithmetic { work: &mut work };
    let a = Z::one() << 4095;
    let b = Z::one() << 4096;
    assert_eq!(arithmetic.mul(&a, &b).unwrap_err(), ReadError::Limit);
    assert_eq!(
        arithmetic.shift(&Z::one(), 8192).unwrap_err(),
        ReadError::Limit
    );
    let mut budget = Budget { remaining_nodes: 0 };
    assert_eq!(
        certify(vec![Z::one(), -Z::one()], &mut arithmetic, &mut budget)
            .err()
            .unwrap(),
        ReadError::Limit
    );
    let c = Z::one() << 100;
    let half = &c / 2;
    let too_close = vec![(&half + 1) * (&half + 2), -&c * (&c + 3), &c * &c];
    let mut budget = Budget::new();
    assert_eq!(
        certify(too_close, &mut arithmetic, &mut budget)
            .err()
            .unwrap(),
        ReadError::Limit
    );
    assert!(arithmetic.work.remaining > 1_000_000);
}

#[test]
fn exact_axis_aligned_metric_and_constant_derivative_have_independent_oracles() {
    let mut work = Work {
        remaining: 16_000_000,
        interrupt: None,
        cancelled: Cell::new(false),
    };
    let mut arithmetic = Arithmetic { work: &mut work };
    let polynomial = arithmetic
        .radius_polynomial(std::array::from_fn(|i| [i as f64, 0., 0.]), 3.)
        .unwrap();
    // Degree9 controls i encode Gamma=(9t,0) exactly, so P is a positive
    // constant multiple of9t²-1. Exact axial zero leaves remain actual zeros.
    assert_eq!(polynomial, vec![Z::from(-1), Z::from(0), Z::from(9)]);
    assert_eq!(
        arithmetic.derivative(&vec![Z::from(7)]).unwrap(),
        vec![Z::from(0)]
    );
    let mut budget = Budget::new();
    let c = certify(polynomial, &mut arithmetic, &mut budget).unwrap();
    assert_eq!(c.roots.len(), 1);
    assert_eq!(
        c.domains.iter().map(|d| d.sign).collect::<Vec<_>>(),
        [-1, 1]
    );
}

#[test]
fn source_parameter_owner_aliases_and_actual_contact_interning_are_exact() {
    let text = actual_h0_text(false, true);
    super::super::super::common_mouth::tests::derived_thread_probe::with_head(&text, &mut |head| {
        let receipt = head
            .contact_thread_carrier(Default::default(), None)
            .unwrap();
        let source = receipt.thread_carrier().source();
        let rims = source.faces[0].outer.len();
        let radial::Curve::Spline { spans, .. } = &source.curves[0] else {
            panic!("actual source rim law");
        };
        // This is the production INTEGER/PARAMETER ownership stage, not a
        // claim that the fixture has geometric roots at these endpoints.
        let (first, start) = root_identity(0, 0, *spans, rims, 0, [Dyad::ZERO; 2]);
        let (last, end_start) =
            root_identity(rims - 1, *spans - 1, *spans, rims, 0, [Dyad::ONE; 2]);
        assert!(start && !end_start);
        assert_eq!(first, last);
        assert_eq!(first, RimContactId::Knot(RimKnot::Vertex(0)));
        let (a, _) = root_identity(0, 0, *spans, rims, 0, [Dyad::ONE; 2]);
        let (b, _) = root_identity(0, 1, *spans, rims, 0, [Dyad::ZERO; 2]);
        assert_eq!(a, b);
        assert_eq!(
            a,
            RimContactId::Knot(RimKnot::SplineKnot { curve: 0, knot: 1 })
        );
        let mut report = Report {
            contacts: Vec::new(),
            domains: Vec::new(),
        };
        let mut work = Work {
            remaining: 100,
            interrupt: None,
            cancelled: Cell::new(false),
        };
        let real = receipt.contacts()[0];
        assert_eq!(intern(&mut report, real, 1, &mut work).unwrap(), 0);
        assert_eq!(intern(&mut report, real, 1, &mut work).unwrap(), 0);
        assert_eq!(report.contacts.len(), 1);
        assert_eq!(
            intern(&mut report, receipt.contacts()[1], 1, &mut work).unwrap_err(),
            ReadError::Limit
        );
        for domain in receipt.sign_domains() {
            assert!(domain.source_span() < *spans);
        }
    });
}

fn offnominal_text(radius: f64, x: f64, y: f64) -> String {
    use crate::{
        head_profile::{Head, Profile},
        helical_thread::Thread,
        phillips_drive::PhillipsScrew,
        profiled_screw::SourceSpec,
    };
    let original = super::super::super::super::tests::fixture(0, false, true, 1);
    let mut thread = original.thread().spec();
    thread.major_radius_mm = radius;
    thread.depth_mm = 0.08;
    thread.length_mm = 0.8;
    thread.crest_width_mm = f64::atan2(y, x) / std::f64::consts::PI * 0.4;
    let head = Head::Crowned {
        radius_mm: 2. * radius,
        height_mm: 1.6 * radius,
        crown_radius_mm: 3.2 * radius,
        top_land_radius_mm: 1.25 * radius,
        neck_radius_mm: 0.1 * radius,
    };
    let source = PhillipsScrew::new(
        Thread::new(thread).unwrap(),
        Profile::new_for_thread_union(head, radius, 5.2).unwrap(),
        SourceSpec {
            head: Some(head),
            shank_length_mm: 5.2,
            minimum_wall_mm: 1e-5,
        },
        original.cavity().clone(),
        100_000,
        &|| false,
    )
    .unwrap();
    crate::step::write_phillips_screw(&source, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7).unwrap()
}
#[test]
fn actual_offnominal_exact_knot_contact_retains_shared_source_identity() {
    use num_rational::BigRational as Q;
    let text = offnominal_text(5., 3., 4.);
    super::super::super::common_mouth::tests::derived_thread_probe::with_head(&text, &mut |head| {
        let receipt = head
            .contact_thread_carrier(
                ThreadContactOptions {
                    maximum_contacts: 4,
                    ..Default::default()
                },
                None,
            )
            .unwrap();
        assert_eq!(receipt.contacts().len(), 4);
        let knots = receipt
            .contacts()
            .iter()
            .enumerate()
            .filter_map(|(i, c)| match c.source_identity() {
                RimContactId::Knot(knot) => Some((i, c, knot)),
                _ => None,
            })
            .collect::<Vec<_>>();
        assert_eq!(knots.len(), 1);
        let (index, contact, knot) = knots[0];
        let input = receipt.thread_carrier().source();
        let point = match knot {
            RimKnot::Vertex(v) => input.vertices[v],
            RimKnot::SplineKnot { curve, knot } => {
                let radial::Curve::Spline { poles, .. } = &input.curves[curve] else {
                    panic!("retained source spline")
                };
                poles[9 * knot]
            }
        };
        let [x, y, _] = point.map(|p| Q::from_float(p).unwrap());
        assert_eq!(&x * &x + &y * &y, Q::from_integer(Z::from(25)));
        assert_eq!(
            receipt
                .sign_domains()
                .iter()
                .flat_map(|d| d.ends())
                .filter(|e| *e == RimDomainEnd::Contact(index))
                .count(),
            2
        );
        assert!(
            contact.isolating_parameter() == [Dyad::ZERO; 2]
                || contact.isolating_parameter() == [Dyad::ONE; 2]
        );
        eprintln!(
            "EXACT KNOT work{} retained{} source{:?}",
            receipt.work_used(),
            receipt.retained_budget_bytes(),
            point
        );
    });
    let larger = offnominal_text(13., 5., 12.);
    super::super::super::common_mouth::tests::derived_thread_probe::with_head(
        &larger,
        &mut |head| {
            assert_eq!(
                head.contact_thread_carrier(Default::default(), None)
                    .err()
                    .unwrap(),
                PhillipsReadError::Limit
            );
        },
    );
}
