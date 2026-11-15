use super::*;
use crate::authored_hanger_bolt::{BoundaryOptions, Hanger};
use num_bigint::BigInt;
use std::cell::Cell;
fn exact_text(s: &str) -> BigRational {
    let (mantissa, exponent) = s.split_once("E-").unwrap_or((s, "0"));
    BigRational::new(
        BigInt::parse_bytes(mantissa.strip_suffix('.').unwrap().as_bytes(), 10).unwrap(),
        BigInt::from(10_u32).pow(exponent.parse().unwrap()),
    )
}
fn kinds(field: &EncodedField<'_>) -> Vec<SideSupportKind> {
    let mut out = vec![
        SideSupportKind::Meridian {
            turn_end: End::Start,
        },
        SideSupportKind::Meridian { turn_end: End::End },
    ];
    for angular_piece in 0..field.angular.coefficients().len() {
        for phase_end in [End::Start, End::End] {
            out.push(SideSupportKind::Helix {
                phase_end,
                angular_piece,
            });
        }
    }
    out
}
#[test]
fn exact_side_controls_and_finite_decimals_reconstruct_the_original_tensor_slices() {
    for (left, nonbinary) in [(false, false), (true, true)] {
        let source = Hanger::new(
            crate::authored_hanger_bolt::tests::fixture(left, nonbinary),
            None,
        )
        .unwrap();
        let boundary = source.boundary(BoundaryOptions::default(), None).unwrap();
        let field = EncodedField::prepare(&boundary, CurveOptions::default(), None).unwrap();
        let mut total = field.work_used();
        let mut peak_query = 0;
        let mut count = 0;
        for meridian in 0..field.meridians.len() {
            let cell = field
                .cells
                .iter()
                .position(|c| c.meridian == meridian)
                .unwrap();
            for kind in kinds(&field) {
                let support = field
                    .side_support(cell, kind, DomainOptions::default(), None)
                    .unwrap();
                assert!(std::ptr::eq(support.field(), &field));
                assert_eq!(support.source_cell(), cell);
                assert_eq!(support.kind(), kind);
                assert!(support.retained_bytes() > size_of::<SideSupport<'_, '_>>());
                assert!(support.position_error_mm() <= CurveOptions::default().position_error_mm);
                assert!(
                    support.derivative_error_mm() <= CurveOptions::default().derivative_error_mm
                );
                total += support.work_used();
                peak_query = peak_query.max(support.work_used());
                count += 1;
                let (fixed, piece, column, parameter) = match kind {
                    SideSupportKind::Meridian {
                        turn_end: End::Start,
                    } => (None, 0, 0, [0., 1.]),
                    SideSupportKind::Meridian { turn_end: End::End } => {
                        (None, field.angular.coefficients().len() - 1, 10, [0., 1.])
                    }
                    SideSupportKind::Helix {
                        phase_end,
                        angular_piece,
                    } => (
                        Some(if phase_end == End::Start {
                            0
                        } else {
                            field.meridians[meridian].degree()
                        }),
                        angular_piece,
                        0,
                        [
                            field.angular.knots()[angular_piece],
                            field.angular.knots()[angular_piece + 1],
                        ],
                    ),
                };
                assert_eq!(support.parameter(), parameter);
                assert_eq!(
                    support.degree(),
                    if fixed.is_some() {
                        10
                    } else {
                        field.meridians[meridian].degree()
                    }
                );
                let mut oracle = Work::new(16_000_000, None).unwrap();
                for (i, pole) in support.poles().iter().enumerate() {
                    let (u, v) = if let Some(u) = fixed {
                        (u, i)
                    } else {
                        (i, column)
                    };
                    let original = field
                        .control(cell, u, piece, v, DomainOptions::default(), None)
                        .unwrap();
                    let weight = pole.weight().number(&mut oracle).unwrap();
                    assert_eq!(exact_text(pole.step_weight().as_str()), weight.0);
                    assert_eq!(
                        weight.0,
                        original.homogeneous[3].number(&mut oracle).unwrap().0
                    );
                    for axis in 0..3 {
                        let coordinate = pole.point()[axis].number(&mut oracle).unwrap();
                        assert_eq!(exact_text(pole.step_point()[axis].as_str()), coordinate.0);
                        assert_eq!(
                            coordinate.mul(&weight, &mut oracle).unwrap().0,
                            original.homogeneous[axis].number(&mut oracle).unwrap().0
                        );
                    }
                }
            }
        }
        eprintln!(
            "exact encoded sides left={left} nonbinary={nonbinary} count={count} standalone_sum_plus_field={total} max_query={peak_query}"
        );
    }
}
#[test]
fn actual_side_support_limits_cancellation_and_invalid_owners_are_atomic() {
    let source = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, false),
        None,
    )
    .unwrap();
    let boundary = source.boundary(BoundaryOptions::default(), None).unwrap();
    let field = EncodedField::prepare(&boundary, CurveOptions::default(), None).unwrap();
    let kind = SideSupportKind::Helix {
        phase_end: End::Start,
        angular_piece: 0,
    };
    for (cell, kind) in [
        (usize::MAX, kind),
        (
            0,
            SideSupportKind::Helix {
                phase_end: End::End,
                angular_piece: usize::MAX,
            },
        ),
    ] {
        assert_eq!(
            field
                .side_support(cell, kind, DomainOptions::default(), None)
                .err(),
            Some(Error::Invalid)
        );
    }
    for maximum_work in [0, 1, 16] {
        assert_eq!(
            field
                .side_support(
                    0,
                    kind,
                    DomainOptions {
                        maximum_work,
                        ..Default::default()
                    },
                    None
                )
                .err(),
            Some(if maximum_work == 0 {
                Error::Invalid
            } else {
                Error::Limit
            })
        );
    }
    for maximum_bytes in [
        0,
        field.retained_bytes(),
        field.retained_bytes()
            + 16 * 1024
            + super::super::super::SCRATCH
            + size_of::<SideSupport<'_, '_>>()
            - 1,
    ] {
        assert_eq!(
            field
                .side_support(
                    0,
                    kind,
                    DomainOptions {
                        maximum_bytes,
                        ..Default::default()
                    },
                    None
                )
                .err(),
            Some(if maximum_bytes == 0 {
                Error::Invalid
            } else {
                Error::Limit
            })
        );
    }
    let calls = Cell::new(0);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    let original = field
        .side_support(0, kind, DomainOptions::default(), Some(&count))
        .unwrap();
    assert!(calls.get() > 64);
    for stop in [0, 1, 8, 64, calls.get() - 2, calls.get() - 1] {
        let n = Cell::new(0);
        let cancel = || {
            let value = n.get();
            n.set(value + 1);
            value >= stop
        };
        assert_eq!(
            field
                .side_support(0, kind, DomainOptions::default(), Some(&cancel))
                .err(),
            Some(Error::Cancelled)
        );
    }
    let repeated = field
        .side_support(0, kind, DomainOptions::default(), None)
        .unwrap();
    assert_eq!(original.poles().len(), repeated.poles().len());
    for (a, b) in original.poles().iter().zip(repeated.poles()) {
        assert_eq!(a.point(), b.point());
        assert_eq!(a.weight(), b.weight());
    }
}
#[test]
fn actual_shared_batch_retains_every_source_owner_and_same_decimal_geometry() {
    for (left, nb) in [(false, false), (true, true)] {
        let source =
            Hanger::new(crate::authored_hanger_bolt::tests::fixture(left, nb), None).unwrap();
        let boundary = source.boundary(BoundaryOptions::default(), None).unwrap();
        let field = EncodedField::prepare(&boundary, CurveOptions::default(), None).unwrap();
        let mut requests = Vec::new();
        for meridian in 0..field.meridians.len() {
            let cell = field
                .cells
                .iter()
                .position(|c| c.meridian == meridian)
                .unwrap();
            for kind in kinds(&field) {
                requests.push((cell, kind));
            }
        }
        let batch = field
            .side_support_batch(&requests, DomainOptions::default(), None)
            .unwrap();
        assert_eq!(batch.supports().len(), 180);
        assert!(batch.retained_bytes() > size_of::<SideSupportBatch<'_, '_>>());
        assert!(batch.work_used() + field.work_used() < 16_000_000);
        for (support, &(cell, kind)) in batch.supports().iter().zip(&requests) {
            assert!(std::ptr::eq(support.field(), &field));
            assert_eq!(support.source_cell(), cell);
            assert_eq!(support.kind(), kind);
            for pole in support.poles() {
                for axis in 0..3 {
                    let mut w = Work::new(16_000_000, None).unwrap();
                    assert_eq!(
                        exact_text(pole.step_point()[axis].as_str()),
                        pole.point()[axis].number(&mut w).unwrap().0
                    );
                }
            }
        }
        let work = batch.work_used();
        drop(batch);
        assert!(
            field
                .side_support_batch(
                    &requests,
                    DomainOptions {
                        maximum_work: work,
                        ..Default::default()
                    },
                    None
                )
                .is_ok()
        );
        assert!(matches!(
            field.side_support_batch(
                &requests,
                DomainOptions {
                    maximum_work: work - 1,
                    ..Default::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
        assert!(matches!(
            field.side_support_batch(
                &requests,
                DomainOptions {
                    maximum_bytes: field.retained_bytes(),
                    ..Default::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
        assert!(matches!(
            field.side_support_batch(&requests, DomainOptions::default(), Some(&|| true)),
            Err(Error::Cancelled)
        ));
        assert!(matches!(
            field.side_support_batch(
                &[(usize::MAX, kinds(&field)[0])],
                DomainOptions::default(),
                None
            ),
            Err(Error::Invalid)
        ));
        let empty = field
            .side_support_batch(&[], DomainOptions::default(), None)
            .unwrap();
        assert!(empty.supports().is_empty());
        eprintln!("source-owned exact batch left={left} nb={nb}work={work}");
    }
}
