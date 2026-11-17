use super::*;
use crate::step_read::npt_stock::full_head::source_tests::with_source;
use num_rational::BigRational as Q;
use num_traits::Zero;
fn q(v: f64) -> Q {
    Q::from_float(v).unwrap()
}
fn contains(b: [f64; 2], v: &Q) {
    assert!(q(b[0]) <= *v && *v <= q(b[1]));
}
fn file(quarter: bool, pose: usize, length: f64) -> String {
    let text = super::super::tests::mapped(quarter, pose);
    if pose == 0 {
        return text;
    }
    let doc = CompoundDocument::parse(&text, Default::default(), None).unwrap();
    let op = doc.represented_members()[0]
        .mapping_link()
        .unwrap()
        .operator_entity();
    let parts = &doc.records()[&op];
    let base = parts
        .iter()
        .find(|e| e.name == "CARTESIAN_TRANSFORMATION_OPERATOR")
        .unwrap();
    let Value::Reference(origin) = base.args[2] else {
        panic!()
    };
    let combined = parts.iter().map(|e| e.source.unwrap()).collect::<String>();
    let operator = combined.replacen(&format!("#{origin},2.)"), &format!("#{origin},1.)"), 1);
    assert_ne!(operator, combined);
    fn change(input: &str, id: u32, body: &str) -> String {
        let start = input.find(&format!("#{id}=")).unwrap();
        let end = start + input[start..].find(';').unwrap() + 1;
        format!("{}#{id}={body};{}", &input[..start], &input[end..])
    }
    let text = change(&text, op, &format!("({operator})"));
    let xy = if pose == 1 { "0.,0." } else { "2.,-3." };
    let literal = crate::step::pipe_hardware::numbers::literal(length).unwrap();
    change(
        &text,
        origin,
        &format!("CARTESIAN_POINT('',({xy},{literal}))"),
    )
}
fn original(pose: usize, length: f64, near: bool) -> LiteralOccurrencePlacement {
    let (matrix, translation) = if pose == 0 {
        ([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]], [0.; 3])
    } else {
        (
            [[1., 0., 0.], [0., -1., 0.], [0., 0., -1.]],
            [0., 0., length],
        )
    };
    let p = LiteralOccurrencePlacement::new(matrix, translation, Default::default(), None).unwrap();
    if pose < 2 {
        p
    } else {
        p.transformed(
            [
                [0.6 + if near { 1e-10 } else { 0. }, -0.8, 0.],
                [0.8, 0.6, 0.],
                [0., 0., 1.],
            ],
            [2., -3., 0.],
            Default::default(),
            None,
        )
        .unwrap()
    }
}
fn intent() -> HexHeadIntent {
    HexHeadIntent {
        across_flats_mm: 18.,
        height_mm: 3.,
    }
}
#[test]
fn actual_mapped_source_preserves_original_law_head_and_intrinsic_plus_requested_factors() {
    for quarter in [false, true] {
        with_source(quarter, |plan| {
            let length = plan.curves().source().source().spec().length_mm;
            for (pose, near) in [(0, false), (1, false), (2, false), (2, true)] {
                let input: Arc<str> = Arc::from(file(quarter, pose, length));
                let owned = MappedActualHead::read(input, Default::default(), None).unwrap();
                let placement = original(pose, length, near);
                let qualified = owned
                    .qualify_original_source(plan, intent(), &placement, Default::default(), None)
                    .unwrap();
                assert!(std::ptr::eq(qualified.actual(), &owned));
                assert!(std::ptr::eq(qualified.source_plan(), plan));
                assert!(std::ptr::eq(qualified.original_placement(), &placement));
                let r = qualified.receipt();
                assert!(r.world_sewn_reference_position_error_mm <= 1e-7);
                assert_eq!(
                    r.thread_axial_transfer_displacement_mm,
                    owned
                        .local_geometry()
                        .axial_sewing_transfer()
                        .domain_displacement_mm()
                );
                assert_eq!(
                    r.aggregate_work,
                    owned.receipt().aggregate_work
                        + plan.work_used()
                        + placement.work_used()
                        + r.head_and_frame_work
                        + r.fresh_thread_correspondence_work
                );
                assert_eq!(
                    r.fresh_thread_correspondence_work,
                    r.original_raw_thread.query_work_used()
                );
                // Fresh old public source proof independently agrees with the
                // reused immutable material's actual raw source correspondence.
                let local = qualify_actual_file(
                    owned.local_geometry().actual_file(),
                    plan,
                    intent(),
                    Default::default(),
                    None,
                )
                .unwrap();
                assert_eq!(
                    local.receipt().original_thread.position_error_mm(),
                    r.original_raw_thread.position_error_mm()
                );
                assert_eq!(
                    local
                        .receipt()
                        .original_thread
                        .original_derivative_error_mm(),
                    r.original_raw_thread.original_derivative_error_mm()
                );
                assert_eq!(
                    local.receipt().head_author_corner_error_mm,
                    r.head_author_corner_error_mm
                );
                let bounds = owned.local_geometry().local_sewn_enclosure_mm();
                let frame = owned.actual_frame();
                let fifth = Q::from_integer(1.into()) / Q::from_integer(5.into());
                let columns = match pose {
                    0 => [
                        [q(1.), q(0.), q(0.)],
                        [q(0.), q(1.), q(0.)],
                        [q(0.), q(0.), q(1.)],
                    ],
                    1 => [
                        [q(1.), q(0.), q(0.)],
                        [q(0.), q(-1.), q(0.)],
                        [q(0.), q(0.), q(-1.)],
                    ],
                    _ => [
                        [&fifth * q(3.), &fifth * q(4.), q(0.)],
                        [&fifth * q(4.), &fifth * q(-3.), q(0.)],
                        [q(0.), q(0.), q(-1.)],
                    ],
                };
                let shift = match pose {
                    0 => [0., 0., 0.],
                    1 => [0., 0., length],
                    _ => [2., -3., length],
                };
                for mask in 0..8 {
                    let point = std::array::from_fn::<_, 3, _>(|j| q(bounds[j][(mask >> j) & 1]));
                    let mut authored = point.clone();
                    // Independent exact rational ORIGINAL factors, evaluated in
                    // append order; no interval-composed columns are the oracle.
                    for factor in placement.transforms() {
                        authored = std::array::from_fn(|i| {
                            (0..3).fold(q(factor.translation_mm()[i]), |s, j| {
                                s + q(factor.matrix()[i][j]) * &authored[j]
                            })
                        });
                    }
                    let imported = std::array::from_fn::<_, 3, _>(|i| {
                        (0..3).fold(q(shift[i]), |s, j| s + &columns[j][i] * &point[j])
                    });
                    let mut delta = Q::zero();
                    for i in 0..3 {
                        let difference = &authored[i] - &imported[i];
                        delta += if difference < Q::zero() {
                            -difference
                        } else {
                            difference
                        };
                        contains(frame.origin_mm()[i], &q(shift[i]));
                    }
                    assert!(delta <= q(r.factor_chain_discrepancy_mm));
                }
                eprintln!("MAPPED_SOURCE quarter={quarter} pose={pose} near={near} {r:?}");
            }
        });
    }
}
#[test]
fn actual_mapped_source_refuses_wrong_head_wrong_datum_wrong_law_and_stricter_unearned_policy() {
    with_source(false, |plan| {
        let length = plan.curves().source().source().spec().length_mm;
        let owned =
            MappedActualHead::read(Arc::from(file(false, 2, length)), Default::default(), None)
                .unwrap();
        let placement = original(2, length, false);
        assert!(matches!(
            owned.qualify_original_source(
                plan,
                HexHeadIntent {
                    across_flats_mm: 20.,
                    height_mm: 3.
                },
                &placement,
                Default::default(),
                None
            ),
            Err(Error::Unsupported)
        ));
        let bad = placement
            .transformed(
                [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
                [1., 0., 0.],
                Default::default(),
                None,
            )
            .unwrap();
        assert!(matches!(
            owned.qualify_original_source(plan, intent(), &bad, Default::default(), None),
            Err(Error::Unsupported)
        ));
        with_source(true, |wrong| {
            assert!(matches!(
                owned.qualify_original_source(
                    wrong,
                    intent(),
                    &placement,
                    Default::default(),
                    None
                ),
                Err(Error::Unsupported)
            ));
        });
        let loose = MappedActualHead::read(
            Arc::from(file(false, 2, length)),
            OwnedHeadOptions {
                material_tolerance_mm: 1e-6,
                ..Default::default()
            },
            None,
        )
        .unwrap();
        assert!(matches!(
            loose.qualify_original_source(plan, intent(), &placement, Default::default(), None),
            Err(Error::Unsupported)
        ));
        loose
            .qualify_original_source(
                plan,
                intent(),
                &placement,
                SourceOptions {
                    material_tolerance_mm: 1e-6,
                    ..Default::default()
                },
                None,
            )
            .unwrap();
        owned
            .qualify_original_source(plan, intent(), &placement, Default::default(), None)
            .unwrap();
    });
}
#[test]
fn actual_mapped_source_work_storage_cancel_and_accuracy_frontiers_are_atomic() {
    use std::cell::Cell;
    with_source(false, |plan| {
        let length = plan.curves().source().source().spec().length_mm;
        let owned =
            MappedActualHead::read(Arc::from(file(false, 2, length)), Default::default(), None)
                .unwrap();
        let placement = original(2, length, false);
        let r = owned
            .qualify_original_source(plan, intent(), &placement, Default::default(), None)
            .unwrap()
            .receipt();
        for (work, bytes, ok) in [
            (r.aggregate_work - 1, 64 * 1024 * 1024, false),
            (r.aggregate_work, 64 * 1024 * 1024, true),
            (272_000_000, r.accounted_peak_bytes - 1, false),
            (272_000_000, r.accounted_peak_bytes, true),
        ] {
            let value = owned.qualify_original_source(
                plan,
                intent(),
                &placement,
                SourceOptions {
                    resources: Options {
                        maximum_work: work,
                        maximum_bytes: bytes,
                        ..Default::default()
                    },
                    ..Default::default()
                },
                None,
            );
            assert_eq!(value.is_ok(), ok);
            if !ok {
                assert!(matches!(value, Err(Error::Limit)));
            }
        }
        let narrow = SourceOptions {
            position_error_mm: r.world_sewn_reference_position_error_mm / 2.,
            ..Default::default()
        };
        assert!(matches!(
            owned.qualify_original_source(plan, intent(), &placement, narrow, None),
            Err(Error::Unsupported)
        ));
        let count = Cell::new(0);
        let observe = || {
            count.set(count.get() + 1);
            false
        };
        owned
            .qualify_original_source(
                plan,
                intent(),
                &placement,
                Default::default(),
                Some(&observe),
            )
            .unwrap();
        for at in [1, count.get() / 2, count.get()] {
            let seen = Cell::new(0);
            let cancel = || {
                seen.set(seen.get() + 1);
                seen.get() >= at
            };
            assert!(matches!(
                owned.qualify_original_source(
                    plan,
                    intent(),
                    &placement,
                    Default::default(),
                    Some(&cancel)
                ),
                Err(Error::Cancelled)
            ));
        }
        let extra = owned
            .qualify_original_source(
                plan,
                intent(),
                &placement,
                SourceOptions {
                    resources: Options {
                        additional_resident_bytes: 4096,
                        ..Default::default()
                    },
                    ..Default::default()
                },
                None,
            )
            .unwrap();
        assert_eq!(
            extra.receipt().accounted_peak_bytes,
            r.accounted_peak_bytes + 4096
        );
        assert_eq!(extra.receipt().aggregate_work, r.aggregate_work);
        // A source/frame may already have been declared live while decoding.
        // Replace the prior external construction scope; count those SAME
        // owners once rather than adding their payload to itself.
        let live = plan.retained_budget_bytes() + placement.retained_budget_bytes();
        let predeclared = MappedActualHead::read(
            Arc::from(file(false, 2, length)),
            OwnedHeadOptions {
                resources: Options {
                    additional_resident_bytes: live,
                    ..Default::default()
                },
                ..Default::default()
            },
            None,
        )
        .unwrap();
        let reused = predeclared
            .qualify_original_source(plan, intent(), &placement, Default::default(), None)
            .unwrap();
        assert_eq!(
            reused.receipt().accounted_peak_bytes,
            r.accounted_peak_bytes
        );
        assert_eq!(reused.receipt().aggregate_work, r.aggregate_work);
        let retry = owned
            .qualify_original_source(plan, intent(), &placement, Default::default(), None)
            .unwrap();
        assert_eq!(
            retry.receipt().world_sewn_reference_position_error_mm,
            r.world_sewn_reference_position_error_mm
        );
    });
}
