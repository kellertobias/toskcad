use super::*;
use std::cell::Cell;
fn source(quarter: bool, pose: usize) -> MappedActualHead {
    MappedActualHead::read(
        Arc::from(super::super::tests::mapped(quarter, pose)),
        Default::default(),
        None,
    )
    .unwrap()
}
fn occurrence() -> LiteralOccurrencePlacement {
    LiteralOccurrencePlacement::new(
        [[0.6, 0., 0.8], [0., 1., 0.], [-0.8, 0., 0.6]],
        [13., -9., 5.],
        Default::default(),
        None,
    )
    .unwrap()
}
#[test]
fn actual_imported_replay_keeps_every_original_local_record_and_earns_positive_scaled_b_then_a() {
    for quarter in [false, true] {
        for pose in 0..3 {
            for placed in [false, true] {
                let source = source(quarter, pose);
                let a = occurrence();
                let a = placed.then_some(&a);
                let output = encode(&source, a, Default::default(), None).unwrap();
                assert!(std::ptr::eq(output.original(), &source));
                assert!(output.occurrence().is_some() == placed);
                assert!(output.actual().receipt().exact_local_material_reused);
                assert!(Arc::ptr_eq(&output.actual.local, &source.local));
                assert_eq!(
                    output.actual.local.actual_file(),
                    source.local.actual_file()
                );
                assert!(output.receipt.maximum_world_position_error_mm < 1e-7);
                rational_operator_oracle(&output, pose, a);
                assert_eq!(
                    output.receipt.inherited_source_work,
                    source.receipt.aggregate_work
                );
                assert_eq!(output.receipt.inherited_occurrence_work, pose_work(a));
                assert_eq!(
                    output.receipt.aggregate_work,
                    source.receipt.aggregate_work
                        + pose_work(a)
                        + output.receipt.fresh_encoding_work
                        + output.receipt.fresh_output_binding_work
                        + output.receipt.operator_comparison_work
                );
                assert!(
                    output.receipt.aggregate_work <= MAX_WORK
                        && output.receipt.accounted_peak_bytes <= MAX_BYTES
                );
                // Fresh public cold read after ALL source/verification owners drop.
                let local = source.local.actual_file().to_string();
                let bytes: Arc<str> = Arc::from(output.as_str());
                let volume = output.actual.world_sewn_volume_mm3();
                let bounds = output.actual.world_sewn_enclosure_mm();
                eprintln!(
                    "IMPORTED_HEAD_REPLAY quarter={quarter} Bpose={pose} A={placed} work={} peak={} C0={}",
                    output.receipt.aggregate_work,
                    output.receipt.accounted_peak_bytes,
                    output.receipt.maximum_world_position_error_mm
                );
                drop(output);
                drop(source);
                let cold = MappedActualHead::read(
                    bytes,
                    OwnedHeadOptions {
                        resources: Options {
                            additional_resident_bytes: local.capacity(),
                            ..Default::default()
                        },
                        ..Default::default()
                    },
                    None,
                )
                .unwrap();
                assert!(!cold.receipt.exact_local_material_reused);
                assert_eq!(cold.local.actual_file(), local);
                assert_eq!(cold.world_sewn_volume_mm3(), volume);
                assert_eq!(cold.world_sewn_enclosure_mm(), bounds);
            }
        }
    }
}
#[test]
fn actual_replay_resources_cancel_and_actual_output_mismatch_refuse_without_changing_source() {
    let source = source(false, 2);
    let a = occurrence();
    let calls = Cell::new(0usize);
    let observe = || {
        calls.set(calls.get() + 1);
        false
    };
    let output = encode(&source, Some(&a), Default::default(), Some(&observe)).unwrap();
    let r = output.receipt();
    let original = source.actual_file().to_string();
    let exact = ReplayOptions {
        resources: Options {
            maximum_work: r.aggregate_work,
            maximum_bytes: r.accounted_peak_bytes,
            ..Default::default()
        },
        ..Default::default()
    };
    assert_eq!(
        encode(&source, Some(&a), exact, None).unwrap().as_str(),
        output.as_str()
    );
    for (work, bytes) in [
        (r.aggregate_work - 1, r.accounted_peak_bytes),
        (r.aggregate_work, r.accounted_peak_bytes - 1),
    ] {
        assert!(matches!(
            encode(
                &source,
                Some(&a),
                ReplayOptions {
                    resources: Options {
                        maximum_work: work,
                        maximum_bytes: bytes,
                        ..Default::default()
                    },
                    ..exact
                },
                None
            ),
            Err(Error::Limit)
        ));
    }
    assert!(matches!(
        encode(&source, Some(&a), exact, Some(&|| true)),
        Err(Error::Cancelled)
    ));
    let total = calls.get();
    let n = Cell::new(0usize);
    assert!(matches!(
        encode(
            &source,
            Some(&a),
            Default::default(),
            Some(&|| {
                n.set(n.get() + 1);
                n.get() >= total
            })
        ),
        Err(Error::Cancelled)
    ));
    assert_eq!(source.actual_file(), original);
    assert_eq!(
        encode(&source, Some(&a), exact, None).unwrap().as_str(),
        output.as_str()
    );
    // Geometry is still a valid occupied file under another actual B, but it
    // is NOT an export of the requested A(B(source)). Material alone is not C0.
    let other = source.actual_file();
    assert!(matches!(
        qualify(
            &source,
            Some(&a),
            Arc::from(other),
            Default::default(),
            None
        ),
        Err(Error::Unsupported)
    ));
    // An output physical coefficient edit must fail exact full-record binding.
    let doc = CompoundDocument::parse(output.as_str(), Default::default(), None).unwrap();
    let surface = doc
        .records()
        .values()
        .flatten()
        .find(|p| p.name == "B_SPLINE_SURFACE_WITH_KNOTS")
        .unwrap();
    let crate::step_read::Value::List(rows) = &surface.args[3] else {
        panic!()
    };
    let crate::step_read::Value::List(row) = &rows[0] else {
        panic!()
    };
    let crate::step_read::Value::Reference(id) = row[0] else {
        panic!()
    };
    let point = &doc.records()[&id][0];
    let crate::step_read::Value::List(xyz) = &point.args[1] else {
        panic!()
    };
    let v: Vec<f64> = xyz
        .iter()
        .map(|x| match x {
            crate::step_read::Value::Number(x) => *x,
            _ => panic!(),
        })
        .collect();
    let changed = output.as_str().replacen(
        point.source.unwrap(),
        &format!("CARTESIAN_POINT('',({},{},{}))", v[0] + 0.001, v[1], v[2]),
        1,
    );
    drop(doc);
    assert_ne!(changed, output.as_str());
    assert!(matches!(
        qualify(
            &source,
            Some(&a),
            Arc::from(changed),
            Default::default(),
            None
        ),
        Err(Error::Unsupported)
    ));
    // Unsupported approximation precision and caller resources fail honestly.
    let mut bad = exact;
    bad.maximum_position_error_mm = f64::NAN;
    assert!(matches!(
        encode(&source, Some(&a), bad, None),
        Err(Error::Invalid)
    ));
    let mut low = exact;
    low.maximum_position_error_mm = 1e-30;
    assert!(matches!(
        encode(&source, Some(&a), low, None),
        Err(Error::Unsupported)
    ));
    assert_eq!(source.actual_file(), original);
}

// Independent exact-rational original operator oracle. B's integer direction
// ratios earn exact3/5,4/5, while A's coefficients are the ORIGINAL binary64
// rationals (not their intended decimal3/5,4/5). Compare every enclosure corner
// to every point in the decoded output interval, not just rounded midpoints.
fn rational_operator_oracle(
    output: &VerifiedReplay<'_, '_>,
    pose: usize,
    a: Option<&LiteralOccurrencePlacement>,
) {
    use num_rational::BigRational as Q;
    use num_traits::{Signed, Zero};
    fn q(x: f64) -> Q {
        Q::from_float(x).unwrap()
    }
    let fifth = Q::from_integer(1.into()) / Q::from_integer(5.into());
    let b = match pose {
        0 => [
            [q(1.), q(0.), q(0.)],
            [q(0.), q(1.), q(0.)],
            [q(0.), q(0.), q(1.)],
        ],
        1 => [
            [q(2.), q(0.), q(0.)],
            [q(0.), q(-2.), q(0.)],
            [q(0.), q(0.), q(-2.)],
        ],
        _ => [
            [&fifth * q(6.), &fifth * q(8.), q(0.)],
            [&fifth * q(8.), &fifth * q(-6.), q(0.)],
            [q(0.), q(0.), q(-2.)],
        ],
    };
    let origin = if pose == 0 {
        [q(0.), q(0.), q(0.)]
    } else {
        [q(2.), q(-3.), q(7.)]
    };
    let bound = output.source.local.local_sewn_enclosure_mm();
    for bits in 0..8 {
        let point = std::array::from_fn::<_, 3, _>(|j| q(bound[j][(bits >> j) & 1]));
        let mut target = std::array::from_fn::<_, 3, _>(|i| {
            (0..3).fold(origin[i].clone(), |s, j| s + &b[i][j] * &point[j])
        });
        if let Some(a) = a {
            for factor in a.transforms() {
                let m = factor.matrix();
                let t = factor.translation_mm();
                target = std::array::from_fn(|i| {
                    (0..3).fold(q(t[i]), |s, j| s + q(m[i][j]) * &target[j])
                });
            }
        }
        // Exact-Q evaluation of the decoded coefficient enclosure avoids
        // introducing new f64 map_box rounding unrelated to operator C0.
        let f = output.actual.actual_frame();
        let scale = f.scale_bounds().map(q);
        let axes = f.axes();
        let origin = f.origin_mm();
        let actual = std::array::from_fn::<_, 3, _>(|i| {
            (0..3).fold(origin[i].map(q), |v, j| {
                let c = axes[j][i].map(q);
                let products = [
                    &scale[0] * &c[0] * &point[j],
                    &scale[0] * &c[1] * &point[j],
                    &scale[1] * &c[0] * &point[j],
                    &scale[1] * &c[1] * &point[j],
                ];
                [
                    v[0].clone() + products.iter().min().unwrap(),
                    v[1].clone() + products.iter().max().unwrap(),
                ]
            })
        });
        let upper = (0..3).fold(Q::zero(), |s, i| {
            s + (&target[i] - actual[i][0].clone())
                .abs()
                .max((&target[i] - actual[i][1].clone()).abs())
        });
        assert!(
            upper <= q(output.receipt.maximum_world_position_error_mm),
            "independent A×B oracle excludes corner {bits}"
        );
    }
}
