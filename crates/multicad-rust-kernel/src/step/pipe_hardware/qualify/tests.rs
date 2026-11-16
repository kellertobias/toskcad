use super::super::test_fixture::engineering_height;
use super::*;

#[test]
fn independent_operator_comparison_checks_two_earned_physical_sources_and_real_limits() {
    let expected = engineering_height(3.);
    let actual_source = engineering_height(4.);
    let file = propose(
        &actual_source,
        Options {
            additional_resident_bytes: expected.receipt().retained_result_bytes_upper_bound,
            ..Default::default()
        },
        None,
    )
    .unwrap();
    let path =
        std::env::temp_dir().join(format!("earned-head4-helper-{}.step", std::process::id()));
    std::fs::write(&path, file.as_str()).unwrap();
    drop(file);
    let actual = std::fs::read_to_string(&path).unwrap();
    let decoded = reader::read(
        &actual,
        1e-7,
        1e-6,
        reader::Options {
            additional_resident_bytes: expected.receipt().retained_result_bytes_upper_bound
                + actual_source.receipt().retained_result_bytes_upper_bound
                + actual.capacity()
                - actual.len(),
            ..Default::default()
        },
        None,
    )
    .unwrap();
    assert_eq!(decoded.material().head_height_mm(), [6., 10.]);
    let calls = std::cell::Cell::new(0usize);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    let mut work = Work {
        remaining: MAX_WORK,
        interrupt: Some(&count),
    };
    compare(
        actual_source.shell().input(),
        decoded.material().shell().input(),
        decoded.parameter_maps(),
        &mut work,
    )
    .unwrap();
    let used = MAX_WORK - work.remaining;
    assert!(used > 0);
    // Helper-only ambiguity contract: duplicate an actual earned plain map
    // value, never its fields or the sealed decoded material. The public
    // decoder emits unique maps and does not expose this ambiguity route.
    let mut ambiguous = decoded.parameter_maps().to_vec();
    for earned in decoded.parameter_maps() {
        ambiguous.push(earned.clone());
        let mut ambiguity = Work {
            remaining: MAX_WORK,
            interrupt: None,
        };
        assert!(matches!(
            compare(
                actual_source.shell().input(),
                decoded.material().shell().input(),
                &ambiguous,
                &mut ambiguity
            ),
            Err(Error::Unsupported)
        ));
        let removed = ambiguous.pop().unwrap();
        assert_eq!(removed.target(), earned.target());
        assert_eq!(removed.actual_entity(), earned.actual_entity());
        assert_eq!(removed.mapping(), earned.mapping());
    }
    drop(ambiguous);
    // The geometry inputs below are genuine and identical; their intended
    // head is3mm high. This immutable ledger belongs to the independently
    // decoded4mm physical file and must not substitute for its operator maps.
    // Public dense binding refuses that owner mismatch earlier; this checks
    // the independent helper contract without fabricating map fields.
    let mut wrong_owner = Work {
        remaining: MAX_WORK,
        interrupt: None,
    };
    assert!(matches!(
        compare(
            expected.shell().input(),
            expected.shell().input(),
            decoded.parameter_maps(),
            &mut wrong_owner
        ),
        Err(Error::Unsupported)
    ));
    // The ledger is a small independent helper contract, not an admitted
    // Source. Every proposed index below comes from these two genuinely earned
    // and successfully compared physical source/operator tables.
    let edge_pairs = actual_source
        .shell()
        .input()
        .faces
        .iter()
        .zip(&decoded.material().shell().input().faces)
        .flat_map(|(a, b)| {
            a.outer
                .iter()
                .zip(&b.outer)
                .map(|(&(x, _), &(y, _))| (x, y))
        })
        .collect::<std::collections::BTreeMap<_, _>>()
        .into_iter()
        .take(3)
        .collect::<Vec<_>>();
    let [(a0, b0), (a1, b1), (a2, b2)]: [(usize, usize); 3] = edge_pairs
        .as_slice()
        .try_into()
        .expect("three real distinct edge mappings");
    assert!(b0 != b1 && b1 != b2 && b0 != b2);
    let mut ledger = Pairs::<128>::new();
    ledger.pair(a0, b0).unwrap();
    ledger.pair(a0, b0).unwrap(); // repeated real face incidence is an alias
    ledger.pair(a1, b1).unwrap();
    let snapshot = (ledger.forward, ledger.reverse);
    assert!(matches!(ledger.pair(a0, b1), Err(Error::Unsupported)));
    assert_eq!((ledger.forward, ledger.reverse), snapshot);
    assert!(matches!(ledger.pair(a2, b0), Err(Error::Unsupported)));
    assert_eq!((ledger.forward, ledger.reverse), snapshot);
    ledger.pair(a2, b2).unwrap();
    ledger.pair(a0, b0).unwrap();
    // Independent helper contract: both arguments below have earned occupied
    // admission. The public qualifier refuses this different height earlier in
    // dense binding; these checks specifically verify its second proof layer.
    let mut mismatch = Work {
        remaining: used,
        interrupt: None,
    };
    assert!(matches!(
        compare(
            expected.shell().input(),
            decoded.material().shell().input(),
            decoded.parameter_maps(),
            &mut mismatch
        ),
        Err(Error::Unsupported)
    ));
    let mut short = Work {
        remaining: used - 1,
        interrupt: None,
    };
    assert!(matches!(
        compare(
            actual_source.shell().input(),
            decoded.material().shell().input(),
            decoded.parameter_maps(),
            &mut short
        ),
        Err(Error::Limit)
    ));
    for stop in [1, calls.get() / 2, calls.get() - 1] {
        let visited = std::cell::Cell::new(0usize);
        let cancel = || {
            visited.set(visited.get() + 1);
            visited.get() == stop
        };
        let mut interrupted = Work {
            remaining: used,
            interrupt: Some(&cancel),
        };
        assert!(matches!(
            compare(
                actual_source.shell().input(),
                decoded.material().shell().input(),
                decoded.parameter_maps(),
                &mut interrupted
            ),
            Err(Error::Cancelled)
        ));
        assert_eq!(visited.get(), stop);
    }
    let mut retry = Work {
        remaining: used,
        interrupt: None,
    };
    compare(
        actual_source.shell().input(),
        decoded.material().shell().input(),
        decoded.parameter_maps(),
        &mut retry,
    )
    .unwrap();
    assert_eq!(retry.remaining, 0);
    assert_eq!(expected.head_height_mm(), [6., 9.]);
    assert_eq!(actual_source.head_height_mm(), [6., 10.]);
    assert_eq!(decoded.original_input(), actual);
}

#[test]
fn actual_exclusive_interior_change_must_reearn_material_before_operator_refusal() {
    use crate::step_read::pipe_hardware::tests::{InteriorChange, exclusive_interior_nextup};
    use num_rational::BigRational;
    let source = engineering_height(3.);
    let candidate = propose(&source, Options::default(), None).unwrap();
    let original_path =
        std::env::temp_dir().join(format!("interior-original-{}.step", std::process::id()));
    let changed_path =
        std::env::temp_dir().join(format!("interior-nextup-{}.step", std::process::id()));
    std::fs::write(&original_path, candidate.as_str()).unwrap();
    drop(candidate);
    let original = std::fs::read_to_string(&original_path).unwrap();
    let qualified = qualify_actual_file(
        &source,
        &original,
        1e-7,
        1e-6,
        QualificationOptions {
            additional_resident_bytes: original.capacity() - original.len() + 4096,
            ..Default::default()
        },
        None,
    )
    .unwrap();
    let snapshot_budget = qualified.material().retained_bytes_upper_bound();
    let before = qualified.material().material().shell().input().clone();
    let InteriorChange {
        text,
        tensor_entity,
        point_entity,
        pole_index,
        old_value,
        new_value,
    } = exclusive_interior_nextup(&original);
    assert!(
        BigRational::from_float(new_value).unwrap() > BigRational::from_float(old_value).unwrap()
    );
    std::fs::write(&changed_path, &text).unwrap();
    drop(text);
    drop(qualified);
    drop(original);
    let changed = std::fs::read_to_string(&changed_path).unwrap();
    // Crucially, this is an independent actual-file material admission FIRST.
    // The original sealed source is never edited, and no manufactured proof
    // stands in for the changed document's real geometry/cap/contact checks.
    let admitted = reader::read(
        &changed,
        1e-7,
        1e-6,
        reader::Options {
            additional_resident_bytes: source.receipt().retained_result_bytes_upper_bound
                + snapshot_budget
                + changed.capacity()
                - changed.len()
                + 4096,
            ..Default::default()
        },
        None,
    )
    .expect("changed actual interior must independently re-earn occupied material");
    assert_eq!(admitted.material().head_height_mm(), [6., 9.]);
    assert_eq!(
        admitted
            .material()
            .capped_thread()
            .classify([0., 0., 3.], None)
            .unwrap(),
        crate::radial_tensor::Location::Inside
    );
    let after = admitted.material().shell().input();
    assert_eq!(before.vertices, after.vertices);
    assert_eq!(before.curves, after.curves);
    assert_eq!(before.edges, after.edges);
    assert_eq!(before.faces, after.faces);
    assert_eq!(before.surfaces.len(), after.surfaces.len());
    let index = admitted
        .parameter_maps()
        .iter()
        .find(|m| m.actual_entity() == tensor_entity)
        .and_then(|m| {
            if let Target::Surface(n) = m.target() {
                Some(n)
            } else {
                None
            }
        })
        .unwrap();
    let mut scalar_changes = 0;
    for (surface_index, (a, b)) in before.surfaces.iter().zip(&after.surfaces).enumerate() {
        match (a, b) {
            (Surface::Tensor(a), Surface::Tensor(b)) => {
                assert_eq!(a.spans(), b.spans());
                assert_eq!(a.poles().len(), b.poles().len());
                if surface_index == index {
                    assert_eq!(a.poles()[pole_index][0].to_bits(), old_value.to_bits());
                    assert_eq!(b.poles()[pole_index][0].to_bits(), new_value.to_bits());
                }
                for (a, b) in a.poles().iter().zip(b.poles()) {
                    scalar_changes += a
                        .iter()
                        .zip(b)
                        .filter(|(a, b)| a.to_bits() != b.to_bits())
                        .count();
                }
            }
            _ => assert_eq!(a, b),
        }
    }
    assert_eq!(
        scalar_changes, 1,
        "all boundary/operator topology stays fixed; exactly one interior scalar changes"
    );
    eprintln!(
        "EARNED_CHANGED_INTERIOR tensor{tensor_entity} point{point_entity} scalar_changes{scalar_changes} preserved_boundarytrue"
    );
    let mut work = Work {
        remaining: MAX_WORK,
        interrupt: None,
    };
    assert!(matches!(
        compare(
            source.shell().input(),
            after,
            admitted.parameter_maps(),
            &mut work
        ),
        Err(Error::Unsupported)
    ));
    drop(admitted);
    drop(before);
    assert!(matches!(
        qualify_actual_file(
            &source,
            &changed,
            1e-7,
            1e-6,
            QualificationOptions {
                additional_resident_bytes: changed.capacity() - changed.len() + 4096,
                ..Default::default()
            },
            None
        ),
        Err(Error::Unsupported)
    ));
    drop(changed);
    let healthy = std::fs::read_to_string(&original_path).unwrap();
    let retry = qualify_actual_file(
        &source,
        &healthy,
        1e-7,
        1e-6,
        QualificationOptions {
            additional_resident_bytes: healthy.capacity() - healthy.len() + 4096,
            ..Default::default()
        },
        None,
    )
    .unwrap();
    assert_eq!(retry.as_str(), healthy);
    assert_eq!(source.head_height_mm(), [6., 9.]);
    std::fs::remove_file(&original_path).unwrap();
    std::fs::remove_file(&changed_path).unwrap();
}
