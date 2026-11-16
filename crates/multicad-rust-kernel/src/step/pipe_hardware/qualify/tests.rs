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
