use super::*;
#[test]
fn immutable_source_jet_cache_refuses_partial_construction_and_retries_cleanly() {
    let body = super::super::tests::cam();
    let mut w = Work {
        remaining: 10_000_000,
        interrupt: None,
    };
    let prepared = Prepared::new(
        &body,
        super::super::tests::options(),
        DEFAULT_LIMITS,
        &mut w,
    )
    .unwrap();
    let tensors = prepared
        .sides
        .iter()
        .filter_map(|s| {
            if let Side::Tensor(t) = s {
                Some(t)
            } else {
                None
            }
        })
        .collect::<Vec<_>>();
    let patch = &tensors[0].cells[0];
    let snapshot = patch
        .normal
        .iter()
        .flatten()
        .flatten()
        .map(|v| [v.lo, v.hi])
        .collect::<Vec<_>>();
    assert!(patch.normal_jet.get().is_none());
    assert_eq!(
        patch
            .taylor_normal_box(
                patch.domain,
                &mut Work {
                    remaining: 0,
                    interrupt: None
                }
            )
            .err(),
        Some(Error::Limit)
    );
    assert!(patch.normal_jet.get().is_none());
    assert_eq!(
        patch
            .taylor_normal_box(
                patch.domain,
                &mut Work {
                    remaining: 1_000_000,
                    interrupt: Some(&|| true)
                }
            )
            .err(),
        Some(Error::Cancelled)
    );
    assert!(patch.normal_jet.get().is_none());
    let bound = patch
        .taylor_normal_box(patch.domain, &mut w)
        .unwrap()
        .map(|v| [v.lo, v.hi]);
    assert!(patch.normal_jet.get().unwrap().is_some());
    assert_eq!(
        patch
            .taylor_normal_box(
                patch.domain,
                &mut Work {
                    remaining: 64,
                    interrupt: Some(&|| true)
                }
            )
            .err(),
        Some(Error::Cancelled)
    );
    let warm = patch
        .taylor_normal_box(
            patch.domain,
            &mut Work {
                remaining: 64,
                interrupt: None,
            },
        )
        .unwrap()
        .map(|v| [v.lo, v.hi]);
    assert_eq!(warm, bound);
    assert_eq!(
        patch
            .normal
            .iter()
            .flatten()
            .flatten()
            .map(|v| [v.lo, v.hi])
            .collect::<Vec<_>>(),
        snapshot
    );

    // The same finite source coefficients in an extreme auxiliary parameter
    // chart cannot make an unresolved fast enclosure authoritative. Remember
    // its numerical refusal; callers must retain their original fallback.
    let narrow = Patch {
        controls: patch.controls.clone(),
        normal: patch.normal.clone(),
        domain: [[0., 1e-300], [0., 1.]],
        normal_derivative: patch.normal_derivative,
        normal_hull: patch.normal_hull,
        curvature: patch.curvature,
        normal_jet: std::sync::OnceLock::new(),
    };
    assert_eq!(
        narrow.taylor_normal_box(narrow.domain, &mut w).err(),
        Some(Error::Unresolved)
    );
    assert!(narrow.normal_jet.get().unwrap().is_none());
    assert_eq!(
        narrow
            .taylor_normal_box(
                narrow.domain,
                &mut Work {
                    remaining: 0,
                    interrupt: None
                }
            )
            .err(),
        Some(Error::Unresolved)
    );
}
