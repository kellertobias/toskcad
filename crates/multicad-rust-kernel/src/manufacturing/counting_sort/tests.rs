use super::*;
fn records(vertices: usize, count: usize) -> Vec<[u32; 3]> {
    let mut seed = 0x715e6755u64;
    (0..count)
        .map(|i| {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            let a = (seed >> 32) as usize % vertices;
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            let b = (seed >> 32) as usize % vertices;
            [a.min(b) as u32, a.max(b) as u32, i as u32]
        })
        .collect()
}
#[test]
fn indexed_endpoint_order_matches_independent_full_lexicographic_sort() {
    for (vertices, count) in [(1, 1), (7, 1024), (53, 4096), (2048, 17)] {
        let original = records(vertices, count);
        let mut expected = original.clone();
        expected.sort_unstable();
        let mut actual = original;
        super::super::sort_counted(
            &mut actual,
            vertices,
            &mut Budget::<_, true>(&mut |_, _| Ok::<_, ()>(())),
        )
        .unwrap();
        assert_eq!(actual, expected);
    }
}
#[test]
fn actual_bounded_buffers_work_frontiers_and_cancellation_preserve_retry() {
    let original = records(53, 4096);
    let mut expected = original.clone();
    expected.sort_unstable();
    let mut used = 0;
    let mut peak = 0;
    let mut calls = 0;
    let mut actual = original.clone();
    sort(
        &mut actual,
        53,
        &mut Budget::<_, true>(&mut |n, bytes| {
            used += n;
            peak = peak.max(bytes.unwrap_or(0));
            calls += 1;
            Ok::<_, ()>(())
        }),
    )
    .unwrap();
    assert_eq!(actual, expected);
    // The sparse-domain fallback is the original comparison sort, applied
    // to the same valid endpoint records and third-key order.
    let mut baseline = original.clone();
    let mut comparison_work = 0;
    super::super::sort_counted(
        &mut baseline,
        original.len() + 1,
        &mut Budget::<_, true>(&mut |n, _| {
            comparison_work += n;
            Ok::<_, ()>(())
        }),
    )
    .unwrap();
    assert_eq!(baseline, expected);
    assert!(used < comparison_work / 2);
    for limit in [used - 1, used] {
        let mut left = limit;
        let result = sort(
            &mut original.clone(),
            53,
            &mut Budget::<_, true>(&mut |n, _| {
                left = left.checked_sub(n).ok_or("work")?;
                Ok(())
            }),
        );
        assert_eq!(
            result,
            if limit == used {
                Ok(())
            } else {
                Err(CountedError::Interrupted("work"))
            }
        );
    }
    for limit in [peak - 1, peak] {
        let result = sort(
            &mut original.clone(),
            53,
            &mut Budget::<_, true>(&mut |_, bytes| {
                if bytes.is_some_and(|n| n > limit) {
                    Err("bytes")
                } else {
                    Ok(())
                }
            }),
        );
        assert_eq!(
            result,
            if limit == peak {
                Ok(())
            } else {
                Err(CountedError::Interrupted("bytes"))
            }
        );
    }
    // Include cancellation after scratch allocation, at its actual-capacity
    // checkpoint, before the first counting pass can mutate any records.
    for stop in [1, 2, 3, calls / 2, calls] {
        let mut seen = 0;
        let result = sort(
            &mut original.clone(),
            53,
            &mut Budget::<_, true>(&mut |_, _| {
                seen += 1;
                if seen == stop { Err("cancel") } else { Ok(()) }
            }),
        );
        assert_eq!(result, Err(CountedError::Interrupted("cancel")));
    }
    let mut retry = original;
    sort(
        &mut retry,
        53,
        &mut Budget::<_, true>(&mut |_, _| Ok::<_, ()>(())),
    )
    .unwrap();
    assert_eq!(retry, expected);
    assert_eq!(storage(usize::MAX, 1, 1), Err(MeshError::BudgetExceeded));
    assert_eq!(storage(0, 0, usize::MAX), Err(MeshError::BudgetExceeded));
}
