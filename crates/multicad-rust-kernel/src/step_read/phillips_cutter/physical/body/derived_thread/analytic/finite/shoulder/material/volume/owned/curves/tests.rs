use super::*;
fn overlap(a: [[f64; 2]; 3], b: [[f64; 2]; 3]) {
    for (a, b) in a.into_iter().zip(b) {
        assert!(a[0] <= b[1] && b[0] <= a[1]);
    }
}
#[test]
fn real_periodic_curve_points_retain_shared_nodes_radii_and_exact_closure() {
    for left in [false, true] {
        let owner = super::super::query::tests::owned(left, 0.03);
        let mut query = owner.source_curve_query(WORK, None).unwrap();
        assert!(std::ptr::eq(query.source(), &owner));
        assert_eq!(query.scratch_upper_bound_bytes(), 65536);
        let r = owner.thread().radius_nodes_mm();
        for (band, expected_radius) in r[..4].iter().copied().enumerate() {
            assert!(query.thread_rail_second_bound_mm(band).unwrap() > 0.);
            assert!(query.thread_rim_second_bound_mm(band).unwrap() > 0.);
            assert!(query.join_rim_second_bound_mm(band).unwrap() > 0.);
            overlap(
                query.join_rim_point_mm(band, 1.).unwrap(),
                query.join_rim_point_mm((band + 1) % 4, 0.).unwrap(),
            );
            for cap in 0..2 {
                overlap(
                    query.thread_rail_point_mm(band, cap as f64).unwrap(),
                    query.thread_rim_point_mm(cap, band, 0.).unwrap(),
                );
                overlap(
                    query.thread_rim_point_mm(cap, band, 1.).unwrap(),
                    query.thread_rim_point_mm(cap, (band + 1) % 4, 0.).unwrap(),
                );
            }
            for t in [0., 0.25, 0.5, 0.75, 1.] {
                let join = query.join_rim_point_mm(band, t).unwrap().map(interval);
                let join_radius_squared = join[0].square() + join[1].square();
                let radius_squared = owner.head_bottom_station_mm()[0].powi(2);
                assert!(
                    join_radius_squared.lo <= radius_squared
                        && join_radius_squared.hi >= radius_squared
                );
                let station = expression(owner.cap_expressions_mm()[0]);
                assert!(join[2].lo <= station.lo && join[2].hi >= station.hi);
                let p = query.thread_rail_point_mm(band, t).unwrap().map(interval);
                let radius = p[0].square() + p[1].square();
                assert!(
                    radius.lo <= expected_radius * expected_radius
                        && radius.hi >= expected_radius * expected_radius
                );
            }
        }
        for station in 0..owner.head().view().segment_count() {
            for band in 0..4 {
                assert!(query.head_rim_second_bound_mm(station, band).unwrap() > 0.);
                overlap(
                    query.head_rim_point_mm(station, band, 1.).unwrap(),
                    query
                        .head_rim_point_mm(station, (band + 1) % 4, 0.)
                        .unwrap(),
                );
            }
        }
        assert_eq!(query.work_used() + query.remaining_work(), WORK);
    }
}
#[test]
fn actual_curve_batch_resources_invalid_inputs_and_cancellation_remain_distinct() {
    let owner = super::super::query::tests::owned(false, 0.03);
    for (work, error) in [
        (0, PhillipsReadError::Limit),
        (8192, PhillipsReadError::Limit),
        (WORK + 1, PhillipsReadError::Invalid),
    ] {
        assert!(matches!(owner.source_curve_query(work,None),Err(e) if e==error));
    }
    assert!(matches!(
        owner.source_curve_query(WORK, Some(&|| true)),
        Err(PhillipsReadError::Cancelled)
    ));
    let mut exact = owner.source_curve_query(8193 + 4352, None).unwrap();
    exact.thread_rail_point_mm(0, 0.5).unwrap();
    assert_eq!(exact.remaining_work(), 0);
    assert!(matches!(
        exact.thread_rail_point_mm(0, 0.5),
        Err(PhillipsReadError::Limit)
    ));
    let mut short = owner.source_curve_query(8193 + 4351, None).unwrap();
    assert!(matches!(
        short.thread_rail_point_mm(0, 0.5),
        Err(PhillipsReadError::Limit)
    ));
    let calls = Cell::new(0);
    let cancel = || {
        let n = calls.get() + 1;
        calls.set(n);
        n >= 5
    };
    let mut late = owner.source_curve_query(WORK, Some(&cancel)).unwrap();
    assert!(matches!(
        late.thread_rail_point_mm(0, 0.5),
        Err(PhillipsReadError::Cancelled)
    ));
    let mut q = owner.source_curve_query(WORK, None).unwrap();
    for t in [f64::NAN, f64::INFINITY, -1., 2.] {
        assert!(matches!(
            q.thread_rail_point_mm(0, t),
            Err(PhillipsReadError::Invalid)
        ));
    }
    assert!(matches!(
        q.thread_rail_point_mm(4, 0.5),
        Err(PhillipsReadError::Invalid)
    ));
    assert!(matches!(
        q.thread_rim_point_mm(2, 0, 0.5),
        Err(PhillipsReadError::Invalid)
    ));
    assert!(matches!(
        q.join_rim_point_mm(4, 0.5),
        Err(PhillipsReadError::Invalid)
    ));
    assert!(matches!(
        q.join_rim_point_mm(0, f64::NAN),
        Err(PhillipsReadError::Invalid)
    ));
    assert!(matches!(
        q.head_rim_point_mm(100, 0, 0.5),
        Err(PhillipsReadError::Invalid)
    ));
    assert!(matches!(
        q.thread_rail_second_bound_mm(4),
        Err(PhillipsReadError::Invalid)
    ));
    assert!(matches!(
        q.thread_rim_second_bound_mm(4),
        Err(PhillipsReadError::Invalid)
    ));
    assert!(matches!(
        q.join_rim_second_bound_mm(4),
        Err(PhillipsReadError::Invalid)
    ));
    assert!(matches!(
        q.head_rim_second_bound_mm(100, 0),
        Err(PhillipsReadError::Invalid)
    ));
}
