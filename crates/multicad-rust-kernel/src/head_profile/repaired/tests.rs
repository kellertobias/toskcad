use super::*;
fn profile() -> Profile {
    Profile::new_for_thread_union(
        super::super::Head::Crowned {
            radius_mm: 2.,
            height_mm: 1.6,
            crown_radius_mm: 3.2,
            top_land_radius_mm: 1.25,
            neck_radius_mm: 0.1,
        },
        1.,
        0.,
    )
    .unwrap()
}
#[test]
fn actual_crowned_source_has_shared_stations_positive_whole_domain_order_and_finite_queries() {
    let p = profile();
    let source = p
        .closed_repaired_meridian(Default::default(), None)
        .unwrap();
    assert!(std::ptr::eq(source.profile_proposal(), &p));
    assert_eq!(source.closed_face_count(), source.segment_count() + 2);
    for i in 0..source.segment_count() {
        let c = source.segment_certificate(i).unwrap();
        assert_eq!(c.stations, [i, i + 1]);
        assert!(c.minimum_radius_mm > 0. && c.minimum_progress_mm > 0.);
        let start = source.evaluate(i, [0., 0.]).unwrap();
        let end = source.evaluate(i, [1., 1.]).unwrap();
        assert_eq!(
            start.position_mm,
            source.station(i).unwrap().map(|x| [x, x])
        );
        assert_eq!(
            end.position_mm,
            source.station(i + 1).unwrap().map(|x| [x, x])
        );
        let all = source.evaluate(i, [0., 1.]).unwrap();
        assert!(
            all.position_mm
                .into_iter()
                .chain(all.derivative_mm)
                .flatten()
                .all(f64::is_finite)
        );
        if i > 0 {
            assert_eq!(
                start.position_mm,
                source.evaluate(i - 1, [1., 1.]).unwrap().position_mm
            )
        }
    }
}
#[test]
fn source_limits_cancellation_and_true_top_cap_overlap_refuse_without_mutation() {
    let p = profile();
    let source = p
        .closed_repaired_meridian(Default::default(), None)
        .unwrap();
    for options in [
        RepairedHeadOptions {
            maximum_transfer_mm: f64::NAN,
            ..Default::default()
        },
        RepairedHeadOptions {
            maximum_transfer_mm: 0.,
            ..Default::default()
        },
        RepairedHeadOptions {
            maximum_transfer_mm: 1e-6,
            ..Default::default()
        },
        RepairedHeadOptions {
            maximum_work: WORK + 1,
            ..Default::default()
        },
        RepairedHeadOptions {
            maximum_bytes: BYTES + 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            p.closed_repaired_meridian(options, None),
            Err(RepairedHeadError::Invalid)
        ));
    }
    for work in [0, 1, source.work_used() - 1] {
        assert!(matches!(
            p.closed_repaired_meridian(
                RepairedHeadOptions {
                    maximum_work: work,
                    ..Default::default()
                },
                None
            ),
            Err(RepairedHeadError::Limit)
        ));
    }
    assert!(
        p.closed_repaired_meridian(
            RepairedHeadOptions {
                maximum_work: source.work_used(),
                ..Default::default()
            },
            None
        )
        .is_ok()
    );
    let frontier = source.retained_budget_bytes() + SCRATCH;
    assert!(matches!(
        p.closed_repaired_meridian(
            RepairedHeadOptions {
                maximum_bytes: frontier - 1,
                ..Default::default()
            },
            None
        ),
        Err(RepairedHeadError::Limit)
    ));
    assert!(
        p.closed_repaired_meridian(
            RepairedHeadOptions {
                maximum_bytes: frontier,
                ..Default::default()
            },
            None
        )
        .is_ok()
    );
    assert!(matches!(
        p.closed_repaired_meridian(
            RepairedHeadOptions {
                maximum_transfer_mm: 1e-30,
                ..Default::default()
            },
            None
        ),
        Err(RepairedHeadError::Unresolved)
    ));
    // This public retained Profile is valid, but its first annulus doubles
    // back over the proposed top disk. A closed-source proof must refuse it.
    let folded = Profile::from_retained_segments(&[
        Segment {
            start: [2., -1.],
            end: [1., -1.],
            curve: Curve::Line,
            role: Role::Bearing,
        },
        Segment {
            start: [1., -1.],
            end: [1., 0.],
            curve: Curve::Line,
            role: Role::Shank,
        },
    ])
    .unwrap();
    assert!(matches!(
        folded.closed_repaired_meridian(Default::default(), None),
        Err(RepairedHeadError::Unresolved)
    ));
    let pole = Profile::new_for_thread_union(
        super::super::Head::Crowned {
            radius_mm: 2.,
            height_mm: 1.6,
            crown_radius_mm: 3.2,
            top_land_radius_mm: 0.,
            neck_radius_mm: 0.1,
        },
        1.,
        0.,
    )
    .unwrap();
    assert!(matches!(
        pole.closed_repaired_meridian(Default::default(), None),
        Err(RepairedHeadError::Unresolved)
    ));
    let calls = std::cell::Cell::new(0usize);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    p.closed_repaired_meridian(Default::default(), Some(&count))
        .unwrap();
    for phase in 1..=calls.get() {
        let polls = std::cell::Cell::new(0);
        let cancel = || {
            polls.set(polls.get() + 1);
            polls.get() == phase
        };
        assert!(matches!(
            p.closed_repaired_meridian(Default::default(), Some(&cancel)),
            Err(RepairedHeadError::Cancelled)
        ));
    }
    assert!(p.closed_repaired_meridian(Default::default(), None).is_ok());
    assert_eq!(source.station(0), Some(p.top()));
}
#[test]
fn source_query_range_resource_and_publication_contracts_preserve_identity() {
    let p = profile();
    let s = p
        .closed_repaired_meridian(Default::default(), None)
        .unwrap();
    assert!(s.station(s.segment_count() + 1).is_none());
    assert!(s.segment_certificate(s.segment_count()).is_none());
    for range in [
        [-1., 0.],
        [1., 0.],
        [0., 2.],
        [f64::NAN, 1.],
        [0., f64::INFINITY],
    ] {
        assert!(matches!(
            s.evaluate(0, range),
            Err(RepairedHeadError::Invalid)
        ));
    }
    assert!(matches!(
        s.evaluate(usize::MAX, [0., 1.]),
        Err(RepairedHeadError::Invalid)
    ));
    for (work, bytes, error) in [
        (WORK + 1, BYTES, RepairedHeadError::Invalid),
        (WORK, BYTES + 1, RepairedHeadError::Invalid),
        (s.work_used() + 256, BYTES, RepairedHeadError::Limit),
        (
            WORK,
            s.retained_budget_bytes() + SCRATCH - 1,
            RepairedHeadError::Limit,
        ),
    ] {
        assert_eq!(
            s.evaluate_with_interrupt(0, [0., 1.], work, bytes, None)
                .unwrap_err(),
            error
        );
    }
    let full = s
        .evaluate_with_interrupt(
            0,
            [0., 1.],
            s.work_used() + 257,
            s.retained_budget_bytes() + SCRATCH,
            None,
        )
        .unwrap();
    assert_eq!(full.parameter, [0., 1.]);
    for phase in 1..=2 {
        let calls = std::cell::Cell::new(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() == phase
        };
        assert!(matches!(
            s.evaluate_with_interrupt(0, [0., 1.], WORK, BYTES, Some(&cancel)),
            Err(RepairedHeadError::Cancelled)
        ));
    }
    assert_eq!(
        s.evaluate(0, [0., 1.]).unwrap().position_mm,
        full.position_mm
    );
}
#[test]
fn actual_straight_conical_and_horizontal_end_bearing_profiles_have_closed_order() {
    for parts in [
        vec![Segment {
            start: [2., -1.],
            end: [1., 0.],
            curve: Curve::Line,
            role: Role::Cone,
        }],
        vec![
            Segment {
                start: [2., -1.],
                end: [2., 0.],
                curve: Curve::Line,
                role: Role::Rim,
            },
            Segment {
                start: [2., 0.],
                end: [1., 0.],
                curve: Curve::Line,
                role: Role::Bearing,
            },
        ],
    ] {
        let p = Profile::from_retained_segments(&parts).unwrap();
        let source = p
            .closed_repaired_meridian(Default::default(), None)
            .unwrap();
        assert_eq!(source.segment_count(), parts.len());
        for (i, p) in parts.iter().enumerate() {
            let cert = source.segment_certificate(i).unwrap();
            assert_eq!(cert.role, p.role);
            assert_eq!(cert.support_transfer_mm, 0.);
            assert_eq!(cert.finite_trim_transfer_mm, 0.);
            let mid = source.evaluate(i, [0.5, 0.5]).unwrap();
            for (axis, b) in mid.position_mm.iter().enumerate() {
                let expected = (p.start[axis] + p.end[axis]) / 2.;
                assert!(b[0] <= expected && expected <= b[1]);
            }
        }
    }
}
#[test]
fn cardinal_neck_positions_and_derivatives_enclose_an_independent_exact_rational_oracle() {
    use num_rational::BigRational as Q;
    let q = |x| Q::from_float(x).unwrap();
    let p = profile();
    let source = p
        .closed_repaired_meridian(Default::default(), None)
        .unwrap();
    let (index, s) = p
        .segments()
        .enumerate()
        .find(|(_, s)| s.role == Role::Neck)
        .unwrap();
    let Curve::Arc {
        center, radius_mm, ..
    } = s.curve
    else {
        panic!("neck is an actual circular support")
    };
    assert_eq!(s.start[0], center[0]);
    assert_eq!(s.end[1], center[1]);
    assert!(s.start[1] < center[1] && s.end[0] < center[0]);
    let (mut lo, mut hi) = (Q::from_integer(1.into()), Q::from_integer(2.into()));
    let two = Q::from_integer(2.into());
    for _ in 0..180 {
        let mid = (&lo + &hi) / &two;
        if &mid * &mid < two {
            lo = mid
        } else {
            hi = mid
        }
    }
    assert!(&lo * &lo < two && &hi * &hi > two);
    let delta0 = [
        q(s.start[0]) - q(center[0]),
        q(s.start[1]) - q(center[1]) + q(radius_mm),
    ];
    let delta1 = [
        q(s.end[0]) - q(center[0]) + q(radius_mm),
        q(s.end[1]) - q(center[1]),
    ];
    let radius = q(radius_mm);
    let mid = source.evaluate(index, [0.5, 0.5]).unwrap();
    for axis in 0..2 {
        let constant = q(center[axis]) + (&delta0[axis] + &delta1[axis]) / &two;
        let expected = [&constant - &radius / &lo, &constant - &radius / &hi];
        assert!(q(mid.position_mm[axis][0]) <= expected[0]);
        assert!(expected[1] <= q(mid.position_mm[axis][1]));
        let constant = (&delta1[axis] - &delta0[axis]) * Q::new(3.into(), 2.into());
        let expected = if axis == 0 {
            [&constant - &radius * &hi, &constant - &radius * &lo]
        } else {
            [&constant + &radius * &lo, &constant + &radius * &hi]
        };
        assert!(q(mid.derivative_mm[axis][0]) <= expected[0]);
        assert!(expected[1] <= q(mid.derivative_mm[axis][1]));
    }
    // The actual normalized-angle winding v then meridian t is outward:
    // its radial normal is 2*pi*r*z', strictly positive here, independently
    // checked against the exact-root oracle rather than certificate metadata.
    assert!(mid.position_mm[0][0] > 0. && mid.derivative_mm[1][0] > 0.);
    let bound = q(source
        .segment_certificate(index)
        .unwrap()
        .support_transfer_mm);
    for delta in [delta0, delta1] {
        assert!(&delta[0] * &delta[0] + &delta[1] * &delta[1] <= &bound * &bound);
    }
}

#[test]
fn actual_countersunk_noncardinal_neck_and_socket_cap_sources_cover_all_finite_parts() {
    for head in [
        super::super::Head::Countersunk {
            radius_mm: 3.,
            height_mm: 1.65,
            included_angle_deg: 90.,
            edge_chamfer_mm: 0.1,
            neck_radius_mm: 0.1,
        },
        super::super::Head::Cylindrical {
            radius_mm: 3.,
            height_mm: 2.,
        },
    ] {
        let p = Profile::new_for_thread_union(head, 1.5, 1.).unwrap();
        let source = p
            .closed_repaired_meridian(Default::default(), None)
            .unwrap();
        let mut neck = false;
        for i in 0..source.segment_count() {
            let c = source.segment_certificate(i).unwrap();
            if c.role == Role::Neck {
                assert_eq!(c.progress, AxialProgress::Strict);
                neck = true;
            }
            let p = source.evaluate(i, [0.25, 0.75]).unwrap();
            assert!(p.position_mm[0][0] > 0.);
            if c.progress != AxialProgress::RadialBearing {
                assert!(p.derivative_mm[1][0] > 0.);
            }
            assert!(c.minimum_progress_mm > 0. && c.minimum_radius_mm > 0.);
        }
        assert_eq!(neck, matches!(head, super::super::Head::Countersunk { .. }));
    }
}
