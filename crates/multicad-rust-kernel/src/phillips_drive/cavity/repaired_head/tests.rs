use super::*;
use crate::head_profile::{Curve, Profile, Role, Segment};
fn cutter() -> Cutter {
    let (s, d) = crate::phillips_drive::test_fixture(0);
    Cutter::new(s, d).unwrap()
}
fn profile(start: f64, end: f64, height: f64) -> Profile {
    Profile::from_retained_segments(&[Segment {
        start: [start, -height / 2.],
        end: [end, height / 2.],
        curve: Curve::Line,
        role: if start == end {
            Role::Shank
        } else if start > end {
            Role::Cone
        } else {
            Role::Chamfer
        },
    }])
    .unwrap()
}
#[test]
fn whole_repaired_source_contains_real_cavity_with_coarse_and_restricted_proofs() {
    let source = cutter();
    let closed = source
        .closed_boundary_certificate(Default::default(), None)
        .unwrap();
    let radius = source.envelope.entry.hi + 1.;
    let height = source.envelope.tip.hi + 1.;
    for end in [radius, 0.05] {
        let p = profile(radius, end, height);
        let h = p
            .closed_repaired_meridian(Default::default(), None)
            .unwrap();
        let proof = closed
            .certify_repaired_head(&h, Default::default(), None)
            .unwrap();
        assert!(std::ptr::eq(proof.source(), &source));
        assert!(std::ptr::eq(proof.head(), &h));
        assert!(proof.minimum_radial_clearance_mm() > 1e-7);
        assert!(proof.axial_floor_clearance_mm() > 0.99);
        assert!(proof.retained_budget_bytes() > std::mem::size_of_val(&proof));
        if end != radius {
            assert!(proof.cells() > 1);
        } else {
            assert_eq!(proof.cells(), 1);
        }
    }
}
#[test]
fn actual_source_refusals_limits_and_every_straight_publication_poll_are_atomic() {
    let source = cutter();
    let closed = source
        .closed_boundary_certificate(Default::default(), None)
        .unwrap();
    let r = source.envelope.entry.hi + 1.;
    let height = source.envelope.tip.hi + 1.;
    let p = profile(r, r, height);
    let h = p
        .closed_repaired_meridian(Default::default(), None)
        .unwrap();
    let proof = closed
        .certify_repaired_head(&h, Default::default(), None)
        .unwrap();
    let peak = std::mem::size_of::<RepairedHeadContainment<'_, '_, '_>>()
        + std::mem::size_of_val(&closed)
        + std::mem::size_of::<Cutter>()
        + h.retained_budget_bytes()
        + 8 * std::mem::size_of::<Cell>()
        + SCRATCH;
    let options = RepairedHeadContainmentOptions {
        maximum_work: proof.work_used(),
        maximum_cells: 1,
        maximum_bytes: peak,
        ..Default::default()
    };
    assert!(closed.certify_repaired_head(&h, options, None).is_ok());
    for options in [
        RepairedHeadContainmentOptions {
            maximum_work: options.maximum_work - 1,
            ..options
        },
        RepairedHeadContainmentOptions {
            maximum_work: 0,
            ..options
        },
        RepairedHeadContainmentOptions {
            maximum_bytes: peak - 1,
            ..options
        },
    ] {
        assert!(matches!(
            closed.certify_repaired_head(&h, options, None),
            Err(Error::Limit)
        ));
    }
    for options in [
        RepairedHeadContainmentOptions {
            minimum_clearance_mm: 0.,
            ..options
        },
        RepairedHeadContainmentOptions {
            minimum_clearance_mm: f64::INFINITY,
            ..options
        },
        RepairedHeadContainmentOptions {
            maximum_cells: 0,
            ..options
        },
        RepairedHeadContainmentOptions {
            maximum_cells: CELLS + 1,
            ..options
        },
        RepairedHeadContainmentOptions {
            maximum_work: WORK + 1,
            ..options
        },
        RepairedHeadContainmentOptions {
            maximum_bytes: BYTES + 1,
            ..options
        },
    ] {
        assert!(matches!(
            closed.certify_repaired_head(&h, options, None),
            Err(Error::Invalid)
        ));
    }
    let n = std::cell::Cell::new(0);
    closed
        .certify_repaired_head(
            &h,
            options,
            Some(&|| {
                n.set(n.get() + 1);
                false
            }),
        )
        .unwrap();
    for at in 0..n.get() {
        let call = std::cell::Cell::new(0);
        assert!(matches!(
            closed.certify_repaired_head(
                &h,
                options,
                Some(&|| {
                    let old = call.get();
                    call.set(old + 1);
                    old == at
                })
            ),
            Err(Error::Cancelled)
        ));
    }
    for p in [
        profile(source.envelope.entry.lo / 2., r, height),
        profile(r, r, source.envelope.tip.lo / 2.),
    ] {
        let h = p
            .closed_repaired_meridian(Default::default(), None)
            .unwrap();
        assert!(matches!(
            closed.certify_repaired_head(&h, Default::default(), None),
            Err(Error::Unresolved)
        ));
    }
    let p = profile(r, 0.05, height);
    let h = p
        .closed_repaired_meridian(Default::default(), None)
        .unwrap();
    assert!(matches!(
        closed.certify_repaired_head(
            &h,
            RepairedHeadContainmentOptions {
                maximum_cells: 1,
                ..Default::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        closed.certify_repaired_head(
            &h,
            RepairedHeadContainmentOptions {
                maximum_work: closed.work_used() + h.work_used() + 511,
                ..Default::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
}
#[test]
fn shallow_cavity_uses_actual_varying_head_instead_of_global_shank_radius() {
    let source = cutter();
    let closed = source
        .closed_boundary_certificate(Default::default(), None)
        .unwrap();
    let r = source.envelope.entry.hi + 1.;
    let height = source.envelope.tip.hi + 1.;
    let p = profile(r, 0.05, height);
    let h = p
        .closed_repaired_meridian(Default::default(), None)
        .unwrap();
    assert!(h.segment_certificate(0).unwrap().minimum_radius_mm < source.envelope.entry.lo);
    let count = std::cell::Cell::new(0);
    let proof = closed
        .certify_repaired_head(
            &h,
            Default::default(),
            Some(&|| {
                count.set(count.get() + 1);
                false
            }),
        )
        .unwrap();
    assert!(proof.cells() > 1);
    for at in [
        2,
        3,
        4,
        5,
        6,
        7,
        count.get() / 2,
        count.get() - 2,
        count.get() - 1,
    ] {
        let calls = std::cell::Cell::new(0);
        assert!(matches!(
            closed.certify_repaired_head(
                &h,
                Default::default(),
                Some(&|| {
                    let n = calls.get();
                    calls.set(n + 1);
                    n == at
                })
            ),
            Err(Error::Cancelled)
        ));
    }
    // A physically too-thin middle skirt, after a valid wide mouth, cannot be
    // accepted merely because the top/bottom stations have enough clearance.
    let narrow = Profile::from_retained_segments(&[
        Segment {
            start: [r, -1.],
            end: [r, -0.99],
            curve: Curve::Line,
            role: Role::Skirt,
        },
        Segment {
            start: [r, -0.99],
            end: [0.05, -0.99],
            curve: Curve::Line,
            role: Role::Bearing,
        },
        Segment {
            start: [0.05, -0.99],
            end: [0.05, height - 1.],
            curve: Curve::Line,
            role: Role::Shank,
        },
    ])
    .unwrap();
    let narrow = narrow
        .closed_repaired_meridian(Default::default(), None)
        .unwrap();
    assert!(matches!(
        closed.certify_repaired_head(&narrow, Default::default(), None),
        Err(Error::Unresolved) | Err(Error::Limit)
    ));
    assert!(
        closed
            .certify_repaired_head(&h, Default::default(), None)
            .is_ok()
    );
}

#[test]
fn frontier_preparation_is_charged_before_query_allocation() {
    let source = cutter();
    let closed = source
        .closed_boundary_certificate(Default::default(), None)
        .unwrap();
    let p = profile(
        source.envelope.entry.hi + 1.,
        source.envelope.entry.hi + 1.,
        source.envelope.tip.hi + 1.,
    );
    let h = p
        .closed_repaired_meridian(Default::default(), None)
        .unwrap();
    assert!(matches!(
        closed.certify_repaired_head(
            &h,
            RepairedHeadContainmentOptions {
                maximum_work: closed.work_used() + h.work_used() + 127,
                ..Default::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
}
