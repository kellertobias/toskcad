use super::*;
use crate::head_profile::{Curve, Role, Segment};
fn straight(radius: f64, height: f64) -> Profile {
    Profile::from_retained_segments(&[Segment {
        start: [radius, -1.],
        end: [radius, height - 1.],
        curve: Curve::Line,
        role: Role::Shank,
    }])
    .unwrap()
}
fn source() -> Cutter {
    let (s, d) = crate::phillips_drive::test_fixture(0);
    Cutter::new(s, d).unwrap()
}
#[test]
fn complete_head_support_envelope_retains_actual_source_and_independent_clearance() {
    let source = source();
    let closed = source
        .closed_boundary_certificate(Default::default(), None)
        .unwrap();
    let radius = source.envelope().entry_radius_mm_bounds()[1] + 1.;
    let height = source.envelope().tip_depth_mm_bounds()[1] + 1.;
    let p = straight(radius, height);
    let receipt = closed
        .certify_head_envelope(&p, Default::default(), None)
        .unwrap();
    assert!(std::ptr::eq(receipt.source(), &source));
    assert!(std::ptr::eq(receipt.profile_proposal(), &p));
    assert_eq!(receipt.cells(), 1);
    assert!(receipt.minimum_radial_clearance_mm() > 0.99);
    assert!(receipt.axial_floor_clearance_mm() > 0.99);
    assert_eq!(receipt.work_used(), closed.work_used() + 257);
    assert!(receipt.retained_budget_bytes() > std::mem::size_of_val(&receipt));
    // A changing actual radius forces subdivision; each accepted whole depth
    // interval is checked against the original narrowing Cutter source.
    let p = Profile::from_retained_segments(&[Segment {
        start: [radius, -1.],
        end: [0.05, height - 1.],
        curve: Curve::Line,
        role: Role::Cone,
    }])
    .unwrap();
    let adaptive = closed
        .certify_head_envelope(&p, Default::default(), None)
        .unwrap();
    assert!(adaptive.cells() > 1);
    assert!(adaptive.minimum_radial_clearance_mm() > 1e-7);
    assert!(matches!(
        closed.certify_head_envelope(
            &p,
            HeadEnvelopeOptions {
                maximum_cells: 1,
                ..Default::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
}
#[test]
fn head_support_limits_invalid_queries_breakthrough_and_cancel_are_atomic() {
    let source = source();
    let closed = source
        .closed_boundary_certificate(Default::default(), None)
        .unwrap();
    let radius = source.envelope().entry_radius_mm_bounds()[1] + 1.;
    let height = source.envelope().tip_depth_mm_bounds()[1] + 1.;
    let p = straight(radius, height);
    for options in [
        HeadEnvelopeOptions {
            minimum_clearance_mm: f64::NAN,
            ..Default::default()
        },
        HeadEnvelopeOptions {
            minimum_clearance_mm: 0.,
            ..Default::default()
        },
        HeadEnvelopeOptions {
            maximum_cells: 0,
            ..Default::default()
        },
        HeadEnvelopeOptions {
            maximum_cells: MAX_CELLS + 1,
            ..Default::default()
        },
        HeadEnvelopeOptions {
            maximum_work: MAX_WORK + 1,
            ..Default::default()
        },
        HeadEnvelopeOptions {
            maximum_bytes: 0,
            ..Default::default()
        },
        HeadEnvelopeOptions {
            maximum_bytes: MAX_BYTES + 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            closed.certify_head_envelope(&p, options, None),
            Err(Error::Invalid)
        ));
    }
    for work in [0, closed.work_used(), closed.work_used() + 256] {
        assert!(matches!(
            closed.certify_head_envelope(
                &p,
                HeadEnvelopeOptions {
                    maximum_work: work,
                    ..Default::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
    }
    assert!(matches!(
        closed.certify_head_envelope(
            &p,
            HeadEnvelopeOptions {
                maximum_bytes: 1,
                ..Default::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    let small = straight(radius, source.envelope().tip_depth_mm_bounds()[0]);
    assert!(matches!(
        closed.certify_head_envelope(&small, Default::default(), None),
        Err(Error::Unresolved)
    ));
    let narrow = straight(source.envelope().entry_radius_mm_bounds()[0], height);
    assert!(matches!(
        closed.certify_head_envelope(&narrow, Default::default(), None),
        Err(Error::Unresolved)
    ));
    let pinched = Profile::from_retained_segments(&[
        Segment {
            start: [radius, -1.],
            end: [0.01, -0.99],
            curve: Curve::Line,
            role: Role::Cone,
        },
        Segment {
            start: [0.01, -0.99],
            end: [0.01, height - 1.],
            curve: Curve::Line,
            role: Role::Shank,
        },
    ])
    .unwrap();
    assert!(matches!(
        closed.certify_head_envelope(&pinched, Default::default(), None),
        Err(Error::Unresolved)
    ));
    for phase in 1..=3 {
        let calls = std::cell::Cell::new(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() == phase
        };
        assert!(matches!(
            closed.certify_head_envelope(&p, Default::default(), Some(&cancel)),
            Err(Error::Cancelled)
        ));
    }
    let retry = closed
        .certify_head_envelope(&p, Default::default(), None)
        .unwrap();
    assert_eq!(retry.cells(), 1);
}

#[test]
fn standalone_head_query_counts_live_closed_atlas_at_the_exact_byte_frontier() {
    let source = source();
    let p = straight(
        source.envelope().entry_radius_mm_bounds()[1] + 1.,
        source.envelope().tip_depth_mm_bounds()[1] + 1.,
    );
    let receipt = {
        let closed = source
            .closed_boundary_certificate(Default::default(), None)
            .unwrap();
        let old_frontier = std::mem::size_of::<HeadEnvelopeClearance<'_, '_>>()
            + std::mem::size_of::<Cutter>()
            + std::mem::size_of::<Profile>()
            + HeadEnvelopeClearance::scratch_budget_bytes(1);
        let frontier = old_frontier + std::mem::size_of_val(&closed);
        for bytes in [old_frontier, frontier - 1] {
            assert!(matches!(
                closed.certify_head_envelope(
                    &p,
                    HeadEnvelopeOptions {
                        maximum_cells: 1,
                        maximum_bytes: bytes,
                        ..Default::default()
                    },
                    None
                ),
                Err(Error::Limit)
            ));
        }
        closed
            .certify_head_envelope(
                &p,
                HeadEnvelopeOptions {
                    maximum_cells: 1,
                    maximum_bytes: frontier,
                    ..Default::default()
                },
                None,
            )
            .unwrap()
    };
    // The input atlas is no longer live: the earned result retains only its
    // immutable Cutter/Profile borrows, as promised by its returned lifetimes.
    assert!(std::ptr::eq(receipt.source(), &source));
    assert_eq!(receipt.cells(), 1);
    assert!(receipt.minimum_radial_clearance_mm() > 0.99);
}
