use super::*;
use std::cell::Cell;
fn profile() -> Profile {
    Profile::new_for_thread_union(
        crate::head_profile::Head::Crowned {
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
fn actual_translated_retained_head_whole_source_bounds_cover_independent_queries() {
    let original = profile();
    let parts: Vec<_> = original
        .segments()
        .map(|mut s| {
            for p in [&mut s.start, &mut s.end] {
                p[0] += 1e-8;
                p[1] += 0.5;
            }
            if let Curve::Arc { ref mut center, .. } = s.curve {
                center[0] += 1e-8;
                center[1] += 0.5;
            }
            s
        })
        .collect();
    let changed = Profile::from_retained_segments(&parts).unwrap();
    let a = original
        .closed_repaired_meridian(Default::default(), None)
        .unwrap();
    let b = changed
        .closed_repaired_meridian(Default::default(), None)
        .unwrap();
    let receipt = a
        .translated_correspondence(&b, 0., 0.5, WORK, None)
        .unwrap();
    assert!(receipt.position_mm >= 1e-8 && receipt.position_mm < 1e-6);
    assert!(receipt.derivative_mm < 1e-6);
    for i in 0..a.segment_count() {
        for j in 0..=128 {
            let t = j as f64 / 128.;
            let x = a.evaluate(i, [t, t]).unwrap();
            let y = b.evaluate(i, [t, t]).unwrap();
            let dp = [0, 1].map(|k| {
                I::new(x.position_mm[k][0], x.position_mm[k][1])
                    - I::new(y.position_mm[k][0], y.position_mm[k][1])
                    - I::point(if k == 1 { -0.5 } else { 0. })
            });
            let dd = [0, 1].map(|k| {
                I::new(x.derivative_mm[k][0], x.derivative_mm[k][1])
                    - I::new(y.derivative_mm[k][0], y.derivative_mm[k][1])
            });
            // Independent point enclosures add their own evaluation roundoff;
            // the analytic whole-domain receipt bounds the underlying sources.
            assert!(magnitude(dp).hi <= receipt.position_mm + 1e-12);
            assert!(magnitude(dd).hi <= receipt.derivative_mm + 1e-12);
        }
    }
    let used = receipt.work_used;
    assert!(matches!(
        a.translated_correspondence(&b, 0., 0.5, used - 1, None),
        Err(RepairedHeadError::Limit)
    ));
    assert!(a.translated_correspondence(&b, 0., 0.5, used, None).is_ok());
    for (s, t, w) in [
        (f64::NAN, 0., WORK),
        (0., f64::INFINITY, WORK),
        (1e8, 0., WORK),
        (0., 0., WORK + 1),
    ] {
        assert!(matches!(
            a.translated_correspondence(&b, s, t, w, None),
            Err(RepairedHeadError::Invalid)
        ));
    }
    let calls = Cell::new(0);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    a.translated_correspondence(&b, 0., 0.5, WORK, Some(&count))
        .unwrap();
    let total = calls.get();
    for stop in 1..=total {
        calls.set(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() == stop
        };
        assert!(matches!(
            a.translated_correspondence(&b, 0., 0.5, WORK, Some(&cancel)),
            Err(RepairedHeadError::Cancelled)
        ));
    }
}

#[test]
fn actual_source_piece_ownership_and_law_changes_refuse() {
    let p = profile();
    let a = p
        .closed_repaired_meridian(Default::default(), None)
        .unwrap();
    let line = Profile::from_retained_segments(&[Segment {
        start: [2., -1.],
        end: [1., 0.],
        curve: Curve::Line,
        role: Role::Cone,
    }])
    .unwrap();
    let b = line
        .closed_repaired_meridian(Default::default(), None)
        .unwrap();
    assert!(matches!(
        a.translated_correspondence(&b, 0., 0., WORK, None),
        Err(RepairedHeadError::Unresolved)
    ));
    let parts: Vec<_> = p
        .segments()
        .map(|mut s| {
            if s.role == Role::Skirt {
                s.role = Role::Shank;
            }
            s
        })
        .collect();
    let changed = Profile::from_retained_segments(&parts).unwrap();
    let b = changed
        .closed_repaired_meridian(Default::default(), None)
        .unwrap();
    assert!(matches!(
        a.translated_correspondence(&b, 0., 0., WORK, None),
        Err(RepairedHeadError::Unresolved)
    ));
}
