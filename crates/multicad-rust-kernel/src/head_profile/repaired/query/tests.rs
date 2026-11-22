use super::*;
use crate::head_profile::Head;
use multicad_precision::PortableMath as _;
fn profile(head: Head) -> Profile {
    Profile::new_for_thread_union(head, 1., 0.).unwrap()
}
fn source(p: &Profile) -> RepairedHeadMeridian<'_> {
    p.closed_repaired_meridian(Default::default(), None)
        .unwrap()
}
fn box_at(p: [f64; 3], side: f64) -> [[f64; 2]; 3] {
    p.map(|v| [v - side / 2., v + side / 2.])
}
#[test]
fn exact_cylindrical_head_queries_keep_finite_caps_bearings_and_whole_box_semantics() {
    use RepairedHeadLocation::*;
    let p = profile(Head::Cylindrical {
        radius_mm: 2.,
        height_mm: 1.6,
    });
    let s = source(&p);
    let opts = RepairedHeadQueryOptions::default();
    for (bounds, expected) in [
        (box_at([0.1, 0.1, -0.8], 0.01), Inside),
        (box_at([1.8, 0.1, -0.8], 0.01), Inside),
        (box_at([2.1, 0.1, -0.8], 0.01), Outside),
        (box_at([0., 0., -2.], 0.01), Outside),
        (box_at([0., 0., 1.], 0.01), Outside),
        (box_at([2., 0., -0.8], 0.01), Unresolved),
        (box_at([0., 0., -1.6], 0.01), Unresolved),
        (box_at([1.5, 0., 0.], 0.01), Unresolved),
        ([[1e300, 1e300], [0., 0.], [-0.8, -0.8]], Outside),
        ([[0., 0.], [-1e300, -1e300], [-0.8, -0.8]], Outside),
        (
            [[-f64::MAX, f64::MAX], [-f64::MAX, f64::MAX], [-0.8, -0.8]],
            Unresolved,
        ),
        ([[1.5, 1.5], [1.5, 1.5], [-0.8, -0.8]], Outside),
    ] {
        let q = s.classify_box(bounds, opts, None).unwrap();
        assert_eq!(q.location(), expected, "{bounds:?}");
        assert!(std::ptr::eq(q.source(), &s));
        assert_eq!(q.bounds_mm(), bounds);
        assert_eq!(q.work_used() + q.remaining_work(), opts.maximum_work);
    }
}
#[test]
fn repaired_crown_and_neck_queries_bound_actual_source_away_from_every_surface() {
    use RepairedHeadLocation::*;
    for head in [
        Head::Crowned {
            radius_mm: 2.,
            height_mm: 1.6,
            crown_radius_mm: 3.2,
            top_land_radius_mm: 1.25,
            neck_radius_mm: 0.1,
        },
        Head::Countersunk {
            radius_mm: 2.,
            height_mm: 1.6,
            included_angle_deg: 90.,
            edge_chamfer_mm: 0.1,
            neck_radius_mm: 0.1,
        },
    ] {
        let p = profile(head);
        let s = source(&p);
        for segment in 0..s.segment_count() {
            if s.segment_certificate(segment).unwrap().progress == AxialProgress::RadialBearing {
                continue;
            }
            for t in [0.125, 0.5, 0.875] {
                // This oracle reads the original corrected source point, not
                // a nominal profile reconstruction or a displayed mesh.
                let point = s
                    .evaluate(segment, [t, t])
                    .unwrap()
                    .position_mm
                    .map(|r| r[0].midpoint(r[1]));
                for theta in [0_f64, 0.7, -2.] {
                    for (offset, expected) in [(-0.001, Inside), (0.001, Outside)] {
                        let radial = point[0] + offset;
                        let bounds = box_at(
                            [
                                radial * theta.portable_cos(),
                                radial * theta.portable_sin(),
                                point[1],
                            ],
                            1e-6,
                        );
                        assert_eq!(
                            s.classify_box(bounds, Default::default(), None)
                                .unwrap()
                                .location(),
                            expected,
                            "segment{segment} t{t} theta{theta} offset{offset} point{point:?} bounds{bounds:?}"
                        );
                    }
                }
            }
        }
    }
}
#[test]
fn invalid_work_cell_and_cancellation_queries_never_publish_partial_results() {
    use std::cell::Cell;
    let p = profile(Head::Crowned {
        radius_mm: 2.,
        height_mm: 1.6,
        crown_radius_mm: 3.2,
        top_land_radius_mm: 1.25,
        neck_radius_mm: 0.1,
    });
    let s = source(&p);
    let b = box_at([1.8, 0.1, -0.8], 1e-6);
    let opts = RepairedHeadQueryOptions::default();
    for options in [
        RepairedHeadQueryOptions {
            maximum_work: WORK + 1,
            ..opts
        },
        RepairedHeadQueryOptions {
            maximum_cells: 0,
            ..opts
        },
        RepairedHeadQueryOptions {
            maximum_cells: 1_000_001,
            ..opts
        },
    ] {
        assert_eq!(
            s.classify_box(b, options, None).err().unwrap(),
            RepairedHeadError::Invalid
        );
    }
    for bounds in [[[f64::NAN, 0.]; 3], [[1., 0.]; 3], [[0., f64::INFINITY]; 3]] {
        assert_eq!(
            s.classify_box(bounds, opts, None).err().unwrap(),
            RepairedHeadError::Invalid
        );
    }
    let count = Cell::new(0);
    let q = s
        .classify_box(
            b,
            opts,
            Some(&|| {
                count.set(count.get() + 1);
                false
            }),
        )
        .unwrap();
    let exact = q.work_used();
    assert_eq!(
        s.classify_box(
            b,
            RepairedHeadQueryOptions {
                maximum_work: exact,
                ..opts
            },
            None
        )
        .unwrap()
        .remaining_work(),
        0
    );
    for maximum_work in [0, 1, exact - 1] {
        assert_eq!(
            s.classify_box(
                b,
                RepairedHeadQueryOptions {
                    maximum_work,
                    ..opts
                },
                None
            )
            .err()
            .unwrap(),
            RepairedHeadError::Limit
        );
    }
    assert_eq!(
        s.classify_box(
            b,
            RepairedHeadQueryOptions {
                maximum_cells: 1,
                ..opts
            },
            None
        )
        .err()
        .unwrap(),
        RepairedHeadError::Limit
    );
    for stop in 0..count.get() {
        let at = Cell::new(0);
        assert_eq!(
            s.classify_box(
                b,
                opts,
                Some(&|| {
                    let n = at.get();
                    at.set(n + 1);
                    n == stop
                })
            )
            .err()
            .unwrap(),
            RepairedHeadError::Cancelled
        );
    }
    assert_eq!(
        s.classify_box(b, opts, None).unwrap().location(),
        q.location()
    );
}
