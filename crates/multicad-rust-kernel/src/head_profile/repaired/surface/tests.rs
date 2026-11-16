use super::*;
fn profile(curved: bool) -> Profile {
    Profile::from_retained_segments(&[Segment {
        role: if curved { Role::Crown } else { Role::Cone },
        start: [if curved { 1. } else { 2. }, -1.],
        end: [if curved { 2. } else { 1. }, 0.],
        curve: if curved {
            Curve::Arc {
                center: [1., 0.],
                radius_mm: 1.,
                sweep: std::f64::consts::FRAC_PI_2,
            }
        } else {
            Curve::Line
        },
    }])
    .unwrap()
}
#[test]
fn actual_repaired_surfaces_enclose_independent_cone_and_quarter_torus_integrals() {
    for curved in [false, true] {
        let p = profile(curved);
        let source = p
            .closed_repaired_meridian(Default::default(), None)
            .unwrap();
        let result = source.surface_properties(Default::default(), None).unwrap();
        assert!(std::ptr::eq(result.source(), &source));
        assert_eq!(result.surfaces().len(), 1);
        let surface = result.surfaces()[0];
        let (area, moment) = if curved {
            (pi() * pi() + I::point(2.) * pi(), -I::point(3.) * pi())
        } else {
            (
                I::point(3.) * pi() * I::point(2.).sqrt(),
                -I::point(5.) * pi() * I::point(2.).sqrt() / I::point(3.),
            )
        };
        for (actual, expected) in [
            (surface.area_mm2, area),
            (surface.first_moment_mm3[2], moment),
            (surface.centroid_mm[2], moment / area),
        ] {
            assert!(
                actual[0] <= expected.lo && expected.hi <= actual[1],
                "{curved} {actual:?} {expected:?}"
            );
            assert!(actual[1] - actual[0] <= 1e-7);
        }
        assert_eq!(surface.first_moment_mm3[0], [0., 0.]);
        assert_eq!(surface.centroid_mm[1], [0., 0.]);
        assert_eq!(surface.segment, 0);
        assert_eq!(result.into_surfaces().len(), 1);
    }
}
#[test]
fn repaired_surface_resources_cancellation_and_retry_preserve_the_same_source() {
    let p = profile(true);
    let source = p
        .closed_repaired_meridian(Default::default(), None)
        .unwrap();
    let o = RepairedHeadSurfaceOptions::default();
    let good = source.surface_properties(o, None).unwrap();
    for options in [
        RepairedHeadSurfaceOptions {
            maximum_work: good.work_used() - 1,
            ..o
        },
        RepairedHeadSurfaceOptions {
            maximum_bytes: good.accounted_peak_bytes() - 1,
            ..o
        },
        RepairedHeadSurfaceOptions {
            maximum_cells: good.cells() - 1,
            ..o
        },
    ] {
        assert!(matches!(
            source.surface_properties(options, None),
            Err(RepairedHeadError::Limit)
        ));
    }
    let polls = std::cell::Cell::new(0usize);
    source
        .surface_properties(
            o,
            Some(&|| {
                polls.set(polls.get() + 1);
                false
            }),
        )
        .unwrap();
    for stop in [1, polls.get() / 2, polls.get()] {
        let n = std::cell::Cell::new(0usize);
        assert!(matches!(
            source.surface_properties(
                o,
                Some(&|| {
                    n.set(n.get() + 1);
                    n.get() >= stop
                })
            ),
            Err(RepairedHeadError::Cancelled)
        ));
    }
    let retry = source.surface_properties(o, None).unwrap();
    assert_eq!(retry.surfaces()[0].area_mm2, good.surfaces()[0].area_mm2);
    for bad in [
        RepairedHeadSurfaceOptions {
            area_tolerance_mm2: f64::NAN,
            ..o
        },
        RepairedHeadSurfaceOptions {
            maximum_cells: 0,
            ..o
        },
        RepairedHeadSurfaceOptions {
            maximum_work: WORK + 1,
            ..o
        },
        RepairedHeadSurfaceOptions {
            maximum_bytes: BYTES + 1,
            ..o
        },
    ] {
        assert!(matches!(
            source.surface_properties(bad, None),
            Err(RepairedHeadError::Invalid)
        ));
    }
    let line_profile = profile(false);
    let line = line_profile
        .closed_repaired_meridian(Default::default(), None)
        .unwrap();
    assert!(matches!(
        line.surface_properties(
            RepairedHeadSurfaceOptions {
                area_tolerance_mm2: 1e-30,
                ..o
            },
            None
        ),
        Err(RepairedHeadError::Unresolved)
    ));
}

#[test]
fn actual_cylindrical_crowned_and_countersunk_heads_measure_every_retained_meridian() {
    use crate::head_profile::Head;
    for head in [
        Head::Cylindrical {
            radius_mm: 2.,
            height_mm: 1.6,
        },
        Head::Crowned {
            radius_mm: 2.,
            height_mm: 1.6,
            crown_radius_mm: 3.2,
            top_land_radius_mm: 1.25,
            neck_radius_mm: 0.1,
        },
        Head::Countersunk {
            radius_mm: 3.,
            height_mm: 1.65,
            included_angle_deg: 90.,
            edge_chamfer_mm: 0.1,
            neck_radius_mm: 0.1,
        },
    ] {
        let shank = if matches!(head, Head::Countersunk { .. }) {
            1.5
        } else {
            1.
        };
        let profile = Profile::new_for_thread_union(head, shank, 0.).unwrap();
        let source = profile
            .closed_repaired_meridian(Default::default(), None)
            .unwrap();
        let actual = source.surface_properties(Default::default(), None).unwrap();
        assert_eq!(actual.surfaces().len(), source.segment_count());
        for surface in actual.surfaces() {
            assert!(surface.area_mm2[0] > 0. && surface.area_mm2[1] - surface.area_mm2[0] <= 1e-7);
            let a = source.station(surface.segment).unwrap()[1];
            let b = source.station(surface.segment + 1).unwrap()[1];
            assert!(
                surface.centroid_mm[2][0] >= a.min(b) - 1e-7
                    && surface.centroid_mm[2][1] <= a.max(b) + 1e-7
            );
        }
    }
}
