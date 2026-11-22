use super::*;
use crate::step_read::phillips_cutter::physical::body::{
    common_mouth::tests::derived_thread_probe::with_head,
    derived_thread::degree11::tests::actual_m3_parts,
};
use num_rational::BigRational as Q;
pub(crate) fn actual_owned(
    label: &str,
    left: bool,
    posed: bool,
    starts: u32,
) -> OwnedDegree11Material {
    let (mut author, mut frame) = actual_m3_parts(true, left);
    if starts != 1 {
        let mut spec = author.thread().spec();
        spec.starts = starts;
        author = crate::phillips_drive::PhillipsScrew::new(
            crate::helical_thread::Thread::new(spec).unwrap(),
            *author.profile(),
            author.spec(),
            author.cavity().clone(),
            100_000,
            &|| false,
        )
        .unwrap();
    }
    if posed {
        frame = [[3., -4., 7.], [0., 1., 0.], [0., 0., 1.]];
    }
    let raw =
        crate::step::write_phillips_screw(&author, frame[0], frame[1], frame[2], 1e-7).unwrap();
    let path = std::env::temp_dir().join(format!(
        "actual-d11-query-{label}-{left}-{posed}-{starts}-{}.step",
        std::process::id()
    ));
    std::fs::write(&path, &raw).unwrap();
    drop(raw);
    let text = std::fs::read_to_string(&path).unwrap();
    std::fs::remove_file(&path).unwrap();
    let mut result = None;
    with_head(&text, &mut |head| {
        let occupied = head
            .degree11_thread_carrier(Default::default(), None)
            .unwrap()
            .certify_occupied(Default::default(), None)
            .unwrap();
        let material = head
            .closed_degree11_material(&occupied, Default::default(), None)
            .unwrap();
        let volume = material.volume_bounds(Default::default(), None).unwrap();
        result = Some(volume.owned_snapshot(Default::default(), None).unwrap());
    });
    drop(text);
    drop(author);
    result.unwrap()
}

#[test]
fn actual_m3_owned_d11_whole_source_face_surface_meets_distinct_requested_goals() {
    let owner = actual_owned("surface-whole-face", false, false, 1);
    // Same source-local extent and literal default linear policy as the
    // Native legacy contract; no Core dependency or inferred authored size.
    let extent = owner
        .local_enclosure_bounds_mm()
        .into_iter()
        .flatten()
        .map(f64::abs)
        .fold(1_f64, f64::max);
    let linear = 1e-7;
    let source = owner.thread_faces()[0].source().clone();
    drop(owner);
    eprintln!("ACTUAL_M3_D11_WHOLE_SURFACE spans{:?}", source.spans());
    let original = source.clone();
    for (label, area, moment, centroid) in [
        ("initial", 1e-3, 1e-3, 1e-3),
        ("default", 1e-5, 1e-5, 1e-5),
        ("native_proxy", 1e-7, 1e-7, 1e-8),
        (
            "native_contract",
            linear * extent / 16.,
            linear * extent * extent / 16.,
            linear / 16.,
        ),
    ] {
        let options = crate::tapered_thread::SurfaceOptions11 {
            area_width_mm2: area,
            first_moment_width_mm3: moment,
            centroid_width_mm: centroid,
            ..Default::default()
        };
        let result = source.surface_properties_with_diagnostics(options, None);
        eprintln!("ACTUAL_M3_D11_WHOLE_SURFACE {label} {result:?}");
        let result = result.unwrap();
        assert!(result.area_mm2[0] > 0.);
        assert!(I::new(result.area_mm2[0], result.area_mm2[1]).width() <= area);
        for b in result.first_moment_mm3 {
            assert!(I::new(b[0], b[1]).width() <= moment);
        }
        for b in result.centroid_mm {
            assert!(I::new(b[0], b[1]).width() <= centroid);
        }
        assert_eq!(result.work_by_stage.iter().sum::<usize>(), result.work_used);
        assert!(result.work_used <= options.max_work);
        assert!(result.accounted_peak_bytes <= options.max_bytes);
        assert_eq!(source, original);
    }
}

#[test]
fn actual_m3_owned_d11_selected_source_cell_surface_survives_all_owner_drop() {
    use crate::tapered_thread::{RetainedTensor11, SurfaceOptions11};
    for left in [false, true] {
        // This is one unchanged actual root tensor net, not a whole-face or
        // whole-screw area claim. All original STEP/parser/author owners died
        // inside actual_owned before these retained coefficients are accessed.
        let owner = actual_owned("surface-cell", left, false, 1);
        let source = owner.thread_faces()[0].source();
        let height = source.dimensions()[1];
        let v = source.spans()[1] / 2;
        let poles = (0..12)
            .flat_map(|i| (0..12).map(move |j| source.poles()[i * height + 11 * v + j]))
            .collect();
        let selected = RetainedTensor11::new([1, 1], poles).unwrap();
        let before = selected.clone();
        drop(owner);
        let options = SurfaceOptions11 {
            area_width_mm2: 1e-5,
            first_moment_width_mm3: 1e-5,
            centroid_width_mm: 1e-5,
            ..Default::default()
        };
        let result = selected.surface_properties(options, None).unwrap();
        eprintln!("ACTUAL_M3_D11_SURFACE left{left} {result:?}");
        assert!(result.area_mm2[0] > 0.);
        assert!(I::new(result.area_mm2[0], result.area_mm2[1]).width() <= options.area_width_mm2);
        for k in 0..3 {
            let min = selected
                .poles()
                .iter()
                .map(|p| p[k])
                .fold(f64::INFINITY, f64::min);
            let max = selected
                .poles()
                .iter()
                .map(|p| p[k])
                .fold(f64::NEG_INFINITY, f64::max);
            assert!(result.centroid_mm[k][0] <= max && result.centroid_mm[k][1] >= min);
            assert!(
                I::new(result.centroid_mm[k][0], result.centroid_mm[k][1]).width()
                    <= options.centroid_width_mm
            );
            assert!(
                I::new(result.first_moment_mm3[k][0], result.first_moment_mm3[k][1]).width()
                    <= options.first_moment_width_mm3
            );
        }
        assert_eq!(selected, before);
        assert_eq!(selected.surface_properties(options, None).unwrap(), result);
    }
}
fn point(p: Point) -> [[f64; 2]; 3] {
    p.map(|x| [x, x])
}
// Original cell corners are EXACT actual surface points. Both radial scalings
// below are independently proved exact binary-rational operations; keeping z
// fixed therefore uses the earned radial graph, not sampled/inferred geometry.
fn witnesses(source: &OwnedDegree11Material) -> [Point; 3] {
    let [low, high] = source.cap_stations_mm();
    for face in source.thread_faces() {
        let [su, sv] = face.source().spans();
        let height = face.source().dimensions()[1];
        for u in 0..=su {
            for v in 1..sv {
                let p = face.source().poles()[11 * u * height + 11 * v];
                if p[2] <= low || p[2] >= high {
                    continue;
                }
                for scale in [0.875, 0.9375, 0.75] {
                    let a = [p[0] * scale, p[1] * scale, p[2]];
                    let b = [2. * a[0], 2. * a[1], p[2]];
                    if !(0..2).all(|k| {
                        Q::from_float(a[k]).unwrap()
                            == Q::from_float(p[k]).unwrap() * Q::from_float(scale).unwrap()
                            && Q::from_float(b[k]).unwrap()
                                == Q::from_float(p[k]).unwrap() * Q::from_float(2. * scale).unwrap()
                    }) {
                        continue;
                    }
                    let r = (I::point(a[0]).square() + I::point(a[1]).square()).sqrt();
                    if r.lo <= source.occupied_core_radius_mm() {
                        continue;
                    }
                    if source
                        .head()
                        .view()
                        .classify_box(point(a), Default::default(), None)
                        .unwrap()
                        .location()
                        != RepairedHeadLocation::Outside
                    {
                        continue;
                    }
                    let cavity = [
                        I::point(a[0]),
                        I::point(a[1]),
                        I::point(a[2]) - I::point(source.chart().outer_top_z),
                    ];
                    if source
                        .cavity()
                        .classify_box_mm(cavity.map(|x| [x.lo, x.hi]))
                        .unwrap()
                        != C::Outside
                    {
                        continue;
                    }
                    eprintln!(
                        "ACTUAL_D11_RAYS core{} radius{} scale{scale} corner{:?}",
                        source.occupied_core_radius_mm(),
                        r.lo,
                        p
                    );
                    return [a, b, p];
                }
            }
        }
    }
    panic!("actual literal corner with exact dyadic rays outside the guaranteed core");
}
#[test]
fn actual_owned_d11_points_and_whole_boxes_survive_owner_drop_and_proper_pose() {
    for (left, posed, starts) in [
        (false, false, 1),
        (true, false, 1),
        (false, true, 1),
        (true, true, 1),
        (false, false, 2),
    ] {
        let source = actual_owned("points", left, posed, starts);
        let [inside, outside, boundary] = witnesses(&source);
        let a = source
            .classify_local_box(point(inside), Default::default(), None)
            .unwrap();
        assert_eq!(
            a.location(),
            C::Inside,
            "left{left} posed{posed} starts{starts}"
        );
        assert!(a.visited_thread_cells() > 0);
        assert!(std::ptr::eq(a.source(), &source));
        assert_eq!(a.bounds_mm(), point(inside));
        let b = source
            .classify_local_box(point(outside), Default::default(), None)
            .unwrap();
        assert_eq!(b.location(), C::Outside);
        assert!(b.visited_thread_cells() > 0);
        let mixed = std::array::from_fn(|k| [inside[k].min(outside[k]), inside[k].max(outside[k])]);
        assert_eq!(
            source
                .classify_local_box(mixed, Default::default(), None)
                .unwrap()
                .location(),
            C::UnresolvedBoundary
        );
        assert_eq!(
            source
                .classify_local_box(point(boundary), Default::default(), None)
                .unwrap()
                .location(),
            C::UnresolvedBoundary
        );
        assert_eq!(
            source
                .classify_local_box([[1e6, 1e6 + 1.]; 3], Default::default(), None)
                .unwrap()
                .location(),
            C::Outside
        );
        // The original left-hand fixture retains its non-cardinal proper
        // [.8,.6,0]/[.36,-.48,-.8] frame when no cardinal override is requested.
        let top = source.head().view().station(0).unwrap();
        assert_eq!(top[1], source.chart().outer_top_z);
        // Inside the complete source enclosure but outside the actual head's
        // radial graph in its own strict slab; thread is axially disjoint.
        let air = [
            (top[0] + source.local_enclosure_bounds_mm()[0][1]) * 0.5,
            0.,
            top[1] + 0.001,
        ];
        assert!(air[2] < source.cap_stations_mm()[0]);
        assert_eq!(
            source
                .head()
                .view()
                .classify_box(point(air), Default::default(), None)
                .unwrap()
                .location(),
            RepairedHeadLocation::Outside
        );
        assert_eq!(
            source
                .classify_local_box(point(air), Default::default(), None)
                .unwrap()
                .location(),
            C::Outside
        );
        let void = [0., 0., top[1] + 0.1];
        let head_point = [0.99 * top[0], 0., top[1] + 0.001];
        let box_at = |p: Point| p.map(|x| [x - 1e-7, x + 1e-7]);
        assert_eq!(
            source
                .classify_local_box(box_at(void), Default::default(), None)
                .unwrap()
                .location(),
            C::Outside
        );
        assert_eq!(
            source
                .classify_local_box(box_at(head_point), Default::default(), None)
                .unwrap()
                .location(),
            C::Inside
        );
        let cavity_crossing = [[0., head_point[0]], [0., 0.], [void[2]; 2]];
        assert_eq!(
            source
                .classify_local_box(cavity_crossing, Default::default(), None)
                .unwrap()
                .location(),
            C::UnresolvedBoundary
        );
        let axis_crossing = [
            [0_f64.min(outside[0]), 0_f64.max(outside[0])],
            [0_f64.min(outside[1]), 0_f64.max(outside[1])],
            [inside[2]; 2],
        ];
        assert_eq!(
            source
                .classify_local_box(axis_crossing, Default::default(), None)
                .unwrap()
                .location(),
            C::UnresolvedBoundary
        );
        let [low, high] = source.cap_stations_mm();
        let delta = (4. * source.transfer_collar().local_total_mm).max(1e-7);
        let small = source.occupied_core_radius_mm() / 4.;
        assert!(2. * delta < high - low);
        for cap in [low, high] {
            let crossing = [[-small, small], [-small, small], [cap - delta, cap + delta]];
            assert_eq!(
                source
                    .classify_local_box(crossing, Default::default(), None)
                    .unwrap()
                    .location(),
                C::UnresolvedBoundary
            );
        }
        assert_eq!(
            source
                .classify_local_box(
                    [
                        [-small, small],
                        [-small, small],
                        [low + delta, low + 2. * delta]
                    ],
                    Default::default(),
                    None
                )
                .unwrap()
                .location(),
            C::Inside
        );
        assert_eq!(
            source
                .classify_local_box(
                    [
                        [-small, small],
                        [-small, small],
                        [high + delta, high + 2. * delta]
                    ],
                    Default::default(),
                    None
                )
                .unwrap()
                .location(),
            C::Outside
        );
        let enclosure = source.local_enclosure_bounds_mm();
        assert_eq!(
            source
                .classify_local_box(enclosure, Default::default(), None)
                .unwrap()
                .location(),
            C::UnresolvedBoundary
        );
        assert!(!source.original_input().is_empty());
        eprintln!(
            "ACTUAL_D11_QUERY left{left} pose{posed} starts{starts} insidework{} outsidework{} cells{} peak{}",
            a.work_used(),
            b.work_used(),
            a.visited_thread_cells(),
            a.accounted_peak_bytes()
        );
    }
}
#[test]
fn actual_owned_d11_query_work_bytes_cells_cancel_and_retry_are_atomic() {
    let source = actual_owned("resources", false, false, 1);
    let [inside, _, boundary] = witnesses(&source);
    let options = Degree11MaterialQueryOptions::default();
    let observed = std::cell::Cell::new(0usize);
    let count = || {
        observed.set(observed.get() + 1);
        false
    };
    let baseline = source
        .classify_local_box(point(inside), options, Some(&count))
        .unwrap();
    let used = baseline.work_used();
    let peak = baseline.accounted_peak_bytes();
    assert_eq!(baseline.remaining_work(), options.maximum_work - used);
    assert!(matches!(
        source.classify_local_box(
            point(inside),
            Degree11MaterialQueryOptions {
                maximum_work: used - 1,
                ..options
            },
            None
        ),
        Err(Degree11MaterialQueryError::Limit)
    ));
    assert_eq!(
        source
            .classify_local_box(
                point(inside),
                Degree11MaterialQueryOptions {
                    maximum_work: used,
                    ..options
                },
                None
            )
            .unwrap()
            .location(),
        C::Inside
    );
    assert!(matches!(
        source.classify_local_box(
            point(inside),
            Degree11MaterialQueryOptions {
                maximum_bytes: peak - 1,
                ..options
            },
            None
        ),
        Err(Degree11MaterialQueryError::Limit)
    ));
    assert_eq!(
        source
            .classify_local_box(
                point(inside),
                Degree11MaterialQueryOptions {
                    maximum_bytes: peak,
                    ..options
                },
                None
            )
            .unwrap()
            .location(),
        C::Inside
    );
    assert!(matches!(
        source.classify_local_box(
            point(inside),
            Degree11MaterialQueryOptions {
                maximum_thread_cells: 1,
                ..options
            },
            None
        ),
        Err(Degree11MaterialQueryError::Limit)
    ));
    assert_eq!(
        source
            .classify_local_box(
                point(boundary),
                Degree11MaterialQueryOptions {
                    maximum_thread_depth: 0,
                    ..options
                },
                None
            )
            .unwrap()
            .location(),
        C::UnresolvedBoundary
    );
    for stop in [
        1,
        3,
        4,
        5,
        observed.get() / 2,
        observed.get() - 1,
        observed.get(),
    ] {
        let visited = std::cell::Cell::new(0usize);
        let cancel = || {
            visited.set(visited.get() + 1);
            visited.get() == stop
        };
        assert!(matches!(
            source.classify_local_box(point(inside), options, Some(&cancel)),
            Err(Degree11MaterialQueryError::Cancelled)
                | Err(Degree11MaterialQueryError::Head(
                    RepairedHeadError::Cancelled
                ))
        ));
        assert_eq!(visited.get(), stop);
        let retry = source
            .classify_local_box(point(inside), options, None)
            .unwrap();
        assert_eq!(retry.location(), C::Inside);
        assert_eq!(retry.work_used(), used);
    }
    // The real head exclusion owns its error until its two checkpoint credits
    // are available; only then can the actual D11 thread consume the balance.
    assert!(matches!(
        source.classify_local_box(
            point(inside),
            Degree11MaterialQueryOptions {
                maximum_work: 8192,
                ..options
            },
            None
        ),
        Err(Degree11MaterialQueryError::Head(RepairedHeadError::Limit))
    ));
    assert!(matches!(
        source.classify_local_box(
            point(inside),
            Degree11MaterialQueryOptions {
                maximum_work: 8194,
                ..options
            },
            None
        ),
        Err(Degree11MaterialQueryError::Limit)
    ));
    for bad in [
        Degree11MaterialQueryOptions {
            maximum_work: MAX_QUERY_WORK + 1,
            ..options
        },
        Degree11MaterialQueryOptions {
            maximum_bytes: MAX_BYTES + 1,
            ..options
        },
        Degree11MaterialQueryOptions {
            maximum_thread_cells: 0,
            ..options
        },
        Degree11MaterialQueryOptions {
            maximum_head_cells: 0,
            ..options
        },
        Degree11MaterialQueryOptions {
            maximum_thread_depth: MAX_THREAD_DEPTH + 1,
            ..options
        },
    ] {
        assert!(matches!(
            source.classify_local_box(point(inside), bad, None),
            Err(Degree11MaterialQueryError::Invalid)
        ));
    }
    for bad in [
        [[f64::NAN, 1.]; 3],
        [[1., 0.]; 3],
        [[0., f64::INFINITY]; 3],
        [[0., 1e8]; 3],
    ] {
        assert!(matches!(
            source.classify_local_box(bad, options, None),
            Err(Degree11MaterialQueryError::Invalid)
        ));
    }
}
