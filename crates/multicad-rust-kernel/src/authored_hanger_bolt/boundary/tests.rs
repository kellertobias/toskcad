use super::*;
#[test]
fn actual_source_incidence_retains_all_components_and_annular_loops() {
    for left in [false, true] {
        let source = Hanger::new(super::super::tests::fixture(left, false), None).unwrap();
        let b = source.boundary(BoundaryOptions::default(), None).unwrap();
        assert_eq!(b.vertices().len(), 18);
        assert_eq!(b.curves().len(), 26);
        for i in 0..source.face_count() {
            let key = source.face_key(i).unwrap();
            let loops = b.face_loops(key).unwrap();
            assert!(!loops[0].is_empty());
            assert_eq!(
                !loops[1].is_empty(),
                matches!(
                    key,
                    FaceKey::PlainWall | FaceKey::WoodShoulder | FaceKey::MetricShoulder
                )
            );
        }
    }
}
fn source(left: bool, starts: u32, nonbinary: bool) -> Hanger {
    let mut spec = super::super::tests::fixture(left, nonbinary);
    let mut m = spec.metric.spec();
    m.starts = starts;
    spec.metric = AuthoredMetricThread::new(m).unwrap();
    let mut wood = spec.wood.spec();
    wood.starts = starts;
    spec.wood = AuthoredThread::new(
        wood,
        crate::authored_flat_thread::SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-6,
        },
    )
    .unwrap();
    Hanger::new(spec, None).unwrap()
}
fn overlap(a: [f64; 2], b: [f64; 2]) -> bool {
    a[0] <= b[1] && b[0] <= a[1]
}
#[test]
fn source_curve_domains_close_at_symbolic_vertices_for_hands_starts_and_nonbinary_stations() {
    for (left, n, nb) in [(false, 1, false), (true, 2, true), (false, 8, false)] {
        let s = source(left, n, nb);
        let b = s.boundary(BoundaryOptions::default(), None).unwrap();
        assert_eq!(b.source(), &s);
        assert!(b.retained_bytes() < 1 << 20);
        for c in b.curves() {
            for (end, t) in [0., 1.].into_iter().enumerate() {
                let p = b.curve_bounds(c.key, [t, t], None).unwrap().point_mm;
                let v = b.vertex_bounds(c.vertices[end], None).unwrap();
                for axis in 0..3 {
                    assert!(
                        overlap(p[axis], v[axis]),
                        "{:?} end{end} axis{axis} {p:?} {v:?}",
                        c.key
                    );
                }
            }
            let whole = b.curve_bounds(c.key, [0., 1.], None).unwrap();
            let point = b.curve_bounds(c.key, [0.37, 0.37], None).unwrap();
            for axis in 0..3 {
                assert!(whole.point_mm[axis][0] <= point.point_mm[axis][0]);
                assert!(whole.point_mm[axis][1] >= point.point_mm[axis][1]);
            }
            // Actual finite differences are enclosed by the analytic first
            // derivative plus the whole-source C2 mean-value remainder.
            let h = 1e-5;
            let a = b.curve_bounds(c.key, [0.37 - h, 0.37 - h], None).unwrap();
            let z = b.curve_bounds(c.key, [0.37 + h, 0.37 + h], None).unwrap();
            for axis in 0..3 {
                let slope = (enclosure_for_test(z.point_mm[axis])
                    - enclosure_for_test(a.point_mm[axis]))
                    / I::point(2. * h);
                let d = enclosure_for_test(point.first_derivative_mm[axis]);
                let reserve =
                    enclosure_for_test(whole.second_derivative_mm[axis]).magnitude() * h + 1e-7;
                assert!(
                    slope.lo >= d.lo - reserve && slope.hi <= d.hi + reserve,
                    "{:?} axis{axis}",
                    c.key
                );
            }
        }
    }
}
fn enclosure_for_test(a: [f64; 2]) -> I {
    I::new(a[0], a[1])
}
#[test]
fn source_curve_queries_refuse_foreign_domains_and_cancel_without_mutating_ownership() {
    let s = source(false, 1, false);
    let b = s.boundary(BoundaryOptions::default(), None).unwrap();
    let valid = b.curves()[0].key;
    for range in [[-0.1, 0.], [1., 1.1], [0.5, 0.4], [f64::NAN, 0.]] {
        assert_eq!(b.curve_bounds(valid, range, None), Err(Error::Invalid));
    }
    let foreign = CurveKey::Rail {
        component: Component::Metric,
        period: 1,
        knot: Knot::CrestEnd,
    };
    assert_eq!(b.curve_bounds(foreign, [0., 1.], None), Err(Error::Invalid));
    assert_eq!(
        b.vertex_bounds(
            VertexKey::Thread {
                component: Component::Metric,
                end: End::Start,
                period: 1,
                knot: Knot::CrestEnd
            },
            None
        ),
        Err(Error::Invalid)
    );
    assert!(
        b.face_loops(FaceKey::MetricThread {
            start: 1,
            band: ThreadBandRole::Root
        })
        .is_none()
    );
    for at in [1, 2] {
        let count = std::cell::Cell::new(0);
        let interrupt = || {
            count.set(count.get() + 1);
            count.get() == at
        };
        assert_eq!(
            b.curve_bounds(valid, [0.1, 0.2], Some(&interrupt)),
            Err(Error::Cancelled)
        );
        count.set(0);
        assert_eq!(
            b.vertex_bounds(b.vertices()[0], Some(&interrupt)),
            Err(Error::Cancelled)
        );
    }
    assert_eq!(b.source(), &s);
}
#[test]
fn construction_memory_work_and_final_publication_have_real_separate_limits() {
    let s = source(false, 2, true);
    let b = s.boundary(BoundaryOptions::default(), None).unwrap();
    let required = b.retained_bytes() + SCRATCH;
    assert!(
        s.boundary(
            BoundaryOptions {
                maximum_bytes: required,
                ..BoundaryOptions::default()
            },
            None
        )
        .is_ok()
    );
    assert!(matches!(
        s.boundary(
            BoundaryOptions {
                maximum_bytes: required - 1,
                ..BoundaryOptions::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        s.boundary(
            BoundaryOptions {
                maximum_work: 1,
                ..BoundaryOptions::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    for options in [
        BoundaryOptions {
            maximum_work: 0,
            ..BoundaryOptions::default()
        },
        BoundaryOptions {
            maximum_work: 16_000_001,
            ..BoundaryOptions::default()
        },
        BoundaryOptions {
            maximum_bytes: 0,
            ..BoundaryOptions::default()
        },
        BoundaryOptions {
            maximum_bytes: (1 << 20) + 1,
            ..BoundaryOptions::default()
        },
    ] {
        assert!(matches!(s.boundary(options, None), Err(Error::Invalid)));
    }
    let count = std::cell::Cell::new(0);
    s.boundary(
        BoundaryOptions::default(),
        Some(&|| {
            count.set(count.get() + 1);
            false
        }),
    )
    .unwrap();
    let total = count.get();
    for at in [1, total] {
        count.set(0);
        assert!(matches!(
            s.boundary(
                BoundaryOptions::default(),
                Some(&|| {
                    count.set(count.get() + 1);
                    count.get() == at
                })
            ),
            Err(Error::Cancelled)
        ));
    }
}
#[test]
fn actual_incidence_receipt_rejects_stale_ownership_orientation_and_pinched_endpoint_maps() {
    let s = source(false, 1, false);
    let make = || s.boundary(BoundaryOptions::default(), None).unwrap();
    let valid = make();
    let (nv, ne, nf) = (valid.vertices.len(), valid.curves.len(), valid.faces.len());
    let check = |b: &Boundary| {
        b.validate(
            nv,
            ne,
            nf,
            &mut Work {
                left: 16_000_000,
                interrupt: None,
            },
        )
    };
    let mut b = make();
    b.vertices[0] = b.vertices[1];
    assert_eq!(check(&b), Err(Error::Unresolved));
    let mut b = make();
    b.curves[0].vertices[0] = VertexKey::Plain { end: End::Start };
    b.vertices
        .retain(|v| *v != VertexKey::Plain { end: End::Start });
    assert_eq!(check(&b), Err(Error::Unresolved));
    // Keep the table counts while making one endpoint refer to a foreign start.
    let mut b = make();
    b.curves[0].vertices[0] = VertexKey::Thread {
        component: Component::Metric,
        end: End::Start,
        period: 1,
        knot: Knot::CrestEnd,
    };
    assert_eq!(check(&b), Err(Error::Unresolved));
    let mut b = make();
    let key = b.curves[0].key;
    let pos = b.coedges.iter().position(|e| e.curve != key).unwrap();
    b.coedges[pos].curve = key;
    assert_eq!(check(&b), Err(Error::Unresolved));
    let mut b = make();
    b.coedges[0].forward = !b.coedges[0].forward;
    assert_eq!(check(&b), Err(Error::Unresolved));
    let mut b = make();
    b.curves[0].vertices.swap(0, 1);
    assert_eq!(check(&b), Err(Error::Unresolved));
    // Exact metadata aliasing is refused even when every face's loop closes.
    let mut b = make();
    let v = b.vertices[0];
    for c in &mut b.curves {
        c.vertices = [v; 2];
    }
    assert_eq!(check(&b), Err(Error::Unresolved));
    let mut b = make();
    let before = b.retained_bytes();
    b.vertices.reserve(1000);
    assert_eq!(b.check_memory(before + SCRATCH), Err(Error::Limit));
}
