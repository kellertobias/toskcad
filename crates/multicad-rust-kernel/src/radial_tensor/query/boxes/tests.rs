use super::*;
use std::cell::Cell;
fn at(point: Point, width: f64) -> [[f64; 2]; 3] {
    point.map(|v| [v - width / 2., v + width / 2.])
}
#[test]
fn analytic_cone_sector_boxes_keep_the_actual_surface_and_finite_collar() {
    let input = crate::radial_tensor::tests::five_bands::graph(
        &[
            [3., -4.],
            [6., -3.],
            [7., -1.],
            [7., 1.],
            [6., 3.],
            [3., 4.],
        ],
        true,
    );
    let body = BoundaryOwnership::new(input, 1e-7, None)
        .unwrap()
        .bind_surfaces(1e-7, None)
        .unwrap()
        .certify_axial(1e-7, 1e-6, None)
        .unwrap()
        .certify_caps(None)
        .unwrap()
        .certify_ordering(None)
        .unwrap()
        .certify_material(1e-7, None)
        .unwrap();
    for (point, location) in [
        ([-5.2, 0.5, 2.], Location::Outside),
        ([-4.8, 0.5, 2.], Location::Inside),
        ([-5., 0., 2.], Location::Unresolved),
    ] {
        assert_eq!(
            body.classify_box(at(point, 1e-6), Default::default(), None)
                .unwrap()
                .location(),
            location
        );
    }
}
#[test]
fn actual_retained_sector_and_core_boxes_preserve_all_coordinates_and_source_identity() {
    let body = crate::radial_tensor::tests::material();
    for (point, location) in [
        ([0.1, 0.2, 4.], Location::Inside),
        ([0.1, 0.2, -1.], Location::Outside),
        ([0., 0., 0.], Location::Unresolved),
    ] {
        let bounds = at(point, 1e-6);
        let q = body.classify_box(bounds, Default::default(), None).unwrap();
        assert_eq!(q.location(), location);
        assert!(std::ptr::eq(q.source(), &body));
        assert_eq!(q.bounds_mm(), bounds);
        assert_eq!(q.work_used() + q.remaining_work(), MAX_WORK);
    }
    // Independently read an ACTUAL retained tensor point. Its open sector
    // has a unique radial crossing, so these finite offsets are opposite
    // material sides. Neither a header nor a nominal thread supplies it.
    let surface = &body.ordered.caps.axial.tensors[1].surface;
    let p = surface
        .point_bounds([0, 0], [[0.37; 2], [0.41; 2]])
        .unwrap()
        .map(|v| v[0].midpoint(v[1]));
    let radius = p[0].hypot(p[1]);
    for (offset, location) in [(-0.01, Location::Inside), (0.01, Location::Outside)] {
        let point = [
            p[0] * (radius + offset) / radius,
            p[1] * (radius + offset) / radius,
            p[2],
        ];
        let bounds = at(point, 1e-6);
        let q = body.classify_box(bounds, Default::default(), None).unwrap();
        assert_eq!(q.location(), location, "{point:?}");
        let cost = q.work_used();
        assert_eq!(
            body.classify_box(bounds, BoxQueryOptions { maximum_work: cost }, None)
                .unwrap()
                .remaining_work(),
            0
        );
        assert_eq!(
            body.classify_box(
                bounds,
                BoxQueryOptions {
                    maximum_work: cost - 1
                },
                None
            )
            .err()
            .unwrap(),
            Error::Limit
        );
        let calls = Cell::new(0);
        body.classify_box(
            bounds,
            Default::default(),
            Some(&|| {
                calls.set(calls.get() + 1);
                false
            }),
        )
        .unwrap();
        for stop in [1, calls.get() / 2, calls.get() - 1] {
            let n = Cell::new(0);
            assert_eq!(
                body.classify_box(
                    bounds,
                    Default::default(),
                    Some(&|| {
                        let at = n.get();
                        n.set(at + 1);
                        at == stop
                    })
                )
                .err()
                .unwrap(),
                Error::Cancelled
            );
        }
        for corner in 0..8 {
            let point = std::array::from_fn(|i| bounds[i][usize::from(corner & (1 << i) != 0)]);
            assert_eq!(body.classify(point, None).unwrap(), location);
        }
    }
}
#[test]
fn whole_box_query_uses_one_exact_balance_and_cancels_before_publication() {
    let body = crate::radial_tensor::tests::material();
    let bounds = at([0.1, 0.2, 4.], 1e-6);
    let calls = Cell::new(0);
    let q = body
        .classify_box(
            bounds,
            Default::default(),
            Some(&|| {
                calls.set(calls.get() + 1);
                false
            }),
        )
        .unwrap();
    let exact = q.work_used();
    assert_eq!(
        body.classify_box(
            bounds,
            BoxQueryOptions {
                maximum_work: exact
            },
            None
        )
        .unwrap()
        .remaining_work(),
        0
    );
    for budget in [0, 1, exact - 1] {
        assert_eq!(
            body.classify_box(
                bounds,
                BoxQueryOptions {
                    maximum_work: budget
                },
                None
            )
            .err()
            .unwrap(),
            Error::Limit
        );
    }
    for stop in 0..calls.get() {
        let n = Cell::new(0);
        assert_eq!(
            body.classify_box(
                bounds,
                Default::default(),
                Some(&|| {
                    let at = n.get();
                    n.set(at + 1);
                    at == stop
                })
            )
            .err()
            .unwrap(),
            Error::Cancelled
        );
    }
    assert_eq!(
        body.classify_box(bounds, Default::default(), None)
            .unwrap()
            .location(),
        Location::Inside
    );
    for bounds in [
        [[1., 0.]; 3],
        [[f64::NAN, 0.]; 3],
        [[0., f64::INFINITY]; 3],
        [[0., 1e7_f64.next_up()]; 3],
    ] {
        assert_eq!(
            body.classify_box(bounds, Default::default(), None)
                .err()
                .unwrap(),
            Error::Invalid
        );
    }
    assert_eq!(
        body.classify_box(
            bounds,
            BoxQueryOptions {
                maximum_work: MAX_WORK + 1
            },
            None
        )
        .err()
        .unwrap(),
        Error::Invalid
    );
}

#[test]
fn source_pose_and_hand_are_retained_in_world_box_classification() {
    for left_hand in [false, true] {
        let thread = crate::tapered_thread::Thread::new(
            crate::helical_thread::Spec {
                major_radius_mm: 4.,
                pitch_mm: 1.25,
                length_mm: 8.,
                depth_mm: 0.5,
                crest_width_mm: 0.25,
                lower_flank_deg: 12.,
                upper_flank_deg: 18.,
                starts: 1,
                left_hand,
                root_radius_mm: None,
            },
            0.03125,
            1e-7,
            1e-6,
        )
        .unwrap();
        let origin = [11., -7., 5.];
        let u = [0.8, 0.6, 0.];
        let v = [-0.36, 0.48, 0.8];
        let n = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        let text = crate::step::write_tapered_thread(&thread, origin, u, v, 1e-7).unwrap();
        let body = crate::step_read::read_radial_boundary(&text, 1e-7)
            .unwrap()
            .bind_surfaces(1e-7, None)
            .unwrap()
            .certify_axial(1e-7, 1e-6, None)
            .unwrap()
            .certify_caps(None)
            .unwrap()
            .certify_ordering(None)
            .unwrap()
            .certify_material(1e-7, None)
            .unwrap();
        for (local, expected) in [
            ([0.1, 0.2, 4.], Location::Inside),
            ([0.1, 0.2, -1.], Location::Outside),
            ([0., 0., 0.], Location::Unresolved),
        ] {
            let point = std::array::from_fn(|i| {
                origin[i] + u[i] * local[0] + v[i] * local[1] + n[i] * local[2]
            });
            assert_eq!(
                body.classify_box(at(point, 1e-6), Default::default(), None)
                    .unwrap()
                    .location(),
                expected
            );
        }
    }
}
