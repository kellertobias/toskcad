use super::*;
fn bits(uv: [f64; 2]) -> [u64; 2] {
    uv.map(f64::to_bits)
}
fn check(surface: &Surface) {
    let prepared = PreparedInverse::new(surface);
    assert!(prepared.grid.get().is_none());
    if let Surface::Nurbs { .. } = surface {
        let old = surface.candidates([0.; 3]);
        let cached = prepared.candidates().unwrap();
        assert_eq!(old.len(), cached.len());
        for (uv, cached) in old.into_iter().zip(cached) {
            assert_eq!(bits(uv), bits(cached.uv));
            assert_eq!(
                surface.evaluate(uv).0.map(f64::to_bits),
                cached.position.map(f64::to_bits)
            );
        }
        assert!(
            std::mem::size_of_val(&prepared)
                + prepared.grid.get().unwrap().as_ref().unwrap().capacity()
                    * std::mem::size_of::<Candidate>()
                <= MAX_BYTES
        );
    }
    for p in [
        [0., 0., 0.],
        [0.37, -0.42, 0.11],
        [20., -3., 9.],
        [f64::NAN, 0., 1.],
        [f64::INFINITY, 1., 0.],
    ] {
        for seed in [None, Some([0., 0.]), Some([0.74, -0.3]), Some([10., 20.])] {
            assert_eq!(
                bits(prepared.invert(p, seed)),
                bits(surface.invert(p, seed)),
                "p={p:?} seed={seed:?}"
            );
        }
    }
    let disabled = PreparedInverse {
        surface,
        grid: OnceCell::new(),
        maximum_bytes: 0,
    };
    assert_eq!(
        bits(disabled.invert([0.3, 0.4, 0.1], None)),
        bits(surface.invert([0.3, 0.4, 0.1], None))
    );
    assert!(disabled.candidates().is_none());
}
#[test]
fn actual_source_spline_surfaces_and_proper_poses_keep_bitwise_inverses() {
    let text = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../multicad-app/tests/fixtures/exchange/corpus/occt-spline-prism.step"
    ));
    let solids = crate::step_read::brep::read_solids(text, 1e-7, None).unwrap();
    for face in &solids[0].boundary().faces {
        check(&face.surface);
    }
    // Real public producer/decoder, without full Solid mass integration.
    let tool = crate::helical_thread::Thread::new(crate::helical_thread::Spec {
        major_radius_mm: 1.,
        pitch_mm: 0.5,
        length_mm: 1.,
        depth_mm: 0.15,
        crest_width_mm: 0.1,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: false,
        root_radius_mm: None,
    })
    .unwrap();
    let written =
        crate::step::write_helical_thread(&tool, None, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7)
            .unwrap();
    let boundaries = crate::step_read::brep::read_boundaries(&written, 1e-7).unwrap();
    let mut count = 0;
    for face in &boundaries[0].faces {
        if !matches!(face.surface, Surface::Nurbs { .. }) {
            continue;
        }
        check(&face.surface);
        for shift in [[1., -2., 3.], [1e3, -1e3, 250.]] {
            let linear = |p: Point| [p[2], p[0], p[1]];
            let map = |p: Point| add(linear(p), shift);
            check(&face.surface.mapped(&map, &linear));
        }
        count += 1;
    }
    assert_eq!(count, 4);
}
#[test]
fn rational_unclamped_repeated_knots_and_closed_grids_keep_exact_proposals() {
    let points = (0..4)
        .flat_map(|u| (0..3).map(move |v| [u as f64, v as f64, 0.1 * (u * v) as f64]))
        .collect::<Vec<_>>();
    let weights = Some((0..12).map(|k| 1. + 0.125 * (k % 3) as f64).collect());
    let s = NurbsSurface::new(
        [2, 2],
        [4, 3],
        points,
        weights,
        [vec![-1., 0., 0.4, 1., 2.], vec![0., 1.]],
        [vec![1, 2, 1, 2, 1], vec![3, 3]],
    )
    .unwrap();
    for closed in [[false, false], [true, false], [true, true]] {
        check(&Surface::Nurbs {
            surface: s.clone(),
            closed,
        });
    }
}
#[test]
fn analytic_and_query_dependent_swept_surfaces_never_cache_point_dependent_seeds() {
    let frame = Frame {
        origin: [0.; 3],
        x: [1., 0., 0.],
        y: [0., 1., 0.],
        z: [0., 0., 1.],
    };
    let curve = Curve::Circle { frame, radius: 2. };
    for surface in [
        Surface::Plane(frame),
        Surface::Extrusion {
            curve: curve.clone(),
            direction: [0., 0., 1.],
        },
        Surface::Revolution {
            curve,
            origin: [0.; 3],
            axis: [0., 1., 0.],
        },
    ] {
        check(&surface);
        assert!(PreparedInverse::new(&surface).candidates().is_none());
    }
}
