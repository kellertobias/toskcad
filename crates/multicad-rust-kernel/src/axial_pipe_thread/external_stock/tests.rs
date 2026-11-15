use super::super::authored::{AuthoredSection, AuthoredSpec, AuthoredThread};
use super::*;
fn source(quarter: bool) -> StockClipPlan {
    let (tpi, e1, l1, l2, od, max) = if quarter {
        (18., 0.49163, 0.2278, 0.4018, 0.540, 0.088)
    } else {
        (27., 0.37360, 0.1615, 0.2639, 0.405, 0.096)
    };
    let p = 25.4 / tpi;
    let k = 1. / 32.;
    let m = 3_f64.sqrt();
    let s = 1. / (m - k) + 1. / (m + k);
    let d = p / s - (0.033 + max) * p;
    AuthoredThread::new(
        AuthoredSpec {
            major_radius_mm: (e1 * 25.4 - 2. * k * l1 * 25.4 + d) / 2.,
            pitch_mm: p,
            length_mm: l2 * 25.4,
            radial_taper: k,
            section: AuthoredSection::NptFlat {
                depth_mm: d,
                flat_width_mm: (p - d * s) / 2.,
            },
            crest_radius_limit_mm: Some(od * 25.4 / 2.),
        },
        1e-7,
        1e-5,
    )
    .unwrap()
    .stock_clip_plan(Default::default(), None)
    .unwrap()
}
fn stock(thread: StockClipPlan, height: f64) -> crate::cylinder::AxialCylinder {
    let r = if thread.stock_radius_limit_mm() > 6. {
        Rational::new(3429, 500).unwrap()
    } else {
        Rational::new(10287, 2000).unwrap()
    };
    let zero = Rational::new(0, 1).unwrap();
    let circle = crate::Circle::new([zero; 2], r).unwrap();
    crate::cylinder::AxialCylinder::new(
        73,
        circle,
        zero,
        Rational::parse_decimal(&height.to_string()).unwrap(),
    )
    .unwrap()
}
#[test]
fn real_external_pipe_cut_keeps_original_owners_and_one_manifold_stock_tail_boundary() {
    for quarter in [false, true] {
        let thread = source(quarter);
        let length = thread.source().spec().length_mm;
        for tail in [0.25, 2., 7.] {
            let cut = ExternalStockCut::new(
                thread,
                stock(thread, length + tail),
                Default::default(),
                None,
            )
            .unwrap();
            assert_eq!(cut.thread().source().spec(), thread.source().spec());
            assert_eq!(cut.stock_length_mm(), length + tail);
            assert_eq!(cut.stock(), &stock(thread, length + tail));
            assert!(
                cut.original_source_position_transfer_mm() > 0.
                    && cut.original_source_position_transfer_mm() < 1e-12
            );
            assert!(cut.original_source_derivative_transfer_mm() < 1e-10);
            assert!(cut.radial_scale_bounds()[0] > 0.);
            assert_eq!(
                cut.stock_radius_mm(),
                cut.stock().circle().radius().as_f64()
            );
            for axis in 0..3 {
                let original = cut.stock().bounds();
                assert!(
                    exact(original[0][axis])
                        >= BigRational::from_float(cut.bounds_mm()[axis][0]).unwrap()
                );
                assert!(
                    exact(original[1][axis])
                        <= BigRational::from_float(cut.bounds_mm()[axis][1]).unwrap()
                );
            }
            let actual = exact(cut.stock().height()) - BigRational::from_float(length).unwrap();
            assert!(
                BigRational::from_float(cut.tail_length_bounds_mm()[0]).unwrap() <= actual
                    && BigRational::from_float(cut.tail_length_bounds_mm()[1]).unwrap() >= actual
            );
            assert!(cut.retained_budget_bytes() > 0);
            let boundary = cut.boundary(Default::default(), None).unwrap();
            assert_eq!(boundary.source().stock_length_mm(), cut.stock_length_mm());
            assert_eq!(boundary.vertices().len(), 18);
            assert_eq!(boundary.curves().len(), 24);
            assert_eq!(boundary.faces().len(), 9);
            assert!(
                boundary
                    .faces()
                    .iter()
                    .all(|f| f.face() != Face::Thread(StockClipFace::EndCap))
            );
            assert_eq!(
                boundary
                    .faces()
                    .iter()
                    .find(|f| f.face() == Face::Shoulder)
                    .unwrap()
                    .loops()[0]
                    .len(),
                6
            );
            assert_eq!(
                boundary
                    .faces()
                    .iter()
                    .find(|f| f.face() == Face::TailCylinder)
                    .unwrap()
                    .loops()
                    .len(),
                2
            );
            assert!(
                thread
                    .physical_curves()
                    .iter()
                    .all(|c| boundary.curves().contains(&Curve::Thread(*c)))
            );
            assert!(
                !boundary
                    .curves()
                    .contains(&Curve::ShoulderOuter(StockClipBand::Crest))
            );
            let crest = Curve::Thread(StockClipCurve::Ring {
                station: StockClipStation::End,
                band: StockClipBand::Crest,
            });
            let owners = boundary
                .faces()
                .iter()
                .filter(|f| f.loops().iter().flatten().any(|c| c.curve == crest))
                .map(FaceLoop::face)
                .collect::<Vec<_>>();
            assert_eq!(
                owners,
                vec![
                    Face::Thread(StockClipFace::StockCylinder),
                    Face::TailCylinder
                ]
            );
            assert!(matches!(
                cut.boundary(
                    StockClipOptions {
                        maximum_work: boundary.work_used() - 1,
                        ..Default::default()
                    },
                    None
                ),
                Err(Error::Limit)
            ));
            assert!(matches!(
                cut.boundary(
                    StockClipOptions {
                        maximum_bytes: 1,
                        ..Default::default()
                    },
                    None
                ),
                Err(Error::Limit)
            ));
            assert!(matches!(
                cut.boundary(Default::default(), Some(&|| true)),
                Err(Error::Cancelled)
            ));
            assert!(boundary.retained_budget_bytes() > 0);
            for &key in boundary.vertices() {
                let point = cut.vertex_bounds(key, Default::default(), None).unwrap();
                assert!(point.position_mm().iter().all(|v| v[0] <= v[1]));
                assert!(point.radius_mm()[0] > 0.);
                assert!(point.work_used() > 0 && point.retained_budget_bytes() > 0);
            }
            for &key in boundary.curves() {
                for t in [0., 0.37, 1.] {
                    let point = cut
                        .curve_bounds(key, [t, t], Default::default(), None)
                        .unwrap();
                    assert!(point.position_mm().iter().all(|v| v[0] <= v[1]));
                    if t == 0. || t == 1. {
                        let v = cut
                            .vertex_bounds(
                                key.endpoints()[usize::from(t == 1.)],
                                Default::default(),
                                None,
                            )
                            .unwrap();
                        for (a, b) in point.position_mm().iter().zip(v.position_mm()) {
                            assert!(
                                a[0] <= b[1] && b[0] <= a[1],
                                "actual curve/owned vertex enclosure {key:?} t={t}"
                            );
                        }
                    }
                }
            }
            use super::super::authored::StockClipClassification::{Inside, Outside, Undetermined};
            for (point, expected) in [
                ([0., 0., length], Inside),
                ([0., 0., length + tail / 2.], Inside),
                (
                    [cut.stock_radius_mm() + 0.1, 0., length + tail / 2.],
                    Outside,
                ),
                ([0., 0., cut.stock_length_mm() + 1.], Outside),
                ([0., 0., cut.stock_length_mm()], Undetermined),
            ] {
                let q = cut
                    .classify_box(point.map(|v| [v, v]), Default::default(), None)
                    .unwrap();
                assert_eq!(q.classification(), expected);
                assert!(q.work_used() > 0 && q.retained_budget_bytes() > 0);
            }
            for z in [length * 0.17, length * 0.63] {
                for radius in [
                    cut.stock_radius_mm() * 0.88,
                    cut.stock_radius_mm() * 0.97,
                    cut.stock_radius_mm() * 1.01,
                ] {
                    for angle in [0.13_f64, 1.72, 4.51] {
                        let point = [radius * angle.cos(), radius * angle.sin(), z].map(|v| [v, v]);
                        let old = thread
                            .classify_box(point, Default::default(), None)
                            .unwrap()
                            .classification();
                        let current = cut
                            .classify_box(point, Default::default(), None)
                            .unwrap()
                            .classification();
                        assert_eq!(current, old);
                    }
                }
            }
            assert!(matches!(
                cut.vertex_bounds(
                    Vertex::ShoulderRoot(StockClipRail::LowerCrest),
                    Default::default(),
                    None
                ),
                Err(Error::Invalid)
            ));
            assert!(matches!(
                cut.curve_bounds(
                    Curve::ShoulderOuter(StockClipBand::Crest),
                    [0., 1.],
                    Default::default(),
                    None
                ),
                Err(Error::Invalid)
            ));
            assert!(matches!(
                cut.curve_bounds(
                    Curve::Thread(StockClipCurve::Ring {
                        station: StockClipStation::Contact,
                        band: StockClipBand::Root
                    }),
                    [0., 1.],
                    Default::default(),
                    None
                ),
                Err(Error::Invalid)
            ));
            assert!(matches!(
                cut.vertex_bounds(
                    boundary.vertices()[0],
                    StockClipOptions {
                        maximum_work: cut.work_used() + 128,
                        ..Default::default()
                    },
                    None
                ),
                Err(Error::Limit)
            ));
            assert!(matches!(
                cut.classify_box(
                    [[f64::NAN, 0.], [0., 0.], [0., 0.]],
                    Default::default(),
                    None
                ),
                Err(Error::Invalid)
            ));
            assert!(matches!(
                cut.classify_box(
                    [[0., 0.], [0., 0.], [length, length]],
                    Default::default(),
                    Some(&|| true)
                ),
                Err(Error::Cancelled)
            ));
            let volume = cut.volume_bounds(1e-6, Default::default(), None).unwrap();
            let old = thread
                .volume_bounds(1e-7, Default::default(), None)
                .unwrap()
                .source_volume_mm3();
            let exact_radius = exact(cut.stock().circle().radius());
            let radial =
                &exact_radius / BigRational::from_float(thread.stock_radius_limit_mm()).unwrap();
            let exact_tail = exact(cut.stock().height()) - BigRational::from_float(length).unwrap();
            let pi_lo = BigRational::from_float(std::f64::consts::PI.next_down()).unwrap();
            let pi_hi = BigRational::from_float(std::f64::consts::PI.next_up()).unwrap();
            let lower = BigRational::from_float(old[0]).unwrap() * &radial * &radial
                + pi_lo * &exact_radius * &exact_radius * &exact_tail;
            let upper = BigRational::from_float(old[1]).unwrap() * &radial * &radial
                + pi_hi * &exact_radius * &exact_radius * &exact_tail;
            assert!(BigRational::from_float(volume.source_volume_mm3()[0]).unwrap() <= lower);
            assert!(BigRational::from_float(volume.source_volume_mm3()[1]).unwrap() >= upper);
            for face in boundary.faces() {
                let m = cut
                    .face_measure(face.face(), Default::default(), None)
                    .unwrap();
                assert_eq!(m.face(), face.face());
                assert!(m.area_mm2()[0] > 0.);
                assert!(m.work_used() > 0 && m.retained_budget_bytes() > 0);
                assert!(m.first_moments_mm3().iter().all(|v| v[0] <= v[1]));
                assert!(m.centroid_mm().iter().all(|v| v[0] <= v[1]));
            }
            assert!(matches!(
                cut.volume_bounds(
                    1e-6,
                    StockClipOptions {
                        maximum_work: volume.work_used() - 1,
                        ..Default::default()
                    },
                    None
                ),
                Err(Error::Limit)
            ));
            assert!(volume.retained_budget_bytes() > 0);
            eprintln!(
                "NPT external tail source quarter={quarter} tail={tail} vertices={} edges={} faces={} work={} bytes={}",
                boundary.vertices().len(),
                boundary.curves().len(),
                boundary.faces().len(),
                boundary.work_used(),
                boundary.retained_budget_bytes()
            );
        }
        for invalid in [0.25, length - 1.] {
            assert!(matches!(
                ExternalStockCut::new(thread, stock(thread, invalid), Default::default(), None),
                Err(Error::Invalid)
            ));
        }
        let cut =
            ExternalStockCut::new(thread, stock(thread, length + 2.), Default::default(), None)
                .unwrap();
        assert!(matches!(
            ExternalStockCut::new(
                thread,
                stock(thread, length + 2.),
                StockClipOptions {
                    maximum_work: cut.work_used() - 1,
                    ..Default::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
        assert!(matches!(
            ExternalStockCut::new(
                thread,
                stock(thread, length + 2.),
                StockClipOptions {
                    maximum_bytes: 1,
                    ..Default::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
        assert!(matches!(
            ExternalStockCut::new(
                thread,
                stock(thread, length + 2.),
                Default::default(),
                Some(&|| true)
            ),
            Err(Error::Cancelled)
        ));
    }
}
#[test]
fn actual_tail_mesh_replaces_the_internal_cap_and_has_one_connected_indexed_shell() {
    for quarter in [false, true] {
        let thread = source(quarter);
        let cut = ExternalStockCut::new(
            thread,
            stock(thread, thread.source().spec().length_mm + 2.),
            Default::default(),
            None,
        )
        .unwrap();
        let q = crate::display::Tessellation {
            deflection: 0.02,
            angle: 0.2,
            triangle_limit: 1_000_000,
        };
        let started = std::time::Instant::now();
        let candidate = cut
            .mesh_candidate(q, 512_000_000, 64 * 1024 * 1024, None)
            .unwrap();
        assert_eq!(candidate.source().stock_length_mm(), cut.stock_length_mm());
        assert!(candidate.connectivity().signed_volume_mm3 > 0.);
        assert!(candidate.source_distance_mm() <= q.deflection);
        assert!(candidate.normal_chord() <= 2. * (q.angle / 2.).sin());
        assert!(candidate.mesh().face_indices.iter().all(|v| *v < 9));
        for face in 0..9 {
            assert!(candidate.mesh().face_indices.contains(&face));
        }
        let h = candidate
            .mesh()
            .vertices
            .as_chunks::<3>()
            .0
            .iter()
            .map(|v| v[2])
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(h <= cut.stock_top_bounds_mm()[1] && h >= cut.stock_top_bounds_mm()[0]);
        let volume = cut
            .volume_bounds(1e-6, Default::default(), None)
            .unwrap()
            .source_volume_mm3();
        eprintln!(
            "NPT tail literal quarter={quarter} seconds={} vertices={} triangles={} work={} bytes={} chord={} normal={} literalvolume={} sourcevolume={volume:?}",
            started.elapsed().as_secs_f64(),
            candidate.mesh().vertices.len() / 3,
            candidate.mesh().indices.len() / 3,
            candidate.work_used(),
            candidate.retained_budget_bytes(),
            candidate.source_distance_mm(),
            candidate.normal_chord(),
            candidate.connectivity().signed_volume_mm3
        );
        assert!(matches!(
            cut.mesh_candidate(q, candidate.work_used() - 1, 64 * 1024 * 1024, None),
            Err(Error::Limit)
        ));
        assert!(matches!(
            cut.mesh_candidate(q, 512_000_000, candidate.retained_budget_bytes() - 1, None),
            Err(Error::Limit)
        ));
        assert!(matches!(
            cut.mesh_candidate(q, 512_000_000, 64 * 1024 * 1024, Some(&|| true)),
            Err(Error::Cancelled)
        ));
        let callbacks = std::cell::Cell::new(0);
        assert!(matches!(
            cut.mesh_candidate(
                q,
                512_000_000,
                64 * 1024 * 1024,
                Some(&|| {
                    let n = callbacks.get() + 1;
                    callbacks.set(n);
                    n >= 50000
                })
            ),
            Err(Error::Cancelled)
        ));
        assert!(callbacks.get() >= 50000);
        let fine = crate::display::Tessellation {
            deflection: 0.005,
            angle: 0.1,
            triangle_limit: 1_000_000,
        };
        let finer = cut
            .mesh_candidate(fine, 512_000_000, 64 * 1024 * 1024, None)
            .unwrap();
        assert!(finer.source_distance_mm() <= fine.deflection);
        assert!(finer.mesh().indices.len() > candidate.mesh().indices.len());
        eprintln!(
            "NPT exact-tail finer quarter={quarter} triangles={} work={} bytes={}",
            finer.mesh().indices.len() / 3,
            finer.work_used(),
            finer.retained_budget_bytes()
        );
        let actual = candidate.connectivity().signed_volume_mm3;
        assert!((actual - volume[0]).abs() < q.deflection * 500.);
    }
}
#[test]
fn exact_stock_datums_and_same_float_aliases_remain_distinct_under_the_declared_map() {
    let thread = source(false);
    let radius = Rational::new(10287, 2000).unwrap();
    let alias = Rational::new(514350000000000001, 100000000000000000).unwrap();
    assert_eq!(radius.as_f64(), alias.as_f64());
    assert_ne!(radius, alias);
    let center = [Rational::new(1, 10).unwrap(), Rational::new(-2, 5).unwrap()];
    let bottom = Rational::new(1, 7).unwrap();
    let height = Rational::new(11, 1).unwrap();
    let a = crate::cylinder::AxialCylinder::new(
        389,
        crate::Circle::new(center, radius).unwrap(),
        bottom,
        height,
    )
    .unwrap();
    let b = crate::cylinder::AxialCylinder::new(
        389,
        crate::Circle::new(center, alias).unwrap(),
        bottom,
        height,
    )
    .unwrap();
    let cut = ExternalStockCut::new(thread, a, Default::default(), None).unwrap();
    let other = ExternalStockCut::new(thread, b, Default::default(), None).unwrap();
    assert_eq!(*cut.stock(), a);
    let axis = cut.stock_axis_base_bounds_mm();
    for (value, bounds) in [center[0], center[1], bottom].into_iter().zip(axis) {
        assert!(
            BigRational::from_float(bounds[0]).unwrap() <= exact(value)
                && BigRational::from_float(bounds[1]).unwrap() >= exact(value)
        );
    }
    assert_eq!(*other.stock(), b);
    assert_ne!(cut.stock(), other.stock());
    assert_eq!(cut.thread().source().spec(), thread.source().spec());
    assert_ne!(
        cut.original_source_position_transfer_mm(),
        other.original_source_position_transfer_mm()
    );
    let original = BigRational::from_float(thread.stock_radius_limit_mm()).unwrap();
    let transfer = (exact(radius) - &original).abs();
    assert!(
        BigRational::from_float(cut.original_source_position_transfer_mm()).unwrap() >= transfer
    );
    assert!(
        BigRational::from_float(cut.original_source_position_transfer_mm().next_down()).unwrap()
            < transfer
    );
    let join = exact(bottom) + BigRational::from_float(thread.source().spec().length_mm).unwrap();
    let tail = exact(height) - BigRational::from_float(thread.source().spec().length_mm).unwrap();
    for (value, bounds) in [
        (exact(a.bottom()), cut.stock_bottom_bounds_mm()),
        (exact(a.top()), cut.stock_top_bounds_mm()),
        (join, cut.thread_join_bounds_mm()),
        (tail, cut.tail_length_bounds_mm()),
    ] {
        assert!(
            BigRational::from_float(bounds[0]).unwrap() <= value
                && BigRational::from_float(bounds[1]).unwrap() >= value
        );
    }
    let major = BigRational::from_float(thread.source().spec().major_radius_mm).unwrap()
        * exact(radius)
        / original;
    let bounds = cut.derived_major_radius_bounds_mm();
    assert!(
        BigRational::from_float(bounds[0]).unwrap() <= major
            && BigRational::from_float(bounds[1]).unwrap() >= major
    );
    let boundary = cut.boundary(Default::default(), None).unwrap();
    for &curve in boundary.curves() {
        for t in [0., 0.31, 1.] {
            let point = cut
                .curve_bounds(curve, [t, t], Default::default(), None)
                .unwrap();
            assert!(point.position_mm().iter().all(|v| v[0] <= v[1]));
        }
    }
    for face in boundary.faces() {
        cut.face_measure(face.face(), Default::default(), None)
            .unwrap();
    }
    let q = crate::display::Tessellation {
        deflection: 0.02,
        angle: 0.2,
        triangle_limit: 1_000_000,
    };
    let mesh = cut
        .mesh_candidate(q, 512_000_000, 64 * 1024 * 1024, None)
        .unwrap();
    assert_eq!(*mesh.source().stock(), a);
    assert!(mesh.source_distance_mm() <= q.deflection);
    assert!(mesh.connectivity().signed_volume_mm3 > 0.);
    let r = interval_for_test(radius);
    let cx = interval_for_test(center[0]);
    let cy = interval_for_test(center[1]);
    for point in mesh.mesh().vertices.as_chunks::<3>().0.iter() {
        let radial =
            ((I::point(point[0]) - cx).square() + (I::point(point[1]) - cy).square()).sqrt();
        assert!(radial.hi <= r.lo);
        assert!(BigRational::from_float(point[2]).unwrap() >= exact(bottom));
        assert!(BigRational::from_float(point[2]).unwrap() <= exact(a.top()));
    }
    eprintln!(
        "exact-stock derived NPT transfer={} d1={} vertices={} triangles={}",
        cut.original_source_position_transfer_mm(),
        cut.original_source_derivative_transfer_mm(),
        mesh.mesh().vertices.len() / 3,
        mesh.mesh().indices.len() / 3
    );
}
fn interval_for_test(value: Rational) -> I {
    enclosure(&exact(value)).unwrap()
}
