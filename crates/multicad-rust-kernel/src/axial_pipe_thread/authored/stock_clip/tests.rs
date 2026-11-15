use super::*;
#[test]
fn canonical_clipped_source_keeps_actual_law_event_and_positive_volume() {
    let mut spec = super::super::tests::source(1. / 32.).spec();
    spec.crest_radius_limit_mm = Some(0.405 * 25.4 / 2.);
    let source = AuthoredThread::new(spec, 1e-7, 1e-5).unwrap();
    let plan = source
        .stock_clip_plan(StockClipOptions::default(), None)
        .unwrap();
    assert_eq!(*plan.source(), source);
    assert_eq!(
        plan.stock_radius_limit_mm(),
        spec.crest_radius_limit_mm.unwrap()
    );
    let event = plan.contact_height_bounds_mm();
    assert!(event[0] <= 5.77370606837215 && event[1] >= 5.773706068372149);
    for patch in plan.patches() {
        let sample = plan
            .sample_bounds(
                patch,
                [0.3, 0.3],
                [0.4, 0.4],
                StockClipOptions::default(),
                None,
            )
            .unwrap();
        assert!(sample.point().radius_mm()[1] <= plan.stock_radius_limit_mm());
        assert!(sample.normal().iter().all(|n| n[0] >= -1. && n[1] <= 1.));
    }
    let volume = plan
        .volume_bounds(1e-6, StockClipOptions::default(), None)
        .unwrap()
        .source_volume_mm3();
    assert!(volume[0] > 0. && volume[1] - volume[0] < 1e-6);
    eprintln!("NPT SOURCE VOLUME {volume:?}");
    let p = plan
        .point_bounds(
            plan.patches()[0],
            [0.3, 0.3],
            [0.4, 0.4],
            Default::default(),
            None,
        )
        .unwrap();
    let n = plan
        .sample_bounds(
            plan.patches()[0],
            [0.3, 0.3],
            [0.4, 0.4],
            Default::default(),
            None,
        )
        .unwrap();
    let v = plan.volume_bounds(1e-6, Default::default(), None).unwrap();
    eprintln!(
        "NPT SOURCE COSTS plan={}/{}B point={}/{}B normal={}/{}B volume={}/{}B",
        plan.work_used(),
        plan.retained_budget_bytes(),
        p.work_used(),
        p.retained_budget_bytes(),
        n.point().work_used(),
        n.point().retained_budget_bytes(),
        v.work_used(),
        v.retained_budget_bytes()
    );
}

pub(in crate::axial_pipe_thread::authored) fn canonical(
    quarter: bool,
    scale: f64,
) -> AuthoredThread {
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
            major_radius_mm: scale * (e1 * 25.4 - 2. * k * l1 * 25.4 + d) / 2.,
            pitch_mm: scale * p,
            length_mm: scale * l2 * 25.4,
            radial_taper: k,
            section: AuthoredSection::NptFlat {
                depth_mm: scale * d,
                flat_width_mm: scale * (p - d * s) / 2.,
            },
            crest_radius_limit_mm: Some(scale * od * 25.4 / 2.),
        },
        1e-7,
        1e-5,
    )
    .unwrap()
}
fn midpoint(v: [f64; 2]) -> f64 {
    v[0] + (v[1] - v[0]) * 0.5
}
fn length(v: [f64; 3]) -> f64 {
    v[0].hypot(v[1]).hypot(v[2])
}
#[test]
fn actual_phase_source_and_independent_differential_normals_match_each_chart() {
    for quarter in [false, true] {
        let source = canonical(quarter, 1.);
        let plan = source.stock_clip_plan(Default::default(), None).unwrap();
        for patch in plan.patches() {
            for (u, v) in [(0.2, 0.3), (0.7, 0.8)] {
                let sample = plan
                    .sample_bounds(patch, [u, u], [v, v], Default::default(), None)
                    .unwrap();
                let point = sample.point();
                let q = midpoint(point.phase_mm()).rem_euclid(source.spec.pitch_mm);
                let z = midpoint(point.axial_mm());
                let radius = source.radial_bounds([q, q], [z, z]).unwrap();
                assert!((midpoint(radius) - midpoint(point.radius_mm())).abs() < 1e-11);
                let h = 1e-6;
                let xyz = |a, b| {
                    plan.point_bounds(patch, [a, a], [b, b], Default::default(), None)
                        .unwrap()
                        .position_mm()
                        .map(midpoint)
                };
                let a = xyz(u - h, v);
                let b = xyz(u + h, v);
                let c = xyz(u, v - h);
                let d = xyz(u, v + h);
                let du = std::array::from_fn::<_, 3, _>(|i| (b[i] - a[i]) / (2. * h));
                let dv = std::array::from_fn::<_, 3, _>(|i| (d[i] - c[i]) / (2. * h));
                let n = [
                    dv[1] * du[2] - dv[2] * du[1],
                    dv[2] * du[0] - dv[0] * du[2],
                    dv[0] * du[1] - dv[1] * du[0],
                ];
                let n = n.map(|v| v / length(n));
                let expected = sample.normal().map(midpoint);
                assert!(length(std::array::from_fn(|i| n[i] - expected[i])) < 1e-6);
                let p = point.position_mm().map(midpoint);
                assert!(p[0] * expected[0] + p[1] * expected[1] > 0.);
            }
        }
    }
}
#[test]
fn owned_event_edges_close_one_oriented_source_atlas_with_manifold_links() {
    use std::collections::{BTreeMap, BTreeSet};
    let source = canonical(false, 1.);
    let plan = source.stock_clip_plan(Default::default(), None).unwrap();
    let vertices = plan.vertices();
    let curves = plan.curves();
    let patches = plan.patches();
    assert_eq!(vertices.len() + patches.len() + 2, curves.len() + 2);
    let mut uses = BTreeMap::<StockClipCurve, Vec<(usize, bool)>>::new();
    let mut wires = Vec::new();
    for patch in patches {
        wires.push(patch.coedges());
    }
    wires.push(plan.cap_coedges(StockClipStation::Start).unwrap());
    wires.push(plan.cap_coedges(StockClipStation::End).unwrap());
    for (face, wire) in wires.iter().enumerate() {
        let endpoints = wire.map(|c| {
            let e = c.curve.endpoints();
            if c.forward { e } else { [e[1], e[0]] }
        });
        for i in 0..4 {
            assert_eq!(endpoints[i][1], endpoints[(i + 1) % 4][0]);
        }
        for c in wire {
            uses.entry(c.curve).or_default().push((face, c.forward));
        }
    }
    assert_eq!(uses.len(), curves.len());
    let mut reached = BTreeSet::new();
    let mut todo = vec![0_usize];
    while let Some(face) = todo.pop() {
        if reached.insert(face) {
            for edge in wires[face] {
                todo.extend(uses[&edge.curve].iter().map(|u| u.0));
            }
        }
    }
    assert_eq!(reached.len(), wires.len());

    for u in uses.values() {
        assert_eq!(u.len(), 2);
        assert_ne!(u[0].1, u[1].1);
    }
    for vertex in vertices {
        let incident: Vec<_> = curves
            .iter()
            .filter(|c| c.endpoints().contains(&vertex))
            .collect();
        let mut adjacency = BTreeMap::<usize, BTreeSet<usize>>::new();
        for e in incident {
            let u = &uses[e];
            adjacency.entry(u[0].0).or_default().insert(u[1].0);
            adjacency.entry(u[1].0).or_default().insert(u[0].0);
        }
        assert!(adjacency.values().all(|neighbors| neighbors.len() == 2));
        let mut reached = BTreeSet::new();
        let mut todo = vec![*adjacency.keys().next().unwrap()];
        while let Some(i) = todo.pop() {
            if reached.insert(i) {
                todo.extend(adjacency[&i].iter().copied());
            }
        }
        assert_eq!(reached.len(), adjacency.len());
    }
    for curve in curves {
        let ends = curve.endpoints();
        for (t, key) in [(0., ends[0]), (1., ends[1])] {
            let a = plan
                .curve_bounds(curve, [t, t], Default::default(), None)
                .unwrap()
                .position_mm();
            let b = plan
                .vertex_bounds(key, Default::default(), None)
                .unwrap()
                .position_mm();
            assert_eq!(a, b, "{curve:?} {key:?}");
        }
        let owners = plan.curve_faces(curve);
        if let StockClipCurve::Ring {
            station: StockClipStation::Contact,
            band,
        } = curve
        {
            assert_eq!(owners[0] != owners[1], band == StockClipBand::Crest);
        }
    }
    for band in BANDS {
        let a = plan
            .point_bounds(
                StockClipPatch {
                    slab: StockClipSlab::Before,
                    band,
                },
                [0.3, 0.3],
                [1., 1.],
                Default::default(),
                None,
            )
            .unwrap();
        let b = plan
            .point_bounds(
                StockClipPatch {
                    slab: StockClipSlab::After,
                    band,
                },
                [0.3, 0.3],
                [0., 0.],
                Default::default(),
                None,
            )
            .unwrap();
        assert_eq!(a.position_mm(), b.position_mm());
        let a = plan
            .sample_bounds(
                StockClipPatch {
                    slab: StockClipSlab::Before,
                    band,
                },
                [0.3, 0.3],
                [1., 1.],
                Default::default(),
                None,
            )
            .unwrap();
        let b = plan
            .sample_bounds(
                StockClipPatch {
                    slab: StockClipSlab::After,
                    band,
                },
                [0.3, 0.3],
                [0., 0.],
                Default::default(),
                None,
            )
            .unwrap();
        if band == StockClipBand::Crest {
            assert_ne!(a.normal(), b.normal());
        } else {
            assert_eq!(a.normal(), b.normal());
        }
    }
    let key = StockClipCurve::Ring {
        station: StockClipStation::Contact,
        band: StockClipBand::Root,
    };
    let a = plan
        .curve_bounds(key, [0., 0.], Default::default(), None)
        .unwrap();
    let b = plan
        .curve_bounds(key, [1., 1.], Default::default(), None)
        .unwrap();
    assert!(b.phase_mm()[0] > a.phase_mm()[1]);
    assert_eq!(
        b.position_mm(),
        plan.vertex_bounds(key.endpoints()[1], Default::default(), None)
            .unwrap()
            .position_mm()
    );
    assert_eq!(*plan.source(), source);
}

fn gauss2(f: impl Fn(f64) -> f64, a: f64, b: f64) -> f64 {
    let m = a + (b - a) * 0.5;
    let h = (b - a) * 0.5;
    let d = h / 3_f64.sqrt();
    h * (f(m - d) + f(m + d))
}
fn independent_volume(source: AuthoredThread) -> f64 {
    let spec = source.spec();
    let limit = spec.crest_radius_limit_mm.unwrap();
    // Independently bisect the ORIGINAL unclipped crest against stock; do not
    // use the plan event/contacts or its positive coefficient volume formula.
    let mut lo = 0.;
    let mut hi = spec.length_mm;
    for _ in 0..64 {
        let mid = lo + (hi - lo) * 0.5;
        if mid == lo || mid == hi {
            break;
        }
        if spec.major_radius_mm + spec.radial_taper * mid < limit {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let event = lo + (hi - lo) * 0.5;
    let section = |z: f64| {
        let radial = |q: f64| {
            midpoint(
                source
                    .radial_bounds([q.rem_euclid(spec.pitch_mm); 2], [z; 2])
                    .unwrap(),
            )
            .powi(2)
        };
        let mut total = 0.;
        for band in 0..3 {
            let (a, ra) = source.endpoint(band, false).unwrap();
            let (b, rb) = source.endpoint(band, true).unwrap();
            let (a, b, ra, rb) = (
                a.midpoint(),
                b.midpoint(),
                ra.midpoint() + spec.radial_taper * z,
                rb.midpoint() + spec.radial_taper * z,
            );
            let t = (limit - ra) / (rb - ra);
            if t > 0. && t < 1. {
                let q = a + (b - a) * t;
                total += gauss2(radial, a, q) + gauss2(radial, q, b);
            } else {
                total += gauss2(radial, a, b);
            }
        }
        let q0 = source.endpoint(0, false).unwrap().0.midpoint();
        let q3 = source.endpoint(2, true).unwrap().0.midpoint();
        total += gauss2(radial, q3, q0 + spec.pitch_mm);
        std::f64::consts::PI / spec.pitch_mm * total
    };
    gauss2(section, 0., event) + gauss2(section, event, spec.length_mm)
}
#[test]
fn positive_source_volume_matches_independent_piecewise_original_law_quadrature() {
    for quarter in [false, true] {
        for scale in [0.1, 1., 3.] {
            let source = canonical(quarter, scale);
            let plan = source.stock_clip_plan(Default::default(), None).unwrap();
            let proof = plan.volume_bounds(1e-6, Default::default(), None).unwrap();
            let actual = independent_volume(source);
            let b = proof.source_volume_mm3();
            assert!(
                b[0] <= actual && actual <= b[1],
                "{quarter}/{scale}: {b:?} vs {actual}"
            );
            assert!(
                proof.work_used() <= 16_000_000
                    && proof.retained_budget_bytes() <= 64 * 1024 * 1024
            );
            let source_before = *plan.source();
            assert!(matches!(
                plan.volume_bounds(f64::from_bits(1), Default::default(), None),
                Err(Error::Unresolved)
            ));
            assert_eq!(*plan.source(), source_before);
            assert_eq!(
                plan.volume_bounds(1e-6, Default::default(), None)
                    .unwrap()
                    .source_volume_mm3(),
                b
            );
        }
    }
}
#[test]
fn source_domain_numerical_resource_and_cancellation_refusals_retain_originals() {
    let source = canonical(false, 1.);
    let plan = source.stock_clip_plan(Default::default(), None).unwrap();
    let patch = plan.patches()[0];
    for options in [
        StockClipOptions {
            maximum_work: 0,
            ..Default::default()
        },
        StockClipOptions {
            maximum_work: 16_000_001,
            ..Default::default()
        },
        StockClipOptions {
            maximum_bytes: 0,
            ..Default::default()
        },
        StockClipOptions {
            maximum_bytes: 64 * 1024 * 1024 + 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            source.stock_clip_plan(options, None),
            Err(Error::Invalid)
        ));
        assert!(matches!(
            plan.sample_bounds(patch, [0., 0.], [0., 0.], options, None),
            Err(Error::Invalid)
        ));
    }
    for options in [
        StockClipOptions {
            maximum_work: 1,
            ..Default::default()
        },
        StockClipOptions {
            maximum_bytes: 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            source.stock_clip_plan(options, None),
            Err(Error::Limit)
        ));
        assert!(matches!(
            plan.sample_bounds(patch, [0., 0.], [0., 0.], options, None),
            Err(Error::Limit)
        ));
        assert!(matches!(
            plan.volume_bounds(1e-6, options, None),
            Err(Error::Limit)
        ));
    }
    for p in [[-0.1, 0.], [0., 1.1], [1., 0.], [f64::NAN, 1.]] {
        assert!(matches!(
            plan.point_bounds(patch, p, [0., 0.], Default::default(), None),
            Err(Error::Invalid)
        ));
        assert!(matches!(
            plan.point_bounds(patch, [0., 0.], p, Default::default(), None),
            Err(Error::Invalid)
        ));
    }
    for width in [0., -1., f64::NAN, f64::INFINITY] {
        assert!(matches!(
            plan.volume_bounds(width, Default::default(), None),
            Err(Error::Invalid)
        ));
    }
    assert!(matches!(
        plan.cap_coedges(StockClipStation::Contact),
        Err(Error::Invalid)
    ));
    let stop = || true;
    assert!(matches!(
        source.stock_clip_plan(Default::default(), Some(&stop)),
        Err(Error::Cancelled)
    ));
    assert!(matches!(
        plan.sample_bounds(patch, [0., 0.], [0., 0.], Default::default(), Some(&stop)),
        Err(Error::Cancelled)
    ));
    assert!(matches!(
        plan.volume_bounds(1e-6, Default::default(), Some(&stop)),
        Err(Error::Cancelled)
    ));
    assert_eq!(*plan.source(), source);
    assert!(
        plan.sample_bounds(patch, [0., 1.], [0., 1.], Default::default(), None)
            .is_ok()
    );
}

#[test]
fn unsupported_full_clips_and_conditioned_events_are_not_redefined_or_snapped() {
    let original = canonical(false, 1.);
    let base = original.spec();
    let mut unsupported = vec![
        AuthoredSpec {
            crest_radius_limit_mm: None,
            ..base
        },
        AuthoredSpec {
            radial_taper: -1. / 32.,
            ..base
        },
        AuthoredSpec {
            crest_radius_limit_mm: Some(base.major_radius_mm - 0.1),
            ..base
        },
        AuthoredSpec {
            crest_radius_limit_mm: Some(
                base.major_radius_mm + base.radial_taper * base.length_mm + 1.,
            ),
            ..base
        },
    ];
    let mut rounded = super::super::tests::bspt_spec(1. / 32.);
    rounded.crest_radius_limit_mm = Some(rounded.major_radius_mm + 0.1);
    unsupported.push(rounded);
    for spec in unsupported {
        let source = AuthoredThread::new(spec, 1e-7, 1e-5).unwrap();
        assert!(matches!(
            source.stock_clip_plan(Default::default(), None),
            Err(Error::Unsupported)
        ));
        assert_eq!(source.spec(), spec);
    }
    let mut found = false;
    let mut length = base.length_mm;
    for _ in 0..64 {
        let spec = AuthoredSpec {
            length_mm: length,
            crest_radius_limit_mm: Some(base.major_radius_mm + base.radial_taper * length),
            ..base
        };
        let source = AuthoredThread::new(spec, 1e-7, 1e-5).unwrap();
        if matches!(
            source.stock_clip_plan(Default::default(), None),
            Err(Error::Unresolved)
        ) {
            found = true;
            assert_eq!(source.spec(), spec);
            break;
        }
        length = length.next_up();
    }
    assert!(
        found,
        "a genuine one-ulp contact/end conditioning frontier is exercised"
    );
    let spec = AuthoredSpec {
        section: AuthoredSection::NptFlat {
            depth_mm: 0.22,
            flat_width_mm: 0.07,
        },
        crest_radius_limit_mm: Some(base.major_radius_mm + 0.005),
        ..base
    };
    let gap = spec.crest_radius_limit_mm.unwrap()
        - (spec.major_radius_mm - 0.22 + spec.radial_taper * spec.length_mm);
    let mut found = false;
    for n in 1..16 {
        let minimum = gap - f64::from(n) * f64::EPSILON * spec.major_radius_mm;
        if let Ok(source) = AuthoredThread::new(spec, 1e-7, minimum)
            && matches!(
                source.stock_clip_plan(Default::default(), None),
                Err(Error::Unresolved)
            )
        {
            found = true;
            assert_eq!(source.spec(), spec);
            break;
        }
    }
    assert!(
        found,
        "actual stock/root clearance must be outward provable"
    );
    let thin = AuthoredSpec {
        crest_radius_limit_mm: Some(base.major_radius_mm.next_up()),
        ..base
    };
    let source = AuthoredThread::new(thin, 1e-7, 1e-5).unwrap();
    let plan = source.stock_clip_plan(Default::default(), None).unwrap();
    assert!(plan.contact_height_bounds_mm()[0] > 0. && plan.contact_height_bounds_mm()[1] < 1e-10);
    assert_eq!(plan.source().spec(), thin);
    assert!(
        plan.sample_bounds(
            plan.patches()[0],
            [0.5, 0.5],
            [0.5, 0.5],
            Default::default(),
            None
        )
        .is_ok()
    );
}
#[test]
fn every_actual_source_query_stage_honors_work_bytes_and_cancellation() {
    use std::cell::Cell;
    let source = canonical(false, 1.);
    let plan = source.stock_clip_plan(Default::default(), None).unwrap();
    let patch = plan.patches()[0];
    let point = plan
        .point_bounds(patch, [0.3, 0.3], [0.4, 0.4], Default::default(), None)
        .unwrap();
    let sample = plan
        .sample_bounds(patch, [0.3, 0.3], [0.4, 0.4], Default::default(), None)
        .unwrap();
    let volume = plan.volume_bounds(1e-6, Default::default(), None).unwrap();
    for work in [64, plan.work_used() - 1] {
        assert!(matches!(
            source.stock_clip_plan(
                StockClipOptions {
                    maximum_work: work,
                    ..Default::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
    }
    for work in [
        plan.work_used(),
        point.work_used() - 1,
        sample.point().work_used() - 1,
    ] {
        assert!(matches!(
            plan.sample_bounds(
                patch,
                [0.3, 0.3],
                [0.4, 0.4],
                StockClipOptions {
                    maximum_work: work,
                    ..Default::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
    }
    assert_eq!(
        plan.sample_bounds(
            patch,
            [0.3, 0.3],
            [0.4, 0.4],
            StockClipOptions {
                maximum_work: sample.point().work_used(),
                ..Default::default()
            },
            None
        )
        .unwrap(),
        sample
    );
    assert!(matches!(
        plan.volume_bounds(
            1e-6,
            StockClipOptions {
                maximum_work: volume.work_used() - 1,
                ..Default::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    for (bytes, which) in [
        (plan.retained_budget_bytes(), 0),
        (point.retained_budget_bytes(), 1),
        (volume.retained_budget_bytes(), 2),
    ] {
        let short = StockClipOptions {
            maximum_bytes: bytes - 1,
            ..Default::default()
        };
        let exact = StockClipOptions {
            maximum_bytes: bytes,
            ..Default::default()
        };
        match which {
            0 => {
                assert!(matches!(
                    source.stock_clip_plan(short, None),
                    Err(Error::Limit)
                ));
                assert!(source.stock_clip_plan(exact, None).is_ok());
            }
            1 => {
                assert!(matches!(
                    plan.point_bounds(patch, [0.3, 0.3], [0.4, 0.4], short, None),
                    Err(Error::Limit)
                ));
                assert_eq!(
                    plan.point_bounds(patch, [0.3, 0.3], [0.4, 0.4], exact, None)
                        .unwrap(),
                    point
                );
            }
            _ => {
                assert!(matches!(
                    plan.volume_bounds(1e-6, short, None),
                    Err(Error::Limit)
                ));
                assert_eq!(plan.volume_bounds(1e-6, exact, None).unwrap(), volume);
            }
        }
    }
    for at in 0..5 {
        let calls = Cell::new(0);
        let stop = || {
            let n = calls.get();
            calls.set(n + 1);
            n >= at
        };
        if at < 2 {
            assert!(matches!(
                source.stock_clip_plan(Default::default(), Some(&stop)),
                Err(Error::Cancelled)
            ));
        }
        calls.set(0);
        assert!(matches!(
            plan.sample_bounds(
                patch,
                [0.3, 0.3],
                [0.4, 0.4],
                Default::default(),
                Some(&stop)
            ),
            Err(Error::Cancelled)
        ));
        if at < 3 {
            calls.set(0);
            assert!(matches!(
                plan.volume_bounds(1e-6, Default::default(), Some(&stop)),
                Err(Error::Cancelled)
            ));
        }
    }
    assert_eq!(*plan.source(), source);
}
#[test]
fn mixed_source_jets_keep_original_chart_and_match_independent_differences() {
    let mut spec = super::super::tests::source(1. / 32.).spec();
    spec.crest_radius_limit_mm = Some(0.405 * 25.4 / 2.);
    let source = AuthoredThread::new(spec, 1e-7, 1e-5).unwrap();
    let plan = source
        .stock_clip_plan(StockClipOptions::default(), None)
        .unwrap();
    let opt = StockClipOptions::default();
    let mid = |p: StockClipPoint| p.position_mm().map(|a| a[0] + (a[1] - a[0]) * 0.5);
    for patch in plan.patches() {
        let law = plan.patch_law(patch, opt, None).unwrap();
        assert!(std::ptr::eq(law.plan(), &plan));
        assert_eq!(law.patch(), patch);
        let (u, v, h) = (0.37, 0.43, 1e-5);
        let jet = law.jets([u, u], [v, v], opt, None).unwrap();
        assert_eq!(
            jet.sample().point().position_mm(),
            plan.point_bounds(patch, [u, u], [v, v], opt, None)
                .unwrap()
                .position_mm()
        );
        let p = |a, b| mid(plan.point_bounds(patch, [a, a], [b, b], opt, None).unwrap());
        let center = p(u, v);
        let up = p(u + h, v);
        let um = p(u - h, v);
        let vp = p(u, v + h);
        let vm = p(u, v - h);
        let pp = p(u + h, v + h);
        let pm = p(u + h, v - h);
        let mp = p(u - h, v + h);
        let mm = p(u - h, v - h);
        let derivatives = [
            (
                jet.du(),
                std::array::from_fn::<_, 3, _>(|i| (up[i] - um[i]) / (2. * h)),
            ),
            (
                jet.dv(),
                std::array::from_fn(|i| (vp[i] - vm[i]) / (2. * h)),
            ),
            (
                jet.duu(),
                std::array::from_fn(|i| (up[i] - 2. * center[i] + um[i]) / (h * h)),
            ),
            (
                jet.duv(),
                std::array::from_fn(|i| (pp[i] - pm[i] - mp[i] + mm[i]) / (4. * h * h)),
            ),
            (
                jet.dvv(),
                std::array::from_fn(|i| (vp[i] - 2. * center[i] + vm[i]) / (h * h)),
            ),
        ];
        for (bound, fd) in derivatives {
            for i in 0..3 {
                let allowance = 1e-3 * (1. + fd[i].abs());
                assert!(
                    fd[i] >= bound[i][0] - allowance && fd[i] <= bound[i][1] + allowance,
                    "{patch:?} {bound:?} {fd:?}"
                );
            }
        }
        let broad = law.jets([u - h, u + h], [v - h, v + h], opt, None).unwrap();
        for (outer, inner) in [
            (broad.du(), jet.du()),
            (broad.dv(), jet.dv()),
            (broad.duu(), jet.duu()),
            (broad.duv(), jet.duv()),
            (broad.dvv(), jet.dvv()),
        ] {
            for i in 0..3 {
                assert!(outer[i][0] <= inner[i][0] && outer[i][1] >= inner[i][1]);
            }
        }
        for (u, v) in [(0., 0.), (1., 0.), (0., 1.), (1., 1.)] {
            let j = law.jets([u, u], [v, v], opt, None).unwrap();
            assert_eq!(
                j.sample().point().position_mm(),
                plan.point_bounds(patch, [u, u], [v, v], opt, None)
                    .unwrap()
                    .position_mm()
            );
        }
    }
}
#[test]
fn mixed_source_jet_queries_keep_work_bytes_cancel_and_invalid_domain_boundaries() {
    let mut spec = super::super::tests::source(1. / 32.).spec();
    spec.crest_radius_limit_mm = Some(0.405 * 25.4 / 2.);
    let source = AuthoredThread::new(spec, 1e-7, 1e-5).unwrap();
    let options = StockClipOptions::default();
    let plan = source.stock_clip_plan(options, None).unwrap();
    let patch = plan.patches()[0];
    let law = plan.patch_law(patch, options, None).unwrap();
    let jet = law.jets([0.2, 0.4], [0.3, 0.5], options, None).unwrap();
    for (work, bytes) in [
        (law.work_used(), law.retained_budget_bytes()),
        (jet.work_used(), jet.retained_budget_bytes()),
    ] {
        assert!(work > 128 && bytes > plan.retained_budget_bytes());
    }
    assert!(matches!(
        plan.patch_law(
            patch,
            StockClipOptions {
                maximum_work: law.work_used() - 1,
                ..options
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(
        plan.patch_law(
            patch,
            StockClipOptions {
                maximum_work: law.work_used(),
                maximum_bytes: law.retained_budget_bytes()
            },
            None
        )
        .is_ok()
    );
    assert!(matches!(
        plan.patch_law(
            patch,
            StockClipOptions {
                maximum_bytes: law.retained_budget_bytes() - 1,
                ..options
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        plan.patch_law(patch, options, Some(&|| true)),
        Err(Error::Cancelled)
    ));
    assert!(matches!(
        law.jets([-0.1, 0.], [0., 1.], options, None),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        law.jets([0., 1.], [0., 1.], options, Some(&|| true)),
        Err(Error::Cancelled)
    ));
    assert!(matches!(
        law.jets(
            [0.2, 0.4],
            [0.3, 0.5],
            StockClipOptions {
                maximum_work: jet.work_used() - 1,
                ..options
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(
        law.jets(
            [0.2, 0.4],
            [0.3, 0.5],
            StockClipOptions {
                maximum_work: jet.work_used(),
                maximum_bytes: jet.retained_budget_bytes()
            },
            None
        )
        .is_ok()
    );
    assert!(matches!(
        law.jets(
            [0.2, 0.4],
            [0.3, 0.5],
            StockClipOptions {
                maximum_bytes: jet.retained_budget_bytes() - 1,
                ..options
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert_eq!(*plan.source(), source);
}
