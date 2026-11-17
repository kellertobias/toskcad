use super::*;
use num_rational::BigRational as Q;
use std::sync::Arc;
fn q(x: f64) -> Q {
    Q::from_float(x).unwrap()
}
fn work() -> Work<'static> {
    Work {
        used: 0,
        remaining: 10_000_000,
        base: 0,
        held: 0,
        peak: 0,
        maximum_bytes: MAX_BYTES,
        interrupt: None,
    }
}
fn cross(a: [Q; 3], b: [Q; 3]) -> [Q; 3] {
    std::array::from_fn(|i| &a[(i + 1) % 3] * &b[(i + 2) % 3] - &a[(i + 2) % 3] * &b[(i + 1) % 3])
}
fn dot(a: &[Q; 3], b: &[Q; 3]) -> Q {
    (0..3).map(|i| &a[i] * &b[i]).sum()
}
#[test]
fn transported_actual_tangent_cone_has_independent_exact_normal_oracles_and_is_not_a_chord() {
    for scale in [0.25, 0.5, 1., 2., 4.] {
        for (f, n) in [
            ([1., 0., 0.], [1., 0.05, 0.025]),
            ([1., 2., 3.], [1.1, 2.05, 3.02]),
        ] {
            let f = f.map(q);
            let n = n.map(q);
            let r = 0.1;
            let c = cross(f.clone(), n.clone());
            let d = dot(&f, &n);
            assert!(d > q(0.) && dot(&c, &c) <= q(r).pow(2) * d.pow(2));
            let k = [q(scale), q(scale), q(1.)];
            let mapped_f = std::array::from_fn(|i| &k[i] * &f[i]);
            let mapped_n = std::array::from_fn(|i| &k[i] * &n[i]);
            let c = cross(mapped_f.clone(), mapped_n.clone());
            let d = dot(&mapped_f, &mapped_n);
            let bound = normal::transported(I::point(scale), r, &mut work()).unwrap();
            assert!(d > q(0.) && dot(&c, &c) <= q(bound).pow(2) * d.pow(2));
        }
    }
    assert!(normal::transported(I::point(1.), 0.5, &mut work()).unwrap() >= 0.5);
    assert!(matches!(
        normal::transported(I::point(0.01), 0.1, &mut work()),
        Err(Error::Unsupported)
    ));
    let a = 0.125;
    let b = 0.25;
    let bound = normal::combine(a, b, &mut work()).unwrap();
    assert!(q(bound) >= (q(a) + q(b)) / (q(1.) - q(a) * q(b)));
    assert!(matches!(
        normal::combine(2., 0.5, &mut work()),
        Err(Error::Unsupported)
    ));
    let triangle = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]];
    let mapped = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]];
    assert!(normal::rounding(triangle, mapped, I::point(2.), &mut work()).unwrap() < 1e-10);
    assert!(matches!(
        normal::rounding(
            triangle,
            [mapped[0], mapped[2], mapped[1]],
            I::point(2.),
            &mut work()
        ),
        Err(Error::Unsupported)
    ));
}
#[test]
fn actual_clipped_head_mesh_probe_both_sizes_and_strict_quality() {
    let mut failed = false;
    for quarter in [false, true] {
        let input: Arc<str> = Arc::from(super::super::source_tests::candidate(quarter, 18.));
        let owner = OwnedActualHead::read(input, Default::default(), None).unwrap();
        for deflection in [0.1, 0.05] {
            let result = owner.mesh_candidate(
                HeadMeshOptions {
                    tessellation: Tessellation {
                        deflection,
                        angle: 0.3,
                        triangle_limit: 500_000,
                    },
                    ..Default::default()
                },
                None,
            );
            match result {
                Ok(mesh) => {
                    let r = mesh.receipt();
                    eprintln!(
                        "ACTUAL_HEAD_MESH quarter={quarter} q={deflection} vertices={} triangles={} {r:?}",
                        mesh.mesh().vertices.len() / 3,
                        mesh.mesh().face_indices.len()
                    );
                    assert!(
                        r.source_position_error_mm <= deflection && r.normal_tangent_cone <= 0.3
                    );
                    assert!(mesh.mesh().face_indices.iter().all(|id| *id < 17));
                    assert_eq!(
                        r.aggregate_work,
                        r.inherited_head_work + r.auxiliary_mesh_work - r.shared_auxiliary_work
                            + r.joining_work
                    );
                    let volume = owner.local_sewn_volume_mm3();
                    let actual = mesh.connectivity().signed_volume_mm3;
                    assert!(actual > volume[0] * 0.95 && actual < volume[1] * 1.05);
                }
                Err(e) => {
                    eprintln!("ACTUAL_HEAD_MESH quarter={quarter} q={deflection} ERROR {e:?}");
                    failed = true;
                }
            }
        }
    }
    assert!(
        !failed,
        "actual mesh controls must pass without caller-quality/quota relaxation"
    );
}
fn decimal(s: &str) -> Q {
    use num_bigint::BigInt;
    let (m, e) = s.split_once('E').unwrap_or((s, "0"));
    let fraction = m.split_once('.').map_or(0, |(_, v)| v.len()) as i32;
    let n = m.replace('.', "").parse::<BigInt>().unwrap();
    let power = e.parse::<i32>().unwrap() - fraction;
    if power >= 0 {
        Q::from_integer(n * BigInt::from(10).pow(power as u32))
    } else {
        Q::new(n, BigInt::from(10).pow((-power) as u32))
    }
}
fn evaluate(mut p: [[Q; 3]; 12], t: Q) -> [Q; 3] {
    for n in (1..12).rev() {
        for i in 0..n {
            p[i] = std::array::from_fn(|k| (q(1.) - &t) * &p[i][k] + &t * &p[i + 1][k]);
        }
    }
    p[0].clone()
}
#[test]
fn actual_cap_trace_chord_bound_contains_independent_original_decimal_curve_and_global_transfer() {
    use num_traits::Signed;
    let owner = OwnedActualHead::read(
        Arc::from(super::super::source_tests::candidate(true, 18.)),
        Default::default(),
        None,
    )
    .unwrap();
    let source = owner.auxiliary_thread();
    let candidate = source.mesh_candidate(Default::default(), None).unwrap();
    let mut heights = source
        .faces()
        .iter()
        .filter_map(|f| match f.support() {
            Support::Plane { height } => Some(decimal(
                source.point(height).unwrap().coordinates()[2].spelling(),
            )),
            _ => None,
        })
        .collect::<Vec<_>>();
    heights.sort();
    let lower = &heights[0];
    let upper = &heights[1];
    let scale = (q(owner.head_stations_mm()[0]) - lower) / (upper - lower);
    let mut checked = 0;
    for cap in candidate.cap_traces() {
        assert_eq!(cap.rim.len(), cap.edges.len());
        for i in [0, cap.edges.len() / 2, cap.edges.len() - 1] {
            let e = cap.edges[i];
            assert_eq!(e.vertices, [cap.rim[i], cap.rim[(i + 1) % cap.rim.len()]]);
            let actual = e.vertices.map(|id| {
                owner
                    .axial_sewing_transfer()
                    .map_raw_box_mm(vertex(candidate.mesh(), id).unwrap().map(|x| [x, x]))
                    .unwrap()
                    .map(|b| I::new(b[0], b[1]).midpoint())
            });
            let bound = caps::chord(
                source,
                owner.axial_sewing_transfer(),
                e.edge,
                e.span,
                e.parameters,
                actual,
                &mut work(),
            )
            .unwrap();
            let curve = &source.curves()[source.edges()[e.edge].curve_index()];
            let poles = std::array::from_fn::<_, 12, _>(|j| {
                source
                    .point(curve.point_indices()[11 * e.span + j])
                    .unwrap()
                    .coordinates()
                    .map(|s| decimal(s.spelling()))
            });
            for sample in 0..=32 {
                let u = Q::new(sample.into(), 32.into());
                let t = q(e.parameters[0]) + (q(e.parameters[1]) - q(e.parameters[0])) * &u;
                let mut point = evaluate(poles.clone(), t);
                point[2] = lower + (&point[2] - lower) * &scale;
                let chord = std::array::from_fn::<_, 3, _>(|k| {
                    (q(1.) - &u) * q(actual[0][k]) + &u * q(actual[1][k])
                });
                let delta: Q = (0..3).map(|k| (&point[k] - &chord[k]).abs()).sum();
                assert!(delta <= q(bound));
            }
            checked += 1;
        }
    }
    assert_eq!(checked, 6);
}
#[test]
fn actual_head_mesh_real_work_storage_and_cancel_frontiers_preserve_source_and_retry() {
    use std::cell::Cell;
    let owner = OwnedActualHead::read(
        Arc::from(super::super::source_tests::candidate(false, 18.)),
        Default::default(),
        None,
    )
    .unwrap();
    let good = owner.mesh_candidate(Default::default(), None).unwrap();
    let r = good.receipt();
    let mesh = good.mesh().clone();
    drop(good);
    for (maximum_work, maximum_bytes, ok) in [
        (r.aggregate_work - 1, MAX_BYTES, false),
        (r.aggregate_work, MAX_BYTES, true),
        (MAX_MESH_WORK, r.accounted_peak_bytes - 1, false),
        (MAX_MESH_WORK, r.accounted_peak_bytes, true),
    ] {
        let result = owner.mesh_candidate(
            HeadMeshOptions {
                maximum_work,
                maximum_bytes,
                ..Default::default()
            },
            None,
        );
        assert_eq!(result.is_ok(), ok);
        if !ok {
            assert!(matches!(result, Err(Error::Limit)));
        }
    }
    let n = Cell::new(0);
    let observe = || {
        n.set(n.get() + 1);
        false
    };
    owner
        .mesh_candidate(Default::default(), Some(&observe))
        .unwrap();
    for at in [1, n.get() / 2, n.get()] {
        let seen = Cell::new(0);
        let cancel = || {
            seen.set(seen.get() + 1);
            seen.get() >= at
        };
        assert!(matches!(
            owner.mesh_candidate(Default::default(), Some(&cancel)),
            Err(Error::Cancelled)
        ));
    }
    let extra = owner
        .mesh_candidate(
            HeadMeshOptions {
                additional_resident_bytes: 4096,
                ..Default::default()
            },
            None,
        )
        .unwrap();
    assert_eq!(extra.receipt().aggregate_work, r.aggregate_work);
    assert_eq!(
        extra.receipt().accounted_peak_bytes,
        r.accounted_peak_bytes + 4096
    );
    let retry = owner.mesh_candidate(Default::default(), None).unwrap();
    assert_eq!(retry.mesh(), &mesh);
    assert_eq!(retry.receipt().aggregate_work, r.aggregate_work);
    assert!(matches!(
        owner.mesh_candidate(
            HeadMeshOptions {
                tessellation: Tessellation {
                    angle: 0.,
                    ..HeadMeshOptions::default().tessellation
                },
                ..Default::default()
            },
            None
        ),
        Err(Error::Invalid)
    ));
}
