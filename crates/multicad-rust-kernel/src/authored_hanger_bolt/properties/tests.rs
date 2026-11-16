use super::*;
use crate::authored_hanger_bolt::{Error, Hanger};
use crate::authored_internal_treatment::{FaceKey as ThreadFace, ThreadBandRole};
use crate::tapered_thread::interval::I;
fn interval(e: Enclosure) -> I {
    I::new(e.lower, e.upper)
}
fn intersects(a: I, b: I) -> bool {
    a.lo <= b.hi && b.lo <= a.hi
}
fn source(left: bool, n: u32, nonbinary: bool) -> Hanger {
    let mut s = super::super::tests::fixture(left, nonbinary);
    let mut m = s.metric.spec();
    m.starts = n;
    s.metric = crate::authored_metric_thread::AuthoredMetricThread::new(m).unwrap();
    let mut w = s.wood.spec();
    w.starts = n;
    s.wood = crate::authored_flat_thread::AuthoredThread::new(
        w,
        crate::authored_flat_thread::SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-6,
        },
    )
    .unwrap();
    Hanger::new(s, None).unwrap()
}
#[test]
fn actual_default_hanger_faces_have_bounded_areas_moments_and_non_axis_caps() {
    let options = FacePropertyOptions::default();
    for (left, n, nb) in [
        (false, 1, false),
        (true, 1, true),
        (false, 2, true),
        (true, 2, false),
    ] {
        let s = source(left, n, nb);
        for i in 0..s.face_count() {
            let key = s.face_key(i).unwrap();
            let p = s
                .local_face_properties(key, options, None)
                .unwrap_or_else(|e| panic!("{left}/{n}/{nb} {key:?}: {e:?}"));
            assert_eq!(p.face, key);
            assert!(p.area_mm2.lower > 0.);
            assert!(interval(p.area_mm2).width() <= options.area_tolerance_mm2);
            for m in p.first_moment_mm3 {
                assert!(interval(m).width() <= options.first_moment_tolerance_mm3);
            }
            for c in p.centroid_mm {
                assert!(interval(c).width() <= options.centroid_tolerance_mm);
            }
            match key {
                FaceKey::PlainWall => assert!(intersects(
                    interval(p.centroid_mm[2]),
                    I::point(-s.spec().plain_length_mm * 0.5)
                )),
                FaceKey::MetricTip => assert!(intersects(
                    interval(p.centroid_mm[2]),
                    I::new(s.sections()[0].start_mm[0], s.sections()[0].start_mm[1])
                )),
                FaceKey::MetricShoulder => assert!(intersects(
                    interval(p.centroid_mm[2]),
                    I::point(-s.spec().plain_length_mm)
                )),
                FaceKey::WoodShoulder => {
                    assert!(intersects(interval(p.centroid_mm[2]), I::point(0.)))
                }
                FaceKey::WoodTip => assert!(intersects(
                    interval(p.centroid_mm[2]),
                    I::point(s.wood().full_length_mm())
                )),
                FaceKey::MetricThread { .. } => assert!(intersects(
                    interval(p.centroid_mm[2]),
                    -I::point(s.spec().plain_length_mm)
                        - I::point(s.metric().full_length_mm()) / I::point(2.)
                )),
                FaceKey::WoodThread { .. } => assert!(intersects(
                    interval(p.centroid_mm[2]),
                    I::point(s.wood().full_length_mm()) / I::point(2.)
                )),
            }
        }
        let v = s
            .volume(
                crate::authored_hanger_bolt::QueryOptions {
                    volume_error_mm3: 1e-7,
                    maximum_cells: 11,
                },
                None,
            )
            .unwrap();
        let metric_area =
            I::new(v.metric_mm3[0], v.metric_mm3[1]) / I::point(s.metric().full_length_mm());
        let wood_area = I::new(v.wood_mm3[0], v.wood_mm3[1]) / I::point(s.wood().full_length_mm());
        let tip = s
            .local_face_properties(FaceKey::MetricTip, options, None)
            .unwrap();
        assert!(intersects(interval(tip.area_mm2), metric_area));
        let tip = s
            .local_face_properties(FaceKey::WoodTip, options, None)
            .unwrap();
        assert!(intersects(interval(tip.area_mm2), wood_area));
        if n == 1 {
            assert!(
                tip.centroid_mm[..2]
                    .iter()
                    .any(|c| c.lower > 0. || c.upper < 0.)
            );
        }
    }
}
#[test]
fn crest_area_and_cap_volume_are_independent_source_oracles() {
    for left in [false, true] {
        let s = source(left, 1, true);
        let options = FacePropertyOptions::default();
        let spec = s.metric().spec();
        let p = s
            .local_face_properties(
                FaceKey::MetricThread {
                    start: 0,
                    band: ThreadBandRole::Crest,
                },
                options,
                None,
            )
            .unwrap();
        let tau = I::new(
            std::f64::consts::TAU.next_down(),
            std::f64::consts::TAU.next_up(),
        );
        let expected = I::point(spec.length_mm)
            * I::point(spec.major_radius_mm)
            * tau
            * I::point(spec.crest_width_mm)
            / I::point(spec.pitch_mm);
        assert!(intersects(interval(p.area_mm2), expected));
        let cap = s
            .metric()
            .local_face_properties(ThreadFace::StartCap, options, None)
            .unwrap();
        let own = s
            .local_face_properties(FaceKey::MetricTip, options, None)
            .unwrap();
        assert!(intersects(interval(own.area_mm2), interval(cap.area_mm2)));
        for i in 0..2 {
            assert!(intersects(
                interval(own.first_moment_mm3[i]),
                interval(cap.first_moment_mm3[i])
            ));
        }
        let wood = s
            .wood()
            .local_face_properties(ThreadFace::EndCap, options, None)
            .unwrap();
        let own = s
            .local_face_properties(FaceKey::WoodTip, options, None)
            .unwrap();
        assert_eq!(own.area_mm2, wood.area_mm2);
        assert_eq!(own.first_moment_mm3, wood.first_moment_mm3);
    }
}
#[test]
fn invalid_precision_resource_and_every_success_publication_checkpoint_are_atomic() {
    let s = source(false, 1, false);
    let before = s;
    let options = FacePropertyOptions::default();
    let root = FaceKey::MetricThread {
        start: 0,
        band: ThreadBandRole::Root,
    };
    for i in 0..s.face_count() {
        let face = s.face_key(i).unwrap();
        if face != FaceKey::PlainWall {
            assert_eq!(
                s.local_face_properties(
                    face,
                    FacePropertyOptions {
                        maximum_cells: 1,
                        ..options
                    },
                    None
                ),
                Err(Error::Limit)
            );
        }
    }

    assert_eq!(
        s.local_face_properties(
            FaceKey::MetricThread {
                start: 1,
                band: ThreadBandRole::Root
            },
            options,
            None
        ),
        Err(Error::Invalid)
    );
    assert_eq!(
        s.metric()
            .local_face_properties(ThreadFace::EntryCone, options, None),
        Err(crate::authored_internal_treatment::Error::Invalid)
    );
    for invalid in [
        FacePropertyOptions {
            area_tolerance_mm2: 0.,
            ..options
        },
        FacePropertyOptions {
            maximum_cells: 0,
            ..options
        },
        FacePropertyOptions {
            first_moment_tolerance_mm3: f64::NAN,
            ..options
        },
    ] {
        assert_eq!(
            s.local_face_properties(root, invalid, None),
            Err(Error::Invalid)
        );
    }
    assert_eq!(
        s.local_face_properties(
            root,
            FacePropertyOptions {
                maximum_cells: 1,
                ..options
            },
            None
        ),
        Err(Error::Limit)
    );
    assert_eq!(
        s.local_face_properties(
            FaceKey::PlainWall,
            FacePropertyOptions {
                area_tolerance_mm2: 1e-30,
                ..options
            },
            None
        ),
        Err(Error::Unresolved)
    );
    for face in [FaceKey::PlainWall, root] {
        let count = std::cell::Cell::new(0);
        s.local_face_properties(
            face,
            options,
            Some(&|| {
                count.set(count.get() + 1);
                false
            }),
        )
        .unwrap();
        let total = count.get();
        for at in [1, total] {
            count.set(0);
            assert_eq!(
                s.local_face_properties(
                    face,
                    options,
                    Some(&|| {
                        count.set(count.get() + 1);
                        count.get() == at
                    })
                ),
                Err(Error::Cancelled)
            );
        }
    }
    assert_eq!(s, before);
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn simpson_weight(i: usize, n: usize) -> f64 {
    if i == 0 || i == n {
        1.
    } else if i.is_multiple_of(2) {
        2.
    } else {
        4.
    }
}
#[test]
fn rounded_and_flat_flank_face_measures_match_independent_two_dimensional_tangent_cross_oracle() {
    use crate::authored_hanger_bolt::{BoundaryOptions, Component, CurveKey, End};
    for (left, n) in [(false, 1), (true, 2)] {
        let s = source(left, n, true);
        let boundary = s.boundary(BoundaryOptions::default(), None).unwrap();
        let spec = s.metric().spec();
        let rate = std::f64::consts::TAU / (spec.pitch_mm * f64::from(spec.starts))
            * if left { -1. } else { 1. };
        for band in [
            ThreadBandRole::Lower,
            ThreadBandRole::Crest,
            ThreadBandRole::Upper,
            ThreadBandRole::Root,
        ] {
            let face = FaceKey::MetricThread { start: 0, band };
            let curve = boundary
                .curves()
                .iter()
                .find(|c| {
                    matches!(
                        c.key,
                        CurveKey::Cap {
                            component: Component::Metric,
                            end: End::Start,
                            ..
                        }
                    ) && c.faces.contains(&face)
                })
                .unwrap();
            let (nu, nv) = (128, 1024);
            let mut sum = [0.; 4];
            for i in 0..=nu {
                let u = i as f64 / nu as f64;
                let base = boundary.curve_bounds(curve.key, [u, u], None).unwrap();
                let p = base.point_mm.map(|b| (b[0] + b[1]) * 0.5);
                let d = base.first_derivative_mm.map(|b| (b[0] + b[1]) * 0.5);
                for j in 0..=nv {
                    let v = j as f64 / nv as f64;
                    let phi = rate * spec.length_mm * v;
                    let (sn, cs) = phi.sin_cos();
                    let point = [
                        cs * p[0] - sn * p[1],
                        sn * p[0] + cs * p[1],
                        p[2] + spec.length_mm * v,
                    ];
                    let du = [cs * d[0] - sn * d[1], sn * d[0] + cs * d[1], 0.];
                    let dv = [
                        -rate * spec.length_mm * point[1],
                        rate * spec.length_mm * point[0],
                        spec.length_mm,
                    ];
                    let c = cross(du, dv);
                    let density = c.into_iter().map(|x| x * x).sum::<f64>().sqrt();
                    let weight = simpson_weight(i, nu) * simpson_weight(j, nv) * density;
                    sum[0] += weight;
                    for axis in 0..3 {
                        sum[axis + 1] += weight * point[axis];
                    }
                }
            }
            let normalizer = 9. * nu as f64 * nv as f64;
            sum = sum.map(|x| x / normalizer);
            let p = s
                .local_face_properties(face, FacePropertyOptions::default(), None)
                .unwrap();
            let actual = [
                interval(p.area_mm2).midpoint(),
                interval(p.first_moment_mm3[0]).midpoint(),
                interval(p.first_moment_mm3[1]).midpoint(),
                interval(p.first_moment_mm3[2]).midpoint(),
            ];
            for axis in 0..4 {
                assert!(
                    (actual[axis] - sum[axis]).abs() < 2e-4,
                    "left{left} n{n} {band:?} axis{axis}: certified{actual:?} oracle{sum:?}"
                );
            }
        }
    }
}
