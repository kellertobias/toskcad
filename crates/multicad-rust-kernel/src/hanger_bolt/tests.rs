use super::*;
use crate::{authored_flat_thread::SourceOptions, helical_thread::Spec as ThreadSpec};
use std::cell::Cell;
pub(super) fn fixture(left: bool, nb: bool) -> Spec {
    // Actual default_custom_connection_hanger_bolt dimensions and historical
    // M3 external6g midpoint construction from the canonical source rows.
    let wood = ThreadSpec {
        major_radius_mm: 3.,
        pitch_mm: 2.,
        length_mm: if nb { 12.3 } else { 12. },
        depth_mm: 1.,
        crest_width_mm: 0.4,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: left,
        root_radius_mm: None,
    };
    let major = (2.874 + 2.980) / 4.;
    let pitch_mid = (2.580 + 2.655) / 2.;
    let pitch = 0.5;
    let crest = pitch / 2. - (2. * major - pitch_mid) / 3_f64.sqrt();
    let radius = pitch / (4. * 3_f64.sqrt());
    let metric = ThreadSpec {
        major_radius_mm: major,
        pitch_mm: pitch,
        length_mm: if nb { 12.7 } else { 12. },
        depth_mm: (pitch - crest) / (2. / 3_f64.sqrt()) - radius,
        crest_width_mm: crest,
        root_radius_mm: Some(radius),
        ..wood
    };
    Spec {
        wood: AuthoredThread::new(
            wood,
            SourceOptions {
                position_limit_mm: 1e-7,
                normal_angle_limit_rad: 1e-7,
            },
        )
        .unwrap(),
        metric: Thread::new(metric).unwrap(),
        plain_length_mm: if nb { 2.3 } else { 2. },
        plain_radius_mm: 1.5,
    }
}
fn options() -> QueryOptions {
    QueryOptions {
        volume_error_mm3: 1e-7,
        maximum_cells: 65536,
    }
}
fn encloses(v: [f64; 2], value: f64) {
    assert!(v[0] <= value && value <= v[1], "{v:?} vs{value}");
}
fn independent_metric_volume(t: Thread) -> f64 {
    let b = t.profile_breaks();
    let root = t.root_circle().unwrap();
    let mut sum = 0.;
    for band in 0..5 {
        let [a, c] = [b[band], b[band + 1]];
        sum += if band == 2 {
            let primitive = |q: f64| {
                let x = q - root.phase;
                let h = (root.radius * root.radius - x * x).sqrt();
                (root.center_radius * root.center_radius + root.radius * root.radius) * x
                    - x.powi(3) / 3.
                    - root.center_radius
                        * (x * h + root.radius * root.radius * (x / root.radius).asin())
            };
            primitive(c) - primitive(a)
        } else {
            let r = |q: f64| {
                if band == 1 || band == 3 {
                    let base = if band == 1 { b[1] } else { b[4] };
                    t.spec().major_radius_mm + t.slope_at_phase(base, band) * (q - base)
                } else {
                    t.spec().major_radius_mm
                }
            };
            let [x, y] = [r(a), r(c)];
            (c - a) * (x * x + x * y + y * y) / 3.
        };
    }
    std::f64::consts::PI * t.spec().length_mm / t.spec().pitch_mm * sum
}
#[test]
fn actual_canonical_hanger_retains_both_sources_and_connected_sections_with_directed_volume() {
    for left in [false, true] {
        for nb in [false, true] {
            let spec = fixture(left, nb);
            let before = spec;
            let h = Hanger::new(spec, None).unwrap();
            assert_eq!(h.spec(), before);
            assert_eq!(*h.wood(), before.wood);
            assert_eq!(*h.metric(), before.metric);
            assert!(h.connection_core_radius_mm() > 1. && h.connection_core_radius_mm() < 1.5);
            let sections = h.sections();
            assert_eq!(
                sections.map(|s| s.kind),
                [SectionKind::Metric, SectionKind::Plain, SectionKind::Wood]
            );
            let exact = -BigRational::from_float(spec.plain_length_mm).unwrap()
                - BigRational::from_float(spec.metric.spec().length_mm).unwrap();
            assert!(
                BigRational::from_float(sections[0].start_mm[0]).unwrap() <= exact
                    && exact <= BigRational::from_float(sections[0].start_mm[1]).unwrap()
            );
            assert_eq!(sections[0].end_mm, sections[1].start_mm);
            assert_eq!(sections[1].end_mm, sections[2].start_mm);
            let volume = h.volume(options(), None).unwrap();
            encloses(volume.metric_mm3, independent_metric_volume(spec.metric));
            let wood = spec
                .wood
                .bands()
                .into_iter()
                .map(|b| {
                    let [a, c] = b.radius_mm;
                    (b.phase_mm[1] - b.phase_mm[0]) * (a * a + a * c + c * c) / 3.
                })
                .sum::<f64>()
                * std::f64::consts::PI
                * spec.wood.full_length_mm()
                / spec.wood.spec().pitch_mm;
            let plain = std::f64::consts::PI * spec.plain_radius_mm.powi(2) * spec.plain_length_mm;
            encloses(volume.wood_mm3, wood);
            encloses(volume.plain_mm3, plain);
            encloses(
                volume.material_mm3,
                wood + plain + independent_metric_volume(spec.metric),
            );
            assert!(volume.material_mm3[1] - volume.material_mm3[0] <= options().volume_error_mm3);
            assert!(volume.metric_phase_error_mm3 > 0.);
            assert_eq!(h.spec(), before);
            assert_eq!(h.bounds_mm()[2][0], sections[0].start_mm[0]);
        }
    }
}
#[test]
fn whole_boxes_respect_real_threads_plain_shoulder_joins_tips_and_original_metric_phase() {
    use Classification::*;
    for left in [false, true] {
        let h = Hanger::new(fixture(left, true), None).unwrap();
        for z in [-10., -2.3, 0., 3.] {
            assert_eq!(
                h.classify_box([[0., 0.], [0., 0.], [z - 0.01, z + 0.01]], None),
                Ok(Inside)
            );
        }
        assert_eq!(
            h.classify_box([[4., 4.], [0., 0.], [1., 1.]], None),
            Ok(Outside)
        );
        assert_eq!(
            h.classify_box([[0., 0.], [0., 0.], [-20., -19.]], None),
            Ok(Outside)
        );
        assert_eq!(
            h.classify_box([[0., 0.], [0., 0.], [20., 21.]], None),
            Ok(Outside)
        );
        assert_eq!(
            h.classify_box([[1.6, 1.6], [0., 0.], [-1., -1.]], None),
            Ok(Outside)
        );
        assert_eq!(
            h.classify_box([[1.5, 1.5], [0., 0.], [-1., -1.]], None),
            Ok(UnresolvedBoundary)
        );
        assert_eq!(
            h.classify_box([[2.5, 2.5], [0., 0.], [0., 0.]], None),
            Ok(UnresolvedBoundary)
        );
        let bottom = h.sections()[0].start_mm;
        assert_eq!(
            h.classify_box([[0., 0.], [0., 0.], bottom], None),
            Ok(UnresolvedBoundary)
        );
        for theta in [0.2, 1.1, 3., 5.7] {
            let local = 5.23;
            let z = local - h.spec.plain_length_mm - h.metric().spec().length_mm;
            let radius = h.metric().radius(theta, local);
            for (offset, expected) in [(-0.005, Inside), (0.005, Outside)] {
                let r = radius + offset;
                let (x, y) = (r * theta.cos(), r * theta.sin());
                assert_eq!(h.classify_box([[x, x], [y, y], [z, z]], None), Ok(expected));
            }
        }
        // The wood root/flank is retained rather than replaced by its envelope.
        let z = 1.;
        assert_eq!(
            h.classify_box([[2.9, 2.9], [0., 0.], [z, z]], None),
            Ok(Outside)
        );
        assert_eq!(
            h.classify_box([[-2.1, 2.1], [-2.1, 2.1], [0.1, 1.9]], None),
            Ok(UnresolvedBoundary)
        );
    }
}
#[test]
fn invalid_unsupported_precision_budget_and_every_query_checkpoint_refuse_atomically() {
    let spec = fixture(false, false);
    let h = Hanger::new(spec, None).unwrap();
    for x in [0., -1., f64::NAN, f64::INFINITY, 1e8] {
        assert_eq!(
            Hanger::new(
                Spec {
                    plain_length_mm: x,
                    ..spec
                },
                None
            ),
            Err(Error::Invalid)
        );
    }
    for r in [1., 2., spec.wood.minimum_radius_mm()] {
        assert_eq!(
            Hanger::new(
                Spec {
                    plain_radius_mm: r,
                    ..spec
                },
                None
            ),
            Err(Error::Unsupported)
        );
    }
    let flat = Thread::new(ThreadSpec {
        root_radius_mm: None,
        depth_mm: 0.2,
        ..spec.metric.spec()
    })
    .unwrap();
    assert_eq!(
        Hanger::new(
            Spec {
                metric: flat,
                ..spec
            },
            None
        ),
        Err(Error::Unsupported)
    );
    assert_eq!(
        h.volume(
            QueryOptions {
                volume_error_mm3: 0.,
                ..options()
            },
            None
        ),
        Err(Error::Invalid)
    );
    for maximum_cells in [0, 1, 65537] {
        assert_eq!(
            h.volume(
                QueryOptions {
                    maximum_cells,
                    ..options()
                },
                None
            ),
            Err(Error::Limit)
        );
    }
    assert_eq!(
        h.volume(
            QueryOptions {
                volume_error_mm3: 1e-30,
                ..options()
            },
            None
        ),
        Err(Error::Unresolved)
    );
    assert_eq!(
        h.classify_box([[f64::NAN, 0.]; 3], None),
        Err(Error::Invalid)
    );
    assert_eq!(
        h.classify_box([[0., 0.], [0., 0.], [f64::MIN, f64::MAX]], None),
        Ok(Classification::UnresolvedBoundary)
    );
    let count = Cell::new(0);
    let poll = || {
        count.set(count.get() + 1);
        false
    };
    h.volume(options(), Some(&poll)).unwrap();
    for at in 1..=count.get() {
        let c = Cell::new(0);
        let cancel = || {
            c.set(c.get() + 1);
            c.get() == at
        };
        assert_eq!(h.volume(options(), Some(&cancel)), Err(Error::Cancelled));
    }
    assert_eq!(Hanger::new(spec, Some(&|| true)), Err(Error::Cancelled));
    assert_eq!(
        h.classify_box([[0., 0.], [0., 0.], [1., 1.]], Some(&|| true)),
        Err(Error::Cancelled)
    );
    assert_eq!(h.spec(), spec);
}
#[test]
fn directed_box_transfers_keep_quadrants_and_plain_and_metric_boundary_ambiguity() {
    use Classification::*;
    for left in [false, true] {
        let h = Hanger::new(fixture(left, false), None).unwrap();
        assert_eq!(
            h.classify_box([[1.3, 1.3], [0., 0.], [-1., -1.]], None),
            Ok(Inside)
        );
        assert_eq!(
            h.classify_box([[f64::MIN, f64::MAX], [f64::MIN, f64::MAX], [1., 1.]], None),
            Ok(UnresolvedBoundary)
        );
        let local = 5.23;
        let z = local - h.spec.plain_length_mm - h.metric().spec().length_mm;
        let r = h.metric().radius(0.2, local);
        assert_eq!(
            h.classify_box(
                [
                    [(r - 0.01) * 0.2_f64.cos(), (r + 0.01) * 0.2_f64.cos()],
                    [(r - 0.01) * 0.2_f64.sin(), (r + 0.01) * 0.2_f64.sin()],
                    [z, z]
                ],
                None
            ),
            Ok(UnresolvedBoundary)
        );
        for [x, y] in [[-1_f64, 0.2], [-1., -0.2], [0., 1.], [0., -1.]] {
            let theta = y.atan2(x);
            let norm = x.hypot(y);
            let z = 1.23;
            let rate =
                (if left { -1. } else { 1. }) * h.wood().spec().pitch_mm / std::f64::consts::TAU;
            let q = (z - rate * theta).rem_euclid(h.wood().spec().pitch_mm);
            let b = h
                .wood()
                .bands()
                .into_iter()
                .find(|b| q <= b.phase_mm[1])
                .unwrap();
            let r = b.radius_mm[0]
                + (b.radius_mm[1] - b.radius_mm[0]) * (q - b.phase_mm[0])
                    / (b.phase_mm[1] - b.phase_mm[0]);
            for (delta, expected) in [(-0.005, Inside), (0.005, Outside)] {
                let xx = (r + delta) * x / norm;
                let yy = (r + delta) * y / norm;
                assert_eq!(
                    h.classify_box([[xx, xx], [yy, yy], [z, z]], None),
                    Ok(expected)
                );
            }
        }
    }
    let mut s = fixture(false, false);
    let mut ws = s.wood.spec();
    ws.major_radius_mm = 1e5;
    s.wood = AuthoredThread::new(
        ws,
        SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-7,
        },
    )
    .unwrap();
    let wide = Hanger::new(s, None).unwrap();
    assert_eq!(wide.volume(options(), None), Err(Error::Unresolved));
}

#[test]
fn tiny_but_real_rounded_root_never_bypasses_the_derivative_precision_floor() {
    let base = fixture(false, false);
    let mut spec = base.metric.spec();
    let radius = (0.25_f64.next_up() - 0.25).next_up();
    spec.root_radius_mm = Some(radius);
    spec.depth_mm = (spec.pitch_mm - spec.crest_width_mm) / (2. / 3_f64.sqrt()) - radius;
    spec.major_radius_mm = spec.depth_mm + 0.2;
    let thread = Thread::new(spec).unwrap();
    let hanger = Hanger::new(
        Spec {
            metric: thread,
            ..base
        },
        None,
    )
    .unwrap();
    assert_eq!(hanger.metric().spec().root_radius_mm, Some(radius));
    assert!(matches!(
        hanger.volume(options(), None),
        Err(Error::Unresolved | Error::Limit)
    ));
    assert_eq!(*hanger.metric(), thread);
}

#[test]
fn exact_periodic_wood_domain_does_not_reuse_binary_connected_end_caps() {
    let original = fixture(false, false);
    assert!(Hanger::new(original, None).is_ok());
    let wood = AuthoredThread::new_with_exact_axial_domain(
        original.wood.spec(),
        crate::Rational::new(123, 10).unwrap(),
        SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-6,
        },
        None,
    )
    .unwrap();
    assert_eq!(
        Hanger::new(Spec { wood, ..original }, None),
        Err(Error::Unsupported)
    );
    assert!(Hanger::new(original, None).is_ok());
}
