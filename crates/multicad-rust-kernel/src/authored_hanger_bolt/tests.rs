use super::*;
use crate::{
    authored_flat_thread::SourceOptions, authored_metric_thread::Spec as MetricSpec,
    helical_thread::Spec as WoodSpec,
};
pub(crate) fn fixture(left: bool, nb: bool) -> Spec {
    let p = 0.5;
    let major = (2.874 + 2.980) / 4.;
    let pd = (2.580 + 2.655) / 2.;
    let metric = AuthoredMetricThread::new(MetricSpec {
        major_radius_mm: major,
        pitch_mm: p,
        length_mm: if nb { 12.7 } else { 12. },
        crest_width_mm: p / 2. - (2. * major - pd) / 3_f64.sqrt(),
        starts: 1,
        left_hand: left,
    })
    .unwrap();
    let wood = AuthoredThread::new(
        WoodSpec {
            major_radius_mm: 3.,
            pitch_mm: 2.,
            length_mm: if nb { 12.3 } else { 12. },
            depth_mm: 1.,
            crest_width_mm: 0.4,
            upper_flank_deg: 30.,
            lower_flank_deg: 30.,
            starts: 1,
            left_hand: left,
            root_radius_mm: None,
        },
        SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-6,
        },
    )
    .unwrap();
    Spec {
        wood,
        metric,
        plain_length_mm: if nb { 2.15 } else { 2. },
        plain_radius_mm: 1.5,
    }
}
fn options() -> QueryOptions {
    QueryOptions {
        volume_error_mm3: 1e-7,
        maximum_cells: 11,
    }
}
#[test]
fn canonical_material_retains_actual_two_laws_sections_and_directed_volume() {
    for (left, nb) in [(false, false), (true, true)] {
        let spec = fixture(left, nb);
        let source = Hanger::new(spec, None).unwrap();
        assert_eq!(source.spec(), spec);
        assert_eq!(*source.wood(), spec.wood);
        assert_eq!(*source.metric(), spec.metric);
        let sections = source.sections();
        assert_eq!(
            sections.map(|s| s.kind),
            [SectionKind::Metric, SectionKind::Plain, SectionKind::Wood]
        );
        let q = |x| num_rational::BigRational::from_float(x).unwrap();
        let bottom = -q(spec.plain_length_mm) - q(spec.metric.full_length_mm());
        assert!(q(sections[0].start_mm[0]) <= bottom && bottom <= q(sections[0].start_mm[1]));
        assert_eq!(sections[0].end_mm, sections[1].start_mm);
        assert_eq!(sections[1].end_mm, sections[2].start_mm);
        assert!(source.connection_core_radius_mm() > 1.);
        assert_eq!(source.bounds_mm()[0], [-3., 3.]);
        let v = source.volume(options(), None).unwrap();
        assert!(v.material_mm3[0] > 0. && v.material_mm3[1] - v.material_mm3[0] <= 1e-7);
        assert!(
            v.plain_mm3[0] <= std::f64::consts::PI * 1.5 * 1.5 * spec.plain_length_mm
                && std::f64::consts::PI * 1.5 * 1.5 * spec.plain_length_mm <= v.plain_mm3[1]
        );
        // Independent dense angular/phase midpoint oracle of the actual circular
        // profile. It is a numerical cross-check, not the directed authority.
        let n = 65536;
        let mut sum = 0.;
        let p = spec.metric.spec().pitch_mm;
        let knots = [
            crate::authored_metric_thread::Knot::Begin,
            crate::authored_metric_thread::Knot::CrestEnd,
            crate::authored_metric_thread::Knot::UpperRoot,
            crate::authored_metric_thread::Knot::LowerRoot,
            crate::authored_metric_thread::Knot::CrestStart,
            crate::authored_metric_thread::Knot::End,
        ]
        .map(|k| spec.metric.knot(k));
        for i in 0..n {
            let phase = (i as f64 + 0.5) / n as f64 * p;
            let band = (0..5).find(|&b| phase <= knots[b + 1].phase_mm[1]).unwrap();
            let a = knots[band].phase_mm[0].midpoint(knots[band].phase_mm[1]);
            let b = knots[band + 1].phase_mm[0].midpoint(knots[band + 1].phase_mm[1]);
            let u = (phase - a) / (b - a);
            let r = spec.metric.profile_bounds(band, [u, u]).unwrap().radius_mm;
            sum += r[0].midpoint(r[1]).powi(2);
        }
        let oracle = std::f64::consts::PI * spec.metric.full_length_mm() * sum / n as f64;
        assert!((v.metric_mm3[0].midpoint(v.metric_mm3[1]) - oracle).abs() < 1e-7);
    }
}
fn point(x: f64, y: f64, z: f64) -> [[f64; 2]; 3] {
    [[x, x], [y, y], [z, z]]
}
#[test]
fn actual_phase_boxes_distinguish_wood_round_root_plain_and_real_join_domains() {
    use Classification::*;
    for left in [false, true] {
        let spec = fixture(left, false);
        let source = Hanger::new(spec, None).unwrap();
        for z in [-13., -2., 0., 3.] {
            assert_eq!(source.classify_box(point(0.5, 0., z), None), Ok(Inside));
        }
        assert_eq!(
            source.classify_box([[1.3, 1.4], [0., 0.], [-0.1, 0.1]], None),
            Ok(Inside)
        );
        assert_eq!(source.classify_box(point(1.4, 0., -1.), None), Ok(Inside));
        assert_eq!(source.classify_box(point(1.6, 0., -1.), None), Ok(Outside));
        assert_eq!(
            source.classify_box(point(1.5, 0., -1.), None),
            Ok(UnresolvedBoundary)
        );
        for z in [-15., 13.] {
            assert_eq!(source.classify_box(point(0., 0., z), None), Ok(Outside));
        }
        assert_eq!(
            source.classify_box(point(0., 0., -14.), None),
            Ok(UnresolvedBoundary)
        );
        assert_eq!(
            source.classify_box(point(0., 0., 12.), None),
            Ok(UnresolvedBoundary)
        );
        assert_eq!(
            source.classify_box([[1.3, 1.4], [0., 0.], [-2.1, -1.9]], None),
            Ok(UnresolvedBoundary)
        );
        // Original full finite metric phase, translated symbolically by P+Lm.
        let shift = spec.plain_length_mm + spec.metric.full_length_mm();
        assert_eq!(
            source.classify_box(point(1.3, 0., 1. - shift), None),
            Ok(Inside)
        );
        assert_eq!(
            source.classify_box(point(1.3, 0., 1.25 - shift), None),
            Ok(Outside)
        );
        assert_eq!(
            source.classify_box(point(0., 1.3, 1.125 - shift), None),
            Ok(if left { Outside } else { Inside })
        );
        assert_eq!(
            source.classify_box(point(0., -1.3, 1.125 - shift), None),
            Ok(if left { Inside } else { Outside })
        );
        assert_eq!(
            source.classify_box([[1.2, 1.4], [0., 0.], [1.1 - shift, 1.2 - shift]], None),
            Ok(UnresolvedBoundary)
        );
        // General phase, including the later assembly's non-quarter nut roll.
        for theta in [0.37, (4.1 + 0.2) / 0.5 * std::f64::consts::TAU] {
            let rate =
                (if left { -1. } else { 1. }) * spec.metric.spec().pitch_mm / std::f64::consts::TAU;
            let (sin, cos) = theta.sin_cos();
            for (phase, expected) in [(6., Inside), (6.25, Outside)] {
                assert_eq!(
                    source.classify_box(
                        point(1.3 * cos, 1.3 * sin, phase + rate * theta - shift),
                        None
                    ),
                    Ok(expected)
                );
            }
        }
        // True shared flat wood source: no cylinder or bbox membership proxy.
        for band in 0..5 {
            for z in [0.7, 1.3] {
                let sample = spec.wood.point_in_band(band, 0, 0.5, z).unwrap();
                let [x, y, z] = sample.point_mm;
                assert_eq!(
                    source.classify_box(point(x * 0.98, y * 0.98, z), None),
                    Ok(Inside)
                );
                assert_eq!(
                    source.classify_box(point(x * 1.02, y * 1.02, z), None),
                    Ok(Outside)
                );
            }
        }
        assert_eq!(
            source.classify_box(point(0., 2.5, 0.5), None),
            Ok(if left { Outside } else { Inside })
        );
        assert_eq!(
            source.classify_box(point(0., -2.5, 0.5), None),
            Ok(if left { Inside } else { Outside })
        );
        assert_eq!(source.classify_box(point(-2.5, 0., 1.), None), Ok(Inside));
        assert_eq!(
            source.classify_box(point(-2.5, -0.001, 1.), None),
            Ok(Inside)
        );
        assert_eq!(
            source.classify_box([[0., 0.], [-2.5, 2.5], [1., 2.]], None),
            Ok(UnresolvedBoundary)
        );
        assert_eq!(
            source.classify_box(
                [[-f64::MAX, f64::MAX], [-f64::MAX, f64::MAX], [1., 2.]],
                None
            ),
            Ok(UnresolvedBoundary)
        );
    }
}
#[test]
fn actual_strict_shoulder_thresholds_preserve_one_ulp_positive_material() {
    let spec = fixture(false, false);
    for radius in [
        spec.metric.maximum_radius_mm(),
        spec.wood.minimum_radius_mm(),
    ] {
        assert_eq!(
            Hanger::new(
                Spec {
                    plain_radius_mm: radius,
                    ..spec
                },
                None
            ),
            Err(Error::Unsupported)
        );
    }
    for radius in [
        spec.metric.maximum_radius_mm().next_up(),
        spec.wood.minimum_radius_mm().next_down(),
    ] {
        let source = Hanger::new(
            Spec {
                plain_radius_mm: radius,
                ..spec
            },
            None,
        )
        .unwrap();
        assert_eq!(source.spec().plain_radius_mm, radius);
        assert!(source.volume(options(), None).unwrap().material_mm3[0] > 0.);
    }
    for value in [0., -1., f64::NAN, f64::INFINITY, 1e7 + 1.] {
        for field in 0..2 {
            let mut invalid = spec;
            if field == 0 {
                invalid.plain_length_mm = value
            } else {
                invalid.plain_radius_mm = value
            };
            assert_eq!(Hanger::new(invalid, None), Err(Error::Invalid));
        }
    }
}
#[test]
fn query_failures_cancelled_publication_and_retries_preserve_original_sources() {
    use std::cell::Cell;
    let spec = fixture(false, true);
    let source = Hanger::new(spec, None).unwrap();
    assert_eq!(Hanger::new(spec, Some(&|| true)), Err(Error::Cancelled));
    let calls = Cell::new(0);
    assert_eq!(
        Hanger::new(
            spec,
            Some(&|| {
                calls.set(calls.get() + 1);
                calls.get() > 1
            })
        ),
        Err(Error::Cancelled)
    );
    for value in [0., -1., f64::NAN, f64::INFINITY] {
        assert_eq!(
            source.volume(
                QueryOptions {
                    volume_error_mm3: value,
                    ..options()
                },
                None
            ),
            Err(Error::Invalid)
        );
    }
    for cells in [0, 1, 5, 10, 65537, usize::MAX] {
        assert_eq!(
            source.volume(
                QueryOptions {
                    maximum_cells: cells,
                    ..options()
                },
                None
            ),
            Err(Error::Limit)
        );
    }
    assert_eq!(
        source.volume(
            QueryOptions {
                volume_error_mm3: 1e-30,
                ..options()
            },
            None
        ),
        Err(Error::Unresolved)
    );
    for stop in [0, 3, 12] {
        let calls = Cell::new(0);
        assert_eq!(
            source.volume(
                options(),
                Some(&|| {
                    calls.set(calls.get() + 1);
                    calls.get() > stop
                })
            ),
            Err(Error::Cancelled)
        );
    }
    assert_eq!(
        source.classify_box(point(f64::NAN, 0., 0.), None),
        Err(Error::Invalid)
    );
    assert_eq!(
        source.classify_box(point(0., 0., 0.), Some(&|| true)),
        Err(Error::Cancelled)
    );
    let calls = Cell::new(0);
    assert_eq!(
        source.classify_box(
            point(0., 0., 0.),
            Some(&|| {
                calls.set(calls.get() + 1);
                calls.get() > 1
            })
        ),
        Err(Error::Cancelled)
    );
    let calls = Cell::new(0);
    assert_eq!(
        source.classify_box(
            point(1.3, 0., -10.),
            Some(&|| {
                calls.set(calls.get() + 1);
                calls.get() > 1
            })
        ),
        Err(Error::Cancelled)
    );
    assert_eq!(source.spec(), spec);
    assert!(source.volume(options(), None).is_ok());
}
#[test]
fn independent_hands_and_starts_change_phase_but_not_complete_angular_volume() {
    let spec = fixture(false, false);
    let baseline = Hanger::new(spec, None)
        .unwrap()
        .volume(options(), None)
        .unwrap();
    let metric = AuthoredMetricThread::new(MetricSpec {
        starts: 2,
        left_hand: true,
        ..spec.metric.spec()
    })
    .unwrap();
    let wood = AuthoredThread::new(
        WoodSpec {
            starts: 2,
            ..spec.wood.spec()
        },
        SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-6,
        },
    )
    .unwrap();
    let source = Hanger::new(
        Spec {
            wood,
            metric,
            ..spec
        },
        None,
    )
    .unwrap();
    assert_eq!(*source.metric(), metric);
    assert_eq!(*source.wood(), wood);
    assert_eq!(source.volume(options(), None).unwrap(), baseline);
    assert_eq!(
        source.classify_box(point(-1.3, 0., -13.), None),
        Ok(Classification::Inside)
    );
}

#[test]
fn actual_canonical_closed_candidate_preserves_connected_source_and_both_volumes() {
    let source = Hanger::new(fixture(false, false), None).unwrap();
    let now = std::time::Instant::now();
    let c = source
        .tessellate_candidate(
            crate::display::Tessellation {
                deflection: 0.1,
                angle: 0.3,
                triangle_limit: 1_000_000,
            },
            16_000_000,
            None,
        )
        .unwrap();
    eprintln!(
        "closedhanger{}facets,{}work,{:?},chord{},normal{},source{:?},literal{:?}",
        c.mesh().face_indices.len(),
        c.work_used(),
        now.elapsed(),
        c.source_chord_bound_mm(),
        c.source_normal_angle_bound_rad(),
        c.source_volume_mm3(),
        c.literal_volume_mm3()
    );
    assert_eq!(c.source(), &source);
    assert!(c.source_chord_bound_mm() <= 0.1 && c.source_normal_angle_bound_rad() <= 0.3);
    assert!(c.literal_volume_mm3()[0] > 0.);
    for &face in &c.mesh().face_indices {
        assert!(source.face_key(face as usize).is_some());
    }
}
#[test]
fn actual_nonbinary_left_hand_closed_candidate_keeps_full_source_and_quota() {
    let source = Hanger::new(fixture(true, true), None).unwrap();
    let now = std::time::Instant::now();
    let c = source
        .tessellate_candidate(
            crate::display::Tessellation {
                deflection: 0.1,
                angle: 0.3,
                triangle_limit: 1_000_000,
            },
            16_000_000,
            None,
        )
        .unwrap();
    eprintln!(
        "closedLHhanger{}facets,{}work,{:?},chord{},normal{}",
        c.mesh().face_indices.len(),
        c.work_used(),
        now.elapsed(),
        c.source_chord_bound_mm(),
        c.source_normal_angle_bound_rad()
    );
    assert_eq!(c.source(), &source);
    assert!(c.work_used() <= 16_000_000);
}

#[test]
fn exact_periodic_wood_domain_requires_a_separate_connected_cap_contract() {
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
