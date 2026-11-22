use super::*;
fn identity() -> Frame {
    Frame {
        origin_mm: [0.; 3],
        axes: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
    }
}
fn cross_f64(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn transported(frame: Frame, outer: Frame) -> Frame {
    let origin_mm = std::array::from_fn(|i| {
        outer.origin_mm[i]
            + (0..3)
                .map(|j| outer.axes[j][i] * frame.origin_mm[j])
                .sum::<f64>()
    });
    let u = std::array::from_fn(|i| (0..3).map(|j| outer.axes[j][i] * frame.axes[0][j]).sum());
    let v = std::array::from_fn(|i| (0..3).map(|j| outer.axes[j][i] * frame.axes[1][j]).sum());
    Frame {
        origin_mm,
        axes: [u, v, cross_f64(u, v)],
    }
}
fn fixture(left: bool, starts: u32, nb: bool) -> (Hanger, Nut, Frame) {
    let mut spec = crate::authored_hanger_bolt::tests::fixture(left, nb);
    let mut ms = spec.metric.spec();
    ms.starts = starts;
    spec.metric = AuthoredMetricThread::new(ms).unwrap();
    let source = Hanger::new(spec, None).unwrap();
    let pitch = 0.5;
    let nut = Nut::new(
        Thread::new(crate::helical_thread::Spec {
            major_radius_mm: 1.5,
            pitch_mm: pitch,
            length_mm: 2.4,
            depth_mm: (3. - (2.459 + 2.599) / 2.) / 2.,
            crest_width_mm: pitch / 2. - (3. - (2.675 + 2.775) / 2.) / 3_f64.sqrt(),
            upper_flank_deg: 30.,
            lower_flank_deg: 30.,
            starts,
            left_hand: left,
            root_radius_mm: None,
        })
        .unwrap(),
        5.5,
    )
    .unwrap();
    // The normal corefixture has grip4.1, washer.2, and nut2.4. These exact
    // actual source placements deliberately produce a non-quarter rotation.
    let nut_z = -4.1 - 0.2 - 2.4;
    let metric_start = -source.spec().plain_length_mm - source.metric().full_length_mm();
    let angle = std::f64::consts::TAU * (nut_z - metric_start) / (pitch * f64::from(starts))
        * if left { -1. } else { 1. };
    let (s, c) = angle.portable_sin_cos();
    let axes = [[c, s, 0.], [-s, c, 0.], [0., 0., c * c + s * s]];
    (
        source,
        nut,
        Frame {
            origin_mm: [0., 0., nut_z],
            axes,
        },
    )
}
fn req() -> HangerMatingRequirements {
    HangerMatingRequirements {
        minimum_engagement_mm: 2.,
        ..Default::default()
    }
}
#[test]
fn actual_nonquarter_nut_mates_with_continuous_metric_for_both_hands_starts_and_proper_far_pose() {
    for (left, starts, nb, pose) in [
        (false, 1, false, 0),
        (true, 1, true, 0),
        (false, 2, true, 1),
        (true, 2, false, 2),
    ] {
        let (source, nut, mut nut_frame) = fixture(left, starts, nb);
        let u = [0.36, 0.48, 0.8];
        let v = [-0.8, 0.6, 0.];
        let frame = if pose == 0 {
            identity()
        } else {
            Frame {
                origin_mm: if pose == 1 {
                    [11., -7., 5.]
                } else {
                    [1e6, -2e6, 3e6]
                },
                axes: [u, v, cross_f64(u, v)],
            }
        };
        if pose != 0 {
            nut_frame = transported(nut_frame, frame);
        }
        let result = certify_hanger_nut(&source, frame, &nut, nut_frame, req(), &|| false).unwrap();
        assert!(result.engagement_mm[0] >= 2.);
        assert!(result.minimum_radial_clearance_lower_mm >= 1e-7);
        assert!(result.proved_phase_cells > 0);
        assert!(result.phase_cells_used >= result.proved_phase_cells);
        assert!(result.phase_cells_used <= 4096);
        assert!(result.relative_angle_error_rad < 1e-6);
        assert!(result.plain_shoulder_clearance_mm[0] > 0.);
    }
}
#[test]
fn actual_retained_round_root_nut_support_is_not_replaced_by_a_flat_profile() {
    let (source, _, frame) = fixture(false, 1, false);
    let e = source.metric().spec();
    let radius = e.pitch_mm / (4. * 3_f64.sqrt());
    let thread = Thread::new(crate::helical_thread::Spec {
        major_radius_mm: e.major_radius_mm + 0.2,
        pitch_mm: e.pitch_mm,
        length_mm: 2.4,
        depth_mm: (e.pitch_mm - e.crest_width_mm) / (2. / 3_f64.sqrt()) - radius,
        crest_width_mm: e.crest_width_mm,
        upper_flank_deg: 30.,
        lower_flank_deg: 30.,
        starts: 1,
        left_hand: false,
        root_radius_mm: Some(radius),
    })
    .unwrap();
    let nut = Nut::new(thread, 5.5).unwrap();
    let mut requirements = req();
    requirements.minimum_radial_clearance_mm = 0.05;
    let proof =
        certify_hanger_nut(&source, identity(), &nut, frame, requirements, &|| false).unwrap();
    assert!(proof.minimum_radial_clearance_lower_mm >= 0.05);
    assert!(nut.thread().root_circle().is_some());
}

#[test]
fn incompatible_laws_real_interference_finite_extent_and_resource_limits_refuse() {
    let (source, nut, frame) = fixture(false, 1, false);
    let mut cases = Vec::new();
    let mut displaced = frame;
    displaced.origin_mm[0] = 0.001;
    cases.push(displaced);
    let mut wrong_roll = frame;
    wrong_roll.axes = identity().axes;
    cases.push(wrong_roll);
    let mut past_tip = frame;
    past_tip.origin_mm[2] = -20.;
    cases.push(past_tip);
    let mut crosses_plain = frame;
    crosses_plain.origin_mm[2] = -3.;
    cases.push(crosses_plain);
    for f in cases {
        assert!(matches!(
            certify_hanger_nut(&source, identity(), &nut, f, req(), &|| false),
            Err(Error::Unsupported | Error::Limit)
        ));
    }
    for (pitch, starts, left_hand) in [(0.6, 1, false), (0.5, 2, false), (0.5, 1, true)] {
        let mut spec = nut.thread().spec();
        spec.pitch_mm = pitch;
        spec.starts = starts;
        spec.left_hand = left_hand;
        let other = Nut::new(Thread::new(spec).unwrap(), 5.5).unwrap();
        assert!(matches!(
            certify_hanger_nut(&source, identity(), &other, frame, req(), &|| false),
            Err(Error::Unsupported)
        ));
    }
    let mut r = req();
    r.minimum_engagement_mm = 3.;
    assert!(matches!(
        certify_hanger_nut(&source, identity(), &nut, frame, r, &|| false),
        Err(Error::Unsupported)
    ));
    r = req();
    r.minimum_radial_clearance_mm = 0.5;
    assert!(matches!(
        certify_hanger_nut(&source, identity(), &nut, frame, r, &|| false),
        Err(Error::Unsupported | Error::Limit)
    ));
    r = req();
    r.maximum_cells = 1;
    assert!(matches!(
        certify_hanger_nut(&source, identity(), &nut, frame, r, &|| false),
        Err(Error::Limit)
    ));
    for n in [0, 4097] {
        r = req();
        r.maximum_cells = n;
        assert!(matches!(
            certify_hanger_nut(&source, identity(), &nut, frame, r, &|| false),
            Err(Error::Invalid)
        ));
    }
    r = req();
    r.minimum_engagement_mm = f64::NAN;
    assert!(matches!(
        certify_hanger_nut(&source, identity(), &nut, frame, r, &|| false),
        Err(Error::Invalid)
    ));
    let mut bad = frame;
    bad.axes[2][2] = -1.;
    assert!(matches!(
        certify_hanger_nut(&source, identity(), &nut, bad, req(), &|| false),
        Err(Error::Unsupported | Error::Invalid)
    ));
}
#[test]
fn every_cancellation_checkpoint_refuses_and_source_retry_is_unchanged() {
    use std::cell::Cell;
    let (source, nut, frame) = fixture(false, 1, false);
    let calls = Cell::new(0);
    let proof = certify_hanger_nut(&source, identity(), &nut, frame, req(), &|| {
        calls.set(calls.get() + 1);
        false
    })
    .unwrap();
    let total = calls.get();
    assert!(total > proof.phase_cells_used);
    for stop in 1..=total {
        let count = Cell::new(0);
        assert!(matches!(
            certify_hanger_nut(&source, identity(), &nut, frame, req(), &|| {
                count.set(count.get() + 1);
                count.get() == stop
            }),
            Err(Error::Cancelled)
        ));
    }
    let again = certify_hanger_nut(&source, identity(), &nut, frame, req(), &|| false).unwrap();
    assert_eq!(again.source_phase_offset_mm, proof.source_phase_offset_mm);
    assert_eq!(
        again.minimum_radial_clearance_lower_mm,
        proof.minimum_radial_clearance_lower_mm
    );
}

#[test]
fn bounded_enclosure_primitives_refuse_uncertainty_without_changing_supports() {
    assert!(matches!(
        ray_angle(I::new(-1., 1.), ip(1.), I::new(1., 2.)),
        Err(Error::Unsupported)
    ));
    assert!(matches!(
        ray_angle(ip(f64::NAN), ip(1.), ip(1.)),
        Err(Error::Unsupported)
    ));
    assert!(matches!(
        phase_midpoint(0.5, 0.5_f64.next_up()),
        Err(Error::Unsupported)
    ));
    assert_eq!(phase_midpoint(0., 0.5).unwrap(), 0.25);
    for (e, want) in [
        (
            crate::authored_metric_thread::Error::Invalid,
            Error::Invalid,
        ),
        (
            crate::authored_metric_thread::Error::Cancelled,
            Error::Cancelled,
        ),
        (crate::authored_metric_thread::Error::Limit, Error::Limit),
        (
            crate::authored_metric_thread::Error::Unresolved,
            Error::Unsupported,
        ),
    ] {
        assert_eq!(metric_error(e), want);
    }
    let (_, nut, _) = fixture(false, 1, false);
    assert!(matches!(
        retained_periodic(nut.thread(), ip(3000.), 0.),
        Err(Error::Unsupported)
    ));
    let (source, nut, frame) = fixture(false, 1, false);
    let mut small = nut.thread().spec();
    small.major_radius_mm = 0.3;
    let small = Nut::new(Thread::new(small).unwrap(), 5.5).unwrap();
    assert!(matches!(
        certify_hanger_nut(&source, identity(), &small, frame, req(), &|| false),
        Err(Error::Unsupported)
    ));
    let proof = certify_hanger_nut(&source, identity(), &nut, frame, req(), &|| false).unwrap();
    for cells in 1..proof.phase_cells_used {
        let mut r = req();
        r.maximum_cells = cells;
        assert!(matches!(
            certify_hanger_nut(&source, identity(), &nut, frame, r, &|| false),
            Err(Error::Limit)
        ));
    }
    let p = nut.thread().spec().pitch_mm;
    let radius = 1e-14;
    let mut rs = nut.thread().spec();
    rs.root_radius_mm = Some(radius);
    rs.depth_mm = (p - rs.crest_width_mm) / (2. / 3_f64.sqrt()) - radius;
    let rounded = Thread::new(rs).unwrap();
    let rounded_nut = Nut::new(rounded, 5.5).unwrap();
    assert!(matches!(
        certify_hanger_nut(&source, identity(), &rounded_nut, frame, req(), &|| false),
        Err(Error::Unsupported)
    ));
    assert!(matches!(
        retained_radius(&rounded, I::new(0., p), 1e-7),
        Err(Error::Unsupported)
    ));
}

#[test]
fn actual_near_axis_metric_requires_a_resolved_angular_annulus_even_with_loose_intent() {
    let (source, nut, mut frame) = fixture(false, 1, false);
    let mut spec = source.spec();
    let mut metric = spec.metric.spec();
    let root = metric.pitch_mm / (4. * 3_f64.sqrt());
    metric.major_radius_mm =
        (metric.pitch_mm - metric.crest_width_mm) * 3_f64.sqrt() / 2. - root + 0.0001;
    spec.metric = AuthoredMetricThread::new(metric).unwrap();
    let source = Hanger::new(spec, None).unwrap();
    frame.origin_mm[0] = 0.00009;
    let mut requirements = req();
    requirements.linear_tolerance_mm = 0.001;
    assert!(source.metric().minimum_radius_mm()[1] < 0.00011);
    assert!(matches!(
        certify_hanger_nut(&source, identity(), &nut, frame, requirements, &|| false),
        Err(Error::Unsupported)
    ));
}

#[test]
fn phase_publication_keeps_all_numerical_resource_guards() {
    use num_bigint::BigInt;
    let huge = Q::from_integer(BigInt::from(1_u8) << 17000);
    assert!(matches!(bounded(huge), Err(Error::Limit)));
    let mut map = Alignment {
        scale: ip(1.),
        angle: 0.,
        angle_error: 0.,
        radial: 0.,
        axial: 0.,
        origin_z: ip(1.),
    };
    let allowance = phase_error::Allowance {
        indexed_mm: 0.,
        evaluator_mm: 0.,
        total_mm: 0.,
    };
    assert!(phase_offset(&map, ip(0.1), ip(0.1), ip(0.5), false, allowance).is_ok());
    map.axial = f64::INFINITY;
    assert!(matches!(
        phase_offset(&map, ip(0.1), ip(0.1), ip(0.5), false, allowance),
        Err(Error::Unsupported)
    ));
    map.axial = 0.;
    map.origin_z = ip(f64::INFINITY);
    assert!(matches!(
        phase_offset(&map, ip(0.1), ip(0.1), ip(0.5), false, allowance),
        Err(Error::Unsupported)
    ));
}
#[test]
fn earned_source_domains_keep_exact_support_and_checked_gap_oracles() {
    let (source, nut, frame) = fixture(false, 1, false);
    let tiny = f64::from_bits(1);
    let u = [1., tiny, tiny];
    let v = [-tiny, 1., tiny];
    let from = Frame {
        origin_mm: [1e7, -1e7, 1e7],
        axes: [u, v, cross_f64(u, v)],
    };
    let to = Frame {
        origin_mm: [-1e7, 1e7, -1e7],
        ..identity()
    };
    scale(from).unwrap();
    scale(to).unwrap();
    let (o, a) = affine(to, from, &|| false).unwrap();
    let large = Nut::new(*nut.thread(), 1e7).unwrap();
    let shift =
        -q(source.spec().plain_length_mm).unwrap() - q(source.metric().full_length_mm()).unwrap();
    for xy in crate::headed_screw::hex_vertices(large.across_flats_mm()) {
        for z in [0., large.thread().spec().length_mm] {
            let xyz = [q(xy[0]).unwrap(), q(xy[1]).unwrap(), q(z).unwrap()];
            let actual = axial_corner(&o[2], [&a[0][2], &a[1][2], &a[2][2]], xyz.clone(), &shift);
            let oracle =
                &o[2] + &a[0][2] * &xyz[0] + &a[1][2] * &xyz[1] + &a[2][2] * &xyz[2] - &shift;
            assert_eq!(actual, oracle);
            assert!(actual.numer().bits() + actual.denom().bits() < 8704);
            assert!(bounded(actual).is_ok());
        }
    }
    for perturb in [false, true] {
        let mut placed = frame;
        if perturb {
            placed.axes[0][0] = placed.axes[0][0].next_up();
        }
        let metric = scale(placed).unwrap();
        let map = align(&source, identity(), placed, &|| false).unwrap();
        let allowance = phase_error::measure(nut.thread(), &|| false).unwrap();
        for k in 0..32 {
            let phase = I::new(f64::from(k) * 0.5 / 32., f64::from(k + 1) * 0.5 / 32.);
            let external = continuous_radius(source.metric(), phase).unwrap();
            let inner = retained_periodic(nut.thread(), phase, allowance.total_mm).unwrap();
            let low =
                finite(ip(inner.lo) - ip(map.scale.hi) * ip(external.hi) - ip(map.radial)).unwrap();
            let high =
                finite(ip(inner.hi) - ip(map.scale.lo) * ip(external.lo) + ip(map.radial)).unwrap();
            let oracle = finite(
                I::new(low.lo, high.hi) * I::new(metric.minimum_scale.lo, metric.maximum_scale.hi),
            )
            .unwrap();
            let actual = physical_gap(inner, external, &map, &metric);
            assert_eq!([actual.lo, actual.hi], [oracle.lo, oracle.hi]);
            assert!(actual.checked().is_ok());
        }
    }
}
