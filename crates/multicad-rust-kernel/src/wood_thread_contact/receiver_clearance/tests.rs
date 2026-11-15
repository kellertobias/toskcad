use super::*;
fn identity() -> Frame {
    Frame {
        origin_mm: [0.; 3],
        axes: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
    }
}
fn profile(width_scale: f64) -> Sweep {
    let mut side = vec![[4., 0.], [16., 0.], [16., 4.5]];
    let flank = [
        [9.25 * width_scale, 4.5],
        [7.75 * width_scale, 8.],
        [5.75 * width_scale, 10.],
        [2.75 * width_scale, 12.25],
    ];
    side.extend(flank.map(|[w, d]| [20. - w, d]));
    side.extend(flank.iter().rev().map(|&[w, d]| [20. + w, d]));
    side.extend([[24., 4.5], [24., 0.], [36., 0.]]);
    let rotate = |[x, y]: [f64; 2]| [40. - y, x];
    let mut center = [36., 4.];
    let mut outer = Vec::new();
    for _ in 0..4 {
        outer.extend(side.windows(2).map(|p| Edge::Line { a: p[0], b: p[1] }));
        outer.push(Edge::Arc {
            center,
            a: side[side.len() - 1],
            b: rotate(side[0]),
            clockwise: false,
        });
        for p in &mut side {
            *p = rotate(*p);
        }
        center = rotate(center);
    }
    Sweep::new(
        vec![
            outer,
            vec![Edge::Circle {
                center: [20.; 2],
                radius_mm: 3.4,
                clockwise: true,
            }],
        ],
        100.,
        1e-7,
        1e-5,
    )
    .unwrap()
}
fn nut() -> crate::line8_nut::Nut {
    use crate::{
        authored_flat_thread::{AuthoredThread, SourceOptions},
        helical_thread::Spec,
        linear_sweep::Sweep,
    };
    let section = vec![
        [-6.9, 0.],
        [6.9, 0.],
        [6.5, 3.],
        [5.4, 5.],
        [3.3, 7.3],
        [-3.3, 7.3],
        [-5.4, 5.],
        [-6.5, 3.],
    ];
    let stock = Sweep::new(vec![section.clone()], vec![section], 23.).unwrap();
    let tap = AuthoredThread::new(
        Spec {
            major_radius_mm: 4.08,
            pitch_mm: 1.25,
            length_mm: 4.1,
            depth_mm: 0.65,
            crest_width_mm: 0.15,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            starts: 1,
            left_hand: false,
            root_radius_mm: None,
        },
        SourceOptions {
            position_limit_mm: 1e-8,
            normal_angle_limit_rad: 1e-8,
        },
    )
    .unwrap();
    crate::line8_nut::Nut::new(
        stock,
        7.5,
        tap,
        [[1., 0.], [0., 1.]],
        crate::line8_nut::Relief {
            start_mm: 4.099,
            radius_mm: 4.1,
            tool_length_mm: 3.22,
        },
        None,
    )
    .unwrap()
}
fn stock_frame() -> Frame {
    Frame {
        origin_mm: [20., 35.5, 32.5],
        axes: [[0., 0., 1.], [-1., 0., 0.], [0., -1., 0.]],
    }
}
fn options() -> ReceiverClearanceOptions {
    ReceiverClearanceOptions {
        maximum_overlap_volume_mm3: 0.001,
        maximum_work: 1_000_000,
    }
}
#[test]
fn original_actual_groove_and_both_members_earn_void_and_outer_source_proofs() {
    let source = nut();
    let member = profile(1.);
    let placed = ClockedLine8Placement {
        source: &source,
        stock_frame: stock_frame(),
        outer_frame: identity(),
        minimum_wall_mm: 0.5,
    };
    let receiver = ReceivingSweep {
        source: &member,
        frame: identity(),
        occurrence: identity(),
        reverse: false,
    };
    let first = certify_line8_receiver_clearance(placed, receiver, options(), &|| false).unwrap();
    assert!(first.used_exterior_cell);
    assert_eq!(first.overlap_volume_upper_mm3, 0.);
    assert!(first.boundary_owners.len() >= 7);
    assert!(std::ptr::eq(first.nut().source, &source));
    assert!(std::ptr::eq(first.receiver().source, &member));
    let other = ReceivingSweep {
        occurrence: Frame {
            origin_mm: [0., 40., 0.],
            ..identity()
        },
        ..receiver
    };
    let second = certify_line8_receiver_clearance(placed, other, options(), &|| false).unwrap();
    assert!(!second.used_exterior_cell);
    assert_eq!(second.overlap_volume_upper_mm3, 0.);
    let altered = profile(0.7);
    assert!(matches!(
        certify_line8_receiver_clearance(
            placed,
            ReceivingSweep {
                source: &altered,
                ..receiver
            },
            options(),
            &|| false
        ),
        Err(Error::Unsupported)
    ));
    let moved = ClockedLine8Placement {
        stock_frame: Frame {
            origin_mm: [20., 35.55, 32.5],
            ..stock_frame()
        },
        ..placed
    };
    assert!(matches!(
        certify_line8_receiver_clearance(moved, receiver, options(), &|| false),
        Err(Error::Unsupported)
    ));
}
#[test]
fn original_occurrence_operators_and_actual_rounding_residual_keep_volume_threshold() {
    let source = nut();
    let member = profile(1.);
    let pose = Frame {
        origin_mm: [37.1, -11.3, 103.7],
        axes: [[0.8, 0.6, 0.], [-0.6, 0.8, 0.], [0., 0., 1.]],
    };
    let placed = ClockedLine8Placement {
        source: &source,
        stock_frame: stock_frame(),
        outer_frame: pose,
        minimum_wall_mm: 0.5,
    };
    let receiver = ReceivingSweep {
        source: &member,
        frame: identity(),
        occurrence: pose,
        reverse: false,
    };
    let exact = certify_line8_receiver_clearance(placed, receiver, options(), &|| false).unwrap();
    assert_eq!(exact.overlap_volume_upper_mm3, 0.);
    let apply = |p: [f64; 3], translation: bool| {
        std::array::from_fn(|i| {
            (if translation { pose.origin_mm[i] } else { 0. })
                + (0..3).map(|j| pose.axes[j][i] * p[j]).sum::<f64>()
        })
    };
    let raw = stock_frame();
    let rounded = Frame {
        origin_mm: apply(raw.origin_mm, true),
        axes: raw.axes.map(|a| apply(a, false)),
    };
    let literal = ClockedLine8Placement {
        stock_frame: rounded,
        outer_frame: identity(),
        ..placed
    };
    let measured =
        certify_line8_receiver_clearance(literal, receiver, options(), &|| false).unwrap();
    eprintln!(
        "Actual rounded occurrence overlap upper={} mm3; charged_work={}; peak_bytes={}",
        measured.overlap_volume_upper_mm3, measured.work_used, measured.accounted_peak_bytes
    );
    assert!(measured.overlap_volume_upper_mm3 > 0.);
    assert!(measured.overlap_volume_upper_mm3 < 1e-8);
    assert!(
        certify_line8_receiver_clearance(
            literal,
            receiver,
            ReceiverClearanceOptions {
                maximum_overlap_volume_mm3: measured.overlap_volume_upper_mm3,
                ..options()
            },
            &|| false
        )
        .is_ok()
    );
    assert!(matches!(
        certify_line8_receiver_clearance(
            literal,
            receiver,
            ReceiverClearanceOptions {
                maximum_overlap_volume_mm3: 0.,
                ..options()
            },
            &|| false
        ),
        Err(Error::Unsupported)
    ));
    assert!(matches!(
        certify_line8_receiver_clearance(placed, receiver, options(), &|| true),
        Err(Error::Cancelled)
    ));
    assert!(matches!(
        certify_line8_receiver_clearance(
            placed,
            receiver,
            ReceiverClearanceOptions {
                maximum_work: 1,
                ..options()
            },
            &|| false
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        certify_line8_receiver_clearance(
            placed,
            receiver,
            ReceiverClearanceOptions {
                maximum_work: 0,
                ..options()
            },
            &|| false
        ),
        Err(Error::Invalid)
    ));
}

#[test]
fn real_arc_incidence_actual_capacity_and_mid_query_cancellation_remain_bounded() {
    let source = nut();
    let member = profile(1.);
    let placed = ClockedLine8Placement {
        source: &source,
        stock_frame: stock_frame(),
        outer_frame: identity(),
        minimum_wall_mm: 0.5,
    };
    let receiver = ReceivingSweep {
        source: &member,
        frame: identity(),
        occurrence: identity(),
        reverse: false,
    };
    let mut rings = member.rings().to_vec();
    let Edge::Arc {
        center,
        a,
        mut b,
        clockwise,
    } = rings[0][13]
    else {
        unreachable!()
    };
    b[0] += 1e-9;
    rings[0][13] = Edge::Arc {
        center,
        a,
        b,
        clockwise,
    };
    let imperfect = Sweep::new(rings, 100., 1e-7, 1e-5).unwrap();
    assert!(matches!(
        certify_line8_receiver_clearance(
            placed,
            ReceivingSweep {
                source: &imperfect,
                ..receiver
            },
            options(),
            &|| false
        ),
        Err(Error::Unsupported)
    ));
    let mut rings = member.rings().to_vec();
    let capacity = 48 * 1024 * 1024 / std::mem::size_of::<Edge>() + 1024;
    let mut original_capacity = Vec::with_capacity(capacity);
    original_capacity.extend_from_slice(&rings[0]);
    rings[0] = original_capacity;
    let large = Sweep::new(rings, 100., 1e-7, 1e-5).unwrap();
    assert!(large.retained_bytes() > 48 * 1024 * 1024);
    assert!(matches!(
        certify_line8_receiver_clearance(
            placed,
            ReceivingSweep {
                source: &large,
                ..receiver
            },
            options(),
            &|| false
        ),
        Err(Error::Limit)
    ));
    let total = std::cell::Cell::new(0usize);
    certify_line8_receiver_clearance(placed, receiver, options(), &|| {
        total.set(total.get() + 1);
        false
    })
    .unwrap();
    assert!(total.get() > 10);
    for stop in [1, total.get() / 4, total.get() / 2, total.get() - 1] {
        let calls = std::cell::Cell::new(0);
        let result = certify_line8_receiver_clearance(placed, receiver, options(), &|| {
            calls.set(calls.get() + 1);
            calls.get() > stop
        });
        assert!(matches!(result, Err(Error::Cancelled)), "checkpoint {stop}");
    }
    let retry = certify_line8_receiver_clearance(placed, receiver, options(), &|| false).unwrap();
    assert_eq!(retry.overlap_volume_upper_mm3, 0.);
    assert!(retry.used_exterior_cell);
    assert!(matches!(
        certify_line8_receiver_clearance(
            placed,
            receiver,
            ReceiverClearanceOptions {
                maximum_overlap_volume_mm3: f64::NAN,
                ..options()
            },
            &|| false
        ),
        Err(Error::Invalid)
    ));
}

#[test]
fn actual_circular_cover_and_complex_retained_section_do_not_use_family_labels() {
    let source = nut();
    let placed = ClockedLine8Placement {
        source: &source,
        stock_frame: stock_frame(),
        outer_frame: identity(),
        minimum_wall_mm: 0.5,
    };
    let circle = Sweep::new(
        vec![vec![Edge::Circle {
            center: [0.; 2],
            radius_mm: 4.,
            clockwise: false,
        }]],
        100.,
        1e-7,
        1e-5,
    )
    .unwrap();
    let receiver = ReceivingSweep {
        source: &circle,
        frame: identity(),
        occurrence: identity(),
        reverse: false,
    };
    let proof = certify_line8_receiver_clearance(placed, receiver, options(), &|| false).unwrap();
    assert_eq!(proof.overlap_volume_upper_mm3, 0.);
    assert!(!proof.used_exterior_cell);
    let points = (0..300)
        .map(|i| {
            let (s, c) = (std::f64::consts::TAU * i as f64 / 300.).sin_cos();
            [20. * c, 20. * s]
        })
        .collect::<Vec<_>>();
    let outer = (0..points.len())
        .map(|i| Edge::Line {
            a: points[i],
            b: points[(i + 1) % points.len()],
        })
        .collect();
    let complex = Sweep::new(vec![outer], 100., 1e-7, 1e-5).unwrap();
    assert!(matches!(
        certify_line8_receiver_clearance(
            placed,
            ReceivingSweep {
                source: &complex,
                ..receiver
            },
            options(),
            &|| false
        ),
        Err(Error::Limit)
    ));
}

#[test]
fn actual_reversed_extent_and_filled_collinear_section_keep_material_semantics() {
    let source = nut();
    let member = profile(1.);
    let placed = ClockedLine8Placement {
        source: &source,
        stock_frame: stock_frame(),
        outer_frame: identity(),
        minimum_wall_mm: 0.5,
    };
    let receiver = ReceivingSweep {
        source: &member,
        frame: identity(),
        occurrence: identity(),
        reverse: true,
    };
    assert_eq!(
        certify_line8_receiver_clearance(placed, receiver, options(), &|| false)
            .unwrap()
            .overlap_volume_upper_mm3,
        0.
    );
    let corners = [[0., 0.], [40., 0.], [40., 40.], [0., 40.]];
    let mut points = Vec::new();
    for i in 0..4 {
        for n in 0..3 {
            points.push(std::array::from_fn(|j| {
                corners[i][j] + (corners[(i + 1) % 4][j] - corners[i][j]) * n as f64 / 3.
            }));
        }
    }
    let edges = (0..points.len())
        .map(|i| Edge::Line {
            a: points[i],
            b: points[(i + 1) % points.len()],
        })
        .collect();
    let filled = Sweep::new(vec![edges], 100., 1e-7, 1e-5).unwrap();
    assert!(matches!(
        certify_line8_receiver_clearance(
            placed,
            ReceivingSweep {
                source: &filled,
                reverse: false,
                ..receiver
            },
            options(),
            &|| false
        ),
        Err(Error::Unsupported)
    ));
}
