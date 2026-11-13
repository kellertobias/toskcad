use super::*;
use crate::{authored_flat_thread::SourceOptions, helical_thread::Spec};
fn identity() -> Frame {
    Frame {
        origin_mm: [0.; 3],
        axes: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
    }
}
fn fixture(left: bool, rounded: bool) -> (Insert, Thread) {
    let major = (2.874 + 2.980) / 4.;
    let pitch_mid = (2.580 + 2.655) / 2.;
    let pitch = 0.5;
    let crest = pitch / 2. - (2. * major - pitch_mid) / 3_f64.sqrt();
    let r = pitch / (4. * 3_f64.sqrt());
    let outer = Spec {
        major_radius_mm: major,
        pitch_mm: pitch,
        length_mm: 4.,
        depth_mm: (pitch - crest) / (2. / 3_f64.sqrt()) - r,
        crest_width_mm: crest,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: left,
        root_radius_mm: rounded.then_some(r),
    };
    let source = |s| {
        AuthoredThread::new(
            s,
            SourceOptions {
                position_limit_mm: 1e-8,
                normal_angle_limit_rad: 1e-6,
            },
        )
        .unwrap()
    };
    let inner = source(Spec {
        major_radius_mm: 1.5,
        depth_mm: (3. - (2.459 + 2.599) / 2.) / 2.,
        crest_width_mm: pitch / 2. - (3. - (2.675 + 2.775) / 2.) / 3_f64.sqrt(),
        root_radius_mm: None,
        ..outer
    });
    let wood = source(Spec {
        major_radius_mm: 5.,
        pitch_mm: 2.,
        depth_mm: 1.,
        crest_width_mm: 0.4,
        root_radius_mm: None,
        ..outer
    });
    (
        Insert::new(wood, inner, 0.5, None).unwrap(),
        Thread::new(outer).unwrap(),
    )
}
fn run(
    insert: &Insert,
    external: &Thread,
    frame: Frame,
    visible: [f64; 2],
    prefix: f64,
) -> R<InsertMating> {
    assert_eq!(visible, [0., external.spec().length_mm]);
    assert_eq!(prefix, 0.);
    let screw = crate::headed_screw::Screw::new(
        *external,
        crate::headed_screw::ScrewSpec {
            head: crate::headed_screw::ScrewHead::SocketCap {
                radius_mm: 2.75,
                height_mm: 3.,
                socket_across_flats_mm: 2.,
                socket_depth_mm: 1.,
            },
            shank_length_mm: 0.,
            shank_radius_mm: 1.5,
        },
    )
    .unwrap();
    certify_insert_mating(
        insert,
        identity(),
        MatingExternal {
            source: MatingScrew::Regular(&screw),
            frame,
        },
        MatingRequirements {
            minimum_engagement_mm: 4.,
            minimum_clearance_mm: 4e-7,
        },
        &|| false,
    )
}

#[test]
fn actual_m3_6g_6h_supports_clear_at_full_four_mm_exact_threshold() {
    for left in [false, true] {
        for rounded in [false, true] {
            let (insert, external) = fixture(left, rounded);
            let proof = run(&insert, &external, identity(), [0., 4.], 0.).unwrap();
            assert!(proof.engagement_mm[0] <= 4. && proof.engagement_mm[1] >= 4.);
            assert!(proof.minimum_radial_clearance_lower_mm > 0.005);
            assert_eq!(insert.inner().spec().major_radius_mm, 1.5);
        }
    }
}
fn regular(thread: Thread) -> crate::headed_screw::Screw {
    crate::headed_screw::Screw::new(
        thread,
        crate::headed_screw::ScrewSpec {
            head: crate::headed_screw::ScrewHead::SocketCap {
                radius_mm: 2.75,
                height_mm: 3.,
                socket_across_flats_mm: 2.,
                socket_depth_mm: 1.,
            },
            shank_length_mm: 0.,
            shank_radius_mm: 1.5,
        },
    )
    .unwrap()
}
fn query(insert: &Insert, screw: MatingScrew<'_>, frame: Frame, minimum: f64) -> R<InsertMating> {
    certify_insert_mating(
        insert,
        identity(),
        MatingExternal {
            source: screw,
            frame,
        },
        MatingRequirements {
            minimum_engagement_mm: minimum,
            minimum_clearance_mm: 4e-7,
        },
        &|| false,
    )
}
#[test]
fn phase_interference_nonexact_axes_prefix_penetration_and_short_overlap_refuse() {
    let (insert, external) = fixture(false, true);
    let screw = regular(external);
    for case in 0..6 {
        let mut frame = identity();
        match case {
            0 => {
                frame.axes[0] = [0., 1., 0.];
                frame.axes[1] = [-1., 0., 0.];
            }
            1 => frame.origin_mm[0] = f64::MIN_POSITIVE,
            2 => frame.origin_mm[2] = 0.001,
            3 => {
                let (s, c) = 0.13_f64.sin_cos();
                frame.axes = [[c, s, 0.], [-s, c, 0.], [0., 0., 1.]];
            }
            4 => {
                frame.axes[1] = [0., -1., 0.];
                frame.axes[2] = [0., 0., -1.];
            }
            _ => frame.origin_mm[2] = -4.,
        }
        assert_eq!(
            query(&insert, MatingScrew::Regular(&screw), frame, 2.).unwrap_err(),
            Error::Unsupported,
            "case{case}"
        );
    }
    assert_eq!(
        query(&insert, MatingScrew::Regular(&screw), identity(), 4.1).unwrap_err(),
        Error::Unsupported
    );
    for field in 0..3 {
        let mut s = external.spec();
        match field {
            0 => s.left_hand = true,
            1 => s.starts = 2,
            _ => s.pitch_mm = 0.6,
        };
        let wrong = regular(Thread::new(s).unwrap());
        assert_eq!(
            query(&insert, MatingScrew::Regular(&wrong), identity(), 2.).unwrap_err(),
            Error::Unsupported
        );
    }
    for length in [3., 4.5, 5.5] {
        let mut s = external.spec();
        s.length_mm = length;
        let actual = regular(Thread::new(s).unwrap());
        let proof = query(&insert, MatingScrew::Regular(&actual), identity(), 2.).unwrap();
        assert_eq!(proof.entire_thread_within_insert, length <= 4.);
        assert!(
            proof.engagement_mm[0] <= length.min(4.) && proof.engagement_mm[1] >= length.min(4.)
        );
        if length > 4. {
            assert!(proof.thread_exit_clearance_mm[1] < 0.);
        }
    }
}
#[test]
fn true_quarter_phase_and_common_world_metric_do_not_rewrite_the_sources() {
    for left in [false, true] {
        let (insert, thread) = fixture(left, true);
        let mut inner = insert.inner().spec();
        inner.starts = 4;
        let inner = AuthoredThread::new(
            inner,
            SourceOptions {
                position_limit_mm: 1e-8,
                normal_angle_limit_rad: 1e-6,
            },
        )
        .unwrap();
        let paired = Insert::new(*insert.outer(), inner, 0.5, None).unwrap();
        let mut external = thread.spec();
        external.starts = 4;
        let external = regular(Thread::new(external).unwrap());
        let world = Frame {
            origin_mm: [105., -20., 59.],
            axes: [[0.8, 0.6, 0.], [0.36, -0.48, -0.8], [-0.48, 0.64, -0.6]],
        };
        for quarter in 0..4 {
            let mut placed = world;
            match quarter {
                1 => {
                    placed.axes[0] = world.axes[1];
                    placed.axes[1] = world.axes[0].map(|v| -v);
                }
                2 => {
                    placed.axes[0] = world.axes[0].map(|v| -v);
                    placed.axes[1] = world.axes[1].map(|v| -v);
                }
                3 => {
                    placed.axes[0] = world.axes[1].map(|v| -v);
                    placed.axes[1] = world.axes[0];
                }
                _ => {}
            }
            let proof = certify_insert_mating(
                &paired,
                world,
                MatingExternal {
                    source: MatingScrew::Regular(&external),
                    frame: placed,
                },
                MatingRequirements {
                    minimum_engagement_mm: 3.9,
                    minimum_clearance_mm: 4e-7,
                },
                &|| false,
            )
            .unwrap();
            assert!(proof.engagement_mm[0] >= 3.9);
            assert!(proof.minimum_radial_clearance_lower_mm > 0.005);
        }
    }
}
#[test]
fn actual_profiled_and_phillips_prefixes_use_complete_retained_meridians() {
    use crate::{
        head_profile::{Head, Profile},
        phillips_drive::{self, Cavity, Cutter, PhillipsScrew},
        profiled_screw::{SourceSpec, Spec},
    };
    let (insert, thread) = fixture(false, true);
    let head = Head::Countersunk {
        radius_mm: 2.75,
        height_mm: 1.65,
        included_angle_deg: 90.,
        edge_chamfer_mm: 0.05,
        neck_radius_mm: 0.1,
    };
    let spec = Spec {
        head,
        shank_length_mm: 1.,
        socket_across_flats_mm: 2.,
        socket_depth_mm: 0.6,
        minimum_wall_mm: 1e-5,
    };
    let profiled = crate::profiled_screw::Screw::new_with_shank_radius(thread, 1.5, spec).unwrap();
    assert_eq!(
        query(&insert, MatingScrew::Profiled(&profiled), identity(), 2.).unwrap_err(),
        Error::Unsupported
    );
    let mut above = identity();
    above.origin_mm[2] = -0.5;
    assert!(query(&insert, MatingScrew::Profiled(&profiled), above, 2.).is_ok());
    let (s, d) = phillips_drive::test_fixture(2);
    let cavity = Cavity::new(Cutter::new(s, d).unwrap(), 1e-5, 100_000).unwrap();
    let actual = PhillipsScrew::new(
        thread,
        Profile::new(head, 1.5, 1.).unwrap(),
        SourceSpec {
            head: Some(head),
            shank_length_mm: 1.,
            minimum_wall_mm: 1e-5,
        },
        cavity,
        100_000,
        &|| false,
    )
    .unwrap();
    assert!(query(&insert, MatingScrew::Phillips(&actual), above, 2.).is_ok());
    assert_eq!(
        query(&insert, MatingScrew::Phillips(&actual), identity(), 2.).unwrap_err(),
        Error::Unsupported
    );
}
#[test]
fn cancelled_invalid_and_outward_period_primitives_keep_honest_boundaries() {
    use std::cell::Cell;
    let (insert, thread) = fixture(false, true);
    let screw = regular(thread);
    let calls = Cell::new(0);
    certify_insert_mating(
        &insert,
        identity(),
        MatingExternal {
            source: MatingScrew::Regular(&screw),
            frame: identity(),
        },
        MatingRequirements {
            minimum_engagement_mm: 4.,
            minimum_clearance_mm: 4e-7,
        },
        &|| {
            calls.set(calls.get() + 1);
            false
        },
    )
    .unwrap();
    for target in 1..=calls.get() {
        let n = Cell::new(0);
        assert_eq!(
            certify_insert_mating(
                &insert,
                identity(),
                MatingExternal {
                    source: MatingScrew::Regular(&screw),
                    frame: identity()
                },
                MatingRequirements {
                    minimum_engagement_mm: 4.,
                    minimum_clearance_mm: 4e-7
                },
                &|| {
                    n.set(n.get() + 1);
                    n.get() >= target
                }
            )
            .unwrap_err(),
            Error::Cancelled
        );
    }
    for value in [f64::NAN, 0., -1., f64::INFINITY, 1e8] {
        assert_eq!(
            query(&insert, MatingScrew::Regular(&screw), identity(), value).unwrap_err(),
            Error::Invalid
        );
    }
    assert_eq!(
        exact_interval(&(q(1.).unwrap() / q(3.).unwrap()))
            .unwrap()
            .lo,
        (1. / 3_f64).next_down()
    );
    assert!(
        periodic_radius(I::new(-1., 1.), 0.5, |p| authored_radius(insert.inner(), p))
            .unwrap()
            .lo
            > 1.
    );
    assert_eq!(
        periodic_radius(ip(3000.), 0.5, |p| authored_radius(insert.inner(), p)).unwrap_err(),
        Error::Unsupported
    );
    assert_eq!(
        authored_radius(insert.inner(), I::new(2., 3.)).unwrap_err(),
        Error::Unsupported
    );
}
