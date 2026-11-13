use super::*;
use crate::{head_profile::Head, helical_thread};
fn thread(left: bool, rounded: bool) -> Thread {
    Thread::new(helical_thread::Spec {
        major_radius_mm: 1.4635,
        pitch_mm: 0.5,
        length_mm: 8.,
        depth_mm: 0.29908756729740676,
        crest_width_mm: 0.07131009168581062,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: left,
        root_radius_mm: rounded.then_some(0.07216878364870323),
    })
    .unwrap()
}
fn head(crowned: bool, plain: f64) -> Profile {
    let h = if crowned {
        Head::Crowned {
            radius_mm: 2.8,
            height_mm: 2.4,
            crown_radius_mm: 3.,
            top_land_radius_mm: 2.1,
            neck_radius_mm: 0.1,
        }
    } else {
        Head::Countersunk {
            radius_mm: 2.75,
            height_mm: 1.65,
            included_angle_deg: 90.,
            edge_chamfer_mm: 0.05,
            neck_radius_mm: 0.1,
        }
    };
    Profile::new_for_thread_union(h, 1.5, plain).unwrap()
}
#[test]
fn actual_h1_counter_and_pan_fit_retained_heads_and_unchanged_metric_core() {
    for crowned in [false, true] {
        let c = super::super::tests::cavity(if crowned { 1 } else { 2 });
        for plain in [0., 0.2, 1.] {
            let p = head(crowned, plain);
            for left in [false, true] {
                for rounded in [false, true] {
                    let t = thread(left, rounded);
                    let original = t;
                    let proof = c
                        .certify_host_envelope(&p, &t, plain, 1e-5, 100_000, &|| false)
                        .unwrap();
                    assert!(proof.cells() > 0 && proof.minimum_radial_clearance_mm() > 1e-5);
                    assert!(proof.axial_floor_clearance_mm() > 1e-5);
                    assert_eq!(t, original);
                    // Actual independent radial probes along every section's
                    // untouched diagonal remain inside the host envelope.
                    for i in 0..100 {
                        let d = c.source.envelope.tip.hi * i as f64 / 100.;
                        let z = p.top()[1] + d;
                        let radius = super::super::tests::oracle(c.source.spec(), d);
                        let host = if z < plain {
                            p.radius_at(z).unwrap()
                        } else {
                            t.spec().major_radius_mm - t.spec().depth_mm
                        };
                        assert!(host > radius + 1e-5, "{crowned}/{plain}/{d}");
                    }
                }
            }
        }
    }
    // A physically narrower common core requires cells wholly below the
    // original thread datum; a mouth or single cross-datum sample cannot fit.
    let c = super::super::tests::cavity(2);
    let mut spec = thread(false, true).spec();
    spec.major_radius_mm = 1.1;
    let t = Thread::new(spec).unwrap();
    assert!(
        c.certify_host_envelope(&head(false, 0.), &t, 0., 0.001, 100_000, &|| false)
            .is_ok()
    );
}
#[test]
fn all_six_actual_metric_heads_fit_their_authored_cavity_envelopes() {
    for (index, (rs, pitch, radius, height, crown, land, neck)) in [
        (1., 0.4, 2., 1.6, 3.2, 1.25, 0.1),
        (1.5, 0.5, 2.8, 2.4, 3., 2.1, 0.1),
        (1.5, 0.5, 2.75, 1.65, 0., 0., 0.1),
        (2., 0.7, 4., 3.1, 6.5, 3., 0.2),
        (3., 1., 6., 4.6, 10., 4.2, 0.25),
        (4., 1.25, 8., 6., 13., 5.4, 0.4),
    ]
    .into_iter()
    .enumerate()
    {
        let source = super::super::tests::cavity(index);
        let h = if index == 2 {
            Head::Countersunk {
                radius_mm: radius,
                height_mm: height,
                included_angle_deg: 90.,
                edge_chamfer_mm: 0.05,
                neck_radius_mm: neck,
            }
        } else {
            Head::Crowned {
                radius_mm: radius,
                height_mm: height,
                crown_radius_mm: crown,
                top_land_radius_mm: land,
                neck_radius_mm: neck,
            }
        };
        let p = Profile::new_for_thread_union(h, rs, 0.).unwrap();
        let mut spec = thread(false, false).spec();
        spec.major_radius_mm = rs - 0.03;
        spec.pitch_mm = pitch;
        spec.depth_mm = 0.6 * pitch;
        spec.crest_width_mm = 0.125 * pitch;
        let t = Thread::new(spec).unwrap();
        source
            .certify_host_envelope(&p, &t, 0., 1e-5, 100_000, &|| false)
            .unwrap();
    }
}
#[test]
fn retained_support_bound_covers_crown_neck_bearing_and_outside_domain() {
    for crowned in [false, true] {
        let p = head(crowned, 1.);
        for s in p.segments() {
            let lo = head_radius_lower(&p, I::new(s.start[1], s.end[1]));
            for i in 0..=30 {
                assert!(lo <= s.point(i as f64 / 30.)[0]);
            }
        }
        assert_eq!(head_radius_lower(&p, I::new(10., 11.)), 0.);
    }
}
#[test]
fn finite_work_cancel_wrong_datums_and_real_wall_breaches_refuse() {
    let c = super::super::tests::cavity(2);
    let p = head(false, 0.);
    let t = thread(false, true);
    let run = |profile: &Profile, thread: &Thread, plain, wall, work, cancel: &dyn Fn() -> bool| {
        c.certify_host_envelope(profile, thread, plain, wall, work, cancel)
    };
    for (plain, wall) in [(f64::NAN, 0.1), (-1., 0.1), (0., 0.), (2., 0.1)] {
        assert_eq!(
            run(&p, &t, plain, wall, 1000, &|| false),
            Err(Error::Invalid)
        );
    }
    let thin_shank = Profile::new_for_thread_union(
        Head::Cylindrical {
            radius_mm: 3.,
            height_mm: 2.4,
        },
        1.,
        0.,
    )
    .unwrap();
    assert_eq!(
        run(&thin_shank, &t, 0., 0.01, 1000, &|| false),
        Err(Error::Invalid)
    );
    for work in [0, 1_000_001, 1] {
        assert_eq!(run(&p, &t, 0., 0.01, work, &|| false), Err(Error::Limit));
    }
    assert_eq!(run(&p, &t, 0., 0.01, 1000, &|| true), Err(Error::Cancelled));
    let checks = std::cell::Cell::new(0);
    assert_eq!(
        run(&p, &t, 0., 0.01, 1000, &|| {
            checks.set(checks.get() + 1);
            checks.get() > 2
        }),
        Err(Error::Cancelled)
    );
    assert_eq!(run(&p, &t, 0., 1., 1000, &|| false), Err(Error::Unresolved));
    let mut short = t.spec();
    short.length_mm = 0.5;
    assert_eq!(
        run(&p, &Thread::new(short).unwrap(), 0., 0.01, 1000, &|| false),
        Err(Error::Unresolved)
    );
    // The mouth fits this shallow cylinder, but the true cut later exits its
    // thin original-thread core. Subdivision cannot turn a breach into proof.
    let shallow = Profile::new_for_thread_union(
        Head::Cylindrical {
            radius_mm: 2.75,
            height_mm: 0.3,
        },
        1.5,
        0.,
    )
    .unwrap();
    assert!(matches!(
        run(&shallow, &t, 0., 0.01, 10_000, &|| false),
        Err(Error::Limit | Error::Unresolved)
    ));
}
