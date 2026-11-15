use super::*;
fn ledger() -> Work<'static> {
    Work {
        remaining: WORK,
        interrupt: None,
        cancelled: Cell::new(false),
    }
}
fn source(crest: f64, depth: f64) -> PeriodicThreadSource {
    use crate::helical_thread::{Spec, Thread};
    let thread = Thread::new(Spec {
        major_radius_mm: 1.,
        pitch_mm: 1.,
        length_mm: 1.,
        depth_mm: depth,
        crest_width_mm: crest,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: false,
        root_radius_mm: None,
    })
    .unwrap();
    PeriodicThreadSource::new(thread, &mut ledger()).unwrap()
}
fn chart() -> CommonMouthChart {
    CommonMouthChart {
        columns: [[0.; 3], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        outer_top_z: 0.,
    }
}
fn raw() -> RawCircle {
    RawCircle::new(
        [I::point(0.); 3],
        [I::point(0.), I::point(0.), I::point(-20.)],
        [I::point(0.), I::point(5.), I::point(0.)],
        I::point(1.),
        &mut ledger(),
    )
    .unwrap()
}
fn quarter_vertices() -> [[I; 3]; 2] {
    let s = I::point(0.5).sqrt();
    [[s, s, I::point(0.)], [s, -s, I::point(0.)]]
}
#[test]
fn exact_algebraic_quarter_arc_oracle_is_invariant_to_raw_reference_axis_scale() {
    // The crest covers pi/4 to-pi/4: true quarter-arc endpoints are exactly
    // (+sqrt(1/2), +/-sqrt(1/2),0). No file/header or sampled trig supplies them.
    let proof = certify(
        &raw(),
        Input {
            source: &source(0.25, 0.2),
            chart: chart(),
            band: 3,
            z: 0.,
            vertices: quarter_vertices(),
            same_sense: true,
        },
        &mut ledger(),
    )
    .unwrap();
    assert!(proof.distance < 1e-9 && proof.derivative < 1e-9 && proof.incidence < 1e-12);
}
#[test]
fn opaque_circle_constructor_rejects_invalid_domains_and_exhausted_work() {
    let o = [I::point(0.); 3];
    let n = [I::point(0.), I::point(0.), I::point(1.)];
    let u = [I::point(1.), I::point(0.), I::point(0.)];
    for r in [
        I::point(0.),
        I::point(-1.),
        I::new(2., 1.),
        I::point(f64::NAN),
    ] {
        assert!(RawCircle::new(o, n, u, r, &mut ledger()).is_err());
    }
    let mut bad = o;
    bad[0] = I::point(f64::INFINITY);
    assert!(RawCircle::new(bad, n, u, I::point(1.), &mut ledger()).is_err());
    assert!(RawCircle::new(o, o, u, I::point(1.), &mut ledger()).is_err());
    assert!(RawCircle::new(o, n, n, I::point(1.), &mut ledger()).is_err());
    let mut work = ledger();
    work.remaining = 0;
    assert!(RawCircle::new(o, n, u, I::point(1.), &mut work).is_err());
}
#[test]
fn untrusted_numeric_correspondence_frontiers_refuse_without_any_head_construction() {
    let s = source(0.25, 0.2);
    let r = raw();
    let c = chart();
    let v = quarter_vertices();
    for band in [0, 4] {
        assert!(
            certify(
                &r,
                Input {
                    source: &s,
                    chart: c,
                    band,
                    z: 0.,
                    vertices: v,
                    same_sense: true
                },
                &mut ledger()
            )
            .is_err()
        );
    }
    for z in [f64::NAN, -1., 2.] {
        assert!(
            certify(
                &r,
                Input {
                    source: &s,
                    chart: c,
                    band: 3,
                    z,
                    vertices: v,
                    same_sense: true
                },
                &mut ledger()
            )
            .is_err()
        );
    }
    let mut bad = c;
    bad.outer_top_z = f64::NAN;
    assert!(
        certify(
            &r,
            Input {
                source: &s,
                chart: bad,
                band: 3,
                z: 0.,
                vertices: v,
                same_sense: true
            },
            &mut ledger()
        )
        .is_err()
    );
    bad = c;
    bad.columns[1][0] = f64::INFINITY;
    assert!(
        certify(
            &r,
            Input {
                source: &s,
                chart: bad,
                band: 3,
                z: 0.,
                vertices: v,
                same_sense: true
            },
            &mut ledger()
        )
        .is_err()
    );
    let mut badv = v;
    badv[0][0] = I::point(f64::NAN);
    assert!(
        certify(
            &r,
            Input {
                source: &s,
                chart: c,
                band: 3,
                z: 0.,
                vertices: badv,
                same_sense: true
            },
            &mut ledger()
        )
        .is_err()
    );
    assert!(
        certify(
            &r,
            Input {
                source: &s,
                chart: c,
                band: 3,
                z: 0.,
                vertices: v,
                same_sense: false
            },
            &mut ledger()
        )
        .is_err()
    );
    badv = v;
    badv[0] = [I::point(0.); 3];
    assert!(
        certify(
            &r,
            Input {
                source: &s,
                chart: c,
                band: 3,
                z: 0.,
                vertices: badv,
                same_sense: true
            },
            &mut ledger()
        )
        .is_err()
    );
    badv = v;
    badv[1] = v[1].map(|x| -x);
    assert!(
        certify(
            &r,
            Input {
                source: &s,
                chart: c,
                band: 3,
                z: 0.,
                vertices: badv,
                same_sense: true
            },
            &mut ledger()
        )
        .is_err()
    );
    let ambiguous = RawCircle::new(
        [I::point(0.); 3],
        [I::new(0.1, 0.2), I::point(0.), I::new(-1., 1.)],
        [I::point(0.), I::point(1.), I::point(0.)],
        I::point(1.),
        &mut ledger(),
    )
    .unwrap();
    assert!(
        certify(
            &ambiguous,
            Input {
                source: &s,
                chart: c,
                band: 3,
                z: 0.,
                vertices: v,
                same_sense: true
            },
            &mut ledger()
        )
        .is_err()
    );
    for budget in [
        0,
        2048,
        2048 + 12544,
        2048 + 2 * 12544,
        2048 + 2 * 12544 + 4096,
    ] {
        let mut work = ledger();
        work.remaining = budget;
        assert!(
            certify(
                &r,
                Input {
                    source: &s,
                    chart: c,
                    band: 3,
                    z: 0.,
                    vertices: v,
                    same_sense: true
                },
                &mut work
            )
            .is_err()
        );
    }
}
#[test]
fn uncertain_endpoint_angles_do_not_choose_near_zero_or_full_turn_arc_branches() {
    for (crest, depth) in [(1e-14, 0.2), (1. - 1e-14, 2e-15)] {
        let s = source(crest, depth);
        let vertices = [0., 1.].map(|u| {
            s.bounds(3, [u; 2], [0.; 2], WORK, None)
                .unwrap()
                .position_mm
                .map(|x| I::new(x[0], x[1]))
        });
        assert!(
            certify(
                &raw(),
                Input {
                    source: &s,
                    chart: chart(),
                    band: 3,
                    z: 0.,
                    vertices,
                    same_sense: true
                },
                &mut ledger()
            )
            .is_err()
        );
    }
}

#[test]
fn finite_input_with_overflowing_world_enclosure_does_not_publish_infinite_receipt() {
    let s = source(0.25, 0.2);
    let mut c = chart();
    c.columns[1][0] = f64::MAX;
    assert!(
        certify(
            &raw(),
            Input {
                source: &s,
                chart: c,
                band: 3,
                z: 0.,
                vertices: quarter_vertices(),
                same_sense: true
            },
            &mut ledger()
        )
        .is_err()
    );
}
