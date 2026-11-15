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
fn input(s: &PeriodicThreadSource) -> Input {
    Input {
        target: ArcTarget::new(*s, [3, 4], 1., 0., &mut ledger()).unwrap(),
        chart: chart(),
        vertices: quarter_vertices(),
        same_sense: true,
    }
}
#[test]
fn exact_algebraic_quarter_arc_oracle_is_invariant_to_raw_reference_axis_scale() {
    // Exact quarter-arc endpoints (+sqrt(1/2), +/-sqrt(1/2),0), independent
    // of file metadata, header fitting or sampled trigonometry.
    let proof = certify(&raw(), input(&source(0.25, 0.2)), &mut ledger()).unwrap();
    assert!(proof.distance < 1e-9 && proof.derivative < 1e-9 && proof.incidence < 1e-12);
}
#[test]
fn exact_algebraic_complement_arc_covers_shared_nodes_zero_to_three_at_other_radius() {
    let s = source(0.25, 0.2);
    let target = ArcTarget::new(s, [0, 3], 2., 0., &mut ledger()).unwrap();
    let root = I::point(2.).sqrt();
    let vertices = [[root, -root, I::point(0.)], [root, root, I::point(0.)]];
    let r = RawCircle::new(
        [I::point(0.); 3],
        [I::point(0.), I::point(0.), I::point(-1.)],
        [I::point(0.), I::point(1.), I::point(0.)],
        I::point(2.),
        &mut ledger(),
    )
    .unwrap();
    let proof = certify(
        &r,
        Input {
            target,
            chart: chart(),
            vertices,
            same_sense: true,
        },
        &mut ledger(),
    )
    .unwrap();
    assert!(proof.distance < 1e-9 && proof.derivative < 1e-9);
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
fn opaque_shared_node_arc_target_validates_phase_radius_axial_and_work_domains() {
    let s = source(0.25, 0.2);
    for (nodes, r, z) in [
        ([4, 5], 1., 0.),
        ([3, 3], 1., 0.),
        ([4, 3], 1., 0.),
        ([0, 5], 1., 0.),
        ([3, 4], 0., 0.),
        ([3, 4], f64::NAN, 0.),
        ([3, 4], 1e8, 0.),
        ([3, 4], 1., -1.),
        ([3, 4], 1., 2.),
        ([3, 4], 1., f64::NAN),
    ] {
        assert!(ArcTarget::new(s, nodes, r, z, &mut ledger()).is_err());
    }
    let mut work = ledger();
    work.remaining = 0;
    assert!(ArcTarget::new(s, [3, 4], 1., 0., &mut work).is_err());
    let target = input(&s).target;
    for budget in [0, 256, 4352, 8448] {
        let mut work = ledger();
        work.remaining = budget;
        assert!(target.endpoint(false, &mut work).is_err());
    }
}
#[test]
fn untrusted_numeric_correspondence_frontiers_refuse_without_any_head_construction() {
    let s = source(0.25, 0.2);
    let r = raw();
    let base = input(&s);
    let mut bad = base;
    bad.chart.outer_top_z = f64::NAN;
    assert!(certify(&r, bad, &mut ledger()).is_err());
    bad = base;
    bad.chart.columns[1][0] = f64::INFINITY;
    assert!(certify(&r, bad, &mut ledger()).is_err());
    bad = base;
    bad.vertices[0][0] = I::point(f64::NAN);
    assert!(certify(&r, bad, &mut ledger()).is_err());
    bad = base;
    bad.same_sense = false;
    assert!(certify(&r, bad, &mut ledger()).is_err());
    bad = base;
    bad.vertices[0] = [I::point(0.); 3];
    assert!(certify(&r, bad, &mut ledger()).is_err());
    bad = base;
    bad.vertices[1] = bad.vertices[1].map(|x| -x);
    assert!(certify(&r, bad, &mut ledger()).is_err());
    let ambiguous = RawCircle::new(
        [I::point(0.); 3],
        [I::new(0.1, 0.2), I::point(0.), I::new(-1., 1.)],
        [I::point(0.), I::point(1.), I::point(0.)],
        I::point(1.),
        &mut ledger(),
    )
    .unwrap();
    assert!(certify(&ambiguous, base, &mut ledger()).is_err());
    for budget in [
        0,
        2048,
        2048 + 12544,
        2048 + 2 * 12544,
        2048 + 2 * 12544 + 4096,
    ] {
        let mut work = ledger();
        work.remaining = budget;
        assert!(certify(&r, base, &mut work).is_err());
    }
}
#[test]
fn uncertain_endpoint_angles_do_not_choose_near_zero_or_full_turn_arc_branches() {
    for (crest, depth) in [(1e-14, 0.2), (1. - 1e-14, 2e-15)] {
        let s = source(crest, depth);
        let target = input(&s).target;
        let vertices = [false, true].map(|end| target.endpoint(end, &mut ledger()).unwrap());
        assert!(
            certify(
                &raw(),
                Input {
                    target,
                    vertices,
                    ..input(&s)
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
    let mut base = input(&s);
    base.chart.columns[1][0] = f64::MAX;
    assert!(certify(&raw(), base, &mut ledger()).is_err());
}
