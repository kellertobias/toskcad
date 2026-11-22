use super::*;
use crate::{Circle, chamfered_bores::Counterbore};
use std::cell::Cell;
fn frame(z: f64) -> Frame {
    Frame {
        origin_mm: [0., 0., z],
        axes: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
    }
}
fn prism(height: Rational, half: i64) -> PolygonPrism {
    let n = |v| Rational::new(v, 1).unwrap();
    PolygonPrism::with_circular_holes(
        vec![
            vec![[-half, -half], [half, -half], [half, half], [-half, half]]
                .into_iter()
                .map(|p| p.map(n))
                .collect(),
        ],
        vec![Circle::new([n(0), n(0)], Rational::new(17, 10).unwrap()).unwrap()],
        height,
    )
    .unwrap()
}
fn cut(p: &PolygonPrism, end: bool) -> ChamferedBores {
    ChamferedBores::counterbored(
        p,
        &[Counterbore {
            hole: 0,
            end,
            radius_mm: 3.5,
            depth_mm: 1.,
        }],
        1e-7,
        1e-5,
    )
    .unwrap()
}
fn requirements() -> InsertClampRequirements {
    InsertClampRequirements {
        minimum_wall_mm: 1.,
        minimum_bearing_land_mm: 1e-5,
        maximum_bearing_gap_mm: 1e-7,
        linear_tolerance_mm: 1e-7,
        minimum_feature_mm: 1e-5,
    }
}
#[allow(clippy::too_many_arguments)]
fn run(
    p: &PolygonPrism,
    c: &ChamferedBores,
    s: &Screw,
    i: &Insert,
    sf: Frame,
    inf: Frame,
    r: InsertClampRequirements,
    cancel: &dyn Fn() -> bool,
) -> R<InsertClampContact> {
    certify_insert_clamp(
        InsertClampSources {
            prism: p,
            cut: c,
            clamp_frame: frame(0.),
            screw: s,
            screw_frame: sf,
            insert: i,
            insert_frame: inf,
        },
        r,
        cancel,
    )
}
#[test]
fn actual_rational_end_plane_and_whole_bearing_are_preserved_for_both_hands() {
    for left in [false, true] {
        let (i, s, mut sf) = super::super::clocked::tests::fixture(left);
        let p = prism(Rational::new(41, 10).unwrap(), 10);
        let c = cut(&p, false);
        sf.origin_mm = [0., 0., 1. - 1e-10];
        let inf = frame(4.1 + 1e-10);
        let proof = run(&p, &c, &s, &i, sf, inf, requirements(), &|| false).unwrap();
        assert!(proof.bearing_face_gap_mm[0] > 0. && proof.bearing_face_gap_mm[1] < 1e-7);
        assert!(proof.insert_to_clamp_gap_mm[0] > 0.);
        assert!(proof.bearing_land_lower_mm > 1.);
        assert!(proof.bore_mouths_mm[1][2][0] <= 4.1 && proof.bore_mouths_mm[1][2][1] >= 4.1);
        assert!(matches!(
            run(&p, &c, &s, &i, sf, frame(4.1), requirements(), &|| false),
            Err(Error::Unsupported)
        ));
        let exact = prism(Rational::new(4, 1).unwrap(), 10);
        let ec = cut(&exact, false);
        sf.origin_mm[2] = 1.;
        let contact = run(&exact, &ec, &s, &i, sf, frame(4.), requirements(), &|| {
            false
        })
        .unwrap();
        assert_eq!(contact.bearing_face_gap_mm, [0.; 2]);
        assert_eq!(contact.insert_to_clamp_gap_mm, [0.; 2]);
    }
}
#[test]
fn actual_intrusion_wrong_sources_land_walls_and_frames_refuse() {
    let (i, s, mut sf) = super::super::clocked::tests::fixture(false);
    sf.origin_mm = [0., 0., 1. - 1e-10];
    let p = prism(Rational::new(41, 10).unwrap(), 10);
    let c = cut(&p, false);
    let inf = frame(4.1 + 1e-10);
    let mut r = requirements();
    r.minimum_wall_mm = f64::NAN;
    assert!(matches!(
        run(&p, &c, &s, &i, sf, inf, r, &|| false),
        Err(Error::Invalid)
    ));
    for wall in [7., 20.] {
        let mut r = requirements();
        r.minimum_wall_mm = wall;
        assert!(run(&p, &c, &s, &i, sf, inf, r, &|| false).is_err());
    }
    let mut r = requirements();
    r.minimum_bearing_land_mm = 3.;
    assert!(run(&p, &c, &s, &i, sf, inf, r, &|| false).is_err());
    let other = prism(Rational::new(5, 1).unwrap(), 10);
    assert!(run(&other, &c, &s, &i, sf, inf, requirements(), &|| false).is_err());
    let reverse = cut(&p, true);
    assert!(run(&p, &reverse, &s, &i, sf, inf, requirements(), &|| false).is_err());
    for z in [1. + 1e-10, 0.9] {
        let mut bad = sf;
        bad.origin_mm[2] = z;
        assert!(run(&p, &c, &s, &i, bad, inf, requirements(), &|| false).is_err());
    }
    let mut bad = sf;
    bad.axes[2][2] = -1.;
    assert!(run(&p, &c, &s, &i, bad, inf, requirements(), &|| false).is_err());
    let mut bad = sf;
    bad.origin_mm[0] = 1.;
    assert!(run(&p, &c, &s, &i, bad, inf, requirements(), &|| false).is_err());
    let mut bad = inf;
    bad.axes[2][2] = -1.;
    assert!(run(&p, &c, &s, &i, sf, bad, requirements(), &|| false).is_err());
}
#[test]
fn every_actual_success_checkpoint_can_cancel_without_publication() {
    let (i, s, mut sf) = super::super::clocked::tests::fixture(false);
    sf.origin_mm = [0., 0., 1. - 1e-10];
    let inf = frame(4.1 + 1e-10);
    let p = prism(Rational::new(41, 10).unwrap(), 10);
    let c = cut(&p, false);
    let calls = Cell::new(0);
    run(&p, &c, &s, &i, sf, inf, requirements(), &|| {
        calls.set(calls.get() + 1);
        false
    })
    .unwrap();
    for stop in 1..=calls.get() {
        let n = Cell::new(0);
        assert!(matches!(
            run(&p, &c, &s, &i, sf, inf, requirements(), &|| {
                n.set(n.get() + 1);
                n.get() == stop
            }),
            Err(Error::Cancelled)
        ));
    }
}

#[test]
fn real_alternative_stock_topologies_and_finite_size_limits_remain_explicit() {
    let (i, s, mut sf) = super::super::clocked::tests::fixture(false);
    sf.origin_mm = [0., 0., 1. - 1e-10];
    let inf = frame(4.1 + 1e-10);
    let n = |v| Rational::new(v, 1).unwrap();
    let circle = |x| Circle::new([n(x), n(0)], Rational::new(17, 10).unwrap()).unwrap();
    let make = |rings: Vec<Vec<[i64; 2]>>, circles| {
        PolygonPrism::with_circular_holes(
            rings
                .into_iter()
                .map(|ring| ring.into_iter().map(|p| p.map(n)).collect())
                .collect(),
            circles,
            Rational::new(41, 10).unwrap(),
        )
        .unwrap()
    };
    for p in [
        make(vec![vec![[-10, -10], [10, -10], [0, 20]]], vec![circle(0)]),
        make(
            vec![vec![[-10, 0], [0, -10], [10, 0], [0, 10]]],
            vec![circle(0)],
        ),
        make(
            vec![
                vec![[-10, -10], [10, -10], [10, 10], [-10, 10]],
                vec![[-8, -8], [-8, -6], [-6, -6], [-6, -8]],
            ],
            vec![circle(0)],
        ),
        make(
            vec![vec![[-20, -10], [20, -10], [20, 10], [-20, 10]]],
            vec![circle(0), circle(10)],
        ),
    ] {
        let c = cut(&p, false);
        assert!(matches!(
            run(&p, &c, &s, &i, sf, inf, requirements(), &|| false),
            Err(Error::Unsupported)
        ));
    }
    let p = prism(Rational::new(41, 10).unwrap(), 10);
    let c = ChamferedBores::counterbored(
        &p,
        &[
            Counterbore {
                hole: 0,
                end: false,
                radius_mm: 3.5,
                depth_mm: 1.,
            },
            Counterbore {
                hole: 0,
                end: true,
                radius_mm: 3.5,
                depth_mm: 1.,
            },
        ],
        1e-7,
        1e-5,
    )
    .unwrap();
    assert!(run(&p, &c, &s, &i, sf, inf, requirements(), &|| false).is_err());
    for p in [
        make(
            vec![vec![
                [9_999_991, -10],
                [10_000_011, -10],
                [10_000_011, 10],
                [9_999_991, 10],
            ]],
            vec![circle(10_000_001)],
        ),
        prism(n(10_000_001), 10),
    ] {
        let c = cut(&p, false);
        assert!(matches!(
            run(&p, &c, &s, &i, sf, inf, requirements(), &|| false),
            Err(Error::Limit)
        ));
    }
}
#[test]
fn actual_plain_shaft_and_proper_reversed_axes_are_not_inferred_from_head_fit() {
    let (i, s, mut sf) = super::super::clocked::tests::fixture(false);
    let p = prism(Rational::new(41, 10).unwrap(), 10);
    let c = cut(&p, false);
    let inf = frame(4.1 + 1e-10);
    for radius in [1.5, 1.8] {
        let mut spec = s.spec();
        spec.shank_length_mm = 0.5;
        spec.shank_radius_mm = radius;
        let actual = Screw::new(*s.thread(), spec).unwrap();
        sf.origin_mm = [0., 0., 1.5 - 1e-10];
        assert_eq!(
            run(&p, &c, &actual, &i, sf, inf, requirements(), &|| false).is_ok(),
            radius == 1.5
        );
    }
    sf = frame(1. - 1e-10);
    sf.axes[0][0] = -1.;
    sf.axes[2][2] = -1.;
    assert!(run(&p, &c, &s, &i, sf, inf, requirements(), &|| false).is_err());
    let mut reverse = inf;
    reverse.axes[0][0] = -1.;
    reverse.axes[2][2] = -1.;
    assert!(
        run(
            &p,
            &c,
            &s,
            &i,
            frame(1. - 1e-10),
            reverse,
            requirements(),
            &|| false
        )
        .is_err()
    );
}

#[test]
fn actual_near_identity_extreme_dyadic_charts_keep_numeric_and_cancel_guards() {
    let (i, s, _) = super::super::clocked::tests::fixture(false);
    let p = prism(Rational::new(41, 10).unwrap(), 10);
    let c = cut(&p, false);
    for exponent in [600, 1000] {
        let e = 2_f64.powi(-exponent);
        let chart = Frame {
            origin_mm: [0.; 3],
            axes: [[1., e, 0.], [0., 1., e], [e, 0., 1.]],
        };
        let mut screw_frame = chart;
        screw_frame.origin_mm = [0., 0., 1. - 1e-10];
        let mut insert_frame = chart;
        insert_frame.origin_mm = [0., 0., 4.1 + 1e-10];
        let query = |interrupt: &dyn Fn() -> bool| {
            certify_insert_clamp(
                InsertClampSources {
                    prism: &p,
                    cut: &c,
                    clamp_frame: chart,
                    screw: &s,
                    screw_frame,
                    insert: &i,
                    insert_frame,
                },
                requirements(),
                interrupt,
            )
        };
        let calls = Cell::new(0);
        let proof = query(&|| {
            calls.set(calls.get() + 1);
            false
        })
        .unwrap();
        assert!(proof.bearing_face_gap_mm[0] > 0. && proof.insert_to_clamp_gap_mm[0] > 0.);
        for stop in 1..=calls.get() {
            let n = Cell::new(0);
            assert!(matches!(
                query(&|| {
                    n.set(n.get() + 1);
                    n.get() == stop
                }),
                Err(Error::Cancelled)
            ));
        }
    }
}
