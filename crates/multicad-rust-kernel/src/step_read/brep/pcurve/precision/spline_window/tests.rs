use super::*;
const BASE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../multicad-app/tests/fixtures/exchange/corpus/occt-cone-frustum.step"
));
fn input(y: &str, natural_end: &str) -> String {
    let added = format!(
        "\n#900001=B_SPLINE_CURVE_WITH_KNOTS('',2,(#900010,#900011,#900012),.UNSPECIFIED.,.F.,.F.,(3,3),(0.,{natural_end}),.UNSPECIFIED.);\n#900010=CARTESIAN_POINT('',(0.,0.,0.));\n#900011=CARTESIAN_POINT('',(0.5,0.,0.));\n#900012=CARTESIAN_POINT('',(1.,1.,0.));\n#900020=CARTESIAN_POINT('',(0.25,{y},0.));\n#900021=CARTESIAN_POINT('',(0.75,0.5625,0.));\n#900030=VERTEX_POINT('',#900020);\n#900031=VERTEX_POINT('',#900021);\n#900100=PCURVE('',#31,#900101);\n#900101=DEFINITIONAL_REPRESENTATION('',(#900102),#900103);\n#900103=(GEOMETRIC_REPRESENTATION_CONTEXT(2) PARAMETRIC_REPRESENTATION_CONTEXT() REPRESENTATION_CONTEXT('',''));\n#900102=B_SPLINE_CURVE_WITH_KNOTS('',1,(#900110,#900111),.UNSPECIFIED.,.F.,.F.,(2,2),(0.25,0.75),.UNSPECIFIED.);\n#900110=CARTESIAN_POINT('',(0.,0.));\n#900111=CARTESIAN_POINT('',(1.,1.));\n"
    );
    let at = BASE.rfind("ENDSEC;").unwrap();
    format!("{}{}{}", &BASE[..at], added, &BASE[at..])
}
fn setup<'s, 'i>(
    s: &'s Source<'i>,
) -> (
    super::super::super::uv::Definition,
    Edge,
    CurveCells<'s, 'i, 'static>,
) {
    let mut b = Budget::new(512_000_000, None);
    let d = super::super::super::uv::definition(s, 900100, &mut b).unwrap();
    let curve = s.curve(900001, 0).unwrap();
    let Curve::Nurbs {
        curve: ref nurbs, ..
    } = curve
    else {
        panic!("actual parsed NURBS");
    };
    let packet =
        super::super::bezier::original_curve(s, 900001, nurbs, Coordinates::Raw(3), &mut b)
            .unwrap();
    let edge = Edge {
        curve,
        vertices: [900030, 900031],
        range: [0.25, 0.75],
    };
    (d, edge, packet)
}
#[test]
fn parabola_global_derivatives_original_window_and_exact_resource_frontiers() {
    let text = input("0.0625", "1.");
    let s = Source::parse(&text, 1e-7).unwrap();
    let (d, e, p) = setup(&s);
    let mut b = Budget::new(512_000_000, None);
    let w = prepare(&s, 900001, &d, &e, &p, &mut b).unwrap();
    // Independently C=(t,t²,0), C'=(1,2t,0), C''=(0,2,0).
    assert!(w.gamma() > 0.99 && w.gamma() <= 1.);
    assert!(w.speed() >= 5f64.sqrt() && w.speed() < 2.237);
    assert!(w.acceleration() >= 2. && w.acceleration() < 2.001);
    assert!(w.closest_delta() < 1e-12 && w.closest_delta_sum() < 2e-12);
    assert!(w.cached_window_delta() < 1e-12 && w.numeric_reserve() < 1e-12);
    assert!(w.raw_uv_width().lo <= 0.5 && w.raw_uv_width().hi >= 0.5);
    assert!(w.raw_uv_domain()[0].lo <= 0.25);
    assert!(w.matches(&s, 900001, &d, &e));
    assert!(!w.matches(&s, 900002, &d, &e));
    let other = Source::parse(&text, 1e-7).unwrap();
    assert!(!w.matches(&other, 900001, &d, &e));
    let used = b.used();
    assert!(matches!(
        prepare(&s, 900001, &d, &e, &p, &mut Budget::new(used - 1, None)),
        Err(ReadError::Limit)
    ));
    assert!(prepare(&s, 900001, &d, &e, &p, &mut Budget::new(used, None)).is_ok());
    let calls = std::cell::Cell::new(0);
    let observe = || {
        calls.set(calls.get() + 1);
        false
    };
    prepare(
        &s,
        900001,
        &d,
        &e,
        &p,
        &mut Budget::new(512_000_000, Some(&observe)),
    )
    .unwrap();
    let last = calls.get();
    calls.set(0);
    let cancel = || {
        calls.set(calls.get() + 1);
        calls.get() == last
    };
    let mut b = Budget::new(512_000_000, Some(&cancel));
    assert!(matches!(
        prepare(&s, 900001, &d, &e, &p, &mut b),
        Err(ReadError::Limit)
    ));
    assert!(b.cancelled.get());
    assert!(prepare(&s, 900001, &d, &e, &p, &mut Budget::new(512_000_000, None)).is_ok());
}
#[test]
fn original_off_curve_vertex_and_same_cached_knot_alias_are_not_silent() {
    let text = input("0.0634765625", "1.");
    let s = Source::parse(&text, 1e-7).unwrap();
    let (d, e, p) = setup(&s);
    let w = prepare(&s, 900001, &d, &e, &p, &mut Budget::new(512_000_000, None)).unwrap();
    // At a=.25, original residual is 1/1024, so 2r/gamma=1/512.
    assert!(w.closest_delta() >= 1. / 512. && w.closest_delta() < 0.002);
    // Exact stationarity polynomial 2t³+(1-2y)t-x is strictly increasing.
    // These exact rational signs bracket the unique global minimizer.
    let y = BigRational::new(65.into(), 1024.into());
    let f = |n: i64, den: i64| {
        let t = BigRational::new(n.into(), den.into());
        BigRational::from_integer(2.into()) * &t * &t * &t
            + (BigRational::from_integer(1.into()) - BigRational::from_integer(2.into()) * &y) * &t
            - BigRational::new(1.into(), 4.into())
    };
    assert!(f(1, 4) < BigRational::from_integer(0.into()));
    assert!(f(129, 512) > BigRational::from_integer(0.into()));
    let text = input("0.0625", "1.0000000000000000001");
    let s = Source::parse(&text, 1e-7).unwrap();
    let (d, e, p) = setup(&s);
    assert!(matches!(
        prepare(&s, 900001, &d, &e, &p, &mut Budget::new(512_000_000, None)),
        Err(ReadError::Unsupported)
    ));
}
#[test]
fn installed_dip_original_truncated_targets_earn_finite_windows_separately() {
    let path = "/Applications/KiCad/KiCad.app/Contents/SharedSupport/3dmodels/Package_DIP.3dshapes/DIP-32_W7.62mm.step";
    let Ok(text) = std::fs::read_to_string(path) else {
        eprintln!("unavailable: external KiCad DIP fixture");
        return;
    };
    let s = Source::parse(&text, 1e-7).unwrap();
    let mut b = Budget::new(512_000_000, None);
    let mut checked = 0;
    for binding in super::super::super::bindings(&s, &mut b).unwrap() {
        let target = super::super::unwrapped(&s, binding.curve, &mut b).unwrap();
        if ![28859, 29471].contains(&target) {
            continue;
        }
        let e = super::super::super::edge(&s, binding, &mut b).unwrap();
        let d = super::super::super::uv::definition(&s, binding.pcurve, &mut b).unwrap();
        let Curve::Nurbs { curve: ref c, .. } = e.curve else {
            panic!("actual target spline");
        };
        let Curve::Nurbs { curve: ref u, .. } = d.geometry else {
            panic!("actual UV spline");
        };
        if u.range()[0] < 0.7 {
            continue;
        }
        let packet =
            super::super::bezier::original_curve(&s, target, c, Coordinates::Raw(3), &mut b)
                .unwrap();
        let w = prepare(&s, target, &d, &e, &packet, &mut b).unwrap();
        assert!(w.matches(&s, target, &d, &e));
        assert!(w.gamma() > 0. && w.speed() > 0. && w.acceleration().is_finite());
        assert!(w.raw_uv_width().lo > w.closest_delta_sum());
        eprintln!(
            "actual DIP target{target} UV{} gamma={} L={} H={} closest={} cached={} numeric={} work={}",
            d.curve_id,
            w.gamma(),
            w.speed(),
            w.acceleration(),
            w.closest_delta(),
            w.cached_window_delta(),
            w.numeric_reserve(),
            b.used()
        );
        checked += 1;
    }
    assert_eq!(checked, 2);
}

#[path = "coverage_tests.rs"]
mod coverage_contracts;
