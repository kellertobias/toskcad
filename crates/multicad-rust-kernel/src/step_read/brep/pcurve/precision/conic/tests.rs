use super::*;
use crate::step_read::brep::pcurve::{Proof, uv};
const TOL: f64 = 1e-7;
const INPUT: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/src/step_read/tests/fixtures/occt-annulus.step"
));
fn fixture(reverse: bool, disturb: bool) -> String {
    let r = 6.125;
    let h = 3_f64.sqrt() * r;
    let mut points = [
        [r, 0.],
        [r, h],
        [-r / 2., h / 2.],
        [-2. * r, 0.],
        [-r / 2., -h / 2.],
        [r, -h],
    ];
    if reverse {
        for p in &mut points {
            p[1] = -p[1];
        }
    }
    if disturb {
        points[1][0] += 0.01;
    }
    let mut records = String::new();
    for (i, p) in points.iter().enumerate() {
        records.push_str(&format!(
            "#{} = CARTESIAN_POINT('',({:.17},{:.17}));\n",
            900000 + i,
            p[0],
            p[1]
        ));
    }
    // Outer multiplicity one is intentional: the active natural domain is
    // [0,3], with three genuine quadratic Bezier spans and shared end poles.
    let curve = "#72 = ( BOUNDED_CURVE() B_SPLINE_CURVE(2,(#900000,#900001,#900002,#900003,#900004,#900005,#900000),.UNSPECIFIED.,.F.,.F.) B_SPLINE_CURVE_WITH_KNOTS((1,2,2,2,2,1),(-1.,0.,1.,2.,3.,4.),.UNSPECIFIED.) CURVE() GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_CURVE((1.,0.5,1.,0.5,1.,0.5,1.)) REPRESENTATION_ITEM('') );";
    let text = INPUT.replace("#72 = CIRCLE('',#73,6.125);", curve);
    assert_ne!(text, INPUT);
    text.replacen(
        "ENDSEC;\nEND-ISO-10303-21;",
        &format!("{records}ENDSEC;\nEND-ISO-10303-21;"),
        1,
    )
}
fn pair(source: &Source<'_>) -> (uv::Definition, Surface, Edge) {
    let d = uv::definition(source, 65, &mut Budget::new(512_000_000, None)).unwrap();
    let curve = source.curve(53, 0).unwrap();
    let p = source.point(25).unwrap();
    let edge = Edge::new(curve, [0, 0], [p, p], TOL).unwrap();
    let surface = source.surface(d.surface).unwrap();
    (d, surface, edge)
}
fn run(source: &Source<'_>, budget: &mut Budget<'_>) -> Result<Option<Proof>> {
    let (d, s, e) = pair(source);
    certify(source, 53, &d, &s, &e, TOL, budget)
}
#[test]
fn actual_unclamped_rational_circle_earns_whole_one_winding_not_samples() {
    let text = fixture(false, false);
    let source = Source::parse(&text, TOL).unwrap();
    let mut budget = Budget::new(512_000_000, None);
    assert!(matches!(run(&source,&mut budget).unwrap(),Some(Proof::Certified(e))if e<=TOL));
    let used = budget.used();
    assert!(used > 0);
    assert!(matches!(
        run(&source, &mut Budget::new(used - 1, None)),
        Err(ReadError::Limit)
    ));
    assert!(matches!(
        run(&source, &mut Budget::new(used, None)).unwrap(),
        Some(Proof::Certified(_))
    ));
    let calls = std::cell::Cell::new(0);
    let cancel = || {
        calls.set(calls.get() + 1);
        calls.get() > 12
    };
    assert!(matches!(
        run(&source, &mut Budget::new(512_000_000, Some(&cancel))),
        Err(ReadError::Limit)
    ));
    assert!(calls.get() > 12);
    assert!(matches!(
        run(&source, &mut Budget::new(512_000_000, None)).unwrap(),
        Some(Proof::Certified(_))
    ));
}
#[test]
fn actual_same_endpoint_wrong_sense_and_changed_control_do_not_gain_certificates() {
    for (reverse, disturb) in [(true, false), (false, true)] {
        let text = fixture(reverse, disturb);
        let source = Source::parse(&text, TOL).unwrap();
        assert!(matches!(
            run(&source, &mut Budget::new(512_000_000, None)),
            Ok(Some(Proof::Unsupported)) | Err(ReadError::Unsupported)
        ));
    }
    let text = fixture(false, false);
    let source = Source::parse(&text, TOL).unwrap();
    assert!(matches!(
        run(&source, &mut Budget::new(512_000_000, None)).unwrap(),
        Some(Proof::Certified(_))
    ));
}
