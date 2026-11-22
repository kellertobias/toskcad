//! A private argument-consistency contract, not a public-path reach claim.
use super::*;
#[test]
fn real_different_curve_record_cannot_supply_this_boundary_degree() {
    let text = input(false);
    let at = text.rfind("ENDSEC;").unwrap();
    let text = format!(
        "{}\n#990001=B_SPLINE_CURVE_WITH_KNOTS('',1,(#900010,#900011,#900012),.UNSPECIFIED.,.F.,.F.,(2,1,2),(0.,0.5,1.),.UNSPECIFIED.);\n{}",
        &text[..at],
        &text[at..]
    );
    let source = Source::parse(&text, 1e-7).unwrap();
    let mut b = Budget::new(512_000_000, None);
    let (d, s, e) = pair(&source, &mut b);
    let Surface::Nurbs { surface, .. } = s else {
        panic!("actual parsed tensor")
    };
    let Curve::Nurbs { curve: own, .. } = e.curve else {
        panic!("actual parsed target")
    };
    let Curve::Nurbs { curve: other, .. } = source.curve(990001, 0).unwrap() else {
        panic!("actual independent curve")
    };
    assert_eq!(own.points().len(), other.points().len());
    assert_ne!(own.degree(), other.degree());
    let iso = original_uv(&source, &d, &mut b).unwrap();
    assert!(matches!(
        coefficient_errors(&source, 900040, &surface, 900001, &other, &iso, &mut b),
        Err(ReadError::Unsupported)
    ));
    let (position, derivative) =
        coefficient_errors(&source, 900040, &surface, 900001, &own, &iso, &mut b).unwrap();
    assert!(position < 1e-12 && derivative < 1e-12);
}
