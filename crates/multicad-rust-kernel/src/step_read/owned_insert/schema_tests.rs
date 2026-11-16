use super::{
    tests::{argument_change, file, record_change},
    *,
};

fn entity(p: &Parsed<'_>, name: &str) -> u32 {
    *p.entities
        .iter()
        .find(|(_, e)| matches!(e.as_slice(), [e] if e.name == name))
        .unwrap()
        .0
}
fn refuses(text: &str, expected: Error) {
    assert_eq!(parse(text, Options::default(), None).err(), Some(expected));
}
#[test]
fn actual_file_rejects_malformed_carrier_and_placement_records() {
    let text = file(false);
    let p = parse(&text, Options::default(), None).unwrap();
    let curve = entity(&p, "B_SPLINE_CURVE_WITH_KNOTS");
    let tensor = entity(&p, "B_SPLINE_SURFACE_WITH_KNOTS");
    let plane = p.planes()[0].entity();
    let axis = p.planes()[0].axis_entity();
    let point = p.points()[0].entity();
    for (id, index, value, expected) in [
        (curve, 4, ".T.", Error::Unsupported),
        (tensor, 5, ".T.", Error::Unsupported),
        (tensor, 3, "((#1),(#1,#1))", Error::Invalid),
        (tensor, 8, "()", Error::Invalid),
        (tensor, 10, "()", Error::Limit),
        (plane, 0, "$", Error::Invalid),
        (axis, 0, "$", Error::Unsupported),
        (point, 1, "(0.)", Error::Unsupported),
    ] {
        refuses(&argument_change(&text, id, index, value), expected);
    }
    let oversized = format!("({})", vec!["0."; 1026].join(","));
    refuses(&argument_change(&text, curve, 7, &oversized), Error::Limit);
    refuses(
        &argument_change(&text, tensor, 11, &oversized),
        Error::Limit,
    );
}
#[test]
fn actual_file_rejects_malformed_face_wire_edge_and_pcurve_records() {
    let text = file(false);
    let p = parse(&text, Options::default(), None).unwrap();
    let face = p.faces()[0].entity();
    let bound = p.faces()[0].loops()[0].bound_entity();
    let wire = p.faces()[0].loops()[0].entity();
    let coedge = entity(&p, "ORIENTED_EDGE");
    let edge = &p.edges()[0];
    let pc = edge.pcurves()[0];
    for (id, index, value, expected) in [
        (face, 0, "$", Error::Invalid),
        (
            face,
            2,
            &format!("#{}", p.points()[0].entity()),
            Error::Unsupported,
        ),
        (bound, 2, ".F.", Error::Unsupported),
        (coedge, 1, "$", Error::Unsupported),
        (wire, 1, "()", Error::Invalid),
        (edge.entity(), 4, ".F.", Error::Unsupported),
        (edge.surface_curve_entity(), 2, "()", Error::Unsupported),
    ] {
        refuses(&argument_change(&text, id, index, value), expected);
    }
    for (id, replacement, expected) in [
        (
            bound,
            format!("#{bound}=UNSUPPORTED_BOUND('',#{wire},.T.);"),
            Error::Unsupported,
        ),
        (wire, format!("#{wire}=EDGE_LOOP('');"), Error::Invalid),
        (
            pc.entity,
            format!("#{}=PCURVE('');", pc.entity),
            Error::Invalid,
        ),
        (
            pc.representation,
            format!("#{}=DEFINITIONAL_REPRESENTATION('');", pc.representation),
            Error::Invalid,
        ),
    ] {
        refuses(&record_change(&text, id, |_| replacement), expected);
    }
}

#[test]
fn actual_part21_nonuniform_knot_records_preserve_interiors_and_refuse_reversed_intervals() {
    let original = file(false);
    let header = &original[..original.find("DATA;").unwrap() + 5];
    for (knots, expected) in [
        ("(0.,0.25,1.)", None),
        ("(0.,1.25,1.)", Some(Error::Invalid)),
    ] {
        let text = format!(
            "{header}#1=B_SPLINE_CURVE_WITH_KNOTS('',1,(#2,#3,#4),.UNSPECIFIED.,.F.,.F.,(2,1,2),{knots},.UNSPECIFIED.);ENDSEC;END-ISO-10303-21;"
        );
        let entities = part21::parse(&text).unwrap();
        let mut work = Ledger {
            used: 0,
            remaining: 16_000_000,
            interrupt: None,
        };
        let curve = schema::curve(1, simple(&entities, 1).unwrap(), &mut work);
        assert_eq!(curve.err(), expected);
    }
}
