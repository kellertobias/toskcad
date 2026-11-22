//! Independent argument contracts: every object comes from a real parser.
use super::*;
fn append_actual(text: &str, record: &str) -> String {
    let at = text.rfind("ENDSEC;").unwrap();
    format!("{}{}{}", &text[..at], record, &text[at..])
}
#[test]
fn actual_line_target_and_uv_line_do_not_satisfy_the_spline_window_operator() {
    let text = input("0.0625", "1.");
    let text = append_actual(
        &text,
        "\n#990100=PCURVE('',#31,#990101);\n#990101=DEFINITIONAL_REPRESENTATION('',(#990102),#900103);\n#990102=LINE('',#990110,#990111);\n#990110=CARTESIAN_POINT('',(0.,0.));\n#990111=VECTOR('',#990112,1.);\n#990112=DIRECTION('',(1.,0.));\n",
    );
    let source = Source::parse(&text, 1e-7).unwrap();
    let (d, e, packet) = setup(&source);
    let mut b = Budget::new(512_000_000, None);
    let line_uv = super::super::super::super::uv::definition(&source, 990100, &mut b).unwrap();
    assert!(matches!(
        prepare(&source, 900001, &line_uv, &e, &packet, &mut b),
        Err(ReadError::Unsupported)
    ));
    let (edge_id, curve_id) = source
        .entities
        .iter()
        .find_map(|(&id, parts)| {
            let [entity] = parts.as_slice() else {
                return None;
            };
            if entity.name != "EDGE_CURVE" {
                return None;
            }
            let curve = reference(&entity.args[3]).unwrap();
            matches!(source.curve(curve, 0).unwrap(), Curve::Line { .. }).then_some((id, curve))
        })
        .unwrap();
    let line_edge = super::super::super::super::edge(
        &source,
        super::super::super::super::Binding {
            edge: edge_id,
            curve: curve_id,
            pcurve: 900100,
        },
        &mut b,
    )
    .unwrap();
    // The packet is genuinely earned; a non-spline target cannot use it.
    assert!(matches!(
        prepare(&source, curve_id, &d, &line_edge, &packet, &mut b),
        Err(ReadError::Unsupported)
    ));
    assert!(prepare(&source, 900001, &d, &e, &packet, &mut b).is_ok());
}
#[test]
fn large_integer_uv_parameter_window_declines_without_a_certificate() {
    let original = input("0.0625", "1.");
    let text = original.replace(
        "(2,2),(0.25,0.75)",
        "(2,2),(9007199254740990.,9007199254740991.)",
    );
    assert_ne!(text, original);
    let source = Source::parse(&text, 1e-7).unwrap();
    let (d, e, packet) = setup(&source);
    let Curve::Nurbs { curve, .. } = &d.geometry else {
        panic!("actual parsed UV")
    };
    assert_eq!(curve.range(), [9007199254740990., 9007199254740991.]);
    // This UV parameter window lies outside the target's natural interval.
    // Check the actual refusal without claiming which numeric guard wins.
    assert!(matches!(
        prepare(
            &source,
            900001,
            &d,
            &e,
            &packet,
            &mut Budget::new(512_000_000, None)
        ),
        Err(ReadError::Unsupported)
    ));
}

#[test]
fn original_adjacent_decimal_uv_parameters_refuse_without_inventing_a_gap() {
    let original = input("0.0625", "1.");
    let text = original.replace("(2,2),(0.25,0.75)", "(2,2),(0.1,0.10000000000000002)");
    assert_ne!(text, original);
    let source = Source::parse(&text, 1e-7).unwrap();
    let (d, e, packet) = setup(&source);
    let Curve::Nurbs { curve, .. } = &d.geometry else {
        panic!("actual parsed UV")
    };
    assert!(curve.range()[1] > curve.range()[0]);
    // A prior numeric transfer may already refuse this non-dyadic case;
    // this control claims refusal, not a hit on a particular later guard.
    assert!(matches!(
        prepare(
            &source,
            900001,
            &d,
            &e,
            &packet,
            &mut Budget::new(512_000_000, None)
        ),
        Err(ReadError::Unsupported)
    ));
}
