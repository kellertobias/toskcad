use super::*;
fn input(posed: bool) -> String {
    let g = crate::rounded_cylinder::RoundedCylinder::new(5., 12., [None, Some(3.)], 1e-6, 1e7)
        .unwrap();
    let (o, u, v) = if posed {
        ([17., -8., 13.], [0.8, 0.6, 0.], [-0.36, 0.48, 0.8])
    } else {
        ([0.; 3], [1., 0., 0.], [0., 1., 0.])
    };
    crate::step::write_rounded_cylinder(&g, o, u, v, 1e-7).unwrap()
}
#[test]
fn actual_upper_fillet_trim_is_regular_while_complete_spindle_body_guard_remains() {
    for posed in [false, true] {
        let text = input(posed);
        let source = Source::parse(&text, 1e-7).unwrap();
        let torus = source
            .entities
            .iter()
            .find(|(_, parts)| parts.iter().any(|e| e.name == "TOROIDAL_SURFACE"))
            .unwrap()
            .0;
        assert!(matches!(
            source.surface(*torus),
            Err(ReadError::Unsupported)
        ));
        let report = validate_pcurves(&text, 1e-7).unwrap();
        assert!(report.pcurves_checked > 0 && report.certified_maximum_deviation_mm <= 1e-7);
        // The support-only path cannot authorize a complete spindle retained B-rep.
        assert!(matches!(
            super::super::super::read_assembly(&text, 1e-7, None),
            Err(crate::brep::Error::Unsupported)
        ));
    }
}
#[test]
fn spindle_trim_support_preserves_resource_and_original_coordinate_guards() {
    let text = input(false);
    let source = Source::parse(&text, 1e-7).unwrap();
    let mut budget = Budget::new(MAX_WORK, None);
    let binding = bindings(&source, &mut budget)
        .unwrap()
        .into_iter()
        .find(|binding| {
            uv::definition(&source, binding.pcurve, &mut Budget::new(MAX_WORK, None)).is_ok_and(
                |d| {
                    source
                        .component(d.surface, "TOROIDAL_SURFACE")
                        .unwrap()
                        .is_some()
                },
            )
        })
        .unwrap();
    let d = uv::definition(&source, binding.pcurve, &mut budget).unwrap();
    assert!(read(&source, &d, &mut budget).is_ok());
    assert!(matches!(
        read(&source, &d, &mut Budget::new(0, None)),
        Err(ReadError::Limit)
    ));
    let cancellation = || true;
    let mut interrupted = Budget::new(MAX_WORK, Some(&cancellation));
    assert!(matches!(
        read(&source, &d, &mut interrupted),
        Err(ReadError::Limit)
    ));
    assert!(interrupted.cancelled.get());
    // Proof is also rejected when the original direction has a nonzero
    // component hidden by binary64 underflow, without rounding it into an iso.
    let [_, _, vector] = source.args(d.curve_id, "LINE").unwrap() else {
        panic!()
    };
    let [_, direction, _] = source.args(reference(vector).unwrap(), "VECTOR").unwrap() else {
        panic!()
    };
    let id = reference(direction).unwrap();
    let statement = source.simple(id).unwrap().source.unwrap();
    assert!(statement.starts_with("DIRECTION("));
    let changed = text.replace(statement, "DIRECTION('',(1E-1000,1.))");
    assert_ne!(changed, text);
    let changed_source = Source::parse(&changed, 1e-7).unwrap();
    assert!(changed_source.args(id, "DIRECTION").is_ok());
    assert!(validate_pcurves(&changed, 1e-7).is_err());
}
