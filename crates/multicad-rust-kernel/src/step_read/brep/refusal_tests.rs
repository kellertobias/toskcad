use super::*;

const PLACED: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../multicad-app/tests/fixtures/exchange/corpus/occt-placed-box.step"
));

#[test]
fn independent_refusal_propagates_an_interrupt_after_inspection_starts() {
    let polls = std::cell::Cell::new(0);
    let interrupt = || {
        polls.set(polls.get() + 1);
        polls.get() >= 3
    };
    assert_eq!(
        retained_geometry_refusal(PLACED, 1e-7, Some(&interrupt)),
        Some(crate::brep::Error::Cancelled)
    );
    assert_eq!(polls.get(), 3);
    // A later fresh request can still inspect the same immutable input.
    assert_eq!(retained_geometry_refusal(PLACED, 1e-7, None), None);
    assert!(read_assembly(PLACED, 1e-7, None).is_ok());
}

#[test]
fn missing_trim_certificate_never_admits_geometry_or_hides_a_proved_bad_solid() {
    // A zero-offset parameter-space curve has the original line geometry,
    // but this representation is outside the current trim certificates.
    let unqualified = PLACED.replace(
        "#38 = LINE('',#39,#40);",
        "#38 = OFFSET_CURVE_2D('',#900001,0.,.F.);\n#900001 = LINE('',#39,#40);",
    );
    assert_ne!(unqualified, PLACED);
    assert_eq!(
        validate_pcurves(&unqualified, 1e-7).unwrap_err().kind,
        ReadError::Unsupported
    );
    assert_eq!(retained_geometry_refusal(&unqualified, 1e-7, None), None);
    assert!(matches!(
        read_assembly(&unqualified, 1e-7, None),
        Err(crate::brep::Error::Unsupported)
    ));
    // Reverse exactly one face. Every entity remains owned, but the shell
    // no longer has consistent orientation, irrespective of its trims.
    let open = unqualified.replace(
        "#17 = ADVANCED_FACE('',(#18),#32,.F.);",
        "#17 = ADVANCED_FACE('',(#18),#32,.T.);",
    );
    assert_ne!(open, unqualified);
    assert_eq!(
        retained_geometry_refusal(&open, 1e-7, None),
        Some(crate::brep::Error::Invalid)
    );
    assert!(matches!(
        read_assembly(&open, 1e-7, None),
        Err(crate::brep::Error::Invalid)
    ));
    assert_eq!(
        retained_geometry_refusal(&unqualified, 1e-7, Some(&|| true)),
        Some(crate::brep::Error::Cancelled)
    );
    assert!(matches!(
        read_assembly(&unqualified, 1e-7, Some(&|| true)),
        Err(crate::brep::Error::Cancelled)
    ));
}
