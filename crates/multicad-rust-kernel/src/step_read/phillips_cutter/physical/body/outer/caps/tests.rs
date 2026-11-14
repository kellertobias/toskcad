use super::*;
#[test]
fn literal_circle_support_enclosure_contains_independent_point_distance_and_bounds_work() {
    let mut w = Work {
        remaining: 1_000_000,
        interrupt: None,
        cancelled: Cell::new(false),
    };
    let circle = Circle {
        frame: support::Frame::literal(&LiteralAxis {
            origin: [0.; 3],
            normal: [0., 0., 1.],
            reference: [1., 0., 0.],
        })
        .unwrap(),
        radius: support::adjacent(1.),
    };
    let delta = 1e-5;
    let poles = vec![[1. + delta, 0., delta]; DEGREE + 1];
    let actual = (((1. + delta) - 1_f64).powi(2) + delta.powi(2)).sqrt();
    let bound = circle.distance(&poles, &mut w).unwrap();
    assert!(bound >= actual && bound < 3. * delta);
    assert!(circle.distance(&[], &mut w).is_err());
    assert!(
        circle
            .distance(&[[f64::INFINITY; 3]; DEGREE + 1], &mut w)
            .is_err()
    );
    w.remaining = 0;
    assert!(circle.distance(&poles, &mut w).is_err());
    assert_eq!(choose(22, 11), 705432.);
}
#[test]
fn malformed_literal_circle_schema_is_refused_before_any_support_receipt() {
    let source = super::super::super::super::tests::fixture(1, false, false, 1);
    let text =
        crate::step::write_phillips_screw(&source, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7)
            .unwrap();
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let shell = Shell::new(&body.entities, 1e-7, true);
    let (_, forward, bounds) = helical::face(&shell, body.tensor_face_entities()[3]).unwrap();
    let wire = helical::orient(&helical::ring(&shell, bounds[0], true).unwrap(), forward);
    let curve = reference(&shell.args(wire[0].0, "EDGE_CURVE").unwrap()[3]).unwrap();
    let old = text
        .lines()
        .find(|l| l.starts_with(&format!("#{curve}=")))
        .unwrap();
    let changed = text.replacen(old, &format!("#{curve}=CIRCLE('');"), 1);
    // Actual public structural extraction refuses this malformed file.
    assert!(read_phillips_body(&changed, 1e-7, 1_000_000).is_err());
    // The shared primitive independently refuses malformed raw schema too;
    // no fabricated expected cutter or inferred gauge participates.
    let parsed = crate::step_read::part21::parse(&changed).unwrap();
    assert!(Caps::new(&body, &Shell::new(&parsed, 1e-7, true)).is_err());
}
