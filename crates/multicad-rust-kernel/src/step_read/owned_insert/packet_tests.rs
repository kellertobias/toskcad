use super::{
    tests::{argument_change, file, record_change},
    *,
};

#[test]
fn actual_packet_requires_coordinate_lists_and_both_finite_cap_supports() {
    let text = file(false);
    let p = parse(&text, Options::default(), None).unwrap();
    let point = p.points()[0].entity();
    let not_a_list = argument_change(&text, point, 1, "$");
    assert_eq!(
        parse(&not_a_list, Options::default(), None).err(),
        Some(Error::Invalid)
    );
    let missing_cap = record_change(&text, p.planes()[0].entity(), |record| {
        record.replacen("PLANE(", "UNSUPPORTED_SURFACE(", 1)
    });
    assert_eq!(
        parse(&missing_cap, Options::default(), None).err(),
        Some(Error::Unsupported)
    );
}

#[test]
fn actual_vertex_must_own_the_original_curve_endpoint_not_another_valid_point() {
    let text = file(false);
    let p = parse(&text, Options::default(), None).unwrap();
    let edge = &p.edges()[0];
    let curve = p.curve(edge.curve_entity()).unwrap();
    let first = curve.pole_entities()[0];
    let last = *curve.pole_entities().last().unwrap();
    assert_ne!(first, last);
    assert!(
        p.point(first)
            .unwrap()
            .coordinates()
            .iter()
            .zip(p.point(last).unwrap().coordinates())
            .any(|(a, b)| a.value() != b.value())
    );
    let wrong = argument_change(&text, edge.vertex_entities()[0], 1, &format!("#{last}"));
    assert_eq!(
        parse(&wrong, Options::default(), None).err(),
        Some(Error::Unsupported)
    );
    assert_eq!(
        p.vertex(edge.vertex_entities()[0]).unwrap().point_entity,
        first
    );
}
