use super::*;

fn policy() -> GeometryPolicy {
    GeometryPolicy::default()
}

fn block(size: [f64; 3], at: [f64; 3]) -> Solid {
    Solid::box_solid(size, policy())
        .unwrap()
        .transformed(RigidTransform {
            translation_mm: at,
            ..Default::default()
        })
        .unwrap()
}

fn compound(members: &[Solid]) -> Solid {
    let members: Vec<&Solid> = members.iter().collect();
    Solid::make_rust_compound(&members, EntityId::new(), policy(), None).unwrap()
}

fn common(a: &Solid, b: &Solid) -> Option<f64> {
    a.compound_intersection(b, BooleanOperation::Intersect, policy())
        .unwrap()
        .unwrap()
        .map(|s| s.volume_mm3().unwrap())
}

#[test]
fn a_solid_meets_a_compound_member_by_member() {
    let posts = compound(&[
        block([2., 2., 10.], [0., 0., 0.]),
        block([2., 2., 10.], [20., 0., 0.]),
    ]);
    // Clear of both, through one, across both.
    let near = |a: Option<f64>, b: f64| assert!((a.unwrap() - b).abs() < 1e-9, "{a:?} != {b}");
    assert_eq!(common(&posts, &block([4., 4., 4.], [8., 0., 0.])), None);
    near(common(&block([4., 4., 4.], [-1., -1., 3.]), &posts), 16.);
    let across = common(&posts, &block([30., 1., 1.], [-5., 0.5, 4.]));
    near(across, 4.);
    // Only intersections distribute; plain pairs stay with other providers.
    let plain = block([1., 1., 1.], [0.; 3]);
    assert!(
        posts
            .compound_intersection(&plain, BooleanOperation::Subtract, policy())
            .is_none()
    );
    assert!(
        plain
            .compound_intersection(&plain, BooleanOperation::Intersect, policy())
            .is_none()
    );
}
