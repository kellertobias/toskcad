use super::*;

fn policy() -> GeometryPolicy {
    GeometryPolicy::default()
}

fn cube(center: [f64; 3], size: f64) -> Solid {
    Solid::box_solid([size; 3], policy())
        .unwrap()
        .transformed(RigidTransform {
            translation_mm: center.map(|c| c - size / 2.),
            ..Default::default()
        })
        .unwrap()
}

/// A 10 x 10 x 5 block with a radius-2 boss standing on its top, axis
/// through (5, 5), reaching z = 8, read back as a shell without a chart.
fn chartless(pose: RigidTransform) -> Solid {
    let stock = Solid::box_solid([10., 10., 5.], policy()).unwrap();
    let boss = Solid::cylinder(
        PlaneFrame {
            origin: [5., 5., 4.],
            ..Default::default()
        },
        2.,
        4.,
        policy(),
    )
    .unwrap();
    let joined = stock
        .boolean(&boss, BooleanOperation::Union, policy())
        .unwrap()
        .unwrap()
        .transformed(pose)
        .unwrap();
    let (body, frame) = joined.warped_union_body().unwrap();
    assert!(body.chart.is_some());
    let solid =
        Solid::from_imported_warped_union(body.geometry.clone(), frame, EntityId::new(), policy())
            .unwrap();
    assert!(solid.warped_union_body().unwrap().0.chart.is_none());
    solid
}

fn relation(shell: &Solid, probe: &Solid) -> Option<Option<f64>> {
    shell
        .warped_shell_probe_boolean(probe, BooleanOperation::Intersect, policy())
        .map(|result| result.unwrap().map(|s| s.volume_mm3().unwrap()))
}

#[test]
fn probes_clear_of_every_support_take_the_section_winding() {
    let d = 0.002;
    let volume = d * d * d;
    let near = |a: f64, b: f64| assert!((a - b).abs() < 1e-15, "{a} != {b}");
    for pose in [
        RigidTransform::default(),
        RigidTransform {
            rotation: [[0.8, -0.36, 0.48], [0.6, 0.48, -0.64], [0., 0.8, 0.6]],
            translation_mm: [13., -7., 42.],
        },
    ] {
        let shell = chartless(pose);
        let place = |p: [f64; 3]| {
            let moved: [f64; 3] =
                std::array::from_fn(|i| dot(pose.rotation[i], p) + pose.translation_mm[i]);
            cube(moved, d)
        };
        // Inside the boss just within its wall: the section is a full
        // circle, followed as an exact ellipse arc.
        near(
            relation(&shell, &place([6.99, 5., 6.5])).unwrap().unwrap(),
            volume,
        );
        // Inside the block under the boss.
        near(
            relation(&shell, &place([1., 1., 2.])).unwrap().unwrap(),
            volume,
        );
        // Just outside the boss wall and above the block top.
        assert_eq!(relation(&shell, &place([7.01, 5., 6.5])), Some(None));
        // Outside everything.
        assert_eq!(relation(&shell, &place([20., 5., 2.])), Some(None));
        // A probe on the block top or the boss wall is unresolved.
        assert_eq!(relation(&shell, &place([1., 1., 5.])), None);
        assert_eq!(relation(&shell, &place([7., 5., 6.5])), None);
    }
}

#[test]
fn proved_probes_share_the_unchanged_operand() {
    let shell = chartless(RigidTransform::default());
    let inside = cube([1., 1., 2.], 0.01);
    let outside = cube([20., 5., 2.], 0.01);
    let run = |a: &Solid, b: &Solid, operation| {
        a.warped_shell_probe_boolean(b, operation, policy())
            .unwrap()
            .unwrap()
            .map(|s| s.volume_mm3().unwrap())
    };
    let whole = shell.volume_mm3().unwrap();
    assert_eq!(run(&shell, &inside, BooleanOperation::Union), Some(whole));
    assert_eq!(run(&inside, &shell, BooleanOperation::Subtract), None);
    assert_eq!(
        run(&shell, &outside, BooleanOperation::Subtract),
        Some(whole)
    );
    assert_eq!(run(&outside, &shell, BooleanOperation::Intersect), None);
    // A cut into the shell, or a disjoint union, changes the boundary.
    assert!(
        shell
            .warped_shell_probe_boolean(&inside, BooleanOperation::Subtract, policy())
            .is_none()
    );
    assert!(
        shell
            .warped_shell_probe_boolean(&outside, BooleanOperation::Union, policy())
            .is_none()
    );
    // Neither operand retained as a warped shell.
    assert!(
        inside
            .warped_shell_probe_boolean(&outside, BooleanOperation::Intersect, policy())
            .is_none()
    );
}
