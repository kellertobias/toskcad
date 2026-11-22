use super::*;
#[test]
fn cylindrical_socket_cap_retains_true_rim_bearing_and_optional_plain_shank() {
    for scale in [0.125, 1., 8.] {
        let (radius, height, shank) = (2.8 * scale, 2.4 * scale, 1.5 * scale);
        for plain in [0., 0.2 * scale, 3. * scale] {
            let head = Head::Cylindrical {
                radius_mm: radius,
                height_mm: height,
            };
            let source = Profile::new(head, shank, plain).unwrap();
            assert_eq!(
                source,
                Profile::new_for_thread_union(head, shank, plain).unwrap()
            );
            assert_eq!(source.top(), [radius, -height]);
            assert_eq!(source.end(), [shank, plain]);
            let walls = source.segments().collect::<Vec<_>>();
            assert_eq!(walls.len(), if plain == 0. { 2 } else { 3 });
            assert_eq!(walls[0].role, Role::Rim);
            assert_eq!(walls[1].role, Role::Bearing);
            assert!(walls.iter().all(|s| s.curve == Curve::Line));
            for pair in walls.windows(2) {
                assert_eq!(pair[0].end, pair[1].start);
            }
            let expected = PI * (radius * radius * height + shank * shank * plain);
            assert!((source.volume_mm3() - expected).abs() < 1e-11 * scale.powi(3));
            assert!((walls[0].area_mm2() - TAU * radius * height).abs() < 1e-11 * scale * scale);
            assert!(
                (walls[1].area_mm2() - PI * (radius * radius - shank * shank)).abs()
                    < 1e-11 * scale * scale
            );
            assert_eq!(source.radius_at(-height / 2.), Some(radius));
            if plain > 0. {
                assert_eq!(source.radius_at(plain / 2.), Some(shank));
            }
            assert_eq!(source.support([0., 0., -1.]), height);
            assert_eq!(source.support([0., 0., 1.]), plain);
            assert_eq!(source.support([1., 0., 0.]), radius);
            assert!(source.admits_socket(1.4 * scale, 1.9 * scale, 0.1 * scale));
            assert!(!source.admits_socket(1.4 * scale, height, 0.1 * scale));
        }
    }
}
#[test]
fn cylindrical_cap_preserves_precision_and_dimension_refusals() {
    for (radius, height, shank, plain) in [
        (1.5, 2.4, 1.5, 0.),
        (2.8, 0., 1.5, 0.),
        (f64::NAN, 2.4, 1.5, 0.),
        (2.8, 2.4, 0., 0.),
        (2.8, 2.4, 1.5, -1.),
        (2.8, 2.4, 1.5, f64::INFINITY),
        (2.8, 2.4, 1.5, 1e-20),
    ] {
        assert_eq!(
            Profile::new(
                Head::Cylindrical {
                    radius_mm: radius,
                    height_mm: height
                },
                shank,
                plain
            ),
            Err(ThreadError::Invalid)
        );
    }
    assert_eq!(
        Profile::new(
            Head::Cylindrical {
                radius_mm: 1e8,
                height_mm: 2.4
            },
            1.5,
            0.
        ),
        Err(ThreadError::Limit)
    );
}
