use super::*;
fn fixture() -> Solid {
    let policy = GeometryPolicy::default();
    let feature = EntityId::new();
    let stock = orthogonal::boolean(
        &[
            Box3 {
                min: [0.; 3],
                max: [40., 40., 5.],
            },
            Box3 {
                min: [0.; 3],
                max: [40., 5., 40.],
            },
        ],
        &[],
        Operation::Union,
        policy.linear_tolerance_mm,
    )
    .unwrap()
    .unwrap();
    let stock = Solid::rust_orthogonal_stock_named(feature, stock, policy).unwrap();
    let tools = [
        [[20., 20., 5.1], [0., 0., -5.2], [1., 0., 0.], [0., -1., 0.]],
        [[20., 5.1, 20.], [0., -5.2, 0.], [1., 0., 0.], [0., 0., 1.]],
    ]
    .into_iter()
    .enumerate()
    .map(
        |(i, [origin, displacement_mm, u, v])| crate::CircularPocketTool {
            frame: PlaneFrame { origin, u, v },
            center_uv: [0.; 2],
            radius_mm: 4.9,
            displacement_mm,
            origin: Origin {
                source: feature,
                role: format!("angle.bore.{i}"),
            },
        },
    )
    .collect::<Vec<_>>();
    let drilled = stock
        .rust_drill_circular_machining(&tools, &[true, true], policy)
        .unwrap();
    drilled
        .rust_countersink_orthogonal(
            &[0, 1].map(|bore| OrthogonalCountersink {
                geometry: Seat {
                    bore,
                    end: false,
                    radius_mm: 8.02,
                    depth_mm: 3.12,
                },
                origin: Origin {
                    source: feature,
                    role: format!("angle.seat.{bore}"),
                },
            }),
            policy,
            None,
        )
        .unwrap()
}
fn actual(text: &str) -> String {
    let path = std::env::temp_dir().join(format!(
        "multicad-angle-{}-{}.step",
        std::process::id(),
        EntityId::new()
    ));
    std::fs::write(&path, text).unwrap();
    let written = std::fs::read_to_string(&path).unwrap();
    std::fs::remove_file(path).unwrap();
    written
}
#[test]
fn qualified_native_standalone_and_compound_files_preserve_both_mass_categories() {
    for (u, v, origin) in [
        ([1., 0., 0.], [0., 1., 0.], [0.; 3]),
        ([0.8, 0.6, 0.], [-0.6, 0.8, 0.], [11., -7., 5.]),
    ] {
        let solid = fixture();
        let (body, _) = solid.orthogonal_seat_body().unwrap();
        let frame = PlaneFrame { origin, u, v };
        let policy = GeometryPolicy::default();
        body.validate_frame(frame, policy).unwrap();
        for compound in [false, true] {
            let mut writer = step::Writer::new();
            if compound {
                writer
                    .append_box(
                        [Rational::new(4, 1).unwrap(); 3],
                        [-100.; 3],
                        [1., 0., 0.],
                        [0., 1., 0.],
                    )
                    .unwrap();
            }
            body.append_step(&mut writer, frame, policy).unwrap();
            let text = actual(&writer.finish().unwrap());
            let document = read::CompoundDocument::parse(&text, Default::default(), None).unwrap();
            let member = document.represented_members()[usize::from(compound)];
            let proof = body
                .verify_step_member(&document, member, frame, policy, None)
                .unwrap();
            let source = body.world_volume_bounds(frame).unwrap();
            assert_eq!(proof.source_volume_bounds_mm3, source);
            assert!(proof.healed_volume_bounds_mm3[0] > 0.);
            assert!(proof.source_to_healed_position_error_mm < policy.linear_tolerance_mm);
            assert!(proof.sewing_correspondence_bound_mm < policy.linear_tolerance_mm);
            assert!(proof.original_coedge_derivative_error_mm < policy.linear_tolerance_mm);
            assert!(
                proof
                    .source_support_derivative_error_mm
                    .into_iter()
                    .all(|v| v < policy.linear_tolerance_mm)
            );
            assert!(proof.matched_entities > 100 && proof.literal_binding_work > 0);
            let bounds = body.world_bounds(frame).unwrap();
            for i in 0..3 {
                assert!(
                    proof.correspondence_enclosure_mm[0][i] <= bounds.min[i]
                        && proof.correspondence_enclosure_mm[1][i] >= bounds.max[i]
                );
            }
        }
    }
}
#[test]
fn warm_native_cache_requires_exact_pose_source_and_policy_and_preserves_cancelled_state() {
    let solid = fixture();
    let (body, frame) = solid.orthogonal_seat_body().unwrap();
    let policy = GeometryPolicy::default();
    assert!(body.qualified_step.borrow().is_none());
    assert_eq!(
        body.qualified_step(frame, policy, Some(&|| true))
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Cancelled
    );
    assert!(body.qualified_step.borrow().is_none());
    let first = body.qualified_step(frame, policy, None).unwrap();
    let start = std::time::Instant::now();
    let warm = body.qualified_step(frame, policy, None).unwrap();
    eprintln!("warm native seated capability={:?}", start.elapsed());
    assert!(Rc::ptr_eq(&first, &warm));
    let moved = PlaneFrame {
        origin: [11., -7., 5.],
        u: [0.8, 0.6, 0.],
        v: [-0.6, 0.8, 0.],
    };
    let moved_file = body.qualified_step(moved, policy, None).unwrap();
    assert!(!Rc::ptr_eq(&first, &moved_file));
    assert!(!first.matches_source(&body.source, transport::placement(moved).unwrap()));
    let changed_policy = GeometryPolicy {
        linear_tolerance_mm: policy.linear_tolerance_mm * 2.,
        ..policy
    };
    let tight = body.qualified_step(moved, changed_policy, None).unwrap();
    assert!(!Rc::ptr_eq(&moved_file, &tight));
    assert_eq!(
        body.qualified_step(moved, changed_policy, Some(&|| true))
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Cancelled
    );
    let stored = Rc::clone(&body.qualified_step.borrow().as_ref().unwrap().2);
    assert!(Rc::ptr_eq(&stored, &tight));
    let cached = body.qualified_step(moved, changed_policy, None).unwrap();
    assert!(Rc::ptr_eq(&cached, &tight));
}
#[test]
fn actual_wrong_source_literal_alias_and_mid_readback_cancellation_do_not_publish() {
    use std::cell::Cell;
    let solid = fixture();
    let (body, frame) = solid.orthogonal_seat_body().unwrap();
    let policy = GeometryPolicy::default();
    let file = body.qualified_step(frame, policy, None).unwrap();
    let original = "4.9000000000000004E0";
    let changed = file
        .as_str()
        .replacen(original, "4.9000000000000004000000000000001E0", 1);
    assert_ne!(changed, file.as_str());
    let document = read::CompoundDocument::parse(&changed, Default::default(), None).unwrap();
    assert_eq!(
        body.verify_step_member(
            &document,
            document.represented_members()[0],
            frame,
            policy,
            None
        )
        .err()
        .unwrap()
        .kind,
        GeometryErrorKind::InvalidResult
    );
    let document = read::CompoundDocument::parse(file.as_str(), Default::default(), None).unwrap();
    let member = document.represented_members()[0];
    let calls = Cell::new(0usize);
    body.verify_step_member(
        &document,
        member,
        frame,
        policy,
        Some(&|| {
            calls.set(calls.get() + 1);
            false
        }),
    )
    .unwrap();
    let stop = calls.get() - 4;
    calls.set(0);
    assert_eq!(
        body.verify_step_member(
            &document,
            member,
            frame,
            policy,
            Some(&|| {
                calls.set(calls.get() + 1);
                calls.get() > stop
            })
        )
        .err()
        .unwrap()
        .kind,
        GeometryErrorKind::Cancelled
    );
    assert!(
        body.verify_step_member(&document, member, frame, policy, None)
            .is_ok()
    );
    assert!(Rc::ptr_eq(
        &file,
        &body.qualified_step(frame, policy, None).unwrap()
    ));
}
