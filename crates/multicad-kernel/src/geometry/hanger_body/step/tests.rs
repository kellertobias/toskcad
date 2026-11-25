use super::*;
fn solid(left: bool, nonbinary: bool) -> Solid {
    Solid::rust_hanger(
        crate::geometry::hanger_material::source_tests::source(left, nonbinary),
        EntityId::new(),
        PlaneFrame::default(),
        GeometryPolicy::default(),
        None,
    )
    .unwrap()
}
#[test]
fn actual_local_hanger_step_exports_cache_qualification_and_keep_healed_mass_distinct() {
    let dir = tempfile::tempdir().unwrap();
    for (left, nb) in [(false, false), (true, true)] {
        let source = solid(left, nb);
        let (body, frame) = source.hanger_body().unwrap();
        let started = std::time::Instant::now();
        assert!(source.supports_step_export());
        let qualified = body.qualified_step(frame, source.policy, None).unwrap();
        let first_seconds = started.elapsed().as_secs_f64();
        let started = std::time::Instant::now();
        assert!(source.supports_step_export());
        assert!(Rc::ptr_eq(
            &qualified,
            &body.qualified_step(frame, source.policy, None).unwrap()
        ));
        eprintln!(
            "Native hanger capability left={left} first={first_seconds} cached={}",
            started.elapsed().as_secs_f64()
        );
        let path = dir
            .path()
            .join(if left { "left.step" } else { "right.step" });
        source.write_step(&path).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        let expected =
            read::CompoundDocument::parse(qualified.as_str(), Default::default(), None).unwrap();
        let actual = read::CompoundDocument::parse(&text, Default::default(), None).unwrap();
        assert!(
            matches!(
                binding::bind_local_members(
                    &expected,
                    expected.represented_members()[0].item_entity(),
                    &actual,
                    actual.represented_members()[0].item_entity(),
                    binding::Options {
                        maximum_work: 64_000_000,
                        ..Default::default()
                    },
                    None,
                ),
                Err(binding::Error::WorkLimit)
            ),
            "explicit smaller callers retain their real refusal frontier"
        );
        let proof = source
            .verify_partitioned_step_source(&text, source.policy.linear_tolerance_mm, None)
            .unwrap()
            .unwrap();
        assert_eq!(proof.members.solid_count, 1);
        assert!(proof.members.retained_members.is_empty());
        assert!(proof.members.decoded_scalar_members.is_empty());
        assert!(proof.source_known_members.is_empty());
        assert!(proof.insert_correspondence.is_empty());
        assert_eq!(proof.healed_members.len(), 1);
        let (index, mass) = &proof.healed_members[0];
        eprintln!(
            "Native hanger written-file binding left={left} work={} entities={} bytes={}",
            mass.literal_binding_work,
            mass.matched_entities,
            text.len()
        );
        assert_eq!(*index, 0);
        assert!(mass.healed_volume_bounds_mm3[0] > 0.);
        assert_eq!(
            mass.healed_volume_bounds_mm3,
            qualified.healed_volume_bounds_mm3()
        );
        assert!(mass.source_to_healed_position_error_mm <= source.policy.linear_tolerance_mm);
        assert!(mass.matched_entities > 100);
        assert!(matches!(
            source
                .verify_partitioned_step_source(
                    &text,
                    source.policy.linear_tolerance_mm,
                    Some(&|| true)
                )
                .unwrap(),
            Err(GeometryError {
                kind: GeometryErrorKind::Cancelled,
                ..
            })
        ));
        let mut strict = source.share().unwrap();
        strict.policy.linear_tolerance_mm = 1e-10;
        assert!(!strict.supports_step_export());
        assert!(strict.write_step(&path).is_err());
        assert_eq!(std::fs::read_to_string(&path).unwrap(), text);
        let placed = source
            .transformed(RigidTransform {
                translation_mm: [1., 2., 3.],
                ..Default::default()
            })
            .unwrap();
        // A placed copy now writes the same local records through one map.
        assert!(placed.supports_step_export());
        placed.write_step(&path).unwrap();
        let moved = std::fs::read_to_string(&path).unwrap();
        assert_eq!(moved.matches("MAPPED_ITEM(").count(), 1);
        assert_eq!(
            moved.matches("ADVANCED_FACE(").count(),
            text.matches("ADVANCED_FACE(").count()
        );
    }
}
#[test]
fn actual_mixed_hanger_exports_keep_member_categories_and_reject_wrong_source() {
    let policy = GeometryPolicy::default();
    let source = solid(false, false);
    let ordinary = Solid::rust_box(
        EntityId::new().uuid_bits(),
        [1, 2, 3].map(|n| Rational::new(n, 1).unwrap()),
        PlaneFrame {
            origin: [31., 0., 0.],
            ..Default::default()
        },
        policy,
    )
    .unwrap();
    let mixed = Solid::make_rust_compound(
        &[&source, &ordinary, &source],
        EntityId::new(),
        policy,
        None,
    )
    .unwrap();
    let file = tempfile::NamedTempFile::new().unwrap();
    mixed.write_step(file.path()).unwrap();
    let text = std::fs::read_to_string(file.path()).unwrap();
    let proof = mixed
        .verify_partitioned_step_source(&text, policy.linear_tolerance_mm, None)
        .unwrap()
        .unwrap();
    assert_eq!(proof.members.solid_count, 3);
    assert!(proof.members.retained_members.is_empty());
    assert_eq!(
        proof
            .healed_members
            .iter()
            .map(|(index, _)| *index)
            .collect::<Vec<_>>(),
        vec![0, 2]
    );
    assert_eq!(proof.members.decoded_scalar_members.len(), 1);
    assert_eq!(proof.members.decoded_scalar_members[0].0, 1);
    assert!((proof.members.decoded_scalar_members[0].1 - 6.).abs() < 1e-10);
    let wrong = solid(true, false);
    let wrong_mixed =
        Solid::make_rust_compound(&[&wrong, &ordinary, &wrong], EntityId::new(), policy, None)
            .unwrap();
    assert!(matches!(
        wrong_mixed
            .verify_partitioned_step_source(&text, policy.linear_tolerance_mm, None)
            .unwrap(),
        Err(GeometryError {
            kind: GeometryErrorKind::InvalidResult,
            ..
        })
    ));
    assert_eq!(std::fs::read_to_string(file.path()).unwrap(), text);
}
#[test]
fn placed_hanger_step_holds_the_local_qualified_records_through_one_map_and_reads_back() {
    let local = solid(false, false);
    let frame = PlaneFrame {
        origin: [20., -7., 9.],
        u: [0., 1., 0.],
        v: [0., 0., 1.],
    };
    let placed = Solid::rust_hanger(
        crate::geometry::hanger_material::source_tests::source(false, false),
        EntityId::new(),
        frame,
        GeometryPolicy::default(),
        None,
    )
    .unwrap();
    assert!(placed.supports_step_export());
    let dir = tempfile::tempdir().unwrap();
    let [local_text, placed_text] = [(&local, "local"), (&placed, "placed")].map(|(s, n)| {
        let path = dir.path().join(format!("{n}.step"));
        s.write_step(&path).unwrap();
        std::fs::read_to_string(path).unwrap()
    });
    assert!(!local_text.contains("MAPPED_ITEM"));
    assert_eq!(placed_text.matches("MAPPED_ITEM(").count(), 1);
    // The placed file's local body is the same qualified records: equal
    // geometry entity counts, and the operator states the frame exactly.
    let count = |t: &str, name: &str| t.matches(name).count();
    for name in ["ADVANCED_FACE(", "EDGE_CURVE(", "VERTEX_POINT(", "PCURVE("] {
        assert_eq!(
            count(&local_text, name),
            count(&placed_text, name),
            "{name}"
        );
    }
    assert!(placed_text.contains("CARTESIAN_TRANSFORMATION_OPERATOR_3D("));
    let document = read::CompoundDocument::parse(&placed_text, Default::default(), None).unwrap();
    let [member] = document.represented_members() else {
        panic!("one placed member")
    };
    assert!(member.mapping_link().is_some());
    // The mapped member reads back through its decoded operator.
    let (body, frame) = placed.hanger_body().unwrap();
    let proof = body
        .verify_step_member(&document, *member, frame, placed.policy, None)
        .unwrap();
    let (local_body, local_frame) = local.hanger_body().unwrap();
    let local_document =
        read::CompoundDocument::parse(&local_text, Default::default(), None).unwrap();
    let local_proof = local_body
        .verify_step_member(
            &local_document,
            local_document.represented_members()[0],
            local_frame,
            local.policy,
            None,
        )
        .unwrap();
    assert!(proof.source_to_healed_position_error_mm <= placed.policy.linear_tolerance_mm);
    assert!(
        proof.source_to_healed_position_error_mm >= local_proof.source_to_healed_position_error_mm
    );
    for k in 0..2 {
        assert!(
            proof.source_support_derivative_error_mm[k]
                >= local_proof.source_support_derivative_error_mm[k]
        );
        assert!(proof.source_support_derivative_error_mm[k] <= 1e-6);
    }
    // A rigid operator keeps the healed volume to rounding.
    let [lo, hi] = proof.healed_volume_bounds_mm3;
    let [llo, lhi] = local_proof.healed_volume_bounds_mm3;
    assert!(lo <= lhi && hi >= llo && hi - lo <= (lhi - llo) + 1e-9);
    // The world enclosure contains the placed source box.
    let world = placed.bounds().unwrap();
    for k in 0..3 {
        assert!(proof.correspondence_enclosure_mm[0][k] <= world.min[k] + 1e-9);
        assert!(proof.correspondence_enclosure_mm[1][k] >= world.max[k] - 1e-9);
    }
    // One micrometre of occurrence shift is refused, and so are the local
    // file for the placed source and the placed file for the local source.
    let moved = PlaneFrame {
        origin: [frame.origin[0] + 1e-3, frame.origin[1], frame.origin[2]],
        ..frame
    };
    let refusal = body
        .verify_step_member(&document, *member, moved, placed.policy, None)
        .unwrap_err();
    assert_eq!(refusal.kind, GeometryErrorKind::InvalidResult);
    assert!(refusal.message.contains("placement"));
    let refusal = body
        .verify_step_member(
            &local_document,
            local_document.represented_members()[0],
            frame,
            placed.policy,
            None,
        )
        .unwrap_err();
    assert_eq!(refusal.kind, GeometryErrorKind::InvalidResult);
    let refusal = local_body
        .verify_step_member(&document, *member, local_frame, local.policy, None)
        .unwrap_err();
    assert_eq!(refusal.kind, GeometryErrorKind::InvalidResult);
    // A reflected operator (negated third axis) is refused, not transported.
    let operator3 = placed_text
        .split("CARTESIAN_TRANSFORMATION_OPERATOR_3D(#")
        .nth(1)
        .and_then(|t| t.split(')').next())
        .unwrap();
    let marker = format!("#{operator3}=DIRECTION('',(");
    let start = placed_text.find(&marker).unwrap() + marker.len();
    let end = start + placed_text[start..].find("))").unwrap();
    let negated = placed_text[start..end]
        .split(',')
        .map(|v| {
            let v = v.trim();
            v.strip_prefix('-')
                .map_or_else(|| format!("-{v}"), str::to_owned)
        })
        .collect::<Vec<_>>()
        .join(",");
    let reflected = format!("{}{negated}{}", &placed_text[..start], &placed_text[end..]);
    let reflected_document =
        read::CompoundDocument::parse(&reflected, Default::default(), None).unwrap();
    assert!(
        body.verify_step_member(
            &reflected_document,
            reflected_document.represented_members()[0],
            frame,
            placed.policy,
            None,
        )
        .is_err()
    );
}
#[test]
fn written_local_and_placed_hangers_read_back_with_their_healed_volume() {
    let dir = tempfile::tempdir().unwrap();
    let source = solid(false, false);
    let (body, frame) = source.hanger_body().unwrap();
    let qualified = body.qualified_step(frame, source.policy, None).unwrap();
    let [lo, hi] = qualified.healed_volume_bounds_mm3();
    let rotation = [[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]];
    for (name, transform) in [
        ("local", None),
        (
            "placed",
            Some(RigidTransform {
                rotation,
                translation_mm: [20., -7., 9.],
            }),
        ),
    ] {
        let shape = match transform {
            Some(t) => source.transformed(t).unwrap(),
            None => source.share().unwrap(),
        };
        let path = dir.path().join(format!("{name}.step"));
        shape.write_step(&path).unwrap();
        let back = Solid::read_step(&path, Default::default()).unwrap();
        assert_eq!(back.solid_count().unwrap(), 1);
        // The independent readback earns the same healed interval the
        // qualified export proved on its own bytes, in any declared pose.
        let (bounds, sewing) = back.healed_import_proof().unwrap();
        assert_eq!(bounds, [lo, hi]);
        assert!(sewing <= source.policy.linear_tolerance_mm);
        let volume = back.volume_mm3().unwrap();
        assert!(lo <= volume && volume <= hi && hi - lo <= 1e-6);
        // The radial envelope encloses the placed source material.
        let expected = shape.bounds().unwrap();
        let actual = back.bounds().unwrap();
        for i in 0..3 {
            assert!(
                actual.min[i] <= expected.min[i] + 1e-9,
                "{name} {actual:?} {expected:?}"
            );
            assert!(
                actual.max[i] >= expected.max[i] - 1e-9,
                "{name} {actual:?} {expected:?}"
            );
            let width = |b: &Bounds| b.max[i] - b.min[i];
            assert!(
                width(&actual) <= width(&expected) * 1.1 + 1e-9,
                "{name} {actual:?} {expected:?}"
            );
        }
        let moved = back
            .transformed(RigidTransform {
                translation_mm: [1., 2., 3.],
                ..Default::default()
            })
            .unwrap();
        assert_eq!(moved.volume_mm3().unwrap(), volume);
        assert!((moved.bounds().unwrap().min[2] - actual.min[2] - 3.).abs() < 1e-9);
        // No boundary is rebuilt: display, faces and re-export stay refusals.
        assert_eq!(
            back.mesh(Default::default()).err().unwrap().kind,
            GeometryErrorKind::Unsupported
        );
        assert!(back.faces().is_err());
        assert!(!back.supports_step_export());
    }
    // Only the declared rigid unit-scale map is read; any other operator is
    // refused with the original trim refusal rather than a guessed pose.
    let placed = std::fs::read_to_string(dir.path().join("placed.step")).unwrap();
    for tampered in [
        placed.replacen(
            ",1.)CARTESIAN_TRANSFORMATION_OPERATOR_3D",
            ",2.)CARTESIAN_TRANSFORMATION_OPERATOR_3D",
            1,
        ),
        placed.replacen("MAPPED_ITEM('',", "MAPPED_ITEM('',#999999999,", 1),
    ] {
        assert_ne!(tampered, placed);
        assert!(Solid::import_rust_step(&tampered, EntityId::new(), Default::default()).is_err());
    }
    let scaled = placed.replacen(
        ",1.)CARTESIAN_TRANSFORMATION_OPERATOR_3D",
        ",2.)CARTESIAN_TRANSFORMATION_OPERATOR_3D",
        1,
    );
    let refusal = Solid::import_rust_step(&scaled, EntityId::new(), Default::default())
        .err()
        .unwrap();
    assert_eq!(refusal.kind, GeometryErrorKind::Unsupported);
    assert!(refusal.message.contains("trim agreement"), "{refusal:?}");
}
