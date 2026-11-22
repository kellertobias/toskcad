use super::*;
pub(in crate::geometry) fn part(left: bool, frame: PlaneFrame) -> Solid {
    Solid::rust_periodic_repaired_phillips_material(
        super::super::periodic::tests::source(left, 0.5),
        EntityId::from_uuid_bits(0xcad352),
        frame,
        Default::default(),
        None,
    )
    .unwrap()
}
#[test]
fn actual_native_periodic_step_writes_freshly_verifies_reopens_and_preserves_atomic_failures() {
    for left in [false, true] {
        let frame = if left {
            PlaneFrame {
                origin: [105., -20., 60.],
                u: [0.8, 0.6, 0.],
                v: [0.36, -0.48, -0.8],
            }
        } else {
            PlaneFrame::default()
        };
        let source = part(left, frame);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("screw.step");
        fs::write(&path, b"existing output").unwrap();
        let polls = std::cell::Cell::new(0);
        let observe = || {
            polls.set(polls.get() + 1);
            false
        };
        let now = std::time::Instant::now();
        let first = source
            .write_periodic_phillips_step_with_verification(&path, 1e-7, Some(&observe))
            .unwrap()
            .unwrap();
        eprintln!("Native fresh physical STEP LH{left} {:?}", now.elapsed());
        assert_eq!(first.represented_member_count, 1);
        assert_eq!(first.member_index, 0);
        assert!(first.decoded_volume_bounds_mm3[0] > 0.);
        assert!(first.maximum_world_position_error_mm <= first.verified_tolerance_mm);
        assert!(first.maximum_world_chart_derivative_error_mm <= first.verified_tolerance_mm);
        assert!(first.correspondence_cavity_cells > 0 && first.correspondence_fixed_work > 0);
        assert!(
            first.decoded_bounds_mm[0]
                .iter()
                .zip(first.decoded_bounds_mm[1])
                .all(|(a, b)| a < &b)
        );
        let text = fs::read_to_string(&path).unwrap();
        let fresh = source
            .verify_periodic_phillips_step_source(&text, 1e-7, None)
            .unwrap()
            .unwrap();
        assert!(first.export_candidate.unwrap().correspondence_cavity_cells > 0);
        let mut fresh_scope = first.clone();
        fresh_scope.export_candidate = None;
        assert_eq!(fresh_scope, fresh);
        let decoded = multicad_rust_kernel::step_read::read_periodic_phillips_material(
            &text,
            Default::default(),
            None,
        )
        .unwrap()
        .into_material();
        let reopened = Solid::rust_periodic_repaired_phillips_material(
            Rc::new(decoded),
            EntityId::from_uuid_bits(0xcad353),
            PlaneFrame::default(),
            Default::default(),
            None,
        )
        .unwrap();
        assert!(reopened.rust_periodic_repaired_phillips_source().is_some());
        assert!(reopened.volume_mm3().unwrap() > 0.);
        let last = polls.get() - 1;
        let count = std::cell::Cell::new(0);
        let stop = || {
            count.set(count.get() + 1);
            count.get() >= last
        };
        assert_eq!(
            source
                .write_periodic_phillips_step_with_verification(&path, 1e-7, Some(&stop))
                .unwrap()
                .unwrap_err()
                .kind,
            GeometryErrorKind::Cancelled
        );
        assert_eq!(fs::read_to_string(&path).unwrap(), text);
        assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 1);
        assert!(
            source
                .write_periodic_phillips_step_with_verification(&path, 1e-16, None)
                .unwrap()
                .is_err()
        );
        assert_eq!(fs::read_to_string(&path).unwrap(), text);
        assert!(
            source
                .verify_periodic_phillips_step_source(&text, 0., None)
                .unwrap()
                .is_err()
        );
        assert!(
            source
                .write_periodic_phillips_step_with_verification(
                    &dir.path().join("missing").join("screw.step"),
                    1e-7,
                    None
                )
                .unwrap()
                .is_err()
        );
    }
}

#[test]
fn actual_native_compound_members_keep_explicit_order_and_fresh_independent_receipts() {
    let a = part(false, PlaneFrame::default());
    let b = part(
        true,
        PlaneFrame {
            origin: [5., 12., 18.],
            u: [0., 1., 0.],
            v: [1., 0., 0.],
        },
    );
    let mut writer = multicad_rust_kernel::step::Writer::new();
    let first = a
        .append_periodic_phillips_step(&mut writer, None)
        .unwrap()
        .unwrap();
    let second = b
        .append_periodic_phillips_step(&mut writer, None)
        .unwrap()
        .unwrap();
    assert!(first.matched_entities > 0 && second.matched_entities > 0);
    let text = writer.finish().unwrap();
    for (index, source) in [(0, &a), (1, &b)] {
        let proof = source
            .verify_periodic_phillips_step_member(&text, index, 1e-7, None)
            .unwrap()
            .unwrap();
        assert_eq!(proof.represented_member_count, 2);
        assert_eq!(proof.member_index, index);
        assert!(proof.selection_work > 0 && proof.selection_query_peak_bytes > 0);
        assert!(proof.decoded_volume_bounds_mm3[0] > 0.);
        assert!(proof.maximum_world_position_error_mm <= 1e-7);
        assert!(
            source
                .verify_periodic_phillips_step_source(&text, 1e-7, None)
                .unwrap()
                .is_err()
        );
        assert!(
            source
                .verify_periodic_phillips_step_member(&text, 2, 1e-7, None)
                .unwrap()
                .is_err()
        );
        assert_eq!(
            source
                .verify_periodic_phillips_step_member(&text, 1 - index, 1e-7, None)
                .unwrap()
                .unwrap_err()
                .kind,
            GeometryErrorKind::Unsupported
        );
    }
}
#[test]
fn actual_written_stale_or_corrupted_geometry_cannot_replace_existing_destination() {
    let expected = part(false, PlaneFrame::default());
    let stale = part(true, PlaneFrame::default());
    let raw = stale.rust_periodic_repaired_phillips_source().unwrap().0;
    let stale_text = exchange::export(
        raw,
        occurrence(PlaneFrame::default()),
        Default::default(),
        None,
    )
    .unwrap();
    let root = expected.rust_periodic_repaired_phillips_source().unwrap().0;
    let prefix = std::cell::Cell::new(0usize);
    let observe = || {
        prefix.set(prefix.get() + 1);
        false
    };
    exchange::export(
        root,
        occurrence(PlaneFrame::default()),
        Default::default(),
        Some(&observe),
    )
    .unwrap();
    let dir = tempfile::tempdir().unwrap();
    let target = dir.path().join("existing.step");
    for (bytes, kind) in [
        (
            stale_text.as_str().as_bytes(),
            GeometryErrorKind::Unsupported,
        ),
        (&b"not STEP"[..], GeometryErrorKind::InvalidResult),
        (&b"\xff"[..], GeometryErrorKind::InvalidResult),
    ] {
        fs::write(&target, b"keep previous output").unwrap();
        let polls = std::cell::Cell::new(0);
        let changed = std::cell::Cell::new(false);
        let external_change = || {
            polls.set(polls.get() + 1);
            // First measure the genuine delegated candidate phase above. Only
            // near its end do we observe the filesystem; no production hooks
            // or fabricated source receipts are needed to reach actual I/O.
            if !changed.get() && polls.get() >= prefix.get().saturating_sub(2) {
                for entry in fs::read_dir(dir.path()).unwrap() {
                    let path = entry.unwrap().path();
                    if path
                        .file_name()
                        .unwrap()
                        .to_string_lossy()
                        .ends_with(".step-tmp")
                    {
                        fs::write(path, bytes).unwrap();
                        changed.set(true);
                        break;
                    }
                }
            }
            false
        };
        assert_eq!(
            expected
                .write_periodic_phillips_step_with_verification(
                    &target,
                    1e-7,
                    Some(&external_change)
                )
                .unwrap()
                .unwrap_err()
                .kind,
            kind
        );
        assert!(changed.get());
        assert_eq!(fs::read(&target).unwrap(), b"keep previous output");
        assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 1);
    }
    let directory = dir.path().join("directory.step");
    fs::create_dir(&directory).unwrap();
    fs::write(directory.join("keep"), b"directory content").unwrap();
    assert_eq!(
        expected
            .write_periodic_phillips_step_with_verification(&directory, 1e-7, None)
            .unwrap()
            .unwrap_err()
            .kind,
        GeometryErrorKind::InvalidResult
    );
    assert_eq!(
        fs::read(directory.join("keep")).unwrap(),
        b"directory content"
    );
    assert!(fs::read_dir(dir.path()).unwrap().all(|e| {
        !e.unwrap()
            .file_name()
            .to_string_lossy()
            .ends_with(".step-tmp")
    }));
}
