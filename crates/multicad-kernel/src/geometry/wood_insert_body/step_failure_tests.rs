//! Independent failures exercise original public source and actual file data.
use super::*;
fn source() -> Solid {
    Solid::rust_wood_insert(
        crate::geometry::wood_insert_world_tests::source(false),
        EntityId::new(),
        PlaneFrame {
            origin: [31., -17., 9.],
            u: [0.8, 0.6, 0.],
            v: [-0.36, 0.48, 0.8],
        },
        GeometryPolicy::default(),
        None,
    )
    .unwrap()
}
fn document(source: &Solid) -> String {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("original-mapped-insert.step");
    source.write_step(&path).unwrap();
    std::fs::read_to_string(path).unwrap()
}
#[test]
fn actual_finite_operator_origin_overflow_is_a_measured_placement_refusal() {
    let solid = source();
    let text = document(&solid);
    solid
        .verify_partitioned_step_source(&text, 1e-7, None)
        .unwrap()
        .unwrap();
    let file = multicad_rust_kernel::step_read::owned_insert::CompoundDocument::parse(
        &text,
        Default::default(),
        None,
    )
    .unwrap();
    let operator = file.represented_members()[0]
        .mapping_link()
        .unwrap()
        .operator_entity();
    let record = text
        .lines()
        .find(|l| l.starts_with(&format!("#{operator}=")))
        .unwrap();
    let args = record
        .split("CARTESIAN_TRANSFORMATION_OPERATOR(")
        .nth(1)
        .unwrap()
        .split(')')
        .next()
        .unwrap();
    let origin = args.split(',').nth(2).unwrap().trim();
    let original = text
        .lines()
        .find(|l| l.starts_with(&format!("{origin}=")))
        .unwrap();
    let changed = text.replacen(
        original,
        &format!("{origin}=CARTESIAN_POINT('',(1.E308,0.,0.));"),
        1,
    );
    let error = solid
        .verify_partitioned_step_source(&changed, 1e-7, None)
        .unwrap()
        .unwrap_err();
    assert_eq!(error.kind, GeometryErrorKind::Unsupported);
    assert!(
        error.message.contains("Placed STEP enclosure overflow"),
        "{error:?}"
    );
    solid
        .verify_partitioned_step_source(&text, 1e-7, None)
        .unwrap()
        .unwrap();
}
#[test]
fn actual_pose_rows_spans_and_final_publication_honor_cancellation() {
    let solid = source();
    let text = document(&solid);
    let (body, _) = solid.wood_insert_body().unwrap();
    // Count ONLY actual retained row/span boundaries, without changing any
    // source, plan or callback protocol. The public observer below independently
    // records the complete verification's callbacks before choosing stops.
    let transport_calls = body
        .with_plan(solid.policy, None, |plan| {
            Ok(1 + plan
                .surface()
                .rows()
                .iter()
                .map(|row| row.knots().len())
                .sum::<usize>())
        })
        .unwrap();
    let calls = std::cell::Cell::new(0usize);
    solid
        .verify_partitioned_step_source(
            &text,
            1e-7,
            Some(&|| {
                calls.set(calls.get() + 1);
                false
            }),
        )
        .unwrap()
        .unwrap();
    let total = calls.get();
    assert!(total > transport_calls + 1);
    let first = total - transport_calls - 1; // Last callback belongs to the coordinator.
    for (stop, message) in [
        (first, "Placed STEP derivative verification superseded"),
        (first + 1, "Placed STEP derivative span superseded"),
        (total - 2, "Placed STEP verification superseded"),
    ] {
        calls.set(0);
        let error = solid
            .verify_partitioned_step_source(
                &text,
                1e-7,
                Some(&|| {
                    let n = calls.get();
                    calls.set(n + 1);
                    n >= stop
                }),
            )
            .unwrap()
            .unwrap_err();
        assert_eq!(error.kind, GeometryErrorKind::Cancelled);
        assert_eq!(error.message, message);
    }
    solid
        .verify_partitioned_step_source(&text, 1e-7, None)
        .unwrap()
        .unwrap();
}
#[test]
fn entry_and_ordinary_compound_keeps_source_known_mass_and_rejects_bad_actual_surface() {
    let policy = GeometryPolicy::default();
    let root = 1. / (4. * 3_f64.sqrt());
    let entry = Solid::rust_helical_thread_entry_chamfer(
        crate::HelicalThread {
            major_radius_mm: 3.,
            pitch_mm: 1.,
            length_mm: 6.6,
            depth_mm: 0.875 * 3_f64.sqrt() / 2. - root,
            crest_width_mm: 0.125,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            radial_taper: 0.,
            starts: 1,
            left_hand: false,
        },
        Some(root),
        0.3,
        None,
        EntityId::new(),
        PlaneFrame::default(),
        policy,
    )
    .unwrap();
    let sphere = Solid::rust_sphere_named(
        1.,
        PlaneFrame {
            origin: [20., 0., 0.],
            ..Default::default()
        },
        EntityId::new(),
        policy,
    )
    .unwrap();
    assert!(
        sphere
            .verify_partitioned_step_source("not STEP", 1e-7, None)
            .is_none()
    );
    let compound = Solid::compound_named(&[&sphere, &entry], EntityId::new(), policy).unwrap();
    assert!(!compound.has_wood_insert_step_source());
    assert!(compound.has_partitioned_step_source());
    let text = document(&compound);
    assert!(
        compound
            .verify_wood_insert_step_compound(&text, 1e-7, None)
            .is_none()
    );
    let proof = compound
        .verify_partitioned_step_source(&text, 1e-7, None)
        .unwrap()
        .unwrap();
    assert_eq!(proof.source_known_members.len(), 1);
    assert_eq!(proof.source_known_members[0].0, 1);
    assert!(proof.members.retained_members.is_empty());
    assert_eq!(proof.members.decoded_scalar_members.len(), 1);
    assert_eq!(proof.members.decoded_scalar_members[0].0, 0);
    let surface = text
        .lines()
        .find(|l| l.contains("=SPHERICAL_SURFACE("))
        .unwrap();
    let prefix = surface.rsplit_once(',').unwrap().0;
    let changed = text.replacen(surface, &format!("{prefix},0.);"), 1);
    let error = compound
        .verify_partitioned_step_source(&changed, 1e-7, None)
        .unwrap()
        .unwrap_err();
    assert_eq!(error.kind, GeometryErrorKind::InvalidResult);
    assert!(
        error.message.contains("Actual member decode failed"),
        "{error:?}"
    );
    compound
        .verify_partitioned_step_source(&text, 1e-7, None)
        .unwrap()
        .unwrap();
}
