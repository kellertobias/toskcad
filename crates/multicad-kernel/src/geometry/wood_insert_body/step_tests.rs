use super::*;

#[test]
fn actual_native_retained_insert_exports_through_the_normal_solid_path() {
    use multicad_rust_kernel::step_read::owned_insert::{File, FluxOptions, Options};
    for left in [false, true] {
        let source = crate::geometry::wood_insert_world_tests::source(left);
        let solid = Solid::rust_wood_insert(
            source,
            EntityId::new(),
            PlaneFrame::default(),
            GeometryPolicy::default(),
            None,
        )
        .unwrap();
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("actual-insert.step");
        solid
            .write_step(&path)
            .expect("normal Native Solid must export its actual retained insert source");
        assert!(solid.supports_step_export());
        let text = std::fs::read_to_string(&path).unwrap();
        let file = File::parse(&text, Options::default(), None).unwrap();
        assert_eq!(file.represented_solids().len(), 1);
        let parsed = file
            .member(file.represented_solids()[0], Options::default(), None)
            .unwrap();
        let bindings = parsed.support_bindings(Options::default(), None).unwrap();
        let material = bindings.material_bounds(Options::default(), None).unwrap();
        let volume = material
            .occupied_volume(FluxOptions::default(), None)
            .unwrap()
            .volume_mm3();
        assert!(volume[0] > 0. && volume[1] >= volume[0]);
        assert!(volume[0] <= 231.37317864323353 && volume[1] >= 231.37317864323353);
        assert_eq!(solid.rust_wood_insert_data().unwrap().source, source);
    }
}

#[test]
fn actual_mixed_insert_export_checks_every_original_member() {
    let policy = GeometryPolicy::default();
    let insert = Solid::rust_wood_insert(
        crate::geometry::wood_insert_world_tests::source(false),
        EntityId::new(),
        PlaneFrame::default(),
        policy,
        None,
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
    let compound = Solid::compound_named(&[&sphere, &insert], EntityId::new(), policy).unwrap();
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("mixed.step");
    compound.write_step(&path).unwrap();
    let text = std::fs::read_to_string(path).unwrap();
    let proof = compound
        .verify_wood_insert_step_compound(&text, 1e-7, None)
        .unwrap()
        .unwrap();
    assert_eq!(proof.members.solid_count, 2);
    assert_eq!(proof.members.retained_members[0].member_index, 1);
    assert_eq!(proof.members.decoded_scalar_members.len(), 1);
    assert_eq!(proof.insert_correspondence[0].0, 1);
    let wrong_sphere = Solid::rust_sphere_named(
        2.,
        PlaneFrame {
            origin: [20., 0., 0.],
            ..Default::default()
        },
        EntityId::new(),
        policy,
    )
    .unwrap();
    let wrong = Solid::compound_named(&[&wrong_sphere, &insert], EntityId::new(), policy).unwrap();
    assert!(
        wrong
            .verify_wood_insert_step_compound(&text, 1e-7, None)
            .unwrap()
            .is_err()
    );
    assert!(
        insert
            .verify_wood_insert_step_compound(&text, 1e-7, None)
            .unwrap()
            .is_err()
    );
}

fn actual_mapped_document(local: &str, frame: PlaneFrame, scale: f64) -> String {
    let file =
        multicad_rust_kernel::step_read::owned_insert::File::parse(local, Default::default(), None)
            .unwrap();
    let solid = file.represented_solids()[0];
    let representation = local
        .lines()
        .find(|l| l.contains("=ADVANCED_BREP_SHAPE_REPRESENTATION("))
        .unwrap();
    let rep_id = representation.split('=').next().unwrap();
    let context = representation
        .rsplit(",#")
        .next()
        .unwrap()
        .trim_end_matches(");");
    let origin = local
        .lines()
        .find(|l| l.contains("=AXIS2_PLACEMENT_3D("))
        .unwrap()
        .split('=')
        .next()
        .unwrap();
    let xyz = |p: [f64; 3]| format!("({},{},{})", p[0], p[1], p[2]);
    let records = format!(
        "#1000000=CARTESIAN_POINT('',{});\n#1000001=DIRECTION('',{});\n#1000002=DIRECTION('',{});\n#1000003=DIRECTION('',{});\n#1000004=(CARTESIAN_TRANSFORMATION_OPERATOR(#1000001,#1000002,#1000000,{scale})CARTESIAN_TRANSFORMATION_OPERATOR_3D(#1000003)FUNCTIONALLY_DEFINED_TRANSFORMATION('','')GEOMETRIC_REPRESENTATION_ITEM()REPRESENTATION_ITEM(''));\n#1000005=SHAPE_REPRESENTATION('',(#{solid},{origin}),#{context});\n#1000006=REPRESENTATION_MAP({origin},#1000005);\n#1000007=MAPPED_ITEM('',#1000006,#1000004);\nENDSEC;\nEND-ISO-10303-21;",
        xyz(frame.origin),
        xyz(frame.u),
        xyz(frame.v),
        xyz(frame.normal())
    );
    local
        .replacen(
            representation,
            &format!("{rep_id}=ADVANCED_BREP_SHAPE_REPRESENTATION('',(#1000007),#{context});"),
            1,
        )
        .replacen("ENDSEC;\nEND-ISO-10303-21;", &records, 1)
}
#[test]
fn actual_mapped_insert_checks_raw_frame_correspondence_and_both_volume_authorities() {
    let policy = GeometryPolicy::default();
    for left in [false, true] {
        let source = crate::geometry::wood_insert_world_tests::source(left);
        let local =
            Solid::rust_wood_insert(source, EntityId::new(), PlaneFrame::default(), policy, None)
                .unwrap();
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("local.step");
        local.write_step(&path).unwrap();
        let text = std::fs::read_to_string(path).unwrap();
        let frame = PlaneFrame {
            origin: [11., -7., 5.],
            u: [0.6, 0.8, 0.],
            v: [0., 0., 1.],
        };
        let placed = Solid::rust_wood_insert(source, EntityId::new(), frame, policy, None).unwrap();
        let mapped = actual_mapped_document(&text, frame, 1.);
        let proof = placed
            .verify_wood_insert_step_compound(&mapped, 1e-7, None)
            .unwrap()
            .unwrap();
        if !left {
            let calls = std::cell::Cell::new(0usize);
            placed
                .verify_wood_insert_step_compound(
                    &mapped,
                    1e-7,
                    Some(&|| {
                        calls.set(calls.get() + 1);
                        false
                    }),
                )
                .unwrap()
                .unwrap();
            let total = calls.get();
            for stop in [0, total / 2, total - 1] {
                calls.set(0);
                let error = placed
                    .verify_wood_insert_step_compound(
                        &mapped,
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
            }
            // Cancellation only abandons a receipt; it never mutates source
            // ownership or the original file arena used by a later request.
            placed
                .verify_wood_insert_step_compound(&mapped, 1e-7, None)
                .unwrap()
                .unwrap();
        }
        let p = &proof.insert_correspondence[0].1;
        assert!(p.maximum_position_error_mm <= 1e-7);
        assert!(p.maximum_row_derivative_error_mm <= 1e-6);
        assert!(p.maximum_ruling_derivative_error_mm <= 1e-6);
        // The normalized STEP basis preserves the independently integrated
        // local file volume. Native's approximate source scalar is a separate
        // contract, not an oracle for the tighter decoded-file interval.
        assert!(p.decoded_volume_bounds_mm3[0] <= 231.37317864323353);
        assert!(p.decoded_volume_bounds_mm3[1] >= 231.37317864323353);
        assert!(p.source_volume_bounds_mm3[0] <= placed.volume_mm3().unwrap());
        assert!(p.source_volume_bounds_mm3[1] >= placed.volume_mm3().unwrap());
        let changed = actual_mapped_document(
            &text,
            PlaneFrame {
                origin: [11.01, -7., 5.],
                ..frame
            },
            1.,
        );
        assert!(
            placed
                .verify_wood_insert_step_compound(&changed, 1e-7, None)
                .unwrap()
                .is_err()
        );
        let scaled = actual_mapped_document(&text, frame, 1.001);
        assert!(
            placed
                .verify_wood_insert_step_compound(&scaled, 1e-7, None)
                .unwrap()
                .is_err()
        );
        let reflected = mapped.replace("DIRECTION('',(0.8,-0.6,0))", "DIRECTION('',(-0.8,0.6,0))");
        assert_ne!(mapped, reflected);
        assert!(
            placed
                .verify_wood_insert_step_compound(&reflected, 1e-7, None)
                .unwrap()
                .is_err()
        );
    }
}
