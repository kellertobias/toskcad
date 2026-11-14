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
