use multicad_core::{
    commands::{Actor, Channel, Session},
    document::{Document, EntityId},
    model::{
        self, MechanicalPartParameters, MechanicalPartRequest, ScrewDrive, ScrewHead, ThreadFit,
        ThreadRepresentation, angle, length,
    },
    project_file::{self, ProjectFile},
};
use multicad_kernel::GeometryPolicy;
use multicad_model::{Compiler, GeometryBackend};

fn actor() -> Actor {
    Actor::local(Channel::Ui)
}
fn compiler() -> Compiler {
    Compiler {
        omit_display: true,
        ..Default::default()
    }
}
fn recipe(kind: usize, fit: ThreadFit, underhead: f64) -> MechanicalPartParameters {
    let mut p = model::default_machine_screw();
    let MechanicalPartParameters::MachineScrew {
        underhead_length,
        shank_diameter,
        head,
        drive,
        thread,
        ..
    } = &mut p
    else {
        unreachable!()
    };
    *underhead_length = Some(length(underhead));
    *shank_diameter = length(3.);
    *head = if kind == 3 {
        ScrewHead::SocketCap {
            diameter: length(6.),
            height: length(2.5),
        }
    } else if kind == 0 {
        ScrewHead::Countersunk {
            profile: Box::new(model::CountersunkScrewHeadParameters {
                diameter: length(6.),
                height: length(1.65),
                included_angle: angle(90.),
                edge_chamfer: length(0.1),
                neck_radius: length(0.1),
            }),
        }
    } else {
        let profile = Box::new(model::CrownedScrewHeadParameters {
            diameter: length(6.),
            height: length(if kind == 1 { 2.5 } else { 1.8 }),
            crown_radius: length(4.),
            top_land_diameter: length(3.),
            neck_radius: length(0.1),
        });
        if kind == 1 {
            ScrewHead::Pan { profile }
        } else {
            ScrewHead::Button { profile }
        }
    };
    *drive = ScrewDrive::Hexalobular {
        profile: Box::new(model::HexalobularDriveParameters {
            major_diameter: length(2.),
            minor_diameter: length(1.5),
            outer_lobe_radius: length(0.2),
            depth: length(0.25),
            bottom_angle: angle(140.),
        }),
    };
    thread.major_diameter = length(3.);
    thread.pitch = length(0.5);
    thread.length = length(2.);
    thread.fit = fit;
    thread.representation = ThreadRepresentation::Modeled;
    p
}

#[test]
fn canonical_hexalobular_heads_create_edit_export_and_reopen_without_native_geometry() {
    for (kind, fit, underhead) in [
        (1, ThreadFit::NominalBasic, 3.),
        (1, ThreadFit::Metric6gExternal, 2.),
        (0, ThreadFit::NominalBasic, 2.),
        (0, ThreadFit::NominalBasic, 3.),
        (2, ThreadFit::NominalBasic, 2.),
        (2, ThreadFit::NominalBasic, 3.),
        (3, ThreadFit::NominalBasic, 2.),
        (3, ThreadFit::NominalBasic, 3.),
        (0, ThreadFit::Metric6gExternal, 2.),
        (2, ThreadFit::Metric6gExternal, 2.),
        (3, ThreadFit::Metric6gExternal, 2.),
    ] {
        let mut session = Session::new(Document::new("Profiled head workflow").unwrap()).unwrap();
        let request =
            MechanicalPartRequest::new("Metric profiled screw", recipe(kind, fit, underhead));
        let candidate = model::prepare_mechanical_part(&session, 0, actor(), &request).unwrap();
        let checked = compiler()
            .compile_with_backend(
                candidate.prepared.document(),
                1,
                GeometryBackend::Rust,
                None,
            )
            .unwrap();
        let source = checked.solid(request.feature).unwrap();
        assert!(source.uses_rust_geometry());
        let (actual, _) = source.rust_profiled_screw_geometry().unwrap();
        assert_eq!(actual.profile().end()[0], 1.5);
        assert!(actual.recess().flat_hex_geometry().is_none());
        let curved = actual.recess().hexalobular_geometry().unwrap();
        assert!(
            (curved.cone().tip_depth_mm - (0.25 + 1. / 70_f64.to_radians().tan())).abs() < 1e-14
        );
        if fit == ThreadFit::Metric6gExternal {
            assert!(actual.thread().spec().major_radius_mm < 1.5);
            assert!(actual.thread().spec().root_radius_mm.is_some());
        }
        let before = source.volume_mm3().unwrap();
        session.commit(candidate.prepared).unwrap();
        let original = session.document().clone();
        let mut edit = recipe(kind, fit, underhead);
        let MechanicalPartParameters::MachineScrew {
            underhead_length, ..
        } = &mut edit
        else {
            unreachable!()
        };
        *underhead_length = Some(length(4.));
        let candidate = model::prepare_mechanical_part_edit(
            &session,
            session.revision(),
            actor(),
            request.feature,
            &edit,
        )
        .unwrap();
        let edited = compiler()
            .compile_with_backend(
                candidate.prepared.document(),
                2,
                GeometryBackend::Rust,
                None,
            )
            .unwrap();
        let solid = edited.solid(request.feature).unwrap();
        assert!(solid.volume_mm3().unwrap() > before);
        assert_eq!(session.document(), &original);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("profiled.step");
        solid.write_step(&path).unwrap();
        let physical = multicad_kernel::Solid::import_rust_step(
            &std::fs::read_to_string(path).unwrap(),
            EntityId::new(),
            GeometryPolicy::default(),
        )
        .unwrap();
        assert!(
            physical
                .rust_profiled_screw_geometry()
                .unwrap()
                .0
                .spec()
                .head
                .is_none()
        );
        assert!((physical.volume_mm3().unwrap() - solid.volume_mm3().unwrap()).abs() < 1e-6);
        session.commit(candidate.prepared).unwrap();
        let mut bytes = Vec::new();
        project_file::write(
            &ProjectFile::new(session.document().clone()).unwrap(),
            &mut bytes,
            Default::default(),
        )
        .unwrap();
        let reopened = project_file::read(
            &mut bytes.as_slice(),
            Default::default(),
            &project_file::supported_capabilities(),
        )
        .unwrap()
        .project;
        let rebuilt = compiler()
            .compile_with_backend(&reopened.document, 0, GeometryBackend::Rust, None)
            .unwrap();
        assert_eq!(
            rebuilt
                .solid(request.feature)
                .unwrap()
                .volume_mm3()
                .unwrap(),
            solid.volume_mm3().unwrap()
        );
    }
}

#[test]
fn deep_metric_recess_retains_plain_and_thread_core_material_and_refuses_breakthrough() {
    for kind in [0, 2] {
        for fit in [ThreadFit::NominalBasic, ThreadFit::Metric6gExternal] {
            let mut session = Session::new(Document::new("Deep metric drive").unwrap()).unwrap();
            let mut parameters = recipe(kind, fit, 3.);
            let MechanicalPartParameters::MachineScrew { drive, .. } = &mut parameters else {
                unreachable!()
            };
            let ScrewDrive::Hexalobular { profile } = drive else {
                unreachable!()
            };
            profile.depth = length(1.8);
            let request = MechanicalPartRequest::new("Deep plain-shank drive", parameters.clone());
            let candidate = model::prepare_mechanical_part(&session, 0, actor(), &request).unwrap();
            let compiled = compiler()
                .compile_with_backend(
                    candidate.prepared.document(),
                    1,
                    GeometryBackend::Rust,
                    None,
                )
                .unwrap();
            let solid = compiled.solid(request.feature).unwrap();
            let (source, _) = solid.rust_profiled_screw_geometry().unwrap();
            let back = source.head_top_z() + source.recess().maximum_depth_mm();
            assert!(back > -1. && back < 0.);
            assert!(!source.contains([0., 0., back - 0.01], 0.));
            assert!(source.contains([0., 0., back + 0.01], 0.));
            assert!(source.contains([1.25, 0., back - 0.01], 0.));
            let mesh = solid
                .manufacturing_mesh_with_receipt(Default::default())
                .unwrap();
            assert_eq!(mesh.edge_pairs as usize * 2, mesh.mesh.indices.len());
            let directory = tempfile::tempdir().unwrap();
            let path = directory.path().join("deep-drive.step");
            solid.write_step(&path).unwrap();
            let physical = multicad_kernel::Solid::import_rust_step(
                &std::fs::read_to_string(path).unwrap(),
                EntityId::new(),
                GeometryPolicy::default(),
            )
            .unwrap();
            assert!((solid.volume_mm3().unwrap() - physical.volume_mm3().unwrap()).abs() < 1e-6);
            let physical_source = physical.rust_profiled_screw_geometry().unwrap().0;
            assert!(!physical_source.contains([0., 0., back - 0.01], 0.));
            assert!(physical_source.contains([0., 0., back + 0.01], 0.));
            session.commit(candidate.prepared).unwrap();
            let unchanged = session.document().clone();
            let MechanicalPartParameters::MachineScrew {
                underhead_length, ..
            } = &mut parameters
            else {
                unreachable!()
            };
            *underhead_length = Some(length(2.));
            let extended = model::prepare_mechanical_part_edit(
                &session,
                session.revision(),
                actor(),
                request.feature,
                &parameters,
            )
            .unwrap();
            let extended_compile = compiler()
                .compile_with_backend(extended.prepared.document(), 2, GeometryBackend::Rust, None)
                .unwrap();
            let extended_solid = extended_compile.solid(request.feature).unwrap();
            let actual = extended_solid.rust_profiled_screw_geometry().unwrap().0;
            let core_floor = actual.head_top_z() + actual.recess().maximum_depth_mm();
            assert!(core_floor > 0.);
            assert!(!actual.contains([0., 0., core_floor - 0.01], 0.));
            assert!(actual.contains([0., 0., core_floor + 0.01], 0.));
            assert!(actual.contains([1.1, 0., core_floor - 0.01], 0.));
            let mesh = extended_solid
                .manufacturing_mesh_with_receipt(Default::default())
                .unwrap();
            assert_eq!(mesh.edge_pairs as usize * 2, mesh.mesh.indices.len());
            let path = directory.path().join("thread-core-drive.step");
            extended_solid.write_step(&path).unwrap();
            let physical = multicad_kernel::Solid::import_rust_step(
                &std::fs::read_to_string(path).unwrap(),
                EntityId::new(),
                GeometryPolicy::default(),
            )
            .unwrap();
            assert!(
                (extended_solid.volume_mm3().unwrap() - physical.volume_mm3().unwrap()).abs()
                    < 1e-6
            );
            let decoded = physical.rust_profiled_screw_geometry().unwrap().0;
            assert!(!decoded.contains([0., 0., core_floor - 0.01], 0.));
            assert!(decoded.contains([0., 0., core_floor + 0.01], 0.));
            assert!(decoded.contains([1.1, 0., core_floor - 0.01], 0.));
            assert_eq!(session.document(), &unchanged);
            let MechanicalPartParameters::MachineScrew { drive, .. } = &mut parameters else {
                unreachable!()
            };
            let ScrewDrive::Hexalobular { profile } = drive else {
                unreachable!()
            };
            profile.depth = length(9.);
            assert!(
                model::prepare_mechanical_part_edit(
                    &session,
                    session.revision(),
                    actor(),
                    request.feature,
                    &parameters
                )
                .is_err()
            );
            assert_eq!(session.document(), &unchanged);
        }
    }
}
