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
    *head = if kind == 0 {
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
    *drive = ScrewDrive::HexSocket {
        across_flats: length(2.),
        depth: length(0.6),
    };
    thread.major_diameter = length(3.);
    thread.pitch = length(0.5);
    thread.length = length(2.);
    thread.fit = fit;
    thread.representation = ThreadRepresentation::Modeled;
    p
}

#[test]
fn canonical_profiled_heads_create_edit_export_and_reopen_without_native_geometry() {
    for (kind, fit, underhead) in [
        (0, ThreadFit::NominalBasic, 3.),
        (0, ThreadFit::NominalBasic, 2.),
        (2, ThreadFit::NominalBasic, 3.),
        (2, ThreadFit::NominalBasic, 2.),
        (0, ThreadFit::Metric6gExternal, 3.),
        (0, ThreadFit::Metric6gExternal, 2.),
        (0, ThreadFit::Metric6gExternal, 2.01),
        (2, ThreadFit::Metric6gExternal, 3.),
        (2, ThreadFit::Metric6gExternal, 2.),
        (2, ThreadFit::Metric6gExternal, 2.01),
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
        if fit == ThreadFit::Metric6gExternal {
            assert!(actual.thread().spec().major_radius_mm < 1.5);
            assert!(actual.thread().spec().root_radius_mm.is_some());
        }
        assert_eq!(actual.thread().spec().length_mm, 2.);
        if underhead <= 2.01 {
            assert!(actual.visible_thread_interval_mm()[0] > 0.);
        }
        let before = source.volume_mm3().unwrap();
        let initial_dir = tempfile::tempdir().unwrap();
        let initial_path = initial_dir.path().join("initial.step");
        source.write_step(&initial_path).unwrap();
        let initial = multicad_kernel::Solid::import_rust_step(
            &std::fs::read_to_string(initial_path).unwrap(),
            EntityId::new(),
            GeometryPolicy::default(),
        )
        .unwrap();
        let (decoded, _) = initial.rust_profiled_screw_geometry().unwrap();
        assert!(decoded.spec().head.is_none());
        assert!((decoded.thread().spec().length_mm - 2.).abs() < 1e-8);
        assert!((initial.volume_mm3().unwrap() - before).abs() < 1e-6);
        let mut initial_bytes = Vec::new();
        project_file::write(
            &ProjectFile::new(candidate.prepared.document().clone()).unwrap(),
            &mut initial_bytes,
            Default::default(),
        )
        .unwrap();
        let initial_reopened = project_file::read(
            &mut initial_bytes.as_slice(),
            Default::default(),
            &project_file::supported_capabilities(),
        )
        .unwrap()
        .project;
        let initial_rebuilt = compiler()
            .compile_with_backend(&initial_reopened.document, 0, GeometryBackend::Rust, None)
            .unwrap();
        assert_eq!(
            initial_rebuilt
                .solid(request.feature)
                .unwrap()
                .volume_mm3()
                .unwrap(),
            before
        );
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
fn modeled_metric_recesses_enter_the_thread_core_and_reopen_without_recipes() {
    for kind in 0..3 {
        let mut parameters = recipe(
            if kind == 0 { 0 } else { 2 },
            ThreadFit::Metric6gExternal,
            2.,
        );
        let MechanicalPartParameters::MachineScrew { head, drive, .. } = &mut parameters else {
            unreachable!()
        };
        let depth = if kind == 1 {
            *head = ScrewHead::SocketCap {
                diameter: length(5.5),
                height: length(2.2),
            };
            2.5
        } else if kind == 0 {
            1.9
        } else {
            2.1
        };
        *drive = ScrewDrive::HexSocket {
            across_flats: length(1.5),
            depth: length(depth),
        };
        let mut session = Session::new(Document::new("Thread core recess").unwrap()).unwrap();
        let request = MechanicalPartRequest::new("Deep metric socket", parameters);
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
        let floor = source.head_top_z() + source.recess().maximum_depth_mm();
        assert!(floor > 0.);
        assert!(!source.contains([0., 0., floor / 2.], 0.));
        assert!(source.contains([0., 0., floor + 0.2], 0.));
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("core.step");
        solid.write_step(&path).unwrap();
        let imported = multicad_kernel::Solid::import_rust_step(
            &std::fs::read_to_string(path).unwrap(),
            EntityId::new(),
            GeometryPolicy::default(),
        )
        .unwrap();
        assert!((solid.volume_mm3().unwrap() - imported.volume_mm3().unwrap()).abs() < 1e-6);
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
