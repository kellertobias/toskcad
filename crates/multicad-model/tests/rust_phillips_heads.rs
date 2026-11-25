//! Actual canonical Phillips sources; STEP admission remains an explicit gate.
#[path = "rust_phillips_heads/export.rs"]
mod export;
use multicad_core::{
    commands::{Actor, Channel, Session},
    document::Document,
    model::{
        self, MechanicalPartParameters, MechanicalPartRequest, PhillipsHParameters, PhillipsHSize,
        PhillipsInnerForm, ScrewDrive, ScrewHead, ThreadFit, ThreadHand, ThreadRepresentation,
        angle, length,
    },
    project_file::{self, ProjectFile},
};
use multicad_kernel::{MeshOptions, PlaneFrame};
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
fn recipe(kind: usize, fit: ThreadFit, left: bool, plain: f64) -> MechanicalPartParameters {
    let (nominal, pitch, head, drive) = match kind {
        0 => (
            2.,
            0.4,
            [4., 1.6, 3.2, 2.5, 0.1],
            (
                PhillipsHSize::H0,
                [1.9, 0.3, 0.61, 0.81],
                7.,
                PhillipsInnerForm::Rounded {
                    radius: length(0.25),
                },
            ),
        ),
        2 => (
            8.,
            1.25,
            [16., 6., 13., 10.8, 0.4],
            (
                PhillipsHSize::H4,
                [9., 1., 3.465, 5.08],
                7.,
                PhillipsInnerForm::ShallowV {
                    angle: angle(153.),
                    e_width: length(2.415),
                },
            ),
        ),
        6 => (
            6.,
            1.,
            [12., 4.6, 10., 8.4, 0.25],
            (
                PhillipsHSize::H3,
                [6.9, 0.8, 2.395, 3.81],
                5.75,
                PhillipsInnerForm::ShallowV {
                    angle: angle(146.),
                    e_width: length(2.005),
                },
            ),
        ),
        _ => (
            3.,
            0.5,
            [5.6, 2.4, 3., 4.2, 0.1],
            (
                PhillipsHSize::H1,
                [3., 0.5, 0.955, 1.27],
                7.,
                PhillipsInnerForm::ShallowV {
                    angle: angle(138.),
                    e_width: length(0.435),
                },
            ),
        ),
    };
    let mut p = model::default_machine_screw();
    let MechanicalPartParameters::MachineScrew {
        underhead_length,
        shank_diameter,
        head: hd,
        drive: dr,
        thread,
        ..
    } = &mut p
    else {
        unreachable!()
    };
    *underhead_length = Some(length(8. + plain));
    *shank_diameter = length(nominal);
    let crowned = Box::new(model::CrownedScrewHeadParameters {
        diameter: length(head[0]),
        height: length(head[1]),
        crown_radius: length(head[2]),
        top_land_diameter: length(head[3]),
        neck_radius: length(head[4]),
    });
    *hd = match kind {
        5 => ScrewHead::Countersunk {
            profile: Box::new(model::CountersunkScrewHeadParameters {
                diameter: length(5.5),
                height: length(1.65),
                included_angle: angle(90.),
                edge_chamfer: length(0.05),
                neck_radius: length(0.1),
            }),
        },
        _ => ScrewHead::Pan { profile: crowned },
    };
    let (size, [m, r, b, g], beta, inner) = drive;
    *dr = ScrewDrive::PhillipsH {
        profile: Box::new(PhillipsHParameters {
            size,
            mouth_diameter: length(m),
            mouth_radius: length(r),
            b_width: length(b),
            g_diameter: length(g),
            tilt_angle: angle(beta),
            outer_flank_angle: angle(26.5),
            floor_angle: angle(28.),
            outer_corner_angle: angle(92.),
            inner,
        }),
    };
    thread.major_diameter = length(nominal);
    thread.pitch = length(pitch);
    thread.length = length(8.);
    thread.fit = fit;
    thread.metric_fit_engagement = if fit == ThreadFit::NominalBasic {
        None
    } else {
        Some(model::MetricFitEngagement::NormalGroupAssumed {})
    };
    thread.hand = if left {
        ThreadHand::Left
    } else {
        ThreadHand::Right
    };
    thread.representation = ThreadRepresentation::Modeled;
    p
}
fn world(frame: PlaneFrame, p: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| {
        frame.origin[i] + frame.u[i] * p[0] + frame.v[i] * p[1] + frame.normal()[i] * p[2]
    })
}
fn occupied(solid: &multicad_kernel::Solid, point: [f64; 3]) -> bool {
    if solid.rust_periodic_repaired_phillips_source().is_some() {
        use multicad_kernel::rust_kernel::phillips_drive::Classification;
        match solid
            .rust_periodic_repaired_phillips_classify_world_box(
                point.map(|x| [x, x]),
                Default::default(),
                None,
            )
            .unwrap()
            .location()
        {
            Classification::Inside => true,
            Classification::Outside => false,
            Classification::UnresolvedBoundary => panic!("separated physical probe must resolve"),
        }
    } else {
        solid.rust_thread_contains(point).unwrap()
    }
}
#[test]
fn actual_metric_fits_and_hands_create_manufacture_and_reopen_without_occ_or_dummy_drive() {
    for (kind, fit, left) in [
        (0, ThreadFit::NominalBasic, false),
        (0, ThreadFit::NominalBasic, true),
        (1, ThreadFit::NominalBasic, false),
        (1, ThreadFit::Metric6gExternal, false),
        (1, ThreadFit::Metric6gExternal, true),
        (2, ThreadFit::NominalBasic, false),
        (6, ThreadFit::Metric6gExternal, true),
        (6, ThreadFit::Metric6gExternal, false),
        (5, ThreadFit::Metric6gExternal, true),
    ] {
        let mut session = Session::new(Document::new("Owned Phillips workflow").unwrap()).unwrap();
        let request =
            MechanicalPartRequest::new("Metric Phillips screw", recipe(kind, fit, left, 1.));
        let candidate = model::prepare_mechanical_part(&session, 0, actor(), &request).unwrap();
        let compiled = Compiler::default()
            .compile_with_backend(
                candidate.prepared.document(),
                1,
                GeometryBackend::Rust,
                None,
            )
            .unwrap();
        let display = &compiled.snapshot().shapes[&request.feature];
        assert!(!display.mesh.indices.is_empty());
        assert!(display.edges.len() >= 76);
        let solid = compiled.solid(request.feature).unwrap();
        assert!(solid.uses_rust_geometry());
        let (frame, head_top) =
            if let Some((source, frame)) = solid.rust_periodic_repaired_phillips_source() {
                assert_eq!(source.thread().left_hand(), left);
                assert_eq!(source.thread().length_mm(), 8.);
                assert_ne!(
                    fit,
                    ThreadFit::Metric6gExternal,
                    "a rounded fitted root must retain its actual source law"
                );
                let top = source.cavity_mouth_expression_mm();
                (frame, top[0] + top[1])
            } else {
                let (source, frame) = solid.rust_phillips_screw_geometry().unwrap();
                assert_eq!(source.thread().spec().left_hand, left);
                assert_eq!(source.thread().spec().length_mm, 8.);
                if fit == ThreadFit::Metric6gExternal {
                    assert!(source.thread().spec().root_radius_mm.is_some());
                }
                (frame, source.head_top_z())
            };
        assert!(!occupied(solid, world(frame, [0., 0., head_top + 0.1])));
        assert!(occupied(solid, world(frame, [0., 0., 4.])));
        let mesh = solid.manufacturing_mesh(MeshOptions::default()).unwrap();
        assert!(!mesh.indices.is_empty());
        let volume = solid.volume_mm3().unwrap();
        assert!(volume > 0.);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("physical-profile.step");
        if solid.rust_periodic_repaired_phillips_source().is_some() {
            assert!(solid.supports_step_export());
            solid.write_step(&path).unwrap();
            let actual = std::fs::read_to_string(&path).unwrap();
            solid
                .verify_periodic_phillips_step_source(&actual, 1e-7, None)
                .unwrap()
                .unwrap();
            assert!(
                multicad_kernel::Solid::import_rust_step(
                    &actual,
                    multicad_core::document::EntityId::new(),
                    Default::default()
                )
                .unwrap()
                .rust_periodic_repaired_phillips_source()
                .is_some()
            );
        } else if solid.rust_authored_phillips_step_source().is_some() {
            assert!(solid.supports_step_export());
            solid.write_step(&path).unwrap();
            let actual = std::fs::read_to_string(&path).unwrap();
            let proof = solid
                .verify_authored_phillips_step_source(&actual, 1e-7, None)
                .unwrap()
                .unwrap();
            assert!(proof.maximum_world_position_error_mm <= 1e-7);
            assert!(proof.maximum_world_chart_derivative_error_mm <= 1e-7);
        } else {
            assert!(!solid.supports_step_export());
            assert!(solid.write_step(&path).is_err());
            assert!(!path.exists());
            std::fs::write(&path, b"existing accepted file").unwrap();
            assert!(solid.write_step(&path).is_err());
            assert_eq!(std::fs::read(&path).unwrap(), b"existing accepted file");
        }
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
            volume
        );
        let edit = model::prepare_mechanical_part_edit(
            &session,
            session.revision(),
            actor(),
            request.feature,
            &recipe(kind, fit, left, 2.),
        )
        .unwrap();
        let edited = compiler()
            .compile_with_backend(
                edit.prepared.document(),
                session.revision() + 1,
                GeometryBackend::Rust,
                None,
            )
            .unwrap();
        assert!(edited.solid(request.feature).unwrap().volume_mm3().unwrap() > volume);
        let original = session.document().clone();
        session.commit(edit.prepared).unwrap();
        session.undo(session.revision(), actor()).unwrap();
        assert_eq!(session.document(), &original);
        session.redo(session.revision(), actor()).unwrap();
        let redone = compiler()
            .compile_with_backend(
                session.document(),
                session.revision(),
                GeometryBackend::Rust,
                None,
            )
            .unwrap();
        assert_eq!(
            redone.solid(request.feature).unwrap().volume_mm3().unwrap(),
            edited.solid(request.feature).unwrap().volume_mm3().unwrap()
        );
    }
}
#[test]
fn invalid_phillips_host_edit_cannot_change_the_accepted_document() {
    let mut session = Session::new(Document::new("Phillips rollback").unwrap()).unwrap();
    let request = MechanicalPartRequest::new(
        "Accepted Phillips screw",
        recipe(1, ThreadFit::Metric6gExternal, false, 1.),
    );
    let candidate = model::prepare_mechanical_part(&session, 0, actor(), &request).unwrap();
    compiler()
        .compile_with_backend(
            candidate.prepared.document(),
            1,
            GeometryBackend::Rust,
            None,
        )
        .unwrap();
    session.commit(candidate.prepared).unwrap();
    let original = session.document().clone();
    let mut invalid = recipe(1, ThreadFit::Metric6gExternal, false, 1.);
    let MechanicalPartParameters::MachineScrew {
        drive: ScrewDrive::PhillipsH { profile },
        ..
    } = &mut invalid
    else {
        unreachable!()
    };
    profile.mouth_diameter = length(100.);
    if let Ok(edit) = model::prepare_mechanical_part_edit(
        &session,
        session.revision(),
        actor(),
        request.feature,
        &invalid,
    ) {
        assert!(
            compiler()
                .compile_with_backend(
                    edit.prepared.document(),
                    session.revision() + 1,
                    GeometryBackend::Rust,
                    None
                )
                .is_err()
        );
    }
    assert_eq!(session.document(), &original);
    let rebuilt = compiler()
        .compile_with_backend(
            session.document(),
            session.revision(),
            GeometryBackend::Rust,
            None,
        )
        .unwrap();
    let restored = rebuilt.solid(request.feature).unwrap();
    assert!(
        restored.rust_phillips_screw_geometry().is_some()
            || restored.rust_periodic_repaired_phillips_source().is_some()
    );
}
