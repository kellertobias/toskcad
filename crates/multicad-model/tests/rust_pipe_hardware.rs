//! Canonical pipe hardware must retain its original thread and neck/mouth datum.
use multicad_core::{
    commands::{Actor, Channel, Session},
    document::Document,
    model::{
        self, MechanicalPartParameters, MechanicalPartRequest, NptGaugeConvention, NptThreadForm,
        NptThreadIntent, NptThreadSize, NptTruncationTarget, PipeGaugeConvention, PipeThreadForm,
        PipeThreadIntent, PipeThreadSize, PrimitiveFrame, ScrewDrive, ScrewHead, ThreadMode,
        ThreadParameters, length,
    },
    parameters,
    project_file::{self, ProjectFile},
};
use multicad_kernel::MeshOptions;
use multicad_model::{Compiler, GeometryBackend};
fn actor() -> Actor {
    Actor::local(Channel::Ui)
}
fn recipe(npt: bool, head: bool, posed: bool) -> MechanicalPartParameters {
    let mut thread: ThreadParameters = if npt {
        NptThreadIntent {
            form: if head {
                NptThreadForm::External
            } else {
                NptThreadForm::Internal
            },
            size: NptThreadSize::OneEighth,
            gauge: NptGaugeConvention::Nbs1959PlainFace,
            truncation: NptTruncationTarget::MidpointConeParallel,
        }
        .parameters()
    } else {
        PipeThreadIntent {
            form: if head {
                PipeThreadForm::R
            } else {
                PipeThreadForm::Rc
            },
            size: PipeThreadSize::OneEighth,
            gauge: PipeGaugeConvention::Iso7Nominal1994,
        }
        .parameters()
    };
    if posed {
        thread.frame = PrimitiveFrame {
            origin: [length(11.), length(-7.), length(5.)],
            u: [0., 1., 0.],
            v: [0., 0., 1.],
        };
    }
    if !head {
        return MechanicalPartParameters::HexNut {
            across_flats: length(18.),
            thread: Box::new(thread),
        };
    }
    let mut source = model::default_machine_screw();
    let MechanicalPartParameters::MachineScrew {
        underhead_length,
        shank_diameter,
        head,
        drive,
        thread: target,
        ..
    } = &mut source
    else {
        unreachable!()
    };
    *underhead_length = Some(thread.length.clone());
    *shank_diameter = thread.major_diameter.clone();
    *head = ScrewHead::ExternalHex {
        across_flats: length(18.),
        height: length(3.),
    };
    *drive = ScrewDrive::ExternalHex;
    **target = thread;
    source
}
fn check_original_pipe_hardware(npt: bool) {
    let compiler = Compiler {
        omit_display: true,
        ..Default::default()
    };
    for head in [false, true] {
        for posed in [false, true] {
            let label = format!("npt={npt} head={head} posed={posed}");
            let mut session = Session::new(Document::new(&label).unwrap()).unwrap();
            let source = recipe(npt, head, posed);
            let request = MechanicalPartRequest::new(&label, source.clone());
            let prepared = model::prepare_mechanical_part(&session, 0, actor(), &request).unwrap();
            session.commit(prepared.prepared).unwrap();
            let model = compiler
                .compile_with_backend(
                    session.document(),
                    session.revision(),
                    GeometryBackend::Rust,
                    None,
                )
                .unwrap_or_else(|e| panic!("{label}: {e}"));
            let solid = model.solid(request.feature).unwrap();
            assert!(solid.uses_rust_geometry());
            assert_eq!(solid.solid_count().unwrap(), 1);
            let volume = solid.volume_mm3().unwrap();
            assert!(volume > 0. && volume.is_finite());
            let thread = match &source {
                MechanicalPartParameters::MachineScrew { thread, .. }
                | MechanicalPartParameters::HexNut { thread, .. } => thread,
                _ => unreachable!(),
            };
            let evaluated = parameters::evaluate(session.document()).unwrap();
            let resolved = thread
                .resolve(
                    &evaluated,
                    if head {
                        ThreadMode::NewExternal
                    } else {
                        ThreadMode::CutInternal
                    },
                )
                .unwrap();
            let bounds = solid.bounds().unwrap();
            let axis = if posed { 0 } else { 2 };
            let origin = if posed { 11. } else { 0. };
            assert!(
                (bounds.min[axis] - (origin - if head { 3. } else { 0. })).abs() < 1e-5,
                "{label}: {bounds:?}"
            );
            assert!(
                (bounds.max[axis] - (origin + resolved.length_mm)).abs() < 1e-5,
                "{label}: {bounds:?}"
            );
            let receipt = solid
                .manufacturing_mesh_with_receipt(MeshOptions {
                    linear_deflection_mm: 0.05,
                    angular_deflection_rad: 0.3,
                    maximum_triangles: 500_000,
                })
                .unwrap_or_else(|e| panic!("{label}: {e}"));
            assert!(receipt.edge_pairs > 0 && receipt.vertex_links > 0);
            let mut edited = source.clone();
            match &mut edited {
                MechanicalPartParameters::MachineScrew {
                    head: ScrewHead::ExternalHex { across_flats, .. },
                    ..
                }
                | MechanicalPartParameters::HexNut { across_flats, .. } => {
                    *across_flats = length(20.)
                }
                _ => unreachable!(),
            }
            let prepared =
                model::prepare_mechanical_part_edit(&session, 1, actor(), request.feature, &edited)
                    .unwrap();
            session.commit(prepared.prepared).unwrap();
            let changed = compiler
                .compile_with_backend(session.document(), 2, GeometryBackend::Rust, None)
                .unwrap();
            let edited_volume = changed
                .solid(request.feature)
                .unwrap()
                .volume_mm3()
                .unwrap();
            assert!(edited_volume > volume);
            session.undo(2, actor()).unwrap();
            assert_eq!(
                model::read_mechanical_part(session.document(), request.feature)
                    .unwrap()
                    .parameters,
                source
            );
            session.redo(3, actor()).unwrap();
            let saved = session.document().clone();
            let mut bytes = Vec::new();
            project_file::write(
                &ProjectFile::new(saved.clone()).unwrap(),
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
            assert_eq!(reopened.document, saved);
            let cold = compiler
                .compile_with_backend(&reopened.document, 4, GeometryBackend::Rust, None)
                .unwrap();
            assert!(
                (cold.solid(request.feature).unwrap().volume_mm3().unwrap() - edited_volume).abs()
                    < 1e-7
            );
        }
    }
}

#[test]
fn original_npt_hex_hardware_keeps_physical_pose_manufacture_history_and_cold_reopen() {
    check_original_pipe_hardware(true);
}
#[test]
fn original_bspt_hex_hardware_keeps_physical_pose_manufacture_history_and_cold_reopen() {
    check_original_pipe_hardware(false);
}

#[test]
fn compiler_owned_prepared_pipe_source_reuses_af_and_frame_edits() {
    use multicad_model::EvaluationCache;
    let compiler = Compiler {
        omit_display: true,
        ..Default::default()
    };
    let mut session = Session::new(Document::new("Owner prepared BSPT").unwrap()).unwrap();
    let original = recipe(false, false, false);
    let request = MechanicalPartRequest::new("BSPT", original);
    session
        .commit(
            model::prepare_mechanical_part(&session, 0, actor(), &request)
                .unwrap()
                .prepared,
        )
        .unwrap();
    let mut cache = EvaluationCache::new(32 * 1024 * 1024, 16);
    let first = compiler
        .compile_cached_with_backend(
            &mut cache,
            session.document(),
            1,
            GeometryBackend::Rust,
            None,
        )
        .unwrap();
    assert_eq!(cache.pipe_thread_factory_stats().prepared, 1);
    assert_eq!(cache.pipe_thread_factory_stats().entries, 1);
    assert!(cache.pipe_thread_factory_stats().retained_bytes < 32 * 1024 * 1024);
    let mut edited = recipe(false, false, true);
    let MechanicalPartParameters::HexNut { across_flats, .. } = &mut edited else {
        unreachable!()
    };
    *across_flats = length(20.);
    session
        .commit(
            model::prepare_mechanical_part_edit(&session, 1, actor(), request.feature, &edited)
                .unwrap()
                .prepared,
        )
        .unwrap();
    let at = std::time::Instant::now();
    let changed = compiler
        .compile_cached_with_backend(
            &mut cache,
            session.document(),
            2,
            GeometryBackend::Rust,
            None,
        )
        .unwrap();
    eprintln!(
        "COMPILER_PREPARED_AF_FRAME ms={:.6} stats={:?}",
        at.elapsed().as_secs_f64() * 1000.,
        cache.pipe_thread_factory_stats()
    );
    assert_eq!(cache.pipe_thread_factory_stats().hits, 1);
    assert_eq!(cache.pipe_thread_factory_stats().prepared, 0);
    let direct = compiler
        .compile_with_backend(session.document(), 2, GeometryBackend::Rust, None)
        .unwrap();
    let a = changed.solid(request.feature).unwrap();
    let b = direct.solid(request.feature).unwrap();
    assert_eq!(a.bounds().unwrap(), b.bounds().unwrap());
    assert_eq!(a.volume_mm3().unwrap(), b.volume_mm3().unwrap());
    assert_eq!(
        a.topology_snapshot(request.feature, 2).unwrap(),
        b.topology_snapshot(request.feature, 2).unwrap()
    );
    assert!(first.solid(request.feature).unwrap().volume_mm3().unwrap() < a.volume_mm3().unwrap());
    cache.clear();
    assert_eq!(cache.pipe_thread_factory_stats().entries, 0);
    assert_eq!(cache.pipe_thread_factory_stats().retained_bytes, 0);
    assert_eq!(cache.estimated_bytes(), 0);
    assert!(a.volume_mm3().unwrap() > 0.);
}

#[test]
fn compiler_late_cancellation_aborts_pending_sources_and_preserves_previous_owner() {
    use multicad_model::EvaluationCache;
    use std::cell::Cell;
    let compiler = Compiler {
        omit_display: true,
        ..Default::default()
    };
    let mut session = Session::new(Document::new("Prepared source cancellation").unwrap()).unwrap();
    let first = MechanicalPartRequest::new("first", recipe(true, true, false));
    session
        .commit(
            model::prepare_mechanical_part(&session, 0, actor(), &first)
                .unwrap()
                .prepared,
        )
        .unwrap();
    let mut observed = EvaluationCache::default();
    let mut cancelled = EvaluationCache::default();
    for cache in [&mut observed, &mut cancelled] {
        compiler
            .compile_cached_with_backend(cache, session.document(), 1, GeometryBackend::Rust, None)
            .unwrap();
        assert_eq!(cache.pipe_thread_factory_stats().entries, 1);
    }
    let second = MechanicalPartRequest::new("second", recipe(true, true, true));
    session
        .commit(
            model::prepare_mechanical_part(&session, 1, actor(), &second)
                .unwrap()
                .prepared,
        )
        .unwrap();
    let ticks = Cell::new(0usize);
    let observe = || {
        ticks.set(ticks.get() + 1);
        false
    };
    compiler
        .compile_cached_with_backend(
            &mut observed,
            session.document(),
            2,
            GeometryBackend::Rust,
            Some(&observe),
        )
        .unwrap();
    assert_eq!(observed.pipe_thread_factory_stats().entries, 2);
    let total = ticks.get();
    ticks.set(0);
    let stop = || {
        ticks.set(ticks.get() + 1);
        ticks.get() == total
    };
    assert!(
        compiler
            .compile_cached_with_backend(
                &mut cancelled,
                session.document(),
                2,
                GeometryBackend::Rust,
                Some(&stop)
            )
            .is_err()
    );
    assert_eq!(cancelled.pipe_thread_factory_stats().entries, 1);
    assert_eq!(cancelled.pipe_thread_factory_stats().pending, 0);
    let retry = compiler
        .compile_cached_with_backend(
            &mut cancelled,
            session.document(),
            2,
            GeometryBackend::Rust,
            None,
        )
        .unwrap();
    assert!(retry.solid(first.feature).unwrap().volume_mm3().unwrap() > 0.);
    assert!(retry.solid(second.feature).unwrap().volume_mm3().unwrap() > 0.);
    assert_eq!(cancelled.pipe_thread_factory_stats().pending, 0);
}
