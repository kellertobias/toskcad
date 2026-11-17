//! Normal typed STEP import keeps actual original source through save/cold reopen.
use multicad_core::{
    commands::{Actor, Channel, Session},
    document::{Document, EntityId},
    model::imports::{StepEncoding, read_step_source},
    project_file::{self, ProjectFile},
};
use multicad_kernel::{GeometryPolicy, Solid, rust_kernel as pure};
use multicad_mcp::{
    Headless, Storage, Workspace,
    tools::{Call, dispatch},
};
use multicad_model::GeometryBackend;
use pure::axial_pipe_thread::authored::{AuthoredSection, AuthoredSpec, AuthoredThread};
use pure::step::npt_stock::{CurveOptions, CurvePlan, SurfacePlan, write_source_candidate};
use pure::step_read::{
    npt_stock::{
        LiteralOccurrencePlacement, Parsed, Support,
        full_head::{OwnedActualHead, transport::encoding},
        head::{self, HeadInput},
    },
    owned_insert::{CompoundDocument, DocumentOptions},
};
use serde_json::{Value, json};
use std::sync::{Arc, atomic::AtomicBool};
fn original_clipped_file(quarter: bool) -> String {
    let (tpi, e1, l1, l2, od, max) = if quarter {
        (18., 0.49163, 0.2278, 0.4018, 0.540, 0.088)
    } else {
        (27., 0.37360, 0.1615, 0.2639, 0.405, 0.096)
    };
    let pitch = 25.4 / tpi;
    let k = 1. / 32.;
    let m = 3f64.sqrt();
    let slopes = 1. / (m - k) + 1. / (m + k);
    let depth = pitch / slopes - (0.033 + max) * pitch;
    let source = AuthoredThread::new(
        AuthoredSpec {
            major_radius_mm: (e1 * 25.4 - 2. * k * l1 * 25.4 + depth) / 2.,
            pitch_mm: pitch,
            length_mm: l2 * 25.4,
            radial_taper: k,
            section: AuthoredSection::NptFlat {
                depth_mm: depth,
                flat_width_mm: (pitch - depth * slopes) / 2.,
            },
            crest_radius_limit_mm: Some(od * 25.4 / 2.),
        },
        1e-7,
        1e-6,
    )
    .unwrap()
    .stock_clip_plan(Default::default(), None)
    .unwrap();
    let options = CurveOptions::default();
    let curves = CurvePlan::new(&source, options, None).unwrap();
    let plan = SurfacePlan::new(curves, options.into(), None).unwrap();
    write_source_candidate(&plan, options.maximum_work, options.maximum_bytes, None)
        .unwrap()
        .into_string()
}
fn physical_head_file(quarter: bool) -> String {
    let text = original_clipped_file(quarter);
    let doc = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let parsed = Parsed::parse(
        &doc,
        doc.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap();
    let cap = parsed
        .faces()
        .iter()
        .filter_map(|f| match f.support() {
            Support::Plane { height } => {
                Some((parsed.points()[height].coordinates()[2].value(), f.entity()))
            }
            _ => None,
        })
        .max_by(|a, b| a.0.total_cmp(&b.0))
        .unwrap()
        .1;
    let head = HeadInput::regular_hex(&parsed, cap, 18., 3.).unwrap();
    let proof = head::certify(&parsed, head, 1e-7, 1e-6, Default::default(), None).unwrap();
    pure::step::npt_stock::full_head::propose(&proof, Default::default(), None)
        .unwrap()
        .as_str()
        .to_owned()
}
fn source_file(quarter: bool, posed: bool) -> String {
    let local = Arc::new(
        OwnedActualHead::read(
            Arc::from(physical_head_file(quarter)),
            Default::default(),
            None,
        )
        .unwrap(),
    );
    let frame = LiteralOccurrencePlacement::new(
        [[1., 0., 0.], [0., -1., 0.], [0., 0., -1.]],
        [0., 0., local.head_stations_mm()[0]],
        Default::default(),
        None,
    )
    .unwrap();
    let frame = if posed {
        frame
            .transformed(
                [[0.6, 0., 0.8], [0., 1., 0.], [-0.8, 0., 0.6]],
                [20., -30., 11.],
                Default::default(),
                None,
            )
            .unwrap()
    } else {
        frame
    };
    encoding::propose_mapped(&local, &frame, Default::default(), None)
        .unwrap()
        .actual_bytes()
        .to_string()
}
struct Store;
impl Storage for Store {
    fn open(&mut self, path: &std::path::Path) -> Result<ProjectFile, String> {
        project_file::read(
            &mut std::fs::File::open(path).map_err(|e| e.to_string())?,
            Default::default(),
            &project_file::supported_capabilities(),
        )
        .map(|read| read.project)
    }
    fn save(&mut self, path: &std::path::Path, project: &ProjectFile) -> Result<(), String> {
        project_file::write(
            project,
            &mut std::fs::File::create(path).map_err(|e| e.to_string())?,
            Default::default(),
        )
        .map_err(|e| e.to_string())
    }
}
fn call(cancel: bool) -> Call {
    Call {
        actor: Actor {
            id: "mcp:original-pipe-hardware".into(),
            channel: Channel::Mcp,
        },
        cancel: Arc::new(AtomicBool::new(cancel)),
    }
}
fn run(host: &Headless, tool: &str, args: Value) -> Value {
    dispatch(host, &call(false), tool, args).unwrap()
}

#[test]
fn public_actual_head_import_typed_update_save_cold_reopen_and_atomic_failures() {
    for quarter in [false, true] {
        for posed in [false, true] {
            let dir = tempfile::tempdir().unwrap();
            let path = dir.path().join("actual-head.tcd");
            Store
                .save(
                    &path,
                    &ProjectFile::new(Document::new("Actual mapped NPT head").unwrap()).unwrap(),
                )
                .unwrap();
            let open = || {
                Headless::with_geometry_backend(
                    Workspace::project(path.clone(), Box::new(Store)).unwrap(),
                    Some(dir.path().to_path_buf()),
                    false,
                    GeometryBackend::Rust,
                )
            };
            let input = source_file(quarter, posed);
            let input_path = dir.path().join("actual-head.step");
            std::fs::write(&input_path, &input).unwrap();
            // Model's public file adapter earns actual geometry before the Core transaction.
            let verified = multicad_model::imports::read_step_file_with_backend(
                &input_path,
                GeometryBackend::Rust,
            )
            .unwrap();
            assert_eq!(verified.encoding, StepEncoding::RustStep);
            assert_eq!(verified.solids.len(), 1);
            assert_eq!(verified.solids[0].1, input);
            let host = open();
            let imported = run(
                &host,
                "import_step",
                json!({"expected_revision":0,"file_name":"actual-head.step"}),
            );
            assert_eq!(imported["revision"], 1);
            assert_eq!(imported["parts"].as_array().unwrap().len(), 1);
            let source = EntityId::parse(imported["source"].as_str().unwrap()).unwrap();
            let initial = host.document();
            let stored = read_step_source(&initial, source).unwrap();
            assert_eq!(stored.encoding, StepEncoding::RustStep);
            assert_eq!(stored.solids[0].brep, input);
            // A real file-operator edit keeps the SAME typed source and Part.
            // An unchanged import is intentionally revision-neutral, so it is
            // not a meaningful edit/reopen fixture.
            let input = source_file(quarter, !posed);
            std::fs::write(&input_path, &input).unwrap();
            let updated = run(
                &host,
                "import_step",
                json!({"expected_revision":1,"source":source,"file_name":"actual-head.step"}),
            );
            assert_eq!(updated["revision"], 2);
            assert_eq!(updated["source"], imported["source"]);
            let saved = host.document();
            let edited = read_step_source(&saved, source).unwrap();
            assert_eq!(edited.solids.len(), 1);
            assert_eq!(edited.solids[0].id, stored.solids[0].id);
            assert_eq!(edited.solids[0].brep, input);
            assert_ne!(edited.solids[0].brep, stored.solids[0].brep);
            assert_eq!(updated["matched"][0]["geometry_changed"], true);
            let model = run(&host, "evaluate_model", json!({}));
            assert_eq!(model["parts"][0]["solid_count"], 1);
            let volume = model["parts"][0]["volume_mm3"].as_f64().unwrap();
            assert!((volume - verified.solids[0].2.volume_mm3).abs() < 1e-7);
            // Malformed actual file has no recipe fallback and cannot replace last source.
            std::fs::write(
                &input_path,
                "ISO-10303-21; HEADER; ENDSEC; DATA; #1=BOGUS(); ENDSEC; END-ISO-10303-21;",
            )
            .unwrap();
            assert!(
                dispatch(
                    &host,
                    &call(false),
                    "import_step",
                    json!({"expected_revision":2,"source":source,"file_name":"actual-head.step"})
                )
                .is_err()
            );
            assert_eq!(host.document(), saved);
            std::fs::write(&input_path, &input).unwrap();
            assert!(
                dispatch(
                    &host,
                    &call(true),
                    "import_step",
                    json!({"expected_revision":2,"source":source,"file_name":"actual-head.step"})
                )
                .is_err()
            );
            assert_eq!(host.document(), saved);
            let output = dir.path().join("held.step");
            std::fs::write(&output, "prior output").unwrap();
            let exported = run(
                &host,
                "export_step",
                json!({"part":imported["parts"][0],"file_name":"held.step"}),
            );
            assert_eq!(
                exported["step_verification"]["method"],
                "actual_mapped_head_reference"
            );
            assert!(
                exported["readback_volume_mm3"].is_null(),
                "sewn reference is not occupied raw-file mass"
            );
            assert!(exported["step_verification"]["decoded_volume_mm3"].is_null());
            let evidence = &exported["step_verification"]["actual_head"];
            assert!(
                evidence["maximum_world_position_error_mm"]
                    .as_f64()
                    .unwrap()
                    <= 1e-7
            );
            assert!(evidence["aggregate_work"].as_u64().unwrap() <= 272_000_000);
            assert!(evidence["accounted_peak_bytes"].as_u64().unwrap() <= 64 * 1024 * 1024);
            let published = std::fs::read_to_string(&output).unwrap();
            assert_ne!(published, "prior output");
            let reopened = Solid::read_step(&output, GeometryPolicy::default()).unwrap();
            assert!((reopened.volume_mm3().unwrap() - volume).abs() < 1e-7);
            assert!(
                dispatch(
                    &host,
                    &call(true),
                    "export_step",
                    json!({"part":imported["parts"][0],"file_name":"held.step"})
                )
                .is_err()
            );
            assert_eq!(std::fs::read_to_string(&output).unwrap(), published);
            assert_eq!(host.document(), saved);
            run(&host, "save_project", json!({"expected_revision":2}));
            drop(host);
            std::fs::remove_file(&input_path).unwrap();
            let cold = open();
            assert_eq!(cold.document(), saved);
            assert_eq!(read_step_source(&cold.document(), source).unwrap(), edited);
            let rebuilt = run(&cold, "evaluate_model", json!({}));
            assert_eq!(
                rebuilt["parts"][0]["volume_mm3"],
                model["parts"][0]["volume_mm3"]
            );
            assert_eq!(rebuilt["parts"][0]["part"], model["parts"][0]["part"]);
            let cold_export = run(
                &cold,
                "export_step",
                json!({"part":imported["parts"][0],"file_name":"cold.step"}),
            );
            assert_eq!(
                cold_export["step_verification"]["method"],
                "actual_mapped_head_reference"
            );
            assert!(cold_export["readback_volume_mm3"].is_null());
            assert_eq!(
                std::fs::read_to_string(dir.path().join("cold.step")).unwrap(),
                published
            );
            // Reopening the canonical embedded source itself independently re-earns the
            // actual original file; the vanished exchange file is never a geometry cache.
            let actual = Solid::import_rust_step(
                &edited.solids[0].brep,
                EntityId::new(),
                GeometryPolicy::default(),
            )
            .unwrap();
            assert_eq!(
                actual
                    .rust_actual_clipped_head_source()
                    .unwrap()
                    .actual_file(),
                input
            );
            assert!(actual.supports_step_export());
            // The same source also survives a real canonical UI transform
            // feature, saved before a new Headless process compiles/exports it.
            let mut session = Session::new(cold.document()).unwrap();
            let body = EntityId::parse(imported["parts"][0].as_str().unwrap()).unwrap();
            let producer =
                EntityId::parse(rebuilt["parts"][0]["feature"].as_str().unwrap()).unwrap();
            let transform = multicad_core::model::TransformRequest::new(
                multicad_core::model::BodyInput { body, producer },
                multicad_core::model::NumericTransform {
                    translation: [3., -4., 5.].map(multicad_core::model::length),
                    axis: [0., 1., 0.],
                    angle: multicad_core::model::angle(37.),
                    ..Default::default()
                },
            );
            session
                .commit(
                    multicad_core::model::prepare_transform(
                        &session,
                        session.revision(),
                        Actor::local(Channel::Ui),
                        &transform,
                    )
                    .unwrap()
                    .prepared,
                )
                .unwrap();
            let moved_document = session.document().clone();
            drop(cold);
            Store
                .save(&path, &ProjectFile::new(moved_document.clone()).unwrap())
                .unwrap();
            let moved = open();
            assert_eq!(moved.document(), moved_document);
            assert_eq!(read_step_source(&moved.document(), source).unwrap(), edited);
            let moved_model = run(&moved, "evaluate_model", json!({}));
            assert_eq!(moved_model["parts"][0]["part"], imported["parts"][0]);
            let moved_export = run(
                &moved,
                "export_step",
                json!({"part":body,"file_name":"moved.step"}),
            );
            assert_eq!(
                moved_export["step_verification"]["method"],
                "actual_mapped_head_reference"
            );
            assert!(moved_export["readback_volume_mm3"].is_null());
            assert_eq!(moved.document(), moved_document);
            let moved_body =
                Solid::read_step(&dir.path().join("moved.step"), GeometryPolicy::default())
                    .unwrap();
            assert_eq!(
                moved_body
                    .rust_actual_clipped_head_source()
                    .unwrap()
                    .local_geometry()
                    .actual_file(),
                actual
                    .rust_actual_clipped_head_source()
                    .unwrap()
                    .local_geometry()
                    .actual_file()
            );

            eprintln!(
                "PUBLIC_HEAD_MCP_COLD quarter={quarter} initial_pose={posed} cold_pose={} revision=2 volume={volume}",
                !posed
            );
        }
    }
}
