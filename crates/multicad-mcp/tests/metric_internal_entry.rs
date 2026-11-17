#[path = "metric_single_internal_support/mod.rs"]
mod support;
#[test]
fn normal_mcp_metric_entry_only_receiving_cut_creates_actual_native_material() {
    support::require_positive(true);
}

#[test]
fn entry_only_normal_mcp_edits_undo_cancel_save_and_cold_reopen_preserve_actual_native_source() {
    use multicad_core::{
        commands::{Actor, Channel, Session},
        document::Document,
        model::{self, PrimitiveParameters, PrimitiveRequest, length},
        project_file::{self, ProjectFile},
    };
    use multicad_kernel::{MeshOptions, RigidTransform, rust_kernel::Rational};
    use multicad_mcp::{
        Headless, Storage, Workspace,
        tools::{Call, dispatch},
    };
    use multicad_model::{Compiler, GeometryBackend};
    use serde_json::{Value, json};
    use std::sync::{Arc, atomic::AtomicBool};
    struct Store;
    impl Storage for Store {
        fn open(&mut self, path: &std::path::Path) -> Result<ProjectFile, String> {
            project_file::read(
                &mut std::fs::File::open(path).map_err(|e| e.to_string())?,
                Default::default(),
                &project_file::supported_capabilities(),
            )
            .map(|p| p.project)
        }
        fn save(&mut self, path: &std::path::Path, p: &ProjectFile) -> Result<(), String> {
            project_file::write(
                p,
                &mut std::fs::File::create(path).map_err(|e| e.to_string())?,
                Default::default(),
            )
            .map_err(|e| e.to_string())
        }
    }
    let mm = |v: f64| json!({"value":v,"unit":"mm"});
    let ends =
        |depth: f64| json!({"start":{"chamfer45":{"radial_depth":mm(depth)}},"end":"flat_cut"});
    let call = |cancel| Call {
        actor: Actor::local(Channel::Mcp),
        cancel: Arc::new(AtomicBool::new(cancel)),
    };
    let run = |host: &Headless, name: &str, args: Value| {
        dispatch(host, &call(false), name, args.clone())
            .unwrap_or_else(|e| panic!("{name} {args}: {e:?}"))
    };
    for fit in ["nominal_basic", "metric_6h_internal"] {
        for hand in ["right", "left"] {
            let dir = tempfile::tempdir().unwrap();
            let path = dir.path().join("entry-only.tcd");
            let mut source =
                Session::new(Document::new("Actual MCP entry receiver").unwrap()).unwrap();
            let stock = PrimitiveRequest::new(
                "Cylinder stock",
                PrimitiveParameters::Cylinder {
                    radius: length(3.5),
                    height: length(3.3000000000000003),
                    frame: Default::default(),
                },
            );
            let prepared = model::prepare_primitive(&source, 0, call(false).actor, &stock).unwrap();
            source.commit(prepared.prepared).unwrap();
            Store
                .save(&path, &ProjectFile::new(source.document().clone()).unwrap())
                .unwrap();
            let open = || {
                Headless::with_geometry_backend(
                    Workspace::project(path.clone(), Box::new(Store)).unwrap(),
                    Some(dir.path().into()),
                    false,
                    GeometryBackend::Rust,
                )
            };
            let host = open();
            let initial = host.document();
            let args = json!({"expected_revision":0,"name":"M3 actual entry receiver","mode":"cut_internal","target_part":stock.body,"family":"iso_metric","fit":fit,"major_diameter":mm(3.),"pitch":mm(0.5),"length":mm(3.3000000000000003),"minimum_wall":mm(0.5),"hand":hand,"end_treatments":ends(0.4)});
            assert_eq!(
                dispatch(&host, &call(true), "set_thread", args.clone())
                    .unwrap_err()
                    .code,
                "cancelled"
            );
            assert_eq!(host.document(), initial);
            let mut invalid = args.clone();
            invalid["minimum_wall"] = mm(3.);
            assert!(dispatch(&host, &call(false), "set_thread", invalid).is_err());
            assert_eq!(host.document(), initial);
            let created = run(&host, "set_thread", args);
            let feature = serde_json::from_value(created["thread"].clone()).unwrap();
            let created_document = host.document();
            let compiled = (Compiler {
                omit_display: true,
                ..Default::default()
            })
            .compile_with_backend(&created_document, 1, GeometryBackend::Rust, None)
            .unwrap();
            let body = compiled.solid(feature).unwrap();
            let data = body.rust_internal_thread_treatment_data().unwrap();
            assert!(data.material.cavity().entry_radius_mm().is_some());
            assert_eq!(data.material.cavity().intent().terminal, None);
            assert_eq!(
                data.material.cavity().source_thread().exact_axial_extent(),
                data.material.cavity().exact_receiving_extent()
            );
            assert_eq!(
                data.material.cavity().source_thread().spec().left_hand,
                hand == "left"
            );
            assert_eq!(
                data.material.stock().height(),
                Rational::new(33_000_000_000_000_003, 10_000_000_000_000_000).unwrap()
            );
            body.manufacturing_mesh_with_receipt(MeshOptions {
                linear_deflection_mm: 0.02,
                angular_deflection_rad: 0.3,
                maximum_triangles: 500_000,
            })
            .unwrap();
            let catalog = run(&host, "thread_profiles", json!({}));
            let scope = &catalog["active_geometry_support"]["internal_entry_treatments"];
            assert_eq!(scope["manufacturing_mesh"], true);
            assert_eq!(scope["step"], false);
            let destination = dir.path().join("entry.step");
            assert!(
                dispatch(
                    &host,
                    &call(false),
                    "export_step",
                    json!({"part":created["part"],"file_name":"entry.step"})
                )
                .is_err()
            );
            assert!(!destination.exists());
            std::fs::write(&destination, b"preserve original manufacturing file").unwrap();
            assert!(
                dispatch(
                    &host,
                    &call(false),
                    "export_step",
                    json!({"part":created["part"],"file_name":"entry.step"})
                )
                .is_err()
            );
            assert_eq!(
                std::fs::read(destination).unwrap(),
                b"preserve original manufacturing file"
            );
            assert_eq!(host.document(), created_document);
            let edit = json!({"expected_revision":1,"thread":created["thread"],"end_treatments":ends(0.5)});
            assert_eq!(
                dispatch(&host, &call(true), "set_thread", edit.clone())
                    .unwrap_err()
                    .code,
                "cancelled"
            );
            assert_eq!(host.document(), created_document);
            let mut stale = edit.clone();
            stale["expected_revision"] = json!(0);
            assert_eq!(
                dispatch(&host, &call(false), "set_thread", stale)
                    .unwrap_err()
                    .code,
                "stale_revision"
            );
            assert_eq!(host.document(), created_document);
            run(&host, "set_thread", edit);
            let edited = host.document();
            run(&host, "undo", json!({"expected_revision":2}));
            assert_eq!(host.document(), created_document);
            run(&host, "redo", json!({"expected_revision":3}));
            assert_eq!(host.document(), edited);
            run(&host, "save_project", json!({"expected_revision":4}));
            drop(compiled);
            drop(host);
            let reopened = open();
            assert_eq!(reopened.document(), edited);
            let rebuilt = (Compiler {
                omit_display: true,
                ..Default::default()
            })
            .compile_with_backend(&reopened.document(), 0, GeometryBackend::Rust, None)
            .unwrap();
            let actual = rebuilt.solid(feature).unwrap();
            let cold = actual.rust_internal_thread_treatment_data().unwrap();
            assert!(cold.material.cavity().entry_radius_mm().is_some());
            assert_eq!(cold.material.cavity().intent().terminal, None);
            assert_eq!(cold.material.cavity().intent().entry_depth_mm, Some(0.5));
            assert_eq!(
                cold.material.cavity().source_thread().exact_axial_extent(),
                cold.material.cavity().exact_receiving_extent()
            );
            assert_eq!(cold.material.stock(), data.material.stock());
            let moved = actual
                .transformed(RigidTransform {
                    rotation: [[0.6, 0., 0.8], [0., 1., 0.], [-0.8, 0., 0.6]],
                    translation_mm: [11., -7., 5.],
                })
                .unwrap();
            for source in [actual, &moved] {
                source
                    .manufacturing_mesh_with_receipt(MeshOptions {
                        linear_deflection_mm: 0.02,
                        angular_deflection_rad: 0.3,
                        maximum_triangles: 500_000,
                    })
                    .unwrap();
            }
        }
    }
}
