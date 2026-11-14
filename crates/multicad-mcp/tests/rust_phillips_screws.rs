use multicad_core::{
    commands::{Actor, Channel},
    document::{Document, EntityId},
    model::{self, MechanicalPartParameters, ScrewDrive},
    project_file::{self, ProjectFile},
};
use multicad_kernel::MeshOptions;
use multicad_mcp::{
    Headless, Storage, Workspace,
    tools::{Call, dispatch},
};
use multicad_model::GeometryBackend;
use serde_json::{Value, json};
use std::sync::{Arc, atomic::AtomicBool};

fn call(cancel: bool) -> Call {
    Call {
        actor: Actor {
            id: "mcp:metric-phillips".into(),
            channel: Channel::Mcp,
        },
        cancel: Arc::new(AtomicBool::new(cancel)),
    }
}
fn run(host: &Headless, tool: &str, args: Value) -> Value {
    dispatch(host, &call(false), tool, args).unwrap()
}
fn mm(value: f64) -> Value {
    json!({"value":value,"unit":"mm"})
}
fn deg(value: f64) -> Value {
    json!({"value":value,"unit":"deg"})
}
struct Store;
impl Storage for Store {
    fn open(&mut self, path: &std::path::Path) -> Result<ProjectFile, String> {
        project_file::read(
            &mut std::fs::File::open(path).map_err(|e| e.to_string())?,
            Default::default(),
            &project_file::supported_capabilities(),
        )
        .map(|r| r.project)
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

#[test]
fn fitted_and_nominal_metric_phillips_mcp_are_transactional_and_reopen_actual_sources() {
    for h0_long in [false, true] {
        for fit in ["nominal_basic", "metric_6g_external"] {
            for hand in ["right", "left"] {
                let dir = tempfile::tempdir().unwrap();
                let path = dir.path().join("metric-phillips.tcd");
                Store
                    .save(
                        &path,
                        &ProjectFile::new(Document::new("Metric Phillips").unwrap()).unwrap(),
                    )
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
                let original = host.document();
                let profiles = run(&host, "mechanical_part_profiles", json!({}));
                assert!(
                    profiles["active_geometry_support"]["head_drive_pairs"]
                        .as_array()
                        .unwrap()
                        .contains(&json!(["pan", "phillips_h"]))
                );
                assert_eq!(
                    profiles["active_geometry_support"]["phillips_h"]["step"],
                    false
                );
                let mut args = json!({"expected_revision":0,"name":"M3 Phillips","recipe":{
                    "kind":"machine_screw","underhead_length":mm(9.),"shank_diameter":mm(3.),
                    "head":{"kind":"pan","diameter":mm(5.6),"height":mm(2.4),"crown_radius":mm(3.),"top_land_diameter":mm(4.2),"neck_radius":mm(0.1)},
                    "drive":{"kind":"phillips_h","profile":{"size":"h1","mouth_diameter":mm(3.),"mouth_radius":mm(0.5),"b_width":mm(0.955),"g_diameter":mm(1.27),"tilt_angle":deg(7.),"outer_flank_angle":deg(26.5),"floor_angle":deg(28.),"outer_corner_angle":deg(92.),"inner":{"kind":"shallow_v","angle":deg(138.),"e_width":mm(0.435)}}},
                    "thread":{"plane":"yz","origin":[mm(10.),mm(20.),mm(30.)],"family":"iso_metric","major_diameter":mm(3.),"pitch":mm(0.5),"length":mm(8.),"representation":"modeled","fit":fit,"hand":hand}
                }});
                if h0_long {
                    args["recipe"]["underhead_length"] = mm(6.);
                    args["recipe"]["shank_diameter"] = mm(2.);
                    args["recipe"]["head"] = json!({"kind":"pan","diameter":mm(4.),"height":mm(1.6),"crown_radius":mm(3.2),"top_land_diameter":mm(2.5),"neck_radius":mm(0.1)});
                    args["recipe"]["drive"]["profile"] = json!({"size":"h0","mouth_diameter":mm(1.9),"mouth_radius":mm(0.3),"b_width":mm(0.61),"g_diameter":mm(0.81),"tilt_angle":deg(7.),"outer_flank_angle":deg(26.5),"floor_angle":deg(28.),"outer_corner_angle":deg(92.),"inner":{"kind":"rounded","radius":mm(0.25)}});
                    args["recipe"]["thread"]["major_diameter"] = mm(2.);
                    args["recipe"]["thread"]["pitch"] = mm(0.4);
                    args["recipe"]["thread"]["length"] = mm(0.8);
                }
                if fit == "metric_6g_external" {
                    args["recipe"]["thread"]["metric_fit_engagement"] =
                        json!({"kind":"normal_group_assumed"});
                }
                // Sourced fit tables currently cover M3/M6. The H0 M2
                // recipe proves the Nominal path and keeps that limit honest.
                if h0_long && fit == "metric_6g_external" {
                    let error =
                        dispatch(&host, &call(false), "set_mechanical_part", args).unwrap_err();
                    assert_eq!(error.code, "invalid_document");
                    assert!(error.message.contains("M3×0.5 and M6×1 only"));
                    assert_eq!(host.document(), original);
                    continue;
                }
                let cancelled =
                    dispatch(&host, &call(true), "set_mechanical_part", args.clone()).unwrap_err();
                assert_eq!(
                    cancelled.code, "cancelled",
                    "h0={h0_long},fit={fit},hand={hand}: {cancelled:?}"
                );
                assert_eq!(host.document(), original);
                let created = run(&host, "set_mechanical_part", args);
                let feature: EntityId = serde_json::from_value(created["feature"].clone()).unwrap();
                assert!(matches!(
                    model::read_mechanical_part(&host.document(), feature)
                        .unwrap()
                        .parameters,
                    MechanicalPartParameters::MachineScrew {
                        drive: ScrewDrive::PhillipsH { .. },
                        ..
                    }
                ));
                let initial = host.document();
                let edited = run(
                    &host,
                    "set_mechanical_part",
                    json!({"expected_revision":1,"feature":feature,"recipe":{"kind":"machine_screw","underhead_length":mm(10.)}}),
                );
                assert_eq!(edited["part"], created["part"]);
                let document = host.document();
                run(&host, "undo", json!({"expected_revision":2}));
                assert_eq!(host.document(), initial);
                run(&host, "redo", json!({"expected_revision":3}));
                assert_eq!(host.document(), document);
                assert!(dispatch(&host,&call(false),"set_mechanical_part",json!({"expected_revision":4,"feature":feature,"recipe":{"kind":"machine_screw","drive":{"kind":"phillips_h","profile":{"mouth_diameter":mm(100.)}}}})).is_err());
                assert_eq!(host.document(), document);
                let export = dir.path().join("unfinished.step");
                assert!(
                    dispatch(
                        &host,
                        &call(false),
                        "export_step",
                        json!({"part":created["part"],"file_name":"unfinished.step"})
                    )
                    .is_err()
                );
                assert!(!export.exists());
                std::fs::write(&export, b"accepted output").unwrap();
                assert!(
                    dispatch(
                        &host,
                        &call(false),
                        "export_step",
                        json!({"part":created["part"],"file_name":"unfinished.step"})
                    )
                    .is_err()
                );
                assert_eq!(std::fs::read(export).unwrap(), b"accepted output");
                run(&host, "save_project", json!({"expected_revision":4}));
                drop(host);
                let reopened = open();
                assert_eq!(reopened.document(), document);
                let compiled = multicad_model::Compiler {
                    omit_display: true,
                    ..Default::default()
                }
                .compile_with_backend(&reopened.document(), 0, GeometryBackend::Rust, None)
                .unwrap();
                let solid = compiled.solid(feature).unwrap();
                let (source, frame) = solid.rust_phillips_screw_geometry().unwrap();
                assert_eq!(
                    source.cavity().source().spec().mouth_diameter_mm,
                    if h0_long { 1.9 } else { 3. }
                );
                assert_eq!(source.thread().spec().left_hand, hand == "left");
                assert_eq!(
                    source.thread().spec().length_mm,
                    if h0_long { 0.8 } else { 8. }
                );
                assert!((frame.normal()[0].abs() - 1.).abs() < 1e-12);
                assert!(
                    !solid
                        .manufacturing_mesh(MeshOptions::default())
                        .unwrap()
                        .indices
                        .is_empty()
                );
            }
        }
    }
}
