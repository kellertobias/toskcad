use multicad_core::{
    commands::{Actor, Channel},
    document::{Document, EntityId},
    model::{self, MechanicalPartParameters, ScrewDrive},
    project_file::{self, ProjectFile},
};
use multicad_kernel::{GeometryPolicy, Solid};
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
            id: "mcp:metric-torx".into(),
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
fn metric_curved_drive_mcp_preserves_real_geometry_history_and_reopen() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("metric-torx.tcd");
    Store
        .save(
            &path,
            &ProjectFile::new(Document::new("Metric Torx").unwrap()).unwrap(),
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
            .contains(&json!(["socket_cap", "hexalobular"]))
    );
    let args = json!({"expected_revision":0,"name":"M3 curved drive","recipe":{
        "kind":"machine_screw",
        "underhead_length":mm(3.),"shank_diameter":mm(3.),
        "head":{"kind":"socket_cap","diameter":mm(6.),"height":mm(2.5)},
        "drive":{"kind":"hexalobular","major_diameter":mm(2.),"minor_diameter":mm(1.5),"outer_lobe_radius":mm(0.2),"depth":mm(0.25),"bottom_angle":deg(140.)},
        "thread":{"plane":"yz","origin":[mm(10.),mm(20.),mm(30.)],"family":"iso_metric","major_diameter":mm(3.),"pitch":mm(0.5),"length":mm(2.),"representation":"modeled","fit":"metric_6g_external"}
    }});
    assert_eq!(
        dispatch(&host, &call(true), "set_mechanical_part", args.clone())
            .unwrap_err()
            .code,
        "cancelled"
    );
    assert_eq!(host.document(), original);
    let created = run(&host, "set_mechanical_part", args);
    let feature: EntityId = serde_json::from_value(created["feature"].clone()).unwrap();
    let definition = model::read_mechanical_part(&host.document(), feature).unwrap();
    assert!(matches!(
        definition.parameters,
        MechanicalPartParameters::MachineScrew {
            drive: ScrewDrive::Hexalobular { .. },
            ..
        }
    ));
    let initial = host.document();
    let edited = run(
        &host,
        "set_mechanical_part",
        json!({"expected_revision":1,"feature":feature,"recipe":{
            "kind":"machine_screw","underhead_length":mm(4.),
            "head":{"kind":"countersunk","diameter":mm(6.),"height":mm(1.65),"included_angle":deg(90.),"edge_chamfer":mm(0.1),"neck_radius":mm(0.1)}
        }}),
    );
    assert_eq!(edited["part"], created["part"]);
    let document = host.document();
    run(&host, "undo", json!({"expected_revision":2}));
    assert_eq!(host.document(), initial);
    run(&host, "redo", json!({"expected_revision":3}));
    assert_eq!(host.document(), document);
    let before_failure = host.document();
    assert!(dispatch(&host,&call(false),"set_mechanical_part",json!({"expected_revision":4,"feature":feature,"recipe":{"kind":"machine_screw","drive":{"kind":"hexalobular","depth":mm(100.)}}})).is_err());
    assert_eq!(host.document(), before_failure);
    run(
        &host,
        "export_step",
        json!({"part":created["part"],"file_name":"metric-torx.step"}),
    );
    let physical = Solid::import_rust_step(
        &std::fs::read_to_string(dir.path().join("metric-torx.step")).unwrap(),
        EntityId::new(),
        GeometryPolicy::default(),
    )
    .unwrap();
    assert!(physical.uses_rust_geometry());
    let (source, frame) = physical.rust_profiled_screw_geometry().unwrap();
    assert!(source.spec().head.is_none());
    assert!(source.recess().hexalobular_geometry().is_some());
    assert!(source.recess().flat_hex_geometry().is_none());
    assert_eq!(source.thread().spec().pitch_mm, 0.5);
    assert!((frame.normal()[0].abs() - 1.).abs() < 1e-12);
    run(&host, "save_project", json!({"expected_revision":4}));
    drop(host);
    let reopened = open();
    assert_eq!(reopened.document(), document);
    run(
        &reopened,
        "export_step",
        json!({"part":created["part"],"file_name":"reopened-torx.step"}),
    );
    let second = Solid::import_rust_step(
        &std::fs::read_to_string(dir.path().join("reopened-torx.step")).unwrap(),
        EntityId::new(),
        GeometryPolicy::default(),
    )
    .unwrap();
    assert!((second.volume_mm3().unwrap() - physical.volume_mm3().unwrap()).abs() < 1e-6);
}

#[test]
fn reported_metric_nut_step_support_has_actual_independent_rust_geometry() {
    let dir = tempfile::tempdir().unwrap();
    let host = Headless::with_geometry_backend(
        Workspace::new(Document::new("Metric nut exchange").unwrap()).unwrap(),
        Some(dir.path().into()),
        false,
        GeometryBackend::Rust,
    );
    let created = run(
        &host,
        "set_mechanical_part",
        json!({"expected_revision":0,"name":"M6 nut","recipe":{
            "kind":"hex_nut","across_flats":mm(10.),"thread":{"major_diameter":mm(6.),"pitch":mm(1.),"length":mm(5.2),"fit":"metric_6h_internal","representation":"modeled","minimum_wall":mm(1.)}
        }}),
    );
    run(
        &host,
        "export_step",
        json!({"part":created["part"],"file_name":"metric-nut.step"}),
    );
    let source = Solid::import_rust_step(
        &std::fs::read_to_string(dir.path().join("metric-nut.step")).unwrap(),
        EntityId::new(),
        GeometryPolicy::default(),
    )
    .unwrap();
    assert!(source.uses_rust_geometry());
    assert!(source.volume_mm3().unwrap() > 0.);
    let profiles = run(&host, "mechanical_part_profiles", json!({}));
    assert!(
        profiles["active_geometry_support"]["step_geometry"]
            .as_array()
            .unwrap()
            .contains(&json!("hex_nut"))
    );
}

#[test]
fn advertised_cosmetic_curved_drive_heads_create_and_reopen_without_manufacturing_claims() {
    for head in [
        json!({"kind":"socket_cap","diameter":mm(5.6),"height":mm(2.4)}),
        json!({"kind":"pan","diameter":mm(5.6),"height":mm(2.4),"crown_radius":mm(5.),"top_land_diameter":mm(4.1),"neck_radius":mm(0.125)}),
        json!({"kind":"button","diameter":mm(5.6),"height":mm(2.4),"crown_radius":mm(5.),"top_land_diameter":mm(4.1),"neck_radius":mm(0.125)}),
        json!({"kind":"countersunk","diameter":mm(5.6),"height":mm(1.8),"included_angle":deg(90.),"edge_chamfer":mm(0.1),"neck_radius":mm(0.)}),
    ] {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("cosmetic-head.tcd");
        Store
            .save(
                &path,
                &ProjectFile::new(Document::new("Cosmetic curved drive").unwrap()).unwrap(),
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
        let profiles = run(&host, "mechanical_part_profiles", json!({}));
        assert!(
            profiles["active_geometry_support"]["cosmetic_head_drive_pairs"]
                .as_array()
                .unwrap()
                .contains(&json!([head["kind"], "hexalobular"]))
        );
        let created = run(
            &host,
            "set_mechanical_part",
            json!({"expected_revision":0,"name":"Cosmetic curved drive","recipe":{
                "kind":"machine_screw","underhead_length":mm(8.),"shank_diameter":mm(3.),"head":head,
                "drive":{"kind":"hexalobular","major_diameter":mm(2.8),"minor_diameter":mm(2.05),"outer_lobe_radius":mm(0.28),"depth":mm(0.5),"bottom_angle":deg(120.)},
                "thread":{"plane":"yz","origin":[mm(10.),mm(20.),mm(30.)],"major_diameter":mm(3.),"pitch":mm(0.5),"length":mm(5.),"representation":"cosmetic_envelope","fit":"nominal_basic"}
            }}),
        );
        let initial = host.document();
        let evaluated = run(&host, "evaluate_model", json!({}));
        assert_eq!(evaluated["parts"].as_array().unwrap().len(), 1);
        let source_hash: [u8; 32] =
            serde_json::from_value(evaluated["source_hash"].clone()).unwrap();
        let error = dispatch(
            &host,
            &call(false),
            "export_step",
            json!({"part":created["part"],"file_name":"before.step"}),
        )
        .unwrap_err();
        assert_eq!(error.code, "unsupported");
        assert!(!dir.path().join("before.step").exists());
        assert_eq!(host.document(), initial);
        run(&host, "save_project", json!({"expected_revision":1}));
        drop(host);
        let reopened = open();
        assert_eq!(reopened.document(), initial);
        assert_eq!(
            run(&reopened, "evaluate_model", json!({}))["source_hash"],
            json!(source_hash)
        );
        assert!(
            dispatch(
                &reopened,
                &call(false),
                "export_step",
                json!({"part":created["part"],"file_name":"reopened.step"})
            )
            .is_err()
        );
        assert!(!dir.path().join("reopened.step").exists());
        assert_eq!(reopened.document(), initial);
    }
}
