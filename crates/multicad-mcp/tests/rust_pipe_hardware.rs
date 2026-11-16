//! Real pipe-hardware commands preserve source intent and transaction state.
use multicad_core::{
    commands::{Actor, Channel},
    document::{Document, EntityId},
    model,
    project_file::{self, ProjectFile},
};
use multicad_mcp::{
    Headless, Storage, Workspace,
    tools::{Call, dispatch},
};
use multicad_model::GeometryBackend;
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
fn mm(value: f64) -> Value {
    json!({"value":value,"unit":"mm"})
}
fn recipe(npt: bool, head: bool) -> Value {
    let (family, fit, field, intent, major, pitch, length) = if npt {
        let intent = model::NptThreadIntent {
            form: if head {
                model::NptThreadForm::External
            } else {
                model::NptThreadForm::Internal
            },
            size: model::NptThreadSize::OneEighth,
            gauge: model::NptGaugeConvention::Nbs1959PlainFace,
            truncation: model::NptTruncationTarget::MidpointConeParallel,
        };
        (
            "npt",
            "npt_historical_midpoint",
            "npt",
            json!(intent),
            intent.base_major_mm(),
            intent.preset().pitch_mm,
            intent.default_length_mm(),
        )
    } else {
        let intent = model::PipeThreadIntent {
            form: if head {
                model::PipeThreadForm::R
            } else {
                model::PipeThreadForm::Rc
            },
            size: model::PipeThreadSize::OneEighth,
            gauge: model::PipeGaugeConvention::Iso7Nominal1994,
        };
        (
            "bspt",
            "bspt_nominal",
            "pipe",
            json!(intent),
            intent.base_major_mm(),
            intent.preset().pitch_mm,
            intent.default_length_mm(),
        )
    };
    let mut thread = json!({
        "family":family,"fit":fit,"representation":"modeled",
        "major_diameter":mm(major),"pitch":mm(pitch),"length":mm(length),
        "plane":"yz","origin":[mm(11.),mm(-7.),mm(5.)]
    });
    thread[field] = intent;
    if head {
        json!({"kind":"machine_screw","underhead_length":mm(length),"shank_diameter":mm(major),
            "head":{"kind":"external_hex","across_flats":mm(18.),"height":mm(3.)},
            "drive":{"kind":"external_hex"},"thread":thread})
    } else {
        json!({"kind":"hex_nut","across_flats":mm(18.),"thread":thread})
    }
}

#[test]
fn original_pipe_hardware_commands_edit_undo_reopen_and_refuse_unsupported_export_atomically() {
    for npt in [true, false] {
        for head in [true, false] {
            let directory = tempfile::tempdir().unwrap();
            let path = directory.path().join("original-pipe-hardware.tcd");
            Store
                .save(
                    &path,
                    &ProjectFile::new(Document::new("Original pipe hardware").unwrap()).unwrap(),
                )
                .unwrap();
            let open = || {
                Headless::with_geometry_backend(
                    Workspace::project(path.clone(), Box::new(Store)).unwrap(),
                    Some(directory.path().into()),
                    false,
                    GeometryBackend::Rust,
                )
            };
            let host = open();
            let support = run(&host, "mechanical_part_profiles", json!({}));
            assert_eq!(
                support["active_geometry_support"]["pipe_hardware"]["manufacturing_mesh"],
                true
            );
            assert_eq!(
                support["active_geometry_support"]["pipe_hardware"]["step"],
                false
            );
            let pipe_support = &support["active_geometry_support"]["pipe_hardware"];
            assert_eq!(pipe_support["hand"], "right");
            assert_eq!(pipe_support["starts"], 1);
            assert_eq!(pipe_support["end_treatments"], false);
            let request =
                json!({"expected_revision":0,"name":"Original pipe law","recipe":recipe(npt,head)});
            let before = host.document();
            assert_eq!(
                dispatch(&host, &call(true), "set_mechanical_part", request.clone())
                    .unwrap_err()
                    .code,
                "cancelled"
            );
            assert_eq!(host.document(), before);
            let created = run(&host, "set_mechanical_part", request);
            let feature: EntityId = serde_json::from_value(created["feature"].clone()).unwrap();
            let original = model::read_mechanical_part(&host.document(), feature)
                .unwrap()
                .parameters;
            let mut edit = json!({"kind":if head {"machine_screw"} else {"hex_nut"}});
            if head {
                edit["head"] = json!({"kind":"external_hex","across_flats":mm(20.)});
            } else {
                edit["across_flats"] = mm(20.);
            }
            let changed = run(
                &host,
                "set_mechanical_part",
                json!({"expected_revision":1,"feature":feature,"recipe":edit}),
            );
            assert_eq!(changed["part"], created["part"]);
            let saved = host.document();
            let source = model::read_mechanical_part(&saved, feature)
                .unwrap()
                .parameters;
            let thread = match &source {
                model::MechanicalPartParameters::MachineScrew { thread, .. }
                | model::MechanicalPartParameters::HexNut { thread, .. } => thread,
                _ => unreachable!(),
            };
            assert_eq!(
                thread.family,
                if npt {
                    model::ThreadFamily::Npt
                } else {
                    model::ThreadFamily::Bspt
                }
            );
            let initial_thread = match &original {
                model::MechanicalPartParameters::MachineScrew { thread, .. }
                | model::MechanicalPartParameters::HexNut { thread, .. } => thread,
                _ => unreachable!(),
            };
            assert_eq!(
                thread, initial_thread,
                "an across-flats edit must preserve the original law and datum"
            );
            run(&host, "undo", json!({"expected_revision":2}));
            assert_eq!(
                model::read_mechanical_part(&host.document(), feature)
                    .unwrap()
                    .parameters,
                original
            );
            run(&host, "redo", json!({"expected_revision":3}));
            assert_eq!(host.document(), saved);
            let mut invalid = json!({"kind":if head {"machine_screw"} else {"hex_nut"}});
            if head {
                invalid["head"] = json!({"kind":"external_hex","across_flats":mm(1.)});
            } else {
                invalid["across_flats"] = mm(1.);
            }
            assert!(
                dispatch(
                    &host,
                    &call(false),
                    "set_mechanical_part",
                    json!({"expected_revision":4,"feature":feature,"recipe":invalid})
                )
                .is_err()
            );
            assert_eq!(host.document(), saved);
            let output = directory.path().join("unavailable.step");
            std::fs::write(&output, b"existing output").unwrap();
            assert!(
                dispatch(
                    &host,
                    &call(false),
                    "export_step",
                    json!({"part":created["part"],"file_name":"unavailable.step"})
                )
                .is_err()
            );
            assert_eq!(std::fs::read(&output).unwrap(), b"existing output");
            assert_eq!(host.document(), saved);
            run(&host, "save_project", json!({"expected_revision":4}));
            drop(host);
            let cold = open();
            assert_eq!(cold.document(), saved);
            assert_eq!(
                model::read_mechanical_part(&cold.document(), feature)
                    .unwrap()
                    .parameters,
                source
            );
            let model = run(&cold, "evaluate_model", json!({}));
            assert_eq!(model["parts"][0]["solid_count"], 1);
        }
    }
}
