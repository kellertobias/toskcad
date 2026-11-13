//! The normal saved-project command path must retain a metric screw's plain zone.
use multicad_core::{
    document::{Document, EntityId},
    project_file::{self, ProjectFile},
};
use multicad_kernel::{GeometryPolicy, Solid};
use serde_json::{Value, json};
use std::{
    io::Write,
    path::Path,
    process::{Command, Stdio},
};
fn invoke(project: &Path, command: &Path, tool: &str, arguments: Value, save: bool) -> Value {
    std::fs::write(
        command,
        serde_json::to_vec(&json!({"tool":tool,"arguments":arguments})).unwrap(),
    )
    .unwrap();
    let mut child = Command::new(env!("CARGO_BIN_EXE_multicad"));
    child.args([
        "--command",
        command.to_str().unwrap(),
        "--project",
        project.to_str().unwrap(),
    ]);
    if save {
        child.arg("--save");
    }
    let output = child.output().unwrap();
    let result: Value = serde_json::from_slice(&output.stdout).unwrap_or_else(|error| {
        panic!(
            "{tool}: {error}; status={}; stdout={}; stderr={}",
            output.status,
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        )
    });
    assert!(
        output.status.success(),
        "{result}; {}",
        String::from_utf8_lossy(&output.stderr)
    );
    result["result"].clone()
}
fn mm(value: f64) -> Value {
    json!({"value":value,"unit":"mm"})
}
fn export_via_stdio(project: &Path, exports: &Path, part: &Value) -> Value {
    let mut child = Command::new(env!("CARGO_BIN_EXE_multicad"))
        .args([
            "--mcp",
            "--project",
            project.to_str().unwrap(),
            "--export-dir",
            exports.to_str().unwrap(),
        ])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let requests = [
        json!({"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2025-06-18","capabilities":{},"clientInfo":{"name":"hex-shank-test","version":"1"}}}),
        json!({"jsonrpc":"2.0","method":"notifications/initialized"}),
        json!({"jsonrpc":"2.0","id":2,"method":"tools/call","params":{"name":"export_step","arguments":{"part":part,"file_name":"actual-plain-shank.step"}}}),
    ];
    {
        let mut input = child.stdin.take().unwrap();
        for request in requests {
            writeln!(input, "{request}").unwrap();
        }
    }
    let output = child.wait_with_output().unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let response = String::from_utf8(output.stdout)
        .unwrap()
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).unwrap())
        .find(|reply| reply["id"] == 2)
        .unwrap();
    assert_eq!(response["result"]["isError"], false, "{response}");
    serde_json::from_str(response["result"]["content"][0]["text"].as_str().unwrap()).unwrap()
}
#[test]
fn saved_metric_hex_head_plain_shank_exports_and_blindly_reimports_through_real_commands() {
    let dir = tempfile::tempdir().unwrap();
    let project = dir.path().join("metric-plain-shank.tcd");
    let command = dir.path().join("command.json");
    let exports = dir.path().join("exports");
    std::fs::create_dir(&exports).unwrap();
    let mut bytes = Vec::new();
    project_file::write(
        &ProjectFile::new(Document::new("Metric plain shank").unwrap()).unwrap(),
        &mut bytes,
        Default::default(),
    )
    .unwrap();
    std::fs::write(&project, bytes).unwrap();
    let source = invoke(
        &project,
        &command,
        "set_mechanical_part",
        json!({
            "expected_revision":0,"name":"M6 hex screw with plain shank","recipe":{
                "kind":"machine_screw","underhead_length":mm(8.),"shank_diameter":mm(6.),
                "head":{"kind":"external_hex","across_flats":mm(10.),"height":mm(4.)},
                "drive":{"kind":"external_hex"},
                "thread":{"family":"iso_metric","fit":"nominal_basic","major_diameter":mm(6.),"pitch":mm(1.),"length":mm(3.7),"representation":"modeled"}
            }
        }),
        true,
    );
    assert_eq!(source["revision"], 1);
    let saved = std::fs::read(&project).unwrap();
    let inspected = invoke(
        &project,
        &command,
        "inspect_mechanical_part",
        json!({"feature":source["feature"]}),
        false,
    );
    assert_eq!(inspected["resolved"]["fully_threaded"], false);
    assert!((inspected["resolved"]["thread_start_mm"].as_f64().unwrap() - 4.3).abs() < 1e-12);
    let receipt = export_via_stdio(&project, &exports, &source["part"]);
    let text = std::fs::read_to_string(exports.join("actual-plain-shank.step")).unwrap();
    let (physical, imported) =
        multicad_kernel::rust_kernel::step_read::read_import_solid(&text, 1e-7).unwrap();
    assert_eq!(physical.solid_count, 1);
    assert_eq!(physical.face_count, 15);
    assert!((physical.volume_mm3 - receipt["readback_volume_mm3"].as_f64().unwrap()).abs() < 1e-6);
    let multicad_kernel::rust_kernel::step_read::ImportProfile::HelicalScrew { screw, .. } =
        imported.profile
    else {
        panic!("physical head and shank required")
    };
    assert!((screw.shank_length_mm - 4.3).abs() < 1e-9);
    assert_eq!(screw.shank_radius_mm, 3.);
    let solid = Solid::import_rust_step(&text, EntityId::new(), GeometryPolicy::default()).unwrap();
    let mesh = solid
        .manufacturing_mesh_with_receipt(Default::default())
        .unwrap();
    assert_eq!(mesh.edge_pairs as usize * 2, mesh.mesh.indices.len());
    assert_eq!(std::fs::read(&project).unwrap(), saved);
}
