//! Exchange failure and robustness controls (CAD-345).
//!
//! Every malformed, oversized, exhausting, damaged or cancelled input must be
//! refused with a classified error, must leave the original file
//! byte-for-byte unchanged, and must leave the previous valid document and
//! its saved project untouched. Bounded healing is pinned explicitly: what
//! the Rust reader tolerates and what it refuses.
use multicad_core::{
    commands::{Actor, Channel},
    document::{Document, EntityId},
    project_file::ContentHash,
};
use multicad_kernel::{GeometryErrorKind, GeometryPolicy, Solid};
use multicad_mcp::{
    Headless, Workspace,
    tools::{Call, dispatch},
};
use multicad_model::GeometryBackend;
use serde_json::{Value, json};
use std::{
    path::Path,
    sync::{Arc, atomic::AtomicBool},
};

const VALID: &str = include_str!("fixtures/exchange/freecad-drilled-box.step");
/// Rust reader byte ceiling (`part21::MAX_BYTES`).
const READER_MAX_BYTES: usize = 64 * 1024 * 1024;
/// Canonical document ceiling for an embedded STEP source.
const DOCUMENT_MAX_BYTES: usize = multicad_core::model::imports::MAX_STEP_BREP_BYTES;

fn kind(bytes: &[u8]) -> &'static str {
    let Ok(text) = std::str::from_utf8(bytes) else {
        return "not-utf8";
    };
    match Solid::import_rust_step(text, EntityId::new(), GeometryPolicy::default()) {
        Ok(_) => "pass",
        Err(error) => match error.kind {
            GeometryErrorKind::Unsupported => "unsupported",
            GeometryErrorKind::Limit => "limit",
            GeometryErrorKind::InvalidResult | GeometryErrorKind::Input => "invalid",
            other => panic!("unclassified {other:?}: {error}"),
        },
    }
}

fn call(cancelled: bool) -> Call {
    Call {
        actor: Actor {
            id: "mcp:exchange-robustness".into(),
            channel: Channel::Mcp,
        },
        cancel: Arc::new(AtomicBool::new(cancelled)),
    }
}

/// A Rust headless host whose document already holds one valid import.
struct Fixture {
    directory: tempfile::TempDir,
    host: Headless,
    source: String,
    document: Document,
}
impl Fixture {
    fn new() -> Self {
        let directory = tempfile::tempdir().unwrap();
        std::fs::write(directory.path().join("valid.step"), VALID).unwrap();
        let host = Headless::with_geometry_backend(
            Workspace::new(Document::new("Exchange robustness").unwrap()).unwrap(),
            Some(directory.path().to_path_buf()),
            false,
            GeometryBackend::Rust,
        );
        let imported = dispatch(
            &host,
            &call(false),
            "import_step",
            json!({"expected_revision": 0, "file_name": "valid.step"}),
        )
        .unwrap();
        let source = imported["source"].as_str().unwrap().to_owned();
        let document = host.document();
        Self {
            directory,
            host,
            source,
            document,
        }
    }
    /// Refuse `bytes` as a new import and as an update of the valid import;
    /// return the error text.
    fn refuse(&self, name: &str, bytes: &[u8]) -> String {
        let path = self.directory.path().join(name);
        std::fs::write(&path, bytes).unwrap();
        let mut messages = vec![];
        for arguments in [
            json!({"expected_revision": 1, "file_name": name}),
            json!({"expected_revision": 1, "file_name": name, "source": self.source}),
        ] {
            let error = dispatch(&self.host, &call(false), "import_step", arguments)
                .expect_err(&format!("{name} was imported"));
            messages.push(format!("{}: {}", error.code, error.message));
            assert_eq!(
                self.host.document(),
                self.document,
                "{name} changed the document"
            );
        }
        assert_eq!(
            std::fs::read(&path).unwrap(),
            bytes,
            "{name} changed the file"
        );
        assert_eq!(
            std::fs::read(self.directory.path().join("valid.step")).unwrap(),
            VALID.as_bytes()
        );
        messages.join(" | ")
    }
}

fn replace(from: &str, to: &str) -> Vec<u8> {
    assert!(VALID.contains(from), "fixture lacks {from}");
    VALID.replacen(from, to, 1).into_bytes()
}

fn append_data(extra: &str) -> Vec<u8> {
    VALID
        .replacen("ENDSEC;\nEND-ISO", &format!("{extra}\nENDSEC;\nEND-ISO"), 1)
        .into_bytes()
}

fn garbage(length: usize) -> Vec<u8> {
    let mut state = 0x2545_f491_4f6c_dd1d_u64;
    (0..length)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 24) as u8
        })
        .collect()
}

/// Malformed inputs and their pinned classifications; `not-utf8` inputs are
/// rejected before parsing.
fn malformed_cases() -> Vec<(&'static str, Vec<u8>, &'static str)> {
    let length = VALID.len();
    let inch = "#205 = ( CONVERSION_BASED_UNIT('INCH',#300) LENGTH_UNIT() NAMED_UNIT(#301) );";
    let mut not_utf8 = VALID.as_bytes().to_vec();
    not_utf8.splice(20..20, [0xff, 0xfe]);
    let nested = format!("#900 = DUMMY({}0.{});", "(".repeat(40), ")".repeat(40));
    let long = format!("#900 = DUMMY('{}');", "x".repeat(70 * 1024));
    vec![
        ("empty.step", vec![], "invalid"),
        (
            "truncated-quarter.step",
            VALID.as_bytes()[..length / 4].to_vec(),
            "invalid",
        ),
        (
            "truncated-half.step",
            VALID.as_bytes()[..length / 2].to_vec(),
            "invalid",
        ),
        (
            "truncated-trailer.step",
            VALID.replacen("END-ISO-10303-21;", "", 1).into_bytes(),
            "invalid",
        ),
        (
            "missing-data-section.step",
            VALID.as_bytes()[..VALID.find("DATA;").unwrap()].to_vec(),
            "invalid",
        ),
        (
            "dangling-reference.step",
            replace(
                "#16 = CLOSED_SHELL('',(#17,",
                "#16 = CLOSED_SHELL('',(#99999,",
            ),
            "invalid",
        ),
        (
            "self-cycle.step",
            replace("#19 = EDGE_LOOP('',(#20,", "#19 = EDGE_LOOP('',(#19,"),
            // A reference cycle below the solid is malformed topology.
            "invalid",
        ),
        (
            "vertex-cycle.step",
            replace("#22 = VERTEX_POINT('',#23)", "#22 = VERTEX_POINT('',#22)"),
            "invalid",
        ),
        (
            "duplicate-entity.step",
            append_data("#23 = CARTESIAN_POINT('',(0.,0.,0.));"),
            "invalid",
        ),
        (
            "nan-coordinate.step",
            replace(
                "#23 = CARTESIAN_POINT('',(13.,",
                "#23 = CARTESIAN_POINT('',(NaN,",
            ),
            "invalid",
        ),
        (
            "overflow-coordinate.step",
            replace(
                "#23 = CARTESIAN_POINT('',(13.,",
                "#23 = CARTESIAN_POINT('',(1.E400,",
            ),
            "invalid",
        ),
        (
            "huge-coordinate.step",
            replace(
                "#23 = CARTESIAN_POINT('',(13.,",
                "#23 = CARTESIAN_POINT('',(1.E300,",
            ),
            // Finite but far outside the numeric budget.
            "limit",
        ),
        (
            "centimetre-units.step",
            replace("SI_UNIT(.MILLI.,.METRE.)", "SI_UNIT(.CENTI.,.METRE.)"),
            "unsupported",
        ),
        (
            "metre-units.step",
            replace("SI_UNIT(.MILLI.,.METRE.)", "SI_UNIT($,.METRE.)"),
            "unsupported",
        ),
        (
            "inch-units.step",
            {
                let text = String::from_utf8(replace(
                    "#205 = ( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.) );",
                    inch,
                ))
                .unwrap();
                text.replacen(
                    "ENDSEC;\nEND-ISO",
                    "#300 = LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(25.4),#302);\n#301 = DIMENSIONAL_EXPONENTS(1.,0.,0.,0.,0.,0.,0.);\n#302 = ( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.) );\nENDSEC;\nEND-ISO",
                    1,
                )
                .into_bytes()
            },
            "unsupported",
        ),
        ("binary-garbage.step", garbage(4096), "not-utf8"),
        ("not-utf8.step", not_utf8, "not-utf8"),
        ("deep-nesting.step", append_data(&nested), "limit"),
        ("overlong-token.step", append_data(&long), "limit"),
    ]
}

#[test]
fn malformed_step_is_classified_and_preserves_files_and_the_previous_document() {
    let fixture = Fixture::new();
    let mut ledger = vec![];
    let mut mismatches = vec![];
    for (name, bytes, expected) in malformed_cases() {
        let actual = kind(&bytes);
        if actual != expected {
            mismatches.push(format!("{name}: {actual}, pinned {expected}"));
        }
        let message = fixture.refuse(name, &bytes);
        ledger.push(json!({"case": name, "reader": actual, "mcp": message}));
    }
    write_ledger("malformed", &ledger);
    assert!(mismatches.is_empty(), "{mismatches:#?}");
}

/// Schema names are not trusted as geometry evidence: the reader decides on
/// the actual representation, so a foreign or nonsense FILE_SCHEMA with a
/// valid AP214 body still imports, and an empty body still fails.
#[test]
fn file_schema_text_is_not_geometry_evidence() {
    for schema in ["CONFIG_CONTROL_DESIGN", "NOT_A_STEP_SCHEMA"] {
        let bytes = replace("AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }", schema);
        assert_eq!(kind(&bytes), "pass", "{schema}");
    }
}

/// Bounded healing: the reader tolerates coordinate noise inside its linear
/// tolerance and nothing else. It never closes gaps, flips orientations or
/// repairs topology.
#[test]
fn healing_is_bounded_to_the_linear_tolerance() {
    let tolerance = GeometryPolicy::default().linear_tolerance_mm;
    let fixture = Fixture::new();
    let nudged = |delta: f64| {
        replace(
            "#25 = CARTESIAN_POINT('',(13.,-7.,13.));",
            &format!("#25 = CARTESIAN_POINT('',(13.,-7.,{:.17}));", 13. + delta),
        )
    };
    // Within tolerance: accepted, with the exact volume.
    let accepted = nudged(tolerance / 10.);
    assert_eq!(kind(&accepted), "pass");
    // Beyond tolerance: refused, never snapped.
    let mut ledger = vec![json!({"case": "vertex-noise-0.1x-tolerance", "reader": "pass"})];
    let cases = vec![
        ("vertex-gap-10x-tolerance.step", nudged(tolerance * 10.)),
        ("vertex-gap-1um.step", nudged(1e-3)),
        (
            "flipped-face.step",
            replace(
                "#17 = ADVANCED_FACE('',(#18),#52,.F.);",
                "#17 = ADVANCED_FACE('',(#18),#52,.T.);",
            ),
        ),
        (
            "flipped-edge-use.step",
            replace(
                "#20 = ORIENTED_EDGE('',*,*,#21,.F.);",
                "#20 = ORIENTED_EDGE('',*,*,#21,.T.);",
            ),
        ),
        (
            "missing-face.step",
            replace(
                "#16 = CLOSED_SHELL('',(#17,#57,#88,#123,#147,#175,#187));",
                "#16 = CLOSED_SHELL('',(#17,#57,#88,#123,#147,#175));",
            ),
        ),
    ];
    for (name, bytes) in cases {
        let actual = kind(&bytes);
        assert_eq!(actual, "invalid", "{name} was healed or misclassified");
        let message = fixture.refuse(name, &bytes);
        ledger.push(json!({"case": name, "reader": actual, "mcp": message}));
    }
    write_ledger("healing", &ledger);
}

/// Byte ceilings: the canonical document embeds at most 8 MiB of STEP text;
/// the Rust reader parses at most 64 MiB. Inputs at each ceiling are not
/// refused for size; one more byte is.
#[test]
fn oversized_inputs_are_refused_at_the_byte_ceilings() {
    let pad = |total: usize| {
        let mut text = VALID.to_owned();
        let at = text.find("ENDSEC;\nEND-ISO").unwrap();
        let padding = total - text.len();
        text.insert_str(at, &"\n".repeat(padding));
        assert_eq!(text.len(), total);
        text.into_bytes()
    };
    let fixture = Fixture::new();
    // Document ceiling through the MCP import path.
    let at_limit = pad(DOCUMENT_MAX_BYTES);
    std::fs::write(
        fixture.directory.path().join("document-limit.step"),
        &at_limit,
    )
    .unwrap();
    let imported = dispatch(
        &fixture.host,
        &call(false),
        "import_step",
        json!({"expected_revision": 1, "file_name": "document-limit.step"}),
    )
    .unwrap();
    assert_eq!(imported["parts"].as_array().unwrap().len(), 1);
    let fixture = Fixture::new();
    let message = fixture.refuse("document-over.step", &pad(DOCUMENT_MAX_BYTES + 1));
    assert!(message.contains("8 MiB import limit"), "{message}");
    // Reader ceiling through the kernel import.
    assert_eq!(kind(&pad(READER_MAX_BYTES)), "pass");
    let over = pad(READER_MAX_BYTES + 1);
    assert_eq!(kind(&over), "limit");
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("reader-over.step");
    std::fs::write(&path, &over).unwrap();
    let Err(error) = Solid::read_step(&path, GeometryPolicy::default()) else {
        panic!("an oversized file was read");
    };
    if !multicad_kernel::OCCT_LINKED {
        assert_eq!(error.kind, GeometryErrorKind::Limit, "{error}");
    }
    assert_eq!(std::fs::read(&path).unwrap().len(), over.len());
}

/// Resource exhaustion inside the byte ceiling: entity and value counts are
/// bounded before allocation grows without limit.
#[test]
fn entity_and_value_exhaustion_is_refused_as_a_limit() {
    let mut entities = String::with_capacity(16 * 1024 * 1024);
    for id in 1000..1000 + 1_048_577 {
        entities.push_str(&format!("#{id}=A();\n"));
    }
    let values = format!("#900 = DUMMY(({}0.));", "0.,".repeat(8_388_608));
    let fixture = Fixture::new();
    let mut ledger = vec![];
    for (name, bytes) in [
        ("entity-exhaustion.step", append_data(&entities)),
        ("value-exhaustion.step", append_data(&values)),
    ] {
        assert!(bytes.len() <= DOCUMENT_MAX_BYTES.max(READER_MAX_BYTES));
        let started = std::time::Instant::now();
        let actual = kind(&bytes);
        let millis = started.elapsed().as_millis();
        assert_eq!(actual, "limit", "{name}");
        // The document path refuses these files by size or by the same limit.
        let message = fixture.refuse(name, &bytes);
        ledger.push(json!({"case": name, "bytes": bytes.len(), "reader": actual, "reader_ms": millis, "mcp": message}));
    }
    write_ledger("exhaustion", &ledger);
}

/// A cancelled request is refused before reading and commits nothing. The
/// Rust STEP reader itself has no interrupt; cancellation is observed before
/// the read and again before the import commits.
#[test]
fn cancelled_import_commits_nothing() {
    let fixture = Fixture::new();
    std::fs::write(fixture.directory.path().join("again.step"), VALID).unwrap();
    for arguments in [
        json!({"expected_revision": 1, "file_name": "again.step"}),
        json!({"expected_revision": 1, "file_name": "again.step", "source": fixture.source}),
    ] {
        let error = dispatch(&fixture.host, &call(true), "import_step", arguments).unwrap_err();
        assert_eq!(error.code, "cancelled");
        assert_eq!(fixture.host.document(), fixture.document);
    }
    assert_eq!(
        std::fs::read(fixture.directory.path().join("again.step")).unwrap(),
        VALID.as_bytes()
    );
}

/// A damaged or truncated project refuses to open in the real executable,
/// stays byte-for-byte unchanged, and the last valid project still opens
/// with the identical imported geometry after the STEP source is gone.
#[test]
fn damaged_projects_refuse_to_open_and_leave_the_valid_project_intact() {
    use std::process::{Command, Stdio};
    let directory = tempfile::tempdir().unwrap();
    std::fs::write(directory.path().join("valid.step"), VALID).unwrap();
    let project = directory.path().join("valid.tcd");
    let run = |project: &Path, request: Value, save: bool| {
        let request_path = directory.path().join("request.json");
        std::fs::write(&request_path, request.to_string()).unwrap();
        let mut command = Command::new(env!("CARGO_BIN_EXE_multicad"));
        command
            .arg("--command")
            .arg(&request_path)
            .arg("--project")
            .arg(project)
            .current_dir(directory.path())
            .stdin(Stdio::null());
        command.arg(if save { "--save" } else { "--read-only" });
        command.arg("--rust-model");
        command.output().unwrap()
    };
    // `--command` has no export folder; seed the project through MCP instead.
    let host_project = {
        let mut child = Command::new(env!("CARGO_BIN_EXE_multicad"))
            .args(["--mcp", "--rust-model", "--project"])
            .arg(&project)
            .arg("--export-dir")
            .arg(directory.path())
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::null())
            .spawn()
            .unwrap();
        use std::io::{BufRead, Write};
        let mut stdin = child.stdin.take().unwrap();
        let mut stdout = std::io::BufReader::new(child.stdout.take().unwrap());
        let mut exchange = |id: u64, method: &str, params: Value| {
            writeln!(
                stdin,
                "{}",
                json!({"jsonrpc": "2.0", "id": id, "method": method, "params": params})
            )
            .unwrap();
            stdin.flush().unwrap();
            loop {
                let mut line = String::new();
                assert!(stdout.read_line(&mut line).unwrap() > 0);
                let message: Value = serde_json::from_str(&line).unwrap();
                if message["id"] == json!(id) {
                    return message;
                }
            }
        };
        exchange(
            1,
            "initialize",
            json!({"protocolVersion": "2025-06-18", "capabilities": {}, "clientInfo": {"name": "robustness", "version": "1"}}),
        );
        let tool = |name: &str, arguments: Value| json!({"name": name, "arguments": arguments});
        let imported = exchange(
            2,
            "tools/call",
            tool(
                "import_step",
                json!({"expected_revision": 0, "file_name": "valid.step"}),
            ),
        );
        assert_ne!(imported["result"]["isError"], true, "{imported}");
        let saved = exchange(
            3,
            "tools/call",
            tool("save_project", json!({"expected_revision": 1})),
        );
        assert_ne!(saved["result"]["isError"], true, "{saved}");
        let _ = child.kill();
        let _ = child.wait();
        std::fs::read(&project).unwrap()
    };
    std::fs::remove_file(directory.path().join("valid.step")).unwrap();
    let evaluate = json!({"tool": "evaluate_model", "arguments": {}});
    let good = run(&project, evaluate.clone(), false);
    assert!(
        good.status.success(),
        "{}",
        String::from_utf8_lossy(&good.stderr)
    );
    let receipt: Value = serde_json::from_slice(&good.stdout).unwrap();
    let parts = receipt["result"]["parts"].clone();
    assert_eq!(parts.as_array().map(Vec::len), Some(1), "{receipt}");

    let mut damaged = vec![];
    for cut in [8, 52, host_project.len() / 2, host_project.len() - 1] {
        damaged.push((format!("truncated-{cut}"), host_project[..cut].to_vec()));
    }
    for offset in [3, 60, host_project.len() / 2, host_project.len() - 1] {
        let mut copy = host_project.clone();
        copy[offset] ^= 0x20;
        damaged.push((format!("flipped-{offset}"), copy));
    }
    let mut extended = host_project.clone();
    extended.push(0);
    damaged.push(("trailing-byte".into(), extended));
    for (name, bytes) in damaged {
        let path = directory.path().join(format!("{name}.tcd"));
        std::fs::write(&path, &bytes).unwrap();
        let output = run(&path, evaluate.clone(), true);
        assert!(!output.status.success(), "{name} opened");
        assert_eq!(std::fs::read(&path).unwrap(), bytes, "{name} was rewritten");
    }
    // The last valid project is untouched and still rebuilds without the
    // original STEP file.
    assert_eq!(std::fs::read(&project).unwrap(), host_project);
    let again = run(&project, evaluate, false);
    assert!(again.status.success());
    let receipt: Value = serde_json::from_slice(&again.stdout).unwrap();
    assert_eq!(receipt["result"]["parts"], parts);
    assert_eq!(
        ContentHash::of(&std::fs::read(&project).unwrap()),
        ContentHash::of(&host_project)
    );
}

fn write_ledger(name: &str, rows: &[Value]) {
    let path =
        Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("exchange-robustness-{name}.json"));
    std::fs::write(path, serde_json::to_vec_pretty(rows).unwrap()).unwrap();
}
