//! Independent cross-CAD STEP exchange corpus (CAD-345).
//!
//! The manifest in `fixtures/exchange/corpus/manifest.json` pins every fixture
//! by producer, units, provenance, redistribution basis, SHA-256, geometry
//! tags and the outcome the Rust importer must produce. Vendored fixtures
//! always run; fixtures referenced from the installed KiCad 3D-model library
//! run only when that library holds identical bytes and otherwise report
//! `unavailable`, never a pass.
//!
//! Every admitted fixture is checked for volume, bounds, topology, an
//! independent Rust readback of the source text, a certified closed
//! manufacturing mesh, a planar section, material persistence, save and
//! source-free reopen, re-export, independent readback of the re-export and
//! the geometric deviation between source and readback meshes. With
//! `--features occt` the re-export is also read back by OpenCascade 7.9.3.
use multicad_core::{
    commands::{Actor, Channel, Session},
    document::{Document, Entity, EntityData, EntityId},
    drawings::{VectorCurve, VectorRole},
    materials::{self, Assignment, HatchStyle, LibraryDefinition},
    model::imports::{self, StepEncoding},
    project_file::{self, ContentHash, FileLimits, ProjectFile},
};
use multicad_kernel::{
    DrawingProjection, DrawingSolid, GeometryErrorKind, GeometryPolicy, MeshOptions, PlaneFrame,
    Solid, project_drawing, rust_kernel::step_read,
};
use multicad_model::{Compiler, GeometryBackend, imports::read_step_file_with_backend};
use serde_json::{Value, json};
use std::path::{Path, PathBuf};

const KICAD_LIBRARY: &str = "/Applications/KiCad/KiCad.app/Contents/SharedSupport/3dmodels";

/// Declared comparison tolerances. They are recorded in the emitted ledger.
const VOLUME_RELATIVE: f64 = 1e-9;
const BOUNDS_MM: f64 = 1e-6;
const ORACLE_VOLUME_RELATIVE: f64 = 1e-9;
const SECTION_LENGTH_RELATIVE: f64 = 1e-9;
/// Mesh-to-mesh deviation between source and readback, both tessellated
/// with the same manufacturing options, must stay within the requested
/// linear deflection.
const MESH_DEFLECTION_MM: f64 = 0.05;
/// A tessellation removes at most this fraction of the exact volume.
const MESH_VOLUME_RELATIVE: f64 = 2e-2;
const OCCT_VOLUME_RELATIVE: f64 = 1e-6;
const OCCT_BOUNDS_MM: f64 = 1e-3;

fn tolerances() -> Value {
    json!({
        "volume_relative": VOLUME_RELATIVE,
        "bounds_mm": BOUNDS_MM,
        "oracle_volume_relative": ORACLE_VOLUME_RELATIVE,
        "section_length_relative": SECTION_LENGTH_RELATIVE,
        "mesh_linear_deflection_mm": MESH_DEFLECTION_MM,
        "mesh_volume_relative": MESH_VOLUME_RELATIVE,
        "occt_volume_relative": OCCT_VOLUME_RELATIVE,
        "occt_bounds_mm": OCCT_BOUNDS_MM,
        "reader_tolerance_mm": GeometryPolicy::default().linear_tolerance_mm,
    })
}

fn kicad_library() -> PathBuf {
    std::env::var_os("MULTICAD_KICAD_3DMODELS")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(KICAD_LIBRARY))
}

fn outcome_of(result: &Result<Solid, multicad_kernel::GeometryError>) -> &'static str {
    match result {
        Ok(_) => "pass",
        Err(error) => match error.kind {
            GeometryErrorKind::Unsupported => "unsupported",
            GeometryErrorKind::Limit => "limit",
            GeometryErrorKind::InvalidResult | GeometryErrorKind::Input => "invalid",
            _ => "fail",
        },
    }
}

fn import(text: &str) -> Result<Solid, multicad_kernel::GeometryError> {
    Solid::import_rust_step(text, EntityId::new(), GeometryPolicy::default())
}

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/exchange")
}

fn repository() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn manifest() -> Value {
    let path = fixtures().join("corpus/manifest.json");
    serde_json::from_slice(&std::fs::read(path).unwrap()).unwrap()
}

/// Resolve an entry to its bytes, or `None` when an external library entry
/// is not installed with the pinned content.
fn entry_bytes(entry: &Value) -> Result<Vec<u8>, String> {
    entry_file(entry).map(|(_, bytes)| bytes)
}

fn entry_file(entry: &Value) -> Result<(PathBuf, Vec<u8>), String> {
    let relative = entry["path"].as_str().unwrap();
    let path = match entry["location"].as_str().unwrap() {
        "vendored" => fixtures().join(relative),
        "repository" => repository().join(relative),
        "kicad-library" => kicad_library().join(relative),
        other => panic!("unknown location {other}"),
    };
    let bytes = std::fs::read(&path).map_err(|e| format!("{}: {e}", path.display()))?;
    let digest = ContentHash::of(&bytes).to_string();
    if digest != entry["sha256"].as_str().unwrap() {
        return Err(format!("{} has sha256 {digest}", path.display()));
    }
    Ok((path, bytes))
}

fn close(actual: f64, expected: f64, relative: f64) -> bool {
    (actual - expected).abs() <= relative * expected.abs().max(1.)
}

fn same_bounds(a: [[f64; 3]; 2], b: [[f64; 3]; 2], tolerance: f64) -> bool {
    a.iter()
        .flatten()
        .zip(b.iter().flatten())
        .all(|(a, b)| (a - b).abs() <= tolerance)
}

fn bounds_of(solid: &Solid) -> [[f64; 3]; 2] {
    let bounds = solid.bounds().unwrap();
    [bounds.min, bounds.max]
}

/// Independent Part 21 readback of every member, without the importer's
/// reconstruction.
fn independent_readback(text: &str) -> Result<Value, step_read::ReadError> {
    let tolerance = GeometryPolicy::default().linear_tolerance_mm;
    let members = step_read::read_analytic_compound(text, tolerance)?;
    let mut bounds = [[f64::INFINITY; 3], [f64::NEG_INFINITY; 3]];
    for member in &members {
        for (axis, (low, high)) in member.bounds_mm[0]
            .iter()
            .zip(member.bounds_mm[1])
            .enumerate()
        {
            bounds[0][axis] = bounds[0][axis].min(*low);
            bounds[1][axis] = bounds[1][axis].max(high);
        }
    }
    Ok(json!({
        "members": members.len(),
        "volume_mm3": members.iter().map(|m| m.volume_mm3).sum::<f64>(),
        "faces": members.iter().map(|m| m.face_count).sum::<usize>(),
        "bounds_mm": bounds,
    }))
}

fn mesh_options() -> MeshOptions {
    MeshOptions {
        linear_deflection_mm: MESH_DEFLECTION_MM,
        ..Default::default()
    }
}

/// Certified manufacturing mesh: every edge is shared by exactly two
/// oppositely oriented triangles, and the enclosed volume agrees with the
/// exact solid.
fn checked_mesh(solid: &Solid, label: &str) -> (multicad_kernel::Mesh, Value) {
    let receipt = solid
        .manufacturing_mesh_with_receipt(mesh_options())
        .unwrap_or_else(|e| panic!("{label}: manufacturing mesh: {e}"));
    let mesh = receipt.mesh.clone();
    let triangles = mesh.indices.len() / 3;
    assert!(triangles > 0, "{label}: empty mesh");
    let mut directed = std::collections::BTreeMap::new();
    for t in mesh.indices.as_chunks::<3>().0 {
        for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            *directed.entry((a, b)).or_insert(0u32) += 1;
        }
    }
    for (&(a, b), &count) in &directed {
        assert_eq!(count, 1, "{label}: repeated directed edge {a}-{b}");
        assert_eq!(
            directed.get(&(b, a)),
            Some(&1),
            "{label}: open or non-manifold edge {a}-{b}"
        );
    }
    let vertex = |i: u32| {
        let i = i as usize * 3;
        [mesh.vertices[i], mesh.vertices[i + 1], mesh.vertices[i + 2]]
    };
    let volume: f64 = mesh
        .indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| {
            let [a, b, c] = [vertex(t[0]), vertex(t[1]), vertex(t[2])];
            (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0]))
                / 6.
        })
        .sum();
    let exact = solid.volume_mm3().unwrap();
    assert!(
        close(volume, exact, MESH_VOLUME_RELATIVE),
        "{label}: mesh volume {volume} differs from {exact}"
    );
    let summary = json!({
        "triangles": triangles,
        "edge_pairs": receipt.edge_pairs,
        "maximum_displacement_mm": receipt.maximum_displacement_mm,
        "volume_mm3": volume,
    });
    (mesh, summary)
}

fn point_triangle_distance(p: [f64; 3], [a, b, c]: [[f64; 3]; 3]) -> f64 {
    let sub = |x: [f64; 3], y: [f64; 3]| [x[0] - y[0], x[1] - y[1], x[2] - y[2]];
    let dot = |x: [f64; 3], y: [f64; 3]| x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
    let (ab, ac, ap) = (sub(b, a), sub(c, a), sub(p, a));
    let (d1, d2) = (dot(ab, ap), dot(ac, ap));
    let point = |q: [f64; 3]| dot(sub(p, q), sub(p, q)).sqrt();
    let along =
        |o: [f64; 3], d: [f64; 3], t: f64| [o[0] + d[0] * t, o[1] + d[1] * t, o[2] + d[2] * t];
    if d1 <= 0. && d2 <= 0. {
        return point(a);
    }
    let bp = sub(p, b);
    let (d3, d4) = (dot(ab, bp), dot(ac, bp));
    if d3 >= 0. && d4 <= d3 {
        return point(b);
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0. && d1 >= 0. && d3 <= 0. {
        return point(along(a, ab, d1 / (d1 - d3)));
    }
    let cp = sub(p, c);
    let (d5, d6) = (dot(ab, cp), dot(ac, cp));
    if d6 >= 0. && d5 <= d6 {
        return point(c);
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0. && d2 >= 0. && d6 <= 0. {
        return point(along(a, ac, d2 / (d2 - d6)));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0. && (d4 - d3) >= 0. && (d5 - d6) >= 0. {
        return point(along(b, sub(c, b), (d4 - d3) / ((d4 - d3) + (d5 - d6))));
    }
    let denominator = 1. / (va + vb + vc);
    let (v, w) = (vb * denominator, vc * denominator);
    point([
        a[0] + ab[0] * v + ac[0] * w,
        a[1] + ab[1] * v + ac[1] * w,
        a[2] + ab[2] * v + ac[2] * w,
    ])
}

/// Largest distance from any vertex of one mesh to the surface of the other,
/// in both directions. Distances beyond `radius` are not searched and report
/// infinity, which fails any tolerance up to that radius.
fn mesh_deviation(a: &multicad_kernel::Mesh, b: &multicad_kernel::Mesh, radius: f64) -> f64 {
    use std::collections::HashMap;
    let vertex = |m: &multicad_kernel::Mesh, i: u32| {
        let i = i as usize * 3;
        [m.vertices[i], m.vertices[i + 1], m.vertices[i + 2]]
    };
    let cell = |p: [f64; 3]| p.map(|x| (x / radius).floor() as i64);
    let one_way = |from: &multicad_kernel::Mesh, to: &multicad_kernel::Mesh| {
        // Bucket each triangle into every grid cell its bounds, expanded by
        // the search radius, overlap; a query then reads one cell.
        let mut grid: HashMap<[i64; 3], Vec<[[f64; 3]; 3]>> = HashMap::new();
        for t in to.indices.as_chunks::<3>().0 {
            let triangle = [vertex(to, t[0]), vertex(to, t[1]), vertex(to, t[2])];
            let low = cell(std::array::from_fn(|axis| {
                triangle
                    .iter()
                    .map(|p| p[axis])
                    .fold(f64::INFINITY, f64::min)
                    - radius
            }));
            let high = cell(std::array::from_fn(|axis| {
                triangle
                    .iter()
                    .map(|p| p[axis])
                    .fold(f64::NEG_INFINITY, f64::max)
                    + radius
            }));
            for x in low[0]..=high[0] {
                for y in low[1]..=high[1] {
                    for z in low[2]..=high[2] {
                        grid.entry([x, y, z]).or_default().push(triangle);
                    }
                }
            }
        }
        from.vertices
            .as_chunks::<3>()
            .0
            .iter()
            .map(|&p| {
                grid.get(&cell(p)).map_or(f64::INFINITY, |triangles| {
                    let nearest = triangles
                        .iter()
                        .map(|&t| point_triangle_distance(p, t))
                        .fold(f64::INFINITY, f64::min);
                    if nearest <= radius {
                        nearest
                    } else {
                        f64::INFINITY
                    }
                })
            })
            .fold(0., f64::max)
    };
    one_way(a, b).max(one_way(b, a))
}

/// Total length and count of the section boundary curves cut by a horizontal
/// plane through the middle of the bounds.
fn section(solid: &Solid, hatch: Option<HatchStyle>) -> Value {
    let [min, max] = bounds_of(solid);
    let height = 0.5 * (min[2] + max[2]);
    let hatch = hatch.map(|h| multicad_kernel::DrawingHatch {
        angle_rad: h.angle_deg.to_radians(),
        spacing_mm: h.spacing_mm,
        line_width_mm: h.line_width_mm,
        pattern: h.pattern,
    });
    let vectors = project_drawing(
        &[DrawingSolid { solid, hatch }],
        DrawingProjection {
            frame: PlaneFrame::default(),
            section_offset_mm: Some(height),
            section_only: true,
            ..Default::default()
        },
    );
    match vectors {
        Ok(vectors) => {
            let boundary: Vec<_> = vectors
                .iter()
                .filter(|v| v.role == VectorRole::Section)
                .collect();
            // Lines and circular arcs have exact lengths; other curves use
            // their published polyline and make the total approximate.
            let mut exact = true;
            let length: f64 = boundary
                .iter()
                .map(|v| match v.curve {
                    VectorCurve::Line { a, b } => (b[0] - a[0]).hypot(b[1] - a[1]),
                    VectorCurve::Ellipse {
                        major_mm,
                        minor_mm,
                        start_rad,
                        end_rad,
                        ..
                    } if major_mm == minor_mm => major_mm * (end_rad - start_rad).abs(),
                    _ => {
                        exact = false;
                        v.preview_mm
                            .windows(2)
                            .map(|w| (w[1][0] - w[0][0]).hypot(w[1][1] - w[0][1]))
                            .sum::<f64>()
                    }
                })
                .sum();
            json!({
                "plane_z_mm": height,
                "curves": boundary.len(),
                "length_mm": length,
                "length_exact": exact,
                "hatch_vectors": vectors.iter().filter(|v| v.role == VectorRole::Hatch).count(),
            })
        }
        Err(error) => json!({ "plane_z_mm": height, "error": error.to_string() }),
    }
}

fn file_name_strings(text: &str) -> Vec<String> {
    let Some(start) = text.find("FILE_NAME") else {
        return vec![];
    };
    let end = text[start..].find(");").map_or(text.len(), |e| start + e);
    text[start..end]
        .split('\'')
        .skip(1)
        .step_by(2)
        .map(String::from)
        .collect()
}

fn actor() -> Actor {
    Actor::local(Channel::Mcp)
}

fn compile_parts(document: &Document, parts: &[EntityId]) -> Vec<Solid> {
    let model = Compiler {
        omit_display: true,
        ..Default::default()
    }
    .compile_with_backend(document, 0, GeometryBackend::Rust, None)
    .unwrap();
    parts
        .iter()
        .map(|&id| {
            model
                .solid(id)
                .expect("compiled imported part")
                .share()
                .unwrap()
        })
        .collect()
}

/// Import through the canonical document, assign a material with a section
/// hatch, save the native project, delete the original file, reopen and
/// rebuild without it. Returns the reopened solid and the persistence receipt.
fn persist_and_reopen(entry: &Value, bytes: &[u8]) -> (Solid, Value) {
    let id = entry["id"].as_str().unwrap();
    let directory = tempfile::tempdir().unwrap();
    let source = directory.path().join(format!("{id}.step"));
    std::fs::write(&source, bytes).unwrap();
    let file = read_step_file_with_backend(&source, GeometryBackend::Rust).unwrap();
    assert_eq!(file.encoding, StepEncoding::RustStep);
    assert_eq!(file.content, ContentHash::of(bytes));
    let mut session = Session::new(Document::new(format!("Corpus {id}")).unwrap()).unwrap();
    let prepared = imports::prepare_step_import(&session, 0, actor(), &file).unwrap();
    let import = prepared.source;
    let parts: Vec<_> = prepared.model.outputs.iter().map(|o| o.body).collect();
    let features: Vec<_> = prepared.model.outputs.iter().map(|o| o.feature).collect();
    session.commit(prepared.model.prepared).unwrap();

    let material = Entity::new(
        "Corpus aluminium",
        EntityData::Material {
            density_kg_m3: Some(2700.),
            color_rgba: [190, 190, 200, 255],
            section_hatch: None,
        },
    );
    let hatch_style = HatchStyle {
        angle_deg: 45.,
        spacing_mm: 0.1,
        line_width_mm: 0.05,
        ..Default::default()
    };
    let hatch = LibraryDefinition::Hatch(hatch_style.clone())
        .entity(EntityId::new(), "Corpus section hatch")
        .unwrap();
    let library = materials::prepare_library(
        &session,
        session.revision(),
        actor(),
        &[material.clone(), hatch.clone()],
    )
    .unwrap();
    session.commit(library).unwrap();
    let assignment = materials::prepare_assignment(
        &session,
        session.revision(),
        actor(),
        &parts,
        &Assignment {
            material: Some(material.id),
            stock: None,
            appearance: None,
            hatch: Some(hatch.id),
        },
    )
    .unwrap();
    session.commit(assignment).unwrap();
    let before = compile_parts(session.document(), &features);

    let project = directory.path().join(format!("{id}.tcd"));
    let file = ProjectFile::new(session.document().clone()).unwrap();
    let mut encoded = vec![];
    project_file::write(&file, &mut encoded, FileLimits::default()).unwrap();
    std::fs::write(&project, &encoded).unwrap();
    std::fs::remove_file(&source).unwrap();
    assert!(!source.exists());

    let opened = project_file::read(
        &mut std::fs::read(&project).unwrap().as_slice(),
        FileLimits::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap();
    let document = opened.project.document;
    assert_eq!(&document, session.document(), "{id}: reopened document");
    let stored = imports::read_step_source(&document, import).unwrap();
    assert_eq!(
        stored.content,
        ContentHash::of(bytes),
        "{id}: embedded source"
    );
    assert_eq!(stored.encoding, StepEncoding::RustStep);
    let after = compile_parts(&document, &features);
    assert_eq!(after.len(), before.len());
    let mut mass = 0.;
    for (part, (old, new)) in parts.iter().zip(before.iter().zip(&after)) {
        assert_eq!(old.volume_mm3().unwrap(), new.volume_mm3().unwrap(), "{id}");
        assert_eq!(bounds_of(old), bounds_of(new), "{id}");
        let resolved = materials::resolve(&document, *part).unwrap();
        assert_eq!(resolved.material, Some(material.id), "{id}: material");
        assert_eq!(resolved.hatch.as_ref(), Some(&hatch_style), "{id}: hatch");
        assert!(
            resolved.diagnostics.is_empty(),
            "{id}: {:?}",
            resolved.diagnostics
        );
        mass += 2700e-9 * new.volume_mm3().unwrap();
    }
    assert_eq!(after.len(), 1, "{id}: one Part per Rust STEP source");
    let solid = after.into_iter().next().unwrap();
    let receipt = json!({
        "project_bytes": encoded.len(),
        "parts": parts.len(),
        "material": "Corpus aluminium 2700 kg/m3",
        "mass_kg": mass,
        "section": section(&solid, Some(hatch_style)),
    });
    (solid, receipt)
}

/// A section either succeeds and matches the analytic perimeter when one is
/// pinned, or fails with the exact error the manifest attributes to its owner.
fn check_section(entry: &Value, section: &Value) {
    let id = entry["id"].as_str().unwrap();
    match (
        section["error"].as_str(),
        entry["section_finding"].as_object(),
    ) {
        (None, None) => {
            if let Some(expected) = entry["oracle"]["section_length_mm"].as_f64() {
                assert_eq!(section["length_exact"], true, "{id}");
                let length = section["length_mm"].as_f64().unwrap();
                assert!(
                    close(length, expected, SECTION_LENGTH_RELATIVE),
                    "{id}: section length {length} differs from {expected}"
                );
            }
        }
        (Some(error), Some(finding)) => assert_eq!(
            error,
            finding["error"].as_str().unwrap(),
            "{id}: section error changed"
        ),
        (error, finding) => panic!("{id}: section error {error:?} but pinned finding {finding:?}"),
    }
}

/// Run one manifest entry, filling its ledger row as checks complete so a
/// failed check still leaves the measured evidence in the row.
fn run_entry(entry: &Value, row: &mut Value) {
    let id = entry["id"].as_str().unwrap();
    *row = json!({
        "id": id,
        "location": entry["location"],
        "path": entry["path"],
        "sha256": entry["sha256"],
        "producer": entry["producer"],
        "tags": entry["tags"],
        "expected": entry["expected"],
        "verdict": entry["verdict"],
        "finding": entry["finding"],
    });
    let (path, bytes) = match entry_file(entry) {
        Ok(file) => file,
        Err(reason) => {
            assert_eq!(
                entry["location"], "kicad-library",
                "{id}: vendored and repository fixtures must be present: {reason}"
            );
            row["outcome"] = json!("unavailable");
            row["reason"] = json!(reason);
            return;
        }
    };
    let text = String::from_utf8(bytes.clone()).expect("Part 21 text");
    let strings = file_name_strings(&text);
    for expected in entry["producer"]["file_name"].as_array().unwrap() {
        assert!(
            strings.iter().any(|s| s == expected.as_str().unwrap()),
            "{id}: FILE_NAME lacks {expected}: {strings:?}"
        );
    }
    assert_eq!(entry["units"], "mm");
    assert!(
        text.contains("SI_UNIT(.MILLI.,.METRE.)"),
        "{id}: not a millimetre file"
    );
    row["source_native"] = native_source(&path);
    if cfg!(feature = "occt") && entry["verdict"] == "fail" {
        // A refusal is only attributed as an incorrect result when the
        // independent reader accepts the same file as a valid solid.
        assert!(
            row["source_native"]["volume_mm3"]
                .as_f64()
                .is_some_and(|v| v > 0.),
            "{id}: OpenCascade does not confirm the failure verdict: {}",
            row["source_native"]
        );
    }
    let started = std::time::Instant::now();
    let result = import(&text);
    let import_ms = started.elapsed().as_secs_f64() * 1e3;
    let outcome = outcome_of(&result);
    row["outcome"] = json!(outcome);
    row["import_ms"] = json!(import_ms);
    let solid = match result {
        Ok(solid) => solid,
        Err(error) => {
            row["message"] = json!(error.to_string());
            return;
        }
    };
    let volume = solid.volume_mm3().unwrap();
    let bounds = bounds_of(&solid);
    if let Some(native) = row["source_native"]["volume_mm3"].as_f64() {
        // Cross-CAD agreement of the imported source with OpenCascade.
        assert!(
            close(native, volume, OCCT_VOLUME_RELATIVE),
            "{id}: OpenCascade source volume {native}"
        );
        let native_bounds: [[f64; 3]; 2] =
            serde_json::from_value(row["source_native"]["bounds_mm"].clone()).unwrap();
        assert!(
            same_bounds(native_bounds, bounds, OCCT_BOUNDS_MM),
            "{id}: OpenCascade source bounds"
        );
    }
    let faces = solid.faces().unwrap().len();
    let solids = solid.solid_count().unwrap();
    row["source"] = json!({
        "volume_mm3": volume,
        "bounds_mm": bounds,
        "faces": faces,
        "solids": solids,
    });
    if let Some(oracle) = entry.get("oracle") {
        let expected = oracle["volume_mm3"].as_f64().unwrap();
        assert!(
            close(volume, expected, ORACLE_VOLUME_RELATIVE),
            "{id}: volume {volume} differs from the analytic {expected}"
        );
        let expected_bounds: [[f64; 3]; 2] =
            serde_json::from_value(oracle["bounds_mm"].clone()).unwrap();
        assert!(
            same_bounds(bounds, expected_bounds, BOUNDS_MM),
            "{id}: bounds {bounds:?} differ from {expected_bounds:?}"
        );
        assert_eq!(
            faces as u64,
            oracle["faces"].as_u64().unwrap(),
            "{id}: faces"
        );
        assert_eq!(
            u64::from(solids),
            oracle["solids"].as_u64().unwrap(),
            "{id}: solids"
        );
    }
    // An independent parse of the source boundary agrees with the
    // reconstructed solid.
    let reader = independent_readback(&text).unwrap();
    assert!(
        close(
            reader["volume_mm3"].as_f64().unwrap(),
            volume,
            VOLUME_RELATIVE
        ),
        "{id}"
    );
    assert_eq!(reader["members"], u64::from(solids), "{id}");
    let reader_bounds: [[f64; 3]; 2] = serde_json::from_value(reader["bounds_mm"].clone()).unwrap();
    assert!(
        same_bounds(reader_bounds, bounds, BOUNDS_MM),
        "{id}: reader bounds"
    );
    row["source_reader"] = reader;
    let (source_mesh, mesh) = checked_mesh(&solid, id);
    row["mesh"] = mesh;
    let source_section = section(&solid, None);
    check_section(entry, &source_section);
    row["section"] = source_section.clone();

    let (reopened, persistence) = persist_and_reopen(entry, &bytes);
    assert_eq!(
        reopened.volume_mm3().unwrap(),
        volume,
        "{id}: reopened volume"
    );
    assert_eq!(bounds_of(&reopened), bounds, "{id}: reopened bounds");
    row["persistence"] = persistence.clone();
    check_section(entry, &persistence["section"]);
    if persistence["section"]["error"].is_null() {
        assert!(
            persistence["section"]["hatch_vectors"].as_u64().unwrap() > 0,
            "{id}: the assigned material hatch is missing from the section"
        );
    }

    let directory = tempfile::tempdir().unwrap();
    let exported = directory.path().join(format!("{id}-export.step"));
    let started = std::time::Instant::now();
    reopened.write_step(&exported).unwrap();
    row["export_ms"] = json!(started.elapsed().as_secs_f64() * 1e3);
    let exported_text = std::fs::read_to_string(&exported).unwrap();
    let readback = independent_readback(&exported_text).unwrap();
    let readback_solid = import(&exported_text).unwrap();
    let readback_volume = readback_solid.volume_mm3().unwrap();
    let readback_bounds = bounds_of(&readback_solid);
    let readback_faces = readback_solid.faces().unwrap().len();
    assert!(
        close(
            readback["volume_mm3"].as_f64().unwrap(),
            volume,
            VOLUME_RELATIVE
        ),
        "{id}: readback volume"
    );
    assert!(
        close(readback_volume, volume, VOLUME_RELATIVE),
        "{id}: reimported volume"
    );
    assert!(
        same_bounds(readback_bounds, bounds, BOUNDS_MM),
        "{id}: readback bounds {readback_bounds:?}"
    );
    assert_eq!(readback_faces, faces, "{id}: readback faces");
    assert_eq!(
        readback_solid.solid_count().unwrap(),
        solids,
        "{id}: readback solids"
    );
    let (readback_mesh, _) = checked_mesh(&readback_solid, id);
    let deviation = mesh_deviation(&source_mesh, &readback_mesh, 4. * MESH_DEFLECTION_MM);
    assert!(
        deviation <= MESH_DEFLECTION_MM,
        "{id}: mesh deviation {deviation}"
    );
    let readback_section = section(&readback_solid, None);
    if let (Some(a), Some(b)) = (
        source_section["length_mm"].as_f64(),
        readback_section["length_mm"].as_f64(),
    ) {
        assert!(
            close(a, b, SECTION_LENGTH_RELATIVE),
            "{id}: section {a} != {b}"
        );
        assert_eq!(source_section["curves"], readback_section["curves"], "{id}");
    }
    row["readback"] = json!({
        "export_bytes": exported_text.len(),
        "reader": readback,
        "volume_mm3": readback_volume,
        "bounds_mm": readback_bounds,
        "faces": readback_faces,
        "volume_deviation_mm3": (readback_volume - volume).abs(),
        "bounds_deviation_mm": readback_bounds.iter().flatten().zip(bounds.iter().flatten()).map(|(a, b)| (a - b).abs()).fold(0., f64::max),
        "mesh_deviation_mm": deviation,
        "section": readback_section,
    });
    row["native_readback"] = native_readback(&exported, volume, bounds, solids, &readback_mesh);
}

#[cfg(not(feature = "occt"))]
fn native_source(_: &Path) -> Value {
    json!("not built; run with --features occt")
}

/// OpenCascade 7.9.3 reads the original source file as an independent
/// validity and mass reference. Its result never changes the Rust outcome.
#[cfg(feature = "occt")]
fn native_source(path: &Path) -> Value {
    match Solid::read_step(path, GeometryPolicy::default()) {
        Ok(native) => match (native.volume_mm3(), native.solid_count()) {
            (Ok(volume), Ok(solids)) => json!({
                "reader": "OpenCascade 7.9.3",
                "volume_mm3": volume,
                "bounds_mm": bounds_of(&native),
                "solids": solids,
                "faces": native.faces().map(|f| f.len()).ok(),
            }),
            (volume, solids) => json!({
                "reader": "OpenCascade 7.9.3",
                "error": format!("{:?} {:?}", volume.err(), solids.err()),
            }),
        },
        Err(error) => json!({ "reader": "OpenCascade 7.9.3", "error": error.to_string() }),
    }
}

#[cfg(not(feature = "occt"))]
fn native_readback(_: &Path, _: f64, _: [[f64; 3]; 2], _: u32, _: &multicad_kernel::Mesh) -> Value {
    json!("not built; run with --features occt")
}

/// OpenCascade 7.9.3 reads the Rust re-export as reference evidence.
#[cfg(feature = "occt")]
fn native_readback(
    path: &Path,
    volume: f64,
    bounds: [[f64; 3]; 2],
    solids: u32,
    rust_mesh: &multicad_kernel::Mesh,
) -> Value {
    let native = Solid::read_step(path, GeometryPolicy::default()).unwrap();
    // Independent tessellations of the same exact boundary: each lies
    // within the requested deflection of the surface, so their mutual
    // deviation is bounded by twice that deflection.
    let native_mesh = native.mesh(mesh_options()).unwrap();
    let deviation = mesh_deviation(rust_mesh, &native_mesh, 4. * MESH_DEFLECTION_MM);
    assert!(
        deviation <= 2. * MESH_DEFLECTION_MM,
        "OpenCascade tessellation deviates {deviation} mm from the Rust readback"
    );
    assert!(!native.uses_rust_geometry());
    let native_volume = native.volume_mm3().unwrap();
    let native_bounds = bounds_of(&native);
    assert!(
        close(native_volume, volume, OCCT_VOLUME_RELATIVE),
        "OpenCascade volume {native_volume} != {volume}"
    );
    assert!(
        same_bounds(native_bounds, bounds, OCCT_BOUNDS_MM),
        "OpenCascade bounds {native_bounds:?}"
    );
    assert_eq!(native.solid_count().unwrap(), solids);
    json!({
        "reader": "OpenCascade 7.9.3",
        "mesh_deviation_mm": deviation,
        "volume_mm3": native_volume,
        "bounds_mm": native_bounds,
        "solids": native.solid_count().unwrap(),
        "faces": native.faces().unwrap().len(),
    })
}

fn git_head() -> Value {
    let output = std::process::Command::new("git")
        .args(["rev-parse", "HEAD"])
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .output();
    let dirty = std::process::Command::new("git")
        .args(["status", "--porcelain", "--untracked-files=no"])
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .output()
        .map(|o| !o.stdout.is_empty())
        .ok();
    json!({
        "commit": output.ok().map(|o| String::from_utf8_lossy(&o.stdout).trim().to_owned()),
        "tracked_changes": dirty,
    })
}

/// The pinned ledger: every entry's import outcome matches the manifest and
/// every admitted entry passes the full exchange workflow. The run writes
/// `exchange-corpus-ledger.json` under the Cargo test temporary directory,
/// or to `MULTICAD_CORPUS_LEDGER` when set.
#[test]
fn corpus_matches_the_pinned_ledger() {
    let manifest = manifest();
    let mut rows = vec![];
    let mut mismatches = vec![];
    for entry in manifest["entries"].as_array().unwrap() {
        let mut row = json!({ "id": entry["id"], "sha256": entry["sha256"] });
        let result =
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| run_entry(entry, &mut row)));
        let row = match result {
            Ok(()) => row,
            Err(panic) => {
                let message = panic
                    .downcast_ref::<String>()
                    .cloned()
                    .or_else(|| panic.downcast_ref::<&str>().map(|s| (*s).to_owned()))
                    .unwrap_or_default();
                mismatches.push(format!("{}: check failed: {message}", entry["id"]));
                row["import_outcome"] = row["outcome"].clone();
                row["outcome"] = json!("check-failed");
                row["check_failure"] = json!(message);
                row
            }
        };
        let outcome = row["outcome"].as_str().unwrap();
        if !matches!(outcome, "unavailable" | "check-failed") && row["outcome"] != entry["expected"]
        {
            mismatches.push(format!(
                "{}: expected {} but import was {outcome}: {}",
                entry["id"], entry["expected"], row["message"]
            ));
        }
        rows.push(row);
    }
    let mut totals = std::collections::BTreeMap::<String, usize>::new();
    let mut verdicts = std::collections::BTreeMap::<String, usize>::new();
    for row in &rows {
        *totals
            .entry(row["outcome"].as_str().unwrap().into())
            .or_default() += 1;
        if row["outcome"] != "unavailable" {
            *verdicts
                .entry(row["verdict"].as_str().unwrap().into())
                .or_default() += 1;
        }
    }
    let ledger = json!({
        "schema": 1,
        "item": "CAD-345",
        "source": git_head(),
        "native_reference": cfg!(feature = "occt"),
        "command": if cfg!(feature = "occt") {
            "CARGO_TARGET_DIR=target/occt cargo test --locked -p multicad-app --features occt --test exchange_corpus corpus_matches_the_pinned_ledger"
        } else {
            "cargo test --locked -p multicad-app --test exchange_corpus corpus_matches_the_pinned_ledger"
        },
        "kicad_library": kicad_library().display().to_string(),
        "tolerances": tolerances(),
        "outcomes": totals,
        "verdicts": verdicts,
        "entries": rows,
    });
    let path = std::env::var_os("MULTICAD_CORPUS_LEDGER")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            Path::new(env!("CARGO_TARGET_TMPDIR")).join("exchange-corpus-ledger.json")
        });
    let mut text = serde_json::to_string_pretty(&ledger).unwrap();
    text.push('\n');
    std::fs::write(&path, text).unwrap();
    println!(
        "outcomes {totals:?} verdicts {verdicts:?}; ledger {}",
        path.display()
    );
    assert!(mismatches.is_empty(), "{mismatches:#?}");
}

/// The manifest records provenance and redistribution for every fixture,
/// covers the required strata with always-available fixtures, and links
/// every refusal to its implementation owner.
#[test]
fn corpus_manifest_pins_provenance_strata_and_owners() {
    let manifest = manifest();
    let entries = manifest["entries"].as_array().unwrap();
    let mut ids = std::collections::BTreeSet::new();
    let mut vendored_tags = std::collections::BTreeSet::new();
    for entry in entries {
        let id = entry["id"].as_str().unwrap();
        assert!(ids.insert(id), "duplicate {id}");
        for field in [
            "location", "path", "sha256", "units", "source", "license", "expected", "verdict",
        ] {
            assert!(
                entry[field].as_str().is_some_and(|s| !s.is_empty()),
                "{id}: {field}"
            );
        }
        assert!(
            entry["producer"]["name"].as_str().is_some(),
            "{id}: producer"
        );
        assert!(
            entry["producer"]["version"].as_str().is_some(),
            "{id}: producer version"
        );
        assert!(
            !entry["producer"]["file_name"]
                .as_array()
                .unwrap()
                .is_empty(),
            "{id}"
        );
        assert!(matches!(
            entry["expected"].as_str().unwrap(),
            "pass" | "unsupported" | "invalid" | "limit"
        ));
        let verdict = entry["verdict"].as_str().unwrap();
        assert!(
            matches!(verdict, "pass" | "unsupported" | "fail"),
            "{id}: {verdict}"
        );
        assert_eq!(verdict == "pass", entry["expected"] == "pass", "{id}");
        if verdict != "pass" {
            let owner = entry["finding"]["owner"].as_str().unwrap_or_default();
            assert!(
                owner.starts_with("CAD-") && owner[4..].parse::<u32>().is_ok(),
                "{id}: refusal without an owner"
            );
            assert!(
                entry["finding"]["summary"]
                    .as_str()
                    .is_some_and(|s| s.len() > 20),
                "{id}"
            );
        }
        if let Some(finding) = entry["section_finding"].as_object() {
            let owner = finding["owner"].as_str().unwrap();
            assert!(
                owner.starts_with("CAD-") && owner[4..].parse::<u32>().is_ok(),
                "{id}"
            );
            assert!(
                finding["error"].is_string() && finding["summary"].is_string(),
                "{id}"
            );
        }
        if entry["location"] != "kicad-library" {
            entry_bytes(entry).unwrap_or_else(|e| panic!("{id}: {e}"));
            for tag in entry["tags"].as_array().unwrap() {
                vendored_tags.insert(tag.as_str().unwrap().to_owned());
            }
        }
    }
    for stratum in [
        "box",
        "cylinder",
        "cone",
        "sphere",
        "torus",
        "revolution",
        "spline",
        "trims-seams",
        "voids",
        "placement",
        "multi-body",
        "board",
    ] {
        assert!(
            vendored_tags.contains(stratum),
            "missing always-available stratum {stratum}"
        );
    }
    // The committed ledger outcomes agree with the manifest.
    let ledger: Value = serde_json::from_slice(
        &std::fs::read(repository().join("docs/architecture/exchange-corpus-ledger.json")).unwrap(),
    )
    .unwrap();
    let rows = ledger["entries"].as_array().unwrap();
    assert_eq!(
        rows.len(),
        entries.len(),
        "ledger and manifest entry counts"
    );
    for (row, entry) in rows.iter().zip(entries) {
        assert_eq!(row["id"], entry["id"]);
        assert_eq!(row["sha256"], entry["sha256"]);
        assert!(
            row["outcome"] == entry["expected"] || row["outcome"] == "unavailable",
            "{}",
            entry["id"]
        );
    }
}

fn step_files(root: &Path, output: &mut Vec<PathBuf>) {
    let Ok(entries) = std::fs::read_dir(root) else {
        return;
    };
    let mut entries: Vec<_> = entries.flatten().map(|e| e.path()).collect();
    entries.sort();
    for path in entries {
        if path.is_dir() {
            step_files(&path, output);
        } else if path
            .extension()
            .is_some_and(|e| e.eq_ignore_ascii_case("step") || e.eq_ignore_ascii_case("stp"))
        {
            output.push(path);
        }
    }
}

/// Producer from the Part 21 FILE_NAME originating system and preprocessor.
fn producer(text: &str) -> String {
    let Some(start) = text.find("FILE_NAME") else {
        return "unknown".into();
    };
    let header = &text[start..text[start..].find(");").map_or(text.len(), |e| start + e)];
    let strings: Vec<&str> = header.split('\'').skip(1).step_by(2).collect();
    let preprocessor = strings.get(strings.len().saturating_sub(3)).copied();
    let system = strings.get(strings.len().saturating_sub(2)).copied();
    let joined = format!("{} | {}", system.unwrap_or(""), preprocessor.unwrap_or(""));
    let lower = joined.to_lowercase();
    if system.unwrap_or("").is_empty() && preprocessor.unwrap_or("").is_empty() {
        "unstated".into()
    } else if lower.contains("freecad") {
        "FreeCAD".into()
    } else if lower.contains("open cascade") || lower.contains("opencascade") {
        "OpenCascade".into()
    } else if lower.contains("solidworks") {
        "SolidWorks".into()
    } else if lower.contains("kicad") {
        "KiCad".into()
    } else {
        joined
    }
}

fn categories(text: &str) -> Vec<&'static str> {
    let count = |needle: &str| text.matches(needle).count();
    let mut tags = vec![];
    if count("MANIFOLD_SOLID_BREP") + count("BREP_WITH_VOIDS") > 1 {
        tags.push("multi-body");
    }
    if count("BREP_WITH_VOIDS") > 0 {
        tags.push("voids");
    }
    if count("B_SPLINE_SURFACE") > 0 {
        tags.push("spline");
    }
    if count("TOROIDAL_SURFACE") > 0 {
        tags.push("torus");
    }
    if count("SPHERICAL_SURFACE") > 0 {
        tags.push("sphere");
    }
    if count("CONICAL_SURFACE") > 0 {
        tags.push("cone");
    }
    if count("SURFACE_OF_REVOLUTION") > 0 {
        tags.push("revolution");
    }
    if count("SHELL_BASED_SURFACE_MODEL") > 0 {
        tags.push("surface-model");
    }
    if count("NEXT_ASSEMBLY_USAGE_OCCURRENCE") > 0 {
        tags.push("assembly");
    }
    if tags.is_empty() {
        if count("CYLINDRICAL_SURFACE") > 0 {
            tags.push("cylinder");
        } else {
            tags.push("planar");
        }
    }
    tags
}

/// Aggregate sweep over the installed KiCad 3D-model library.
///
/// `cargo test --locked -p multicad-app --test exchange_corpus -- --ignored
/// kicad_library_sweep --nocapture` writes `kicad-sweep.json` under the Cargo
/// test temporary directory and prints pass/unsupported/fail counts per
/// producer and category. Set `MULTICAD_KICAD_3DMODELS` to another library
/// root, and `MULTICAD_KICAD_SWEEP_LIMIT` to sample the first N files.
#[test]
#[ignore = "broad sweep of an external library; run explicitly"]
fn kicad_library_sweep() {
    use std::collections::BTreeMap;
    let root = kicad_library();
    let mut files = vec![];
    step_files(&root, &mut files);
    if files.is_empty() {
        println!("unavailable: no STEP files under {}", root.display());
        return;
    }
    if let Some(limit) = std::env::var("MULTICAD_KICAD_SWEEP_LIMIT")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
    {
        files.truncate(limit);
    }
    let mut rows = vec![];
    let mut by_producer: BTreeMap<(String, &'static str), usize> = BTreeMap::new();
    let mut by_category: BTreeMap<(&'static str, &'static str), usize> = BTreeMap::new();
    for path in &files {
        let bytes = std::fs::read(path).unwrap();
        let relative = path.strip_prefix(&root).unwrap().display().to_string();
        let text = String::from_utf8_lossy(&bytes);
        let started = std::time::Instant::now();
        let result = std::panic::catch_unwind(|| import(&text));
        let millis = started.elapsed().as_millis();
        let (outcome, message) = match &result {
            Ok(result) => (
                outcome_of(result),
                result.as_ref().err().map(ToString::to_string),
            ),
            Err(_) => ("fail", Some("panic".into())),
        };
        let producer = producer(&text);
        let tags = categories(&text);
        *by_producer.entry((producer.clone(), outcome)).or_default() += 1;
        for tag in &tags {
            *by_category.entry((tag, outcome)).or_default() += 1;
        }
        rows.push(serde_json::json!({
            "path": relative,
            "sha256": ContentHash::of(&bytes).to_string(),
            "producer": producer,
            "tags": tags,
            "outcome": outcome,
            "message": message,
            "millis": millis,
        }));
    }
    let mut totals: BTreeMap<&str, usize> = BTreeMap::new();
    for row in &rows {
        *totals.entry(row["outcome"].as_str().unwrap()).or_default() += 1;
    }
    let summary = serde_json::json!({
        "library": root.display().to_string(),
        "files": rows.len(),
        "totals": totals,
        "by_producer": by_producer.iter().map(|((p, o), n)| serde_json::json!([p, o, n])).collect::<Vec<_>>(),
        "by_category": by_category.iter().map(|((c, o), n)| serde_json::json!([c, o, n])).collect::<Vec<_>>(),
        "rows": rows,
    });
    let out = Path::new(env!("CARGO_TARGET_TMPDIR")).join("kicad-sweep.json");
    std::fs::write(&out, serde_json::to_vec_pretty(&summary).unwrap()).unwrap();
    println!("totals {totals:?}");
    for ((producer, outcome), n) in &by_producer {
        println!("producer {producer:?} {outcome} {n}");
    }
    for ((category, outcome), n) in &by_category {
        println!("category {category} {outcome} {n}");
    }
    println!("wrote {}", out.display());
    // Every outcome must be a classified refusal or a pass; panics and
    // unclassified errors fail the sweep.
    assert_eq!(
        totals.get("fail"),
        None,
        "unclassified failures in the sweep"
    );
}

/// Developer probe: `MULTICAD_STEP_PROBE=a.step:b.step cargo test --locked -p
/// multicad-app --test exchange_corpus -- --ignored probe --nocapture` prints
/// the Rust import outcome, volume and bounds of each file.
#[test]
#[ignore = "developer probe for explicit files"]
fn probe() {
    let Some(paths) = std::env::var_os("MULTICAD_STEP_PROBE") else {
        return;
    };
    for path in std::env::split_paths(&paths) {
        let text = std::fs::read_to_string(&path).unwrap();
        let result = import(&text);
        let outcome = outcome_of(&result);
        match result {
            Ok(solid) => println!(
                "{}: {outcome} volume {:?} bounds {:?} solids {:?} faces {:?}",
                path.display(),
                solid.volume_mm3(),
                solid.bounds(),
                solid.solid_count(),
                solid.faces().map(|f| f.len()),
            ),
            Err(error) => println!("{}: {outcome} {error}", path.display()),
        }
    }
}

/// Regenerate the vendored OpenCascade 7.9.3 producer fixtures. Run only with
/// the optional native reference build:
///
/// `CARGO_TARGET_DIR=target/occt cargo test --locked -p multicad-app
/// --features occt --test exchange_corpus -- --ignored
/// regenerate_opencascade_fixtures`
///
/// The files are committed; the default build never runs OpenCascade.
#[cfg(feature = "occt")]
#[test]
#[ignore = "rewrites vendored fixtures with the optional OpenCascade build"]
fn regenerate_opencascade_fixtures() {
    use multicad_kernel::{
        BSpline2d, BooleanOperation, CurvedEdge2d, CurvedProfile, PlaneFrame, RevolutionAxis,
        RigidTransform,
    };
    let policy = GeometryPolicy::default();
    let directory = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/exchange/corpus");
    std::fs::create_dir_all(&directory).unwrap();
    let write = |name: &str, solid: &Solid| {
        solid.write_step(&directory.join(name)).unwrap();
    };
    let boxed = |size: [f64; 3], at: [f64; 3]| {
        Solid::box_solid(size, policy)
            .unwrap()
            .transformed(RigidTransform {
                translation_mm: at,
                ..Default::default()
            })
            .unwrap()
    };
    let subtract = |a: &Solid, b: &Solid| {
        a.boolean(b, BooleanOperation::Subtract, policy)
            .unwrap()
            .unwrap()
    };
    let (sin, cos) = 30_f64.to_radians().sin_cos();
    write(
        "occt-placed-box.step",
        &Solid::box_solid([24., 16., 9.], policy)
            .unwrap()
            .transformed(RigidTransform {
                rotation: [[cos, -sin, 0.], [sin, cos, 0.], [0., 0., 1.]],
                translation_mm: [50., -20., 7.],
            })
            .unwrap(),
    );
    write(
        "occt-cylinder.step",
        &Solid::cylinder(
            PlaneFrame {
                origin: [3., 4., 5.],
                u: [1., 0., 0.],
                v: [0., 1., 0.],
            },
            6.,
            15.,
            policy,
        )
        .unwrap(),
    );
    let xz = PlaneFrame {
        origin: [0.; 3],
        u: [1., 0., 0.],
        v: [0., 0., 1.],
    };
    let z_axis = RevolutionAxis {
        origin_mm: [0.; 3],
        direction: [0., 0., 1.],
    };
    let line = |a: [f64; 2], b: [f64; 2]| CurvedEdge2d::Line { a, b };
    let revolve = |outer: Vec<CurvedEdge2d>| {
        Solid::revolve_curved_profile(
            &CurvedProfile {
                frame: xz,
                outer,
                holes: vec![],
            },
            z_axis,
            std::f64::consts::TAU,
            policy,
        )
        .unwrap()
    };
    write(
        "occt-cone-frustum.step",
        &revolve(vec![
            line([0., 0.], [8., 0.]),
            line([8., 0.], [4., 10.]),
            line([4., 10.], [0., 10.]),
            line([0., 10.], [0., 0.]),
        ]),
    );
    write(
        "occt-sphere.step",
        &revolve(vec![
            line([0., -5.], [0., 5.]),
            CurvedEdge2d::Arc {
                center: [0., 0.],
                a: [0., 5.],
                b: [0., -5.],
                clockwise: true,
            },
        ]),
    );
    write(
        "occt-torus.step",
        &revolve(vec![CurvedEdge2d::Circle {
            center: [12., 0.],
            radius_mm: 3.,
        }]),
    );
    write(
        "occt-revolved-vase.step",
        &revolve(vec![
            line([0., 0.], [10., 0.]),
            line([10., 0.], [10., 5.]),
            CurvedEdge2d::Arc {
                center: [6., 5.],
                a: [10., 5.],
                b: [6., 9.],
                clockwise: false,
            },
            line([6., 9.], [6., 14.]),
            line([6., 14.], [0., 14.]),
            line([0., 14.], [0., 0.]),
        ]),
    );
    write(
        "occt-spline-prism.step",
        &Solid::extrude_curved_profile(
            &CurvedProfile {
                frame: PlaneFrame::default(),
                outer: vec![
                    line([0., 0.], [20., 0.]),
                    CurvedEdge2d::BSpline(BSpline2d {
                        degree: 3,
                        poles: vec![[20., 0.], [20., 12.], [0., 12.], [0., 0.]],
                        knots: vec![0., 1.],
                        multiplicities: vec![4, 4],
                        weights: None,
                    }),
                ],
                holes: vec![],
            },
            6.,
            policy,
        )
        .unwrap(),
    );
    write(
        "occt-box-void.step",
        &subtract(
            &boxed([30., 20., 10.], [0.; 3]),
            &boxed([10., 6., 4.], [10., 7., 3.]),
        ),
    );
    let cylinder = Solid::cylinder(
        PlaneFrame {
            origin: [30., 0., 0.],
            ..Default::default()
        },
        4.,
        8.,
        policy,
    )
    .unwrap();
    let block = boxed([10., 10., 5.], [0.; 3]);
    write(
        "occt-two-bodies.step",
        &Solid::compound_named(&[&block, &cylinder], EntityId::new(), policy).unwrap(),
    );
    let bore = Solid::cylinder(
        PlaneFrame {
            origin: [10., 10., 3.],
            ..Default::default()
        },
        3.,
        6.,
        policy,
    )
    .unwrap();
    write(
        "occt-pocketed-plate.step",
        &subtract(
            &subtract(&boxed([40., 30., 8.], [0.; 3]), &bore),
            &boxed([12., 8., 4.], [22., 15., 5.]),
        ),
    );
}

/// Minimal newline-delimited JSON-RPC client for the real executable.
struct Client {
    child: std::process::Child,
    stdin: std::process::ChildStdin,
    stdout: std::io::BufReader<std::process::ChildStdout>,
    next: u64,
}
impl Client {
    fn start(arguments: &[&str]) -> Self {
        use std::process::{Command, Stdio};
        let mut child = Command::new(env!("CARGO_BIN_EXE_multicad"))
            .arg("--mcp")
            .args(arguments)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::inherit())
            .spawn()
            .expect("start multicad --mcp");
        let stdin = child.stdin.take().unwrap();
        let stdout = std::io::BufReader::new(child.stdout.take().unwrap());
        let mut client = Self {
            child,
            stdin,
            stdout,
            next: 1,
        };
        let init = client.request(
            "initialize",
            json!({
                "protocolVersion": "2025-06-18",
                "capabilities": {},
                "clientInfo": { "name": "exchange-corpus", "version": "1.0" }
            }),
        );
        assert_eq!(init["result"]["serverInfo"]["name"], "multicad");
        client.send(json!({"jsonrpc": "2.0", "method": "notifications/initialized", "params": {}}));
        client
    }
    fn send(&mut self, message: Value) {
        use std::io::Write;
        writeln!(self.stdin, "{message}").unwrap();
        self.stdin.flush().unwrap();
    }
    fn request(&mut self, method: &str, params: Value) -> Value {
        use std::io::BufRead;
        let id = self.next;
        self.next += 1;
        self.send(json!({ "jsonrpc": "2.0", "id": id, "method": method, "params": params }));
        loop {
            let mut line = String::new();
            assert!(
                self.stdout.read_line(&mut line).unwrap() > 0,
                "server closed stdout"
            );
            let message: Value = serde_json::from_str(&line).unwrap();
            if message["id"] == json!(id) {
                return message;
            }
        }
    }
    fn call(&mut self, name: &str, arguments: Value) -> (bool, Value) {
        let response = self.request(
            "tools/call",
            json!({ "name": name, "arguments": arguments }),
        );
        assert!(response["error"].is_null(), "{response}");
        let result = &response["result"];
        (
            result["isError"].as_bool().unwrap_or(false),
            result["structuredContent"].clone(),
        )
    }
    fn ok(&mut self, name: &str, arguments: Value) -> Value {
        let (error, value) = self.call(name, arguments);
        assert!(!error, "{name} failed: {value}");
        value
    }
    fn fail(&mut self, name: &str, arguments: Value) -> Value {
        let (error, value) = self.call(name, arguments);
        assert!(error, "{name} unexpectedly succeeded: {value}");
        value
    }
    fn revision(&mut self) -> u64 {
        self.ok("document_info", json!({}))["revision"]
            .as_u64()
            .unwrap()
    }
}
impl Drop for Client {
    fn drop(&mut self) {
        let _ = self.child.kill();
        let _ = self.child.wait();
    }
}

fn manifest_entry(id: &str) -> Value {
    manifest()["entries"]
        .as_array()
        .unwrap()
        .iter()
        .find(|e| e["id"] == id)
        .unwrap_or_else(|| panic!("{id} is not in the manifest"))
        .clone()
}

/// The real executable over stdio MCP imports representative admitted
/// fixtures from every producer, saves, reopens after the original file is
/// deleted, re-exports with a readback receipt, and an independent reader
/// confirms the exported mass and bounds. A refused fixture in the same
/// session leaves the revision, the saved project and the source file
/// unchanged.
#[test]
fn rust_stdio_corpus_workflow_survives_source_removal() {
    for id in [
        "occt-placed-box",
        "occt-box-void",
        "occt-torus",
        "occt-pocketed-plate",
        "kicad-vendored-L_Wuerth_HCM-7070",
        "kicad-vendored-Mounting_Wuerth_WA-SMST-3.3mm_H1mm_9774010960",
    ] {
        let entry = manifest_entry(id);
        let bytes = entry_bytes(&entry).unwrap();
        let expected = import(std::str::from_utf8(&bytes).unwrap()).unwrap();
        let volume = expected.volume_mm3().unwrap();
        let bounds = bounds_of(&expected);
        let directory = tempfile::tempdir().unwrap();
        let source = directory.path().join("source.step");
        std::fs::write(&source, &bytes).unwrap();
        let refused = manifest_entry("occt-sphere");
        let refused_bytes = entry_bytes(&refused).unwrap();
        let refused_path = directory.path().join("refused.step");
        std::fs::write(&refused_path, &refused_bytes).unwrap();
        let project = directory.path().join("corpus.tcd");
        let arguments = [
            "--rust-model",
            "--project",
            project.to_str().unwrap(),
            "--export-dir",
            directory.path().to_str().unwrap(),
        ];
        let mut client = Client::start(&arguments);
        let imported = client.ok(
            "import_step",
            json!({"expected_revision": 0, "file_name": "source.step"}),
        );
        let [part] = imported["parts"].as_array().unwrap().as_slice() else {
            panic!("{id}: {imported}");
        };
        let model = client.ok("evaluate_model", json!({}));
        let evaluated = &model["parts"][0];
        assert_eq!(evaluated["part"], *part);
        assert!(
            close(
                evaluated["volume_mm3"].as_f64().unwrap(),
                volume,
                VOLUME_RELATIVE
            ),
            "{id}"
        );
        let min: [f64; 3] = serde_json::from_value(evaluated["bounds_mm"]["min"].clone()).unwrap();
        let max: [f64; 3] = serde_json::from_value(evaluated["bounds_mm"]["max"].clone()).unwrap();
        assert!(same_bounds([min, max], bounds, BOUNDS_MM), "{id}");
        client.ok("save_project", json!({"expected_revision": 1}));
        let saved = std::fs::read(&project).unwrap();
        // A refused import changes nothing.
        client.fail(
            "import_step",
            json!({"expected_revision": 1, "file_name": "refused.step"}),
        );
        assert_eq!(client.revision(), 1, "{id}");
        assert_eq!(std::fs::read(&project).unwrap(), saved, "{id}");
        assert_eq!(std::fs::read(&refused_path).unwrap(), refused_bytes, "{id}");
        drop(client);

        std::fs::remove_file(&source).unwrap();
        let mut reopened = Client::start(&arguments);
        assert_eq!(
            reopened.ok("evaluate_model", json!({}))["parts"],
            model["parts"],
            "{id}"
        );
        let receipt = reopened.ok(
            "export_step",
            json!({"part": part, "file_name": "exported.step"}),
        );
        assert!(
            close(
                receipt["readback_volume_mm3"].as_f64().unwrap(),
                volume,
                VOLUME_RELATIVE
            ),
            "{id}: {receipt}"
        );
        let exported = std::fs::read_to_string(directory.path().join("exported.step")).unwrap();
        let readback = independent_readback(&exported).unwrap();
        assert!(
            close(
                readback["volume_mm3"].as_f64().unwrap(),
                volume,
                VOLUME_RELATIVE
            ),
            "{id}"
        );
        let readback_bounds: [[f64; 3]; 2] =
            serde_json::from_value(readback["bounds_mm"].clone()).unwrap();
        assert!(same_bounds(readback_bounds, bounds, BOUNDS_MM), "{id}");
        assert_eq!(
            std::fs::read(&project).unwrap(),
            saved,
            "{id}: export changed the project"
        );
    }
}
