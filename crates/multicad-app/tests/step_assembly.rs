//! STEP product structure through the application import path (CAD-361).
//!
//! Assemblies import as one Part whose members are the placed occurrences of
//! shared definitions. These tests follow an assembly through the canonical
//! document: placements and names, save and source-free reopen, structured
//! re-export and its independent readback, stable references to occurrence
//! faces, and atomic refusals of malformed structure. The comparison of every
//! occurrence with OpenCascade 7.9.3's own transfer was removed with the
//! native provider; its last result (b0b2ea52) is recorded in
//! `docs/architecture/exchange-corpus-ledger.md` as historical evidence.
use multicad_core::{
    commands::{Actor, Channel, Session},
    document::{Document, EntityId},
    model::imports::{self, StepEncoding},
    project_file::{self, ContentHash, FileLimits, ProjectFile},
    references::{Origin, TopologyKind},
};
use multicad_kernel::{
    AssemblyOccurrence, GeometryErrorKind, GeometryPolicy, Solid, rust_kernel::step_read,
};
use multicad_model::{Compiler, GeometryBackend, imports::read_step_file_with_backend};
use std::path::{Path, PathBuf};

const TWO: &str = include_str!("fixtures/exchange/corpus/occt-two-bodies.step");
const BOARD: &str = include_str!("fixtures/exchange/corpus/kicad-board-with-part.step");

fn import(text: &str) -> Result<Solid, multicad_kernel::GeometryError> {
    Solid::import_rust_step(text, EntityId::new(), GeometryPolicy::default())
}

fn occurrences(solid: &Solid) -> Vec<AssemblyOccurrence> {
    solid.assembly_occurrences().expect("imported assembly")
}

fn swap(text: &str, from: &str, to: &str) -> String {
    assert!(text.contains(from), "{from}");
    text.replacen(from, to, 1)
}

fn with(text: &str, extra: &str) -> String {
    let end = text.rfind("ENDSEC;").unwrap();
    format!("{}{extra}\n{}", &text[..end], &text[end..])
}

/// The two-body assembly with a second box placed 90 degrees about z at
/// (100, 0, 0), sharing the first box's definition.
fn repeated() -> String {
    let text = with(
        TWO,
        "#20000=AXIS2_PLACEMENT_3D('',#20001,#20002,#20003);\
         #20001=CARTESIAN_POINT('',(100.,0.,0.));#20002=DIRECTION('',(0.,0.,1.));\
         #20003=DIRECTION('',(0.,1.,0.));\
         #20010=NEXT_ASSEMBLY_USAGE_OCCURRENCE('3','again','',#5,#31,$);\
         #20011=PRODUCT_DEFINITION_SHAPE('','',#20010);\
         #20012=ITEM_DEFINED_TRANSFORMATION('','',#11,#20000);\
         #20013=(REPRESENTATION_RELATIONSHIP('','',#36,#10)REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION(#20012)SHAPE_REPRESENTATION_RELATIONSHIP());\
         #20014=CONTEXT_DEPENDENT_SHAPE_REPRESENTATION(#20013,#20011);",
    );
    swap(
        &text,
        "#10 = SHAPE_REPRESENTATION('',(#11,#15,#19),#23);",
        "#10 = SHAPE_REPRESENTATION('',(#11,#15,#19,#20000),#23);",
    )
}

/// The cylinder nested in the box, which the root offsets by (5, 6, 7) and
/// which turns the cylinder 90 degrees about z.
fn nested() -> String {
    let text = swap(
        TWO,
        "NEXT_ASSEMBLY_USAGE_OCCURRENCE('2','','',#5,#380,$)",
        "NEXT_ASSEMBLY_USAGE_OCCURRENCE('2','','',#31,#380,$)",
    );
    let text = swap(&text, "('','',#385,#10)", "('','',#385,#36)");
    let text = swap(
        &text,
        "#36 = ADVANCED_BREP_SHAPE_REPRESENTATION('',(#11,#37),#367);",
        "#36 = ADVANCED_BREP_SHAPE_REPRESENTATION('',(#11,#37,#19),#367);",
    );
    let text = swap(
        &text,
        "#16 = CARTESIAN_POINT('',(0.,0.,0.));",
        "#16 = CARTESIAN_POINT('',(5.,6.,7.));",
    );
    swap(
        &text,
        "#22 = DIRECTION('',(1.,0.,-0.));",
        "#22 = DIRECTION('',(0.,1.,0.));",
    )
}

/// The cylinder placed by a rigid Cartesian transformation operator: 90
/// degrees about z at (100, 0, 0).
fn operated() -> String {
    with(
        &swap(
            TWO,
            "REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION(#491)",
            "REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION(#20030)",
        ),
        "#20031=DIRECTION('',(0.,1.,0.));#20032=DIRECTION('',(-1.,0.,0.));\
         #20033=DIRECTION('',(0.,0.,1.));#20034=CARTESIAN_POINT('',(100.,0.,0.));\
         #20030=CARTESIAN_TRANSFORMATION_OPERATOR_3D('','','',#20031,#20032,#20034,1.,#20033);",
    )
}

/// The installed KiCad 3D-model library, or `MULTICAD_KICAD_3DMODELS`.
fn library() -> PathBuf {
    std::env::var_os("MULTICAD_KICAD_3DMODELS")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            PathBuf::from("/Applications/KiCad/KiCad.app/Contents/SharedSupport/3dmodels")
        })
}

/// The MeanWell IRM-60 converter: an assembly whose component solid `#110`
/// has horn-torus corner blends (major = minor = 0.125 mm).
const CONVERTER: &str = "Converter_ACDC.3dshapes/Converter_ACDC_MeanWell_IRM-60-xx_THT.step";

fn close(a: f64, b: f64, relative: f64) -> bool {
    (a - b).abs() <= relative * b.abs().max(1.)
}

fn near(a: [f64; 3], b: [f64; 3]) -> bool {
    a.iter().zip(b).all(|(x, y)| (x - y).abs() <= 1e-9)
}

fn actor() -> Actor {
    Actor::local(Channel::Mcp)
}

#[test]
fn a_board_places_its_component_at_the_stated_occurrence() {
    let solid = import(BOARD).unwrap();
    let found = occurrences(&solid);
    assert_eq!(found.len(), 2);
    let component = &found[0];
    assert_eq!(component.path, ["board-with-part 1", "U1"]);
    assert_eq!(component.definition_name, "occt-cylinder");
    assert!(near(component.origin_mm, [130., -120., 1.615]));
    // KiCad turns U1 by 90 degrees: the component's x runs along y.
    assert!(near(component.axes[0], [1.110223024625E-16, 1., 0.]));
    let board = &found[1];
    assert_eq!(board.path, ["board-with-part 1", "=>[0:1:1:3]"]);
    assert_eq!(board.definition_name, "board-with-part_PCB");
    let total: f64 = found.iter().map(|o| o.volume_mm3).sum();
    assert!(close(total, solid.volume_mm3().unwrap(), 1e-12));
    assert_eq!(solid.solid_count().unwrap(), 2);
}

#[test]
fn repeated_components_share_one_definition_at_distinct_placements() {
    let solid = import(&repeated()).unwrap();
    let found = occurrences(&solid);
    assert_eq!(found.len(), 3);
    assert_eq!(found[0].definition, found[2].definition);
    assert_ne!(found[0].definition, found[1].definition);
    assert_eq!(found[2].path[1], "again");
    assert!(near(found[2].origin_mm, [100., 0., 0.]));
    assert!(close(found[2].volume_mm3, found[0].volume_mm3, 1e-12));
    assert!(near(found[2].bounds.min, [90., 0., 0.]));
    assert!(near(found[2].bounds.max, [100., 10., 5.]));
    // Each occurrence has its own references; the shared definition gives
    // corresponding faces the same role.
    let topology = solid
        .topology_with_display(EntityId::new(), 0, None)
        .unwrap()
        .topology;
    topology.validate().unwrap();
    let faces: Vec<_> = topology
        .elements
        .iter()
        .filter(|e| e.key.kind == TopologyKind::Face)
        .collect();
    let first: Vec<&Origin> = faces[..6].iter().flat_map(|e| &e.origins).collect();
    let again: Vec<&Origin> = faces[faces.len() - 6..]
        .iter()
        .flat_map(|e| &e.origins)
        .collect();
    for (a, b) in first.iter().zip(&again) {
        assert_eq!(a.role, b.role);
        assert_ne!(a.source, b.source);
    }
}

#[test]
fn nested_and_operator_placements_compose() {
    let found = occurrences(&import(&nested()).unwrap());
    assert_eq!(found[1].path.len(), 3);
    assert!(near(found[1].origin_mm, [5., 6., 7.]));
    assert!(near(found[1].axes[0], [0., 1., 0.]));
    let found = occurrences(&import(&operated()).unwrap());
    assert!(near(found[1].origin_mm, [100., 0., 0.]));
    assert!(near(found[1].axes[1], [-1., 0., 0.]));
}

/// Commit an import of `text` and return the session with the import, part
/// and feature ids.
fn imported(text: &str, directory: &Path) -> (Session, EntityId, EntityId, EntityId) {
    let source = directory.join("assembly.step");
    std::fs::write(&source, text).unwrap();
    let file = read_step_file_with_backend(&source, GeometryBackend::Rust).unwrap();
    assert_eq!(file.encoding, StepEncoding::RustStep);
    let mut session = Session::new(Document::new("Assembly").unwrap()).unwrap();
    let prepared = imports::prepare_step_import(&session, 0, actor(), &file).unwrap();
    let [output] = prepared.model.outputs.as_slice() else {
        panic!("one Part per STEP source");
    };
    let (part, feature) = (output.body, output.feature);
    let import = prepared.source;
    session.commit(prepared.model.prepared).unwrap();
    std::fs::remove_file(&source).unwrap();
    (session, import, part, feature)
}

fn compiled(document: &Document, feature: EntityId) -> Solid {
    let model = Compiler {
        omit_display: true,
        ..Default::default()
    }
    .compile_with_backend(document, 0, GeometryBackend::Rust, None)
    .unwrap();
    model.solid(feature).unwrap().share().unwrap()
}

fn face_origins(solid: &Solid) -> Vec<Vec<Origin>> {
    let topology = solid
        .topology_with_display(EntityId::new(), 0, None)
        .unwrap()
        .topology;
    topology
        .elements
        .into_iter()
        .filter(|e| e.key.kind == TopologyKind::Face)
        .map(|e| e.origins.into_iter().collect())
        .collect()
}

#[test]
fn an_assembly_survives_reopen_without_its_source_and_reexports_its_structure() {
    let directory = tempfile::tempdir().unwrap();
    let text = repeated();
    let (session, source, _, feature) = imported(&text, directory.path());
    let before = compiled(session.document(), feature);
    let project = directory.path().join("assembly.tcd");
    let mut encoded = vec![];
    project_file::write(
        &ProjectFile::new(session.document().clone()).unwrap(),
        &mut encoded,
        FileLimits::default(),
    )
    .unwrap();
    std::fs::write(&project, &encoded).unwrap();
    let opened = project_file::read(
        &mut std::fs::read(&project).unwrap().as_slice(),
        FileLimits::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap();
    let document = opened.project.document;
    let stored = imports::read_step_source(&document, source).unwrap();
    assert_eq!(stored.content, ContentHash::of(text.as_bytes()));
    let after = compiled(&document, feature);
    assert_eq!(occurrences(&after), occurrences(&before));
    // References to occurrence faces resolve to the same faces after reopen.
    assert_eq!(face_origins(&after), face_origins(&before));

    // The re-export keeps one definition per product and one occurrence per
    // instance, and reads back to the same placed solids.
    let exported = directory.path().join("export.step");
    after.write_step(&exported).unwrap();
    let written = std::fs::read_to_string(&exported).unwrap();
    assert_eq!(
        written.matches("NEXT_ASSEMBLY_USAGE_OCCURRENCE(").count(),
        3
    );
    assert_eq!(written.matches("MANIFOLD_SOLID_BREP(").count(), 2);
    let tolerance = GeometryPolicy::default().linear_tolerance_mm;
    let structure = step_read::brep::read_assembly(&written, tolerance, None).unwrap();
    assert_eq!(structure.nodes.len(), 3);
    assert_eq!(structure.occurrences().len(), 3);
    let volume = after.verify_imported_brep_step(&written).unwrap().unwrap();
    assert!(close(volume, after.volume_mm3().unwrap(), 1e-9));
    let back = import(&written).unwrap();
    let (a, b) = (occurrences(&back), occurrences(&after));
    assert_eq!(a.len(), b.len());
    for (a, b) in a.iter().zip(&b) {
        assert_eq!((&a.path, a.definition), (&b.path, b.definition));
        assert!(close(a.volume_mm3, b.volume_mm3, 1e-9));
        assert!(near(a.origin_mm, b.origin_mm) && near(a.bounds.min, b.bounds.min));
    }
    // Reimporting the same source with the same identity keeps references.
    let again = Solid::import_rust_step(&text, feature, GeometryPolicy::default()).unwrap();
    assert_eq!(face_origins(&again), face_origins(&after));
}

#[test]
fn malformed_structure_refuses_atomically_with_its_kind() {
    let cases = [
        // A product placed inside itself.
        (
            with(
                TWO,
                "#20000=NEXT_ASSEMBLY_USAGE_OCCURRENCE('3','','',#31,#31,$);\
                 #20001=PRODUCT_DEFINITION_SHAPE('','',#20000);\
                 #20002=ITEM_DEFINED_TRANSFORMATION('','',#11,#11);\
                 #20003=(REPRESENTATION_RELATIONSHIP('','',#36,#36)REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION(#20002)SHAPE_REPRESENTATION_RELATIONSHIP());\
                 #20004=CONTEXT_DEPENDENT_SHAPE_REPRESENTATION(#20003,#20001);",
            ),
            GeometryErrorKind::InvalidResult,
        ),
        // A dangling component.
        (
            swap(TWO, "('1','','',#5,#31,$)", "('1','','',#5,#99999,$)"),
            GeometryErrorKind::InvalidResult,
        ),
        // A scaled and a mirroring placement.
        (
            swap(&operated(), ",#20034,1.,", ",#20034,2.,"),
            GeometryErrorKind::Unsupported,
        ),
        (
            swap(
                &operated(),
                "#20032=DIRECTION('',(-1.,0.,0.))",
                "#20032=DIRECTION('',(1.,0.,0.))",
            ),
            GeometryErrorKind::Unsupported,
        ),
        // A component in metres beside millimetres.
        (
            swap(
                TWO,
                "#485 = ( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.) );",
                "#485 = ( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT($,.METRE.) );",
            ),
            GeometryErrorKind::Unsupported,
        ),
    ];
    let directory = tempfile::tempdir().unwrap();
    let (session, ..) = imported(&repeated(), directory.path());
    let document = session.document().clone();
    let revision = session.revision();
    for (index, (text, kind)) in cases.iter().enumerate() {
        assert_eq!(import(text).err().unwrap().kind, *kind, "case {index}");
        let path = directory.path().join(format!("refused-{index}.step"));
        std::fs::write(&path, text).unwrap();
        assert!(read_step_file_with_backend(&path, GeometryBackend::Rust).is_err());
        assert_eq!(session.revision(), revision);
        assert_eq!(session.document(), &document);
    }
}

/// Every placed solid of the converter imports on its own terms: exact
/// volume and bounds per occurrence, the horn-torus component among them,
/// and a certified closed manufacturing mesh of each solid whose enclosed
/// volume agrees with the exact one.
#[test]
fn every_converter_solid_meshes_closed_including_the_horn_torus_blends() {
    use multicad_kernel::rust_kernel::{brep, manufacturing};
    let path = library().join(CONVERTER);
    if !path.is_file() {
        eprintln!("{} is not installed; skipped", path.display());
        return;
    }
    let text = std::fs::read_to_string(&path).unwrap();
    let tolerance = GeometryPolicy::default().linear_tolerance_mm;
    let assembly = step_read::brep::read_assembly(&text, tolerance, None).unwrap();
    let occurrences = assembly.occurrences();
    assert_eq!(occurrences.len(), 6, "OpenCascade transfers six solids");
    let mut rows = vec![];
    let mut horn = 0;
    let mut total = 0.;
    for (index, occurrence) in occurrences.iter().enumerate() {
        let solid = assembly.place(occurrence).unwrap();
        let horns = solid
            .boundary()
            .faces
            .iter()
            .filter(|f| matches!(f.surface, brep::Surface::Torus { major, minor, .. } if major == minor))
            .count();
        horn += horns;
        let options = brep::mesh::Options {
            deflection_mm: 0.05,
            angle_rad: 0.3,
            triangle_limit: 5_000_000,
        };
        let mesh = brep::mesh::tessellate(&solid, options, None)
            .unwrap_or_else(|e| panic!("solid {index}: mesh: {e:?}"));
        assert!(
            mesh.maximum_deviation_mm <= 0.05,
            "solid {index}: deviation {}",
            mesh.maximum_deviation_mm
        );
        let proof = manufacturing::certify(&mesh.display())
            .unwrap_or_else(|e| panic!("solid {index}: certification: {e:?}"));
        let exact = solid.volume_mm3();
        let relative = (proof.signed_volume_mm3 - exact).abs() / exact;
        assert!(relative <= 0.02, "solid {index}: mesh volume {relative}");
        total += exact;
        rows.push(serde_json::json!({
            "solid": index,
            "node": assembly.nodes[occurrence.node].name,
            "volume_mm3": exact,
            "volume_error_mm3": solid.volume_error_mm3(),
            "bounds_mm": solid.bounds_mm(),
            "faces": solid.faces().len(),
            "horn_torus_faces": horns,
            "triangles": mesh.triangles.len(),
            "edge_pairs": proof.edge_pairs,
            "mesh_deviation_mm": mesh.maximum_deviation_mm,
            "mesh_volume_relative": relative,
        }));
    }
    assert_eq!(horn, 4, "four horn-torus corner blends");
    assert!(close(total, 133_479.602_204_707, 1e-9), "total {total}");
    println!("{}", serde_json::to_string_pretty(&rows).unwrap());
}
