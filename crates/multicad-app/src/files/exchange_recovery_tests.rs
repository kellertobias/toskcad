//! Exchange recovery controls (CAD-345): a project holding independently
//! produced STEP imports and pinned library definitions survives recovery,
//! damaged recovery copies and failed saves, and rebuilds without the
//! original foreign files.
use super::*;
use multicad_core::{
    commands::{Actor, Channel, Session},
    document::{Entity, EntityData, Mutation},
    materials::{self, Assignment, DiagnosticKind, HatchStyle, LibraryDefinition},
    model::imports::{self, StepEncoding},
};
use multicad_kernel::Solid;
use multicad_model::{Compiler, GeometryBackend, imports::read_step_file_with_backend};

const FREECAD: &[u8] = include_bytes!("../../tests/fixtures/exchange/freecad-drilled-box.step");
const OPENCASCADE: &[u8] =
    include_bytes!("../../tests/fixtures/exchange/corpus/occt-box-void.step");
const KICAD_LIBRARY: &[u8] =
    include_bytes!("../../tests/fixtures/exchange/corpus/kicad-library/L_Wuerth_HCM-7070.step");

fn actor() -> Actor {
    Actor::local(Channel::Ui)
}

struct Imported {
    session: Session,
    features: Vec<EntityId>,
    parts: Vec<EntityId>,
    hatch: Entity,
}

/// Import three foreign producers and assign a pinned material and hatch.
fn imported(directory: &Path) -> Imported {
    let mut session = Session::new(Document::new("Foreign exchange recovery").unwrap()).unwrap();
    let mut features = vec![];
    let mut parts = vec![];
    for (name, bytes) in [
        ("freecad.step", FREECAD),
        ("opencascade.step", OPENCASCADE),
        ("kicad-library.step", KICAD_LIBRARY),
    ] {
        let path = directory.join(name);
        std::fs::write(&path, bytes).unwrap();
        let file = read_step_file_with_backend(&path, GeometryBackend::Rust).unwrap();
        assert_eq!(file.encoding, StepEncoding::RustStep);
        let prepared =
            imports::prepare_step_import(&session, session.revision(), actor(), &file).unwrap();
        features.extend(prepared.model.outputs.iter().map(|o| o.feature));
        parts.extend(prepared.model.outputs.iter().map(|o| o.body));
        session.commit(prepared.model.prepared).unwrap();
    }
    let material = Entity::new(
        "Recovery steel",
        EntityData::Material {
            density_kg_m3: Some(7850.),
            color_rgba: [80, 80, 90, 255],
            section_hatch: None,
        },
    );
    let hatch = LibraryDefinition::Hatch(HatchStyle::default())
        .entity(EntityId::new(), "Recovery hatch")
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
    Imported {
        session,
        features,
        parts,
        hatch,
    }
}

fn rebuild(document: &Document, features: &[EntityId]) -> Vec<(u32, f64, [[f64; 3]; 2])> {
    let model = Compiler {
        omit_display: true,
        ..Default::default()
    }
    .compile_with_backend(document, 0, GeometryBackend::Rust, None)
    .unwrap();
    features
        .iter()
        .map(|&id| {
            let solid: &Solid = model.solid(id).unwrap();
            let bounds = solid.bounds().unwrap();
            (
                solid.solid_count().unwrap(),
                solid.volume_mm3().unwrap(),
                [bounds.min, bounds.max],
            )
        })
        .collect()
}

#[test]
fn foreign_step_recovery_copies_rebuild_without_sources_and_never_replace_the_original() {
    let sources = tempfile::tempdir().unwrap();
    let imported = imported(sources.path());
    let before = rebuild(imported.session.document(), &imported.features);
    assert_eq!(before.len(), 3);
    assert!(
        before
            .iter()
            .all(|(count, volume, _)| *count == 1 && *volume > 0.)
    );

    let scratch = tempfile::tempdir().unwrap();
    let path = scratch.path().join("foreign.tcd");
    let mut store = ProjectStore::new(scratch.path().join("recovery"));
    let project = ProjectFile::new(imported.session.document().clone()).unwrap();
    let saved = store
        .save_project(&path, &project, SaveMode::Create)
        .unwrap();
    let original = std::fs::read(&path).unwrap();
    // The foreign files are gone; only the embedded sources remain.
    sources.close().unwrap();

    let entry = store
        .write_recovery(&project, imported.session.revision(), Some(&path))
        .unwrap();
    let restored = store.restore_recovery(entry.info.id).unwrap();
    assert!(restored.stamp.is_none());
    assert_eq!(restored.file.document, project.document);
    assert_eq!(rebuild(&restored.file.document, &imported.features), before);
    for part in &imported.parts {
        let effective = materials::resolve(&restored.file.document, *part).unwrap();
        assert!(effective.material.is_some());
        assert_eq!(effective.hatch, Some(HatchStyle::default()));
        assert!(
            effective.diagnostics.is_empty(),
            "{:?}",
            effective.diagnostics
        );
    }
    for source in imports::step_sources(&restored.file.document) {
        let stored = imports::read_step_source(&restored.file.document, source.id).unwrap();
        assert!(
            [FREECAD, OPENCASCADE, KICAD_LIBRARY]
                .iter()
                .any(|bytes| ContentHash::of(bytes) == stored.content),
            "embedded source differs from the foreign file"
        );
    }
    assert_eq!(std::fs::read(&path).unwrap(), original);

    // A damaged recovery copy is refused, listed with a diagnostic, and
    // leaves the original project byte-for-byte unchanged.
    let recovery = store.recovery_path(entry.info.id);
    let mut damaged = std::fs::read(&recovery).unwrap();
    let middle = damaged.len() / 2;
    damaged[middle] ^= 0x20;
    std::fs::write(&recovery, &damaged).unwrap();
    // Listing reads only the checked manifest summary; restoring verifies
    // every chunk and refuses the damaged source payload.
    assert!(store.restore_recovery(entry.info.id).is_err());
    let mut manifest_damage = std::fs::read(&recovery).unwrap();
    manifest_damage[60] ^= 0x20;
    std::fs::write(&recovery, &manifest_damage).unwrap();
    let listing = store.recoveries().unwrap();
    assert!(
        listing.entries.is_empty(),
        "a damaged recovery manifest was listed"
    );
    assert!(
        !listing.diagnostics.is_empty(),
        "a damaged recovery copy was hidden"
    );
    assert!(store.restore_recovery(entry.info.id).is_err());
    std::fs::write(&recovery, &damaged[..damaged.len() / 3]).unwrap();
    assert!(store.restore_recovery(entry.info.id).is_err());
    assert_eq!(
        std::fs::read(&recovery).unwrap(),
        &damaged[..damaged.len() / 3]
    );
    assert_eq!(std::fs::read(&path).unwrap(), original);

    // Failed saves keep the last valid project: a stale replace stamp and a
    // non-project extension are both refused before writing.
    let mut edited = project.clone();
    edited.document = edited
        .document
        .edited(&[Mutation::RenameDocument("Edited elsewhere".into())])
        .unwrap();
    assert!(
        store
            .save_project(&path, &edited, SaveMode::Replace(saved.stamp))
            .is_ok()
    );
    let current = std::fs::read(&path).unwrap();
    assert!(
        store
            .save_project(&path, &project, SaveMode::Replace(saved.stamp))
            .is_err(),
        "a stale stamp replaced newer bytes"
    );
    assert!(
        store
            .save_project(
                &scratch.path().join("foreign.step"),
                &project,
                SaveMode::Create
            )
            .is_err()
    );
    assert_eq!(std::fs::read(&path).unwrap(), current);
    let reopened = store.open_project(&path).unwrap();
    assert_eq!(rebuild(&reopened.file.document, &imported.features), before);
}

#[test]
fn pinned_library_definitions_survive_recovery_and_report_changed_content() {
    let sources = tempfile::tempdir().unwrap();
    let mut imported = imported(sources.path());
    let scratch = tempfile::tempdir().unwrap();
    let mut store = ProjectStore::new(scratch.path().join("recovery"));
    let project = ProjectFile::new(imported.session.document().clone()).unwrap();
    let entry = store
        .write_recovery(&project, imported.session.revision(), None)
        .unwrap();

    // Change the pinned hatch definition in the live session.
    let mut changed = imported.hatch.clone();
    changed.data = LibraryDefinition::Hatch(HatchStyle {
        angle_deg: 15.,
        ..Default::default()
    })
    .entity(changed.id, changed.name.clone())
    .unwrap()
    .data;
    imported
        .session
        .apply(multicad_core::commands::Request {
            expected_revision: imported.session.revision(),
            actor: actor(),
            label: "Edit pinned hatch".into(),
            edits: vec![Mutation::Replace(changed)],
        })
        .unwrap();
    for part in &imported.parts {
        let effective = materials::resolve(imported.session.document(), *part).unwrap();
        assert!(
            effective
                .diagnostics
                .iter()
                .any(|d| d.kind == DiagnosticKind::DefinitionChanged
                    && d.target == imported.hatch.id),
            "a changed pinned definition was accepted silently"
        );
    }
    // The recovery copy still holds the pinned original definition.
    let restored = store.restore_recovery(entry.info.id).unwrap();
    assert_eq!(restored.file.document, project.document);
    for part in &imported.parts {
        let effective = materials::resolve(&restored.file.document, *part).unwrap();
        assert!(effective.diagnostics.is_empty());
        assert_eq!(effective.hatch, Some(HatchStyle::default()));
    }
}
