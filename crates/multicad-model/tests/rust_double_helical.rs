//! Gapped double-helical gears through the normal Rust model provider.
use multicad_core::{
    commands::{Actor, Channel, Session},
    document::Document,
    model::{self, GearFamily, GearHand, GearHub, GearParameters, GearRequest, angle, length},
    project_file::{self, FileLimits, ProjectFile},
    references::{Origin, TopologyKind},
};
use multicad_model::{Compiler, GeometryBackend};
use multicad_precision::PortableMath as _;

fn actor() -> Actor {
    Actor::local(Channel::Ui)
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected} ±{tolerance}"
    );
}

fn compiler() -> Compiler {
    Compiler {
        omit_display: true,
        ..Default::default()
    }
}

fn request() -> GearRequest {
    GearRequest::new(
        "Double helical",
        GearParameters {
            family: GearFamily::DoubleHelical,
            hand: GearHand::Right,
            helix_angle: angle(20.),
            face_width: length(12.),
            center_gap: length(3.),
            bore_diameter: length(8.),
            hub: Some(GearHub {
                diameter: length(16.),
                axial_start: length(-3.),
                axial_end: length(5.),
            }),
            ..Default::default()
        },
    )
}

#[test]
fn rust_gap_edits_history_reopen_and_refusals_keep_the_actual_core() {
    let mut session = Session::new(Document::new("Double helical").unwrap()).unwrap();
    let request = request();
    let prepared = model::prepare_gear(&session, 0, actor(), &request).unwrap();
    session.commit(prepared.prepared).unwrap();
    let compile = |document: &Document, revision| {
        compiler()
            .compile_with_backend(document, revision, GeometryBackend::Rust, None)
            .unwrap()
    };
    let first = compile(session.document(), 1);
    let solid = first.solid(request.feature).unwrap();
    assert!(solid.uses_rust_geometry());
    let volume = solid.volume_mm3().unwrap();
    let shape = &first.snapshot().shapes[&request.feature];
    let role = |role: &str| {
        shape
            .topology
            .elements
            .iter()
            .filter(|element| {
                element.key.kind == TopologyKind::Face
                    && element.origins.contains(&Origin {
                        source: request.feature,
                        role: role.into(),
                    })
            })
            .count()
    };
    assert_eq!(role("gear/root.core"), 1);
    let teeth = role("gear/gap.lower");
    assert!(teeth > 4);
    assert_eq!(role("gear/gap.upper"), teeth);
    assert_eq!(role("gear/hub.side"), 1);

    // Widening the gap by 1 mm replaces 1 mm of tooth section with root
    // annulus; the difference is the tooth area outside the root circle.
    let resolved = request
        .parameters
        .resolve_geometry(session.parameters())
        .unwrap();
    let model::ResolvedGearGeometry::Rotary(gear) = resolved else {
        panic!("Expected a rotary gear");
    };
    assert_eq!(gear.center_gap_mm, 3.);
    let shoulders = solid
        .faces()
        .unwrap()
        .iter()
        .filter(|face| {
            face.plane.is_some()
                && (face.center_mm[2] - (12. - 3.) / 2.).abs() < 1e-9
                && face.center_mm[0].portable_hypot(face.center_mm[1]) > gear.root_radius_mm
        })
        .map(|face| face.area_mm2)
        .sum::<f64>();
    let widened = model::GearParameters {
        center_gap: length(4.),
        ..request.parameters.clone()
    };
    let prepared =
        model::prepare_gear_edit(&session, 1, actor(), request.feature, &widened).unwrap();
    let candidate = compile(prepared.prepared.document(), 2);
    let changed = candidate
        .solid(request.feature)
        .unwrap()
        .volume_mm3()
        .unwrap();
    near(volume - changed, shoulders, 1e-6);
    session.commit(prepared.prepared).unwrap();
    let edited = session.document().clone();
    session.undo(2, actor()).unwrap();
    near(
        compile(session.document(), 3)
            .solid(request.feature)
            .unwrap()
            .volume_mm3()
            .unwrap(),
        volume,
        1e-9,
    );
    session.redo(3, actor()).unwrap();
    assert_eq!(session.document(), &edited);

    // A gap that leaves no tooth band is refused without changing the model.
    let before = session.document().clone();
    let invalid = model::GearParameters {
        center_gap: length(12.),
        ..request.parameters.clone()
    };
    if let Ok(prepared) = model::prepare_gear_edit(
        &session,
        session.revision(),
        actor(),
        request.feature,
        &invalid,
    ) {
        assert!(
            compiler()
                .compile_with_backend(prepared.prepared.document(), 5, GeometryBackend::Rust, None)
                .is_err()
        );
    }
    assert_eq!(session.document(), &before);

    let mut bytes = Vec::new();
    project_file::write(
        &ProjectFile::new(session.document().clone()).unwrap(),
        &mut bytes,
        FileLimits::default(),
    )
    .unwrap();
    let reopened = project_file::read(
        &mut bytes.as_slice(),
        FileLimits::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap();
    let readback = compile(&reopened.project.document, 0);
    near(
        readback
            .solid(request.feature)
            .unwrap()
            .volume_mm3()
            .unwrap(),
        changed,
        1e-9,
    );
    assert_eq!(
        readback.snapshot().shapes[&request.feature]
            .topology
            .elements,
        candidate.snapshot().shapes[&request.feature]
            .topology
            .elements
    );

    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("double-helical.step");
    let solid = readback.solid(request.feature).unwrap();
    assert!(solid.supports_step_export());
    solid.write_step(&path).unwrap();
    near(
        multicad_model::export::verify_rust_step_shape(
            &path,
            solid,
            &readback.snapshot().shapes[&request.feature],
            compiler().policy.linear_tolerance_mm,
        )
        .unwrap(),
        changed,
        1e-5,
    );
    if multicad_kernel::OCCT_LINKED {
        let native = compiler()
            .compile_with_backend(
                &reopened.project.document,
                0,
                GeometryBackend::OpenCascade,
                None,
            )
            .unwrap();
        let native = native.solid(request.feature).unwrap();
        near(native.volume_mm3().unwrap(), changed, 2e-4);
        let [a, b] = [native.bounds().unwrap(), solid.bounds().unwrap()];
        for (x, y) in a
            .min
            .into_iter()
            .chain(a.max)
            .zip(b.min.into_iter().chain(b.max))
        {
            near(x, y, 1e-5);
        }
        let imported = Solid::read_step(&path, Default::default()).unwrap();
        near(imported.volume_mm3().unwrap(), changed, 2e-4);
    }
}

use multicad_kernel::Solid;
