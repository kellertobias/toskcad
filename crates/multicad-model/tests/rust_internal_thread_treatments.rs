//! Actual canonical combined receiving cuts; no retained native or display proxy.
use multicad_core::{
    commands::{Actor, Channel, Session},
    document::Document,
    model::{
        self, BodyInput, PrimitiveFrame, PrimitiveParameters, PrimitiveRequest, ThreadEndTreatment,
        ThreadEndTreatments, ThreadFit, ThreadHand, ThreadMode, ThreadParameters, ThreadRequest,
        ThreadStartTreatment, length,
    },
    project_file::{self, ProjectFile},
    references::{self, GeometryIntent, Resolution, TopologyKind},
};
use multicad_kernel::PlaneFrame;
use multicad_model::{Compiler, GeometryBackend};
fn actor() -> Actor {
    Actor::local(Channel::Ui)
}
fn compiler() -> Compiler {
    Compiler {
        omit_display: true,
        ..Default::default()
    }
}
fn ends(chamfer: f64, width: f64) -> ThreadEndTreatments {
    ThreadEndTreatments {
        start: ThreadStartTreatment::Chamfer45 {
            radial_depth: length(chamfer),
        },
        end: ThreadEndTreatment::CylindricalRelief {
            length: length(width),
            radial_depth: length(0.125),
        },
    }
}
fn fixture(hand: ThreadHand, starts: u32, posed: bool, height: f64) -> (Session, ThreadRequest) {
    let mut session = Session::new(Document::new("Combined internal treatments").unwrap()).unwrap();
    let frame = if posed {
        PrimitiveFrame {
            origin: [length(10.), length(20.), length(30.)],
            u: [0.8, 0.6, 0.],
            v: [-0.36, 0.48, 0.8],
        }
    } else {
        Default::default()
    };
    let blank = PrimitiveRequest::new(
        "Cylinder stock",
        PrimitiveParameters::Cylinder {
            radius: length(2.5),
            height: length(height),
            frame: frame.clone(),
        },
    );
    let p = model::prepare_primitive(&session, 0, actor(), &blank).unwrap();
    session.commit(p.prepared).unwrap();
    let thread = ThreadRequest::receiving(
        "Treated M3 internal",
        ThreadMode::CutInternal,
        BodyInput {
            body: blank.body,
            producer: blank.feature,
        },
        ThreadParameters {
            major_diameter: length(3.),
            pitch: length(0.5),
            length: length(4.),
            fit: ThreadFit::Metric6HInternal,
            hand,
            starts,
            minimum_wall: length(0.5),
            end_treatments: Some(ends(0.5, 0.5)),
            frame,
            ..Default::default()
        },
    )
    .unwrap();
    (session, thread)
}
#[test]
fn combined_cylindrical_internal_treatments_preserve_actual_source_stock_and_material() {
    for (hand, starts, posed) in [
        (ThreadHand::Right, 1, false),
        (ThreadHand::Left, 1, false),
        (ThreadHand::Right, 1, true),
        (ThreadHand::Left, 1, true),
    ] {
        let (session, thread) = fixture(hand, starts, posed, 4.);
        let stock = compiler()
            .compile_with_backend(session.document(), 1, GeometryBackend::Rust, None)
            .unwrap();
        let input = thread.target.as_ref().unwrap();
        let original = &stock.snapshot().shapes[&input.producer].topology;
        let intents: Vec<_> = original
            .elements
            .iter()
            .filter(|e| e.key.kind == TopologyKind::Face)
            .map(|e| GeometryIntent::capture(original, e.key).unwrap())
            .collect();
        let before = session.document().clone();
        let candidate = model::prepare_thread(&session, 1, actor(), &thread).unwrap();
        let compiled = Compiler::default()
            .compile_with_backend(
                candidate.prepared.document(),
                2,
                GeometryBackend::Rust,
                None,
            )
            .unwrap();
        let solid = compiled.solid(thread.feature).unwrap();
        let data = solid.rust_internal_thread_treatment_data().unwrap();
        let resolved = thread
            .parameters
            .resolve(session.parameters(), ThreadMode::CutInternal)
            .unwrap();
        let spec = data.material.cavity().source_thread().spec();
        assert_eq!(
            (spec.major_radius_mm, spec.pitch_mm, spec.length_mm),
            (resolved.major_radius_mm, resolved.pitch_mm, 4.)
        );
        assert_eq!(
            (
                spec.depth_mm,
                spec.crest_width_mm,
                spec.lower_flank_deg,
                spec.upper_flank_deg
            ),
            (
                resolved.profile.depth_mm,
                resolved.profile.crest_width_mm,
                resolved.profile.lower_flank_deg,
                resolved.profile.upper_flank_deg
            )
        );
        assert_eq!(
            (spec.starts, spec.left_hand),
            (starts, hand == ThreadHand::Left)
        );
        assert_eq!(
            data.frame,
            PlaneFrame {
                origin: resolved.frame.origin,
                u: resolved.frame.u,
                v: resolved.frame.v
            }
        );
        assert_eq!(data.material.stock().feature(), input.producer.uuid_bits());
        assert!(data.material.minimum_wall_mm() >= resolved.minimum_wall_mm);
        assert!(data.material.minimum_wall_mm() - resolved.minimum_wall_mm < 1e-10);
        assert_eq!(data.material.cavity().intent().entry_depth_mm, Some(0.5));
        assert_eq!(data.material.cavity().terminal_station_mm(), Some(3.5));
        let cavity = data.material.cavity();
        let (legacy, transfer) = cavity.source_thread().legacy_approximation();
        assert!(
            transfer.maximum_position_error_mm
                <= multicad_kernel::GeometryPolicy::default().linear_tolerance_mm / 16.
        );
        for theta in [0., 0.43, 2.1] {
            assert!(
                (cavity.radius_mm(theta, 1.7) - legacy.radius(theta, 1.7)).abs()
                    <= transfer.maximum_position_error_mm + 1e-14
            );
        }

        for theta in [0., 0.43, 2.1] {
            assert!(
                (cavity.radius_mm(theta, 0.01)
                    - (resolved.major_radius_mm - resolved.profile.depth_mm + 0.49))
                    .abs()
                    < 1e-12
            );
            assert!(
                (cavity.radius_mm(theta, 3.9) - (resolved.major_radius_mm + 0.125)).abs() < 1e-12
            );
        }

        for z in [0.01, 1.7, 3.9] {
            for theta in [0., 0.43, 2.1] {
                let r = cavity.radius_mm(theta, z);
                for (d, inside) in [(-0.001, true), (0.001, false)] {
                    assert_eq!(
                        cavity.contains_void([(r + d) * theta.cos(), (r + d) * theta.sin(), z], 0.),
                        inside
                    );
                }
            }
        }
        let snapshot = &compiled.snapshot().shapes[&thread.feature];
        assert!(!snapshot.mesh.indices.is_empty());
        assert!(!snapshot.edges.is_empty() && !snapshot.vertices.is_empty());
        assert_eq!(snapshot.topology.revision, 2);
        for intent in intents {
            assert!(matches!(
                references::resolve(&intent, &snapshot.topology, 2),
                Resolution::Resolved { .. }
            ));
        }
        let receipt = solid
            .manufacturing_mesh_with_receipt(Default::default())
            .unwrap();
        assert_eq!(receipt.edge_pairs as usize * 2, receipt.mesh.indices.len());
        assert!(
            receipt
                .mesh
                .face_indices
                .iter()
                .all(|&f| f < data.material.face_count())
        );
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("not-qualified.step");
        std::fs::write(&path, b"keep existing destination").unwrap();
        assert!(solid.write_step(&path).is_err());
        assert_eq!(std::fs::read(&path).unwrap(), b"keep existing destination");
        assert_eq!(session.document(), &before);
    }
}
#[test]
fn canonical_edit_history_and_reopen_keep_the_same_authored_source() {
    let (mut session, request) = fixture(ThreadHand::Left, 1, false, 4.);
    let p = model::prepare_thread(&session, 1, actor(), &request).unwrap();
    session.commit(p.prepared).unwrap();
    let original = session.document().clone();
    let mut edited = request.parameters.clone();
    edited.end_treatments = Some(ends(0.625, 0.625));
    let p = model::prepare_thread_edit(&session, 2, actor(), request.feature, &edited).unwrap();
    let compiled = compiler()
        .compile_with_backend(p.prepared.document(), 3, GeometryBackend::Rust, None)
        .unwrap();
    let actual = compiled
        .solid(request.feature)
        .unwrap()
        .rust_internal_thread_treatment_data()
        .unwrap();
    assert_eq!(actual.material.cavity().full_length_mm(), 4.);
    assert_eq!(
        actual.material.cavity().intent().entry_depth_mm,
        Some(0.625)
    );
    assert_eq!(actual.material.cavity().terminal_station_mm(), Some(3.375));
    session.commit(p.prepared).unwrap();
    let saved = session.document().clone();
    session.undo(3, actor()).unwrap();
    assert_eq!(session.document(), &original);
    session.redo(4, actor()).unwrap();
    assert_eq!(session.document(), &saved);
    let mut bytes = vec![];
    project_file::write(
        &ProjectFile::new(saved).unwrap(),
        &mut bytes,
        Default::default(),
    )
    .unwrap();
    let reopened = project_file::read(
        &mut bytes.as_slice(),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap()
    .project;
    assert_eq!(
        model::read_thread(&reopened.document, request.feature)
            .unwrap()
            .parameters,
        edited
    );
    let rebuilt = compiler()
        .compile_with_backend(&reopened.document, 0, GeometryBackend::Rust, None)
        .unwrap();
    assert_eq!(
        rebuilt
            .solid(request.feature)
            .unwrap()
            .rust_internal_thread_treatment_data()
            .unwrap(),
        actual
    );
}
#[test]
fn partial_shifted_undersized_and_unresolved_sources_refuse_without_committing() {
    for case in 0..6 {
        let (session, mut request) =
            fixture(ThreadHand::Right, 1, false, if case == 0 { 5. } else { 4. });
        match case {
            1 => request.parameters.frame.origin[0] = length(0.01),
            2 => request.parameters.minimum_wall = length(1.),
            3 => request.parameters.end_treatments = Some(ends(0.1, 0.5)),
            4 => {
                request.parameters.representation =
                    multicad_core::model::ThreadRepresentation::CosmeticEnvelope
            }
            5 => {
                request.parameters.frame.u = [0., 1., 0.];
                request.parameters.frame.v = [-1., 0., 0.];
            }
            _ => {}
        }
        let before = session.document().clone();
        let Ok(p) = model::prepare_thread(&session, 1, actor(), &request) else {
            assert_eq!(session.document(), &before);
            continue;
        };
        assert!(
            compiler()
                .compile_with_backend(p.prepared.document(), 2, GeometryBackend::Rust, None)
                .is_err(),
            "case {case}"
        );
        assert_eq!(session.document(), &before);
    }
    let (session, request) = fixture(ThreadHand::Right, 1, false, 4.);
    let before = session.document().clone();
    let p = model::prepare_thread(&session, 1, actor(), &request).unwrap();
    assert!(
        compiler()
            .compile_with_backend(
                p.prepared.document(),
                2,
                GeometryBackend::Rust,
                Some(&|| true)
            )
            .is_err()
    );
    assert_eq!(session.document(), &before);
}
#[test]
fn nominal_basic_source_default_limit_remains_a_real_diagnostic() {
    use multicad_kernel::{
        GeometryErrorKind, GeometryPolicy, HelicalThread, InternalThreadTreatmentCut,
        rust_kernel::authored_internal_treatment as source,
    };
    let (session, mut request) = fixture(ThreadHand::Right, 1, false, 4.);
    let six = request
        .parameters
        .resolve(session.parameters(), ThreadMode::CutInternal)
        .unwrap();
    eprintln!(
        "Paired actual6H R={} p={} L={} depth={} crest={} lower={} upper={} starts={} hand={:?} frame={:?}",
        six.major_radius_mm,
        six.pitch_mm,
        six.length_mm,
        six.profile.depth_mm,
        six.profile.crest_width_mm,
        six.profile.lower_flank_deg,
        six.profile.upper_flank_deg,
        six.starts,
        six.hand,
        six.frame
    );
    request.parameters.fit = ThreadFit::NominalBasic;
    let r = request
        .parameters
        .resolve(session.parameters(), ThreadMode::CutInternal)
        .unwrap();
    let profile = HelicalThread {
        major_radius_mm: r.major_radius_mm,
        pitch_mm: r.pitch_mm,
        length_mm: r.length_mm,
        depth_mm: r.profile.depth_mm,
        crest_width_mm: r.profile.crest_width_mm,
        lower_flank_deg: r.profile.lower_flank_deg,
        upper_flank_deg: r.profile.upper_flank_deg,
        radial_taper: r.profile.radial_taper,
        starts: r.starts,
        left_hand: false,
    };
    let frame = PlaneFrame {
        origin: r.frame.origin,
        u: r.frame.u,
        v: r.frame.v,
    };
    eprintln!(
        "NominalBasic actual profile={profile:?} frame={frame:?} display={:?}",
        Compiler::default().mesh
    );
    let model = compiler()
        .compile_with_backend(session.document(), 1, GeometryBackend::Rust, None)
        .unwrap();
    let stock = model
        .solid(request.target.as_ref().unwrap().producer)
        .unwrap();
    let native = stock.rust_cut_authored_internal_thread_treatments_with_interrupt(
        InternalThreadTreatmentCut {
            profile,
            frame,
            minimum_wall_mm: r.minimum_wall_mm,
            intent: source::Intent {
                entry_depth_mm: Some(0.5),
                terminal: Some(source::Terminal {
                    width_mm: 0.5,
                    radial_depth_mm: 0.125,
                }),
            },
        },
        request.feature,
        GeometryPolicy::default(),
        None,
    );
    let native = match native {
        Ok(s) => s,
        Err(e) => {
            eprintln!("NominalBasic native construction: {e:?}");
            assert_eq!(e.kind, GeometryErrorKind::Limit);
            return;
        }
    };
    let mut limits = 0;
    macro_rules! stage {
        ($name:literal,$call:expr) => {
            match $call {
                Ok(_) => eprintln!(concat!("NominalBasic ", $name, ": PASS")),
                Err(e) => {
                    eprintln!(concat!("NominalBasic ", $name, ": {:?}"), e);
                    assert_eq!(e.kind, GeometryErrorKind::Limit);
                    limits += 1;
                }
            }
        };
    }
    stage!("volume", native.volume_mm3());
    stage!("faces", native.faces());
    stage!("topology", native.topology_snapshot(request.feature, 2));
    stage!(
        "display edges",
        native.display_edges(Compiler::default().mesh.linear_deflection_mm)
    );
    stage!("display vertices", native.display_vertices());
    stage!("display mesh", native.mesh(Compiler::default().mesh));
    stage!(
        "manufacturing",
        native.manufacturing_mesh_with_receipt(Default::default())
    );
    // Keep this real failure control until the source owner independently
    // qualifies the same default geometry; no quality/work-limit waiver.
    assert!(limits > 0);
    let before = session.document().clone();
    let candidate = model::prepare_thread(&session, 1, actor(), &request).unwrap();
    assert!(
        compiler()
            .compile_with_backend(
                candidate.prepared.document(),
                2,
                GeometryBackend::Rust,
                None
            )
            .is_err()
    );
    assert_eq!(session.document(), &before);
}
