use multicad_core::{
    commands::{Actor, Channel, Session},
    document::{Document, EntityId},
    model::{
        self, CustomConnectionHardwareParameters, MechanicalPartParameters, MechanicalPartRequest,
        PrimitiveFrame, ResolvedCustomConnectionHardware, ResolvedMechanicalPart, ThreadHand,
        length,
    },
    project_file::{self, ProjectFile},
};
use multicad_kernel::{GeometryPolicy, Solid};
use multicad_model::{Compiler, GeometryBackend};
use std::f64::consts::PI;

#[test]
fn custom_authored_wood_threads_keep_their_exact_start_frame_and_physical_geometry() {
    let compiler = Compiler {
        omit_display: true,
        ..Default::default()
    };
    for (left, plain, posed) in [
        (false, 0., false),
        (true, 0., true),
        (false, 2.3, true),
        (true, 2.3, false),
    ] {
        let mut session =
            Session::new(Document::new("Actual authored wood thread").unwrap()).unwrap();
        let mut parameters = model::default_custom_connection_wood_screw();
        let MechanicalPartParameters::CustomConnection { parameters: custom } = &mut parameters
        else {
            unreachable!()
        };
        let CustomConnectionHardwareParameters::WoodScrew {
            wood_thread,
            plain_length,
            frame,
            ..
        } = custom.as_mut()
        else {
            unreachable!()
        };
        wood_thread.hand = if left {
            ThreadHand::Left
        } else {
            ThreadHand::Right
        };
        *plain_length = length(plain);
        if posed {
            *frame = PrimitiveFrame {
                origin: [length(1000.3), length(-2000.7), length(0.1)],
                u: [0., 1., 0.],
                v: [0., 0., 1.],
            };
        }
        let ResolvedMechanicalPart::CustomConnection {
            parameters: resolved,
        } = parameters.resolve(session.parameters()).unwrap()
        else {
            unreachable!()
        };
        let ResolvedCustomConnectionHardware::WoodScrew {
            wood_thread: source,
            ..
        } = resolved.as_ref()
        else {
            unreachable!()
        };
        let request = MechanicalPartRequest::new("Explicit custom profile", parameters);
        session
            .commit(
                model::prepare_mechanical_part(&session, 0, Actor::local(Channel::Ui), &request)
                    .unwrap()
                    .prepared,
            )
            .unwrap();
        let built = compiler
            .compile_with_backend(
                session.document(),
                session.revision(),
                GeometryBackend::Rust,
                None,
            )
            .unwrap();
        let solid = built.solid(request.feature).unwrap();
        let (actual, frame) = solid.rust_machine_screw_geometry().unwrap();
        assert_eq!(
            frame.origin.map(f64::to_bits),
            source.frame.origin.map(f64::to_bits)
        );
        assert_eq!(frame.u, source.frame.u);
        assert_eq!(frame.v, source.frame.v);
        assert_eq!(actual.thread().spec().left_hand, left);
        let r = source.major_radius_mm;
        let root = r - source.profile.depth_mm;
        let slopes = source.profile.lower_flank_deg.to_radians().tan()
            + source.profile.upper_flank_deg.to_radians().tan();
        let crest = source.profile.crest_width_mm;
        let root_width = source.pitch_mm - crest - (r - root) * slopes;
        let thread_volume = PI
            * source.length_mm
            * (crest * r * r + slopes * (r.powi(3) - root.powi(3)) / 3. + root_width * root * root)
            / source.pitch_mm;
        let expected = thread_volume + PI * r * r * plain + 3_f64.sqrt() * 50. * 4.;
        assert!((solid.volume_mm3().unwrap() - expected).abs() < expected * 1e-9);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("actual-profile.step");
        let mesh = solid
            .manufacturing_mesh_with_receipt(Default::default())
            .unwrap();
        assert_eq!(mesh.edge_pairs as usize * 2, mesh.mesh.indices.len());
        solid.write_step(&path).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        let (physical, _) =
            multicad_kernel::rust_kernel::step_read::read_import_solid(&text, 1e-7).unwrap();
        assert_eq!(physical.solid_count, 1);
        assert!((physical.volume_mm3 - expected).abs() < expected * 1e-8);
        let decoded =
            Solid::import_rust_step(&text, EntityId::new(), GeometryPolicy::default()).unwrap();
        let mesh = decoded
            .manufacturing_mesh_with_receipt(Default::default())
            .unwrap();
        assert_eq!(mesh.edge_pairs as usize * 2, mesh.mesh.indices.len());
        let mut saved = Vec::new();
        project_file::write(
            &ProjectFile::new(session.document().clone()).unwrap(),
            &mut saved,
            Default::default(),
        )
        .unwrap();
        let reopened = project_file::read(
            &mut saved.as_slice(),
            Default::default(),
            &project_file::supported_capabilities(),
        )
        .unwrap();
        let fresh = compiler
            .compile_with_backend(&reopened.project.document, 0, GeometryBackend::Rust, None)
            .unwrap();
        assert_eq!(
            fresh
                .solid(request.feature)
                .unwrap()
                .rust_machine_screw_geometry()
                .unwrap()
                .1,
            frame
        );
        assert!(
            model::manufacturing_limitation_in_history(session.document(), request.feature)
                .unwrap()
                .is_none()
        );
    }
}

#[test]
fn unsupported_custom_hanger_keeps_the_previous_document() {
    let compiler = Compiler {
        omit_display: true,
        ..Default::default()
    };
    let session = Session::new(Document::new("Atomic custom hardware refusal").unwrap()).unwrap();
    let before = session.document().clone();
    {
        let source = model::default_custom_connection_hanger_bolt();
        let request = MechanicalPartRequest::new("Needs connected material", source);
        let candidate =
            model::prepare_mechanical_part(&session, 0, Actor::local(Channel::Ui), &request)
                .unwrap();
        let failure = compiler
            .compile_with_backend(
                candidate.prepared.document(),
                1,
                GeometryBackend::Rust,
                None,
            )
            .err()
            .unwrap();
        assert!(
            failure.message.contains("hanger bolts require their own connected-material provider"),
            "{failure:?}"
        );
        assert_eq!(session.document(), &before);
        assert_eq!(session.revision(), 0);
    }
}
