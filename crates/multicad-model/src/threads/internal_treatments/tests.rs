use super::*;
use multicad_core::{
    commands::Session,
    document::Document,
    model::{
        self, ThreadEndTreatment, ThreadEndTreatments, ThreadFit, ThreadParameters,
        ThreadStartTreatment, length,
    },
};
use multicad_kernel::rust_kernel::{Circle, Rational, cylinder::AxialCylinder};
#[test]
fn real_resolved_roots_and_native_cancellation_keep_refusal_precedence() {
    let session = Session::new(Document::new("Source routing controls").unwrap()).unwrap();
    let parameters = ThreadParameters {
        major_diameter: length(3.),
        pitch: length(0.5),
        length: length(4.),
        minimum_wall: length(0.5),
        fit: ThreadFit::Metric6gExternal,
        end_treatments: Some(ThreadEndTreatments {
            start: ThreadStartTreatment::Chamfer45 {
                radial_depth: length(0.5),
            },
            end: ThreadEndTreatment::CylindricalRelief {
                length: length(0.5),
                radial_depth: length(0.125),
            },
        }),
        ..Default::default()
    };
    let external = parameters
        .resolve(session.parameters(), ThreadMode::NewExternal)
        .unwrap();
    assert!(
        external
            .metric_fit
            .as_ref()
            .unwrap()
            .external_root_radius_mm
            .is_some()
    );
    let q = |n, d| Rational::new(n, d).unwrap();
    let stock = Solid::from_axial_cylinder(
        AxialCylinder::new(
            1,
            Circle::new([q(0, 1); 2], q(5, 2)).unwrap(),
            q(0, 1),
            q(4, 1),
        )
        .unwrap(),
        PlaneFrame::default(),
        GeometryPolicy::default(),
    )
    .unwrap();
    assert!(
        build(
            EntityId::from_uuid_bits(2),
            &external,
            q(4, 1),
            &stock,
            GeometryPolicy::default(),
            None
        )
        .err()
        .unwrap()
        .message
        .contains("rounded roots")
    );
    let internal_parameters = ThreadParameters {
        fit: ThreadFit::Metric6HInternal,
        ..parameters.clone()
    };
    let missing = model::ThreadDefinition {
        output: model::BodyOutput {
            body: EntityId::from_uuid_bits(4),
        },
        target: None,
        mode: ThreadMode::CutInternal,
        parameters: internal_parameters,
    };
    assert!(missing.validate().is_err());
    assert!(
        super::super::build_rust(
            EntityId::from_uuid_bits(2),
            &missing,
            session.parameters(),
            &BTreeMap::new(),
            GeometryPolicy::default(),
            None
        )
        .err()
        .unwrap()
        .message
        .contains("no target predecessor")
    );
    let internal = ThreadParameters {
        fit: ThreadFit::Metric6HInternal,
        ..parameters
    }
    .resolve(session.parameters(), ThreadMode::CutInternal)
    .unwrap();
    assert!(
        build(
            EntityId::from_uuid_bits(2),
            &internal,
            q(4, 1),
            &stock,
            GeometryPolicy::default(),
            Some(&|| true)
        )
        .err()
        .unwrap()
        .message
        .contains("superseded")
    );
    let success = build(
        EntityId::from_uuid_bits(2),
        &internal,
        q(4, 1),
        &stock,
        GeometryPolicy::default(),
        None,
    )
    .unwrap();
    assert!(success.rust_internal_thread_treatment_data().is_some());
    let unsupported_stock = Solid::box_solid([5.; 3], GeometryPolicy::default()).unwrap();
    assert!(
        build(
            EntityId::from_uuid_bits(2),
            &internal,
            q(4, 1),
            &unsupported_stock,
            GeometryPolicy::default(),
            None
        )
        .is_err()
    );
    // Canonical custom treatment validation already rejects cosmetic input;
    // the geometry helper must also refuse that genuinely resolved source.
    let cosmetic = model::ThreadParameters {
        representation: ThreadRepresentation::CosmeticEnvelope,
        end_treatments: None,
        ..Default::default()
    }
    .resolve(session.parameters(), ThreadMode::CutInternal)
    .unwrap();
    assert!(
        build(
            EntityId::from_uuid_bits(2),
            &cosmetic,
            q(4, 1),
            &stock,
            GeometryPolicy::default(),
            None
        )
        .is_err()
    );
}
