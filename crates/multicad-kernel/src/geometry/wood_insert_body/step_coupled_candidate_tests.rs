//! Candidate-only Native source integration. Independent file material/mass
//! is not source correspondence and does not enable the normal export route.
use super::*;
use multicad_rust_kernel::{
    authored_flat_thread::{AuthoredThread, SourceOptions},
    step::insert::{
        CurveOptions,
        coupled::{CoupledPlan, topology::CoupledTopology},
    },
    step_read::owned_insert::{CompoundDocument, DocumentOptions, FluxOptions, Options},
};
fn actual_source(length: f64, left: bool) -> source::Insert {
    let original = crate::geometry::wood_insert_world_tests::source(left);
    let quality = SourceOptions {
        position_limit_mm: 1e-7,
        normal_angle_limit_rad: 1e-7,
    };
    let with_length = |law: &AuthoredThread| {
        let mut spec = law.spec();
        spec.length_mm = length;
        AuthoredThread::new(spec, quality).unwrap()
    };
    source::Insert::new(
        with_length(original.outer()),
        with_length(original.inner()),
        original.minimum_wall_mm(),
        None,
    )
    .unwrap()
}
#[test]
fn native_actual_tiny_trim_sources_produce_mapped_material_candidates_without_admission() {
    let posed = PlaneFrame {
        origin: [40., -18., 7.],
        u: [0.8, 0.6, 0.],
        v: [-0.36, 0.48, 0.8],
    };
    for (length, left, frame) in [
        (4_f64.next_up(), false, PlaneFrame::default()),
        (4_f64.next_up(), true, posed),
        (4_f64.next_down(), false, posed),
        (8., true, PlaneFrame::default()),
    ] {
        let source = actual_source(length, left);
        let policy = GeometryPolicy::default();
        let solid = Solid::rust_wood_insert(source, EntityId::new(), frame, policy, None).unwrap();
        let (body, actual_frame) = solid.wood_insert_body().unwrap();
        assert_eq!(actual_frame, frame);
        assert_eq!(body.source, source);
        let plan = CoupledPlan::prepare(
            body.boundary(),
            CurveOptions {
                position_error_mm: policy.linear_tolerance_mm * 0.25,
                ..Default::default()
            },
            None,
        )
        .unwrap();
        let topology = CoupledTopology::prepare(&plan, 16_000_000, 64 * 1024 * 1024, None).unwrap();
        let text = topology
            .write_mapped_candidate(
                frame.origin,
                frame.u,
                frame.v,
                16_000_000,
                64 * 1024 * 1024,
                None,
            )
            .unwrap();
        let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
        let [member] = document.represented_members() else {
            panic!("one actual mapped candidate required");
        };
        let packet = document
            .member_packet(member.item_entity(), Options::default(), None)
            .unwrap();
        let parsed = packet.parsed();
        let mapping = parsed
            .physical_mapping(
                member.mapping_link().unwrap().operator_entity(),
                Options::default(),
                None,
            )
            .unwrap();
        assert_eq!(mapping.handedness().unwrap(), 1);
        let unit = plan.coefficient_unit_mm();
        let scale = mapping.scale_bounds();
        assert!(scale[0] <= unit && unit <= scale[1]);
        assert!(scale[1] < 1e-20);
        let binding = parsed.support_bindings(Options::default(), None).unwrap();
        let material = binding.material_bounds(Options::default(), None).unwrap();
        let determinant = mapping.volume_scale_bounds().unwrap();
        let width = (B::point(1e-6) / B::point(determinant[1])).lo;
        let raw = material
            .occupied_volume(
                FluxOptions {
                    maximum_width_mm3: width,
                    ..Default::default()
                },
                None,
            )
            .unwrap()
            .volume_mm3();
        let decoded = B::new(raw[0], raw[1]) * B::new(determinant[0], determinant[1]);
        assert!(decoded.checked().is_some() && decoded.lo > 0.);
        assert!((B::point(decoded.hi) - B::point(decoded.lo)).hi <= 1e-6);
        // The source query and decoded boundary integral are independent.
        // The candidate is an approximation, so interval equality is not expected.
        // This regression comparison is NOT a shape-correspondence certificate.
        let native = WoodInsertMaterial::new(source, frame, policy, None)
            .unwrap()
            .volume_bounds(
                source::QueryOptions {
                    volume_error_mm3: 1e-7,
                    maximum_cells: 65536,
                },
                None,
            )
            .unwrap();
        let difference = decoded - B::new(native[0], native[1]);
        assert!(
            difference.magnitude() <= 1e-6,
            "length={length:?}, left={left}, decoded={decoded:?}, source={native:?}"
        );
        let mapped = mapping.map_box_mm(material.bounds_mm()).unwrap();
        assert!(
            mapped
                .into_iter()
                .all(|b| b[0].is_finite() && b[1].is_finite() && b[0] < b[1])
        );
        assert_eq!(solid.rust_wood_insert_data().unwrap().source, source);
        assert_eq!(source.full_length_mm().to_bits(), length.to_bits());
    }
}

#[test]
fn actual_tiny_trim_mixed_export_checks_mapping_source_and_atomic_cancellation() {
    let policy = GeometryPolicy::default();
    let frame = PlaneFrame {
        origin: [40., -18., 7.],
        u: [0.8, 0.6, 0.],
        v: [-0.36, 0.48, 0.8],
    };
    let source = actual_source(4_f64.next_up(), true);
    let insert = Solid::rust_wood_insert(source, EntityId::new(), frame, policy, None).unwrap();
    let sphere = Solid::rust_sphere_named(
        1.,
        PlaneFrame {
            origin: [20., 0., 0.],
            ..Default::default()
        },
        EntityId::new(),
        policy,
    )
    .unwrap();
    let mixed = Solid::compound_named(&[&sphere, &insert], EntityId::new(), policy).unwrap();
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("actual-mixed-tiny-insert.step");
    mixed.write_step(&path).unwrap();
    let text = std::fs::read_to_string(&path).unwrap();
    let proof = mixed
        .verify_wood_insert_step_compound(&text, 1e-7, None)
        .unwrap()
        .unwrap();
    assert_eq!(proof.members.solid_count, 2);
    assert_eq!(proof.members.retained_members.len(), 1);
    assert_eq!(proof.members.decoded_scalar_members.len(), 1);
    assert_eq!(proof.insert_correspondence[0].0, 1);
    assert!(proof.insert_correspondence[0].1.maximum_position_error_mm <= 1e-7);
    let cancelled = mixed
        .verify_wood_insert_step_compound(&text, 1e-7, Some(&|| true))
        .unwrap()
        .err()
        .unwrap();
    assert_eq!(cancelled.kind, GeometryErrorKind::Cancelled);
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let member = document.represented_members()[1];
    let op = member.mapping_link().unwrap().operator_entity();
    let prefix = format!("#{op}=");
    let record = text.lines().find(|l| l.starts_with(&prefix)).unwrap();
    let attributes = record
        .split("CARTESIAN_TRANSFORMATION_OPERATOR(")
        .nth(1)
        .unwrap()
        .split(')')
        .next()
        .unwrap();
    let scale = attributes.rsplit(',').next().unwrap();
    let changed_record = record.replacen(&format!(",{scale})"), ",1.)", 1);
    assert_ne!(changed_record, record);
    let changed = text.replacen(record, &changed_record, 1);
    assert!(
        mixed
            .verify_wood_insert_step_compound(&changed, 1e-7, None)
            .unwrap()
            .is_err()
    );
    // Matching local geometry and scale do not authorize another occurrence.
    let origin_id: u32 = attributes
        .split(',')
        .nth(2)
        .unwrap()
        .trim()
        .trim_start_matches('#')
        .parse()
        .unwrap();
    let origin_prefix = format!("#{origin_id}=");
    let origin_record = text
        .lines()
        .find(|l| l.starts_with(&origin_prefix))
        .unwrap();
    let shifted_record = format!("#{origin_id}=CARTESIAN_POINT('',(41.,-18.,7.));");
    let shifted = text.replacen(origin_record, &shifted_record, 1);
    let error = mixed
        .verify_wood_insert_step_compound(&shifted, 1e-7, None)
        .unwrap()
        .err()
        .unwrap();
    assert_eq!(error.kind, GeometryErrorKind::InvalidResult);
    assert_eq!(std::fs::read_to_string(path).unwrap(), text);
    assert_eq!(insert.rust_wood_insert_data().unwrap().source, source);
}
