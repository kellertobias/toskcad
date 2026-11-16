//! Public actual-file mesh plus independent Native literal embedding gate.
//! Normal public singleton import; imported STEP export remains closed.
use multicad_core::{document::EntityId, references::TopologyKind};
use multicad_kernel::{GeometryPolicy, MeshOptions, Solid, rust_kernel as pure};
use pure::axial_pipe_thread::authored::{AuthoredSection, AuthoredSpec, AuthoredThread};
use pure::step::npt_stock::{CurveOptions, CurvePlan, SurfacePlan, write_source_candidate};
use pure::step_read::{
    npt_stock::{
        LiteralOccurrencePlacement, Parsed, Support,
        full_head::{OwnedActualHead, transport::encoding},
        head::{self, HeadInput},
    },
    owned_insert::{CompoundDocument, DocumentOptions},
};
use std::sync::Arc;
fn original_clipped_file(quarter: bool) -> String {
    let (tpi, e1, l1, l2, od, max) = if quarter {
        (18., 0.49163, 0.2278, 0.4018, 0.540, 0.088)
    } else {
        (27., 0.37360, 0.1615, 0.2639, 0.405, 0.096)
    };
    let pitch = 25.4 / tpi;
    let k = 1. / 32.;
    let m = 3f64.sqrt();
    let slopes = 1. / (m - k) + 1. / (m + k);
    let depth = pitch / slopes - (0.033 + max) * pitch;
    let source = AuthoredThread::new(
        AuthoredSpec {
            major_radius_mm: (e1 * 25.4 - 2. * k * l1 * 25.4 + depth) / 2.,
            pitch_mm: pitch,
            length_mm: l2 * 25.4,
            radial_taper: k,
            section: AuthoredSection::NptFlat {
                depth_mm: depth,
                flat_width_mm: (pitch - depth * slopes) / 2.,
            },
            crest_radius_limit_mm: Some(od * 25.4 / 2.),
        },
        1e-7,
        1e-6,
    )
    .unwrap()
    .stock_clip_plan(Default::default(), None)
    .unwrap();
    let options = CurveOptions::default();
    let curves = CurvePlan::new(&source, options, None).unwrap();
    let plan = SurfacePlan::new(curves, options.into(), None).unwrap();
    write_source_candidate(&plan, options.maximum_work, options.maximum_bytes, None)
        .unwrap()
        .into_string()
}
fn physical_head_file(quarter: bool) -> String {
    let text = original_clipped_file(quarter);
    let doc = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let parsed = Parsed::parse(
        &doc,
        doc.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap();
    let cap = parsed
        .faces()
        .iter()
        .filter_map(|f| match f.support() {
            Support::Plane { height } => {
                Some((parsed.points()[height].coordinates()[2].value(), f.entity()))
            }
            _ => None,
        })
        .max_by(|a, b| a.0.total_cmp(&b.0))
        .unwrap()
        .1;
    let head = HeadInput::regular_hex(&parsed, cap, 18., 3.).unwrap();
    let proof = head::certify(&parsed, head, 1e-7, 1e-6, Default::default(), None).unwrap();
    pure::step::npt_stock::full_head::propose(&proof, Default::default(), None)
        .unwrap()
        .as_str()
        .to_owned()
}
fn source_file(quarter: bool, posed: bool) -> String {
    let local = Arc::new(
        OwnedActualHead::read(
            Arc::from(physical_head_file(quarter)),
            Default::default(),
            None,
        )
        .unwrap(),
    );
    let frame = LiteralOccurrencePlacement::new(
        [[1., 0., 0.], [0., -1., 0.], [0., 0., -1.]],
        [0., 0., local.head_stations_mm()[0]],
        Default::default(),
        None,
    )
    .unwrap();
    let frame = if posed {
        frame
            .transformed(
                [[0.6, 0., 0.8], [0., 1., 0.], [-0.8, 0., 0.6]],
                [20., -30., 11.],
                Default::default(),
                None,
            )
            .unwrap()
    } else {
        frame
    };
    encoding::propose_mapped(&local, &frame, Default::default(), None)
        .unwrap()
        .actual_bytes()
        .to_string()
}
#[test]
fn normal_public_singleton_import_retains_original_file_and_placement_and_fresh_manufacturing() {
    for quarter in [false, true] {
        for posed in [false, true] {
            let input = source_file(quarter, posed);
            let feature = EntityId::from_uuid_bits(0xfeed);
            let source = pure::step_read::read_source_compound(&input, 1e-7).unwrap();
            assert!(matches!(
                source.as_slice(),
                [pure::step_read::SourceImport::ActualMappedHead(_)]
            ));
            drop(source);
            let solid =
                Solid::import_rust_step(&input, feature, GeometryPolicy::default()).unwrap();
            assert_eq!(
                solid
                    .rust_actual_clipped_head_source()
                    .unwrap()
                    .actual_file(),
                input
            );
            assert!(solid.rust_actual_clipped_head_occurrence().is_none());
            assert!(!solid.supports_step_export());
            assert_eq!(
                solid
                    .rust_actual_clipped_head_runtime_evidence()
                    .unwrap()
                    .face_measure_builds,
                0
            );
            let topology = solid.topology_snapshot(feature, 0).unwrap();
            assert_eq!(
                topology
                    .elements
                    .iter()
                    .filter(|v| v.key.kind == TopologyKind::Face)
                    .count(),
                17
            );
            let options = MeshOptions {
                linear_deflection_mm: 0.05,
                angular_deflection_rad: 0.3,
                maximum_triangles: 500_000,
            };
            let m = solid.manufacturing_mesh_with_receipt(options).unwrap();
            assert!(m.edge_pairs > 0 && m.vertex_links > 0);
            let dir = tempfile::tempdir().unwrap();
            let path = dir.path().join("actual.step");
            std::fs::write(&path, &input).unwrap();
            let reopened = Solid::read_step(&path, GeometryPolicy::default()).unwrap();
            std::fs::remove_file(&path).unwrap();
            assert_eq!(
                reopened
                    .rust_actual_clipped_head_source()
                    .unwrap()
                    .actual_file(),
                input
            );
            assert!((reopened.volume_mm3().unwrap() - solid.volume_mm3().unwrap()).abs() < 1e-8);
            assert!(reopened.manufacturing_mesh(options).is_ok());
            eprintln!(
                "NORMAL_PUBLIC_HEAD_IMPORT quarter={quarter} posed={posed} facets={}",
                m.mesh.face_indices.len()
            );
        }
    }
}
