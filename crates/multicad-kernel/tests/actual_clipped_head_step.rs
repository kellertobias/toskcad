//! Public actual-file mesh plus independent Native literal embedding gate.
//! Normal public singleton import; imported STEP export remains closed.
use multicad_core::document::EntityId;
use multicad_kernel::{
    GeometryErrorKind, GeometryPolicy, MeshOptions, RigidTransform, Solid, rust_kernel as pure,
};
use pure::axial_pipe_thread::authored::{AuthoredSection, AuthoredSpec, AuthoredThread};
use pure::step::npt_stock::{CurveOptions, CurvePlan, SurfacePlan, write_source_candidate};
use pure::step_read::npt_stock::full_head::transport::replay::ReplayOptions;
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
fn imported_actual_head_normal_step_reexport_retains_whole_local_geometry_and_b_then_a() {
    for quarter in [false, true] {
        for posed in [false, true] {
            for extra in [false, true] {
                let input = source_file(quarter, posed);
                let base =
                    Solid::import_rust_step(&input, EntityId::new(), GeometryPolicy::default())
                        .unwrap();
                let solid = if extra {
                    base.transformed(RigidTransform {
                        rotation: [[0.6, 0., 0.8], [0., 1., 0.], [-0.8, 0., 0.6]],
                        translation_mm: [-13., 21., 4.],
                    })
                    .unwrap()
                } else {
                    base
                };
                assert!(solid.supports_step_export());
                let local = solid
                    .rust_actual_clipped_head_source()
                    .unwrap()
                    .local_geometry()
                    .actual_file()
                    .to_owned();
                let dir = tempfile::tempdir().unwrap();
                let path = dir.path().join("actual-output.step");
                let receipt = solid
                    .write_actual_clipped_head_step_with_verification(
                        &path,
                        Default::default(),
                        None,
                    )
                    .unwrap()
                    .unwrap();
                assert!(
                    receipt.aggregate_work <= 272_000_000
                        && receipt.accounted_peak_bytes <= 64 * 1024 * 1024
                );
                assert!(receipt.maximum_world_position_error_mm <= 1e-7);
                assert!(
                    receipt.candidate_preparation_work > 0
                        && receipt.actual_output_binding_work > 0
                        && receipt.native_and_io_work > 0
                );
                let actual: Arc<str> = Arc::from(std::fs::read_to_string(&path).unwrap());
                let checked = solid
                    .verify_actual_clipped_head_step(actual.clone(), Default::default(), None)
                    .unwrap()
                    .unwrap();
                assert!(checked.maximum_world_position_error_mm <= 1e-7);
                let cold = Solid::read_step(&path, GeometryPolicy::default()).unwrap();
                std::fs::remove_file(&path).unwrap();
                assert_eq!(
                    cold.rust_actual_clipped_head_source()
                        .unwrap()
                        .local_geometry()
                        .actual_file(),
                    local
                );
                let owner = EntityId::new();
                let vertices = |s: &Solid| {
                    s.topology_snapshot(owner, 0)
                        .unwrap()
                        .elements
                        .into_iter()
                        .filter_map(|e| match e.geometry {
                            Some(multicad_core::references::AnalyticGeometry::Point {
                                position_mm,
                            }) => Some((e.key.ordinal, position_mm)),
                            _ => None,
                        })
                        .collect::<Vec<_>>()
                };
                let before = vertices(&solid);
                let after = vertices(&cold);
                assert_eq!(before.len(), 24);
                assert_eq!(after.len(), 24);
                for ((i, a), (j, b)) in before.into_iter().zip(after) {
                    assert_eq!(i, j);
                    for (a, b) in a.into_iter().zip(b) {
                        assert!((a - b).abs() < 1e-7);
                    }
                }
                assert!((solid.volume_mm3().unwrap() - cold.volume_mm3().unwrap()).abs() < 1e-7);
                assert_eq!(
                    cold.rust_actual_clipped_head_runtime_evidence()
                        .unwrap()
                        .face_measure_builds,
                    0
                );
                assert!(
                    cold.rust_actual_clipped_head_occurrence().is_none(),
                    "A is encoded in new actual B exactly once"
                );
                cold.manufacturing_mesh(MeshOptions {
                    linear_deflection_mm: 0.05,
                    angular_deflection_rad: 0.3,
                    maximum_triangles: 500_000,
                })
                .unwrap();
                // The ordinary writer routes this provider too; no private export API is required.
                solid.write_step(&path).unwrap();
                assert_eq!(std::fs::read_to_string(&path).unwrap(), actual.as_ref());
                eprintln!(
                    "NATIVE_ACTUAL_HEAD_REEXPORT quarter={quarter} Bpose={posed} A={extra} work={} peak={} C0={}",
                    receipt.aggregate_work,
                    receipt.accounted_peak_bytes,
                    receipt.maximum_world_position_error_mm
                );
            }
        }
    }
}
#[test]
fn actual_head_step_limits_and_last_cancel_preserve_existing_output_and_source() {
    use std::cell::Cell;
    let input = source_file(true, true);
    let solid =
        Solid::import_rust_step(&input, EntityId::new(), GeometryPolicy::default()).unwrap();
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("keep.step");
    let n = Cell::new(0usize);
    let observe = || {
        n.set(n.get() + 1);
        false
    };
    let receipt = solid
        .write_actual_clipped_head_step_with_verification(&path, Default::default(), Some(&observe))
        .unwrap()
        .unwrap();
    let total = n.get();
    assert!(total > 2);
    let exact = ReplayOptions {
        resources: pure::step_read::npt_stock::full_head::Options {
            maximum_work: receipt.aggregate_work,
            maximum_bytes: receipt.accounted_peak_bytes,
            ..Default::default()
        },
        ..Default::default()
    };
    solid
        .write_actual_clipped_head_step_with_verification(&path, exact, None)
        .unwrap()
        .unwrap();
    for (w, b) in [
        (receipt.aggregate_work - 1, receipt.accounted_peak_bytes),
        (receipt.aggregate_work, receipt.accounted_peak_bytes - 1),
        (1, receipt.accounted_peak_bytes),
    ] {
        std::fs::write(&path, "prior output").unwrap();
        let mut low = exact;
        low.resources.maximum_work = w;
        low.resources.maximum_bytes = b;
        let error = solid
            .write_actual_clipped_head_step_with_verification(&path, low, None)
            .unwrap()
            .unwrap_err();
        assert_eq!(error.kind, GeometryErrorKind::Limit);
        assert_eq!(std::fs::read_to_string(&path).unwrap(), "prior output");
    }
    let mut tiny = exact;
    tiny.maximum_position_error_mm = 1e-30;
    assert_eq!(
        solid
            .write_actual_clipped_head_step_with_verification(&path, tiny, None)
            .unwrap()
            .unwrap_err()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(std::fs::read_to_string(&path).unwrap(), "prior output");
    for late in [false, true] {
        let count = Cell::new(0usize);
        let stop = || {
            count.set(count.get() + 1);
            !late || count.get() >= total
        };
        assert_eq!(
            solid
                .write_actual_clipped_head_step_with_verification(
                    &path,
                    Default::default(),
                    Some(&stop)
                )
                .unwrap()
                .unwrap_err()
                .kind,
            GeometryErrorKind::Cancelled
        );
        assert_eq!(std::fs::read_to_string(&path).unwrap(), "prior output");
    }
    assert_eq!(
        solid
            .rust_actual_clipped_head_source()
            .unwrap()
            .actual_file(),
        input
    );
    assert_eq!(
        std::fs::read_dir(dir.path()).unwrap().count(),
        1,
        "temporary artifacts cleaned after every refusal"
    );
    solid
        .write_actual_clipped_head_step_with_verification(&path, exact, None)
        .unwrap()
        .unwrap();
    assert!(Solid::read_step(&path, GeometryPolicy::default()).is_ok());
}
