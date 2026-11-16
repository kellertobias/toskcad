//! Public actual-file mesh plus independent Native literal embedding gate.
//! This does not enable a Solid provider or a global STEP capability.
use multicad_kernel::{Mesh, rust_kernel as pure};
use pure::axial_pipe_thread::authored::{AuthoredSection, AuthoredSpec, AuthoredThread};
use pure::step::npt_stock::{CurveOptions, CurvePlan, SurfacePlan, write_source_candidate};
use pure::step_read::{
    npt_stock::{
        LiteralOccurrencePlacement, Parsed, Support,
        full_head::{
            OwnedActualHead,
            mesh::HeadMeshOptions,
            transport::{MappedActualHead, encoding},
        },
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
#[test]
fn actual_core_head_strict_world_facets_earn_fresh_embedding_within_original_remaining_balance() {
    for quarter in [false, true] {
        let local = Arc::new(
            OwnedActualHead::read(
                Arc::from(physical_head_file(quarter)),
                Default::default(),
                None,
            )
            .unwrap(),
        );
        let intrinsic = LiteralOccurrencePlacement::new(
            [[1., 0., 0.], [0., -1., 0.], [0., 0., -1.]],
            [0., 0., local.head_stations_mm()[0]],
            Default::default(),
            None,
        )
        .unwrap();
        let encoded =
            encoding::propose_mapped(&local, &intrinsic, Default::default(), None).unwrap();
        let actual = MappedActualHead::read_reusing_local(
            encoded.actual_bytes(),
            local,
            Default::default(),
            None,
        )
        .unwrap();
        let occurrence = LiteralOccurrencePlacement::new(
            [[0.6, -0.8, 0.], [0.8, 0.6, 0.], [0., 0., 1.]],
            [20., -30., 11.],
            Default::default(),
            None,
        )
        .unwrap();
        for placed in [false, true] {
            let options = HeadMeshOptions {
                tessellation: pure::display::Tessellation {
                    deflection: 0.05,
                    angle: 0.3,
                    triangle_limit: 500_000,
                },
                ..Default::default()
            };
            let candidate = if placed {
                actual.placed_mesh_candidate(&occurrence, options, None)
            } else {
                actual.mesh_candidate(options, None)
            }
            .unwrap();
            let receipt = candidate.receipt();
            let raw = candidate.into_mesh();
            let mesh = Mesh {
                vertices: raw.vertices,
                indices: raw.indices,
                face_indices: raw.face_indices,
            };
            let proof = mesh
                .prove_literal_embedding(receipt.remaining_work, None)
                .unwrap();
            assert!(std::ptr::eq(proof.mesh(), &mesh));
            assert_eq!(
                receipt.aggregate_work + proof.work_used() + proof.remaining_work(),
                1_600_000_000
            );
            eprintln!(
                "ACTUAL_HEAD_NATIVE_EMBEDDING quarter={quarter} placed={placed} facets={} inherited={} fresh={} pairs={} remaining={}",
                mesh.face_indices.len(),
                receipt.aggregate_work,
                proof.work_used(),
                proof.pairs(),
                proof.remaining_work()
            );
            let used = proof.work_used();
            assert!(mesh.prove_literal_embedding(used - 1, None).is_err());
            assert_eq!(
                mesh.prove_literal_embedding(used, None)
                    .unwrap()
                    .remaining_work(),
                0
            );
            let calls = std::cell::Cell::new(0usize);
            let observe = || {
                calls.set(calls.get() + 1);
                false
            };
            mesh.prove_literal_embedding(used, Some(&observe)).unwrap();
            let at = calls.get();
            let seen = std::cell::Cell::new(0usize);
            let late_cancel = || {
                seen.set(seen.get() + 1);
                seen.get() == at
            };
            assert!(
                mesh.prove_literal_embedding(used, Some(&late_cancel))
                    .is_err()
            );
            assert_eq!(
                mesh.prove_literal_embedding(used, None)
                    .unwrap()
                    .work_used(),
                used
            );
            assert!(
                mesh.prove_literal_embedding(receipt.remaining_work, Some(&|| true))
                    .is_err()
            );
        }
    }
}
