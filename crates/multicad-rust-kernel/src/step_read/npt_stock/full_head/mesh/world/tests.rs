use super::*;
use crate::step_read::{
    npt_stock::full_head::transport::encoding,
    owned_insert::{CompoundDocument, DocumentOptions},
};
use num_bigint::BigInt;
use num_rational::BigRational as Q;
use num_traits::Signed;
use std::sync::Arc;
fn q(x: f64) -> Q {
    Q::from_float(x).unwrap()
}
fn decimal(s: &str) -> Q {
    let (m, e) = s.split_once('E').unwrap_or((s, "0"));
    let f = m.split_once('.').map_or(0, |(_, v)| v.len()) as i32;
    let n = m.replace('.', "").parse::<BigInt>().unwrap();
    let p = e.parse::<i32>().unwrap() - f;
    if p >= 0 {
        Q::from_integer(n * BigInt::from(10).pow(p as u32))
    } else {
        Q::new(n, BigInt::from(10).pow((-p) as u32))
    }
}
fn source(quarter: bool) -> MappedActualHead {
    let local = Arc::new(
        OwnedActualHead::read(
            Arc::from(super::super::super::source_tests::candidate(quarter, 18.)),
            Default::default(),
            None,
        )
        .unwrap(),
    );
    let f = LiteralOccurrencePlacement::new(
        [[1., 0., 0.], [0., -1., 0.], [0., 0., -1.]],
        [0., 0., local.head_stations_mm()[0]],
        Default::default(),
        None,
    )
    .unwrap();
    let candidate = encoding::propose_mapped(&local, &f, Default::default(), None).unwrap();
    MappedActualHead::read_reusing_local(candidate.actual_bytes(), local, Default::default(), None)
        .unwrap()
}
fn placement() -> LiteralOccurrencePlacement {
    LiteralOccurrencePlacement::new(
        [[0.6, -0.8, 0.], [0.8, 0.6, 0.], [0., 0., 1.]],
        [20., -30., 11.],
        Default::default(),
        None,
    )
    .unwrap()
}
fn actual_file_origin(input: &str) -> ([Q; 3], Q) {
    use crate::step_read::part21::Value;
    let doc = CompoundDocument::parse(input, DocumentOptions::default(), None).unwrap();
    let link = doc.represented_members()[0].mapping_link().unwrap();
    let parts = &doc.records()[&link.operator_entity()];
    let base = parts
        .iter()
        .find(|e| e.name == "CARTESIAN_TRANSFORMATION_OPERATOR")
        .unwrap();
    let Value::Reference(origin) = base.args[2] else {
        panic!()
    };
    let raw = crate::step_read::part21::numeric_spellings(&doc.records()[&origin][0]).unwrap();
    let origin = std::array::from_fn(|i| decimal(raw[i]));
    // Actual encoded direction leaves independently establish this proper F,
    // rather than inferring it from a feature name or a regenerated recipe.
    let Value::Reference(x) = base.args[0] else {
        panic!()
    };
    let Value::Reference(y) = base.args[1] else {
        panic!()
    };
    let dimensional = parts
        .iter()
        .find(|e| e.name == "CARTESIAN_TRANSFORMATION_OPERATOR_3D")
        .unwrap();
    let Value::Reference(z) = dimensional.args[0] else {
        panic!()
    };
    for (id, expected) in [(x, [1., 0., 0.]), (y, [0., -1., 0.]), (z, [0., 0., -1.])] {
        let raw = crate::step_read::part21::numeric_spellings(&doc.records()[&id][0]).unwrap();
        let got: [Q; 3] = std::array::from_fn(|i| decimal(raw[i]));
        for i in 0..3 {
            if expected[i] == 0. {
                assert_eq!(got[i], q(0.));
            } else {
                assert!(got[i].clone() * q(expected[i]) > q(0.));
            }
        }
    }
    let scale = crate::step_read::part21::numeric_spellings(base).unwrap();
    assert_eq!(scale.len(), 1);
    (origin, decimal(scale[0]))
}
#[test]
fn actual_file_intrinsic_pose_then_only_extra_native_occurrence_has_independent_world_corner_and_flux_oracles()
 {
    for quarter in [false, true] {
        let source = source(quarter);
        let (original_origin, original_scale) = actual_file_origin(source.actual_file());
        let extra = placement();
        let near_rigid = LiteralOccurrencePlacement::new(
            [[0.60000001, -0.8, 0.], [0.8, 0.6, 0.], [0., 0., 1.]],
            [20., -30., 11.],
            Default::default(),
            None,
        )
        .unwrap();
        for occurrence in [None, Some(&extra), Some(&near_rigid)] {
            let placed = occurrence.is_some();
            let options = HeadMeshOptions {
                tessellation: Tessellation {
                    deflection: 0.05,
                    angle: 0.3,
                    triangle_limit: 500_000,
                },
                ..Default::default()
            };
            let candidate = if let Some(occurrence) = occurrence {
                source.placed_mesh_candidate(occurrence, options, None)
            } else {
                source.mesh_candidate(options, None)
            }
            .unwrap();
            assert!(std::ptr::eq(candidate.source(), &source));
            assert_eq!(candidate.occurrence().is_some(), placed);
            let r = candidate.receipt();
            assert!(r.source_position_error_mm <= 0.05 && r.normal_tangent_cone <= 0.3);
            assert_eq!(r.aggregate_work + r.remaining_work, MAX_MESH_WORK);
            assert_eq!(
                r.aggregate_work,
                r.inherited_mapped_material_work
                    + r.inherited_occurrence_work
                    + r.local.aggregate_work
                    - source.local_geometry().receipt().aggregate_work
                    + r.world_transfer_work
            );
            let independent = manufacturing::certify(candidate.mesh()).unwrap();
            assert_eq!(independent.edge_pairs, candidate.connectivity().edge_pairs);
            assert_eq!(
                independent.vertex_links,
                candidate.connectivity().vertex_links
            );
            let flux = candidate.connectivity().signed_volume_mm3;
            assert!((independent.signed_volume_mm3 - flux).abs() < 1e-10 * flux);
            for &p in &source.local_geometry().head().vertices {
                let mut expected = [
                    &original_scale * q(p[0]) + &original_origin[0],
                    -&original_scale * q(p[1]) + &original_origin[1],
                    -&original_scale * q(p[2]) + &original_origin[2],
                ];
                if let Some(occurrence) = occurrence {
                    for factor in occurrence.transforms() {
                        expected = std::array::from_fn(|i| {
                            (0..3).fold(q(factor.translation_mm()[i]), |s, j| {
                                s + q(factor.matrix()[i][j]) * &expected[j]
                            })
                        });
                    }
                }
                let closest = candidate
                    .mesh()
                    .vertices
                    .as_chunks::<3>()
                    .0
                    .iter()
                    .min_by(|a, b| {
                        let distance = |v: &[f64; 3]| {
                            (0..3).map(|i| (q(v[i]) - &expected[i]).abs()).sum::<Q>()
                        };
                        distance(a).cmp(&distance(b))
                    })
                    .unwrap();
                let distance: Q = (0..3).map(|i| (q(closest[i]) - &expected[i]).abs()).sum();
                assert!(distance <= q(r.literal_world_rounding_mm));
            }
            eprintln!(
                "ACTUAL_WORLD_HEAD quarter={quarter} placed={placed} facets={} {r:?}",
                candidate.mesh().face_indices.len()
            );
        }
    }
}
#[test]
fn actual_world_mesh_work_storage_cancel_and_retry_are_atomic_and_bounded() {
    use std::cell::Cell;
    let source = source(false);
    let extra = placement();
    let good = source
        .placed_mesh_candidate(&extra, Default::default(), None)
        .unwrap();
    let r = good.receipt();
    let expected = good.mesh().clone();
    drop(good);
    for (maximum_work, maximum_bytes, ok) in [
        (r.aggregate_work - 1, MAX_BYTES, false),
        (r.aggregate_work, MAX_BYTES, true),
        (MAX_MESH_WORK, r.accounted_peak_bytes - 1, false),
        (MAX_MESH_WORK, r.accounted_peak_bytes, true),
    ] {
        let result = source.placed_mesh_candidate(
            &extra,
            HeadMeshOptions {
                maximum_work,
                maximum_bytes,
                ..Default::default()
            },
            None,
        );
        assert_eq!(result.is_ok(), ok);
        if !ok {
            assert!(matches!(result, Err(Error::Limit)));
        }
    }
    let n = Cell::new(0);
    let observe = || {
        n.set(n.get() + 1);
        false
    };
    source
        .placed_mesh_candidate(&extra, Default::default(), Some(&observe))
        .unwrap();
    for at in [1, n.get() / 2, n.get()] {
        let seen = Cell::new(0);
        let cancel = || {
            seen.set(seen.get() + 1);
            seen.get() >= at
        };
        assert!(matches!(
            source.placed_mesh_candidate(&extra, Default::default(), Some(&cancel)),
            Err(Error::Cancelled)
        ));
    }
    let more = source
        .placed_mesh_candidate(
            &extra,
            HeadMeshOptions {
                additional_resident_bytes: 4096,
                ..Default::default()
            },
            None,
        )
        .unwrap();
    assert_eq!(more.receipt().aggregate_work, r.aggregate_work);
    assert_eq!(
        more.receipt().accounted_peak_bytes,
        r.accounted_peak_bytes + 4096
    );
    let retry = source
        .placed_mesh_candidate(&extra, Default::default(), None)
        .unwrap();
    assert_eq!(retry.mesh(), &expected);
    assert_eq!(retry.receipt().aggregate_work, r.aggregate_work);
}

#[test]
fn actual_cold_uniform_scale_changes_whole_geometry_without_reapplying_intrinsic_frame() {
    use crate::step_read::part21::Value;
    let source = source(false);
    let doc =
        CompoundDocument::parse(source.actual_file(), DocumentOptions::default(), None).unwrap();
    let operator = doc.represented_members()[0]
        .mapping_link()
        .unwrap()
        .operator_entity();
    let base = doc.records()[&operator]
        .iter()
        .find(|e| e.name == "CARTESIAN_TRANSFORMATION_OPERATOR")
        .unwrap();
    assert!(matches!(base.args[3], Value::Number(_)));
    let raw = crate::step_read::part21::numeric_spellings(base).unwrap();
    let old = raw[0];
    let record = base.source.unwrap();
    let replacement = record.replacen(old, "2.", 1);
    assert_ne!(record, replacement);
    let input = source.actual_file().replacen(record, &replacement, 1);
    drop(doc);
    let scaled = MappedActualHead::read(Arc::from(input), Default::default(), None).unwrap();
    let (_, actual_scale) = actual_file_origin(scaled.actual_file());
    assert_eq!(actual_scale, q(2.));
    let options = HeadMeshOptions {
        tessellation: Tessellation {
            deflection: 0.1,
            angle: 0.3,
            triangle_limit: 500_000,
        },
        ..Default::default()
    };
    let mesh = scaled.mesh_candidate(options, None).unwrap();
    assert!(mesh.receipt().source_position_error_mm <= 0.1);
    let independent = manufacturing::certify(mesh.mesh()).unwrap();
    assert!(
        (independent.signed_volume_mm3 - mesh.connectivity().signed_volume_mm3).abs()
            < 1e-10 * independent.signed_volume_mm3
    );
    let local_volume = source.local_geometry().local_sewn_volume_mm3();
    let world_volume = scaled.world_sewn_volume_mm3();
    assert!(world_volume[0] <= 8. * local_volume[1]);
    assert!(world_volume[1] >= 8. * local_volume[0]);
}
