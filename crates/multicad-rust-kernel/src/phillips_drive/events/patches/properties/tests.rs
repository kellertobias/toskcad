use super::*;
use std::cell::Cell;
fn cavity(metric: usize) -> Cavity {
    let (spec, domain) = crate::phillips_drive::test_fixture(metric);
    Cavity::new(Cutter::new(spec, domain).unwrap(), 1e-5, 100_000).unwrap()
}
fn overlap(a: Bounds, b: Bounds) {
    assert!(a[0] <= b[1] && b[0] <= a[1], "{a:?} / {b:?}");
}
fn negate(a: Bounds) -> Bounds {
    [-a[1], -a[0]]
}
#[test]
fn every_actual_metric_face_has_bounded_area_moments_centroid_and_exact_quarter_ownership() {
    for metric in 0..6 {
        let source = cavity(metric);
        let options = FacePropertyOptions {
            area_tolerance_mm2: 1e-8,
            first_moment_tolerance_mm3: 1e-8,
            centroid_tolerance_mm: 1e-8,
            ..FacePropertyOptions::default()
        };
        let properties = source
            .face_properties_with_cancel(options, &|| false)
            .unwrap();
        let topology = source.source().boundary_topology().unwrap();
        assert_eq!(properties.len(), 29);
        let bbox = source.source().bounding_box_mm();
        for (property, face) in properties.iter().zip(topology.faces()) {
            assert_eq!(property.face, face.source_identity());
            assert!(property.area_mm2[0] > 0.);
            assert!(property.area_mm2[1] - property.area_mm2[0] <= options.area_tolerance_mm2);
            for (axis, &centroid) in property.centroid_mm.iter().enumerate() {
                let moment = property.first_moment_mm3[axis];
                assert!(moment[1] - moment[0] <= options.first_moment_tolerance_mm3);
                assert!(centroid[1] - centroid[0] <= options.centroid_tolerance_mm);
                assert!(centroid[0] >= bbox[0][axis] - options.centroid_tolerance_mm);
                assert!(centroid[1] <= bbox[1][axis] + options.centroid_tolerance_mm);
            }
            let rotated_key = match property.face {
                FaceKey::Mouth => continue,
                FaceKey::Envelope { role, quarter } => FaceKey::Envelope {
                    role,
                    quarter: (quarter + 1) % 4,
                },
                FaceKey::Wall {
                    kind,
                    quarter,
                    lower,
                } => FaceKey::Wall {
                    kind,
                    quarter: (quarter + 1) % 4,
                    lower,
                },
            };
            let rotated = properties.iter().find(|p| p.face == rotated_key).unwrap();
            overlap(property.area_mm2, rotated.area_mm2);
            overlap(negate(property.centroid_mm[1]), rotated.centroid_mm[0]);
            overlap(property.centroid_mm[0], rotated.centroid_mm[1]);
            overlap(property.centroid_mm[2], rotated.centroid_mm[2]);
        }
        let mouth = properties
            .iter()
            .find(|p| p.face == FaceKey::Mouth)
            .unwrap();
        for c in mouth.centroid_mm {
            assert!(c[0] <= 0. && c[1] >= 0.);
        }
        let local = source
            .local_face_properties_with_cancel(properties[1].face, options, &|| false)
            .unwrap();
        overlap(local.area_mm2, properties[1].area_mm2);
        for i in 0..3 {
            overlap(local.centroid_mm[i], properties[1].centroid_mm[i]);
        }
    }
}
// Independent two-dimensional quadrature uses the original chart's actual
// tangent cross product. It never uses the analytic v-integrated density.
fn oracle(atlas: &PatchAtlas<'_>, subdivisions: usize) -> Vec<[f64; 4]> {
    const X: [f64; 4] = [
        -0.8611363115940526,
        -0.3399810435848563,
        0.3399810435848563,
        0.8611363115940526,
    ];
    const W: [f64; 4] = [
        0.34785484513745385,
        0.6521451548625461,
        0.6521451548625461,
        0.34785484513745385,
    ];
    let mut result = vec![[0.; 4]; atlas.topology.faces().len()];
    let mut compensation = result.clone();
    let midpoint = |b: Bounds| b[0].midpoint(b[1]);
    for (index, patch) in atlas.patches().iter().enumerate() {
        for row in 0..subdivisions {
            for col in 0..subdivisions {
                for i in 0..4 {
                    for j in 0..4 {
                        let u = (row as f64 + (1. + X[i]) / 2.) / subdivisions as f64;
                        let v = (col as f64 + (1. + X[j]) / 2.) / subdivisions as f64;
                        let state = atlas.patch_bounds(index, [[u, u], [v, v]]).unwrap();
                        let a = state.du_mm.map(midpoint);
                        let b = state.dv_mm.map(midpoint);
                        let cross = [
                            a[1] * b[2] - a[2] * b[1],
                            a[2] * b[0] - a[0] * b[2],
                            a[0] * b[1] - a[1] * b[0],
                        ];
                        let area =
                            cross.into_iter().map(|v| v * v).sum::<f64>().sqrt() * W[i] * W[j]
                                / (4. * (subdivisions * subdivisions) as f64);
                        let point = state.point_mm.map(midpoint);
                        let raw = [area, area * point[0], area * point[1], area * point[2]];
                        for (k, value) in raw.into_iter().enumerate() {
                            let y = value - compensation[patch.face_index()][k];
                            let next = result[patch.face_index()][k] + y;
                            compensation[patch.face_index()][k] =
                                (next - result[patch.face_index()][k]) - y;
                            result[patch.face_index()][k] = next;
                        }
                    }
                }
            }
        }
    }
    result
}
#[test]
fn actual_tangent_cross_product_integrals_agree_with_reduced_source_integrals() {
    for metric in [0, 1, 5] {
        let source = cavity(metric);
        let properties = source
            .face_properties_with_cancel(FacePropertyOptions::default(), &|| false)
            .unwrap();
        let atlas = source.source().surface_patches().unwrap();
        let coarse = oracle(&atlas, 2);
        let fine = oracle(&atlas, 4);
        for ((property, coarse), fine) in properties.iter().zip(coarse).zip(fine) {
            let intervals = [
                property.area_mm2,
                property.first_moment_mm3[0],
                property.first_moment_mm3[1],
                property.first_moment_mm3[2],
            ];
            for i in 0..4 {
                let scale = fine[i].abs().max(1.);
                let delta = (fine[i] - coarse[i]).abs();
                assert!(
                    delta <= 2e-5 * scale,
                    "metric{metric} {:?} quantity{i} oracle unconverged {delta}",
                    property.face
                );
                let allowance = 5. * delta + 2e-10 * scale;
                assert!(
                    fine[i] >= intervals[i][0] - allowance
                        && fine[i] <= intervals[i][1] + allowance,
                    "metric{metric} {:?} quantity{i}: {:?} oracle {} allowance{}",
                    property.face,
                    intervals[i],
                    fine[i],
                    allowance
                );
            }
        }
    }
}
#[test]
fn cancellation_resource_precision_and_invalid_queries_publish_no_partial_properties() {
    let source = cavity(1);
    let options = FacePropertyOptions::default();
    let count = Cell::new(0usize);
    let reference = source
        .face_properties_with_cancel(options, &|| {
            count.set(count.get() + 1);
            false
        })
        .unwrap();
    let calls = count.get();
    for target in [0, 1, 4, 16, calls / 4, calls / 2, calls - 2, calls - 1] {
        let count = Cell::new(0usize);
        let result = source.face_properties_with_cancel(options, &|| {
            let n = count.get();
            count.set(n + 1);
            n == target
        });
        assert!(
            matches!(result, Err(Error::Cancelled)),
            "target{target} calls{calls}: {result:?}"
        );
    }
    assert!(matches!(
        source.face_properties_with_cancel(
            FacePropertyOptions {
                maximum_cells: 1,
                ..options
            },
            &|| false
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        source.face_properties_with_cancel(
            FacePropertyOptions {
                area_tolerance_mm2: 1e-30,
                first_moment_tolerance_mm3: 1e-30,
                centroid_tolerance_mm: 1e-30,
                maximum_cells: 128
            },
            &|| false
        ),
        Err(Error::Limit) | Err(Error::Unresolved)
    ));
    for invalid in [
        FacePropertyOptions {
            area_tolerance_mm2: 0.,
            ..options
        },
        FacePropertyOptions {
            first_moment_tolerance_mm3: f64::NAN,
            ..options
        },
        FacePropertyOptions {
            centroid_tolerance_mm: -1.,
            ..options
        },
        FacePropertyOptions {
            maximum_cells: 0,
            ..options
        },
        FacePropertyOptions {
            maximum_cells: MAXIMUM_CELLS + 1,
            ..options
        },
    ] {
        assert!(matches!(
            source.face_properties_with_cancel(invalid, &|| false),
            Err(Error::Invalid)
        ));
        assert!(matches!(
            source.local_face_properties_with_cancel(FaceKey::Mouth, invalid, &|| false),
            Err(Error::Invalid)
        ));
    }
    assert!(matches!(
        source.local_face_properties_with_cancel(
            FaceKey::Envelope {
                role: EnvelopeRole::Floor,
                quarter: 255
            },
            options,
            &|| false
        ),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        source.local_face_properties_with_cancel(FaceKey::Mouth, options, &|| true),
        Err(Error::Cancelled)
    ));
    let after = source
        .face_properties_with_cancel(options, &|| false)
        .unwrap();
    for (a, b) in reference.iter().zip(after) {
        assert_eq!(a.face, b.face);
        assert_eq!(a.area_mm2, b.area_mm2);
        assert_eq!(a.centroid_mm, b.centroid_mm);
    }
}
