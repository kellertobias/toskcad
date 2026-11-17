//! Public actual-file mesh plus independent Native literal embedding gate.
//! Explicit actual provider workflow; global import/STEP flags remain closed.
use multicad_core::{document::EntityId, references::TopologyKind};
use multicad_kernel::{
    GeometryPolicy, MeshOptions, PlaneFrame, RigidTransform, Solid, rust_kernel as pure,
};
use pure::axial_pipe_thread::authored::{AuthoredSection, AuthoredSpec, AuthoredThread};
use pure::step::npt_stock::{CurveOptions, CurvePlan, SurfacePlan, write_source_candidate};
use pure::step_read::{
    npt_stock::{
        LiteralOccurrencePlacement, Parsed, Support,
        full_head::{
            OwnedActualHead,
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
fn actual(quarter: bool) -> Arc<MappedActualHead> {
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
    let output = encoding::propose_mapped(&local, &intrinsic, Default::default(), None).unwrap();
    // Independent public cold file decoding, not a retained-source recipe.
    Arc::new(MappedActualHead::read(output.actual_bytes(), Default::default(), None).unwrap())
}
fn pose() -> RigidTransform {
    RigidTransform {
        rotation: [[0.6, 0., 0.8], [0., 1., 0.], [-0.8, 0., 0.6]],
        translation_mm: [20., -30., 11.],
    }
}
fn place(p: [f64; 3], placed: bool) -> [f64; 3] {
    if placed {
        [
            0.6 * p[0] + 0.8 * p[2] + 20.,
            p[1] - 30.,
            -0.8 * p[0] + 0.6 * p[2] + 11.,
        ]
    } else {
        p
    }
}
#[test]
fn actual_native_head_cold_source_planes_selected_properties_queries_and_manufacturing_both_sizes_and_poses()
 {
    for quarter in [false, true] {
        let source = actual(quarter);
        let length = source.local_geometry().head_stations_mm()[0];
        for placed in [false, true] {
            let feature = EntityId::from_uuid_bits(0xabc);
            let initial = Solid::from_actual_clipped_head_source(
                Arc::clone(&source),
                feature,
                GeometryPolicy::default(),
                None,
            )
            .unwrap();
            let solid = if placed {
                initial.transformed(pose()).unwrap()
            } else {
                initial
            };
            assert!(solid.uses_rust_geometry());
            assert!(solid.supports_step_export());
            assert!(std::ptr::eq(
                solid.rust_actual_clipped_head_source().unwrap(),
                source.as_ref()
            ));
            assert_eq!(
                solid.rust_actual_clipped_head_occurrence().is_some(),
                placed
            );
            assert_eq!(
                solid
                    .rust_actual_clipped_head_runtime_evidence()
                    .unwrap()
                    .measured_faces,
                0
            );
            for (index, local_origin) in [(8, [0., 0., length]), (9, [0., 0., 0.])] {
                let (origin, normal) = solid.face_plane(index).unwrap().unwrap();
                let expected_origin = place(local_origin, placed);
                let expected_normal = if placed { [0.8, 0., 0.6] } else { [0., 0., 1.] };
                assert!(
                    origin
                        .into_iter()
                        .zip(expected_origin)
                        .map(|(a, b)| (a - b).abs())
                        .sum::<f64>()
                        < 1e-7
                );
                assert!(
                    normal
                        .into_iter()
                        .zip(expected_normal)
                        .map(|(a, b)| (a - b).abs())
                        .sum::<f64>()
                        < 1e-10
                );
            }
            let mut planes = 0;
            for i in 0..17 {
                if solid.face_plane(i).unwrap().is_some() {
                    planes += 1;
                }
            }
            assert_eq!(planes, 9);
            assert_eq!(
                solid
                    .rust_actual_clipped_head_runtime_evidence()
                    .unwrap()
                    .face_measure_builds,
                0
            );
            let topology = solid.topology_snapshot(feature, 7).unwrap();
            assert_eq!(
                topology
                    .elements
                    .iter()
                    .filter(|e| e.key.kind == TopologyKind::Face)
                    .count(),
                17
            );
            assert_eq!(
                topology
                    .elements
                    .iter()
                    .filter(|e| e.key.kind == TopologyKind::Edge)
                    .count(),
                38
            );
            assert_eq!(
                topology
                    .elements
                    .iter()
                    .filter(|e| e.key.kind == TopologyKind::Vertex)
                    .count(),
                24
            );
            assert_eq!(
                solid
                    .rust_actual_clipped_head_runtime_evidence()
                    .unwrap()
                    .face_measure_builds,
                0
            );

            let top = 9 + source
                .local_geometry()
                .head()
                .faces
                .iter()
                .position(|f| f.outer.len() == 6 && f.inner.is_empty())
                .unwrap() as u32;
            let frame = topology
                .elements
                .iter()
                .find(|e| e.key.kind == TopologyKind::Face && e.key.ordinal == top)
                .unwrap()
                .frames[0]
                .frame;
            let frame = PlaneFrame {
                origin: frame.origin_mm,
                u: frame.u,
                v: frame.v,
            };
            let contour = solid.planar_face_contours_with_sources(top, frame).unwrap();
            assert_eq!(contour.contours.loops.len(), 1);
            assert_eq!(contour.contours.loops[0].edges.len(), 6);
            assert_eq!(
                solid
                    .rust_actual_clipped_head_runtime_evidence()
                    .unwrap()
                    .face_measure_builds,
                0
            );
            assert!(solid.planar_face_contours(0, frame).is_err());
            assert!(solid.planar_face_contours(17, frame).is_err());
            assert_eq!(
                solid
                    .rust_actual_clipped_head_runtime_evidence()
                    .unwrap()
                    .face_measure_builds,
                0
            );
            let overlay = solid.topology_with_display(feature, 7, Some(0.05)).unwrap();
            assert_eq!(overlay.edges.len(), 38);
            assert_eq!(overlay.vertices.len(), 24);
            assert_eq!(solid.display_edges(0.05).unwrap().len(), 38);
            assert_eq!(solid.display_vertices().unwrap().len(), 24);
            assert_eq!(
                solid
                    .rust_actual_clipped_head_runtime_evidence()
                    .unwrap()
                    .face_measure_builds,
                0
            );

            for (p, expected) in [
                ([0., 0., -1.], true),
                ([0., 0., length / 2.], true),
                ([10., -8., -1.], false),
                ([30., 0., length / 2.], false),
            ] {
                assert_eq!(
                    solid
                        .rust_actual_clipped_head_contains(place(p, placed))
                        .unwrap(),
                    expected
                );
            }
            let bounds = solid.bounds().unwrap();
            assert!(bounds.min.into_iter().chain(bounds.max).all(f64::is_finite));
            assert!(solid.volume_mm3().unwrap() > 0.);
            let planar = solid.face(9).unwrap();
            assert!(planar.area_mm2 > 0. && planar.plane.is_some());
            let curved = solid.face(0).unwrap();
            assert!(curved.area_mm2 > 0. && curved.plane.is_none());
            let before = solid.rust_actual_clipped_head_runtime_evidence().unwrap();
            assert_eq!(before.measured_faces, 2);
            assert_eq!(before.face_measure_builds, 2);
            assert_eq!(solid.face(0).unwrap().area_mm2, curved.area_mm2);
            assert_eq!(
                solid
                    .rust_actual_clipped_head_runtime_evidence()
                    .unwrap()
                    .face_measure_builds,
                2
            );
            let strict = MeshOptions {
                linear_deflection_mm: 0.05,
                angular_deflection_rad: 0.3,
                maximum_triangles: 500_000,
            };
            let m = solid.manufacturing_mesh_with_receipt(strict).unwrap();
            assert_eq!(m.edge_pairs as usize * 2, m.mesh.indices.len());
            assert_eq!(m.vertex_links as usize, m.mesh.vertices.len() / 3);
            let r = solid.rust_actual_clipped_head_runtime_evidence().unwrap();
            assert!(
                r.last_mesh_work.unwrap() <= 1_600_000_000
                    && r.last_mesh_peak_bytes.unwrap() <= 64 * 1024 * 1024
            );
            assert!(
                r.last_mesh_source_error_mm.unwrap() <= 0.05
                    && r.last_mesh_normal_tangent_cone.unwrap() <= 0.3
            );
            eprintln!(
                "ACTUAL_NATIVE_PROVIDER quarter={quarter} placed={placed} facets={} work={} peak={} pairs={}",
                m.mesh.face_indices.len(),
                r.last_mesh_work.unwrap(),
                r.last_mesh_peak_bytes.unwrap(),
                r.last_embedding_pairs.unwrap()
            );
            let display = solid.mesh(MeshOptions::default()).unwrap();
            let builds = solid
                .rust_actual_clipped_head_runtime_evidence()
                .unwrap()
                .mesh_builds;
            assert_eq!(solid.mesh(MeshOptions::default()).unwrap(), display);
            assert_eq!(
                solid
                    .rust_actual_clipped_head_runtime_evidence()
                    .unwrap()
                    .mesh_builds,
                builds
            );
            let previous = solid
                .rust_actual_clipped_head_runtime_evidence()
                .unwrap()
                .last_mesh_work;
            assert!(
                solid
                    .manufacturing_mesh(MeshOptions {
                        maximum_triangles: 1,
                        ..strict
                    })
                    .is_err()
            );
            assert_eq!(
                solid
                    .rust_actual_clipped_head_runtime_evidence()
                    .unwrap()
                    .last_mesh_work,
                previous
            );
            // Existing display cache stays resident through fresh manufacturing.
            solid.manufacturing_mesh(strict).unwrap();
            assert_eq!(solid.mesh(MeshOptions::default()).unwrap(), display);
            let temp = tempfile::tempdir().unwrap();
            let path = temp.path().join("must-stay.step");
            std::fs::write(&path, b"prior bytes").unwrap();
            let mut low =
                pure::step_read::npt_stock::full_head::transport::replay::ReplayOptions::default();
            low.resources.maximum_work = 1;
            assert!(
                solid
                    .write_actual_clipped_head_step_with_verification(&path, low, None)
                    .unwrap()
                    .is_err()
            );
            assert_eq!(std::fs::read(path).unwrap(), b"prior bytes");
            assert!(solid.face(17).is_err());
            assert!(solid.face_plane(17).is_err());
        }
    }
}
#[test]
fn actual_provider_cancellation_and_invalid_pose_preserve_original_source() {
    let source = actual(false);
    let feature = EntityId::from_uuid_bits(0xabc);
    assert!(
        Solid::from_actual_clipped_head_source(
            Arc::clone(&source),
            feature,
            GeometryPolicy::default(),
            Some(&|| true)
        )
        .is_err()
    );

    let count = std::cell::Cell::new(0usize);
    let observe = || {
        count.set(count.get() + 1);
        false
    };
    Solid::from_actual_clipped_head_source(
        Arc::clone(&source),
        feature,
        GeometryPolicy::default(),
        Some(&observe),
    )
    .unwrap();
    let at = count.get();
    let seen = std::cell::Cell::new(0usize);
    let cancel = || {
        seen.set(seen.get() + 1);
        seen.get() == at
    };
    assert!(
        Solid::from_actual_clipped_head_source(
            Arc::clone(&source),
            feature,
            GeometryPolicy::default(),
            Some(&cancel)
        )
        .is_err()
    );
    let solid = Solid::from_actual_clipped_head_source(
        Arc::clone(&source),
        feature,
        GeometryPolicy::default(),
        None,
    )
    .unwrap();
    let mut invalid = pose();
    invalid.rotation[0][0] += 0.001;
    assert!(solid.transformed(invalid).is_err());
    assert!(std::ptr::eq(
        solid.rust_actual_clipped_head_source().unwrap(),
        source.as_ref()
    ));
    assert_eq!(
        solid
            .rust_actual_clipped_head_runtime_evidence()
            .unwrap()
            .measured_faces,
        0
    );
}
