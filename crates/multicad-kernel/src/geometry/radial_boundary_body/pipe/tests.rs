use super::*;
fn spec() -> AuthoredSpec {
    let p = 25.4 / 28.;
    AuthoredSpec {
        major_radius_mm: 9.728 / 2. - 4. / 32.,
        pitch_mm: p,
        length_mm: 6.5,
        radial_taper: 1. / 32.,
        section: AuthoredSection::BsptRounded {
            depth_mm: p / (3. * 27.5_f64.to_radians().portable_tan()),
        },
        crest_radius_limit_mm: None,
    }
}
#[test]
fn authored_pipe_has_actual_owned_geometry_and_distinct_source_provenance() {
    let feature = EntityId::new();
    let law = spec();
    let body = Solid::rust_authored_axial_pipe_thread_with_interrupt(
        law,
        feature,
        PlaneFrame::default(),
        GeometryPolicy::default(),
        None,
    )
    .unwrap();
    let (source, pose) = body.rust_authored_axial_pipe_thread_source().unwrap();
    assert_eq!(source.transport.spec(), law);
    assert_eq!(pose, PlaneFrame::default());
    assert!(body.rust_retained_helical_thread_source().is_none());
    let top = body.topology_snapshot(EntityId::new(), 0).unwrap();
    assert!(
        top.elements
            .iter()
            .any(|e| e.origins.iter().any(|o| o.role == "cap.end"))
    );
    assert_eq!(source.band_faces.len(), 5);
    assert_eq!(source.root_faces.len(), 1);
    let copy = body.share().unwrap();
    assert!(std::ptr::eq(
        source,
        copy.rust_authored_axial_pipe_thread_source().unwrap().0
    ));
    let posed = body
        .transformed(RigidTransform {
            rotation: [[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]],
            translation_mm: [10., 20., 30.],
        })
        .unwrap();
    let (moved, occurrence) = posed.rust_authored_axial_pipe_thread_source().unwrap();
    assert!(std::ptr::eq(source, moved));
    assert_eq!(occurrence.origin, [10., 20., 30.]);
    assert_eq!(source.retained_radius_upper_bound_mm, None);
    let v = body.volume_mm3().unwrap();
    assert!((v - 423.7556150148513).abs() < 1e-3);
}
#[test]
fn invalid_cancelled_and_active_stock_clips_do_not_become_pipe_bodies() {
    let p = GeometryPolicy::default();
    let frame = PlaneFrame::default();
    let feature = EntityId::new();
    assert_eq!(
        Solid::rust_authored_axial_pipe_thread_with_interrupt(
            spec(),
            EntityId::nil(),
            frame,
            p,
            None
        )
        .err()
        .unwrap()
        .kind,
        GeometryErrorKind::Input
    );
    assert_eq!(
        Solid::rust_authored_axial_pipe_thread_with_interrupt(
            spec(),
            feature,
            frame,
            p,
            Some(&|| true)
        )
        .err()
        .unwrap()
        .kind,
        GeometryErrorKind::Cancelled
    );
    let calls = std::cell::Cell::new(0);
    let cancel = || {
        let n = calls.get() + 1;
        calls.set(n);
        n > 1
    };
    assert_eq!(
        Solid::rust_authored_axial_pipe_thread_with_interrupt(
            spec(),
            feature,
            frame,
            p,
            Some(&cancel)
        )
        .err()
        .unwrap()
        .kind,
        GeometryErrorKind::Cancelled
    );
    let mut clipped = spec();
    clipped.crest_radius_limit_mm = Some(4.5);
    assert_eq!(
        Solid::rust_authored_axial_pipe_thread_with_interrupt(clipped, feature, frame, p, None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Unsupported
    );
}

#[test]
fn pipe_manufacturing_and_step_readback_keep_actual_geometry_without_guessing_authorship() {
    let body = Solid::rust_authored_axial_pipe_thread_with_interrupt(
        spec(),
        EntityId::new(),
        PlaneFrame::default(),
        GeometryPolicy::default(),
        None,
    )
    .unwrap();
    let mesh = body
        .manufacturing_mesh_with_receipt(MeshOptions {
            linear_deflection_mm: 0.1,
            angular_deflection_rad: 0.3,
            maximum_triangles: 1_000_000,
        })
        .unwrap();
    assert!(mesh.edge_pairs > 0 && mesh.vertex_links > 0);
    assert_eq!(mesh.removed_exact_zero_triangles, 0);
    let file = tempfile::NamedTempFile::new().unwrap();
    body.write_step(file.path()).unwrap();
    let text = std::fs::read_to_string(file.path()).unwrap();
    assert_eq!(text.matches("ADVANCED_FACE(").count(), 7);
    let proof = body
        .verify_retained_step_source(&text, 1e-7)
        .unwrap()
        .unwrap();
    assert_eq!(proof.retained_members.len(), 1);
    assert!(proof.retained_members[0].decoded_volume_bounds_mm3[0] > 0.);
    let decoded =
        Solid::import_rust_step(&text, EntityId::new(), GeometryPolicy::default()).unwrap();
    assert!(decoded.rust_radial_boundary().is_some());
    assert!(decoded.rust_authored_axial_pipe_thread_source().is_none());
    assert!(decoded.rust_retained_helical_thread_source().is_none());
    assert!((decoded.volume_mm3().unwrap() - body.volume_mm3().unwrap()).abs() < 1e-3);
}
#[test]
fn pipe_semantic_refs_survive_authored_dimension_and_occurrence_edits() {
    use multicad_core::references::{GeometryIntent, Resolution, resolve};
    let feature = EntityId::new();
    let owner = EntityId::new();
    let before = Solid::rust_authored_axial_pipe_thread_with_interrupt(
        spec(),
        feature,
        PlaneFrame::default(),
        GeometryPolicy::default(),
        None,
    )
    .unwrap();
    let topology = before.topology_snapshot(owner, 0).unwrap();
    let cap = topology
        .elements
        .iter()
        .find(|e| e.key.kind == TopologyKind::Face && e.origins.iter().any(|o| o.role == "cap.end"))
        .unwrap();
    let intent = GeometryIntent::capture(&topology, cap.key).unwrap();
    let mut edited = spec();
    edited.length_mm = 6.7;
    let frame = PlaneFrame {
        origin: [17., -8., 13.],
        u: [0.8, 0.6, 0.],
        v: [-0.36, 0.48, 0.8],
    };
    let after = Solid::rust_authored_axial_pipe_thread_with_interrupt(
        edited,
        feature,
        frame,
        GeometryPolicy::default(),
        None,
    )
    .unwrap();
    let current = after.topology_snapshot(owner, 1).unwrap();
    assert!(matches!(
        resolve(&intent, &current, 1),
        Resolution::Resolved { .. }
    ));
    assert_eq!(
        after
            .rust_authored_axial_pipe_thread_source()
            .unwrap()
            .0
            .construction_frame,
        frame
    );
}

#[test]
fn stock_admission_proves_the_whole_retained_body_and_preserves_the_limit() {
    let mut law = spec();
    law.crest_radius_limit_mm = Some(10.);
    let body = Solid::rust_authored_axial_pipe_thread_with_interrupt(
        law,
        EntityId::new(),
        PlaneFrame::default(),
        GeometryPolicy::default(),
        None,
    )
    .unwrap();
    let source = body.rust_authored_axial_pipe_thread_source().unwrap().0;
    assert_eq!(source.transport.spec(), law);
    let upper = source.retained_radius_upper_bound_mm.unwrap();
    assert!(upper > law.major_radius_mm + law.radial_taper * law.length_mm && upper < 10.);
    let mut near = law;
    near.crest_radius_limit_mm =
        Some(law.major_radius_mm + law.radial_taper * law.length_mm + 1e-10);
    let rejected = Solid::rust_authored_axial_pipe_thread_with_interrupt(
        near,
        EntityId::new(),
        PlaneFrame::default(),
        GeometryPolicy::default(),
        None,
    )
    .err()
    .unwrap();
    assert_eq!(rejected.kind, GeometryErrorKind::Unsupported);
    assert!(rejected.to_string().contains("stock OD"));
    let posed = Solid::rust_authored_axial_pipe_thread_with_interrupt(
        law,
        EntityId::new(),
        PlaneFrame {
            origin: [1., 0., 0.],
            ..PlaneFrame::default()
        },
        GeometryPolicy::default(),
        None,
    )
    .err()
    .unwrap();
    assert_eq!(posed.kind, GeometryErrorKind::Unsupported);
    assert!(posed.to_string().contains("posed pipe"));
}

#[test]
fn npt_source_is_not_a_bspt_or_custom_alias() {
    let p = 1.25;
    let k = 1. / 32.;
    let depth = 0.5;
    let flank = 3_f64.sqrt();
    let law = AuthoredSpec {
        major_radius_mm: 4.,
        pitch_mm: p,
        length_mm: 2.5,
        radial_taper: k,
        section: AuthoredSection::NptFlat {
            depth_mm: depth,
            flat_width_mm: (p - depth * (1. / (flank - k) + 1. / (flank + k))) / 2.,
        },
        crest_radius_limit_mm: Some(10.),
    };
    let body = Solid::rust_authored_axial_pipe_thread_with_interrupt(
        law,
        EntityId::new(),
        PlaneFrame::default(),
        GeometryPolicy::default(),
        None,
    )
    .unwrap();
    let source = body.rust_authored_axial_pipe_thread_source().unwrap().0;
    assert_eq!(source.transport.spec(), law);
    assert_eq!(source.band_faces.len(), 3);
    assert_eq!(source.root_faces.len(), 1);
    assert!(body.rust_retained_helical_thread_source().is_none());
    assert!(body.volume_mm3().unwrap() > 100.);
}

#[test]
fn semantic_provenance_refuses_inconsistent_actual_patch_roles() {
    let body = Solid::rust_authored_axial_pipe_thread_with_interrupt(
        spec(),
        EntityId::new(),
        PlaneFrame::default(),
        GeometryPolicy::default(),
        None,
    )
    .unwrap();
    let source = body.rust_authored_axial_pipe_thread_source().unwrap().0;
    let boundary = body.rust_radial_boundary().unwrap().0;
    let mut wrong = source.clone();
    wrong.cap_faces[0] = wrong.band_faces[0];
    assert!(AuthoredPipe::new(wrong, boundary.input()).is_err());
    let mut wrong = source.clone();
    wrong.band_faces.pop();
    assert!(AuthoredPipe::new(wrong, boundary.input()).is_err());
    let mut wrong = source.clone();
    wrong.band_faces[0] = wrong.cap_faces[0];
    assert!(AuthoredPipe::new(wrong, boundary.input()).is_err());
    let mut wrong = source.clone();
    wrong.band_faces[0] = wrong.band_faces[1];
    assert!(AuthoredPipe::new(wrong, boundary.input()).is_err());
    let mut wrong_kind = source.clone();
    let k = 1. / 32.;
    let flank = 3_f64.sqrt();
    wrong_kind.transport = AuthoredThread::new(
        AuthoredSpec {
            major_radius_mm: 4.,
            pitch_mm: 1.25,
            length_mm: 2.5,
            radial_taper: k,
            section: AuthoredSection::NptFlat {
                depth_mm: 0.5,
                flat_width_mm: (1.25 - 0.5 * (1. / (flank - k) + 1. / (flank + k))) / 2.,
            },
            crest_radius_limit_mm: None,
        },
        1e-7,
        1e-5,
    )
    .unwrap();
    wrong_kind.band_faces.truncate(3);
    wrong_kind.root_faces = vec![source.band_faces[4]];
    assert!(AuthoredPipe::new(wrong_kind, boundary.input()).is_err());
    let mut wrong = source.clone();
    wrong.root_faces.clear();
    assert!(AuthoredPipe::new(wrong, boundary.input()).is_err());
    let mut wrong = source.clone();
    wrong.root_faces[0] = wrong.band_faces[1];
    assert!(AuthoredPipe::new(wrong, boundary.input()).is_err());
}

#[test]
fn semantic_provenance_requires_actual_canonical_incidence_not_just_equal_volume() {
    use multicad_rust_kernel::radial_tensor::{BoundaryOwnership, Curve};
    let p = 1.25;
    let k = 1. / 32.;
    let d = 0.5;
    let f = 3_f64.sqrt();
    let law = AuthoredSpec {
        major_radius_mm: 4.,
        pitch_mm: p,
        length_mm: 2.5,
        radial_taper: k,
        section: AuthoredSection::NptFlat {
            depth_mm: d,
            flat_width_mm: (p - d * (1. / (f - k) + 1. / (f + k))) / 2.,
        },
        crest_radius_limit_mm: None,
    };
    let solid = Solid::rust_authored_axial_pipe_thread_with_interrupt(
        law,
        EntityId::new(),
        PlaneFrame::default(),
        GeometryPolicy::default(),
        None,
    )
    .unwrap();
    let source = solid.rust_authored_axial_pipe_thread_source().unwrap().0;
    let boundary = solid.rust_radial_boundary().unwrap().0;
    let mut input = boundary.input().clone();
    let index = input.faces[source.cap_faces[0]]
        .outer
        .iter()
        .map(|x| x.0)
        .find(|&i| matches!(input.curves[input.edges[i].curve], Curve::Circle { .. }))
        .unwrap();
    let edge = input.edges[index];
    let Curve::Circle { axis, radius_mm } = input.curves[edge.curve] else {
        panic!()
    };
    let a = input.vertices[edge.vertices[0]];
    let b = input.vertices[edge.vertices[1]];
    let v = [
        a[0] + b[0] - 2. * axis.origin[0],
        a[1] + b[1] - 2. * axis.origin[1],
        a[2] + b[2] - 2. * axis.origin[2],
    ];
    let norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
    let mid = std::array::from_fn(|i| axis.origin[i] + radius_mm * v[i] / norm);
    let vertex = input.vertices.len();
    input.vertices.push(mid);
    input.edges[index].vertices[1] = vertex;
    let added = input.edges.len();
    input.edges.push(multicad_rust_kernel::radial_tensor::Edge {
        vertices: [vertex, edge.vertices[1]],
        ..edge
    });
    for face in &mut input.faces {
        let mut ring = vec![];
        for &(i, forward) in &face.outer {
            if i != index {
                ring.push((i, forward));
            } else if forward {
                ring.extend([(index, true), (added, true)]);
            } else {
                ring.extend([(added, false), (index, false)]);
            }
        }
        face.outer = ring;
    }
    let refined = BoundaryOwnership::new(input, 1e-7, None).unwrap();
    // A valid finite manifold can contain extra rim vertices. Metadata must
    // reject ambiguous canonical roles before it could reach any body factory;
    // this control does not claim occupied-material certification for it.
    assert!(AuthoredPipe::new(source.clone(), refined.input()).is_err());
    let mut wrong = source.clone();
    wrong.root_faces[0] = wrong.band_faces[0];
    // Switching the private semantic declaration cannot make an uncovered
    // physical face disappear from the complete finite incidence inventory.
    assert!(AuthoredPipe::new(wrong, boundary.input()).is_err());
    let mut wrong = source.clone();
    wrong.root_faces[0] = source.cap_faces[0];
    assert!(AuthoredPipe::new(wrong, boundary.input()).is_err());
}

#[test]
fn actual_canonical_bspt_eight_mm_has_complete_owned_boundary() {
    use multicad_core::{document::Document, model, parameters};
    let document = Document::new("actual BSPT request").unwrap();
    let values = parameters::evaluate(&document).unwrap();
    let intent = model::PipeThreadIntent {
        form: model::PipeThreadForm::R,
        size: model::PipeThreadSize::OneEighth,
        gauge: model::PipeGaugeConvention::Iso7Nominal1994,
    };
    let mut parameters = intent.parameters();
    parameters.length = model::length(8.);
    let resolved = parameters
        .resolve(&values, model::ThreadMode::NewExternal)
        .unwrap();
    let spec = AuthoredSpec {
        major_radius_mm: resolved.major_radius_mm,
        pitch_mm: resolved.pitch_mm,
        length_mm: resolved.length_mm,
        radial_taper: resolved.profile.radial_taper,
        section: AuthoredSection::BsptRounded {
            depth_mm: resolved.profile.depth_mm,
        },
        crest_radius_limit_mm: None,
    };
    let source = AuthoredThread::new(spec, 1e-7, 1e-6).unwrap();
    let raw = source
        .raw_boundary_ownership(1e-7, None)
        .expect("actual raw boundary");
    let body = raw
        .ownership
        .bind_surfaces(1e-7, None)
        .expect("actual support binding");
    let body = body
        .certify_axial(1e-7, 1e-6, None)
        .expect("actual axial chart");
    let body = body.certify_caps(None).expect("actual caps");
    let body = body.certify_ordering(None).expect("actual radial order");
    let body = body.certify_material(1e-7, None).expect("actual material");
    assert!(body.volume_with_interrupt(1e-5, None).unwrap().volume_mm3[0] > 0.);
    let solid = Solid::rust_authored_axial_pipe_thread_with_interrupt(
        spec,
        EntityId::new(),
        PlaneFrame::default(),
        GeometryPolicy::default(),
        None,
    )
    .unwrap();
    let mesh = solid
        .manufacturing_mesh_with_receipt(MeshOptions {
            linear_deflection_mm: 0.1,
            angular_deflection_rad: 0.3,
            maximum_triangles: 1_000_000,
        })
        .unwrap();
    assert!(mesh.edge_pairs > 0 && mesh.vertex_links > 0);
    assert_eq!(mesh.removed_exact_zero_triangles, 0);
}
