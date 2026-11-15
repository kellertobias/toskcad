use super::*;
use multicad_core::references::TopologyKind;
use multicad_rust_kernel::{
    display::DisplayMesh,
    head_profile::Head,
    hexalobular_drive::{Cone, Drive, TangentSpec},
    manufacturing::certify,
    profiled_screw::SourceSpec,
    screw_recess::Recess,
};
fn source(head: Head, frame: PlaneFrame, left: bool) -> Solid {
    let drive = Drive::from_authored_tangent(
        TangentSpec {
            major_diameter_mm: 2.8,
            minor_diameter_mm: 2.05,
            outer_lobe_radius_mm: 0.28,
        },
        Cone {
            tip_depth_mm: 0.5 + 1.4 / 3f64.sqrt(),
            depth_per_radius: 1. / 3f64.sqrt(),
        },
    )
    .unwrap();
    Solid::rust_profiled_machine_screw_with_recess(
        crate::HelicalThread {
            major_radius_mm: 1.5,
            pitch_mm: 0.5,
            length_mm: 5.,
            depth_mm: 0.2,
            crest_width_mm: 0.0625,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            radial_taper: 0.,
            starts: 1,
            left_hand: left,
        },
        None,
        1.5,
        SourceSpec {
            head: Some(head),
            shank_length_mm: 3.,
            minimum_wall_mm: 0.1,
        },
        Recess::conical_hexalobular(drive, 1e-7, 100_000).unwrap(),
        EntityId::new(),
        frame,
        GeometryPolicy::default(),
    )
    .unwrap()
}
#[test]
fn cosmetic_profiled_native_preserves_actual_head_faces_refs_and_closed_display_material() {
    for head in [
        Head::Cylindrical {
            radius_mm: 2.8,
            height_mm: 2.4,
        },
        Head::Countersunk {
            radius_mm: 2.8,
            height_mm: 1.8,
            included_angle_deg: 90.,
            edge_chamfer_mm: 0.1,
            neck_radius_mm: 0.,
        },
        Head::Crowned {
            radius_mm: 2.8,
            height_mm: 2.4,
            crown_radius_mm: 5.,
            top_land_radius_mm: 2.05,
            neck_radius_mm: 0.125,
        },
    ] {
        for left in [false, true] {
            let frame = PlaneFrame {
                origin: [17., -8., 13.],
                u: [0.8, 0.6, 0.],
                v: [-0.36, 0.48, 0.8],
            };
            let original = source(head, frame, left);
            let volume = original.volume_mm3().unwrap();
            let cosmetic = original.cosmetic_profiled_screw_envelope().unwrap();
            assert_eq!(original.volume_mm3().unwrap(), volume);
            assert!(cosmetic.volume_mm3().unwrap() > volume);
            assert_eq!(cosmetic.solid_count().unwrap(), 1);
            let before = original.topology_snapshot(EntityId::new(), 0).unwrap();
            let after = cosmetic.topology_snapshot(EntityId::new(), 0).unwrap();
            for element in before.elements.iter().filter(|e| {
                e.origins.iter().any(|o| {
                    (o.role.starts_with("screw.head/") || o.role.starts_with("screw.drive/"))
                        && !o.role.starts_with("screw.head/thread.junction/")
                })
            }) {
                let new = after
                    .elements
                    .iter()
                    .find(|e| e.key.kind == element.key.kind && e.origins == element.origins)
                    .unwrap_or_else(|| panic!("lost head reference {:?}", element.origins));
                assert_eq!(new.geometry, element.geometry);
                assert_eq!(new.frames, element.frames);
            }
            let mesh = cosmetic.mesh(Default::default()).unwrap();
            assert!(
                certify(&DisplayMesh {
                    vertices: mesh.vertices,
                    indices: mesh.indices,
                    face_indices: mesh.face_indices
                })
                .is_ok()
            );
            let edges = cosmetic.display_edges(0.01).unwrap();
            assert_eq!(
                edges.len(),
                after
                    .elements
                    .iter()
                    .filter(|e| e.key.kind == TopologyKind::Edge)
                    .count()
            );
            for face in after
                .elements
                .iter()
                .filter(|e| e.key.kind == TopologyKind::Face && e.geometry.is_some())
            {
                let Some(AnalyticGeometry::Plane { frame: plane }) = face.geometry else {
                    unreachable!()
                };
                let plane = PlaneFrame {
                    origin: plane.origin_mm,
                    u: plane.u,
                    v: plane.v,
                };
                assert!(
                    !cosmetic
                        .planar_face_contours(face.key.ordinal, plane)
                        .unwrap()
                        .loops
                        .is_empty()
                );
            }
            assert!(!cosmetic.supports_step_export());
            assert_eq!(
                cosmetic
                    .manufacturing_mesh(Default::default())
                    .err()
                    .unwrap()
                    .kind,
                GeometryErrorKind::Unsupported
            );
            let shared = cosmetic.share().unwrap();
            assert_eq!(shared.volume_mm3().unwrap(), cosmetic.volume_mm3().unwrap());
        }
    }
}

#[test]
fn cosmetic_profiled_refuses_unrelated_sources_and_unresolvable_edge_requests() {
    use multicad_rust_kernel::display::{MeshError, Tessellation};
    let unrelated = Solid::box_solid([3., 4., 5.], GeometryPolicy::default()).unwrap();
    assert_eq!(
        unrelated
            .cosmetic_profiled_screw_envelope()
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Unsupported
    );
    let original = source(
        Head::Cylindrical {
            radius_mm: 2.8,
            height_mm: 2.4,
        },
        PlaneFrame::default(),
        false,
    );
    let cosmetic = original.cosmetic_profiled_screw_envelope().unwrap();
    let (body, frame) = cosmetic.cosmetic_profiled_body_ref().unwrap();
    let options = Tessellation {
        deflection: 0.05,
        angle: 0.3,
        triangle_limit: 100_000,
    };
    let floor = body
        .boundaries()
        .iter()
        .position(|b| matches!(b.edge, Edge::Floor(_)))
        .unwrap() as u32;
    assert!(body.planar_loops(u32::MAX).is_empty());
    assert!(body.sample_edge(frame, u32::MAX, options).is_none());
    assert!(body.sample_edge(frame, 0, options).is_none());
    assert_eq!(
        body.sample_edge(
            frame,
            floor,
            Tessellation {
                deflection: 0.,
                ..options
            }
        )
        .unwrap(),
        Err(MeshError::InvalidOptions)
    );
    assert_eq!(
        body.sample_edge(
            frame,
            floor,
            Tessellation {
                triangle_limit: 1,
                ..options
            }
        )
        .unwrap(),
        Err(MeshError::BudgetExceeded)
    );
    let remote = PlaneFrame {
        origin: [1e16, 0., 0.],
        ..frame
    };
    assert_eq!(
        body.sample_edge(remote, floor, options).unwrap(),
        Err(MeshError::UnresolvedBoundary)
    );
    assert_eq!(
        original.volume_mm3().unwrap(),
        source(
            Head::Cylindrical {
                radius_mm: 2.8,
                height_mm: 2.4
            },
            PlaneFrame::default(),
            false
        )
        .volume_mm3()
        .unwrap()
    );
}
