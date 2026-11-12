use super::*;
use crate::{
    head_profile::Head,
    helical_thread::{Spec, Thread},
    hexalobular_drive::{Cone, TangentSpec},
    manufacturing::certify,
    profiled_screw::{Face, Screw},
    screw_recess::Recess,
};

fn drive(scale: f64) -> Drive {
    Drive::from_authored_tangent(
        TangentSpec {
            major_diameter_mm: 2.8 * scale,
            minor_diameter_mm: 2.05 * scale,
            outer_lobe_radius_mm: 0.28 * scale,
        },
        Cone {
            tip_depth_mm: (1.1 + 1.4 / 3_f64.sqrt()) * scale,
            depth_per_radius: 1. / 3_f64.sqrt(),
        },
    )
    .unwrap()
}

fn options(scale: f64) -> Tessellation {
    Tessellation {
        deflection: 0.02 * scale,
        angle: 0.25,
        triangle_limit: 500_000,
    }
}

fn screw(scale: f64, left: bool, union: bool, fitted: bool) -> Screw {
    let source = drive(scale);
    let recess = Recess::conical_hexalobular(source, 1e-7 * scale.powi(3), 100_000).unwrap();
    let thread = Thread::new(Spec {
        major_radius_mm: if fitted { 1.4635 } else { 1.5 } * scale,
        pitch_mm: 0.5 * scale,
        length_mm: 3. * scale,
        depth_mm: if fitted { 0.29908756729740676 } else { 0.2 } * scale,
        crest_width_mm: if fitted { 0.07131009168581062 } else { 0.0625 } * scale,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: if left { 2 } else { 1 },
        left_hand: left,
        root_radius_mm: fitted.then_some(0.07216878364870323 * scale),
    })
    .unwrap();
    Screw::new_with_recess(
        thread,
        1.5 * scale,
        Head::Crowned {
            radius_mm: 2.8 * scale,
            height_mm: 2.4 * scale,
            crown_radius_mm: 5. * scale,
            top_land_radius_mm: 2.05 * scale,
            neck_radius_mm: 0.125 * scale,
        },
        if union { 0. } else { 3. * scale },
        0.1 * scale,
        recess,
    )
    .unwrap()
}

#[test]
fn metric_torx_mesh_retains_every_circle_contact_and_actual_cone_across_thread_variants() {
    for scale in [0.125, 1., 8.] {
        for left in [false, true] {
            for union in [false, true] {
                for fitted in [false, true] {
                    let screw = screw(scale, left, union, fitted);
                    let thread = *screw.thread();
                    let source = screw.recess().hexalobular_geometry().unwrap();
                    let opts = options(scale);
                    let mesh = tessellate_profiled_screw(&screw, opts).unwrap();
                    let proof = certify(&mesh).unwrap();
                    let bound = profiled_screw_volume_error_bound(&screw, opts).unwrap();
                    assert!((proof.signed_volume_mm3 - screw.volume_mm3()).abs() <= bound);
                    assert_eq!(*screw.thread(), thread);
                    let vertices = mesh.vertices.as_chunks::<3>().0;
                    let top = screw.head_top_z();
                    for arc in source.arcs() {
                        let [x, y] = arc.start_mm;
                        assert!(vertices.contains(&[x, y, top]));
                        let z = top + source.cone().tip_depth_mm
                            - x.hypot(y) * source.cone().depth_per_radius;
                        assert!(vertices.contains(&[x, y, z]));
                    }
                    let faces = screw.faces();
                    let wall = faces
                        .iter()
                        .position(|&f| f == Face::SocketWall(0))
                        .unwrap() as u32;
                    let floor = faces.iter().position(|&f| f == Face::SocketFloor).unwrap() as u32;
                    for index in 0..12 {
                        assert!(mesh.face_indices.contains(&(wall + index)));
                    }
                    for (triangle, &face) in mesh
                        .indices
                        .as_chunks::<3>()
                        .0
                        .iter()
                        .zip(&mesh.face_indices)
                    {
                        let p = triangle.map(|i| vertices[i as usize]);
                        if face == floor {
                            // The varying rim and central apex lie on the actual cone.
                            for [x, y, z] in p {
                                let expected = top + source.cone().tip_depth_mm
                                    - x.hypot(y) * source.cone().depth_per_radius;
                                assert!((z - expected).abs() <= 4. * screw.source_error_mm());
                            }
                            let center = std::array::from_fn::<_, 3, _>(|i| {
                                (p[0][i] + p[1][i] + p[2][i]) / 3.
                            });
                            let cone_z = top + source.cone().tip_depth_mm
                                - center[0].hypot(center[1]) * source.cone().depth_per_radius;
                            assert!((center[2] - cone_z).abs() < opts.deflection);
                        } else if (wall..wall + 12).contains(&face) {
                            let arc = source.arcs()[(face - wall) as usize];
                            for point in p {
                                let radius = (point[0] - arc.center_mm[0])
                                    .hypot(point[1] - arc.center_mm[1]);
                                assert!(
                                    (radius - arc.radius_mm).abs() <= 4. * source.source_error_mm()
                                );
                            }
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn torx_budget_precision_and_unresolved_geometry_refuse_without_changing_source() {
    let screw = screw(1., false, true, true);
    let source = screw.recess().hexalobular_geometry().unwrap();
    let snapshot = source.arcs().to_owned();
    let opts = options(1.);
    assert_eq!(
        stations(
            source,
            Tessellation {
                triangle_limit: 1,
                ..opts
            }
        ),
        Err(MeshError::BudgetExceeded)
    );
    assert_eq!(
        stations(
            source,
            Tessellation {
                deflection: f64::MIN_POSITIVE,
                ..opts
            }
        ),
        Err(MeshError::BudgetExceeded)
    );
    assert_eq!(
        tessellate_profiled_screw(
            &screw,
            Tessellation {
                deflection: screw.source_error_mm(),
                ..opts
            }
        ),
        Err(MeshError::UnresolvedBoundary)
    );
    assert_eq!(
        tessellate_profiled_screw(
            &screw,
            Tessellation {
                triangle_limit: 1,
                ..opts
            }
        ),
        Err(MeshError::BudgetExceeded)
    );
    let bad = [[1., 0., 0.], [1., 0., 0.], [0., 1., 0.]];
    assert_eq!(
        validate_ring(&bad, source, opts),
        Err(MeshError::UnresolvedBoundary)
    );
    let invalid = [[f64::NAN, 0., 0.]; 3];
    assert_eq!(
        validate_ring(&invalid, source, opts),
        Err(MeshError::UnresolvedBoundary)
    );
    assert_eq!(
        volume_error(source, opts, f64::MAX, 1.),
        Err(MeshError::UnresolvedBoundary)
    );
    assert_eq!(
        profiled_screw_volume_error_bound(
            &screw,
            Tessellation {
                deflection: f64::MAX,
                ..opts
            }
        ),
        Err(MeshError::UnresolvedBoundary)
    );
    assert_eq!(source.arc_index_at_angle(f64::NAN), None);
    assert_eq!(source.boundary_volume_error_bound_mm3(-1., 0.), None);
    assert_eq!(
        source.boundary_volume_error_bound_mm3(0., f64::INFINITY),
        None
    );
    assert_eq!(
        source.boundary_volume_error_bound_mm3(f64::MAX, f64::MAX),
        None
    );
    assert!(source.boundary_volume_error_bound_mm3(0., 0.).unwrap() >= 0.);
    assert_eq!(source.arcs(), &snapshot);
    assert!(certify(&tessellate_profiled_screw(&screw, opts).unwrap()).is_ok());
}

#[test]
fn existing_socket_wall_budget_failure_keeps_the_original_source_reusable() {
    use crate::headed_screw::{Screw as Legacy, ScrewFace, ScrewHead, ScrewSpec};
    let thread = *screw(1., false, false, false).thread();
    let source = Legacy::new(
        thread,
        ScrewSpec {
            head: ScrewHead::SocketCap {
                radius_mm: 3.,
                height_mm: 3.,
                socket_across_flats_mm: 2.,
                socket_depth_mm: 1.,
            },
            shank_length_mm: 1.,
            shank_radius_mm: 1.5,
        },
    )
    .unwrap();
    let opts = options(1.);
    let mesh = super::super::tessellate_screw(&source, opts).unwrap();
    let first_wall = source
        .faces()
        .iter()
        .position(|&f| f == ScrewFace::SocketWall(0))
        .unwrap() as u32;
    let first_triangle = mesh
        .face_indices
        .iter()
        .position(|&f| f >= first_wall)
        .unwrap();
    assert_eq!(
        super::super::tessellate_screw(
            &source,
            Tessellation {
                triangle_limit: first_triangle + 1,
                ..opts
            }
        ),
        Err(MeshError::BudgetExceeded)
    );
    assert_eq!(super::super::tessellate_screw(&source, opts).unwrap(), mesh);
}
