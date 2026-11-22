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
                            - x.portable_hypot(y) * source.cone().depth_per_radius;
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
                                    - x.portable_hypot(y) * source.cone().depth_per_radius;
                                assert!((z - expected).abs() <= 4. * screw.source_error_mm());
                            }
                            let center = std::array::from_fn::<_, 3, _>(|i| {
                                (p[0][i] + p[1][i] + p[2][i]) / 3.
                            });
                            let cone_z = top + source.cone().tip_depth_mm
                                - center[0].portable_hypot(center[1])
                                    * source.cone().depth_per_radius;
                            assert!((center[2] - cone_z).abs() < opts.deflection);
                        } else if (wall..wall + 12).contains(&face) {
                            let arc = source.arcs()[(face - wall) as usize];
                            for point in p {
                                let radius = (point[0] - arc.center_mm[0])
                                    .portable_hypot(point[1] - arc.center_mm[1]);
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

#[test]
fn cylindrical_socket_caps_keep_true_torx_cavity_and_full_original_threads() {
    for left in [false, true] {
        for plain in [0., 3.] {
            let prototype = screw(1., left, false, true);
            let actual = Screw::new_with_recess(
                *prototype.thread(),
                1.5,
                Head::Cylindrical {
                    radius_mm: 2.8,
                    height_mm: 2.4,
                },
                plain,
                0.1,
                prototype.recess().clone(),
            )
            .unwrap();
            assert!(!actual.has_head_thread_overlap());
            assert_eq!(
                actual.visible_thread_interval_mm(),
                [0., actual.thread().spec().length_mm]
            );
            let mesh = tessellate_profiled_screw(&actual, options(1.)).unwrap();
            let receipt = certify(&mesh).unwrap();
            let bound = profiled_screw_volume_error_bound(&actual, options(1.)).unwrap();
            assert!((receipt.signed_volume_mm3 - actual.volume_mm3()).abs() <= bound);
            let vertices = mesh.vertices.as_chunks::<3>().0;
            let drive = actual.recess().hexalobular_geometry().unwrap();
            for arc in drive.arcs() {
                let [x, y] = arc.start_mm;
                assert!(vertices.contains(&[x, y, actual.head_top_z()]));
                assert!(vertices.contains(&[
                    x,
                    y,
                    actual.head_top_z() + drive.cone().tip_depth_mm
                        - x.portable_hypot(y) * drive.cone().depth_per_radius
                ]));
            }
            let expected = actual.thread().volume_mm3()
                + PI * (2.8_f64.powi(2) * 2.4 + 1.5_f64.powi(2) * plain)
                - actual.recess().nominal_volume_mm3();
            assert!((actual.volume_mm3() - expected).abs() < 1e-10);
        }
    }
}

#[test]
fn cosmetic_profiled_torx_retains_real_head_and_conical_drive_with_one_smooth_shaft() {
    for scale in [0.125, 1., 8.] {
        for left in [false, true] {
            for union in [false, true] {
                let source = screw(scale, left, union, false);
                let before = source.clone();
                let mesh = tessellate_profiled_cosmetic_screw(&source, options(scale)).unwrap();
                let proof = certify(&mesh).unwrap();
                let r = source.thread().spec().major_radius_mm;
                let interval = source.visible_thread_interval_mm();
                let expected = source.volume_mm3() - source.visible_thread_volume_mm3()
                    + PI * r * r * (interval[1] - interval[0]);
                assert!(
                    (proof.signed_volume_mm3 - expected).abs() < 0.1 * scale.powi(3),
                    "{} != {}",
                    proof.signed_volume_mm3,
                    expected
                );
                let faces = source.faces();
                let body = faces
                    .iter()
                    .position(|f| *f == Face::ThreadWall(0))
                    .unwrap() as u32;
                for (triangle, face) in mesh
                    .indices
                    .as_chunks::<3>()
                    .0
                    .iter()
                    .zip(&mesh.face_indices)
                {
                    if *face == body {
                        for vertex in triangle {
                            let [x, y, z] = mesh.vertices.as_chunks::<3>().0[*vertex as usize];
                            assert!((x.portable_hypot(y) - r).abs() < 1e-12 * scale);
                            assert!(z == interval[0] || z == interval[1]);
                        }
                    }
                }
                for arc in source.recess().hexalobular_geometry().unwrap().arcs() {
                    assert!(mesh.vertices.as_chunks::<3>().0.contains(&[
                        arc.start_mm[0],
                        arc.start_mm[1],
                        source.head_top_z()
                    ]));
                }
                assert_eq!(source, before);
            }
        }
    }
}

#[test]
fn cosmetic_profiled_source_has_its_own_cylinder_volume_faces_and_queries() {
    use crate::cosmetic_profiled_screw::Envelope;
    for union in [false, true] {
        let source = screw(1., false, union, false);
        let before = source.clone();
        let body = Envelope::new(source);
        assert_eq!(body.head_source(), &before);
        assert_eq!(body, body.clone());
        let interval = body.shaft_interval_mm();
        let r = body.shaft_radius_mm();
        assert!(
            !body
                .faces()
                .iter()
                .any(|f| matches!(f,Face::ThreadWall(i) if *i!=0))
        );
        assert_eq!(
            body.faces().contains(&Face::ThreadShoulder),
            before.profile().end()[0] != r
        );
        let mesh = body.mesh(options(1.)).unwrap();
        let proof = certify(&mesh).unwrap();
        assert!((proof.signed_volume_mm3 - body.volume_mm3()).abs() < 0.1);
        assert!((mesh.face_indices.iter().copied().max().unwrap() as usize) < body.faces().len());
        assert_eq!(body.support([0., 0., 1.]), interval[1]);
        assert_eq!(body.support([0., 0., -1.]), -before.head_top_z());
        assert_eq!(body.support([0.; 3]), 0.);
        assert!(body.support([f64::NAN, 0., 0.]).is_nan());
        assert!(!body.contains([f64::NAN, 0., 0.], 0.));
        assert!(!body.contains([0.; 3], -1.));
        assert!(!body.contains([0.; 3], f64::NAN));
        assert!(!body.contains([0., 0., before.head_top_z()], 0.));
        assert!(!body.contains([r * 2., 0., interval[1]], 0.));
        assert!(!body.contains([0., 0., interval[1] + 1.], 0.));
        let z = interval[0].midpoint(interval[1]);
        let mut fills_old_groove = false;
        for i in 0..128 {
            let a = TAU * i as f64 / 128.;
            let p = [r * 0.99 * a.portable_cos(), r * 0.99 * a.portable_sin(), z];
            assert!(body.contains(p, 0.));
            fills_old_groove |= !before.contains(p, 0.);
        }
        assert!(fills_old_groove);
        assert!(
            body.mesh(Tessellation {
                triangle_limit: 2,
                ..options(1.)
            })
            .is_err()
        );
        assert!(
            body.mesh(Tessellation {
                angle: 0.,
                ..options(1.)
            })
            .is_err()
        );
    }
}
