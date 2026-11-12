use super::*;
use crate::{
    display::Tessellation,
    helical_thread,
    helical_thread_mesh::{profiled_screw_volume_error_bound, tessellate_profiled_screw},
    manufacturing::certify,
};

fn thread(starts: u32, left_hand: bool) -> Thread {
    Thread::new(helical_thread::Spec {
        major_radius_mm: 1.5,
        pitch_mm: 0.5,
        length_mm: 7.5,
        depth_mm: 0.2,
        crest_width_mm: 0.0625,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts,
        left_hand,
        root_radius_mm: None,
    })
    .unwrap()
}

fn spec(crowned: bool) -> Spec {
    Spec {
        head: if crowned {
            Head::Crowned {
                radius_mm: 3.,
                height_mm: 2.5,
                crown_radius_mm: 4.,
                top_land_radius_mm: 1.5,
                neck_radius_mm: 0.1,
            }
        } else {
            Head::Countersunk {
                radius_mm: 3.,
                height_mm: 1.65,
                included_angle_deg: 90.,
                edge_chamfer_mm: 0.1,
                neck_radius_mm: 0.1,
            }
        },
        shank_length_mm: 1.,
        socket_across_flats_mm: 2.,
        socket_depth_mm: 0.6,
        minimum_wall_mm: 0.1,
    }
}

fn options() -> Tessellation {
    Tessellation {
        deflection: 0.02,
        angle: 0.25,
        triangle_limit: 500_000,
    }
}

#[test]
fn real_head_neck_drive_and_unchanged_thread_have_one_material_domain() {
    for crowned in [false, true] {
        let thread = thread(1, false);
        let spec = spec(crowned);
        let screw = Screw::new(thread, spec).unwrap();
        assert_eq!(*screw.thread(), thread);
        assert_eq!(screw.spec(), spec);
        assert_eq!(screw.profile().end(), [1.5, 1.]);
        let top = screw.head_top_z();
        assert!(!screw.contains([0., 0., top + 0.1], 0.));
        assert!(screw.contains([0., 0., top + 0.6], 0.));
        assert!(screw.contains([1.4, 0., top + 0.1], 0.));
        assert!(!screw.contains([3.1, 0., -1.5], 0.));
        assert!(screw.contains([1.4, 0., -0.5], 0.));
        assert!(screw.contains([0., 0., 0.5], 0.));
        assert!(!screw.contains([0., 0., top - 0.1], 0.));
        assert!(!screw.contains([0., 0., 8.], 0.));
        assert!(screw.contains([1.5001, 0., -0.5], 0.001));
        for t in [-1., f64::NAN, f64::INFINITY] {
            assert!(!screw.contains([0.; 3], t));
        }
        assert!(!screw.contains([f64::NAN, 0., 0.], 0.));
        assert_eq!(screw.support([0.; 3]), 0.);
        assert!(screw.support([f64::NAN, 0., 0.]).is_nan());
        assert!((screw.support([0., 0., -1.]) + top).abs() < 1e-12);
        assert_eq!(screw.support([0., 0., 1.]), 7.5);
        assert_eq!(screw.support([1., 0., 0.]), 3.);
        assert_eq!(screw.faces().len(), 19);
        assert!(screw.volume_mm3() > thread.volume_mm3());
    }
    let t = thread(1, false);
    for s in [
        Spec {
            socket_across_flats_mm: 1e-300,
            ..spec(false)
        },
        Spec {
            socket_depth_mm: 1e-300,
            ..spec(false)
        },
        Spec {
            socket_across_flats_mm: 8.,
            ..spec(false)
        },
        Spec {
            socket_depth_mm: 1.7,
            ..spec(false)
        },
        Spec {
            shank_length_mm: 0.,
            ..spec(false)
        },
    ] {
        assert_eq!(Screw::new(t, s), Err(ThreadError::Invalid));
    }
}

#[test]
fn sampled_real_profiles_share_the_original_shoulder_and_a_closed_oriented_shell() {
    for crowned in [false, true] {
        for (starts, left) in [(1, false), (1, true), (2, false)] {
            let screw = Screw::new(thread(starts, left), spec(crowned)).unwrap();
            let mesh = tessellate_profiled_screw(&screw, options()).unwrap();
            let proof = certify(&mesh).unwrap();
            assert!(proof.signed_volume_mm3 > 0.);
            let bound = profiled_screw_volume_error_bound(&screw, options()).unwrap();
            assert!((proof.signed_volume_mm3 - screw.volume_mm3()).abs() < 0.5);
            assert!((proof.signed_volume_mm3 - screw.volume_mm3()).abs() <= bound);
            let roles: std::collections::BTreeSet<_> = mesh.face_indices.iter().copied().collect();
            assert_eq!(roles.len(), screw.faces().len());
        }
    }
}

#[test]
fn malformed_options_and_excess_profile_refinement_fail_without_poisoning_source() {
    use crate::display::MeshError;
    let screw = Screw::new(thread(1, false), spec(false)).unwrap();
    assert_eq!(
        profiled_screw_volume_error_bound(
            &screw,
            Tessellation {
                deflection: 0.,
                ..options()
            }
        ),
        Err(MeshError::InvalidOptions)
    );
    assert_eq!(
        profiled_screw_volume_error_bound(
            &screw,
            Tessellation {
                triangle_limit: 1,
                ..options()
            }
        ),
        Err(MeshError::BudgetExceeded)
    );
    assert_eq!(
        tessellate_profiled_screw(
            &screw,
            Tessellation {
                deflection: 0.,
                ..options()
            }
        ),
        Err(MeshError::InvalidOptions)
    );
    assert_eq!(
        tessellate_profiled_screw(
            &screw,
            Tessellation {
                triangle_limit: 1,
                ..options()
            }
        ),
        Err(MeshError::BudgetExceeded)
    );
    assert_eq!(
        tessellate_profiled_screw(
            &screw,
            Tessellation {
                triangle_limit: 300,
                ..options()
            }
        ),
        Err(MeshError::BudgetExceeded)
    );
    let t = Thread::new(helical_thread::Spec {
        length_mm: 0.5,
        ..thread(1, false).spec()
    })
    .unwrap();
    let crown = Screw::new(t, spec(true)).unwrap();
    assert_eq!(
        profiled_screw_volume_error_bound(
            &crown,
            Tessellation {
                angle: 0.08,
                triangle_limit: 12_000,
                ..options()
            }
        ),
        Err(MeshError::BudgetExceeded)
    );
    assert_eq!(
        tessellate_profiled_screw(
            &crown,
            Tessellation {
                angle: 0.08,
                triangle_limit: 12_000,
                ..options()
            }
        ),
        Err(MeshError::BudgetExceeded)
    );
    let mesh = tessellate_profiled_screw(&screw, options()).unwrap();
    // Exercise exhaustion throughout the thread, head, drive and floor;
    // every failed attempt leaves the same retained source reusable.
    for triangle_limit in (128..mesh.indices.len() / 3).step_by(127) {
        assert_eq!(
            tessellate_profiled_screw(
                &screw,
                Tessellation {
                    triangle_limit,
                    ..options()
                }
            ),
            Err(MeshError::BudgetExceeded)
        );
    }
    assert!(certify(&tessellate_profiled_screw(&screw, options()).unwrap()).is_ok());
}

#[test]
fn profile_triangles_track_real_cone_crown_neck_and_bearing_geometry() {
    use crate::head_profile::Curve;
    let distance = |a: [f64; 2], b: [f64; 2]| (a[0] - b[0]).hypot(a[1] - b[1]);
    for crowned in [false, true] {
        let screw = Screw::new(thread(1, false), spec(crowned)).unwrap();
        let mesh = tessellate_profiled_screw(&screw, options()).unwrap();
        let vertices = mesh.vertices.as_chunks::<3>().0;
        let faces = screw.faces();
        let segments: Vec<_> = screw.profile().segments().collect();
        for (triangle, &face) in mesh
            .indices
            .as_chunks::<3>()
            .0
            .iter()
            .zip(&mesh.face_indices)
        {
            let Face::ProfileWall(index) = faces[face as usize] else {
                continue;
            };
            let segment = segments[index as usize];
            let [a, b, c] = triangle.map(|i| vertices[i as usize]);
            let center: [f64; 3] = std::array::from_fn(|i| (a[i] + b[i] + c[i]) / 3.);
            let p = [
                center[0].hypot(center[1]),
                center[2] + screw.spec().shank_length_mm,
            ];
            let d = match segment.curve {
                Curve::Line => {
                    let delta: [f64; 2] =
                        std::array::from_fn(|i| segment.end[i] - segment.start[i]);
                    let t = ((p[0] - segment.start[0]) * delta[0]
                        + (p[1] - segment.start[1]) * delta[1])
                        / (delta[0].powi(2) + delta[1].powi(2));
                    distance(p, segment.point(t.clamp(0., 1.)))
                }
                Curve::Arc { center, sweep } => {
                    let angle = (p[1] - center[1]).atan2(p[0] - center[0]);
                    let start = (segment.start[1] - center[1]).atan2(segment.start[0] - center[0]);
                    let advance = if sweep > 0. {
                        (angle - start).rem_euclid(std::f64::consts::TAU)
                    } else {
                        (start - angle).rem_euclid(std::f64::consts::TAU)
                    };
                    if advance <= sweep.abs() {
                        (distance(p, center) - distance(segment.start, center)).abs()
                    } else {
                        distance(p, segment.start).min(distance(p, segment.end))
                    }
                }
            };
            assert!(d <= options().deflection, "{:?}: chord {d}", segment.role);
            let u: [f64; 3] = std::array::from_fn(|i| b[i] - a[i]);
            let v: [f64; 3] = std::array::from_fn(|i| c[i] - a[i]);
            let normal = [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ];
            let norm = normal[0].hypot(normal[1]).hypot(normal[2]);
            for vertex in [a, b, c] {
                let [radial, axial] = match segment.curve {
                    Curve::Line => [
                        segment.end[1] - segment.start[1],
                        segment.start[0] - segment.end[0],
                    ],
                    Curve::Arc { center, sweep } => [
                        (vertex[0].hypot(vertex[1]) - center[0]) * sweep.signum(),
                        (vertex[2] + screw.spec().shank_length_mm - center[1]) * sweep.signum(),
                    ],
                };
                let theta = vertex[1].atan2(vertex[0]);
                let target = [radial * theta.cos(), radial * theta.sin(), axial];
                let dot: f64 = normal.into_iter().zip(target).map(|(a, b)| a * b).sum();
                let angle = (dot / (norm * radial.hypot(axial))).clamp(-1., 1.).acos();
                assert!(
                    angle <= options().angle,
                    "{:?}: normal {angle}",
                    segment.role
                );
            }
        }
    }
}
