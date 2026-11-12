use super::*;
use crate::headed_screw::hex_area;
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
fn retained_drive_orientation_controls_material_and_conforming_manufacturing_mesh() {
    for crowned in [false, true] {
        for (starts, left) in [(1, false), (2, true)] {
            for union in [false, true] {
                let mut authored = spec(crowned);
                if union {
                    authored.shank_length_mm = 0.;
                }
                let original = if union {
                    Screw::new_with_head_thread_union(thread(starts, left), authored)
                } else {
                    Screw::new(thread(starts, left), authored)
                }
                .unwrap();
                let mut physical = original.spec();
                physical.head = None;
                for angle in [
                    0.37,
                    -std::f64::consts::PI / 6.,
                    -0.03125 / original.thread().helical_rate() - std::f64::consts::PI / 6.,
                ] {
                    let (s, c) = angle.sin_cos();
                    let points = crate::headed_screw::hex_vertices(2.)
                        .map(|[x, y]| [c * x - s * y + 1e-14, s * x + c * y]);
                    let retained = if union {
                        Screw::from_source_with_union_and_drive(
                            *original.thread(),
                            *original.profile(),
                            physical,
                            points,
                        )
                    } else {
                        Screw::from_source_with_drive(
                            *original.thread(),
                            *original.profile(),
                            physical,
                            points,
                        )
                    }
                    .unwrap();
                    assert_eq!(retained.drive().corners_mm(), points);
                    assert_eq!(*retained.thread(), *original.thread());
                    assert_eq!(
                        retained.visible_thread_interval_mm(),
                        original.visible_thread_interval_mm()
                    );
                    assert!(!retained.drive().is_canonical());
                    let mesh = tessellate_profiled_screw(&retained, options()).unwrap();
                    assert_eq!(
                        tessellate_profiled_screw(
                            &retained,
                            Tessellation {
                                deflection: 1.,
                                angle: 1.,
                                triangle_limit: 48
                            }
                        ),
                        Err(crate::display::MeshError::BudgetExceeded)
                    );
                    let receipt = certify(&mesh).unwrap();
                    assert!(
                        (receipt.signed_volume_mm3 - retained.volume_mm3()).abs()
                            <= profiled_screw_volume_error_bound(&retained, options()).unwrap()
                    );
                    for point in points {
                        for z in [
                            retained.head_top_z(),
                            retained.head_top_z() + retained.spec().socket_depth_mm,
                        ] {
                            assert!(
                                mesh.vertices
                                    .as_chunks::<3>()
                                    .0
                                    .contains(&[point[0], point[1], z])
                            );
                        }
                    }
                    let floor = retained.head_top_z() + retained.spec().socket_depth_mm;
                    assert!(!retained.contains([0., 0., floor - 0.1], 0.));
                    assert!(retained.contains([0., 0., floor + 0.1], 0.));
                    let mut bad = points;
                    bad[0][0] += 0.01;
                    assert!(
                        Screw::from_source_with_drive(
                            *original.thread(),
                            *original.profile(),
                            physical,
                            bad
                        )
                        .is_err()
                    );
                }
                let corners = original.drive().corners_mm();
                assert!(
                    Screw::from_source_with_drive(
                        *original.thread(),
                        *original.profile(),
                        original.spec(),
                        corners
                    )
                    .is_err()
                );
                assert!(
                    Screw::from_source_with_union_and_drive(
                        *original.thread(),
                        *original.profile(),
                        original.spec(),
                        corners
                    )
                    .is_err()
                );
            }
        }
    }
}

#[test]
fn real_head_neck_drive_and_unchanged_thread_have_one_material_domain() {
    for crowned in [false, true] {
        let thread = thread(1, false);
        let spec = spec(crowned);
        let screw = Screw::new(thread, spec).unwrap();
        assert_eq!(*screw.thread(), thread);
        assert_eq!(screw.spec(), spec.into());
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
    for deflection in [
        screw.profile().rounding_error_mm() / 2.,
        4. * screw.profile().rounding_error_mm(),
    ] {
        let tight = Tessellation {
            deflection,
            ..options()
        };
        assert_eq!(
            tessellate_profiled_screw(&screw, tight),
            Err(MeshError::UnresolvedBoundary)
        );
        assert_eq!(
            profiled_screw_volume_error_bound(&screw, tight),
            Err(MeshError::UnresolvedBoundary)
        );
    }
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
                Curve::Arc {
                    center,
                    radius_mm: radius,
                    sweep,
                } => {
                    let angle = (p[1] - center[1]).atan2(p[0] - center[0]);
                    let start = (segment.start[1] - center[1]).atan2(segment.start[0] - center[0]);
                    let advance = if sweep > 0. {
                        (angle - start).rem_euclid(std::f64::consts::TAU)
                    } else {
                        (start - angle).rem_euclid(std::f64::consts::TAU)
                    };
                    if advance <= sweep.abs() {
                        (distance(p, center) - radius).abs()
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
                    Curve::Arc { center, sweep, .. } => [
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

#[test]
fn physical_profile_and_drive_need_no_authored_head_recipe() {
    let original = Screw::new(thread(1, false), spec(false)).unwrap();
    let profile =
        Profile::from_retained_segments(&original.profile().segments().collect::<Vec<_>>())
            .unwrap();
    let mut data = original.spec();
    assert_eq!(
        Screw::from_source(*original.thread(), profile, data),
        Err(ThreadError::Invalid)
    );
    data.head = None;
    let physical = Screw::from_source(*original.thread(), profile, data).unwrap();
    assert_eq!(physical.spec().head, None);
    assert_eq!(physical.volume_mm3(), original.volume_mm3());
    assert_eq!(physical.faces(), original.faces());
    assert_eq!(
        physical.profile().segments().collect::<Vec<_>>(),
        original.profile().segments().collect::<Vec<_>>()
    );
    assert!(certify(&tessellate_profiled_screw(&physical, options()).unwrap()).is_ok());
    for bad in [
        SourceSpec {
            shank_length_mm: 0.,
            ..data
        },
        SourceSpec {
            socket_depth_mm: 100.,
            ..data
        },
    ] {
        assert_eq!(
            Screw::from_source(*original.thread(), profile, bad),
            Err(ThreadError::Invalid)
        );
    }
    let different = Thread::new(helical_thread::Spec {
        major_radius_mm: 1.6,
        ..thread(1, false).spec()
    })
    .unwrap();
    assert_eq!(
        Screw::from_source(different, profile, data),
        Err(ThreadError::Invalid)
    );
}

#[test]
fn opt_in_union_retains_full_thread_phase_and_actual_neck_material() {
    for major in [1.5, 3.] {
        for crowned in [false, true] {
            for (starts, left) in [(1, false), (1, true), (2, false), (2, true)] {
                for plain in [0., 0.01] {
                    let t = Thread::new(helical_thread::Spec {
                        major_radius_mm: major,
                        ..thread(starts, left).spec()
                    })
                    .unwrap();
                    let mut s = spec(crowned);
                    s.head = match s.head {
                        Head::Countersunk {
                            included_angle_deg,
                            edge_chamfer_mm,
                            neck_radius_mm,
                            ..
                        } => Head::Countersunk {
                            radius_mm: 2. * major,
                            height_mm: major + 0.15,
                            included_angle_deg,
                            edge_chamfer_mm,
                            neck_radius_mm,
                        },
                        Head::Crowned { neck_radius_mm, .. } => Head::Crowned {
                            radius_mm: 2. * major,
                            height_mm: major + 1.,
                            crown_radius_mm: 3. * major,
                            top_land_radius_mm: major,
                            neck_radius_mm,
                        },
                    };
                    s.shank_length_mm = plain;
                    assert!(Screw::new(t, s).is_err());
                    let screw = Screw::new_with_head_thread_union(t, s).unwrap();
                    assert_eq!(*screw.thread(), t);
                    assert_eq!(screw.spec(), s.into());
                    assert!(screw.has_head_thread_overlap());
                    let [lo, hi] = screw.visible_thread_interval_mm();
                    assert!(lo > 0. && hi == t.spec().length_mm);
                    assert!(screw.contains([major, 0., lo / 2.], 0.));
                    assert!(screw.visible_thread_volume_mm3() < t.volume_mm3());
                    let expected = screw.profile().volume_mm3() + t.volume_mm3() * (hi - lo) / hi
                        - hex_area(s.socket_across_flats_mm) * s.socket_depth_mm;
                    assert!((screw.volume_mm3() - expected).abs() < 1e-12);
                    let mesh = tessellate_profiled_screw(&screw, options()).unwrap();
                    let proof = certify(&mesh).unwrap();
                    assert!(
                        (proof.signed_volume_mm3 - screw.volume_mm3()).abs()
                            <= profiled_screw_volume_error_bound(&screw, options()).unwrap()
                    );
                    let faces = screw.faces();
                    let vertices = mesh.vertices.as_chunks::<3>().0;
                    let mut phase_changed = false;
                    for (triangle, &face) in mesh
                        .indices
                        .as_chunks::<3>()
                        .0
                        .iter()
                        .zip(&mesh.face_indices)
                    {
                        if let Face::ThreadWall(_) = faces[face as usize] {
                            for &id in triangle {
                                let [x, y, z] = vertices[id as usize];
                                assert!(z >= lo);
                                if z == lo {
                                    let theta = y.atan2(x);
                                    assert!((x.hypot(y) - t.radius(theta, lo)).abs() < 1e-12);
                                    phase_changed |=
                                        (t.radius(theta, lo) - t.radius(theta, 0.)).abs() > 1e-4;
                                }
                            }
                        }
                    }
                    assert!(phase_changed);
                    let mut source = screw.spec();
                    source.head = None;
                    assert!(Screw::from_source(t, *screw.profile(), source).is_err());
                    let imported =
                        Screw::from_source_with_union(t, *screw.profile(), source).unwrap();
                    assert_eq!(imported.visible_thread_interval_mm(), [lo, hi]);
                    assert_eq!(imported.volume_mm3(), screw.volume_mm3());
                    assert!(Screw::from_source_with_union(t, *screw.profile(), s.into()).is_err());
                }
            }
        }
    }
}

#[test]
fn union_refuses_real_crescents_and_exhausted_thread_without_mutating_retry() {
    let t = thread(1, false);
    let mut s = spec(true);
    s.shank_length_mm = 0.;
    s.head = Head::Crowned {
        radius_mm: 3.,
        height_mm: 2.5,
        crown_radius_mm: 4.,
        top_land_radius_mm: 1.5,
        neck_radius_mm: 0.2,
    };
    // Literal 1.7 - 0.2 is strictly smaller than literal 1.5.
    assert!(Screw::new_with_head_thread_union(t, s).is_err());
    assert!(
        Screw::new(
            t,
            Spec {
                shank_length_mm: 1.,
                ..s
            }
        )
        .is_ok()
    );
    s.head = Head::Crowned {
        radius_mm: 3.,
        height_mm: 2.5,
        crown_radius_mm: 4.,
        top_land_radius_mm: 1.5,
        neck_radius_mm: 0.125,
    };
    let good = Screw::new_with_head_thread_union(t, s).unwrap();
    let short_spec = helical_thread::Spec {
        pitch_mm: 0.05,
        depth_mm: 0.02,
        crest_width_mm: 0.00625,
        ..t.spec()
    };
    for length_mm in [0.1, 0.125, 0.125 + 1e-14] {
        let short = Thread::new(helical_thread::Spec {
            length_mm,
            ..short_spec
        })
        .unwrap();
        assert!(Screw::new_with_head_thread_union(short, s).is_err());
    }
    let short = Thread::new(helical_thread::Spec {
        length_mm: 0.1251,
        ..short_spec
    })
    .unwrap();
    assert!(Screw::new_with_head_thread_union(short, s).is_ok());
    let mut source = good.spec();
    source.head = None;
    for plain in [-1., f64::NAN, 0.25] {
        source.shank_length_mm = plain;
        assert!(Screw::from_source_with_union(t, *good.profile(), source).is_err());
    }
    let old = Screw::new(t, spec(false)).unwrap();
    assert_eq!(
        old,
        Screw::new_with_head_thread_union(t, spec(false)).unwrap()
    );
    assert!(!old.has_head_thread_overlap());
    assert_eq!(old.visible_thread_volume_mm3(), t.volume_mm3());
    assert!(certify(&tessellate_profiled_screw(&good, options()).unwrap()).is_ok());
}

#[test]
fn retained_union_proves_actual_cylindrical_overlap_and_bounds_complete_head() {
    let t = thread(1, false);
    let old = Screw::new(t, spec(false)).unwrap();
    let source = SourceSpec {
        head: None,
        shank_length_mm: 0.,
        ..old.spec()
    };
    let union = Screw::from_source_with_union(t, *old.profile(), source).unwrap();
    assert_eq!(union.visible_thread_interval_mm(), [1., 7.5]);
    assert!(union.contains([1.5, 0., 0.9], 0.));
    assert!(certify(&tessellate_profiled_screw(&union, options()).unwrap()).is_ok());
    assert!(
        Screw::new_with_head_thread_union(
            t,
            Spec {
                shank_length_mm: -1.,
                ..spec(false)
            }
        )
        .is_err()
    );
    assert_eq!(
        Profile::new_for_thread_union(
            Head::Crowned {
                radius_mm: 4.,
                height_mm: 9_999_999.,
                crown_radius_mm: 5.,
                top_land_radius_mm: 1.,
                neck_radius_mm: 2.,
            },
            1.,
            0.
        ),
        Err(ThreadError::Limit)
    );
}

#[test]
fn nominal_shank_preserves_fitted_thread_shoulders_round_roots_and_head_unions() {
    for crowned in [false, true] {
        for left in [false, true] {
            for starts in [1, 2] {
                for rounded in [false, true] {
                    for plain in [0., 1.] {
                        let mut data = thread(starts, left).spec();
                        data.major_radius_mm = 1.48;
                        if rounded {
                            data.root_radius_mm = Some(0.05);
                            data.depth_mm =
                                (data.pitch_mm - data.crest_width_mm) * 3_f64.sqrt() / 2. - 0.05;
                        }
                        let t = Thread::new(data).unwrap();
                        let mut s = spec(crowned);
                        s.shank_length_mm = plain;
                        let screw = Screw::new_with_shank_thread_union(t, 1.5, s).unwrap();
                        assert_eq!(*screw.thread(), t);
                        assert_eq!(screw.profile().end()[0], 1.5);
                        let [start, end] = screw.visible_thread_interval_mm();
                        assert_eq!(end, data.length_mm);
                        assert!(!screw.contains([1.49, 0., start + data.pitch_mm], 0.));
                        assert!(screw.contains([1.49, 0., start - 0.0001], 0.));
                        if plain > 0. {
                            assert_eq!(Screw::new_with_shank_radius(t, 1.5, s).unwrap(), screw);
                        }
                        let mesh = tessellate_profiled_screw(&screw, options()).unwrap();
                        let receipt = certify(&mesh).unwrap();
                        assert!(
                            (receipt.signed_volume_mm3 - screw.volume_mm3()).abs()
                                <= profiled_screw_volume_error_bound(&screw, options()).unwrap()
                        );
                        let mut physical = screw.spec();
                        physical.head = None;
                        assert_eq!(
                            Screw::from_source_with_union(t, *screw.profile(), physical)
                                .unwrap()
                                .profile(),
                            screw.profile()
                        );
                        for r in [1.47, f64::NAN, f64::INFINITY] {
                            assert!(Screw::new_with_shank_thread_union(t, r, s).is_err());
                        }
                    }
                }
            }
        }
    }
}
