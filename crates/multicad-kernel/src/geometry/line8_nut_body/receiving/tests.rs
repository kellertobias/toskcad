use super::*;
use multicad_core::{
    materials::{ProfileShape, Stock, StockKind},
    model::{
        self, PathSegment, ProfileGrooveSide, ProfileSection, ResolvedMember, SectionEdge,
        SectionLoop,
    },
};
fn resolved(start: Point, end: Point) -> ResolvedMember {
    let stock = Stock {
        material: EntityId::new(),
        specification: StockKind::Profile {
            designation: "item 0.0.026.03".into(),
            length_mm: 6000.,
            width_mm: 40.,
            height_mm: 40.,
            wall_mm: None,
            shape: ProfileShape::ItemLine8_40x40,
            section: None,
            corner_radius_mm: None,
            web_mm: None,
        },
    };
    let delta = std::array::from_fn::<_, 3, _>(|i| end[i] - start[i]);
    let length = delta[0].portable_hypot(delta[1]).portable_hypot(delta[2]);
    ResolvedMember {
        stock: EntityId::new(),
        section: ProfileSection::from_stock(&stock).unwrap(),
        stock_spec: stock,
        path: PathSegment {
            start,
            end,
            reference: None,
        },
        path_length_mm: length,
        frame: model::member_frame(delta.map(|v| v / length), None, 0.),
        anchor_mm: [20.; 2],
        rotation_rad: 0.,
        extend_mm: [0.; 2],
    }
}
fn edges(wire: &SectionLoop) -> Vec<crate::CurvedEdge2d> {
    match wire {
        SectionLoop::Curved(edges) => edges
            .iter()
            .map(|edge| match *edge {
                SectionEdge::Line { a, b } => crate::CurvedEdge2d::Line { a, b },
                SectionEdge::Arc {
                    center,
                    a,
                    b,
                    clockwise,
                } => crate::CurvedEdge2d::Arc {
                    center,
                    a,
                    b,
                    clockwise,
                },
            })
            .collect(),
        SectionLoop::Circle { center, radius_mm } => vec![crate::CurvedEdge2d::Circle {
            center: *center,
            radius_mm: *radius_mm,
        }],
        SectionLoop::Polygon(points) => (0..points.len())
            .map(|i| crate::CurvedEdge2d::Line {
                a: points[i],
                b: points[(i + 1) % points.len()],
            })
            .collect(),
    }
}
fn actual_member(member: &ResolvedMember) -> Solid {
    actual_member_mouth(member, 0.)
}
fn actual_member_mouth(member: &ResolvedMember, inward_lip_mm: f64) -> Solid {
    let mut profile = crate::CurvedProfile {
        frame: PlaneFrame {
            origin: member.section_origin(0.),
            u: member.frame.x,
            v: member.frame.y,
        },
        outer: edges(&member.section.outer(0.)),
        holes: member.section.holes().iter().map(edges).collect(),
    };
    if inward_lip_mm != 0. {
        for edge in &mut profile.outer {
            if let crate::CurvedEdge2d::Line { a, b } = edge {
                for p in [a, b] {
                    if [0., 4.5, 35.5, 40.].contains(&p[1]) {
                        if p[0] == 16. {
                            p[0] += inward_lip_mm;
                        }
                        if p[0] == 24. {
                            p[0] -= inward_lip_mm;
                        }
                    }
                }
            }
        }
    }
    let owner = EntityId::new();
    let origins = (0..profile.outer.len() + profile.holes.iter().map(Vec::len).sum::<usize>())
        .map(|i| Origin {
            source: owner,
            role: format!("actual-member.{i}"),
        })
        .collect::<Vec<_>>();
    Solid::rust_extrude_normal_profile_named(
        &profile,
        member.length_mm(),
        owner,
        &origins,
        GeometryPolicy::default(),
    )
    .unwrap()
}
fn actual_nut(member: &ResolvedMember) -> Result<Solid, GeometryError> {
    let surface = member.groove_frame(ProfileGrooveSide::MaxY, 40.).unwrap();
    let inward = PlaneFrame {
        origin: surface.origin_mm,
        u: surface.u,
        v: surface.v.map(|v| -v),
    };
    let normal = inward.normal();
    let seat = std::array::from_fn(|i| inward.origin[i] + normal[i] * 4.5);
    let section = [
        [-6.9, 0.],
        [6.9, 0.],
        [6.5, 3.],
        [5.4, 5.],
        [3.3, 7.3],
        [-3.3, 7.3],
        [-5.4, 5.],
        [-6.5, 3.],
    ];
    let outer = (0..8)
        .map(|i| crate::CurvedEdge2d::Line {
            a: section[i],
            b: section[(i + 1) % 8],
        })
        .collect();
    let blank = Solid::extrude_curved_profile(
        &crate::CurvedProfile {
            frame: PlaneFrame {
                origin: std::array::from_fn(|i| seat[i] - inward.u[i] * 7.5),
                u: inward.v,
                v: normal,
            },
            outer,
            holes: vec![],
        },
        23.,
        GeometryPolicy::default(),
    )
    .unwrap();
    let pitch = 1.25;
    let height = 4.4;
    let external_start = -5. + height + 1.25;
    let station: f64 = 4.5 - external_start;
    let angle = std::f64::consts::TAU * station.rem_euclid(pitch) / pitch;
    let (s, c) = angle.portable_sin_cos();
    let raw = s.portable_atan2(c);
    let roll = (raw + std::f64::consts::TAU / 2.).rem_euclid(std::f64::consts::TAU)
        - std::f64::consts::TAU / 2.;
    let (s, c) = roll.portable_sin_cos();
    let u = std::array::from_fn(|i| inward.u[i] * c + inward.v[i] * s);
    let v = std::array::from_fn(|i| inward.u[i] * -s + inward.v[i] * c);
    blank.rust_cut_line8_nut(
        Line8NutCut {
            profile: HelicalThread {
                major_radius_mm: 4. + 0.08,
                pitch_mm: pitch,
                length_mm: 4.1,
                depth_mm: 5. * 3_f64.sqrt() / 16. * pitch,
                crest_width_mm: pitch / 8. + 0.04,
                lower_flank_deg: 30.,
                upper_flank_deg: 30.,
                radial_taper: 0.,
                starts: 1,
                left_hand: false,
            },
            frame: PlaneFrame { origin: seat, u, v },
            bore_from_start_mm: 7.5,
            relief: Relief {
                start_mm: 4.1 - 0.001,
                radius_mm: 4.10,
                tool_length_mm: 7.3 - 4.1 + 0.02,
            },
            minimum_wall_mm: 0.5,
        },
        EntityId::new(),
        GeometryPolicy::default(),
        None,
    )
}
#[test]
fn original_model_member_sections_and_both_occurrences_earn_separate_nut_clearance() {
    for rotated in [false, true] {
        let rotate = |[x, y, z]: Point| {
            if rotated {
                [0.8 * x - 0.6 * y, 0.6 * x + 0.8 * y, z]
            } else {
                [x, y, z]
            }
        };
        let members = [
            ([0., 0., 0.], [100., 0., 0.]),
            ([0., 40., 0.], [100., 40., 0.]),
        ]
        .map(|(a, b)| resolved(rotate(a), rotate(b)));
        let solids = members.each_ref().map(actual_member);
        let direct = actual_nut(&members[0]);
        let nut = if rotated {
            // The generic vector extrusion retains its genuine tiny shear.
            // Never alias that unsupported source to a straight prism.
            assert_eq!(direct.err().unwrap().kind, GeometryErrorKind::Unsupported);
            actual_nut(&resolved([0., 0., 0.], [100., 0., 0.]))
                .unwrap()
                .transformed(RigidTransform {
                    rotation: [[0.8, -0.6, 0.], [0.6, 0.8, 0.], [0., 0., 1.]],
                    ..RigidTransform::default()
                })
                .unwrap()
        } else {
            direct.unwrap()
        };
        let started = std::time::Instant::now();
        let a = nut
            .rust_line8_receiver_clearance(&solids[0], RigidTransform::default(), 0.001, None)
            .unwrap();
        let b = nut
            .rust_line8_receiver_clearance(&solids[1], RigidTransform::default(), 0.001, None)
            .unwrap();
        eprintln!(
            "Line8 receiver rotated={rotated}: first+second={:?}, charged_work={}/{}, peak_bytes={}/{}, overlap_upper={}/{} mm3",
            started.elapsed(),
            a.work_used,
            b.work_used,
            a.accounted_peak_bytes,
            b.accounted_peak_bytes,
            a.overlap_volume_upper_mm3,
            b.overlap_volume_upper_mm3
        );
        assert!(a.used_exterior_cell);
        assert!(!b.used_exterior_cell);
        assert!(a.overlap_volume_upper_mm3 < 1e-8);
        assert_eq!(b.overlap_volume_upper_mm3, 0.);
        assert!(std::ptr::eq(
            a.nut().source,
            nut.rust_line8_nut_source().unwrap()
        ));
        assert_eq!(a.receiver().source.rings()[0].len(), 56);
        assert_eq!(a.receiver().source.rings().len(), 2);
        assert_eq!(
            solids[0]
                .rust_line8_receiver_clearance(&solids[1], RigidTransform::default(), 0.001, None)
                .err()
                .unwrap()
                .kind,
            GeometryErrorKind::Unsupported
        );
        assert_eq!(
            nut.rust_line8_receiver_clearance(&nut, RigidTransform::default(), 0.001, None)
                .err()
                .unwrap()
                .kind,
            GeometryErrorKind::Unsupported
        );
        assert_eq!(
            nut.rust_line8_receiver_clearance(
                &solids[0],
                RigidTransform::default(),
                f64::NAN,
                None
            )
            .err()
            .unwrap()
            .kind,
            GeometryErrorKind::Input
        );
        assert_eq!(
            nut.rust_line8_receiver_clearance(
                &solids[0],
                RigidTransform {
                    translation_mm: [f64::NAN, 0., 0.],
                    ..RigidTransform::default()
                },
                0.001,
                None
            )
            .err()
            .unwrap()
            .kind,
            GeometryErrorKind::Input
        );
        let other = Solid::box_solid(
            [1.; 3],
            GeometryPolicy {
                linear_tolerance_mm: 2e-7,
                ..GeometryPolicy::default()
            },
        )
        .unwrap();
        assert_eq!(
            nut.rust_line8_receiver_clearance(&other, RigidTransform::default(), 0.001, None)
                .err()
                .unwrap()
                .kind,
            GeometryErrorKind::Input
        );
        if !rotated {
            let points = (0..300)
                .map(|i| {
                    let (s, c) = (std::f64::consts::TAU * i as f64 / 300.).portable_sin_cos();
                    [20. * c, 20. * s]
                })
                .collect::<Vec<_>>();
            let outer = (0..points.len())
                .map(|i| crate::CurvedEdge2d::Line {
                    a: points[i],
                    b: points[(i + 1) % points.len()],
                })
                .collect();
            let owner = EntityId::new();
            let origins = (0..300)
                .map(|i| Origin {
                    source: owner,
                    role: format!("actual-complex.{i}"),
                })
                .collect::<Vec<_>>();
            let complex = Solid::rust_extrude_normal_profile_named(
                &crate::CurvedProfile {
                    frame: PlaneFrame::default(),
                    outer,
                    holes: vec![],
                },
                100.,
                owner,
                &origins,
                GeometryPolicy::default(),
            )
            .unwrap();
            assert_eq!(
                nut.rust_line8_receiver_clearance(&complex, RigidTransform::default(), 0.001, None)
                    .err()
                    .unwrap()
                    .kind,
                GeometryErrorKind::Limit
            );
        }
        let pose = RigidTransform {
            rotation: [[0.8, -0.6, 0.], [0.6, 0.8, 0.], [0., 0., 1.]],
            translation_mm: [37.1, -11.3, 103.7],
        };
        let moved = nut.transformed(pose).unwrap();
        for receiver in &solids {
            assert!(
                moved
                    .rust_line8_receiver_clearance(receiver, pose, 0.001, None)
                    .unwrap()
                    .overlap_volume_upper_mm3
                    < 1e-8
            );
        }
        let surface = members[0]
            .groove_frame(ProfileGrooveSide::MaxY, 40.)
            .unwrap();
        let penetrating = nut
            .transformed(RigidTransform {
                translation_mm: surface.normal().map(|v| v * 0.05),
                ..RigidTransform::default()
            })
            .unwrap();
        assert_eq!(
            penetrating
                .rust_line8_receiver_clearance(&solids[0], RigidTransform::default(), 0.001, None)
                .err()
                .unwrap()
                .kind,
            GeometryErrorKind::Unsupported
        );
        assert_eq!(
            nut.rust_line8_receiver_clearance(
                &solids[0],
                RigidTransform::default(),
                0.001,
                Some(&|| true)
            )
            .err()
            .unwrap()
            .kind,
            GeometryErrorKind::Cancelled
        );
    }
}
#[test]
fn actual_profiled_screw_receiver_guard_preserves_errors_and_inward_collision() {
    use multicad_rust_kernel::{
        head_profile::Head,
        helical_thread::{Spec, Thread},
        profiled_screw::{Screw, Spec as ScrewSpec},
    };
    let member = resolved([0., 0., 0.], [100., 0., 0.]);
    let receiver = actual_member(&member);
    let surface = member.groove_frame(ProfileGrooveSide::MaxY, 40.).unwrap();
    let mut frame = PlaneFrame {
        origin: surface.origin_mm,
        u: surface.u,
        v: surface.v.map(|v| -v),
    };
    let normal = frame.normal();
    frame.origin = std::array::from_fn(|i| frame.origin[i] + normal[i] * 0.65);
    let source = Screw::new_with_shank_radius(
        Thread::new(Spec {
            major_radius_mm: 3.92,
            pitch_mm: 1.25,
            length_mm: 8.35,
            depth_mm: 0.65,
            crest_width_mm: 0.15,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            starts: 1,
            left_hand: false,
            root_radius_mm: None,
        })
        .unwrap(),
        4.,
        ScrewSpec {
            head: Head::Countersunk {
                radius_mm: 8.,
                height_mm: 4.4,
                included_angle_deg: 90.,
                edge_chamfer_mm: 0.05,
                neck_radius_mm: 0.4,
            },
            shank_length_mm: 1.25,
            socket_across_flats_mm: 5.,
            socket_depth_mm: 3.2,
            minimum_wall_mm: 0.5,
        },
    )
    .unwrap();
    let screw = Solid::rust_profiled_screw_from_source(
        source,
        EntityId::new(),
        frame,
        GeometryPolicy::default(),
    )
    .unwrap();
    let proof = screw
        .rust_profiled_screw_receiver_clearance(&receiver, RigidTransform::default(), 0.001, None)
        .unwrap();
    assert!(proof.overlap_volume_upper_mm3 < 1e-8);
    assert!(std::ptr::eq(
        proof.screw().source,
        screw.rust_profiled_screw_geometry().unwrap().0
    ));
    assert_eq!(
        screw
            .rust_profiled_screw_receiver_clearance(
                &receiver,
                RigidTransform::default(),
                f64::NAN,
                None
            )
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Input
    );
    let narrowed = actual_member_mouth(&member, 0.1);
    assert!(narrowed.volume_mm3().unwrap() > receiver.volume_mm3().unwrap());
    assert_eq!(
        screw
            .rust_profiled_screw_receiver_clearance(
                &narrowed,
                RigidTransform::default(),
                0.001,
                None
            )
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Unsupported
    );
    let actual = screw.rust_profiled_screw_geometry().unwrap().0;
    let mut thread = actual.thread().spec();
    thread.length_mm = 14.;
    let long_source = Screw::from_source(
        Thread::new(thread).unwrap(),
        *actual.profile(),
        multicad_rust_kernel::profiled_screw::SourceSpec {
            head: None,
            ..actual.spec()
        },
        actual.recess().clone(),
    )
    .unwrap();
    let long = Solid::rust_profiled_screw_from_source(
        long_source,
        EntityId::new(),
        frame,
        GeometryPolicy::default(),
    )
    .unwrap();
    assert_eq!(
        long.rust_profiled_screw_receiver_clearance(
            &receiver,
            RigidTransform::default(),
            0.001,
            None
        )
        .err()
        .unwrap()
        .kind,
        GeometryErrorKind::Unsupported
    );
    let polls = std::cell::Cell::new(0);
    assert_eq!(
        screw
            .rust_profiled_screw_receiver_clearance(
                &receiver,
                RigidTransform::default(),
                0.001,
                Some(&|| {
                    let n = polls.get() + 1;
                    polls.set(n);
                    n > 20
                })
            )
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Cancelled
    );
    assert!(polls.get() > 20);
    let inward = screw
        .transformed(RigidTransform {
            translation_mm: normal.map(|v| v * 0.8),
            ..Default::default()
        })
        .unwrap();
    assert_eq!(
        inward
            .rust_profiled_screw_receiver_clearance(
                &receiver,
                RigidTransform::default(),
                0.001,
                None
            )
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(
        receiver
            .rust_profiled_screw_receiver_clearance(
                &receiver,
                RigidTransform::default(),
                0.001,
                None
            )
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(
        screw
            .rust_profiled_screw_receiver_clearance(&screw, RigidTransform::default(), 0.001, None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Unsupported
    );
}
