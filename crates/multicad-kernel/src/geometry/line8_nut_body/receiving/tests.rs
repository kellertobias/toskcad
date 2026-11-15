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
        },
    };
    let delta = std::array::from_fn::<_, 3, _>(|i| end[i] - start[i]);
    let length = delta[0].hypot(delta[1]).hypot(delta[2]);
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
                SectionEdge::Arc { center, a, b } => crate::CurvedEdge2d::Arc {
                    center,
                    a,
                    b,
                    clockwise: false,
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
    let profile = crate::CurvedProfile {
        frame: PlaneFrame {
            origin: member.section_origin(0.),
            u: member.frame.x,
            v: member.frame.y,
        },
        outer: edges(&member.section.outer(0.)),
        holes: member.section.holes().iter().map(edges).collect(),
    };
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
    let (s, c) = angle.sin_cos();
    let raw = s.atan2(c);
    let roll = (raw + std::f64::consts::TAU / 2.).rem_euclid(std::f64::consts::TAU)
        - std::f64::consts::TAU / 2.;
    let (s, c) = roll.sin_cos();
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
