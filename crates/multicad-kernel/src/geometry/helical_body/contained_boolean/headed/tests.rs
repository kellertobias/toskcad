use super::*;
use multicad_rust_kernel::headed_screw::ScrewSpec;

fn thread(left: bool) -> HelicalThread {
    HelicalThread {
        major_radius_mm: 2.,
        pitch_mm: 0.7,
        length_mm: 4.2,
        depth_mm: 0.4,
        crest_width_mm: 0.1,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        radial_taper: 0.,
        starts: 1,
        left_hand: left,
    }
}
fn probe(frame: PlaneFrame, center: Point) -> Solid {
    let point = position(frame, center);
    Solid::box_solid([0.01; 3], Default::default())
        .unwrap()
        .transformed(RigidTransform {
            translation_mm: point.map(|v| v - 0.005),
            ..Default::default()
        })
        .unwrap()
}
#[test]
fn actual_heads_sockets_shanks_and_nuts_retain_only_proved_whole_operands() {
    let policy = GeometryPolicy::default();
    for left in [false, true] {
        for posed in [false, true] {
            let frame = if posed {
                PlaneFrame {
                    origin: [11., -7., 5.],
                    u: [0.8, 0.6, 0.],
                    v: [-0.36, 0.48, 0.8],
                }
            } else {
                PlaneFrame::default()
            };
            for head in [
                ScrewHead::SocketCap {
                    radius_mm: 4.5,
                    height_mm: 2.5,
                    socket_across_flats_mm: 2.,
                    socket_depth_mm: 1.,
                },
                ScrewHead::ExternalHex {
                    across_flats_mm: 8.,
                    height_mm: 2.5,
                },
            ] {
                let solid = Solid::rust_machine_screw_at_thread_frame(
                    thread(left),
                    None,
                    ScrewSpec {
                        head,
                        shank_length_mm: 1.3,
                        shank_radius_mm: 2.1,
                    },
                    EntityId::new(),
                    frame,
                    policy,
                )
                .unwrap();
                let socket = head.socket().is_some();
                for (point, occupied) in [
                    ([0.1, 0.1, -3.3], !socket),
                    ([3.5, 0.1, -3.3], true),
                    ([0.1, 0.1, -2.2], true),
                    ([0.1, 0.1, -0.6], true),
                    ([0.1, 0.1, 2.], true),
                    ([5., 0.1, -3.3], false),
                    ([0.1, 0.1, -4.], false),
                    ([0.1, 0.1, 5.], false),
                ] {
                    assert_eq!(
                        solid.rust_thread_contains(position(frame, point)).unwrap(),
                        occupied
                    );
                    let small = probe(frame, point);
                    for reversed in [false, true] {
                        let result = if reversed {
                            small.boolean(&solid, BooleanOperation::Intersect, policy)
                        } else {
                            solid.boolean(&small, BooleanOperation::Intersect, policy)
                        }
                        .unwrap();
                        assert_eq!(result.is_some(), occupied, "{left}/{posed}/{point:?}");
                        if let Some(result) = result {
                            assert_eq!(result.volume_mm3().unwrap(), small.volume_mm3().unwrap());
                        }
                    }
                }
                let boundary = probe(
                    frame,
                    if socket {
                        [1., 0., -3.3]
                    } else {
                        [4., 0., -3.3]
                    },
                );
                assert_eq!(
                    solid
                        .boolean(&boundary, BooleanOperation::Intersect, policy)
                        .err()
                        .unwrap()
                        .kind,
                    GeometryErrorKind::Unsupported
                );
                let (body, _) = solid.helical_body().unwrap();
                let HelicalForm::Screw(source) = &body.form else {
                    panic!()
                };
                assert_eq!(
                    screw_relation(source, [[f64::NAN, 1.]; 3]),
                    UnresolvedBoundary
                );
                assert_eq!(screw_relation(source, [[1., 0.]; 3]), UnresolvedBoundary);
            }
            let nut =
                Solid::rust_hex_nut(thread(left), 8., 1., EntityId::new(), frame, policy).unwrap();
            for (point, occupied) in [
                ([0.1, 0.1, 2.], false),
                ([3.9, 0., 2.], true),
                ([0., 4.5, 2.], true),
                ([4.1, 0., 2.], false),
                ([0., 4.7, 2.], false),
                ([3., 0., -0.1], false),
                ([3., 0., 4.3], false),
            ] {
                assert_eq!(
                    nut.rust_thread_contains(position(frame, point)).unwrap(),
                    occupied
                );
                let small = probe(frame, point);
                let result = nut
                    .boolean(&small, BooleanOperation::Intersect, policy)
                    .unwrap();
                assert_eq!(result.is_some(), occupied, "nut:{left}/{posed}/{point:?}");
                if let Some(result) = result {
                    assert_eq!(result.volume_mm3().unwrap(), small.volume_mm3().unwrap());
                }
            }
            let boundary = probe(frame, [4., 0., 2.]);
            assert_eq!(
                nut.boolean(&boundary, BooleanOperation::Intersect, policy)
                    .err()
                    .unwrap()
                    .kind,
                GeometryErrorKind::Unsupported
            );
            let (body, _) = nut.helical_body().unwrap();
            let HelicalForm::Nut(source) = &body.form else {
                panic!()
            };
            assert_eq!(
                nut_relation(source, [[f64::NAN, 1.]; 3]),
                UnresolvedBoundary
            );
        }
    }
}

#[test]
fn numeric_enclosures_and_empty_shanks_cannot_publish_uncertain_material() {
    assert_eq!(slab(B::new(f64::NAN, 1.), 0., 1.), UnresolvedBoundary);
    assert_eq!(
        circle([B::new(-f64::MAX, f64::MAX); 2], 2.),
        UnresolvedBoundary
    );
    assert_eq!(
        circle([B::new(1.9, 2.1), B::point(0.)], 2.),
        UnresolvedBoundary
    );
    assert_eq!(
        hexagon([B::new(-f64::MAX, f64::MAX); 2], 8.),
        UnresolvedBoundary
    );
    let source = Screw::new(
        Thread::new(Spec {
            major_radius_mm: 2.,
            pitch_mm: 0.7,
            length_mm: 4.2,
            depth_mm: 0.4,
            crest_width_mm: 0.1,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            root_radius_mm: None,
            starts: 1,
            left_hand: false,
        })
        .unwrap(),
        ScrewSpec {
            head: ScrewHead::ExternalHex {
                across_flats_mm: 8.,
                height_mm: 2.5,
            },
            shank_length_mm: 0.,
            shank_radius_mm: 2.,
        },
    )
    .unwrap();
    assert_eq!(
        screw_relation(&source, [[0.1, 0.2], [0.1, 0.2], [5., 6.]]),
        Outside
    );
    assert_eq!(
        screw_relation(&source, [[-f64::MAX, f64::MAX]; 3]),
        UnresolvedBoundary
    );
    let inside = [[0.1, 0.2], [0.1, 0.2], [-1., -0.9]];
    assert_eq!(screw_relation(&source, inside), Inside);
    let cavity = Nut::new(*source.thread(), 8.).unwrap();
    assert_eq!(nut_relation(&cavity, [[1., 0.]; 3]), UnresolvedBoundary);
    assert_eq!(
        nut_relation(&cavity, [[-f64::MAX, f64::MAX]; 3]),
        UnresolvedBoundary
    );
}
