use super::*;

#[test]
fn reversed_receiving_stock_keeps_its_actual_third_axis_and_numeric_refusals() {
    let policy = GeometryPolicy::default();
    for posed in [false, true] {
        let frame = frame(posed);
        let stock_frame = PlaneFrame {
            origin: position(frame, [0., 0., 6.]),
            u: frame.u,
            v: frame.v.map(|v| -v),
        };
        let receiver = Solid::cylinder(stock_frame, 5., 6., policy).unwrap();
        let bore = receiver
            .rust_cut_internal_thread(spec(false), 1., EntityId::new(), frame, policy)
            .unwrap();
        for (point, occupied) in [
            ([0.1, 0.1, 1.], false),
            ([0.1, 0.1, 5.], true),
            ([4.5, 0., 5.], true),
            ([4.5, 0., 7.], false),
        ] {
            let small = probe(frame, point, 0.01, true);
            let result = bore
                .boolean(&small, BooleanOperation::Intersect, policy)
                .unwrap();
            assert_eq!(result.is_some(), occupied);
            if let Some(result) = result {
                assert_eq!(result.volume_mm3().unwrap(), small.volume_mm3().unwrap());
            }
        }
        let (body, _) = bore.helical_body().unwrap();
        let HelicalForm::Bore(stock) = &body.form else {
            panic!("actual receiving bore");
        };
        assert!(bore_relation(body, stock, [[-f64::MAX, f64::MAX]; 3]).is_none());
        assert!(bore_relation(body, stock, [[f64::NAN, 1.]; 3]).is_none());
    }
}

#[test]
fn receiving_stock_probes_cover_void_wall_floor_and_real_boundary_crossings() {
    let policy = GeometryPolicy::default();
    for left in [false, true] {
        for posed in [false, true] {
            let frame = frame(posed);
            let receiver = Solid::cylinder(frame, 5., 6., policy).unwrap();
            for length in [4.2, 7.] {
                let bore = receiver
                    .rust_cut_internal_thread(
                        HelicalThread {
                            length_mm: length,
                            ..spec(left)
                        },
                        1.,
                        EntityId::new(),
                        frame,
                        policy,
                    )
                    .unwrap();
                for (point, occupied) in [
                    ([0.1, 0.1, 1.], false),
                    ([4.5, 0.1, 1.], true),
                    ([5.5, 0.1, 1.], false),
                    ([0.1, 0.1, -1.], false),
                    ([0.1, 0.1, 7.], false),
                    ([0.1, 0.1, 5.], length < 6.),
                ] {
                    let small = probe(frame, point, 0.01, true);
                    // Independent point witnesses agree, while admission uses
                    // the entire closed box, never only these point samples.
                    for corner in 0..8 {
                        let p = std::array::from_fn(|i| {
                            point[i]
                                + if corner & (1 << i) == 0 {
                                    -0.005
                                } else {
                                    0.005
                                }
                        });
                        assert_eq!(
                            bore.rust_thread_contains(position(frame, p)).unwrap(),
                            occupied
                        );
                    }
                    for (a, b) in [(&bore, &small), (&small, &bore)] {
                        let result = a.boolean(b, BooleanOperation::Intersect, policy).unwrap();
                        assert_eq!(result.is_some(), occupied);
                        if let Some(result) = result {
                            assert_eq!(result.volume_mm3().unwrap(), small.volume_mm3().unwrap());
                        }
                        if !occupied {
                            assert_eq!(
                                a.boolean(b, BooleanOperation::Subtract, policy)
                                    .unwrap()
                                    .unwrap()
                                    .volume_mm3()
                                    .unwrap(),
                                a.volume_mm3().unwrap()
                            );
                        }
                    }
                }
                // Finite cap, outer wall and actual helical wall crossings
                // must not be converted into an unchanged operand.
                let (body, _) = bore.helical_body().unwrap();
                for point in [
                    [5., 0., 1.],
                    [4.5, 0., 0.],
                    [body.geometry.radius(0., 1.), 0., 1.],
                ] {
                    let crossing = probe(frame, point, 0.01, true);
                    assert!(
                        bore.contained_thread_boolean(
                            &crossing,
                            BooleanOperation::Intersect,
                            policy
                        )
                        .is_none()
                    );
                }
            }
        }
    }
}

#[test]
fn actual_same_chart_shafts_clear_finite_cavities_and_nuts_without_waiving_floors_or_phase() {
    let policy = GeometryPolicy::default();
    for left in [false, true] {
        for posed in [false, true] {
            let frame = frame(posed);
            let shaft =
                Solid::rust_helical_thread(spec(left), None, EntityId::new(), frame, policy)
                    .unwrap();
            let receiver = Solid::cylinder(frame, 5., 6., policy).unwrap();
            let cavity_spec = HelicalThread {
                major_radius_mm: 2.2,
                ..spec(left)
            };
            let bore = receiver
                .rust_cut_internal_thread(cavity_spec, 1., EntityId::new(), frame, policy)
                .unwrap();
            let nut = Solid::rust_hex_nut(
                HelicalThread {
                    length_mm: 2.1,
                    ..cavity_spec
                },
                8.,
                0.5,
                EntityId::new(),
                frame,
                policy,
            )
            .unwrap();
            let prefix = receiver
                .rust_cut_internal_thread(
                    HelicalThread {
                        length_mm: 7.,
                        ..cavity_spec
                    },
                    1.,
                    EntityId::new(),
                    frame,
                    policy,
                )
                .unwrap();
            for cavity in [&bore, &nut, &prefix] {
                for (a, b) in [(&shaft, cavity), (cavity, &shaft)] {
                    assert!(
                        a.boolean(b, BooleanOperation::Intersect, policy)
                            .unwrap()
                            .is_none()
                    );
                    let remainder = a
                        .boolean(b, BooleanOperation::Subtract, policy)
                        .unwrap()
                        .unwrap();
                    assert_eq!(remainder.volume_mm3().unwrap(), a.volume_mm3().unwrap());
                    assert!(
                        a.disjoint_thread_mating(b, BooleanOperation::Union)
                            .is_none()
                    );
                }
            }
            // The actual shaft extends into a genuine retained blind floor.
            let long = Solid::rust_helical_thread(
                HelicalThread {
                    length_mm: 5.6,
                    ..spec(left)
                },
                None,
                EntityId::new(),
                frame,
                policy,
            )
            .unwrap();
            let witness = position(frame, [0.1, 0., 4.5]);
            assert!(long.rust_thread_contains(witness).unwrap());
            assert!(bore.rust_thread_contains(witness).unwrap());
            for (a, b) in [(&long, &bore), (&bore, &long)] {
                assert!(
                    a.disjoint_thread_mating(b, BooleanOperation::Intersect)
                        .is_none()
                );
            }
            let wrong = receiver
                .rust_cut_internal_thread(
                    HelicalThread {
                        left_hand: !left,
                        ..cavity_spec
                    },
                    1.,
                    EntityId::new(),
                    frame,
                    policy,
                )
                .unwrap();
            assert!(
                shaft
                    .disjoint_thread_mating(&wrong, BooleanOperation::Intersect)
                    .is_none()
            );
            let shifted = bore
                .transformed(RigidTransform {
                    translation_mm: [0.01, 0., 0.],
                    ..Default::default()
                })
                .unwrap();
            assert!(
                shaft
                    .disjoint_thread_mating(&shifted, BooleanOperation::Intersect)
                    .is_none()
            );
            let interference = receiver
                .rust_cut_internal_thread(
                    HelicalThread {
                        major_radius_mm: 1.8,
                        ..cavity_spec
                    },
                    1.,
                    EntityId::new(),
                    frame,
                    policy,
                )
                .unwrap();
            assert!(
                shaft
                    .disjoint_thread_mating(&interference, BooleanOperation::Intersect)
                    .is_none()
            );
        }
    }
}

fn spec(left: bool) -> HelicalThread {
    HelicalThread {
        major_radius_mm: 2.,
        pitch_mm: 0.7,
        length_mm: 4.2,
        depth_mm: 0.35,
        crest_width_mm: 0.1,
        lower_flank_deg: 30.,
        upper_flank_deg: 35.,
        radial_taper: 0.,
        starts: 1,
        left_hand: left,
    }
}
fn polygon(side: f64, reverse: bool) -> Solid {
    let d = Rational::parse_decimal(&side.to_string()).unwrap();
    let z = Rational::new(0, 1).unwrap();
    let source = multicad_rust_kernel::polygon::PolygonPrism::new(
        vec![vec![[z, z], [d, z], [d, d], [z, d]]],
        d,
    )
    .unwrap();
    let feature = EntityId::new();
    let origins: Vec<_> = (0..4)
        .map(|i| Origin {
            source: feature,
            role: format!("actual.edge.{i}"),
        })
        .collect();
    Solid::from_polygon_prism_sweep(
        source,
        PlaneFrame::default(),
        reverse,
        feature,
        &origins,
        &[false; 4],
        GeometryPolicy::default(),
    )
    .unwrap()
}
fn probe(frame: PlaneFrame, centre: Point, side: f64, prism: bool) -> Solid {
    let policy = GeometryPolicy::default();
    let solid = if prism {
        let d = Rational::parse_decimal(&side.to_string()).unwrap();
        Solid::rust_box(
            EntityId::new().uuid_bits(),
            [d; 3],
            PlaneFrame::default(),
            policy,
        )
        .unwrap()
    } else {
        polygon(side, false)
    };
    let translation_mm = position(frame, centre.map(|v| v - side / 2.));
    let n = frame.normal();
    solid
        .transformed(RigidTransform {
            rotation: std::array::from_fn(|i| [frame.u[i], frame.v[i], n[i]]),
            translation_mm,
        })
        .unwrap()
}
fn frame(posed: bool) -> PlaneFrame {
    if posed {
        PlaneFrame {
            origin: [11., -7., 5.],
            u: [0.8, 0.6, 0.],
            v: [-0.36, 0.48, 0.8],
        }
    } else {
        PlaneFrame::default()
    }
}

#[test]
fn whole_source_relations_return_only_exact_unchanged_operands_or_empty_material() {
    let policy = GeometryPolicy::default();
    for left in [false, true] {
        for posed in [false, true] {
            let frame = frame(posed);
            let thread =
                Solid::rust_helical_thread(spec(left), None, EntityId::new(), frame, policy)
                    .unwrap();
            let before = thread.volume_mm3().unwrap();
            for prism in [false, true] {
                for inside in [false, true] {
                    let small = probe(
                        frame,
                        [if inside { 0.4 } else { 2.5 }, 0.1, 1.],
                        0.01,
                        prism,
                    );
                    let small_volume = small.volume_mm3().unwrap();
                    for reverse in [false, true] {
                        let (a, b) = if reverse {
                            (&small, &thread)
                        } else {
                            (&thread, &small)
                        };
                        for operation in [
                            BooleanOperation::Union,
                            BooleanOperation::Subtract,
                            BooleanOperation::Intersect,
                        ] {
                            let expected = match (inside, operation, reverse) {
                                (true, BooleanOperation::Union, _) => Some(Some(before)),
                                (true, BooleanOperation::Intersect, _) => Some(Some(small_volume)),
                                (true, BooleanOperation::Subtract, true)
                                | (false, BooleanOperation::Intersect, _) => Some(None),
                                (false, BooleanOperation::Subtract, _) => {
                                    Some(Some(if reverse { small_volume } else { before }))
                                }
                                _ => None,
                            };
                            if let Some(expected) = expected {
                                let actual = a.boolean(b, operation, policy).unwrap();
                                assert_eq!(actual.is_some(), expected.is_some());
                                if let Some(expected) = expected {
                                    assert_eq!(actual.unwrap().volume_mm3().unwrap(), expected);
                                }
                            } else {
                                assert!(a.contained_thread_boolean(b, operation, policy).is_none());
                            }
                        }
                    }
                    assert_eq!(thread.volume_mm3().unwrap(), before);
                    assert_eq!(small.volume_mm3().unwrap(), small_volume);
                }
            }
        }
    }
}

#[test]
fn reverse_sweeps_keep_the_actual_negative_axial_extent() {
    let policy = GeometryPolicy::default();
    for posed in [false, true] {
        let frame = frame(posed);
        let thread =
            Solid::rust_helical_thread(spec(false), None, EntityId::new(), frame, policy).unwrap();
        let n = frame.normal();
        for x in [0.4, 2.5] {
            let probe = polygon(0.01, true)
                .transformed(RigidTransform {
                    rotation: std::array::from_fn(|i| [frame.u[i], frame.v[i], n[i]]),
                    translation_mm: position(frame, [x, 0.1, 1.]),
                })
                .unwrap();
            let intersection = thread
                .boolean(&probe, BooleanOperation::Intersect, policy)
                .unwrap();
            assert_eq!(intersection.is_some(), x < 2.);
            if let Some(intersection) = intersection {
                assert_eq!(
                    intersection.volume_mm3().unwrap(),
                    probe.volume_mm3().unwrap()
                );
            }
        }
    }
}
#[test]
fn actual_flank_probes_use_all_box_points_and_boundary_crossings_stay_unresolved() {
    let policy = GeometryPolicy::default();
    for left in [false, true] {
        for posed in [false, true] {
            let frame = frame(posed);
            let thread =
                Solid::rust_helical_thread(spec(left), None, EntityId::new(), frame, policy)
                    .unwrap();
            let (source, _) = thread.helical_body().unwrap();
            for z in [0.23, 0.47, 1.13] {
                let radius = source.geometry.radius(0., z);
                for (offset, occupied) in [(-0.001, true), (0.001, false)] {
                    let small = probe(frame, [radius + offset, 0., z], 0.00002, false);
                    let result = thread
                        .boolean(&small, BooleanOperation::Intersect, policy)
                        .unwrap();
                    assert_eq!(result.is_some(), occupied);
                    if let Some(result) = result {
                        assert_eq!(result.volume_mm3().unwrap(), small.volume_mm3().unwrap());
                    }
                }
                let crossing = probe(frame, [radius, 0., z], 0.01, false);
                assert!(
                    thread
                        .contained_thread_boolean(&crossing, BooleanOperation::Intersect, policy)
                        .is_none()
                );
            }
            let cap = probe(frame, [0.4, 0., 0.], 0.01, false);
            assert!(
                thread
                    .contained_thread_boolean(&cap, BooleanOperation::Intersect, policy)
                    .is_none()
            );
        }
    }
}

#[test]
fn rounded_root_probes_retain_the_actual_circle_for_both_hands_and_poses() {
    let policy = GeometryPolicy::default();
    for left in [false, true] {
        for posed in [false, true] {
            let frame = frame(posed);
            let mut parameters = spec(left);
            parameters.upper_flank_deg = 30.;
            parameters.depth_mm = (0.7 - 0.1) / (2. * 30_f64.to_radians().tan())
                - 0.06 * (1. / 30_f64.to_radians().sin() - 1.);
            let thread =
                Solid::rust_helical_thread(parameters, Some(0.06), EntityId::new(), frame, policy)
                    .unwrap();
            let (source, _) = thread.helical_body().unwrap();
            for z in [0.32, 0.35, 0.38] {
                let radius = source.geometry.radius(0., z);
                for (offset, occupied) in [(-0.001, true), (0.001, false)] {
                    let small = probe(frame, [radius + offset, 0., z], 0.00002, false);
                    let result = thread
                        .boolean(&small, BooleanOperation::Intersect, policy)
                        .unwrap();
                    assert_eq!(result.is_some(), occupied);
                    if let Some(result) = result {
                        assert_eq!(result.volume_mm3().unwrap(), small.volume_mm3().unwrap());
                    }
                }
            }
        }
    }
}
#[test]
fn unrelated_families_curved_proposals_and_fuzzy_or_mismatched_policies_keep_their_dispatch() {
    let policy = GeometryPolicy::default();
    let thread =
        Solid::rust_helical_thread(spec(false), None, EntityId::new(), frame(false), policy)
            .unwrap();
    let small = probe(frame(false), [0.4, 0., 1.], 0.01, false);
    let sphere = Solid::rust_sphere_named(0.01, frame(false), EntityId::new(), policy).unwrap();
    let q = |n| Rational::new(n, 1).unwrap();
    let curved = Solid::from_rust_prism(
        multicad_rust_kernel::prism::Prism::extrude(42, q(10), q(2), q(2)).unwrap(),
        policy,
    )
    .unwrap();
    let nut =
        Solid::rust_hex_nut(spec(false), 8., 0.5, EntityId::new(), frame(false), policy).unwrap();
    for (a, b) in [(&small, &sphere), (&thread, &sphere), (&thread, &curved)] {
        assert!(
            a.contained_thread_boolean(b, BooleanOperation::Intersect, policy)
                .is_none()
        );
    }
    assert!(
        nut.contained_thread_boolean(&small, BooleanOperation::Intersect, policy)
            .expect("a probe wholly inside the retained threaded cavity is proved empty")
            .unwrap()
            .is_none()
    );
    for changed in [
        GeometryPolicy {
            boolean_fuzzy_mm: 1e-6,
            ..policy
        },
        GeometryPolicy {
            linear_tolerance_mm: 1e-6,
            ..policy
        },
    ] {
        assert!(
            thread
                .contained_thread_boolean(&small, BooleanOperation::Intersect, changed)
                .is_none()
        );
    }
}
