use super::*;

fn square(left: bool, wraps: usize) -> Vec<Degree11ThreadFace> {
    let corners = [[1., -1.], [1., 1.], [-1., 1.], [-1., -1.]];
    let count = 4 * wraps;
    (0..count)
        .map(|face| {
            let index = if left { (4 - face % 4) % 4 } else { face % 4 };
            let next = if left {
                (index + 3) % 4
            } else {
                (index + 1) % 4
            };
            let poles = (0..12)
                .flat_map(|u| {
                    (0..12).map(move |v| {
                        let t = u as f64 / 11.;
                        [
                            corners[index][0] * (1. - t) + corners[next][0] * t,
                            corners[index][1] * (1. - t) + corners[next][1] * t,
                            2. * v as f64 / 11.,
                        ]
                    })
                })
                .collect();
            Degree11ThreadFace {
                source: RetainedTensor11::new([1, 1], poles).unwrap(),
                original_controls: (0..144).collect(),
                forward: !left,
                transfer: ThreadCarrierTransfer {
                    face_entity: face as u32 + 1,
                    surface_entity: face as u32 + 101,
                    original_rail_entities: [face as u32, ((face + 1) % count) as u32],
                    original_cells: [1, 1],
                    source_cells: [1, 1],
                    ..Default::default()
                },
            }
        })
        .collect()
}
fn run(faces: &[Degree11ThreadFace], o: Degree11OccupiedOptions) -> Result<Proof> {
    let mut work = Work {
        remaining: o.maximum_work,
        interrupt: None,
        cancelled: Cell::new(false),
    };
    certify(faces, [0., 2.], o, &mut work)
}
fn change(face: &mut Degree11ThreadFace, f: impl FnOnce(&mut Vec<Point>)) {
    let mut poles = face.source.poles().to_vec();
    f(&mut poles);
    face.source = RetainedTensor11::new(face.source.spans(), poles).unwrap();
}
fn polynomial() -> [[Point; 12]; 12] {
    // Degree-elevated x=2+u²+uv,y=-1+u+v²,z=2v.
    std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            let u = i as f64 / 11.;
            let v = j as f64 / 11.;
            let u2 = if i < 2 {
                0.
            } else {
                (i * (i - 1)) as f64 / 110.
            };
            let v2 = if j < 2 {
                0.
            } else {
                (j * (j - 1)) as f64 / 110.
            };
            [2. + u2 + u * v, -1. + u + v2, 2. * v]
        })
    })
}
#[test]
fn antisymmetric_coefficients_match_independent_full_exact_product_oracle() {
    use num_rational::BigRational as Q;
    let mut work = Work {
        remaining: MAX_WORK,
        interrupt: None,
        cancelled: Cell::new(false),
    };
    let weights = Weights::new(&mut work).unwrap();
    let points = polynomial();
    let actual = angular_coefficients(&points, &weights, &mut work).unwrap();
    let q = |x| Q::from_float(x).unwrap();
    let ratio = |a, b| Q::new((a as u64).into(), (b as u64).into());
    let p: Vec<_> = points
        .iter()
        .map(|row| row.iter().map(|p| [q(p[0]), q(p[1])]).collect::<Vec<_>>())
        .collect();
    let mut expected = vec![Q::from_integer(0.into()); 22 * 23];
    // Independent generic P×P_u product, degree21×22, using exact rational
    // representations of the ACTUAL binary64 coefficients. No antisymmetric
    // weights or real-polynomial rounding assumptions enter this oracle.
    for i in 0..12 {
        for k in 0..11 {
            for j in 0..12 {
                for b in 0..12 {
                    let dx = (&p[k + 1][b][0] - &p[k][b][0]) * Q::from_integer(11.into());
                    let dy = (&p[k + 1][b][1] - &p[k][b][1]) * Q::from_integer(11.into());
                    let c = &p[i][j][0] * dy - &p[i][j][1] * dx;
                    let weight = ratio(choose(11, i) * choose(10, k), choose(21, i + k))
                        * ratio(choose(11, j) * choose(11, b), choose(22, j + b));
                    expected[(i + k) * 23 + j + b] += c * weight;
                }
            }
        }
    }
    for i in 0..22 {
        for j in 0..23 {
            let mut raised = I::point(0.);
            if i > 0 {
                raised = raised + actual[(i - 1) * 23 + j] * (I::point(i as f64) / I::point(21.));
            }
            if i < 21 {
                raised = raised + actual[i * 23 + j] * (I::point((21 - i) as f64) / I::point(21.));
            }
            assert!(
                q(raised.lo) <= expected[i * 23 + j] && q(raised.hi) >= expected[i * 23 + j],
                "coefficient{i},{j}"
            );
        }
    }
}
#[test]
fn nontrivial_polynomial_leaf_jacobian_keeps_both_original_chart_factors() {
    let points = polynomial();
    let tensor = RetainedTensor11::new([1, 1], points.into_iter().flatten().collect()).unwrap();
    let u = tensor.cell([0, 0]).unwrap().subdivide(0).unwrap()[0].clone();
    let cell = u.subdivide(1).unwrap()[1].clone();
    assert_eq!(cell.domain(), [[0., 0.5], [0.5, 1.]]);
    let d = cell.derivative_bounds_at([[0.5, 0.5]; 2]).unwrap();
    let p = d.position.map(|x| I::new(x[0], x[1]));
    let du = d.du.map(|x| I::new(x[0], x[1]));
    let dv = d.dv.map(|x| I::new(x[0], x[1]));
    let j =
        cross([p[0], p[1]], [du[0], du[1]]) * dv[2] - cross([p[0], p[1]], [dv[0], dv[1]]) * du[2];
    // Analytic J at original(.25,.75), times BOTH dyadic chart factors.
    let expected = 4.96875 * 0.5 * 0.5;
    assert!(j.lo <= expected + 1e-12 && j.hi >= expected - 1e-12);
    assert!(j.lo > expected - 1e-10 && j.hi < expected + 1e-10);
}
#[test]
fn directed_normal_norm_stays_positive_finite_through_square_underflow() {
    for value in [0., f64::from_bits(1), 1e-300, 1e-100, 1., 5e17] {
        let length = (I::point(value).square() + I::point(0.).square()).sqrt();
        assert!(length.valid() && length.hi > 0. && length.hi.is_finite());
        assert!(length.lo <= value && length.hi >= value);
    }
    // A downward Cauchy quotient may still underflow. That remains an
    // unearned lower bound, never an authorization to claim positive radius.
    assert!(!positive(I::point(f64::from_bits(1)) / I::point(1e18)));
}
#[test]
fn independent_exact_prism_normals_jacobian_caps_and_occupied_core_both_hands() {
    // Each exact degree-elevated bilinear wall bounds the square prism
    // [-1,1]^2 ×[0,2]. Independent analytic du,dv give outward normal area4,
    // cylindrical determinant4 and an occupied unit-radius core. No mesh or
    // inferred thread law supplies these expectations.
    for left in [false, true] {
        let faces = square(left, 1);
        for face in &faces {
            let d = face
                .source
                .cell([0, 0])
                .unwrap()
                .first_derivative_bounds()
                .unwrap();
            let a = face.source.poles()[0];
            let b = face.source.poles()[132];
            let du = [b[0] - a[0], b[1] - a[1], 0.];
            for (axis, &value) in du.iter().enumerate() {
                assert!(d.du[axis][0] <= value && d.du[axis][1] >= value);
                let dv = if axis == 2 { 2. } else { 0. };
                assert!(d.dv[axis][0] <= dv && d.dv[axis][1] >= dv);
            }
            let normal = [2. * du[1], -2. * du[0]];
            let sign = if left { -1. } else { 1. };
            assert_eq!(sign * (a[0] * normal[0] + a[1] * normal[1]), 4.);
        }
        let proof = run(&faces, Default::default()).unwrap();
        assert!(proof.radius > 0.999999999999 && proof.radius <= 1.);
        assert!(proof.jacobian > 3.99999999999 && proof.jacobian <= 4.);
        assert!(proof.axial > 1.99999999999 && proof.axial <= 2.);
        assert_eq!(proof.enclosure, [[-1., 1.], [-1., 1.], [0., 2.]]);
        for angle in proof.winding {
            assert!(angle[0] <= 2. * std::f64::consts::PI);
            assert!(angle[1] >= 2. * std::f64::consts::PI);
        }
        assert_eq!(proof.cells, 4);
    }
}
#[test]
fn actual_graph_caps_orientation_wrap_and_core_failures_never_publish_occupancy() {
    let o = Degree11OccupiedOptions::default();
    let mut cases = Vec::new();
    let mut face = square(false, 1);
    face[0].transfer.original_rail_entities[1] = 99;
    cases.push(face);
    let mut face = square(false, 1);
    change(&mut face[0], |p| p[143][0] += 0.001);
    cases.push(face);
    let mut face = square(false, 1);
    face[0].forward = false;
    cases.push(face);
    let mut face = square(false, 1);
    for f in &mut face {
        f.forward = false;
    }
    cases.push(face);
    let mut face = square(false, 1);
    change(&mut face[0], |p| p[12][2] = 0.001);
    cases.push(face);
    let mut face = square(false, 1);
    change(&mut face[0], |p| {
        // Endpoints alone look like the original cap edge. Its complete
        // polynomial bows across the opposite half-plane and must refuse.
        p[5 * 12][0] = -10.;
        p[5 * 12 + 11][0] = -10.;
    });
    cases.push(face);
    let mut face = square(false, 1);
    change(&mut face[0], |p| {
        // A true interior axial fold, with both original cap stations fixed.
        for row in p.as_chunks_mut::<12>().0 {
            row[5][2] = 3.;
        }
    });
    cases.push(face);
    // The same locally regular annulus wrapping twice must be rejected by the
    // complete cap winding certificate, despite valid local Jacobians.
    cases.push(square(false, 2));
    for faces in cases {
        assert!(run(&faces, o).is_err());
    }
    let mut o = o;
    o.minimum_radius_mm = 2.;
    assert!(run(&square(false, 1), o).is_err());
    let mut thin = square(false, 1);
    for face in &mut thin {
        change(face, |p| {
            for x in p {
                x[0] *= 1e-18;
            }
        });
    }
    // A genuine thin prism, whose near-pi cap angle cannot be resolved by
    // the declared directed atan enclosure. Refuse rather than sample it.
    assert!(
        run(
            &thin,
            Degree11OccupiedOptions {
                minimum_radius_mm: 1e-20,
                ..o
            }
        )
        .is_err()
    );
    assert!(
        certify(
            &square(false, 1),
            [2., 0.],
            o,
            &mut Work {
                remaining: MAX_WORK,
                interrupt: None,
                cancelled: Cell::new(false)
            }
        )
        .is_err()
    );
}
#[test]
fn degree11_occupied_query_limits_and_cancellation_are_bounded_and_retryable() {
    let faces = square(false, 1);
    let o = Degree11OccupiedOptions::default();
    assert!(valid(o));
    assert!(!valid(Degree11OccupiedOptions {
        maximum_depth: MAX_DEPTH + 1,
        ..o
    }));
    assert!(!valid(Degree11OccupiedOptions {
        minimum_radius_mm: f64::NAN,
        ..o
    }));
    assert!(matches!(
        run(
            &faces,
            Degree11OccupiedOptions {
                maximum_work: 1,
                ..o
            }
        ),
        Err(ReadError::Limit)
    ));
    assert!(matches!(
        run(
            &faces,
            Degree11OccupiedOptions {
                maximum_cells: 1,
                ..o
            }
        ),
        Err(ReadError::Limit)
    ));
    let stop = || true;
    let mut w = Work {
        remaining: MAX_WORK,
        interrupt: Some(&stop),
        cancelled: Cell::new(false),
    };
    assert!(certify(&faces, [0., 2.], o, &mut w).is_err());
    assert!(w.cancelled.get());
    assert!(run(&faces, o).is_ok());
}

#[test]
fn genuinely_positive_nonuniform_wall_uses_bounded_dyadic_fallback() {
    let mut faces = square(false, 1);
    change(&mut faces[0], |poles| {
        // Same exact square-prism wall, with a regular nonuniform cubic
        // parametrization. y'=2*((u-.5)^2+.01)/(1/3-.5+.26)>0, yet its
        // elevated degree20 derivative has negative Bernstein coefficients.
        // Whole-cell coefficient positivity must refuse; restriction earns
        // the genuine strict sign without changing the physical source.
        for (i, row) in poles.as_chunks_mut::<12>().0.iter_mut().enumerate() {
            let t = i as f64 / 11.;
            let t2 = if i < 2 {
                0.
            } else {
                (i * (i - 1)) as f64 / 110.
            };
            let t3 = if i < 3 {
                0.
            } else {
                (i * (i - 1) * (i - 2)) as f64 / 990.
            };
            let y = if i == 0 {
                -1.
            } else if i == 11 {
                1.
            } else {
                -1. + 2. * (t3 / 3. - t2 / 2. + 0.26 * t) / (1. / 3. - 0.5 + 0.26)
            };
            for p in row {
                p[1] = y;
            }
        }
    });
    let o = Degree11OccupiedOptions::default();
    let proof = run(&faces, o).unwrap();
    assert!(proof.cells > 4);
    assert!(proof.radius > 0.99);
    assert!(
        run(
            &faces,
            Degree11OccupiedOptions {
                maximum_depth: 0,
                ..o
            }
        )
        .is_err()
    );
}

fn actual_m3(fitted: bool, left: bool) -> String {
    use crate::{
        head_profile::{Head, Profile},
        helical_thread::{Spec, Thread},
        phillips_drive::{Cavity, Cutter, PhillipsScrew},
        profiled_screw::SourceSpec,
    };
    // Exact public Core::ThreadParameters.resolve result captured in
    // cc-m3-retained-receipt.json (CAD352,65e2321b); the flat control is NOT6g.
    // Core remains independent of this crate; no dependency-cycle/graph change.
    let thread = Thread::new(Spec {
        major_radius_mm: if fitted { 1.4635 } else { 1.47 },
        pitch_mm: 0.5,
        length_mm: 8.,
        depth_mm: if fitted { 0.29908756729740676 } else { 0.3 },
        crest_width_mm: if fitted { 0.07131009168581062 } else { 0.0625 },
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: left,
        root_radius_mm: fitted.then_some(0.07216878364870323),
    })
    .unwrap();
    let head = Head::Crowned {
        radius_mm: 2.8,
        height_mm: 2.4,
        crown_radius_mm: 3.,
        top_land_radius_mm: 2.1,
        neck_radius_mm: 0.1,
    };
    let (spec, domain) = crate::phillips_drive::test_fixture(1);
    let source = PhillipsScrew::new(
        thread,
        Profile::new_for_thread_union(head, 1.5, 1.).unwrap(),
        SourceSpec {
            head: Some(head),
            shank_length_mm: 1.,
            minimum_wall_mm: 1e-5,
        },
        Cavity::new(Cutter::new(spec, domain).unwrap(), 1e-5, 100_000).unwrap(),
        100_000,
        &|| false,
    )
    .unwrap();
    let (origin, u, v) = if left {
        ([105., -20., 60.], [0.8, 0.6, 0.], [0.36, -0.48, -0.8])
    } else {
        ([0.; 3], [1., 0., 0.], [0., 1., 0.])
    };
    crate::step::write_phillips_screw(&source, origin, u, v, 1e-7).unwrap()
}
#[test]
fn actual_m3_flat_annulus_and_retained_6g_d11_nets_earn_occupancy_both_hands() {
    use crate::step_read::phillips_cutter::physical::body::common_mouth::tests::derived_thread_probe::with_head;
    for fitted in [false, true] {
        for left in [false, true] {
            let text = actual_m3(fitted, left);
            if !left && let Ok(directory) = std::env::var("MULTICAD_RETAINED_PROBE_OUTPUT") {
                let original = std::fs::read_to_string(
                    std::path::Path::new(&directory)
                        .join(format!("cc-m3-retained-fitted-{fitted}.step")),
                )
                .unwrap();
                assert_eq!(
                    text, original,
                    "same saved actual fit-engine-derived raw file"
                );
            }
            with_head(&text, &mut |head| {
                let carrier = head
                    .degree11_thread_carrier(Default::default(), None)
                    .unwrap();
                let original: Vec<_> = carrier.faces().iter().map(|f| f.source().clone()).collect();
                let inherited = carrier.total_work_used();
                let original_caps = carrier.cap_stations_mm();
                let occupied = carrier.certify_occupied(Default::default(), None).unwrap();
                eprintln!(
                    "D11_OCCUPIED fitted={fitted} left={left} new={} inherited={inherited} total={} peak={} retained={} cells={} core={} jacobian={}",
                    occupied.new_work_used(),
                    occupied.total_work_used(),
                    occupied.accounted_peak_bytes(),
                    occupied.retained_budget_bytes(),
                    occupied.visited_cells(),
                    occupied.occupied_core_radius_mm(),
                    occupied.minimum_leaf_jacobian_mm3()
                );
                assert!(std::ptr::eq(occupied.body(), head.body()));
                assert_eq!(occupied.chart(), head.chart());
                assert_eq!(occupied.cap_stations_mm(), original_caps);
                assert!(original_caps[0].abs() < 1e-7 && (original_caps[1] - 8.).abs() < 1e-7);
                assert_eq!(occupied.body().literal_input(), text);
                assert_eq!(
                    occupied
                        .body()
                        .inferred_thread_spec()
                        .root_radius_mm
                        .is_some(),
                    fitted
                );
                assert!(occupied.total_work_used() <= MAX_TOTAL);
                assert!(occupied.accounted_peak_bytes() <= BYTES);
                for (face, original) in occupied.carrier().faces().iter().zip(&original) {
                    assert_eq!(face.source(), original);
                }
            });
        }
    }
}

#[test]
fn public_d11_occupied_constructor_preflights_exact_limits_cancels_and_retries() {
    use crate::step_read::phillips_cutter::physical::body::common_mouth::tests::derived_thread_probe::with_head;
    let text = actual_m3(true, false);
    with_head(&text, &mut |head| {
        let make = || {
            head.degree11_thread_carrier(Default::default(), None)
                .unwrap()
        };
        let ticks = Cell::new(0usize);
        let observe = || {
            ticks.set(ticks.get() + 1);
            false
        };
        let reference = make()
            .certify_occupied(Default::default(), Some(&observe))
            .unwrap();
        let calls = ticks.get();
        let o = Degree11OccupiedOptions::default();
        for limited in [
            Degree11OccupiedOptions {
                maximum_work: reference.new_work_used() - 1,
                ..o
            },
            Degree11OccupiedOptions {
                maximum_total_work: reference.total_work_used() - 1,
                ..o
            },
            Degree11OccupiedOptions {
                maximum_total_work: make().total_work_used() - 1,
                ..o
            },
            Degree11OccupiedOptions {
                maximum_bytes: reference.accounted_peak_bytes() - 1,
                ..o
            },
            Degree11OccupiedOptions {
                maximum_cells: reference.visited_cells() - 1,
                ..o
            },
        ] {
            assert!(matches!(
                make().certify_occupied(limited, None),
                Err(PhillipsReadError::Limit)
            ));
        }
        for at in [1, calls / 2, calls] {
            let seen = Cell::new(0usize);
            let stop = || {
                seen.set(seen.get() + 1);
                seen.get() == at
            };
            assert!(matches!(
                make().certify_occupied(o, Some(&stop)),
                Err(PhillipsReadError::Cancelled)
            ));
        }
        assert!(matches!(
            make().certify_occupied(
                Degree11OccupiedOptions {
                    maximum_bytes: 0,
                    ..o
                },
                None
            ),
            Err(PhillipsReadError::Invalid)
        ));
        let retry = make()
            .certify_occupied(
                Degree11OccupiedOptions {
                    maximum_work: reference.new_work_used(),
                    maximum_total_work: reference.total_work_used(),
                    maximum_bytes: reference.accounted_peak_bytes(),
                    ..o
                },
                None,
            )
            .unwrap();
        assert_eq!(
            reference.occupied_core_radius_mm(),
            retry.occupied_core_radius_mm()
        );
        assert_eq!(reference.new_work_used(), retry.new_work_used());
        assert_eq!(
            reference.local_enclosure_bounds_mm(),
            retry.local_enclosure_bounds_mm()
        );
        assert!(retry.minimum_leaf_axial_derivative_mm() > 0.);
        assert_eq!(
            reference.cap_winding_angle_bounds_rad(),
            retry.cap_winding_angle_bounds_rad()
        );
        assert_eq!(head.body().literal_input(), text);
    });
}
