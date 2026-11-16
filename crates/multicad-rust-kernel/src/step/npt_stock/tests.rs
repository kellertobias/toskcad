use super::*;
use crate::axial_pipe_thread::authored::{AuthoredSection, AuthoredSpec, AuthoredThread};
pub(super) fn source(quarter: bool) -> StockClipPlan {
    let (tpi, e1, l1, l2, od, max) = if quarter {
        (18., 0.49163, 0.2278, 0.4018, 0.540, 0.088)
    } else {
        (27., 0.37360, 0.1615, 0.2639, 0.405, 0.096)
    };
    let pitch = 25.4 / tpi;
    let k = 1. / 32.;
    let m = 3f64.sqrt();
    let slopes = 1. / (m - k) + 1. / (m + k);
    let depth = pitch / slopes - (0.033 + max) * pitch;
    AuthoredThread::new(
        AuthoredSpec {
            major_radius_mm: (e1 * 25.4 - 2. * k * l1 * 25.4 + depth) / 2.,
            pitch_mm: pitch,
            length_mm: l2 * 25.4,
            radial_taper: k,
            section: AuthoredSection::NptFlat {
                depth_mm: depth,
                flat_width_mm: (pitch - depth * slopes) / 2.,
            },
            crest_radius_limit_mm: Some(od * 25.4 / 2.),
        },
        1e-7,
        1e-6,
    )
    .unwrap()
    .stock_clip_plan(Default::default(), None)
    .unwrap()
}
pub(super) fn evaluate<const N: usize>(mut poles: [[f64; 3]; N], t: f64) -> [f64; 3] {
    for level in (1..N).rev() {
        for i in 0..level {
            let next = poles[i + 1];
            for (j, value) in poles[i].iter_mut().enumerate() {
                *value = (1. - t) * *value + t * next[j];
            }
        }
    }
    poles[0]
}
#[test]
fn actual_twenty_edge_degree_eleven_candidates_share_vertices_and_bound_source_c0_c1() {
    for quarter in [false, true] {
        let source = source(quarter);
        let before = *source.source();
        let contact = source.contact_height_bounds_mm();
        let options = CurveOptions::default();
        let plan = CurvePlan::new(&source, options, None).unwrap();
        eprintln!(
            "NPT CURVE CANDIDATE quarter={quarter} spans={} work={} bytes={}",
            plan.curves().iter().map(|c| c.spans().len()).sum::<usize>(),
            plan.work_used(),
            plan.retained_budget_bytes()
        );
        assert!(std::ptr::eq(plan.source(), &source));
        assert_eq!(plan.curves().len(), 20);
        assert!(plan.contact_plane_error_mm() > 0.);
        assert!(plan.contact_plane_error_mm() < options.position_error_mm);
        assert!(
            contact[0] <= plan.contact_plane_height_mm()
                && plan.contact_plane_height_mm() <= contact[1]
        );
        for (curve, key) in plan.curves().iter().zip(source.curves()) {
            assert_eq!(curve.source_curve(), key);
            let law = source.curve_law(key, Default::default(), None).unwrap();
            let spans = curve.spans();
            assert_eq!(spans[0].domain()[0], 0.);
            assert_eq!(spans.last().unwrap().domain()[1], 1.);
            for ends in spans.windows(2) {
                assert_eq!(ends[0].domain()[1], ends[1].domain()[0]);
                assert_eq!(ends[0].endpoint_owners()[1], ends[1].endpoint_owners()[0]);
                assert_eq!(ends[0].poles_mm()[11], ends[1].poles_mm()[0]);
            }
            for (i, key) in key.endpoints().into_iter().enumerate() {
                let vertex = plan.vertices().iter().find(|v| v.owner() == key).unwrap();
                let span = if i == 0 {
                    &spans[0]
                } else {
                    spans.last().unwrap()
                };
                assert_eq!(span.endpoint_owners()[i], CurvePointOwner::Vertex(key));
                assert_eq!(span.poles_mm()[i * 11], vertex.position_mm());
            }
            for span in spans {
                assert!(span.position_error_mm() <= options.position_error_mm);
                assert!(span.original_derivative_error_mm() <= options.derivative_error_mm);
                if let StockClipCurve::Ring { station, .. } = key {
                    let z = match station {
                        crate::axial_pipe_thread::authored::StockClipStation::Start => 0.,
                        crate::axial_pipe_thread::authored::StockClipStation::Contact => {
                            plan.contact_plane_height_mm()
                        }
                        crate::axial_pipe_thread::authored::StockClipStation::End => {
                            before.spec().length_mm
                        }
                    };
                    assert!(span.poles_mm().iter().all(|p| p[2] == z));
                }
                let [a, b] = span.domain();
                let h = b - a;
                let derivative: [[f64; 3]; 11] = std::array::from_fn(|i| {
                    std::array::from_fn(|j| {
                        11. * (span.poles_mm()[i + 1][j] - span.poles_mm()[i][j]) / h
                    })
                });
                for local in [0., 0.13, 0.5, 0.91, 1.] {
                    let t = a + h * local;
                    let source = law
                        .taylor([t, t], Default::default(), None)
                        .unwrap()
                        .coefficients_mm();
                    let point = evaluate(*span.poles_mm(), local);
                    let d = evaluate(derivative, local);
                    for j in 0..3 {
                        let p = source[0][j][0].midpoint(source[0][j][1]);
                        let slope = source[1][j][0].midpoint(source[1][j][1]);
                        assert!((point[j] - p).abs() <= span.position_error_mm() + 1e-11);
                        assert!(
                            (d[j] - slope).abs() <= span.original_derivative_error_mm() + 1e-10
                        );
                    }
                }
            }
        }
        assert_eq!(*source.source(), before);
        assert_eq!(source.contact_height_bounds_mm(), contact);
    }
}

#[test]
fn candidate_errors_remain_source_only_and_preserve_whole_request_limits() {
    let source = source(false);
    let options = CurveOptions::default();
    let result = CurvePlan::new(&source, options, None).unwrap();
    // Exhaust exactly the first real source-vertex query, so the next query
    // sees zero remaining work rather than a renewed/invalid source budget.
    let first_vertex_work = source
        .vertex_bounds(source.vertices()[0], Default::default(), None)
        .unwrap()
        .work_used();
    for bad in [
        CurveOptions {
            position_error_mm: 0.,
            ..options
        },
        CurveOptions {
            derivative_error_mm: f64::NAN,
            ..options
        },
        CurveOptions {
            maximum_spans_per_curve: 3,
            ..options
        },
        CurveOptions {
            maximum_spans_per_curve: 512,
            ..options
        },
        CurveOptions {
            maximum_work: 0,
            ..options
        },
        CurveOptions {
            maximum_work: options.maximum_work + 1,
            ..options
        },
        CurveOptions {
            maximum_bytes: 0,
            ..options
        },
        CurveOptions {
            maximum_bytes: options.maximum_bytes + 1,
            ..options
        },
    ] {
        assert!(matches!(
            CurvePlan::new(&source, bad, None),
            Err(Error::Invalid)
        ));
    }
    for bad in [
        CurveOptions {
            maximum_work: 1,
            ..options
        },
        CurveOptions {
            maximum_work: first_vertex_work,
            ..options
        },
        CurveOptions {
            maximum_work: 25_000,
            ..options
        },
        CurveOptions {
            maximum_work: result.work_used() - 1,
            ..options
        },
        CurveOptions {
            maximum_bytes: 1,
            ..options
        },
        CurveOptions {
            maximum_bytes: result.retained_budget_bytes() - 1,
            ..options
        },
    ] {
        assert!(matches!(
            CurvePlan::new(&source, bad, None),
            Err(Error::Limit)
        ));
    }
    assert_eq!(
        CurvePlan::new(
            &source,
            CurveOptions {
                maximum_work: result.work_used(),
                maximum_bytes: result.retained_budget_bytes(),
                ..options
            },
            None
        )
        .unwrap()
        .work_used(),
        result.work_used()
    );
    assert!(matches!(
        CurvePlan::new(
            &source,
            CurveOptions {
                maximum_spans_per_curve: 1,
                ..options
            },
            None
        ),
        Err(Error::Unresolved)
    ));
    assert!(matches!(
        CurvePlan::new(&source, options, Some(&|| true)),
        Err(Error::Cancelled)
    ));
    let polls = std::cell::Cell::new(0usize);
    let count = || {
        polls.set(polls.get() + 1);
        false
    };
    CurvePlan::new(&source, options, Some(&count)).unwrap();
    let seen = std::cell::Cell::new(0usize);
    let cancel = || {
        let n = seen.get();
        seen.set(n + 1);
        n + 1 == polls.get()
    };
    assert!(matches!(
        CurvePlan::new(&source, options, Some(&cancel)),
        Err(Error::Cancelled)
    ));
}

#[test]
fn actual_shared_surface_candidates_keep_original_source_and_whole_work_quota() {
    for quarter in [false, true] {
        let source = source(quarter);
        let original = *source.source();
        let contact = source.contact_height_bounds_mm();
        let options = CurveOptions::default();
        let curves = CurvePlan::new(&source, options, None).unwrap();
        let result = SurfacePlan::new(curves, options.into(), None).unwrap();
        eprintln!(
            "NPT SURFACE COMPLETE quarter={quarter} work={} bytes={} controls={}",
            result.work_used(),
            result.retained_budget_bytes(),
            result
                .patches()
                .iter()
                .map(|p| p.poles_mm().len())
                .sum::<usize>()
        );
        assert_eq!(result.patches().len(), 8);
        assert!(result.work_used() <= options.maximum_work);
        for patch in result.patches() {
            let law = source
                .patch_law(patch.source_patch(), Default::default(), None)
                .unwrap();
            let dims = patch.spans().map(|n| 11 * n + 1);
            assert_eq!(patch.poles_mm().len(), dims[0] * dims[1]);
            assert_eq!(patch.cells().len(), patch.spans()[0] * patch.spans()[1]);
            for (side, edge) in patch.source_patch().coedges().into_iter().enumerate() {
                let curve = result
                    .curves()
                    .curves()
                    .iter()
                    .find(|c| c.source_curve() == edge.curve)
                    .unwrap();
                let count = if side % 2 == 0 { dims[1] } else { dims[0] };
                for index in 0..count {
                    let span = (index / 11).min(curve.spans().len() - 1);
                    let expected = curve.spans()[span].poles_mm()[index - 11 * span];
                    let [u, v] = match side {
                        0 => [0, index],
                        1 => [index, dims[1] - 1],
                        2 => [dims[0] - 1, index],
                        _ => [index, 0],
                    };
                    assert_eq!(patch.poles_mm()[u * dims[1] + v], expected);
                }
            }
            for (index, cell) in patch.cells().iter().enumerate() {
                assert!(cell.position_error_mm() <= options.position_error_mm);
                assert!(
                    cell.original_derivative_error_mm()
                        .iter()
                        .all(|e| *e <= options.derivative_error_mm)
                );
                let i = index / patch.spans()[1];
                let j = index % patch.spans()[1];
                for [u, v] in [[0., 0.], [0.13, 0.29], [0.5, 0.5], [0.91, 0.83], [1., 1.]] {
                    let rows: [[f64; 3]; 12] = std::array::from_fn(|a| {
                        evaluate(
                            std::array::from_fn::<_, 12, _>(|b| {
                                patch.poles_mm()[(11 * i + a) * dims[1] + 11 * j + b]
                            }),
                            v,
                        )
                    });
                    let actual = evaluate(rows, u);
                    let domain = cell.domain();
                    let uv = [
                        domain[0][0] + u * (domain[0][1] - domain[0][0]),
                        domain[1][0] + v * (domain[1][1] - domain[1][0]),
                    ];
                    let jets = law
                        .jets([uv[0], uv[0]], [uv[1], uv[1]], Default::default(), None)
                        .unwrap();
                    let expected = jets.sample().point().position_mm();
                    let du = evaluate(
                        std::array::from_fn::<_, 11, _>(|a| {
                            std::array::from_fn(|axis| {
                                11. * (rows[a + 1][axis] - rows[a][axis])
                                    / (domain[0][1] - domain[0][0])
                            })
                        }),
                        u,
                    );
                    let dv_rows = std::array::from_fn::<_, 12, _>(|a| {
                        evaluate(
                            std::array::from_fn::<_, 11, _>(|b| {
                                std::array::from_fn(|axis| {
                                    11. * (patch.poles_mm()
                                        [(11 * i + a) * dims[1] + 11 * j + b + 1][axis]
                                        - patch.poles_mm()[(11 * i + a) * dims[1] + 11 * j + b]
                                            [axis])
                                        / (domain[1][1] - domain[1][0])
                                })
                            }),
                            v,
                        )
                    });
                    let dv = evaluate(dv_rows, u);
                    for axis in 0..3 {
                        let middle = expected[axis][0].midpoint(expected[axis][1]);
                        assert!((actual[axis] - middle).abs() <= cell.position_error_mm() + 1e-10);
                        for (direction, (observed, expected)) in
                            [(du, jets.du()), (dv, jets.dv())].into_iter().enumerate()
                        {
                            let middle = expected[axis][0].midpoint(expected[axis][1]);
                            assert!(
                                (observed[axis] - middle).abs()
                                    <= cell.original_derivative_error_mm()[direction] + 1e-10
                            );
                        }
                    }
                }
            }
        }
        assert_eq!(*source.source(), original);
        assert_eq!(source.contact_height_bounds_mm(), contact);
    }
}

#[test]
fn shared_surface_failures_keep_real_source_and_one_aggregate_resource_boundary() {
    let source = source(false);
    let original = *source.source();
    let options = CurveOptions::default();
    let make = || CurvePlan::new(&source, options, None).unwrap();
    let result = SurfacePlan::new(make(), options.into(), None).unwrap();
    assert!(std::ptr::eq(result.curves().source(), &source));
    for bad in [
        CurveOptions {
            maximum_spans_per_curve: 1,
            ..options
        },
        CurveOptions {
            maximum_work: result.work_used() - 1,
            ..options
        },
        CurveOptions {
            maximum_bytes: result.retained_budget_bytes() - 1,
            ..options
        },
        CurveOptions {
            maximum_work: 1,
            ..options
        },
        CurveOptions {
            maximum_bytes: 1,
            ..options
        },
    ] {
        assert!(matches!(
            SurfacePlan::new(make(), bad.into(), None),
            Err(Error::Limit)
        ));
    }
    assert_eq!(
        SurfacePlan::new(
            make(),
            CurveOptions {
                maximum_work: result.work_used(),
                maximum_bytes: result.retained_budget_bytes(),
                ..options
            }
            .into(),
            None
        )
        .unwrap()
        .work_used(),
        result.work_used()
    );
    assert!(matches!(
        SurfacePlan::new(
            make(),
            CurveOptions {
                position_error_mm: 0.,
                ..options
            }
            .into(),
            None
        ),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        SurfacePlan::new(
            make(),
            CurveOptions {
                position_error_mm: 1e-14,
                derivative_error_mm: 1e-14,
                ..options
            }
            .into(),
            None
        ),
        Err(Error::Unresolved)
    ));
    assert!(matches!(
        SurfacePlan::new(make(), options.into(), Some(&|| true)),
        Err(Error::Cancelled)
    ));
    let polls = std::cell::Cell::new(0usize);
    let count = || {
        polls.set(polls.get() + 1);
        false
    };
    SurfacePlan::new(make(), options.into(), Some(&count)).unwrap();
    let seen = std::cell::Cell::new(0usize);
    let cancel = || {
        seen.set(seen.get() + 1);
        seen.get() == polls.get()
    };
    assert!(matches!(
        SurfacePlan::new(make(), options.into(), Some(&cancel)),
        Err(Error::Cancelled)
    ));
    let controls = result
        .patches()
        .iter()
        .map(|p| p.poles_mm().len())
        .sum::<usize>();
    assert!(matches!(
        SurfacePlan::new(
            make(),
            SurfaceOptions {
                query: options,
                maximum_control_points: controls - 1
            },
            None
        ),
        Err(Error::Limit)
    ));
    for maximum_control_points in [0, 262_145] {
        assert!(matches!(
            SurfacePlan::new(
                make(),
                SurfaceOptions {
                    query: options,
                    maximum_control_points
                },
                None
            ),
            Err(Error::Invalid)
        ));
    }
    assert_eq!(
        SurfacePlan::new(
            make(),
            SurfaceOptions {
                query: options,
                maximum_control_points: controls
            },
            None
        )
        .unwrap()
        .patches()
        .len(),
        8
    );
    assert_eq!(SurfaceOptions::default().maximum_control_points, 262_144);
    assert_eq!(*source.source(), original);
}

#[test]
fn actual_source_known_closed_step_packet_keeps_one_registry_and_original_total_limits() {
    use crate::step_read::owned_insert::{CompoundDocument, DocumentOptions};
    for quarter in [false, true] {
        let source = source(quarter);
        let options = CurveOptions::default();
        let curves = CurvePlan::new(&source, options, None).unwrap();
        let plan = SurfacePlan::new(curves, options.into(), None).unwrap();
        let file = write_source_candidate(&plan, options.maximum_work, options.maximum_bytes, None)
            .unwrap();
        eprintln!(
            "NPT FILE quarter={quarter} bytes={} work={} peak={}",
            file.as_str().len(),
            file.receipt().work_used(),
            file.receipt().accounted_peak_bytes()
        );
        assert_eq!(file.as_str().matches("=CLOSED_SHELL(").count(), 1);
        assert_eq!(file.as_str().matches("=MANIFOLD_SOLID_BREP(").count(), 1);
        assert_eq!(file.as_str().matches("=ADVANCED_FACE(").count(), 10);
        assert_eq!(file.as_str().matches("=EDGE_CURVE(").count(), 20);
        assert_eq!(file.as_str().matches("=PCURVE(").count(), 40);
        let document =
            CompoundDocument::parse(file.as_str(), DocumentOptions::default(), None).unwrap();
        assert_eq!(document.represented_members().len(), 1);
        let again =
            write_source_candidate(&plan, options.maximum_work, options.maximum_bytes, None)
                .unwrap()
                .into_string();
        assert_eq!(file.as_str(), again);
    }
}

#[test]
fn source_known_step_output_refuses_real_budgets_and_cancellation_atomically() {
    let source = source(false);
    let original = *source.source();
    let options = CurveOptions::default();
    let curves = CurvePlan::new(&source, options, None).unwrap();
    let plan = SurfacePlan::new(curves, options.into(), None).unwrap();
    let file =
        write_source_candidate(&plan, options.maximum_work, options.maximum_bytes, None).unwrap();
    let receipt = file.receipt();
    for (work, bytes) in [
        (0, options.maximum_bytes),
        (options.maximum_work + 1, options.maximum_bytes),
        (options.maximum_work, 0),
        (options.maximum_work, options.maximum_bytes + 1),
    ] {
        assert!(matches!(
            write_source_candidate(&plan, work, bytes, None),
            Err(Error::Invalid)
        ));
    }
    for (work, bytes) in [
        (1, options.maximum_bytes),
        (receipt.work_used() - 1, options.maximum_bytes),
        (options.maximum_work, 1),
    ] {
        assert!(matches!(
            write_source_candidate(&plan, work, bytes, None),
            Err(Error::Limit)
        ));
    }
    assert_eq!(
        write_source_candidate(&plan, receipt.work_used(), options.maximum_bytes, None)
            .unwrap()
            .as_str(),
        file.as_str()
    );
    assert!(matches!(
        write_source_candidate(
            &plan,
            options.maximum_work,
            options.maximum_bytes,
            Some(&|| true)
        ),
        Err(Error::Cancelled)
    ));
    let polls = std::cell::Cell::new(0usize);
    let count = || {
        polls.set(polls.get() + 1);
        false
    };
    write_source_candidate(
        &plan,
        options.maximum_work,
        options.maximum_bytes,
        Some(&count),
    )
    .unwrap();
    let seen = std::cell::Cell::new(0usize);
    let cancel = || {
        seen.set(seen.get() + 1);
        seen.get() == polls.get()
    };
    assert!(matches!(
        write_source_candidate(
            &plan,
            options.maximum_work,
            options.maximum_bytes,
            Some(&cancel)
        ),
        Err(Error::Cancelled)
    ));
    // Find an ACTUAL neighboring combined byte/work frontier. Allocation
    // growth can adapt to a smaller quota, so the default peak-minus-one is
    // deliberately not assumed to be a failing request.
    let mut bad = 1;
    let mut good = receipt.accounted_peak_bytes();
    while good - bad > 1 {
        let middle = bad + (good - bad) / 2;
        match write_source_candidate(&plan, options.maximum_work, middle, None) {
            Ok(candidate) => {
                assert_eq!(candidate.as_str(), file.as_str());
                good = middle
            }
            Err(Error::Limit) => bad = middle,
            Err(other) => panic!("unexpected resource result {other:?}"),
        }
    }
    assert!(matches!(
        write_source_candidate(&plan, options.maximum_work, bad, None),
        Err(Error::Limit)
    ));
    assert_eq!(
        write_source_candidate(&plan, options.maximum_work, good, None)
            .unwrap()
            .as_str(),
        file.as_str()
    );
    let existing = || {
        let mut writer = crate::step::Writer::new();
        writer
            .append_box(
                [5, 7, 9].map(|n| crate::Rational::new(n, 1).unwrap()),
                [30., 0., 0.],
                [1., 0., 0.],
                [0., 1., 0.],
            )
            .unwrap();
        writer
    };
    let baseline = existing().finish().unwrap();
    let mut writer = existing();
    assert!(matches!(
        writer.append_source_npt_candidate(
            &plan,
            plan.work_used() + 8192,
            options.maximum_bytes,
            None
        ),
        Err(Error::Limit)
    ));
    assert_eq!(writer.finish().unwrap(), baseline);
    assert_eq!(*source.source(), original);
}
