use super::*;
use crate::axial_pipe_thread::authored::{AuthoredSection, AuthoredSpec, AuthoredThread};
fn source(quarter: bool) -> StockClipPlan {
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
fn evaluate<const N: usize>(mut poles: [[f64; 3]; N], t: f64) -> [f64; 3] {
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
