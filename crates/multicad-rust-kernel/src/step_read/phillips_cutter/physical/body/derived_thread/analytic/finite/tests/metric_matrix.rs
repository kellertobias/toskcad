//! Existing six REAL metric writer packets, followed through every production
//! prerequisite. A refusal is recorded at its actual stage, never hidden by an
//! unwrap or represented as supported manufacturing/import.
use super::*;
#[derive(Debug)]
struct Failure {
    stage: &'static str,
    error: PhillipsReadError,
}
#[derive(Debug)]
struct Outcome {
    known_work: usize,
    finite_work: usize,
    retained: usize,
    field_bound: [f64; 3],
    curve_bound: [f64; 2],
    kinds: [usize; 2],
}
fn failure(stage: &'static str) -> impl FnOnce(PhillipsReadError) -> Failure {
    move |error| Failure { stage, error }
}
fn qualify(text: &str) -> std::result::Result<Outcome, Failure> {
    let body = read_phillips_body(text, 1e-7, 1_000_000).map_err(failure("parse"))?;
    let outer = body
        .outer_carrier_bounds(
            OuterCarrierOptions {
                maximum_distance_mm: 1e-7,
                work_limit: 10_000_000,
            },
            None,
        )
        .map_err(failure("outer"))?;
    let cavity = body.cavity();
    let ends = cavity
        .endpoint_coupling_bounds(
            EndpointCouplingOptions {
                maximum_distance_mm: 1e-7,
                work_limit: 1_000_000,
            },
            None,
        )
        .map_err(failure("cavity-endpoints"))?;
    let analytic = cavity
        .analytic_support_bounds(
            &ends,
            AnalyticSupportOptions {
                maximum_distance_mm: 1e-7,
                work_limit: 10_000_000,
            },
            None,
        )
        .map_err(failure("cavity-supports"))?;
    let lifts = analytic
        .finite_edge_lifts(ends, AnalyticEdgeLiftOptions::default(), None)
        .map_err(failure("cavity-lifts"))?;
    let cubics = cavity
        .cubic_support_bounds(
            CubicSupportOptions {
                maximum_distance_mm: 1e-7,
                work_limit: 10_000_000,
            },
            None,
        )
        .map_err(failure("cavity-cubics"))?;
    let sewn = lifts
        .sewn_edge_sources(&body, &cubics, Default::default(), None)
        .map_err(failure("cavity-sewn"))?;
    let rep = body
        .derived_cavity_representative(&sewn, Default::default(), None)
        .map_err(failure("cavity-fit"))?;
    let t = rep.transfer_bounds();
    assert_eq!(
        (t.vertex_count(), t.edge_count(), t.face_count()),
        (49, 76, 29)
    );
    let cavity_bound = [
        t.maximum_vertex_error_mm(),
        t.maximum_analytic_support_error_mm(),
        t.maximum_polynomial_error_mm(),
        t.maximum_polynomial_support_error_mm(),
    ];
    assert!(cavity_bound.into_iter().all(|x| x <= 1e-7));
    let profile = body.retained_meridian_header();
    eprintln!(
        "METRIC CAVITY FIT inner{:?} actual-transfer{:?} head-height{} cavity-tip{:?} work{}",
        rep.derived_cavity_source().spec().inner,
        cavity_bound,
        profile.end()[1] - profile.top()[1],
        rep.derived_cavity_source().envelope().tip_depth_mm_bounds(),
        rep.work_used()
    );
    let closed = rep
        .closed_boundary_certificate(Default::default(), None)
        .map_err(failure("cavity-closed"))?;
    let common = closed
        .common_mouth_head_envelope(&outer, Default::default(), None)
        .map_err(failure("common-chart"))?;
    let head = common
        .closed_repaired_head(Default::default(), None)
        .map_err(failure("head-closed"))?;
    let field = head
        .common_chart_periodic_thread_repair(Default::default(), None)
        .map_err(failure("periodic-field"))?;
    let finite = head
        .periodic_thread_finite_curves(&field, Default::default(), None)
        .map_err(failure("finite-curves"))?;
    let field_bound = field.transfers().iter().fold([0_f64; 3], |a, b| {
        [
            a[0].max(b.distance_mm),
            a[1].max(b.derivative_mm[0]),
            a[2].max(b.derivative_mm[1]),
        ]
    });
    let curve_bound = finite.curves().iter().fold([0_f64; 2], |a, b| {
        [a[0].max(b.distance_mm), a[1].max(b.cell_derivative_mm)]
    });
    let kinds = finite.curves().iter().fold([0; 2], |mut a, c| {
        a[usize::from(c.kind == PeriodicFiniteKind::Circle)] += 1;
        a
    });
    Ok(Outcome {
        known_work: finite.total_work_used(),
        finite_work: finite.new_work_used(),
        retained: finite.retained_budget_bytes(),
        field_bound,
        curve_bound,
        kinds,
    })
}
#[test]
fn existing_six_metric_file_support_matrix_records_actual_supported_and_refused_stages() {
    let mut supported = 0;
    let mut refused = 0;
    for index in 0..6 {
        for left in [false, true] {
            for nominal in [false, true] {
                let source = super::super::super::super::super::super::tests::fixture(
                    index, left, nominal, 1,
                );
                let spec = source.thread().spec();
                assert_eq!(spec.starts, 1);
                assert_eq!(spec.root_radius_mm, None);
                let (o, u, v) = if left {
                    ([105., -20., 60.], [0.8, 0.6, 0.], [0.36, -0.48, -0.8])
                } else {
                    ([0.; 3], [1., 0., 0.], [0., 1., 0.])
                };
                let text = crate::step::write_phillips_screw(&source, o, u, v, 1e-7).unwrap();
                let start = std::time::Instant::now();
                match qualify(&text) {
                    Ok(outcome) => {
                        assert_ne!(
                            index, 2,
                            "countersunk must use a future combined host qualifier"
                        );
                        supported += 1;
                        assert!(
                            outcome
                                .field_bound
                                .into_iter()
                                .chain(outcome.curve_bound)
                                .all(|x| x <= 1e-7)
                        );
                        assert!(
                            outcome.known_work <= FINITE_TOTAL
                                && outcome.finite_work <= WORK
                                && outcome.retained < BYTES
                        );
                        assert_eq!(outcome.kinds.iter().sum::<usize>(), 16);
                        eprintln!(
                            "METRIC MATRIX index{index} left{left} nominal{nominal} thread{spec:?} head{:?} bytes{} SUPPORTED {outcome:?} elapsed{:?}",
                            source.spec().head,
                            text.len(),
                            start.elapsed()
                        );
                    }
                    Err(f) => {
                        assert_eq!(index, 2, "unexpected unsupported metric source: {f:?}");
                        assert_eq!(f.stage, "common-chart");
                        assert_eq!(f.error, PhillipsReadError::Unsupported);
                        refused += 1;
                        eprintln!(
                            "METRIC MATRIX index{index} left{left} nominal{nominal} thread{spec:?} head{:?} bytes{} REFUSED stage{} {:?} elapsed{:?}",
                            source.spec().head,
                            text.len(),
                            f.stage,
                            f.error,
                            start.elapsed()
                        );
                    }
                }
            }
        }
    }
    eprintln!("METRIC MATRIX TOTAL supported{supported} refused{refused}");
    assert_eq!((supported, refused), (20, 4));
}
