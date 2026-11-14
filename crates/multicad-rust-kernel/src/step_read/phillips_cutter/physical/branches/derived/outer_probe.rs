//! Diagnostic only: a fitted common mouth is not the original exact raw plane.
use super::*;
use num_rational::BigRational as Q;
use num_traits::ToPrimitive;
fn decimal(token: &str) -> Q {
    let (m, e) = token.split_once('E').unwrap();
    let (a, b) = m.split_once('.').unwrap();
    let n: num_bigint::BigInt = format!("{a}{b}").parse().unwrap();
    let exponent = e.parse::<i32>().unwrap() - b.len() as i32;
    if exponent >= 0 {
        Q::from_integer(n * num_bigint::BigInt::from(10).pow(exponent as u32))
    } else {
        Q::new(n, num_bigint::BigInt::from(10).pow((-exponent) as u32))
    }
}
fn exact(p: Point) -> [Q; 3] {
    p.map(|x| Q::from_float(x).unwrap())
}
fn dot(a: &[Q; 3], b: &[Q; 3]) -> Q {
    (0..3).fold(Q::from_integer(0.into()), |s, i| s + &a[i] * &b[i])
}
fn cross(a: &[Q; 3], b: &[Q; 3]) -> [Q; 3] {
    std::array::from_fn(|i| &a[(i + 1) % 3] * &b[(i + 2) % 3] - &a[(i + 2) % 3] * &b[(i + 1) % 3])
}
#[test]
fn actual_outer_and_cavity_mouth_frames_require_explicit_derived_authority() {
    for scale in [1., 1.137] {
        for posed in [false, true] {
            let text = super::super::representative::scaled_text(scale, posed);
            let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
            let outer = body
                .outer_carrier_bounds(
                    OuterCarrierOptions {
                        maximum_distance_mm: 1e-7,
                        work_limit: 10_000_000,
                    },
                    None,
                )
                .unwrap();
            assert!(std::ptr::eq(outer.body(), &body));
            let b = body.cavity();
            let endpoints = b
                .endpoint_coupling_bounds(
                    EndpointCouplingOptions {
                        maximum_distance_mm: 1e-7,
                        work_limit: 1_000_000,
                    },
                    None,
                )
                .unwrap();
            let analytic = b
                .analytic_support_bounds(
                    &endpoints,
                    AnalyticSupportOptions {
                        maximum_distance_mm: 1e-7,
                        work_limit: 10_000_000,
                    },
                    None,
                )
                .unwrap();
            let lifts = analytic
                .finite_edge_lifts(endpoints, LiftOptions::default(), None)
                .unwrap();
            let cubics = b
                .cubic_support_bounds(
                    CubicSupportOptions {
                        maximum_distance_mm: 1e-7,
                        work_limit: 10_000_000,
                    },
                    None,
                )
                .unwrap();
            let sewn = lifts
                .sewn_edge_sources(&body, &cubics, SewnOptions::default(), None)
                .unwrap();
            let derived = body
                .derived_cavity_representative(&sewn, Default::default(), None)
                .unwrap();
            let closed = derived
                .closed_boundary_certificate(Default::default(), None)
                .unwrap();
            let common = closed
                .common_mouth_head_envelope(&outer, Default::default(), None)
                .unwrap();
            let columns = derived.derived_affine_columns();
            let [o, u, v, n] = columns.map(exact);
            let face = b
                .faces()
                .iter()
                .find(|f| f.material_entity.is_none())
                .unwrap();
            let entities = body.source_entities();
            let surface = &entities[&face.surface_entity][0];
            let axis = reference(&surface.args[1]).unwrap();
            let axis = &entities[&axis][0];
            let point = body
                .exact_point_leaf(reference(&axis.args[1]).unwrap(), 1_000_000, None)
                .unwrap();
            let po = point.coordinate_spellings().map(decimal);
            let direction = &entities[&reference(&axis.args[2]).unwrap()][0];
            let words = part21::numeric_spellings(direction).unwrap();
            let pn = std::array::from_fn(|i| decimal(words[i]));
            let delta = std::array::from_fn(|i| &o[i] - &po[i]);
            let raw_plane = dot(&pn, &delta);
            let raw_u = dot(&pn, &u);
            let raw_v = dot(&pn, &v);
            if scale == 1. && !posed {
                assert_ne!(raw_plane, Q::from_integer(0.into()));
            }
            // Independent exact-decimal oracle over all four corners of the
            // enclosing source mouth square. Compare squares, without a rounded
            // normal normalization or floating square root in the oracle.
            let bound = Q::from_float(common.raw_mouth_plane_distance_mm()).unwrap();
            let radius = Q::from_float(
                derived
                    .derived_cavity_source()
                    .envelope()
                    .entry_radius_mm_bounds()[1],
            )
            .unwrap();
            let norm_squared = dot(&pn, &pn);
            for a in [-1, 1] {
                for b in [-1, 1] {
                    let residual = &raw_plane
                        + Q::from_integer(a.into()) * &radius * &raw_u
                        + Q::from_integer(b.into()) * &radius * &raw_v;
                    assert!(&residual * &residual <= &bound * &bound * &norm_squared);
                }
            }
            let [oo, ou, ov] = body.inferred_thread_frame().map(exact);
            let on = cross(&ou, &ov);
            let top = Q::from_float(body.retained_meridian_header().top()[1]).unwrap();
            let proposal_delta: [Q; 3] = std::array::from_fn(|i| &o[i] - (&oo[i] + &on[i] * &top));
            let frame_delta = std::array::from_fn::<_, 3, _>(|i| &n[i] - &on[i]);
            let scalar = |q: &Q| q.to_f64().unwrap();
            eprintln!(
                "scale{scale} posed{posed} rawPlaneDot={} rawDU={} rawDV={} bodyProposalMouthDelta={:?} Ndelta={:?} outerMax={} profileReserve={}",
                scalar(&raw_plane),
                scalar(&raw_u),
                scalar(&raw_v),
                proposal_delta.map(|q| scalar(&q)),
                frame_delta.map(|q| scalar(&q)),
                outer
                    .tensors()
                    .iter()
                    .map(|t| t.control_deviation_mm)
                    .fold(0., f64::max),
                outer.meridian_reconstruction_error_mm()
            );
            assert_eq!(derived.body().literal_input(), text);
        }
    }
}
#[test]
fn decimal_oracle_preserves_valid_real_spellings_across_exponent_signs() {
    assert_eq!(decimal("1.25E+3"), Q::from_integer(1250.into()));
    assert_eq!(decimal("1.250E+3"), Q::from_integer(1250.into()));
    assert_eq!(decimal("-1.25E-3"), Q::new((-1).into(), 800.into()));
    assert_eq!(decimal("3.125E-2"), Q::new(1.into(), 32.into()));
    assert_eq!(decimal("0.0E+4"), Q::from_integer(0.into()));
}
