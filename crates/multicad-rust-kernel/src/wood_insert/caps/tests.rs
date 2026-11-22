use super::*;
use multicad_precision::PortableMath as _;
use std::cell::Cell;
fn options() -> CapOptions {
    CapOptions {
        maximum_work: 4_000_000,
        maximum_bytes: 1024 * 1024,
    }
}
fn fixture(left: bool, starts: u32, nb: bool) -> Insert {
    super::super::tests::fixture(left, starts, nb)
}
fn nominal(source: &AuthoredThread, turn: f64, z: f64) -> f64 {
    let lead = source.spec().pitch_mm * f64::from(source.spec().starts);
    let phase = (z - (if source.spec().left_hand { -1. } else { 1. }) * lead * turn)
        .rem_euclid(source.spec().pitch_mm);
    let bands = source.bands();
    let band = bands.iter().find(|b| phase <= b.phase_mm[1]).unwrap();
    let u = (phase - band.phase_mm[0]) / (band.phase_mm[1] - band.phase_mm[0]);
    band.radius_mm[0] + u * (band.radius_mm[1] - band.radius_mm[0])
}
#[test]
fn actual_source_knots_remain_independent_and_cells_cover_complete_annulus() {
    for left in [false, true] {
        for starts in [1, 2, 8] {
            for nb in [false, true] {
                for end in [CapEnd::Start, CapEnd::End] {
                    let body = fixture(left, starts, nb);
                    let cap = body.cap_domains(end, options(), None).unwrap();
                    assert_eq!(cap.owner(), &body);
                    assert_eq!(cap.end(), end);
                    assert_eq!(cap.knots().len(), 10 * starts as usize);
                    assert!(cap.retained_bytes() < 128 * 1024);
                    assert!(cap.construction_work() > 0);
                    let mut w = work(options(), None).unwrap();
                    assert_eq!(cap.angles[0].decode(&mut w).unwrap(), Q::zero());
                    assert_eq!(cap.angles.last().unwrap().decode(&mut w).unwrap(), Q::one());
                    for pair in cap.angles.windows(2) {
                        assert!(pair[0].decode(&mut w).unwrap() < pair[1].decode(&mut w).unwrap());
                    }
                    let z = if end == CapEnd::Start {
                        0.
                    } else {
                        body.full_length_mm()
                    };
                    for knot in cap.knots() {
                        let source = if knot.side == Side::Outer {
                            body.outer()
                        } else {
                            body.inner()
                        };
                        let q = cap.angles[knot.event].decode(&mut w).unwrap();
                        let turn = q.to_f64().unwrap();
                        let point = source
                            .point_in_band(
                                usize::from(knot.source_boundary),
                                knot.source_start,
                                0.,
                                z,
                            )
                            .unwrap();
                        let r = source.bands()[usize::from(knot.source_boundary)].radius_mm[0];
                        let theta = std::f64::consts::TAU * turn;
                        assert!(
                            (point.point_mm[0] - r * theta.portable_cos()).abs()
                                + (point.point_mm[1] - r * theta.portable_sin()).abs()
                                < point.error_mm + 5e-13
                        );
                    }
                    let mut area = 0.;
                    for (i, cell) in cap.cells().iter().enumerate() {
                        assert_eq!(cell.events, [i, i + 1]);
                        let b = cap.cell_bounds(i, options(), None).unwrap();
                        for k in 0..2 {
                            let turn = cap.angles[cell.events[k]]
                                .decode(&mut w)
                                .unwrap()
                                .to_f64()
                                .unwrap();
                            for (source, bounds) in [
                                (body.outer(), b.outer_radius_mm[k]),
                                (body.inner(), b.inner_radius_mm[k]),
                            ] {
                                let expected = nominal(source, turn, z);
                                assert!(
                                    bounds[0] - 1e-12 <= expected && expected <= bounds[1] + 1e-12
                                );
                            }
                            assert!(
                                b.outer_radius_mm[k][0] - b.inner_radius_mm[k][1]
                                    >= body.minimum_wall_mm()
                            );
                        }
                        let a = cap.angles[i].decode(&mut w).unwrap().to_f64().unwrap();
                        let c = cap.angles[i + 1].decode(&mut w).unwrap().to_f64().unwrap();
                        let mean = |v: [[f64; 2]; 2]| {
                            let [x, y] = v.map(|b| b[0].midpoint(b[1]));
                            (x * x + x * y + y * y) / 3.
                        };
                        area += std::f64::consts::PI
                            * (c - a)
                            * (mean(b.outer_radius_mm) - mean(b.inner_radius_mm));
                    }
                    let v = body
                        .volume(
                            QueryOptions {
                                volume_error_mm3: 1e-7,
                                maximum_cells: 10,
                            },
                            None,
                        )
                        .unwrap();
                    assert!((area - v.material_mm3[0] / body.full_length_mm()).abs() < 1e-10);
                }
            }
        }
    }
}
#[test]
fn real_positive_distinct_angle_events_are_never_proximity_merged() {
    let body = fixture(false, 1, true);
    let cap = body.cap_domains(CapEnd::End, options(), None).unwrap();
    let mut w = work(options(), None).unwrap();
    let mut closest = f64::INFINITY;
    for pair in cap.angles.windows(2) {
        let a = pair[0].decode(&mut w).unwrap();
        let b = pair[1].decode(&mut w).unwrap();
        assert!(b > a);
        closest = closest.min((b - a).to_f64().unwrap());
    }
    assert!(closest > 0. && closest < 0.1);
    assert!(cap.cells().len() >= cap.knots().len() / 2);
}
#[test]
fn exact_graph_work_storage_cancel_and_retry_are_atomic() {
    let body = fixture(false, 1, true);
    let cap = body.cap_domains(CapEnd::End, options(), None).unwrap();
    assert!(matches!(
        body.cap_domains(
            CapEnd::End,
            CapOptions {
                maximum_work: 0,
                ..options()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        body.cap_domains(
            CapEnd::End,
            CapOptions {
                maximum_bytes: SCRATCH_BYTES - 1,
                ..options()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        body.cap_domains(
            CapEnd::End,
            CapOptions {
                maximum_bytes: SCRATCH_BYTES + 1,
                ..options()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert_eq!(
        cap.cell_bounds(usize::MAX, options(), None),
        Err(Error::Invalid)
    );
    assert_eq!(
        cap.cell_bounds(
            0,
            CapOptions {
                maximum_bytes: SCRATCH_BYTES,
                ..options()
            },
            None
        ),
        Err(Error::Limit)
    );
    let calls = Cell::new(0);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    body.cap_domains(CapEnd::End, options(), Some(&count))
        .unwrap();
    let total = calls.get();
    for stop in 1..=total {
        let calls = Cell::new(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() == stop
        };
        assert!(matches!(
            body.cap_domains(CapEnd::End, options(), Some(&cancel)),
            Err(Error::Cancelled)
        ));
    }
    let calls = Cell::new(0);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    cap.cell_bounds(0, options(), Some(&count)).unwrap();
    let total = calls.get();
    for stop in 1..=total {
        let calls = Cell::new(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() == stop
        };
        assert_eq!(
            cap.cell_bounds(0, options(), Some(&cancel)),
            Err(Error::Cancelled)
        );
    }
    assert_eq!(
        body.cap_domains(CapEnd::End, options(), None)
            .unwrap()
            .knots(),
        cap.knots()
    );
}
#[test]
fn exact_scalar_bit_preflight_refuses_before_unbounded_arithmetic() {
    let mut w = work(options(), None).unwrap();
    let large = Q::from_integer(BigInt::one() << 8192);
    assert_eq!(op(&large, &large, Op::Mul, &mut w), Err(Error::Limit));
    let numerator = Q::from_integer(BigInt::one() << 100);
    assert_eq!(modulo(&numerator, &Q::one(), &mut w), Err(Error::Limit));
}

#[test]
fn one_ulp_source_knot_change_keeps_exact_distinct_proof_events() {
    let body = fixture(false, 1, false);
    let mut inner = body.inner().spec();
    inner.crest_width_mm = 0.075_f64.next_up();
    let inner = AuthoredThread::new(
        inner,
        crate::authored_flat_thread::SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-7,
        },
    )
    .unwrap();
    let body = Insert::new(*body.outer(), inner, 0.5, None).unwrap();
    let cap = body.cap_domains(CapEnd::Start, options(), None).unwrap();
    let mut w = work(options(), None).unwrap();
    let mut found = false;
    let mut overlapping_enclosures = false;
    for pair in cap.angles.windows(2) {
        let a = pair[0].decode(&mut w).unwrap();
        let b = pair[1].decode(&mut w).unwrap();
        let width = (&b - &a).to_f64().unwrap();
        if width < 1e-16 {
            assert!(width > 0.);
            let ai = enclosure(&a, &mut w).unwrap();
            let bi = enclosure(&b, &mut w).unwrap();
            overlapping_enclosures |= ai[1] >= bi[0];
            found = true;
        }
    }
    assert!(found);
    assert!(
        overlapping_enclosures,
        "distinct exact events may have overlapping f64 enclosures"
    );
}

#[test]
fn canonical_default_wood_insert_dimensions_keep_both_original_source_laws() {
    // Current canonical default_custom_connection_wood_insert has outerD10/p2/L4
    // with depth1/crest.4/30deg flanks, and M3x.5 internal6H on the same frame.
    let original = fixture(false, 1, false);
    let mut outer = original.outer().spec();
    outer.crest_width_mm = 0.4;
    outer.lower_flank_deg = 30.;
    outer.upper_flank_deg = 30.;
    let mut inner = original.inner().spec();
    inner.left_hand = false;
    let options_source = crate::authored_flat_thread::SourceOptions {
        position_limit_mm: 1e-7,
        normal_angle_limit_rad: 1e-7,
    };
    let body = Insert::new(
        AuthoredThread::new(outer, options_source).unwrap(),
        AuthoredThread::new(inner, options_source).unwrap(),
        0.5,
        None,
    )
    .unwrap();
    for end in [CapEnd::Start, CapEnd::End] {
        let cap = body.cap_domains(end, options(), None).unwrap();
        assert_eq!(cap.owner().outer().spec(), outer);
        assert_eq!(cap.owner().inner().spec(), inner);
        for index in 0..cap.cells().len() {
            let b = cap.cell_bounds(index, options(), None).unwrap();
            assert!(b.outer_radius_mm[0][0] > b.inner_radius_mm[0][1]);
        }
    }
}

#[test]
fn checked_numeric_projection_cannot_publish_infinite_quotients() {
    let mut w = work(options(), None).unwrap();
    let huge = Q::from_integer(BigInt::one() << 2000);
    assert_eq!(enclosure(&huge, &mut w), Err(Error::Unresolved));
    let tiny = Q::new(BigInt::one(), BigInt::one() << 2000);
    let b = enclosure(&tiny, &mut w).unwrap();
    assert!(binary(b[0]) <= tiny && tiny <= binary(b[1]));
}
