use super::*;
use crate::authored_hanger_bolt::{BoundaryOptions, Hanger};
use num_bigint::BigInt;
use num_traits::{One, Zero};
type Q = BigRational;
pub(in crate::step::hanger::encoded) fn number(v: &EncodedScalar) -> Q {
    let sign = match v.signum() {
        -1 => Sign::Minus,
        0 => Sign::NoSign,
        _ => Sign::Plus,
    };
    Q::new(
        BigInt::from_bytes_le(sign, v.numerator_magnitude_le()),
        BigInt::from_bytes_le(Sign::Plus, v.denominator_magnitude_le()),
    )
}
pub(in crate::step::hanger::encoded) fn bezier<const N: usize>(
    mut c: Vec<[Q; N]>,
    t: &Q,
) -> [Q; N] {
    let other = Q::one() - t;
    for n in (1..c.len()).rev() {
        for i in 0..n {
            c[i] = std::array::from_fn(|k| &c[i][k] * &other + &c[i + 1][k] * t);
        }
    }
    c.remove(0)
}
#[test]
fn exact_encoded_cap_pcurves_remain_on_the_same_scalar_owner_planes() {
    for (left, nb) in [(false, false), (true, true)] {
        let s = Hanger::new(crate::authored_hanger_bolt::tests::fixture(left, nb), None).unwrap();
        let b = s.boundary(BoundaryOptions::default(), None).unwrap();
        let f = EncodedField::prepare(&b, CurveOptions::default(), None).unwrap();
        let mut count = 0;
        for i in 0..f.cells().len() {
            for end in [End::Start, End::End] {
                if let Some(trace) = f.cap_trace(i, end, DomainOptions::default(), None).unwrap() {
                    count += 1;
                    assert!(
                        trace.endpoints[0].meridian_parameter[1]
                            < trace.endpoints[1].meridian_parameter[0]
                    );
                    assert!((2..=3).contains(&trace.degree()));
                    let cell = f.cells()[i];
                    let meridian = &f.meridians()[cell.meridian];
                    let expected_plane = match (meridian.component(), end) {
                        (Component::Metric, End::Start) => {
                            -Q::from_float(s.metric().full_length_mm()).unwrap()
                                - Q::from_float(s.spec().plain_length_mm).unwrap()
                        }
                        (Component::Metric, End::End) => {
                            -Q::from_float(s.spec().plain_length_mm).unwrap()
                        }
                        (Component::Wood, End::Start) => Q::zero(),
                        (Component::Wood, End::End) => {
                            Q::from_float(s.wood().full_length_mm()).unwrap()
                        }
                    };
                    assert_eq!(number(&trace.plane_z), expected_plane);
                    let u = (Q::from_float(trace.endpoints[0].meridian_parameter[1]).unwrap()
                        + Q::from_float(trace.endpoints[1].meridian_parameter[0]).unwrap())
                        / Q::from_integer(2.into());
                    let uv = bezier(
                        trace
                            .controls()
                            .iter()
                            .map(|c| c.each_ref().map(number))
                            .collect(),
                        &u,
                    );
                    assert_eq!(&uv[0] / &uv[2], u, "the SAME original meridian parameter");
                    let v = &uv[1] / &uv[2];
                    assert!(v >= Q::zero() && v <= Q::one());
                    let piece = f
                        .angular()
                        .knots()
                        .windows(2)
                        .position(|k| {
                            v >= Q::from_float(k[0]).unwrap() && v <= Q::from_float(k[1]).unwrap()
                        })
                        .unwrap();
                    let a = Q::from_float(f.angular().knots()[piece]).unwrap();
                    let width = Q::from_float(f.angular().knots()[piece + 1]).unwrap() - &a;
                    let local = (&v - &a) / width;
                    let mut columns = Vec::new();
                    for j in 0..=10 {
                        let mut rows = Vec::new();
                        for r in 0..=meridian.degree() {
                            let c = f
                                .control(i, r, piece, j, DomainOptions::default(), None)
                                .unwrap();
                            rows.push(c.homogeneous.each_ref().map(number));
                        }
                        columns.push(bezier(rows, &u));
                    }
                    // Independent exact rational evaluation of ACTUAL tensor
                    // controls, not a repeat of the source affine-z formula.
                    let point = bezier(columns.clone(), &local);
                    assert_eq!(&point[2] / &point[3], expected_plane);
                    // A real encoded-coefficient mutation breaks that identity.
                    columns[0][2] += Q::one();
                    let wrong = bezier(columns, &local);
                    if local < Q::one() {
                        assert_ne!(&wrong[2] / &wrong[3], expected_plane);
                    }
                    let original = &trace.endpoints[0].original_parameter;
                    assert!(number(original) >= Q::zero() && number(original) <= Q::one());
                    assert!(b.curves().iter().any(|c| c.key == trace.source_curve
                        && c.faces.contains(&trace.cap_face)
                        && c.faces.contains(&cell.source_key.face)));
                }
            }
        }
        eprintln!("cap traces left={left} nonbinary={nb} count={count}");
        assert!(count >= 20);
    }
}
#[test]
fn positive_algebraic_trim_sliver_is_refused_without_welding_or_source_change() {
    let mut spec = crate::authored_hanger_bolt::tests::fixture(false, false);
    spec.metric = AuthoredMetricThread::new(crate::authored_metric_thread::Spec {
        major_radius_mm: 1.46,
        pitch_mm: 0.9,
        length_mm: 1.2375,
        crest_width_mm: 0.09,
        starts: 1,
        left_hand: false,
    })
    .unwrap();
    let s = Hanger::new(spec, None).unwrap();
    let b = s.boundary(BoundaryOptions::default(), None).unwrap();
    let field = EncodedField::prepare(&b, CurveOptions::default(), None).unwrap();
    let index = field
        .cells()
        .iter()
        .position(|c| {
            c.source_key.total_period == 1
                && matches!(
                    c.source_key.face,
                    FaceKey::MetricThread {
                        band: ThreadBandRole::Root,
                        ..
                    }
                )
        })
        .unwrap();
    assert!(matches!(
        field.cap_trace(index, End::End, DomainOptions::default(), None),
        Err(Error::Unresolved)
    ));
    assert_eq!(s.spec(), b.source().spec());
}

#[test]
fn shared_trace_ledger_counts_empty_selections_and_retains_identical_source_tokens() {
    let s = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(true, true),
        None,
    )
    .unwrap();
    let b = s.boundary(BoundaryOptions::default(), None).unwrap();
    let f = EncodedField::prepare(&b, CurveOptions::default(), None).unwrap();
    let mut w = Work::new(16_000_000, None).unwrap();
    w.charge(f.work_used()).unwrap();
    let mut m = Memory::new(f.retained_bytes(), 64 * 1024 * 1024).unwrap();
    m.add(SCRATCH + size_of::<CapTrace>()).unwrap();
    let mut empty = 0;
    let mut positive = 0;
    for index in 0..f.cells().len() {
        for end in [End::Start, End::End] {
            let old = f
                .cap_trace(index, end, DomainOptions::default(), None)
                .unwrap();
            let before_work = w.used;
            let before_memory = m.used;
            let shared = f
                .cap_trace_with_ledgers(index, end, &mut m, &mut w)
                .unwrap();
            assert!(w.used > before_work);
            match (old, shared) {
                (None, None) => {
                    empty += 1;
                    assert_eq!(m.used, before_memory);
                }
                (Some(old), Some(shared)) => {
                    positive += 1;
                    assert_eq!(shared.work_used, w.used);
                    assert_eq!(old.source_curve, shared.source_curve);
                    assert_eq!(old.source_key, shared.source_key);
                    assert_eq!(old.cap_face, shared.cap_face);
                    assert_eq!(number(&old.plane_z), number(&shared.plane_z));
                    for (a, b) in old.controls().iter().zip(shared.controls()) {
                        assert_eq!(a.each_ref().map(number), b.each_ref().map(number));
                    }
                    for (a, b) in old.endpoints.iter().zip(&shared.endpoints) {
                        assert_eq!(number(&a.phase), number(&b.phase));
                        assert_eq!(number(&a.original_parameter), number(&b.original_parameter));
                        assert_eq!(number(&a.physical_turn), number(&b.physical_turn));
                        assert_eq!(a.meridian_parameter, b.meridian_parameter);
                    }
                    drop((old, shared));
                    m.used = before_memory;
                }
                _ => panic!("shared ledger changed the actual source selection"),
            }
        }
    }
    assert!(empty > 0 && positive > 0);
    let calls = std::cell::Cell::new(0);
    let cancel = || {
        calls.set(calls.get() + 1);
        calls.get() > 1
    };
    let mut w = Work::new(16_000_000, Some(&cancel)).unwrap();
    assert!(matches!(
        f.cap_trace_with_ledgers(0, End::Start, &mut m, &mut w),
        Err(Error::Cancelled)
    ));
    let mut w = Work::new(1, None).unwrap();
    assert!(matches!(
        f.cap_trace_with_ledgers(0, End::Start, &mut m, &mut w),
        Err(Error::Limit)
    ));
}
