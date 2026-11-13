use super::*;
use std::cell::Cell;
fn options() -> CapOptions {
    CapOptions {
        maximum_work: 16_000_000,
        maximum_bytes: 1024 * 1024,
    }
}
fn body() -> Insert {
    super::super::super::tests::fixture(false, 1, false)
}
fn proposal(body: &Insert, w: &mut Work<'_>) -> ([[Q; 2]; 3], [[f64; 3]; 3]) {
    let uv = [
        [Q::new(5.into(), 64.into()), Q::one()],
        [Q::new(5.into(), 64.into()), Q::zero()],
        [Q::new(6.into(), 64.into()), Q::one()],
    ];
    let p = uv
        .each_ref()
        .map(|v| super::super::affine_probe::literal(body, &Q::zero(), &v[0], &v[1], w));
    (uv, p)
}
#[test]
fn directed_numeric_primitives_refuse_bad_proposals_overflow_and_degenerate_domains() {
    let mut w = work(options(), None).unwrap();
    for bad in [f64::NAN, f64::INFINITY] {
        assert_eq!(binary_work(bad, &mut w), Err(Error::Unresolved));
    }
    for x in [
        0.,
        -0.,
        f64::from_bits(1),
        -f64::from_bits(1),
        0.7,
        -12345.,
        f64::MAX,
    ] {
        assert_eq!(binary_work(x, &mut w).unwrap(), binary(x));
    }
    assert!(matches!(
        interval(&binary(f64::MAX), &mut w),
        Err(Error::Unresolved)
    ));
    assert!(interval_proposal(&binary(0.5), 0.5, &mut w).is_ok());
    assert!(matches!(
        interval_proposal(&binary(0.5), 0., &mut w),
        Err(Error::Unresolved)
    ));
    let c = [[I::point(f64::MAX), I::point(0.)]; 2];
    assert!(matches!(
        field(c, I::point(0.5), I::point(0.5), I::point(0.), &mut w),
        Err(Error::Unresolved)
    ));
    let bounds = [[Q::zero(), Q::one()], [Q::zero(), Q::one()]];
    let c = [[I::point(1.), I::point(0.)], [I::point(2.), I::point(0.)]];
    let a = [[I::point(f64::MAX); 3]; 3];
    assert_eq!(
        residual(c, &bounds, I::point(0.), a, &mut w),
        Err(Error::Unresolved)
    );
    let uv = std::array::from_fn(|_| [Q::zero(), Q::zero()]);
    assert!(matches!(
        affine(&uv, [[0.; 3]; 3], &mut w),
        Err(Error::Unresolved)
    ));
}
#[test]
fn actual_cap_proof_rejects_changed_literal_geometry_and_bounds_its_proof_resources() {
    let body = body();
    let cap = body.cap_domains(CapEnd::Start, options(), None).unwrap();
    let mut w = work(options(), None).unwrap();
    let proof = CapProof::new(&cap, &mut w).unwrap();
    let (uv, p) = proposal(&body, &mut w);
    assert!(proof.certify(&uv, p, 0.1, &mut w).unwrap() <= 0.1);
    assert_eq!(
        proof.certify_bounded(&uv, p, 0.1, 0, 12, &mut w),
        Err(Error::Limit)
    );
    let mut wrong = p;
    wrong[0][0] += 1.;
    assert_eq!(
        proof.certify_bounded(&uv, wrong, 0.1, 4096, 0, &mut w),
        Err(Error::Unresolved)
    );
    let mut foreign = uv.clone();
    foreign[0][0] = Q::from_integer((-1).into());
    assert_eq!(
        proof.certify(&foreign, p, 0.1, &mut w),
        Err(Error::Unresolved)
    );
    let mut lifted = uv.clone();
    for v in &mut lifted {
        v[0] += Q::one();
    }
    assert!(proof.certify(&lifted, p, 0.1, &mut w).unwrap() <= 0.1);
}
#[test]
fn cap_coefficient_construction_cancels_at_every_actual_checkpoint_without_source_mutation() {
    let body = body();
    let before = body;
    let cap = body.cap_domains(CapEnd::Start, options(), None).unwrap();
    let count = Cell::new(0);
    let callback = || {
        count.set(count.get() + 1);
        false
    };
    CapProof::new(&cap, &mut work(options(), Some(&callback)).unwrap()).unwrap();
    let total = count.get();
    for stop in 1..=total {
        let count = Cell::new(0);
        let cancel = || {
            count.set(count.get() + 1);
            count.get() == stop
        };
        let result = work(options(), Some(&cancel)).and_then(|mut w| CapProof::new(&cap, &mut w));
        assert!(matches!(result, Err(Error::Cancelled)));
    }
    assert_eq!(body, before);
}
