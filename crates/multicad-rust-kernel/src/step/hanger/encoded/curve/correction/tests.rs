use super::super::restrict::tests::jet;
use super::*;
use crate::authored_hanger_bolt::{BoundaryOptions, Hanger};
use crate::step::hanger::encoded::cap::tests::number;
use num_traits::{One, Signed};
type Q = BigRational;
#[test]
fn actual_full_root_basis_endpoint_correction_is_exact_on_finite_interval() {
    for (left, nb, period, end, piece) in [
        (false, false, -1, End::Start, 3),
        (true, true, 25, End::End, 0),
    ] {
        let s = Hanger::new(crate::authored_hanger_bolt::tests::fixture(left, nb), None).unwrap();
        let b = s.boundary(BoundaryOptions::default(), None).unwrap();
        let f = EncodedField::prepare(&b, CurveOptions::default(), None).unwrap();
        let index = f
            .cells()
            .iter()
            .position(|c| {
                c.source_key.total_period == period
                    && matches!(
                        c.source_key.face,
                        FaceKey::MetricThread {
                            band: ThreadBandRole::Root,
                            ..
                        }
                    )
            })
            .unwrap();
        let curve = f
            .cap_curve_piece(index, end, piece, DomainOptions::default(), None)
            .unwrap()
            .unwrap();
        let mut w = Work::new(16_000_000, None).unwrap();
        let mut m = Memory::new(0, 64 * 1024 * 1024).unwrap();
        let a = binary(
            (curve.parameter_bounds[0][0] + curve.parameter_bounds[0][1]) / 2.,
            &mut w,
        )
        .unwrap();
        let b = binary(
            (curve.parameter_bounds[1][0] + curve.parameter_bounds[1][1]) / 2.,
            &mut w,
        )
        .unwrap();
        let h = curve
            .controls()
            .iter()
            .map(|p| p.each_ref().map(|x| N(number(x))))
            .collect::<Vec<_>>();
        let small = N(Q::new(1.into(), num_bigint::BigInt::one() << 80));
        let delta = [
            [small.clone(), N::integer(0), N::integer(0)],
            [N::integer(0), small, N::integer(0)],
        ];
        let corrected = correct(&h, &delta, &a, &b, &mut m, &mut w).unwrap();
        assert_eq!(corrected.controls.len(), 24);
        let old = h.iter().map(|p| p.clone().map(|x| x.0)).collect::<Vec<_>>();
        let new = corrected
            .controls
            .iter()
            .map(|p| p.clone().map(|x| x.0))
            .collect::<Vec<_>>();
        for k in 0..=8 {
            let t = Q::new(k.into(), 8.into());
            let u = &a.0 + (&b.0 - &a.0) * &t;
            let (point, d) = jet(old.clone(), &u);
            let (actual, ad) = jet(new.clone(), &u);
            let correction = std::array::from_fn::<_, 3, _>(|axis| {
                &delta[0][axis].0 * (Q::one() - &t) + &delta[1][axis].0 * &t
            });
            for axis in 0..3 {
                assert_eq!(&actual[axis] - &point[axis], correction[axis]);
                assert_eq!(
                    &ad[axis] - &d[axis],
                    (&delta[1][axis].0 - &delta[0][axis].0) / (&b.0 - &a.0)
                );
            }
            assert!(
                (0..3).map(|i| (&actual[i] - &point[i]).abs()).sum::<Q>()
                    <= Q::from_float(corrected.position_error.hi).unwrap()
            );
            assert!(
                (0..3).map(|i| (&ad[i] - &d[i]).abs()).sum::<Q>()
                    <= Q::from_float(corrected.derivative_error.hi).unwrap()
            );
        }
        eprintln!(
            "full basis correction left={left} nb={nb} work={} max_weight_bits={}",
            w.used,
            corrected
                .controls
                .iter()
                .map(|p| p[3].0.denom().bits())
                .max()
                .unwrap()
        );
    }
}
