use super::*;
// Actual polynomial patch P(u,v)=(u,v,u²v+v²), whose oriented normal is
// (-2uv,-u²-2v,1). Degree elevation retains its exact rational coefficients.
fn field() -> Grid {
    (0_usize..=17)
        .map(|i| {
            (0..=17)
                .map(|j| {
                    let u = I::point(i as f64) / I::point(17.);
                    let v = I::point(j as f64) / I::point(17.);
                    let u2 = I::point((i * i.saturating_sub(1)) as f64) / I::point(272.);
                    [-I::point(2.) * u * v, -u2 - I::point(2.) * v, I::point(1.)]
                })
                .collect()
        })
        .collect()
}
fn work() -> Work<'static> {
    Work {
        remaining: 1_000_000,
        interrupt: None,
    }
}
fn contains(outer: I, inner: I) {
    assert!(
        outer.lo <= inner.lo && inner.hi <= outer.hi,
        "{outer:?} excludes {inner:?}"
    );
}
#[test]
fn actual_quadratic_patch_mixed_normal_and_global_chart_scalings_are_enclosed() {
    let source = field();
    for domain in [[[0., 1.]; 2], [[2., 5.], [-4., -2.]]] {
        let mut w = work();
        let jet = NormalJet::new(&source, domain, &mut w).unwrap();
        let width = domain.map(|d| I::point(d[1]) - I::point(d[0]));
        contains(jet.hessian[0][1], -I::point(2.) / width[0].square());
        contains(jet.hessian[1][0], -I::point(2.) / (width[0] * width[1]));
        for u in [0.125, 0.5, 0.875] {
            for v in [0.125, 0.5, 0.875] {
                let at = [u, v];
                let point = std::array::from_fn::<_, 2, _>(|k| {
                    domain[k][0] + (domain[k][1] - domain[k][0]) * at[k]
                });
                let bound = jet.enclose(point.map(|x| [x; 2]), &mut w).unwrap();
                let u = I::point(u);
                let v = I::point(v);
                let exact = [
                    -I::point(2.) * u * v,
                    -u.square() - I::point(2.) * v,
                    I::point(1.),
                ];
                let whole = jet.enclose(domain, &mut w).unwrap();
                for k in 0..3 {
                    contains(bound[k], exact[k]);
                    contains(whole[k], exact[k]);
                }
            }
        }
    }
}

#[test]
fn bounded_numeric_refusals_cancellation_and_retries_do_not_change_source() {
    let source = field();
    let before = source.clone();
    let domain = [[0., 1.]; 2];
    let cancel = || true;
    let stop: Option<&dyn Fn() -> bool> = Some(&cancel);
    for (remaining, interrupt, expected) in
        [(0, None, Error::Limit), (1_000_000, stop, Error::Cancelled)]
    {
        let mut w = Work {
            remaining,
            interrupt,
        };
        assert_eq!(
            NormalJet::new(&source, domain, &mut w).err(),
            Some(expected)
        );
    }
    let mut w = work();
    assert!(NormalJet::new(&source, [[0., 0.], [0., 1.]], &mut w).is_err());
    assert!(NormalJet::new(&source, [[0., 1e-300], [0., 1.]], &mut w).is_err());
    let mut broken = source.clone();
    broken[8][8][0] = I::point(f64::NAN);
    assert!(NormalJet::new(&broken, domain, &mut w).is_err());
    let jet = NormalJet::new(&source, domain, &mut w).unwrap();
    let range = [[0.25, 0.75]; 2];
    let bound = jet.enclose(range, &mut w).unwrap();
    assert_eq!(
        jet.enclose(
            range,
            &mut Work {
                remaining: 0,
                interrupt: None
            }
        )
        .err(),
        Some(Error::Limit)
    );
    assert_eq!(
        jet.enclose(
            range,
            &mut Work {
                remaining: 1_000_000,
                interrupt: Some(&|| true)
            }
        )
        .err(),
        Some(Error::Cancelled)
    );
    assert!(jet.enclose([[2., 3.], [0., 1.]], &mut w).is_err());
    assert_eq!(
        jet.enclose([[f64::NAN; 2], [0., 1.]], &mut w).err(),
        Some(Error::Invalid)
    );
    assert_eq!(
        jet.enclose(range, &mut w).unwrap().map(|v| [v.lo, v.hi]),
        bound.map(|v| [v.lo, v.hi])
    );
    assert_eq!(
        source
            .iter()
            .flatten()
            .flatten()
            .map(|v| [v.lo, v.hi])
            .collect::<Vec<_>>(),
        before
            .iter()
            .flatten()
            .flatten()
            .map(|v| [v.lo, v.hi])
            .collect::<Vec<_>>()
    );
}
