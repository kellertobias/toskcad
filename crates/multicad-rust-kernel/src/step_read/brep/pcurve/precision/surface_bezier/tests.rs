use super::*;
const BASE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../multicad-app/tests/fixtures/exchange/corpus/occt-cone-frustum.step"
));
fn input(weighted: bool) -> String {
    let mut extra = String::new();
    let mut rows = Vec::new();
    let mut weights = Vec::new();
    for u in 0i64..4 {
        let mut row = Vec::new();
        let mut wr = Vec::new();
        for v in 0i64..4 {
            let id = 900010 + u * 4 + v;
            row.push(format!("#{id}"));
            wr.push(format!("{}.", 1 + (u + v) % 3));
            extra += &format!(
                "#{id}=CARTESIAN_POINT('',({u}.,{v}.,{}.));\n",
                (u - v) * (u - v)
            );
        }
        rows.push(format!("({})", row.join(",")));
        weights.push(format!("({})", wr.join(",")));
    }
    let grid = rows.join(",");
    extra += &if weighted {
        format!(
            "#900001=(BOUNDED_SURFACE() B_SPLINE_SURFACE(2,2,({grid}),.UNSPECIFIED.,.F.,.F.,.F.) B_SPLINE_SURFACE_WITH_KNOTS((3,1,3),(3,1,3),(0.,0.5,1.),(0.,0.5,1.),.UNSPECIFIED.) GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_SURFACE(({})) REPRESENTATION_ITEM('') SURFACE());\n",
            weights.join(",")
        )
    } else {
        format!(
            "#900001=B_SPLINE_SURFACE_WITH_KNOTS('',2,2,({grid}),.UNSPECIFIED.,.F.,.F.,.F.,(3,1,3),(3,1,3),(0.,0.5,1.),(0.,0.5,1.),.UNSPECIFIED.);\n"
        )
    };
    let at = BASE.rfind("ENDSEC;").unwrap();
    format!("{}{}{}", &BASE[..at], extra, &BASE[at..])
}
fn q(n: i64, d: i64) -> BigRational {
    BigRational::new(n.into(), d.into())
}
// Independent exact Cox-de Boor oracle, not the extraction recursion.
fn basis(i: usize, p: usize, t: &BigRational, k: &[BigRational]) -> BigRational {
    if p == 0 {
        return if k[i] <= *t && *t < k[i + 1] {
            q(1, 1)
        } else {
            q(0, 1)
        };
    }
    let a = if k[i + p] > k[i] {
        (t - &k[i]) / (&k[i + p] - &k[i]) * basis(i, p - 1, t, k)
    } else {
        q(0, 1)
    };
    let b = if k[i + p + 1] > k[i + 1] {
        (&k[i + p + 1] - t) / (&k[i + p + 1] - &k[i + 1]) * basis(i + 1, p - 1, t, k)
    } else {
        q(0, 1)
    };
    a + b
}
fn oracle(u: &BigRational, v: &BigRational, weighted: bool) -> [BigRational; 3] {
    let knots = vec![
        q(0, 1),
        q(0, 1),
        q(0, 1),
        q(1, 2),
        q(1, 1),
        q(1, 1),
        q(1, 1),
    ];
    let mut h: [BigRational; 4] = std::array::from_fn(|_| q(0, 1));
    for i in 0..4 {
        for j in 0..4 {
            let at = |n, t: &BigRational| {
                if *t == q(1, 1) {
                    if n == 3 { q(1, 1) } else { q(0, 1) }
                } else {
                    basis(n, 2, t, &knots)
                }
            };
            let weight = if weighted { 1 + (i + j) % 3 } else { 1 };
            let b = at(i, u) * at(j, v) * q(weight as i64, 1);
            let values = [i as i64, j as i64, (i as i64 - j as i64).pow(2), 1];
            for axis in 0..4 {
                h[axis] += &b * q(values[axis], 1);
            }
        }
    }
    std::array::from_fn(|i| &h[i] / &h[3])
}
fn evaluate(p: &Patch, u: f64, v: f64) -> V {
    let mut rows = Vec::new();
    for row in p.controls().chunks_exact(p.shape()[1]) {
        let mut net = row.to_vec();
        for n in (1..net.len()).rev() {
            for i in 0..n {
                net[i] = std::array::from_fn(|j| {
                    I::point(1. - v) * net[i][j] + I::point(v) * net[i + 1][j]
                });
            }
        }
        rows.push(net[0]);
    }
    for n in (1..rows.len()).rev() {
        for i in 0..n {
            rows[i] = std::array::from_fn(|j| {
                I::point(1. - u) * rows[i][j] + I::point(u) * rows[i + 1][j]
            });
        }
    }
    std::array::from_fn(|i| rows[0][i] / rows[0][3])
}
#[test]
fn both_original_tensor_axes_match_an_independent_exact_oracle() {
    for weighted in [false, true] {
        let text = input(weighted);
        let source = Source::parse(&text, 1e-7).unwrap();
        let Surface::Nurbs { surface, .. } = source.surface(900001).unwrap() else {
            panic!()
        };
        let mut budget = Budget::new(512_000_000, None);
        let data = original_surface(&source, 900001, &surface, &mut budget).unwrap();
        assert!(data.matches(&source, 900001));
        assert!(!data.matches(&source, 900002));
        let other = Source::parse(&text, 1e-7).unwrap();
        assert!(!data.matches(&other, 900001));
        assert_eq!(data.degrees(), [2, 2]);
        assert_eq!(data.counts(), [2, 2]);
        assert_eq!(data.patches().len(), 4);
        assert!(data.retained_bytes() <= data.peak_bytes() && data.peak_bytes() <= MAX_BYTES);
        for (i, p) in data.patches().iter().enumerate() {
            let [pu, pv] = [i / 2, i % 2];
            assert_eq!(p.shape(), [3, 3]);
            assert_eq!(p.knot_indices(), [[pu, pu + 1], [pv, pv + 1]]);
            for a in 0..=4 {
                for b in 0..=4 {
                    let value = evaluate(p, a as f64 / 4., b as f64 / 4.);
                    let exact = oracle(
                        &q((pu * 4 + a) as i64, 8),
                        &q((pv * 4 + b) as i64, 8),
                        weighted,
                    );
                    for axis in 0..3 {
                        assert!(BigRational::from_float(value[axis].lo).unwrap() <= exact[axis]);
                        assert!(BigRational::from_float(value[axis].hi).unwrap() >= exact[axis]);
                    }
                }
            }
            if !weighted {
                assert!(p.controls().iter().all(|h| h[3].lo == 1. && h[3].hi == 1.));
            }
        }
    }
}
#[test]
fn actual_constructor_has_exact_work_and_late_cancellation_frontiers() {
    let text = input(true);
    let source = Source::parse(&text, 1e-7).unwrap();
    let Surface::Nurbs { surface, .. } = source.surface(900001).unwrap() else {
        panic!()
    };
    let mut b = Budget::new(512_000_000, None);
    original_surface(&source, 900001, &surface, &mut b).unwrap();
    let used = b.used();
    assert!(matches!(
        original_surface(&source, 900001, &surface, &mut Budget::new(used - 1, None)),
        Err(ReadError::Limit)
    ));
    original_surface(&source, 900001, &surface, &mut Budget::new(used, None)).unwrap();
    let calls = std::cell::Cell::new(0);
    let observe = || {
        calls.set(calls.get() + 1);
        false
    };
    original_surface(
        &source,
        900001,
        &surface,
        &mut Budget::new(512_000_000, Some(&observe)),
    )
    .unwrap();
    let last = calls.get();
    calls.set(0);
    let cancel = || {
        calls.set(calls.get() + 1);
        calls.get() == last
    };
    let mut cancelled = Budget::new(512_000_000, Some(&cancel));
    assert!(matches!(
        original_surface(&source, 900001, &surface, &mut cancelled),
        Err(ReadError::Limit)
    ));
    assert!(cancelled.cancelled.get());
    original_surface(
        &source,
        900001,
        &surface,
        &mut Budget::new(512_000_000, None),
    )
    .unwrap();
}
