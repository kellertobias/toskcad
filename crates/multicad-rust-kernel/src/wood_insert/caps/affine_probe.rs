//! Test-only actual cap correspondence diagnostic; no runtime admission.
use super::*;
fn add(a: &Q, b: &Q, w: &mut Work<'_>) -> Q {
    op(a, b, Op::Add, w).unwrap()
}
fn sub(a: &Q, b: &Q, w: &mut Work<'_>) -> Q {
    op(a, b, Op::Sub, w).unwrap()
}
fn mul(a: &Q, b: &Q, w: &mut Work<'_>) -> Q {
    op(a, b, Op::Mul, w).unwrap()
}
fn div(a: &Q, b: &Q, w: &mut Work<'_>) -> Q {
    op(a, b, Op::Div, w).unwrap()
}
fn ij(q: &Q, w: &mut Work<'_>) -> I {
    let [a, b] = enclosure(q, w).unwrap();
    I::new(a, b)
}
fn coeff(source: &AuthoredThread, z: &Q, band: u8, k: i64, w: &mut Work<'_>) -> [Q; 2] {
    let a = radius(source, z, &Q::zero(), band, k, w).unwrap();
    let b = radius(source, z, &Q::one(), band, k, w).unwrap();
    let slope = sub(&b, &a, w);
    [a, slope]
}
fn literal(body: &Insert, z: &Q, u: &Q, t: &Q, w: &mut Work<'_>) -> [f64; 3] {
    let r = |source: &AuthoredThread, w: &mut Work<'_>| {
        let (b, k) = band_at(source, z, u, w).unwrap();
        radius(source, z, u, b, k, w).unwrap()
    };
    let inner = r(body.inner(), w);
    let outer = r(body.outer(), w);
    let radius = add(&inner, &mul(&sub(&outer, &inner, w), t, w), w);
    let tau = I::new(
        std::f64::consts::PI.next_down(),
        std::f64::consts::PI.next_up(),
    ) * I::point(2.);
    let angle = ij(u, w) * tau;
    let r = ij(&radius, w);
    let p = [r * angle.cos(), r * angle.sin(), ij(z, w)];
    p.map(|p| p.midpoint())
}
fn affine(uv: &[[Q; 2]; 3], p: [[f64; 3]; 3], w: &mut Work<'_>) -> [[I; 3]; 3] {
    let du = [sub(&uv[1][0], &uv[0][0], w), sub(&uv[2][0], &uv[0][0], w)];
    let dt = [sub(&uv[1][1], &uv[0][1], w), sub(&uv[2][1], &uv[0][1], w)];
    let determinant = sub(&mul(&du[0], &dt[1], w), &mul(&du[1], &dt[0], w), w);
    assert!(determinant > Q::zero());
    std::array::from_fn(|axis| {
        let dp = [
            sub(&binary(p[1][axis]), &binary(p[0][axis]), w),
            sub(&binary(p[2][axis]), &binary(p[0][axis]), w),
        ];
        let a = div(
            &sub(&mul(&dp[0], &dt[1], w), &mul(&dp[1], &dt[0], w), w),
            &determinant,
            w,
        );
        let b = div(
            &sub(&mul(&du[0], &dp[1], w), &mul(&du[1], &dp[0], w), w),
            &determinant,
            w,
        );
        let c = sub(
            &sub(&binary(p[0][axis]), &mul(&a, &uv[0][0], w), w),
            &mul(&b, &uv[0][1], w),
            w,
        );
        [ij(&a, w), ij(&b, w), ij(&c, w)]
    })
}
fn cross(a: &[Q; 2], b: &[Q; 2], p: &[Q; 2], w: &mut Work<'_>) -> Q {
    sub(
        &mul(&sub(&b[0], &a[0], w), &sub(&p[1], &a[1], w), w),
        &mul(&sub(&b[1], &a[1], w), &sub(&p[0], &a[0], w), w),
        w,
    )
}
fn outside(uv: &[[Q; 2]; 3], bounds: &[[Q; 2]; 2], w: &mut Work<'_>) -> bool {
    let corners = [
        [bounds[0][0].clone(), bounds[1][0].clone()],
        [bounds[0][1].clone(), bounds[1][0].clone()],
        [bounds[0][0].clone(), bounds[1][1].clone()],
        [bounds[0][1].clone(), bounds[1][1].clone()],
    ];
    for e in 0..3 {
        if corners
            .iter()
            .all(|p| cross(&uv[e], &uv[(e + 1) % 3], p, w) < Q::zero())
        {
            return true;
        }
    }
    false
}
// Actual C2 piece: r(w,t)=ri(w)+t*(ro(w)−ri(w)), each ri/ro affine.
// Columns per Cartesian component: value,dw,dt,dww,dwt (dtt exactlyzero).
fn field(c: [[I; 2]; 2], u: I, t: I, z: I) -> [[I; 5]; 3] {
    let tau = I::new(
        std::f64::consts::PI.next_down(),
        std::f64::consts::PI.next_up(),
    ) * I::point(2.);
    let angle = u * tau;
    let co = angle.cos();
    let si = angle.sin();
    let inner = c[0][0] + c[0][1] * u;
    let gap = c[1][0] - c[0][0] + (c[1][1] - c[0][1]) * u;
    let r = inner + t * gap;
    let rw = c[0][1] + t * (c[1][1] - c[0][1]);
    let slope = c[1][1] - c[0][1];
    [
        [
            r * co,
            rw * co - tau * r * si,
            gap * co,
            -I::point(2.) * tau * rw * si - tau.square() * r * co,
            slope * co - tau * gap * si,
        ],
        [
            r * si,
            rw * si + tau * r * co,
            gap * si,
            I::point(2.) * tau * rw * co - tau.square() * r * si,
            slope * si + tau * gap * co,
        ],
        [z, I::point(0.), I::point(0.), I::point(0.), I::point(0.)],
    ]
}
fn residual(c: [[I; 2]; 2], bounds: &[[Q; 2]; 2], z: I, a: [[I; 3]; 3], w: &mut Work<'_>) -> f64 {
    let center: [Q; 2] = std::array::from_fn(|i| {
        div(
            &add(&bounds[i][0], &bounds[i][1], w),
            &Q::from_integer(2.into()),
            w,
        )
    });
    let center_i = center.each_ref().map(|q| ij(q, w));
    let range: [I; 2] = std::array::from_fn(|i| {
        let lo = ij(&bounds[i][0], w);
        let hi = ij(&bounds[i][1], w);
        I::new(lo.lo, hi.hi)
    });
    let delta: [I; 2] = std::array::from_fn(|i| {
        let lo = ij(&sub(&bounds[i][0], &center[i], w), w);
        let hi = ij(&sub(&bounds[i][1], &center[i], w), w);
        I::new(lo.lo, hi.hi)
    });
    let f0 = field(c, center_i[0], center_i[1], z);
    let whole = field(c, range[0], range[1], z);
    let mut squares = I::point(0.);
    for axis in 0..3 {
        let l = a[axis][0] * center_i[0] + a[axis][1] * center_i[1] + a[axis][2];
        let r0 = f0[axis][0] - l;
        let rw = f0[axis][1] - a[axis][0];
        let rt = f0[axis][2] - a[axis][1];
        let bound = I::point(r0.magnitude())
            + I::point(rw.magnitude()) * I::point(delta[0].magnitude())
            + I::point(rt.magnitude()) * I::point(delta[1].magnitude())
            + I::point(0.5)
                * I::point(whole[axis][3].magnitude())
                * I::point(delta[0].magnitude()).square()
            + I::point(whole[axis][4].magnitude())
                * I::point(delta[0].magnitude())
                * I::point(delta[1].magnitude());
        squares = squares + bound.square();
    }
    squares.sqrt().hi
}
fn proof(
    cap: &CapDomains,
    uv: &[[Q; 2]; 3],
    p: [[f64; 3]; 3],
    subdivisions: usize,
    w: &mut Work<'_>,
) -> (f64, usize) {
    let affine = affine(uv, p, w);
    let z = station(cap.owner(), cap.end());
    let z_i = ij(&z, w);
    let mut maximum = 0_f64;
    let mut boxes = 0;
    let min = uv.iter().map(|p| p[0].clone()).min().unwrap();
    let max = uv.iter().map(|p| p[0].clone()).max().unwrap();
    for cell in cap.cells() {
        let lo = cap.angles[cell.events[0]]
            .decode(w)
            .unwrap()
            .max(min.clone());
        let hi = cap.angles[cell.events[1]]
            .decode(w)
            .unwrap()
            .min(max.clone());
        if lo >= hi {
            continue;
        }
        let inner = coeff(
            cap.owner().inner(),
            &z,
            cell.inner_band,
            cell.inner_period,
            w,
        )
        .map(|q| ij(&q, w));
        let outer = coeff(
            cap.owner().outer(),
            &z,
            cell.outer_band,
            cell.outer_period,
            w,
        )
        .map(|q| ij(&q, w));
        for i in 0..subdivisions {
            for j in 0..subdivisions {
                let endpoint = |k: usize, w: &mut Work<'_>| {
                    add(
                        &lo,
                        &mul(&sub(&hi, &lo, w), &Q::new(k.into(), subdivisions.into()), w),
                        w,
                    )
                };
                let bounds = [
                    [endpoint(i, w), endpoint(i + 1, w)],
                    [
                        Q::new(j.into(), subdivisions.into()),
                        Q::new((j + 1).into(), subdivisions.into()),
                    ],
                ];
                if outside(uv, &bounds, w) {
                    continue;
                }
                boxes += 1;
                maximum = maximum.max(residual([inner, outer], &bounds, z_i, affine, w));
            }
        }
    }
    (maximum, boxes)
}
fn adaptive(
    cap: &CapDomains,
    uv: &[[Q; 2]; 3],
    p: [[f64; 3]; 3],
    limit: f64,
    w: &mut Work<'_>,
) -> Result<(f64, usize, usize), Error> {
    let a = affine(uv, p, w);
    let z = station(cap.owner(), cap.end());
    let z_i = ij(&z, w);
    let mut maximum = 0_f64;
    let mut checked = 0;
    let mut accepted = 0;
    let min = uv.iter().map(|p| p[0].clone()).min().unwrap();
    let max = uv.iter().map(|p| p[0].clone()).max().unwrap();
    for cell in cap.cells() {
        let lo = cap.angles[cell.events[0]].decode(w)?.max(min.clone());
        let hi = cap.angles[cell.events[1]].decode(w)?.min(max.clone());
        if lo >= hi {
            continue;
        }
        let inner = coeff(
            cap.owner().inner(),
            &z,
            cell.inner_band,
            cell.inner_period,
            w,
        )
        .map(|q| ij(&q, w));
        let outer = coeff(
            cap.owner().outer(),
            &z,
            cell.outer_band,
            cell.outer_period,
            w,
        )
        .map(|q| ij(&q, w));
        let mut stack = vec![([[lo, hi], [Q::zero(), Q::one()]], 0)];
        while let Some((bounds, depth)) = stack.pop() {
            if outside(uv, &bounds, w) {
                continue;
            }
            checked += 1;
            if checked > 4096 {
                return Err(Error::Limit);
            }
            let e = residual([inner, outer], &bounds, z_i, a, w);
            if e <= limit {
                maximum = maximum.max(e);
                accepted += 1;
                continue;
            }
            if depth == 12 {
                return Err(Error::Unresolved);
            }
            let mid: [Q; 2] = std::array::from_fn(|i| {
                div(
                    &add(&bounds[i][0], &bounds[i][1], w),
                    &Q::from_integer(2.into()),
                    w,
                )
            });
            for i in 0..2 {
                for j in 0..2 {
                    let half = |axis: usize, upper: usize| {
                        if upper == 0 {
                            [bounds[axis][0].clone(), mid[axis].clone()]
                        } else {
                            [mid[axis].clone(), bounds[axis][1].clone()]
                        }
                    };
                    stack.push(([half(0, i), half(1, j)], depth + 1));
                }
            }
        }
    }
    Ok((maximum, checked, accepted))
}
#[test]
fn actual_piecewise_cap_source_to_literal_affine_chord_is_provable_without_new_vertices() {
    let body = super::super::tests::fixture(false, 1, false);
    let options = CapOptions {
        maximum_work: 16_000_000,
        maximum_bytes: 1024 * 1024,
    };
    let cap = body.cap_domains(CapEnd::Start, options, None).unwrap();
    let mut work = work(options, None).unwrap();
    let uv = [
        [Q::new(5.into(), 64.into()), Q::one()],
        [Q::new(5.into(), 64.into()), Q::zero()],
        [Q::new(6.into(), 64.into()), Q::one()],
    ];
    let z = Q::zero();
    let p = uv
        .each_ref()
        .map(|uv| literal(&body, &z, &uv[0], &uv[1], &mut work));
    let mut results = Vec::new();
    for grid in [1, 2, 4, 8, 16, 32] {
        let (bound, boxes) = proof(&cap, &uv, p, grid, &mut work);
        results.push((grid, bound, boxes));
    }
    eprintln!(
        "actual piecewise cap chord proof {results:?} {}work",
        work.used
    );
    assert!(results.last().unwrap().1 < 0.1);
    let mut wrong = p;
    wrong[0][0] += 1.;
    assert!(proof(&cap, &uv, wrong, 8, &mut work).0 > 0.1);
}

#[test]
fn adaptive_exact_proof_partition_keeps_geometry_and_needs_fewer_boxes() {
    let body = super::super::tests::fixture(false, 1, false);
    let options = CapOptions {
        maximum_work: 16_000_000,
        maximum_bytes: 1024 * 1024,
    };
    let cap = body.cap_domains(CapEnd::Start, options, None).unwrap();
    let mut work = work(options, None).unwrap();
    let uv = [
        [Q::new(5.into(), 64.into()), Q::one()],
        [Q::new(5.into(), 64.into()), Q::zero()],
        [Q::new(6.into(), 64.into()), Q::one()],
    ];
    let p = uv
        .each_ref()
        .map(|uv| literal(&body, &Q::zero(), &uv[0], &uv[1], &mut work));
    let result = adaptive(&cap, &uv, p, 0.1, &mut work).unwrap();
    eprintln!("adaptive actualcap {result:?} {}weightedQwork", work.used);
    assert!(result.0 <= 0.1);
    assert!(result.1 < 726);
    let mut wrong = p;
    wrong[0][0] += 1.;
    assert_eq!(
        adaptive(&cap, &uv, wrong, 0.1, &mut work),
        Err(Error::Unresolved)
    );
}
