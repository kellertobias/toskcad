use super::*;
use multicad_precision::PortableMath as _;
fn frame() -> RawFrame {
    RawFrame {
        origin: cached([1., -2., 0.5]),
        axes: [
            cached([1., 0., 0.]),
            cached([0., 1., 0.]),
            cached([0., 0., 1.]),
        ],
    }
}
/// Degree-`n` Bernstein controls interpolating `f` at Chebyshev nodes.
/// The interpolant is only a TEST INPUT; the proof below never trusts it.
fn interpolate(n: usize, f: impl Fn(f64) -> f64) -> Vec<f64> {
    let nodes: Vec<f64> = (0..=n)
        .map(|k| {
            0.5 - 0.5
                * ((2 * k + 1) as f64 * std::f64::consts::PI / (2 * n + 2) as f64).portable_cos()
        })
        .collect();
    let choose = |n: usize, k: usize| (1..=k).fold(1., |a, i| a * (n + 1 - i) as f64 / i as f64);
    let mut m: Vec<Vec<f64>> = nodes
        .iter()
        .map(|&t| {
            let mut row: Vec<f64> = (0..=n)
                .map(|i| choose(n, i) * t.powi(i as i32) * (1. - t).powi((n - i) as i32))
                .collect();
            row.push(f(t));
            row
        })
        .collect();
    for c in 0..=n {
        let pivot = (c..=n)
            .max_by(|&a, &b| m[a][c].abs().total_cmp(&m[b][c].abs()))
            .unwrap();
        m.swap(c, pivot);
        for r in 0..=n {
            if r != c {
                let factor = m[r][c] / m[c][c];
                let pivot = m[c].clone();
                for (x, y) in m[r][c..].iter_mut().zip(&pivot[c..]) {
                    *x -= factor * y;
                }
            }
        }
    }
    (0..=n).map(|i| m[i][n + 1] / m[i][i]).collect()
}
fn poly(values: &[f64], b: &mut Budget<'_>) -> P {
    P::new(&values.iter().map(|&v| I::point(v)).collect::<Vec<_>>(), b).unwrap()
}
/// theta(t)=theta0+span*t and v(t)=height*t as degree-7 Bernstein curves on
/// a cylinder (slope 0) or cone, with an interpolated degree-7 target.
fn scene(span: f64, slope: f64, b: &mut Budget<'_>) -> (Support, [P; 2], [P; 3]) {
    let (radius, theta0, height) = (3.165, 0.3, 0.35);
    let linear = |a: f64, c: f64| -> Vec<f64> { (0..=7).map(|i| a + c * i as f64 / 7.).collect() };
    let uv = [poly(&linear(theta0, span), b), poly(&linear(0., height), b)];
    let f = frame();
    let point = |t: f64, axis: usize| {
        let theta = theta0 + span * t;
        let v = height * t;
        let rho = radius + slope * v;
        let local = [rho * theta.portable_cos(), rho * theta.portable_sin(), v];
        f.origin[axis].lo + (0..3).map(|k| f.axes[k][axis].lo * local[k]).sum::<f64>()
    };
    let target = std::array::from_fn(|axis| poly(&interpolate(7, |t| point(t, axis)), b));
    (
        Support {
            frame: f,
            radius: I::point(radius),
            slope: I::point(slope),
        },
        uv,
        target,
    )
}
fn strict(tol: f64) -> impl FnMut(&CellProof, &mut Budget<'_>) -> Result<bool> {
    move |p: &CellProof, _: &mut Budget<'_>| Ok(p.position <= tol && p.derivative < p.speed)
}
#[test]
fn short_cylinder_and_cone_arcs_certify_as_one_whole_cell() {
    for slope in [0., std::f64::consts::FRAC_PI_3.portable_tan()] {
        let mut b = Budget::new(512_000_000, None);
        let (support, uv, target) = scene(0.25, slope, &mut b);
        let mut pieces = 0;
        let mut accept = |p: &CellProof, _: &mut Budget<'_>| {
            pieces += 1;
            Ok(p.position <= 1e-7 && p.derivative < p.speed)
        };
        assert!(certify(&support, &uv, &target, 1e-7, &mut accept, &mut b).unwrap());
        assert_eq!(pieces, 1, "slope {slope}");
    }
}
#[test]
fn a_wide_arc_is_proved_only_through_exact_dyadic_pieces() {
    let mut b = Budget::new(512_000_000, None);
    let (support, uv, target) = scene(1.2, std::f64::consts::FRAC_PI_3.portable_tan(), &mut b);
    // The whole cell exceeds the degree-80 algebra at the required order.
    let [theta, v] = [lift(&uv[0]).unwrap(), lift(&uv[1]).unwrap()];
    let target80 = target.each_ref().map(|p| lift(p).unwrap());
    assert!(matches!(
        proof(&support, [theta, v], target80, 1e-7, &mut b).unwrap(),
        Cell::Undecided
    ));
    let mut depths = Vec::new();
    let mut accept = |p: &CellProof, _: &mut Budget<'_>| {
        depths.push(p.speed);
        Ok(p.position <= 1e-7 && p.derivative < p.speed)
    };
    assert!(certify(&support, &uv, &target, 1e-7, &mut accept, &mut b).unwrap());
    assert!(depths.len() >= 2);
}
#[test]
fn a_displaced_target_or_a_wrong_radius_is_never_certified() {
    let mut b = Budget::new(512_000_000, None);
    let (support, uv, target) = scene(0.25, 0., &mut b);
    let mut moved = target;
    let mut controls = target[2].controls().to_vec();
    controls[4] = controls[4] + I::point(1e-5);
    moved[2] = P::new(&controls, &mut b).unwrap();
    assert!(!certify(&support, &uv, &moved, 1e-7, &mut strict(1e-7), &mut b).unwrap());
    let wrong = Support {
        radius: I::point(3.165 + 1e-6),
        ..support
    };
    assert!(!certify(&wrong, &uv, &target, 1e-7, &mut strict(1e-7), &mut b).unwrap());
}
#[test]
fn reversed_target_direction_fails_the_independent_derivative_check() {
    let mut b = Budget::new(512_000_000, None);
    let (support, uv, target) = scene(0.25, 0., &mut b);
    // Same image set traversed backwards: positions differ, and even a
    // permissive position limit cannot admit the opposite tangent.
    let reversed = target.map(|p| {
        let mut c = p.controls().to_vec();
        c.reverse();
        P::new(&c, &mut b).unwrap()
    });
    let mut accept = |p: &CellProof, _: &mut Budget<'_>| Ok(p.derivative < p.speed);
    assert!(!certify(&support, &uv, &reversed, 1e-7, &mut accept, &mut b).unwrap());
}
#[test]
fn degrees_above_seven_and_non_positive_cone_radius_refuse() {
    let mut b = Budget::new(512_000_000, None);
    let (support, uv, target) = scene(0.25, 0., &mut b);
    let mut high = uv;
    high[0] = uv[0].elevate(8, &mut b).unwrap();
    assert_eq!(
        certify(&support, &high, &target, 1e-7, &mut strict(1e-7), &mut b).err(),
        Some(ReadError::Unsupported)
    );
    let collapsed = Support {
        radius: I::point(-1.),
        ..support
    };
    assert_eq!(
        certify(&collapsed, &uv, &target, 1e-7, &mut strict(1e-7), &mut b).err(),
        Some(ReadError::Unsupported)
    );
}
#[test]
fn every_charge_is_metered_and_exhaustion_is_a_limit() {
    let mut setup = Budget::new(512_000_000, None);
    let (support, uv, target) = scene(1.2, 0., &mut setup);
    let mut full = Budget::new(512_000_000, None);
    assert!(certify(&support, &uv, &target, 1e-7, &mut strict(1e-7), &mut full).unwrap());
    let used = full.used();
    let mut short = Budget::new(used - 1, None);
    assert_eq!(
        certify(&support, &uv, &target, 1e-7, &mut strict(1e-7), &mut short).err(),
        Some(ReadError::Limit)
    );
    let mut exact = Budget::new(used, None);
    assert!(certify(&support, &uv, &target, 1e-7, &mut strict(1e-7), &mut exact).unwrap());
}
#[test]
fn an_end_beyond_tolerance_refuses_without_bisecting() {
    let mut b = Budget::new(512_000_000, None);
    let (support, uv, target) = scene(1.2, 0., &mut b);
    // Move the whole target 2e-7 along the axis: every END is now farther
    // than tol from the exact image, so the first piece decides.
    let shifted = [
        target[0],
        target[1],
        target[2]
            .add(&P::constant(I::point(2e-7), &mut b).unwrap(), &mut b)
            .unwrap(),
    ];
    let mut pieces = 0;
    let mut accept = |_: &CellProof, _: &mut Budget<'_>| {
        pieces += 1;
        Ok(true)
    };
    let before = b.used();
    assert!(!certify(&support, &uv, &shifted, 1e-7, &mut accept, &mut b).unwrap());
    assert_eq!(pieces, 0);
    let refusal = b.used() - before;
    let mut fresh = Budget::new(512_000_000, None);
    assert!(certify(&support, &uv, &target, 1e-7, &mut strict(1e-7), &mut fresh).unwrap());
    assert!(refusal < fresh.used());
}
#[test]
fn the_refusal_scan_rejects_a_displaced_target_and_never_rejects_a_true_one() {
    let mut b = Budget::new(512_000_000, None);
    for span in [0.25, 1.2] {
        let (support, uv, target) = scene(span, std::f64::consts::FRAC_PI_3.portable_tan(), &mut b);
        assert!(scan(&support, &uv, &target, 1e-7, &mut b).unwrap());
        let moved = [
            target[0],
            target[1]
                .add(&P::constant(I::point(3e-7), &mut b).unwrap(), &mut b)
                .unwrap(),
            target[2],
        ];
        assert!(!scan(&support, &uv, &moved, 1e-7, &mut b).unwrap());
    }
    let (support, uv, target) = scene(0.25, 0., &mut b);
    let mut high = target;
    high[0] = target[0].elevate(8, &mut b).unwrap();
    assert_eq!(
        scan(&support, &uv, &high, 1e-7, &mut b).err(),
        Some(ReadError::Unsupported)
    );
}
