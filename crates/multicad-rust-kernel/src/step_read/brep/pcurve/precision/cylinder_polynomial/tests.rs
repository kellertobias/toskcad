use super::*;
fn q(n: i64, d: i64) -> BigRational {
    BigRational::new(n.into(), d.into())
}
fn choose(n: usize, k: usize) -> BigRational {
    let fact = |n| (1..=n).fold(BigInt::from(1), |v, i| v * BigInt::from(i));
    BigRational::new(fact(n), fact(k) * fact(n - k))
}
fn scene(b: &mut Budget<'_>) -> Result<(RawFrame, [P; 2], [P; 3])> {
    let h = q(1, 1000);
    let uv = [
        P::new(&[I::point(0.), enclose(&h, b)?], b)?,
        P::new(&[I::point(0.), I::point(1.)], b)?,
    ];
    let h2 = &h * &h;
    let h3 = &h2 * &h;
    let h4 = &h3 * &h;
    let powers = [
        vec![q(2, 1), q(0, 1), -&h2, q(0, 1), &h4 / q(12, 1)],
        vec![q(0, 1), q(2, 1) * &h, q(0, 1), -&h3 / q(3, 1), q(0, 1)],
        vec![q(0, 1), q(1, 1), q(0, 1), q(0, 1), q(0, 1)],
    ];
    let mut target = [P::constant(I::point(0.), b)?; 3];
    for axis in 0..3 {
        let mut controls = [I::point(0.); 5];
        for (i, c) in controls.iter_mut().enumerate() {
            let exact = (0..=i).fold(q(0, 1), |s, j| {
                s + &powers[axis][j] * choose(i, j) / choose(4, j)
            });
            *c = enclose(&exact, b)?;
        }
        target[axis] = P::new(&controls, b)?;
    }
    let frame = RawFrame {
        origin: cached([0.; 3]),
        axes: [
            cached([1., 0., 0.]),
            cached([0., 1., 0.]),
            cached([0., 0., 1.]),
        ],
    };
    Ok((frame, uv, target))
}
#[test]
fn exact_power_series_target_is_bounded_over_the_whole_cell() {
    let mut b = Budget::new(512_000_000, None);
    let (f, uv, target) = scene(&mut b).unwrap();
    let result = cell(&f, I::point(2.), &uv, &target, &mut b).unwrap();
    assert!(result.position < 1e-7 && result.derivative < result.speed && result.speed > 0.);
    let mut changed = target;
    let mut controls = target[2].controls().to_vec();
    controls[2] = controls[2] + I::point(0.1);
    changed[2] = P::new(&controls, &mut b).unwrap();
    let result = cell(&f, I::point(2.), &uv, &changed, &mut b).unwrap();
    assert!(result.position > 1e-3);
}
#[test]
fn actual_math_has_exact_work_and_late_cancellation_frontiers() {
    let mut setup = Budget::new(512_000_000, None);
    let (f, uv, target) = scene(&mut setup).unwrap();
    let mut b = Budget::new(512_000_000, None);
    cell(&f, I::point(2.), &uv, &target, &mut b).unwrap();
    let used = b.used();
    assert!(matches!(
        cell(
            &f,
            I::point(2.),
            &uv,
            &target,
            &mut Budget::new(used - 1, None)
        ),
        Err(ReadError::Limit)
    ));
    cell(&f, I::point(2.), &uv, &target, &mut Budget::new(used, None)).unwrap();
    let calls = std::cell::Cell::new(0);
    let observe = || {
        calls.set(calls.get() + 1);
        false
    };
    cell(
        &f,
        I::point(2.),
        &uv,
        &target,
        &mut Budget::new(512_000_000, Some(&observe)),
    )
    .unwrap();
    let last = calls.get();
    calls.set(0);
    let cancel = || {
        calls.set(calls.get() + 1);
        calls.get() == last
    };
    let mut b = Budget::new(512_000_000, Some(&cancel));
    assert!(matches!(
        cell(&f, I::point(2.), &uv, &target, &mut b),
        Err(ReadError::Limit)
    ));
    assert!(b.cancelled.get());
    cell(
        &f,
        I::point(2.),
        &uv,
        &target,
        &mut Budget::new(512_000_000, None),
    )
    .unwrap();
}

#[test]
fn installed_dip_truncated_spline_curves_and_whole_validator_are_qualified() {
    let path = "/Applications/KiCad/KiCad.app/Contents/SharedSupport/3dmodels/Package_DIP.3dshapes/DIP-32_W7.62mm.step";
    let Ok(text) = std::fs::read_to_string(path) else {
        eprintln!("unavailable: external KiCad DIP fixture");
        return;
    };
    let source = Source::parse(&text, 1e-7).unwrap();
    let mut setup = Budget::new(512_000_000, None);
    let mut count = 0;
    for binding in super::super::super::bindings(&source, &mut setup).unwrap() {
        if ![29066, 29678].contains(&binding.pcurve) {
            continue;
        }
        let edge = super::super::super::edge(&source, binding, &mut setup).unwrap();
        let definition =
            super::super::super::uv::definition(&source, binding.pcurve, &mut setup).unwrap();
        let surface = source.surface(definition.surface).unwrap();
        let mut budget = Budget::new(512_000_000, None);
        let proof = certify(
            &source,
            binding.curve,
            &definition,
            &surface,
            &edge,
            1e-7,
            &mut budget,
        )
        .unwrap();
        match proof {
            Some(super::super::super::Proof::Certified(error)) => eprintln!(
                "DIP spline edge{} pc{} C0={error} work={}",
                binding.edge,
                binding.pcurve,
                budget.used()
            ),
            _ => panic!(
                "unqualified original DIP spline edge{}pc{}",
                binding.edge, binding.pcurve
            ),
        }
        count += 1;
    }
    assert_eq!(count, 2);
    let report = super::super::super::validate_pcurves(&text, 1e-7).unwrap();
    eprintln!(
        "whole DIP edges{} pcurves{} C0={} work={}",
        report.edges_checked,
        report.pcurves_checked,
        report.certified_maximum_deviation_mm,
        report.work_used
    );
}
