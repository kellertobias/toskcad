use super::super::c0_trace;
use super::*;
use crate::step_read::brep::pcurve::{bindings, edge, trim_surface, uv};

const TOL: f64 = 1e-7;

/// The kernel's own conical-entry export: its thread walls are degree-11
/// tensor supports with C0 interior knots, trimmed by cubic UV traces.
fn exported() -> String {
    let tool = crate::helical_thread::Thread::new(crate::helical_thread::Spec {
        major_radius_mm: 2.942,
        pitch_mm: 1.,
        length_mm: 6.6,
        depth_mm: 0.5,
        crest_width_mm: 0.125,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: false,
        root_radius_mm: None,
    })
    .unwrap();
    let body = crate::thread_entry_chamfer::Chamfer::new(tool, 0.4, None).unwrap();
    crate::step::write_thread_entry_chamfer(&body, [0.; 3], [1., 0., 0.], [0., 1., 0.], TOL)
        .unwrap()
}
/// The cheapest curved trace on a C0 tensor support: its outcome.
fn first_curved(text: &str, tol: f64) -> Option<Result<Option<Proof>>> {
    let source = Source::parse(text, TOL).unwrap();
    let mut setup = Budget::new(512_000_000, None);
    let mut best: Option<(usize, _)> = None;
    for binding in bindings(&source, &mut setup).unwrap() {
        let d = uv::definition(&source, binding.pcurve, &mut setup).unwrap();
        let s = trim_surface::read(&source, &d, &mut setup).unwrap();
        let (Surface::Nurbs { surface: net, .. }, Curve::Nurbs { curve, .. }) = (&s, &d.geometry)
        else {
            continue;
        };
        let c0 = (0..2).any(|axis| {
            let (_, m) = net.knots(axis);
            m.len() > 2 && m[1..m.len() - 1].contains(&net.degrees()[axis])
        });
        if c0 && curve.degree() >= 2 && best.as_ref().is_none_or(|(n, _)| curve.points().len() < *n)
        {
            best = Some((curve.points().len(), (binding, d, s)));
        }
    }
    let (_, (binding, d, s)) = best?;
    let e = edge(&source, binding, &mut setup).ok()?;
    let mut b = Budget::new(512_000_000, None);
    Some(c0_trace::certify(
        &source,
        binding.curve,
        &d,
        &s,
        &e,
        tol,
        &mut b,
    ))
}

#[test]
fn a_curved_wall_trace_on_a_c0_support_certifies_and_a_tighter_tolerance_refuses() {
    let text = exported();
    let Some(Ok(Some(Proof::Certified(error)))) = first_curved(&text, TOL) else {
        panic!("curved C0 trace was not certified");
    };
    assert!(error <= TOL);
    // The same trace cannot be certified to 1e-12 mm: never a false pass.
    assert!(matches!(
        first_curved(&text, 1e-12),
        Some(Ok(Some(Proof::Unsupported)))
    ));
}

#[test]
fn jets_and_restriction_agree_with_exact_polynomials() {
    let mut b = Budget::new(1_000_000, None);
    // p(t) = (1-t)^2 * 1 + 2t(1-t) * 3 + t^2 * 2 => 1 + 4t - 3t^2.
    let p = [[I::point(1.)], [I::point(3.)], [I::point(2.)]];
    let [j] = jets(&p, I::point(0.5), &mut b).unwrap();
    for (got, want) in j.iter().zip([2.25, 1., -3., 0., 0.]) {
        assert!(got.lo <= want && want <= got.hi);
    }
    let half = restrict(&p, I::point(0.), I::point(0.5), &mut b).unwrap();
    let [h] = jets(&half, I::point(1.), &mut b).unwrap();
    assert!(h[0].lo <= 2.25 && 2.25 <= h[0].hi);
    assert_eq!(
        restrict::<1>(&[], I::point(0.), I::point(1.), &mut b).err(),
        Some(ReadError::Unsupported)
    );
}

#[test]
fn one_casteljau_run_gives_every_low_derivative_and_zero_above_the_degree() {
    let mut b = Budget::new(10_000_000, None);
    // A degree-11 curve: compare with re-differenced control nets.
    let points: Vec<[I; 2]> = (0..12)
        .map(|i| {
            let x = f64::from(i);
            [I::point((0.3 * x).sin() + 0.1 * x * x), I::point(x.sqrt())]
        })
        .collect();
    for t in [0., 0.25, 0.5, 0.875, 1.] {
        let t = I::point(t);
        let fast = derivatives(&points, t, 4, &mut b).unwrap();
        let mut net = points.clone();
        for value in &fast {
            let slow = evaluate(&net, t, &mut b).unwrap();
            for axis in 0..2 {
                let (a, c) = (value[axis], slow[axis]);
                let scale = 1. + a.magnitude().max(c.magnitude());
                assert!(a.valid() && (a.midpoint() - c.midpoint()).abs() <= 1e-9 * scale);
                // Both are enclosures of the same exact value.
                assert!(a.lo <= c.hi + 1e-9 * scale && c.lo <= a.hi + 1e-9 * scale);
            }
            net = derivative(&net);
        }
    }
    // A quadratic has no third derivative.
    let quadratic = [[I::point(1.)], [I::point(3.)], [I::point(2.)]];
    let all = derivatives(&quadratic, I::point(0.5), 4, &mut b).unwrap();
    for (got, want) in all.iter().zip([2.25, 1., -6., 0.]) {
        assert!(got[0].lo <= want && want <= got[0].hi, "{got:?} {want}");
    }
    assert_eq!(
        derivatives::<1>(&[], I::point(0.5), 4, &mut b).err(),
        Some(ReadError::Unsupported)
    );
}
