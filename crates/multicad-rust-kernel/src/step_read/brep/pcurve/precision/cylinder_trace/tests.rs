use super::super::super::{bindings, edge, uv};
use super::*;
fn input(posed: bool) -> String {
    let cut = crate::helical_tube_cut::Cut::new(
        crate::helical_tube_cut::Spec {
            stock_radius_mm: 8.,
            stock_height_mm: 30.,
            tube_radius_mm: 0.6,
            pitch_mm: 2.,
            path_height_mm: 15.,
            start_height_mm: 16.2,
            start_angle_rad: 0.,
        },
        1e-7,
    )
    .unwrap();
    let (o, u, v) = if posed {
        ([17., 4., 0.], [0.6, 0., -0.8], [0., 1., 0.])
    } else {
        ([0.; 3], [1., 0., 0.], [0., 1., 0.])
    };
    crate::step::write_helical_tube_cut(cut, o, u, v, 1e-7).unwrap()
}
#[test]
fn actual_original_cylinder_rails_earn_whole_c0_and_independent_positive_tangent() {
    for posed in [false, true] {
        let text = input(posed);
        let source = Source::parse(&text, 1e-7).unwrap();
        let mut budget = Budget::new(512_000_000, None);
        let mut checked = 0;
        for binding in bindings(&source, &mut budget).unwrap() {
            let e = edge(&source, binding, &mut budget).unwrap();
            let d = uv::definition(&source, binding.pcurve, &mut budget).unwrap();
            let s = source.surface(d.surface).unwrap();
            match certify(&source, binding.curve, &d, &s, &e, 1e-7, &mut budget).unwrap() {
                Some(Proof::Certified(error)) => {
                    assert!(error <= 1e-7);
                    checked += 1;
                }
                None => {}
                _ => panic!(
                    "unearned rail posed={posed} edge={} pcurve={}",
                    binding.edge, binding.pcurve
                ),
            }
        }
        assert_eq!(checked, 4, "both seam uses and the two finite source rails");
        assert!(budget.used() < 512_000_000);
    }
}
fn wave() -> Wave {
    Wave {
        origin: [I::point(0.); 3],
        axes: [
            cached([1., 0., 0.]),
            cached([0., 1., 0.]),
            cached([0., 0., 1.]),
        ],
        radius: I::point(8.),
        theta: I::point(0.3),
        height: I::point(1.),
        angle_step: I::point(1.),
        height_step: I::point(0.5),
    }
}
#[test]
fn actual_polynomial_perturbation_and_independent_derivative_refuse_without_budget_escape() {
    let w = wave();
    let mut b = Budget::new(512_000_000, None);
    let h = hermite11(&w, &mut b).unwrap();
    let mut controls = h.map(|p| [p[0], p[1], p[2], I::point(1.)]);
    assert!(matches!(
        cell(&w, &controls, 1e-7, &mut b).unwrap(),
        Proof::Certified(_)
    ));
    controls[5][0] = controls[5][0] + I::point(0.01);
    assert!(matches!(
        cell(&w, &controls, 1e-7, &mut b).unwrap(),
        Proof::Unsupported
    ));
    assert!(matches!(
        cell(&w, &controls, 1e-7, &mut Budget::new(0, None)),
        Err(ReadError::Limit)
    ));
    let cancel = || true;
    let mut interrupted = Budget::new(512_000_000, Some(&cancel));
    assert!(matches!(
        cell(&w, &controls, 1e-7, &mut interrupted),
        Err(ReadError::Limit)
    ));
    assert!(interrupted.cancelled.get());
    let degenerate = Wave {
        angle_step: I::point(0.),
        height_step: I::point(0.),
        ..wave()
    };
    let flat = hermite11(&degenerate, &mut b)
        .unwrap()
        .map(|p| [p[0], p[1], p[2], I::point(1.)]);
    assert!(matches!(
        cell(&degenerate, &flat, 1e-7, &mut b).unwrap(),
        Proof::Unsupported
    ));
    // A position error below tolerance must NOT hide a reversed/uncertain
    // parameter flow on a very slow source arc.
    let slow = Wave {
        angle_step: I::point(1e-12),
        height_step: I::point(0.),
        ..wave()
    };
    let mut almost = hermite11(&slow, &mut b)
        .unwrap()
        .map(|p| [p[0], p[1], p[2], I::point(1.)]);
    almost[5][0] = almost[5][0] + I::point(5e-9);
    assert!(matches!(
        cell(&slow, &almost, 1e-7, &mut b).unwrap(),
        Proof::Unsupported
    ));
}
#[test]
fn original_constructor_work_frontier_and_late_cancel_preserve_same_source_retry() {
    let text = input(true);
    let source = Source::parse(&text, 1e-7).unwrap();
    let mut preparation = Budget::new(512_000_000, None);
    let (id, definition, surface, edge) = bindings(&source, &mut preparation)
        .unwrap()
        .into_iter()
        .find_map(|binding| {
            let edge = edge(&source, binding, &mut preparation).unwrap();
            let definition = uv::definition(&source, binding.pcurve, &mut preparation).unwrap();
            let surface = source.surface(definition.surface).unwrap();
            if matches!(surface, Surface::Cylinder { .. })
                && matches!(&definition.geometry, Curve::Nurbs { curve, closed: false } if curve.degree()==1 && curve.points().len()==2)
                && matches!(&edge.curve, Curve::Nurbs { curve, closed: false } if curve.degree()==11)
            {
                Some((binding.curve, definition, surface, edge))
            } else {
                None
            }
        })
        .expect("the unchanged public writer supplies actual cylinder rails");
    let run =
        |budget: &mut Budget<'_>| certify(&source, id, &definition, &surface, &edge, 1e-7, budget);
    let mut full = Budget::new(512_000_000, None);
    assert!(matches!(run(&mut full).unwrap(), Some(Proof::Certified(_))));
    let used = full.used();
    assert!(matches!(
        run(&mut Budget::new(used - 1, None)),
        Err(ReadError::Limit)
    ));
    assert!(matches!(
        run(&mut Budget::new(used, None)).unwrap(),
        Some(Proof::Certified(_))
    ));
    let calls = std::cell::Cell::new(0usize);
    let observe = || {
        calls.set(calls.get() + 1);
        false
    };
    assert!(matches!(
        run(&mut Budget::new(used, Some(&observe))).unwrap(),
        Some(Proof::Certified(_))
    ));
    let total = calls.get();
    calls.set(0);
    let cancel = || {
        calls.set(calls.get() + 1);
        calls.get() == total
    };
    let mut interrupted = Budget::new(used, Some(&cancel));
    assert!(matches!(run(&mut interrupted), Err(ReadError::Limit)));
    assert!(interrupted.cancelled.get());
    assert!(matches!(
        run(&mut Budget::new(used, None)).unwrap(),
        Some(Proof::Certified(_))
    ));
}
