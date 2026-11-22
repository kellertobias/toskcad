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
fn run(text: &str, resource: usize, interrupt: Option<&dyn Fn() -> bool>) -> Result<(f64, usize)> {
    let source = Source::parse(text, 1e-7)?;
    let mut setup = Budget::new(512_000_000, None);
    let mut count = 0;
    let mut maximum = 0f64;
    let mut used = 0;
    for binding in bindings(&source, &mut setup)? {
        let e = edge(&source, binding, &mut setup)?;
        let d = uv::definition(&source, binding.pcurve, &mut setup)?;
        if !matches!(&e.curve, Curve::Ellipse { .. }) {
            continue;
        }
        let surface = source.surface(d.surface)?;
        let mut budget = Budget::new(resource, interrupt);
        if let Some(proof) = certify(&source, binding.curve, &d, &surface, &e, 1e-7, &mut budget)? {
            count += 1;
            match proof {
                Proof::Certified(error) => maximum = maximum.max(error),
                _ => return Err(ReadError::Unsupported),
            }
            used = used.max(budget.used());
        }
    }
    check(count == 1, ReadError::Unsupported)?;
    Ok((maximum, used))
}
#[test]
fn actual_original_floor_ellipse_is_certified_without_a_cut_recipe() {
    for posed in [false, true] {
        let text = input(posed);
        let (error, used) = run(&text, 512_000_000, None).unwrap();
        assert!(error <= 1e-7);
        assert!(used < 512_000_000);
        assert!(matches!(run(&text, used - 1, None), Err(ReadError::Limit)));
        assert!(run(&text, used, None).is_ok());
        assert!(matches!(
            run(&text, 512_000_000, Some(&|| true)),
            Err(ReadError::Limit)
        ));
        let callbacks = std::cell::Cell::new(0usize);
        let observe = || {
            callbacks.set(callbacks.get() + 1);
            false
        };
        run(&text, 512_000_000, Some(&observe)).unwrap();
        let final_call = callbacks.get();
        callbacks.set(0);
        let cancel = || {
            callbacks.set(callbacks.get() + 1);
            callbacks.get() == final_call
        };
        assert!(matches!(
            run(&text, 512_000_000, Some(&cancel)),
            Err(ReadError::Limit)
        ));
        assert!(run(&text, 512_000_000, None).is_ok());
    }
}
#[test]
fn a_changed_original_uv_control_does_not_earn_a_curve_certificate() {
    let text = input(false);
    let source = Source::parse(&text, 1e-7).unwrap();
    let mut budget = Budget::new(512_000_000, None);
    let binding = bindings(&source, &mut budget)
        .unwrap()
        .into_iter()
        .find(|b| {
            matches!(
                edge(&source, *b, &mut budget).unwrap().curve,
                Curve::Ellipse { .. }
            )
        })
        .unwrap();
    let d = uv::definition(&source, binding.pcurve, &mut budget).unwrap();
    let fields = fields(&source, d.curve_id, 1, &mut budget).unwrap();
    let controls = list(fields.controls).unwrap();
    let id = reference(&controls[5]).unwrap();
    let statement = source.simple(id).unwrap().source.unwrap();
    let [_, values] = source.args(id, "CARTESIAN_POINT").unwrap() else {
        panic!()
    };
    let values = list(values).unwrap();
    let Value::Number(u) = values[0] else {
        panic!()
    };
    let Value::Number(v) = values[1] else {
        panic!()
    };
    let changed = text.replace(
        statement,
        &format!("CARTESIAN_POINT('',({u:.17E},{:.17E}))", v + 0.1),
    );
    assert_ne!(text, changed);
    Source::parse(&changed, 1e-7).unwrap();
    assert!(matches!(
        run(&changed, 512_000_000, None),
        Err(ReadError::Unsupported)
    ));
}
