use super::super::super::{Binding, bindings, edge, uv};
use super::*;
use multicad_precision::PortableMath as _;
const TOL: f64 = 1e-7;
fn text(left: bool, posed: bool) -> String {
    let tool = crate::helical_thread::Thread::new(crate::helical_thread::Spec {
        major_radius_mm: 2.942,
        pitch_mm: 1.,
        length_mm: 6.6,
        depth_mm: 0.5,
        crest_width_mm: 0.125,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: left,
        root_radius_mm: None,
    })
    .unwrap();
    let body = crate::thread_entry_chamfer::Chamfer::new(tool, 0.4, None).unwrap();
    let (o, u, v) = if posed {
        ([17., -8., 13.], [0.8, 0.6, 0.], [-0.36, 0.48, 0.8])
    } else {
        ([0.; 3], [1., 0., 0.], [0., 1., 0.])
    };
    // This must be the corrected actual AXIAL writer, never a repair of
    // archived bad bytes. Its own full graph/source-known verifier also runs.
    let text = crate::step::write_thread_entry_chamfer(&body, o, u, v, TOL).unwrap();
    crate::step_read::verify_thread_entry_chamfer(&text, &body, o, u, v, TOL).unwrap();
    text
}
struct Pair {
    binding: Binding,
    curve: u32,
    definition: uv::Definition,
    surface: Surface,
    edge: Edge,
}
fn pairs(source: &Source<'_>) -> Vec<Pair> {
    let mut budget = Budget::new(512_000_000, None);
    let mut out = Vec::new();
    for binding in bindings(source, &mut budget).unwrap() {
        let edge = edge(source, binding, &mut budget).unwrap();
        let definition = uv::definition(source, binding.pcurve, &mut budget).unwrap();
        let surface = source.surface(definition.surface).unwrap();
        let (
            Surface::Cone { .. },
            Curve::Nurbs { curve: uv, .. },
            Curve::Nurbs { curve: target, .. },
        ) = (&surface, &definition.geometry, &edge.curve)
        else {
            continue;
        };
        if matches!((uv.degree(), target.degree()), (3, 3) | (1, 11)) {
            out.push(Pair {
                binding,
                curve: binding.curve,
                definition,
                surface,
                edge,
            });
        }
    }
    out
}
#[test]
fn actual_corrected_writer_cubic_cone_traces_and_degree_eleven_rails_earn_full_fields() {
    for left in [false, true] {
        for posed in [false, true] {
            let text = text(left, posed);
            let source = Source::parse(&text, TOL).unwrap();
            let mut budget = Budget::new(512_000_000, None);
            let mut cubic = 0;
            let mut rail = 0;
            for pair in pairs(&source) {
                let actual = certify(
                    &source,
                    pair.curve,
                    &pair.definition,
                    &pair.surface,
                    &pair.edge,
                    TOL,
                    &mut budget,
                )
                .unwrap();
                assert!(
                    matches!(actual,Some(Proof::Certified(error))if error<=TOL),
                    "unearned left={left} posed={posed}"
                );
                if let Curve::Nurbs { curve, .. } = &pair.definition.geometry {
                    if curve.degree() == 1 {
                        rail += 1
                    } else {
                        cubic += 1
                    }
                }
            }
            assert!(cubic >= 2 && rail >= 1);
            assert!(budget.used() < 512_000_000);
        }
    }
}
#[test]
fn changed_actual_cone_field_is_not_hidden_by_parameter_or_tolerance_guesses() {
    let text = text(false, false);
    let source = Source::parse(&text, TOL).unwrap();
    let pair = pairs(&source)
        .into_iter()
        .find(|p| matches!(&p.definition.geometry,Curve::Nurbs{curve,..}if curve.degree()==3))
        .unwrap();
    let id = pair.definition.surface;
    let record = source.simple(id).unwrap();
    let [_, axis, Value::Number(radius), angle] = record.args.as_slice() else {
        panic!()
    };
    let axis = reference(axis).unwrap();
    let Value::Number(angle) = angle else {
        panic!()
    };
    let changed = text.replacen(
        record.source.unwrap(),
        &format!(
            "CONICAL_SURFACE('',#{axis},{:.16E},{angle:.16E})",
            radius + 0.2
        ),
        1,
    );
    assert_ne!(changed, text);
    let actual = Source::parse(&changed, TOL).unwrap();
    let pair = pairs(&actual)
        .into_iter()
        .find(|p| {
            p.definition.surface == id
                && matches!(&p.definition.geometry,Curve::Nurbs{curve,..}if curve.degree()==3)
        })
        .unwrap();
    assert!(matches!(
        certify(
            &actual,
            pair.curve,
            &pair.definition,
            &pair.surface,
            &pair.edge,
            TOL,
            &mut Budget::new(512_000_000, None)
        ),
        Ok(Some(Proof::Unsupported)) | Err(ReadError::Unsupported)
    ));
}
#[test]
fn whole_target_endpoint_identity_and_original_chart_are_not_inferred_from_cached_ranges() {
    let text = text(false, false);
    let source = Source::parse(&text, TOL).unwrap();
    let pair = pairs(&source)
        .into_iter()
        .find(|p| matches!(&p.definition.geometry,Curve::Nurbs{curve,..}if curve.degree()==3))
        .unwrap();
    let target = unwrapped(&source, pair.curve, &mut Budget::new(512_000_000, None)).unwrap();
    let f = fields(&source, target, 1, &mut Budget::new(512_000_000, None)).unwrap();
    let controls = list(f.controls).unwrap();
    let interior = reference(&controls[3]).unwrap();
    let edge_record = source.simple(pair.binding.edge).unwrap();
    let [_, _, old_end, _, _] = edge_record.args.as_slice() else {
        panic!()
    };
    let old_end = reference(old_end).unwrap();
    let original = edge_record.source.unwrap();
    let replacement = original.replacen(&format!(",#{old_end},"), ",#999999,", 1);
    assert_ne!(replacement, original);
    let partial = text.replacen(original, &replacement, 1);
    let at = partial.rfind("ENDSEC;").unwrap();
    let partial = format!(
        "{}#999999=VERTEX_POINT('',#{interior});\n{}",
        &partial[..at],
        &partial[at..]
    );
    let actual = Source::parse(&partial, TOL).unwrap();
    // The actual parsed EDGE_CURVE now ends at the exact first cubic knot
    // join's original point. Its positive finite window is earned by the
    // normal edge parser; no positions or cached ranges are fabricated.
    let segment = edge(&actual, pair.binding, &mut Budget::new(512_000_000, None)).unwrap();
    assert!(matches!(
        full_target(
            &actual,
            target,
            &segment,
            &mut Budget::new(512_000_000, None)
        ),
        Err(ReadError::Unsupported)
    ));
}
#[test]
fn actual_cubic_source_work_frontier_late_cancel_and_retry_preserve_operands() {
    let text = text(false, true);
    let source = Source::parse(&text, TOL).unwrap();
    let pair = pairs(&source)
        .into_iter()
        .find(|p| matches!(&p.definition.geometry,Curve::Nurbs{curve,..}if curve.degree()==3))
        .unwrap();
    let run = |b: &mut Budget<'_>| {
        certify(
            &source,
            pair.curve,
            &pair.definition,
            &pair.surface,
            &pair.edge,
            TOL,
            b,
        )
    };
    let mut full = Budget::new(512_000_000, None);
    let Some(Proof::Certified(first)) = run(&mut full).unwrap() else {
        panic!()
    };
    let used = full.used();
    let mut exact = Budget::new(used, None);
    assert!(matches!(run(&mut exact).unwrap(),Some(Proof::Certified(error))if error==first));
    assert_eq!(exact.remaining, 0);
    assert!(matches!(
        run(&mut Budget::new(used - 1, None)),
        Err(ReadError::Limit)
    ));
    let count = std::cell::Cell::new(0);
    let stop = || {
        count.set(count.get() + 1);
        count.get() == 40
    };
    let mut cancelled = Budget::new(512_000_000, Some(&stop));
    assert!(run(&mut cancelled).is_err());
    assert!(cancelled.cancelled.get());
    let mut retry = Budget::new(512_000_000, None);
    assert!(matches!(run(&mut retry).unwrap(),Some(Proof::Certified(error))if error==first));
    assert_eq!(retry.used(), used);
}

#[test]
fn hull_tangent_proof_accepts_agreeing_cells_and_never_reversed_ones() {
    let mut b = Budget::new(10_000_000, None);
    let poly = |values: &[f64], b: &mut Budget<'_>| {
        let values: Vec<I> = values.iter().map(|&x| I::point(x)).collect();
        P::new(&values, b).unwrap()
    };
    // A radius-2 band at z = 0 swept by 0.1 rad, and the cubic Bernstein
    // controls of its exact circular arc's Taylor cubic.
    let rho = poly(&[2.], &mut b);
    let theta = poly(&[0., 0.1 / 3., 0.2 / 3., 0.1], &mut b);
    let z = poly(&[0.], &mut b);
    let arc = |reverse: bool, b: &mut Budget<'_>| -> VecP {
        let point = |t: f64| {
            let a = 0.1 * if reverse { 1. - t } else { t };
            [2. * a.portable_cos(), 2. * a.portable_sin(), 0.]
        };
        // Interpolating cubic through four arc points, in Bernstein form.
        let p: Vec<[f64; 3]> = [0., 1. / 3., 2. / 3., 1.].map(point).to_vec();
        let control = |k: usize| -> [f64; 3] {
            std::array::from_fn(|i| match k {
                0 => p[0][i],
                3 => p[3][i],
                1 => (-5. * p[0][i] + 18. * p[1][i] - 9. * p[2][i] + 2. * p[3][i]) / 6.,
                _ => (2. * p[0][i] - 9. * p[1][i] + 18. * p[2][i] - 5. * p[3][i]) / 6.,
            })
        };
        std::array::from_fn(|i| {
            poly(
                &[control(0)[i], control(1)[i], control(2)[i], control(3)[i]],
                b,
            )
        })
    };
    let forward = arc(false, &mut b);
    assert!(tangent_by_hulls(&rho, &theta, &z, &forward, &mut b).unwrap());
    assert!(tangent_agrees(&rho, &theta, &z, &forward, 1, &mut b).unwrap());
    // The same image run backwards has the opposite tangent everywhere.
    let backward = arc(true, &mut b);
    assert!(!tangent_by_hulls(&rho, &theta, &z, &backward, &mut b).unwrap());
    assert!(!tangent_agrees(&rho, &theta, &z, &backward, 3, &mut b).unwrap());
    // The position bound holds for the agreeing cell and not a shifted one.
    let bound = position_error(&rho, &theta, &z, &forward, 4, &mut b).unwrap();
    assert!(bound < 1e-5, "{bound}");
    let mut shifted = forward;
    shifted[2] = poly(&[1e-3; 4], &mut b);
    assert!(position_error(&rho, &theta, &z, &shifted, 4, &mut b).unwrap() >= 1e-3);
}
