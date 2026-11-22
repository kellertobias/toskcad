//! Independent generic cone semantics beside source-known full readback.
use super::*;
const TOL: f64 = 1e-7;
fn body(left: bool) -> crate::thread_entry_chamfer::Chamfer {
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
    crate::thread_entry_chamfer::Chamfer::new(tool, 0.4, None).unwrap()
}
fn placement(posed: bool) -> [[f64; 3]; 3] {
    if posed {
        [[17., -8., 13.], [0.8, 0.6, 0.], [-0.36, 0.48, 0.8]]
    } else {
        [[0.; 3], [1., 0., 0.], [0., 1., 0.]]
    }
}
struct Group {
    uv_ids: Vec<u32>,
    space_ids: Vec<u32>,
    support: u32,
    height: Vec<f64>,
    axial: Vec<f64>,
}
fn uv_point(value: &Value<'_>) -> [f64; 2] {
    let values = list(value).unwrap();
    assert_eq!(values.len(), 2);
    [number(&values[0]).unwrap(), number(&values[1]).unwrap()]
}
fn groups(source: &Source<'_>) -> Vec<Group> {
    let mut output = Vec::new();
    for parts in source.entities.values() {
        let [e] = parts.as_slice() else { continue };
        if !matches!(e.name, "SURFACE_CURVE" | "SEAM_CURVE") {
            continue;
        }
        let [_, curve, associations, _] = e.args.as_slice() else {
            panic!()
        };
        let Curve::Nurbs { curve: actual, .. } =
            source.curve(reference(curve).unwrap(), 0).unwrap()
        else {
            continue;
        };
        for pc in refs(associations).unwrap() {
            let [_, support, representation] = source.args(pc, "PCURVE").unwrap() else {
                panic!()
            };
            let support = reference(support).unwrap();
            let Ok([_, axis, _, _]) = source.args(support, "CONICAL_SURFACE") else {
                continue;
            };
            let frame = source.frame(reference(axis).unwrap()).unwrap();
            let [_, items, _] = source
                .args(
                    reference(representation).unwrap(),
                    "DEFINITIONAL_REPRESENTATION",
                )
                .unwrap()
            else {
                panic!()
            };
            let uv_items = refs(items).unwrap();
            let [uv] = uv_items.as_slice() else { panic!() };
            let [_, degree, controls, _, _, _, multiplicity, knots, _] =
                source.args(*uv, "B_SPLINE_CURVE_WITH_KNOTS").unwrap()
            else {
                panic!()
            };
            let Value::Number(degree) = degree else {
                panic!()
            };
            if *degree != 3. || actual.degree() != 3 {
                continue;
            }
            let ids = refs(controls).unwrap();
            let heights = ids
                .iter()
                .map(|&id| {
                    let [_, values] = source.args(id, "CARTESIAN_POINT").unwrap() else {
                        panic!()
                    };
                    uv_point(values)[1]
                })
                .collect::<Vec<_>>();
            assert_eq!(heights.len(), actual.points().len());
            let declared_knots = list(knots)
                .unwrap()
                .iter()
                .map(number)
                .collect::<Result<Vec<_>>>()
                .unwrap();
            let declared_m = list(multiplicity)
                .unwrap()
                .iter()
                .map(|v| number(v).unwrap() as usize)
                .collect::<Vec<_>>();
            assert_eq!(
                actual.knots(),
                (declared_knots.as_slice(), declared_m.as_slice())
            );
            let axial = actual
                .points()
                .iter()
                .map(|p| (0..3).map(|i| (p[i] - frame.origin[i]) * frame.z[i]).sum())
                .collect::<Vec<_>>();
            // Independent generic support uses AXIAL v, including under pose.
            let (at, _) = source.surface(support).unwrap().evaluate([0., 0.4]);
            let z: f64 = (0..3).map(|i| (at[i] - frame.origin[i]) * frame.z[i]).sum();
            assert!((z - 0.4).abs() < 1e-12);
            let [_, _, space, _, _, _, _, _, _] = source
                .args(reference(curve).unwrap(), "B_SPLINE_CURVE_WITH_KNOTS")
                .unwrap()
            else {
                panic!()
            };
            output.push(Group {
                uv_ids: ids,
                space_ids: refs(space).unwrap(),
                support,
                height: heights,
                axial,
            });
        }
    }
    output
}
#[test]
fn emitted_height_coefficients_follow_generic_axial_cone_for_both_hands_and_pose() {
    for left in [false, true] {
        for posed in [false, true] {
            let body = body(left);
            let [o, u, v] = placement(posed);
            let text = crate::step::write_thread_entry_chamfer(&body, o, u, v, TOL).unwrap();
            crate::step_read::verify_thread_entry_chamfer(&text, &body, o, u, v, TOL).unwrap();
            let source = Source::parse(&text, TOL).unwrap();
            let groups = groups(&source);
            assert!(!groups.is_empty());
            for group in groups {
                // Matching complete positive B-spline bases make the axial
                // residual a convex combination of THESE actual coefficients.
                // This bounds the whole field, rather than sampled curve points.
                assert!(
                    group
                        .height
                        .iter()
                        .zip(group.axial)
                        .all(|(v, z)| (v - z).abs() < 1e-11)
                );
            }
        }
    }
}
#[test]
fn former_slant_height_cannot_be_repaired_by_reparameterizing_a_constant_height_trace() {
    let body = body(false);
    let [o, u, v] = placement(false);
    let text = crate::step::write_thread_entry_chamfer(&body, o, u, v, TOL).unwrap();
    let source = Source::parse(&text, TOL).unwrap();
    let group = groups(&source)
        .into_iter()
        .find(|g| {
            g.height[0] > 0.1
                && g.height.iter().all(|v| *v == g.height[0])
                && g.axial.iter().all(|z| (*z - g.axial[0]).abs() < 1e-12)
        })
        .unwrap();
    let mut edits = Vec::new();
    for &id in &group.uv_ids {
        let entity = source.simple(id).unwrap();
        let original = entity.source.unwrap();
        let [_, values] = entity.args.as_slice() else {
            panic!()
        };
        let [x, z] = uv_point(values);
        let old_slant = z / std::f64::consts::FRAC_PI_4.cos();
        let replacement = format!("CARTESIAN_POINT('',({x:.16E},{old_slant:.16E}))");
        let offset = original.as_ptr() as usize - text.as_ptr() as usize;
        edits.push((offset, original.len(), replacement));
    }
    edits.sort_by_key(|e| e.0);
    edits.dedup_by_key(|e| e.0);
    let mut bad = text.clone();
    for (at, len, replacement) in edits.into_iter().rev() {
        bad.replace_range(at..at + len, &replacement);
    }
    assert_ne!(text, bad);
    assert!(crate::step_read::verify_thread_entry_chamfer(&bad, &body, o, u, v, TOL).is_err());
    let actual = Source::parse(&bad, TOL).unwrap();
    let [_, axis, _, _] = actual.args(group.support, "CONICAL_SURFACE").unwrap() else {
        panic!()
    };
    let [_, origin, normal, _] = actual
        .args(reference(axis).unwrap(), "AXIS2_PLACEMENT_3D")
        .unwrap()
    else {
        panic!()
    };
    // This negative has actual exact identity axial placement, proved from
    // original numeric spellings, not decoded-coordinate coincidence.
    assert!(
        raw_values(&actual, reference(origin).unwrap())
            .iter()
            .all(|value| value == &q("0.0"))
    );
    assert_eq!(
        raw_values(&actual, reference(normal).unwrap()),
        vec![q("0.0"), q("0.0"), q("1.0")]
    );
    let min_v = group
        .uv_ids
        .iter()
        .map(|&id| raw_values(&actual, id)[1].clone())
        .min()
        .unwrap();
    let max_z = group
        .space_ids
        .iter()
        .map(|&id| raw_values(&actual, id)[2].clone())
        .max()
        .unwrap();
    // Every original constant-v curve point and every original 3D spline
    // point lies in its positive-basis coefficient slab. This exact rational
    // gap is a pointset lower bound, independent of parameterization/window.
    assert!(min_v - max_z > q("0.16"));
}

fn q(text: &str) -> num_rational::BigRational {
    let (n, d) = multicad_precision::Rational::parse_decimal(text)
        .unwrap()
        .ratio();
    num_rational::BigRational::new(n.into(), d.into())
}
fn raw_values(source: &Source<'_>, id: u32) -> Vec<num_rational::BigRational> {
    let record = source.simple(id).unwrap();
    part21::numeric_spellings(record)
        .unwrap()
        .into_iter()
        .map(q)
        .collect()
}
