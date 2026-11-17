use super::*;
fn text() -> String {
    use crate::{
        head_profile::Head,
        helical_thread::{Spec, Thread},
        hexalobular_drive::{Cone, Drive, TangentSpec},
        profiled_screw::Screw,
        screw_recess::Recess,
        step::Writer,
    };
    let drive = Drive::from_authored_tangent(
        TangentSpec {
            major_diameter_mm: 2.8,
            minor_diameter_mm: 2.05,
            outer_lobe_radius_mm: 0.28,
        },
        Cone {
            tip_depth_mm: 1.1 + 1.4 / 3f64.sqrt(),
            depth_per_radius: 1. / 3f64.sqrt(),
        },
    )
    .unwrap();
    let recess = Recess::conical_hexalobular(drive, 1e-7, 100_000).unwrap();
    let thread = Thread::new(Spec {
        major_radius_mm: 1.5,
        pitch_mm: 0.5,
        length_mm: 3.,
        depth_mm: 0.2,
        crest_width_mm: 0.0625,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: false,
        root_radius_mm: None,
    })
    .unwrap();
    let torx = Screw::new_with_recess(
        thread,
        1.5,
        Head::Crowned {
            radius_mm: 2.8,
            height_mm: 2.4,
            crown_radius_mm: 5.,
            top_land_radius_mm: 2.05,
            neck_radius_mm: 0.125,
        },
        0.,
        0.1,
        recess,
    )
    .unwrap();
    let mut file = Writer::new();
    file.append_profiled_screw(&torx, [-20., 0., 0.], [1., 0., 0.], [0., 1., 0.], 1e-7)
        .unwrap();
    file.finish().unwrap()
}
struct Pair {
    curve: u32,
    definition: uv::Definition,
    surface: Surface,
    edge: Edge,
}
fn pairs(source: &Source<'_>) -> Vec<Pair> {
    let mut output = Vec::new();
    let mut budget = Budget::new(MAX_WORK, None);
    for parts in source.entities.values() {
        let [entity] = parts.as_slice() else { continue };
        if entity.name != "EDGE_CURVE" {
            continue;
        }
        let [_, a, b, curve, same] = entity.args.as_slice() else {
            panic!()
        };
        let curve = reference(curve).unwrap();
        let record = source.simple(curve).unwrap();
        if !matches!(
            record.name,
            "SURFACE_CURVE" | "SEAM_CURVE" | "INTERSECTION_CURVE"
        ) {
            continue;
        }
        let [_, _, associations, _] = record.args.as_slice() else {
            panic!()
        };
        let mut vertices = [reference(a).unwrap(), reference(b).unwrap()];
        if !sense(same).unwrap() {
            vertices.swap(0, 1)
        }
        let points = vertices.map(|id| {
            let [_, point] = source.args(id, "VERTEX_POINT").unwrap() else {
                panic!()
            };
            let point = reference(point).unwrap();
            equivalent_vertex(source, point, &mut budget)
                .unwrap()
                .unwrap_or_else(|| source.point(point).unwrap())
        });
        let edge = Edge::new(
            source.curve(curve, 0).unwrap(),
            vertices.map(|v| v as usize),
            points,
            1e-7,
        )
        .unwrap();
        for pc in refs(associations).unwrap() {
            if source.simple(pc).unwrap().name != "PCURVE" {
                continue;
            }
            let definition = uv::definition(source, pc, &mut budget).unwrap();
            let surface = source.surface(definition.surface).unwrap();
            if matches!(surface, Surface::Sphere { .. } | Surface::Cone { .. }) {
                output.push(Pair {
                    curve,
                    definition,
                    surface,
                    edge: edge.clone(),
                });
            }
        }
    }
    output
}
#[test]
fn actual_crowned_torx_sphere_and_apex_finite_traces_keep_original_geometry() {
    let input = text();
    let source = Source::parse(&input, 1e-7).unwrap();
    let pairs = pairs(&source);
    assert_eq!(pairs.len(), 4);
    let mut budget = Budget::new(MAX_WORK, None);
    for p in pairs {
        let proof = certify(
            &source,
            p.curve,
            &p.definition,
            &p.surface,
            &p.edge,
            1e-7,
            &mut budget,
        )
        .unwrap()
        .unwrap();
        assert!(matches!(proof,Proof::Certified(error) if (0. ..=1e-7).contains(&error)));
    }
}
#[test]
fn actual_spherical_curve_wrong_half_and_overlong_arc_are_not_full_circle_admissions() {
    let input = text();
    let source = Source::parse(&input, 1e-7).unwrap();
    let pair = pairs(&source)
        .into_iter()
        .find(|p| matches!(p.surface, Surface::Sphere { .. }))
        .unwrap();
    // Same ACTUAL curve operator with different finite subarcs. These are
    // mathematical domain controls, not forged STEP/topology certificates.
    for range in [
        [std::f64::consts::PI, std::f64::consts::PI + 0.2],
        [0., 2. * std::f64::consts::TAU],
    ] {
        let edge = Edge {
            range,
            ..pair.edge.clone()
        };
        let mut budget = Budget::new(MAX_WORK, None);
        assert!(matches!(
            certify(
                &source,
                pair.curve,
                &pair.definition,
                &pair.surface,
                &edge,
                1e-7,
                &mut budget
            )
            .unwrap()
            .unwrap(),
            Proof::Unsupported
        ));
    }
}
fn change_polyline_v(input: &str, source: &Source<'_>, point: u32, v: &str) -> String {
    let old = source.simple(point).unwrap().source.unwrap();
    let u = coordinate(source, point, 0).unwrap();
    let Value::Number(u) = u else { panic!() };
    input.replacen(old, &format!("CARTESIAN_POINT('',({u:.17E},{v}))"), 1)
}
#[test]
fn actual_undersized_cone_trace_misses_material_edge_and_is_refused() {
    let input = text();
    let source = Source::parse(&input, 1e-7).unwrap();
    let pair = pairs(&source)
        .into_iter()
        .find(|p| matches!(p.surface, Surface::Cone { .. }))
        .unwrap();
    let ids = polyline(
        &source,
        pair.definition.curve_id,
        &mut Budget::new(MAX_WORK, None),
    )
    .unwrap()
    .unwrap();
    let changed = change_polyline_v(&input, &source, ids[0], "0.1");
    let source = Source::parse(&changed, 1e-7).unwrap();
    let pair = pairs(&source)
        .into_iter()
        .find(|p| matches!(p.surface, Surface::Cone { .. }))
        .unwrap();
    let mut budget = Budget::new(MAX_WORK, None);
    assert!(
        matches!(certify(&source,pair.curve,&pair.definition,&pair.surface,&pair.edge,1e-7,&mut budget).unwrap().unwrap(),Proof::Disagreement(lower) if lower>0.5)
    );
}
#[test]
fn entire_collapsed_row_not_just_equivalent_point_is_required() {
    let input = text();
    let source = Source::parse(&input, 1e-7).unwrap();
    let (&id, parts) = source
        .entities
        .iter()
        .find(|(_, parts)| {
            parts
                .iter()
                .any(|p| p.name == "EVALUATED_DEGENERATE_PCURVE")
        })
        .unwrap();
    let entity = &parts[0];
    let [_, _, representation, equivalent] = entity.args.as_slice() else {
        panic!()
    };
    let representation = reference(representation).unwrap();
    let equivalent = reference(equivalent).unwrap();
    let mut budget = Budget::new(MAX_WORK, None);
    assert_eq!(
        equivalent_vertex(&source, id, &mut budget)
            .unwrap()
            .unwrap(),
        source.point(equivalent).unwrap()
    );
    let [_, items, _] = source
        .args(representation, "DEFINITIONAL_REPRESENTATION")
        .unwrap()
    else {
        panic!()
    };
    let [curve] = list(items).unwrap() else {
        panic!()
    };
    let points = polyline(&source, reference(curve).unwrap(), &mut budget)
        .unwrap()
        .unwrap();
    let wrong = change_polyline_v(&input, &source, points[1], "0.01");
    let bad = Source::parse(&wrong, 1e-7).unwrap();
    assert!(equivalent_vertex(&bad, id, &mut Budget::new(MAX_WORK, None)).is_err());
    let old = entity.source.unwrap();
    let new = old.replace(&format!(",#{equivalent})"), ",#900001)");
    assert_ne!(old, new);
    let shifted = input.replacen(old, &new, 1);
    let at = shifted.rfind("ENDSEC;").unwrap();
    let shifted = format!(
        "{}#900001=CARTESIAN_POINT('',(-19.,0.,0.));\n{}",
        &shifted[..at],
        &shifted[at..]
    );
    let bad = Source::parse(&shifted, 1e-7).unwrap();
    assert!(equivalent_vertex(&bad, id, &mut Budget::new(MAX_WORK, None)).is_err());
}
#[test]
fn finite_trace_work_frontier_cancel_and_retry_preserve_original_source() {
    let input = text();
    let source = Source::parse(&input, 1e-7).unwrap();
    let pair = pairs(&source)
        .into_iter()
        .find(|p| matches!(p.surface, Surface::Sphere { .. }))
        .unwrap();
    let initial = MAX_WORK;
    let mut first = Budget::new(initial, None);
    assert!(matches!(
        certify(
            &source,
            pair.curve,
            &pair.definition,
            &pair.surface,
            &pair.edge,
            1e-7,
            &mut first
        )
        .unwrap()
        .unwrap(),
        Proof::Certified(_)
    ));
    let used = initial - first.remaining;
    let mut exact = Budget::new(used, None);
    assert!(matches!(
        certify(
            &source,
            pair.curve,
            &pair.definition,
            &pair.surface,
            &pair.edge,
            1e-7,
            &mut exact
        )
        .unwrap()
        .unwrap(),
        Proof::Certified(_)
    ));
    let mut small = Budget::new(used - 1, None);
    assert!(matches!(
        certify(
            &source,
            pair.curve,
            &pair.definition,
            &pair.surface,
            &pair.edge,
            1e-7,
            &mut small
        ),
        Err(ReadError::Limit)
    ));
    let count = std::cell::Cell::new(0);
    let stop = || {
        count.set(count.get() + 1);
        count.get() == 20
    };
    let mut cancelled = Budget::new(initial, Some(&stop));
    assert!(
        certify(
            &source,
            pair.curve,
            &pair.definition,
            &pair.surface,
            &pair.edge,
            1e-7,
            &mut cancelled
        )
        .is_err()
    );
    assert!(cancelled.cancelled.get());
    let mut retry = Budget::new(initial, None);
    assert!(matches!(
        certify(
            &source,
            pair.curve,
            &pair.definition,
            &pair.surface,
            &pair.edge,
            1e-7,
            &mut retry
        )
        .unwrap()
        .unwrap(),
        Proof::Certified(_)
    ));
    assert_eq!(first.remaining, retry.remaining);
}
#[test]
fn binary_tau_used_as_legal_bound_is_below_true_tau_by_independent_exact_series() {
    use num_bigint::BigInt;
    fn sum(d: u32, n: usize) -> Q {
        (0..n).fold(Q::from_integer(0.into()), |s, k| {
            let term = Q::new(
                1.into(),
                BigInt::from(2 * k + 1) * BigInt::from(d).pow((2 * k + 1) as u32),
            );
            if k % 2 == 0 { s + term } else { s - term }
        })
    }
    // Even alternating partial sums are lower bounds; odd ones are upper.
    // Machin: pi = 16 atan(1/5) - 4 atan(1/239).
    let lower = Q::from_integer(2.into())
        * (Q::from_integer(16.into()) * sum(5, 12) - Q::from_integer(4.into()) * sum(239, 5));
    assert!(Q::from_float(std::f64::consts::TAU).unwrap() < lower);
}
