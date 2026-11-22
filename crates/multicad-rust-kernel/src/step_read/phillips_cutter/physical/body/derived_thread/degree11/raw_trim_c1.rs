//! Real saved M3 data/decimal correspondence experiment; no file admission.
use super::*;
use crate::phillips_drive::{
    CurveKey, TRIM_SPAN_COMPARE_WORK, TRIM_SPAN_SOURCE_WORK, TrimSpanCorrespondence,
};
use crate::step_read::phillips_cutter::physical::{
    Curve, body::common_mouth::tests::derived_thread_probe::with_head,
};
use std::collections::BTreeMap;
fn interval([lo, hi]: [f64; 2]) -> I {
    I::new(lo, hi)
}
#[test]
fn actual_m3_raw_trim_c1_uses_original_decimal_points_knots_and_all_corner_orbits() {
    for left in [false, true] {
        let text = tests::actual_m3(true, left);
        with_head(&text, &mut |head| {
            let representative = head.common_mouth().closed_cavity().representative();
            let body = head.body();
            let events = representative
                .derived_cavity_source()
                .boundary_events()
                .unwrap();
            let [origin, u, v, n] = representative.derived_affine_columns();
            let mut prepared: Vec<TrimSpanCorrespondence<'_, '_>> = Vec::new();
            let mut keys: BTreeMap<(usize, [[u64; 2]; 2]), usize> = BTreeMap::new();
            let mut points: BTreeMap<u32, [[f64; 2]; 3]> = BTreeMap::new();
            let mut work = 0usize;
            let mut spans = 0usize;
            let mut curves = 0usize;
            let mut position = 0_f64;
            let mut derivative = 0_f64;
            let mut whole_trim_derivative = 0_f64;
            let mut worst = None;
            for mapping in representative.edges() {
                let CurveKey::Trim {
                    envelope,
                    wall,
                    quarter,
                    lower,
                } = mapping.source_role
                else {
                    continue;
                };
                let edge = body
                    .cavity()
                    .edges()
                    .iter()
                    .find(|e| e.geometry == mapping.raw_curve_entity)
                    .unwrap();
                let Curve::CubicIntersection {
                    spline_entity,
                    pole_entities,
                    spans: raw,
                    ..
                } = &edge.curve
                else {
                    continue;
                };
                assert!(mapping.raw_curve_forward ^ mapping.raw_edge_reversed);
                let index = events
                    .trims()
                    .iter()
                    .position(|t| t.envelope_role() == envelope && t.wall_kind() == wall)
                    .unwrap();
                let axes = match quarter {
                    0 => [u, v, n],
                    1 => [v, u.map(|x| -x), n],
                    2 => [u.map(|x| -x), v.map(|x| -x), n],
                    3 => [v.map(|x| -x), u, n],
                    _ => unreachable!(),
                };
                let axes = if lower {
                    [axes[0], axes[1].map(|x| -x), axes[2]]
                } else {
                    axes
                };
                let numbers = body
                    .original_numeric_leaves(*spline_entity, 16_000_000 - work, None)
                    .unwrap();
                assert!(std::ptr::eq(numbers.body(), body));
                work += numbers.work_used();
                let count = raw.len() + 1;
                assert_eq!(numbers.leaves().len(), 1 + 2 * count);
                let knots = &numbers.leaves()[1 + count..];
                let first = interval(knots[0].bounds());
                let last = interval(knots[count - 1].bounds());
                let length = last - first;
                assert!(length.lo > 0.);
                let normalized = |i: usize| {
                    // SAME original scalar expression identities, not cached
                    // float==0/1 as a claim about the original decimal token.
                    if i == 0 {
                        I::point(0.)
                    } else if i + 1 == count {
                        I::point(1.)
                    } else {
                        (interval(knots[i].bounds()) - first) / length
                    }
                };
                curves += 1;
                for s in 0..raw.len() {
                    let parameters = [normalized(s), normalized(s + 1)].map(|x| [x.lo, x.hi]);
                    let key = (index, parameters.map(|p| p.map(f64::to_bits)));
                    // <=4096 canonical cache entries, height/11key node search,
                    // key comparison+index/check bookkeeping fits512/lookup.
                    work += 512;
                    assert!(
                        work <= 16_000_000,
                        "fresh scope work{work} spans{spans} curves{curves} canonical{} c0{position} c1{derivative}",
                        prepared.len()
                    );
                    let cache = if let Some(&i) = keys.get(&key) {
                        i
                    } else {
                        assert!(prepared.len() < 4096);
                        work += TRIM_SPAN_SOURCE_WORK;
                        assert!(
                            work <= 16_000_000,
                            "fresh scope work{work} spans{spans} curves{curves} canonical{} c0{position} c1{derivative}",
                            prepared.len()
                        );
                        let value = events
                            .prepare_trim_span_correspondence(
                                index,
                                parameters,
                                TRIM_SPAN_SOURCE_WORK,
                                64 * 1024 * 1024,
                                None,
                            )
                            .unwrap();
                        let i = prepared.len();
                        prepared.push(value);
                        keys.insert(key, i);
                        i
                    };
                    let mut controls = [[[0.; 2]; 3]; 4];
                    for k in 0..4 {
                        work += 512;
                        assert!(
                            work <= 16_000_000,
                            "fresh scope work{work} spans{spans} curves{curves} canonical{} c0{position} c1{derivative}",
                            prepared.len()
                        );
                        let id = pole_entities[3 * s + k];
                        let p = if let Some(p) = points.get(&id) {
                            *p
                        } else {
                            let point = body.exact_point_leaf(id, 16_000_000 - work, None).unwrap();
                            assert!(std::ptr::eq(point.body(), body));
                            assert_eq!(point.entity(), id);
                            work += point.work_used();
                            assert!(
                                work <= 16_000_000,
                                "fresh scope work{work} spans{spans} curves{curves} canonical{} c0{position} c1{derivative}",
                                prepared.len()
                            );
                            let p = point.coordinate_bounds_mm();
                            points.insert(id, p);
                            p
                        };
                        controls[k] = p;
                    }
                    work += TRIM_SPAN_COMPARE_WORK;
                    assert!(
                        work <= 16_000_000,
                        "fresh scope work{work} spans{spans} curves{curves} canonical{} c0{position} c1{derivative}",
                        prepared.len()
                    );
                    let compared = prepared[cache]
                        .compare_controls(
                            &controls,
                            origin,
                            axes,
                            TRIM_SPAN_COMPARE_WORK,
                            64 * 1024 * 1024,
                            None,
                        )
                        .unwrap();
                    let p = compared.maximum_position_mm();
                    let d = compared.maximum_derivative_mm();
                    let width = interval(parameters[1]) - interval(parameters[0]);
                    assert!(width.lo > 0.);
                    whole_trim_derivative =
                        whole_trim_derivative.max((I::point(d) / I::point(width.lo)).hi);
                    position = position.max(p);
                    if d > derivative {
                        derivative = d;
                        worst = Some((
                            mapping.raw_curve_entity,
                            s,
                            index,
                            quarter,
                            lower,
                            parameters,
                        ));
                    }
                    spans += 1;
                }
            }
            eprintln!(
                "RAW_TRIM_C1 left={left} curves={curves} spans={spans} canonical={} original_points={} fresh_work={work} inherited_head_work={} maximum_c0={position:.17E} maximum_c1={derivative:.17E} whole_trim_c1={whole_trim_derivative:.17E} worst={worst:?}",
                prepared.len(),
                points.len(),
                head.work_used()
            );
            assert!(curves > 0 && spans > 0 && prepared.len() < spans);
            assert!(
                whole_trim_derivative <= 1e-7,
                "actual raw whole-trim C1 {whole_trim_derivative}"
            );
            assert!(position <= 1e-7, "actual raw C0 policy {position}");
            assert!(
                derivative <= 1e-7,
                "actual raw ORIGINAL-span C1 policy {derivative}"
            );
        });
    }
}
