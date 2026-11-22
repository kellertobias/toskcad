use super::*;
use num_bigint::BigInt;
use num_rational::BigRational as Q;
use num_traits::{Signed, Zero};
use std::{cell::Cell, sync::Arc};
fn q(x: f64) -> Q {
    Q::from_float(x).unwrap()
}
fn decimal(s: &str) -> Q {
    let (m, e) = s.split_once('E').unwrap_or((s, "0"));
    let places = m.split_once('.').map_or(0, |(_, v)| v.len()) as i32;
    let n = m.replace('.', "").parse::<BigInt>().unwrap();
    let power = e.parse::<i32>().unwrap() - places;
    if power >= 0 {
        Q::from_integer(n * BigInt::from(10).pow(power as u32))
    } else {
        Q::new(n, BigInt::from(10).pow((-power) as u32))
    }
}
fn choose(n: usize, k: usize) -> usize {
    (0..k).fold(1, |p, i| p * (n - i) / (i + 1))
}
fn product(a: &[Q], b: &[Q]) -> Vec<Q> {
    let n = a.len() - 1;
    let m = b.len() - 1;
    (0..=n + m)
        .map(|k| {
            (0..=n)
                .filter(|&i| k >= i && k - i <= m)
                .fold(Q::zero(), |sum, i| {
                    sum + &a[i]
                        * &b[k - i]
                        * Q::new(
                            BigInt::from(choose(n, i) * choose(m, k - i)),
                            BigInt::from(choose(n + m, k)),
                        )
                })
        })
        .collect()
}
fn integrate(a: &[Q], b: &[Q]) -> Q {
    let p = product(a, b);
    p.iter().sum::<Q>() / Q::from_integer(BigInt::from(p.len()))
}
/// Independent exact polynomial Green integrals of ORIGINAL decimal leaves;
/// no production CapWeights or reported measure fields authorize this oracle.
fn cap(source: &OwnedActualHead, index: usize) -> (Q, [Q; 2]) {
    let aux = source.auxiliary_thread();
    let owner = &aux.faces()[index];
    let mut area = Q::zero();
    let mut moment = [Q::zero(), Q::zero()];
    for co in owner.coedges() {
        let curve = &aux.curves()[aux.edges()[co.edge_index()].curve_index()];
        for span in 0..curve.spans() {
            let p: Vec<[Q; 2]> = (0..12)
                .map(|i| {
                    aux.point(curve.point_indices()[11 * span + i])
                        .unwrap()
                        .coordinates()
                        .map(|v| decimal(v.spelling()))
                })
                .map(|p| [p[0].clone(), p[1].clone()])
                .collect();
            let x: Vec<Q> = p.iter().map(|v| v[0].clone()).collect();
            let y: Vec<Q> = p.iter().map(|v| v[1].clone()).collect();
            let dx: Vec<Q> = x.windows(2).map(|p| (&p[1] - &p[0]) * q(11.)).collect();
            let dy: Vec<Q> = y.windows(2).map(|p| (&p[1] - &p[0]) * q(11.)).collect();
            let a = (integrate(&x, &dy) - integrate(&y, &dx)) / q(2.);
            let mx = integrate(&product(&x, &x), &dy) / q(2.);
            let my = -integrate(&product(&y, &y), &dx) / q(2.);
            if co.forward() {
                area += a;
                moment[0] += &mx;
                moment[1] += &my;
            } else {
                area -= a;
                moment[0] -= &mx;
                moment[1] -= &my;
            }
        }
    }
    if area.is_negative() {
        area = -area;
        moment = moment.map(|v| -v);
    }
    (area, moment)
}
fn enclosed(bounds: [f64; 2], value: &Q) {
    assert!(
        q(bounds[0]) <= *value && *value <= q(bounds[1]),
        "{bounds:?} misses {value}"
    );
}
#[test]
fn selected_actual_planar_faces_enclose_independent_original_green_and_literal_polygon_integrals() {
    for quarter in [false, true] {
        let source = OwnedActualHead::read(
            Arc::from(super::super::source_tests::candidate(quarter, 18.)),
            Default::default(),
            None,
        )
        .unwrap();
        let polygon = source.head_polygon_mm();
        let mut twice = Q::zero();
        let mut outer_m = [Q::zero(), Q::zero()];
        for i in 0..6 {
            let a = polygon[i].map(q);
            let b = polygon[(i + 1) % 6].map(q);
            let c = &a[0] * &b[1] - &a[1] * &b[0];
            twice += &c;
            for k in 0..2 {
                outer_m[k] += (&a[k] + &b[k]) * &c / q(6.);
            }
        }
        let outer = twice / q(2.);
        let [low, high] = source.head_stations_mm();
        let roles = source.physical_roles();
        let mut planar = 0;
        let upper = source
            .auxiliary_thread()
            .faces()
            .iter()
            .position(|f| f.entity() == source.head().auxiliary_cap_entity)
            .unwrap();
        let (cap_area, cap_m) = cap(&source, upper);
        for ordinal in 0..17 {
            let value = source.planar_face_measure(ordinal, Default::default(), None);
            if ordinal < 9 {
                let entity = roles.thread_face_entities[ordinal];
                let index = source
                    .auxiliary_thread()
                    .faces()
                    .iter()
                    .position(|f| f.entity() == entity)
                    .unwrap();
                if matches!(
                    source.auxiliary_thread().faces()[index].support(),
                    Support::Tensor { .. }
                ) {
                    assert!(matches!(value, Err(Error::Unsupported)));
                    continue;
                }
                let value = value.unwrap();
                let (area, m) = cap(&source, index);
                enclosed(value.area_mm2(), &area);
                for (k, mean) in m.iter().enumerate() {
                    enclosed(value.centroid_mm()[k], &(mean / &area));
                }
            } else {
                let value = value.unwrap();
                assert!(std::ptr::eq(value.source(), &source));
                assert_eq!(value.physical_face_ordinal(), ordinal);
                assert_eq!(
                    value.physical_face_entity(),
                    roles.head_face_entities[ordinal - 9]
                );
                let face = &source.head().faces[ordinal - 9];
                if face.outer.len() == 6 {
                    let (area, m, z) = if face.inner.is_empty() {
                        (outer.clone(), outer_m.clone(), high)
                    } else {
                        (
                            &outer - &cap_area,
                            std::array::from_fn(|k| &outer_m[k] - &cap_m[k]),
                            low,
                        )
                    };
                    enclosed(value.area_mm2(), &area);
                    for (k, m) in m.iter().enumerate() {
                        enclosed(value.centroid_mm()[k], &(m / &area));
                    }
                    enclosed(value.centroid_mm()[2], &q(z));
                } else {
                    let ring: Vec<_> = face
                        .outer
                        .iter()
                        .map(|&(e, forward)| {
                            source.head().vertices[source.head().edges[e][usize::from(!forward)]]
                        })
                        .collect();
                    let mut cross = [Q::zero(), Q::zero(), Q::zero()];
                    for i in 1..ring.len() - 1 {
                        let a: [Q; 3] = std::array::from_fn(|k| q(ring[i][k]) - q(ring[0][k]));
                        let b: [Q; 3] = std::array::from_fn(|k| q(ring[i + 1][k]) - q(ring[0][k]));
                        for k in 0..3 {
                            cross[k] += &a[(k + 1) % 3] * &b[(k + 2) % 3]
                                - &a[(k + 2) % 3] * &b[(k + 1) % 3];
                        }
                    }
                    let square = cross.iter().map(|v| v * v).sum::<Q>() / q(4.);
                    let bounds = value.area_mm2();
                    assert!(
                        q(bounds[0]) * q(bounds[0]) <= square
                            && square <= q(bounds[1]) * q(bounds[1])
                    );
                    for k in 0..3 {
                        let mean = ring.iter().map(|p| q(p[k])).sum::<Q>() / q(4.);
                        enclosed(value.centroid_mm()[k], &mean);
                    }
                }
            }
            planar += 1;
        }
        assert_eq!(planar, 9);
    }
}
#[test]
fn actual_selected_planar_measure_frontiers_cancel_precision_and_retry() {
    let source = OwnedActualHead::read(
        Arc::from(super::super::source_tests::candidate(false, 18.)),
        Default::default(),
        None,
    )
    .unwrap();
    for ordinal in 8..17 {
        let good = source
            .planar_face_measure(ordinal, Default::default(), None)
            .unwrap();
        let work = good.work_used();
        let bytes = good.accounted_peak_bytes();
        let area = good.area_mm2();
        for (maximum_work, maximum_bytes, ok) in [
            (work - 1, MAX_BYTES, false),
            (work, MAX_BYTES, true),
            (MEASURE_WORK, bytes - 1, false),
            (MEASURE_WORK, bytes, true),
        ] {
            let result = source.planar_face_measure(
                ordinal,
                HeadMeasureOptions {
                    measure: FaceMeasureOptions {
                        maximum_work,
                        maximum_bytes,
                        ..Default::default()
                    },
                    ..Default::default()
                },
                None,
            );
            assert_eq!(result.is_ok(), ok);
            if !ok {
                assert!(matches!(result, Err(Error::Limit)));
            }
        }
        let n = Cell::new(0);
        let observe = || {
            n.set(n.get() + 1);
            false
        };
        source
            .planar_face_measure(ordinal, Default::default(), Some(&observe))
            .unwrap();
        let at = n.get();
        let seen = Cell::new(0);
        let cancel = || {
            seen.set(seen.get() + 1);
            seen.get() == at
        };
        assert!(matches!(
            source.planar_face_measure(ordinal, Default::default(), Some(&cancel)),
            Err(Error::Cancelled)
        ));
        assert_eq!(
            source
                .planar_face_measure(ordinal, Default::default(), None)
                .unwrap()
                .area_mm2(),
            area
        );
        assert_eq!(
            source
                .planar_face_measure(
                    ordinal,
                    HeadMeasureOptions {
                        additional_resident_bytes: 4096,
                        ..Default::default()
                    },
                    None
                )
                .unwrap()
                .accounted_peak_bytes(),
            bytes + 4096
        );
    }
    assert!(matches!(
        source.planar_face_measure(17, Default::default(), None),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        source.planar_face_measure(
            8,
            HeadMeasureOptions {
                measure: FaceMeasureOptions {
                    area_width_mm2: f64::from_bits(1),
                    ..Default::default()
                },
                ..Default::default()
            },
            None
        ),
        Err(Error::Unsupported)
    ));
}
