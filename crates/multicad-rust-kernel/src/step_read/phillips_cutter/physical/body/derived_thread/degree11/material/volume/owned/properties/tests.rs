use super::*;
use num_rational::BigRational as Q;
fn contains(actual: I, n: i64, d: i64) {
    let expected = Q::new(n.into(), d.into());
    assert!(
        Q::from_float(actual.lo).unwrap() <= expected
            && expected <= Q::from_float(actual.hi).unwrap(),
        "{actual:?} misses {expected}"
    );
}
#[test]
fn exact_cap_green_integrals_match_square_and_curved_polynomial_wedge() {
    let weights = planar::Weights::new();
    let mut square = [I::point(0.); 3];
    for side in 0..4 {
        let p = std::array::from_fn(|i| {
            let i = i as f64;
            let p = match side {
                0 => [i, 0.],
                1 => [11., i],
                2 => [11. - i, 11.],
                _ => [0., 11. - i],
            };
            p.map(I::point)
        });
        let value = weights.curve(p);
        square = std::array::from_fn(|i| square[i] + value[i]);
    }
    contains(square[0], 121, 1);
    contains(square[1], 1331, 2);
    contains(square[2], 1331, 2);
    let curve = weights.curve(std::array::from_fn(|i| {
        [I::point((i * i) as f64), I::point(i as f64)]
    }));
    let closing = weights.curve(std::array::from_fn(|i| {
        [I::point((11 * (11 - i)) as f64), I::point((11 - i) as f64)]
    }));
    let wedge = std::array::from_fn::<_, 3, _>(|i| curve[i] + closing[i]);
    contains(wedge[0], -605, 3);
    contains(wedge[1], -19965, 2);
    contains(wedge[2], -6655, 6);
}
#[test]
fn actual_owned_d11_end_disks_survive_owner_drop_and_resource_refusals() {
    for left in [false, true] {
        let source = super::super::query::tests::actual_owned("end-disks", left, left, 1);
        let o = Degree11EndDiskOptions::default();
        let good = source.end_disk_properties(o, None).unwrap();
        assert!(std::ptr::eq(good.source(), &source));
        for (end, value) in good.ends().iter().enumerate() {
            assert_eq!(value.station_mm, source.cap_stations_mm()[end]);
            assert!(
                value.area_mm2[0] > 0.
                    && value.area_mm2[1] - value.area_mm2[0] <= o.maximum_area_width_mm2
            );
            assert_eq!(value.centroid_mm[2], [value.station_mm; 2]);
            // A single-start helical section need not have an axial centroid.
            // Its disk is bounded by the actual cap-control convex hull.
            for axis in 0..2 {
                let bounds = source
                    .thread_faces()
                    .iter()
                    .flat_map(|f| {
                        let height = f.source().dimensions()[1];
                        f.source()
                            .poles()
                            .chunks_exact(height)
                            .map(move |r| r[if end == 0 { 0 } else { height - 1 }][axis])
                    })
                    .fold([f64::INFINITY, f64::NEG_INFINITY], |b, x| {
                        [b[0].min(x), b[1].max(x)]
                    });
                assert!(
                    value.centroid_mm[axis][0] >= bounds[0] - 1e-7
                        && value.centroid_mm[axis][1] <= bounds[1] + 1e-7
                );
            }
        }
        for bad in [
            Degree11EndDiskOptions {
                maximum_work: good.work_used() - 1,
                ..o
            },
            Degree11EndDiskOptions {
                maximum_bytes: good.accounted_peak_bytes() - 1,
                ..o
            },
        ] {
            assert!(matches!(
                source.end_disk_properties(bad, None),
                Err(PhillipsReadError::Limit)
            ));
        }
        let polls = Cell::new(0usize);
        source
            .end_disk_properties(
                o,
                Some(&|| {
                    polls.set(polls.get() + 1);
                    false
                }),
            )
            .unwrap();
        for stop in [1, polls.get() / 2, polls.get()] {
            let count = Cell::new(0usize);
            assert!(matches!(
                source.end_disk_properties(
                    o,
                    Some(&|| {
                        count.set(count.get() + 1);
                        count.get() >= stop
                    })
                ),
                Err(PhillipsReadError::Cancelled)
            ));
        }
        assert!(matches!(
            source.end_disk_properties(
                Degree11EndDiskOptions {
                    maximum_work: END_WORK + 1,
                    ..o
                },
                None
            ),
            Err(PhillipsReadError::Invalid)
        ));
        assert!(matches!(
            source.end_disk_properties(
                Degree11EndDiskOptions {
                    maximum_area_width_mm2: f64::NAN,
                    ..o
                },
                None
            ),
            Err(PhillipsReadError::Invalid)
        ));
        assert!(matches!(
            source.end_disk_properties(
                Degree11EndDiskOptions {
                    maximum_area_width_mm2: 1e-30,
                    ..o
                },
                None
            ),
            Err(PhillipsReadError::Unsupported)
        ));
        assert_eq!(
            source.end_disk_properties(o, None).unwrap().ends()[1].area_mm2,
            good.ends()[1].area_mm2
        );
    }
}
