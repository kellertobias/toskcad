use super::*;
use crate::{authored_flat_thread::SourceOptions, helical_thread::Spec};
use num_bigint::{BigInt, Sign};
use num_rational::BigRational;
fn number(v: ExactBoundaryRational<'_>) -> BigRational {
    let sign = match v.signum() {
        -1 => Sign::Minus,
        0 => Sign::NoSign,
        _ => Sign::Plus,
    };
    BigRational::new(
        BigInt::from_bytes_le(sign, v.numerator_magnitude_le()),
        BigInt::from_bytes_le(Sign::Plus, v.denominator_magnitude_le()),
    )
}
fn thread(left: bool, starts: u32) -> AuthoredThread {
    AuthoredThread::new(
        Spec {
            major_radius_mm: 5.,
            pitch_mm: 1.9,
            length_mm: 4.7,
            depth_mm: 1.,
            crest_width_mm: 0.4,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            starts,
            left_hand: left,
            root_radius_mm: None,
        },
        SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-7,
        },
    )
    .unwrap()
}
#[test]
fn original_rational_cells_restore_periodic_aliases_for_both_hands() {
    for (left, starts) in [(false, 1), (true, 2)] {
        let original = thread(left, starts);
        let body = original
            .boundary_with_interrupt(BoundaryOptions::default(), None)
            .unwrap();
        let mut zero = false;
        let mut one = false;
        for (index, cell) in body.cells().iter().enumerate() {
            let values = body
                .exact_cell_parameters(index, 16_000_000, 0, None)
                .unwrap();
            assert!(std::ptr::eq(values.boundary(), &body));
            assert_eq!(values.domain(), ExactBoundaryDomain::Cell(index));
            assert_eq!(values.row_count(), cell.count);
            assert_eq!(values.column_count(), 2);
            assert!(values.owned_bytes() < 64 * 1024);
            assert!(values.work_used() > 0);
            let mut w = Work {
                remaining: 16_000_000,
                used: 0,
                interrupt: None,
            };
            let patch = Patch::from_chart(
                body.graph.source.chart(),
                cell.period,
                cell.band,
                1,
                false,
                &mut w,
            )
            .unwrap();
            let polygon = patch.cell(0, 0, 1, Kind::Wall, &mut w).unwrap();
            assert_eq!(polygon.points.len(), values.row_count());
            for (i, p) in polygon.points.iter().enumerate() {
                let q = number(values.coordinate(i, 0).unwrap());
                let v = number(values.coordinate(i, 1).unwrap());
                assert_eq!(q, p.q.rational(&mut w).unwrap());
                assert_eq!(v, p.v.rational(&mut w).unwrap());
                zero |= v == BigRational::from_integer(0.into());
                one |= v == BigRational::from_integer(1.into());
            }
            assert!(values.coordinate(0, 2).is_none());
            assert!(values.coordinate(usize::MAX, 0).is_none());
        }
        assert!(zero && one);
        for (index, curve) in body.source_curves() {
            let values = body
                .exact_curve_coefficients(index, 16_000_000, 0, None)
                .unwrap();
            assert_eq!(values.domain(), ExactBoundaryDomain::Curve(index));
            assert_eq!(values.row_count(), 2);
            assert_eq!(values.column_count(), 3);
            let mut w = Work {
                remaining: 16_000_000,
                used: 0,
                interrupt: None,
            };
            for row in 0..2 {
                for column in 0..3 {
                    let expected = close::scalar(&curve.exact, row * 3 + column, &mut w)
                        .unwrap()
                        .rational(&mut w)
                        .unwrap();
                    assert_eq!(number(values.coordinate(row, column).unwrap()), expected);
                }
            }
        }
    }
}
#[test]
fn exact_query_limits_cancellation_invalid_keys_and_unsupported_cones_fail_closed() {
    let body = thread(false, 1)
        .boundary_with_interrupt(BoundaryOptions::default(), None)
        .unwrap();
    assert!(matches!(
        body.exact_cell_parameters(usize::MAX, 16_000_000, 0, None),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        body.exact_curve_coefficients(u32::MAX, 16_000_000, 0, None),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        body.exact_cell_parameters(0, 16_000_001, 0, None),
        Err(Error::Invalid)
    ));
    for reserved in [usize::MAX, 64 * 1024 * 1024] {
        assert!(matches!(
            body.exact_cell_parameters(0, 16_000_000, reserved, None),
            Err(Error::Limit)
        ));
        assert!(matches!(
            body.exact_curve_coefficients(0, 16_000_000, reserved, None),
            Err(Error::Limit)
        ));
    }
    assert!(matches!(
        body.exact_cell_parameters(0, 1, 0, None),
        Err(Error::Limit)
    ));
    assert!(matches!(
        body.exact_cell_parameters(0, 16_000_000, 0, Some(&|| true)),
        Err(Error::Cancelled)
    ));
    let polls = std::cell::Cell::new(0);
    let stop = || {
        polls.set(polls.get() + 1);
        polls.get() > 10
    };
    assert!(matches!(
        body.exact_curve_coefficients(0, 16_000_000, 0, Some(&stop)),
        Err(Error::Cancelled)
    ));
    let source =
        crate::authored_internal_treatment::mesh::tests::source(false, 1, Some(0.4), false);
    let treated = source
        .boundary_with_interrupt(BoundaryOptions::default(), None)
        .unwrap();
    let wall = treated
        .cells()
        .iter()
        .position(|c| !c.is_entry_cone())
        .unwrap();
    let wall_parameters = treated
        .exact_cell_parameters(wall, 16_000_000, 0, None)
        .unwrap();
    assert_eq!(
        wall_parameters.row_count(),
        treated.cells()[wall].points().len()
    );
    assert_eq!(wall_parameters.domain(), ExactBoundaryDomain::Cell(wall));
    let mut w = Work {
        remaining: 16_000_000,
        used: 0,
        interrupt: None,
    };
    let cell = &treated.cells()[wall];
    let patch = Patch::from_chart(
        treated.graph.source.chart(),
        cell.period,
        cell.band,
        1,
        false,
        &mut w,
    )
    .unwrap();
    let polygon = patch.cell(0, 0, 1, Kind::Wall, &mut w).unwrap();
    for (i, p) in polygon.points.iter().enumerate() {
        assert_eq!(
            number(wall_parameters.coordinate(i, 0).unwrap()),
            p.q.rational(&mut w).unwrap()
        );
        assert_eq!(
            number(wall_parameters.coordinate(i, 1).unwrap()),
            p.v.rational(&mut w).unwrap()
        );
    }
    for small in [5, 6] {
        assert!(matches!(
            body.exact_cell_parameters(0, small, 0, None),
            Err(Error::Limit)
        ));
        assert!(matches!(
            treated.exact_cell_parameters(wall, small, 0, None),
            Err(Error::Limit)
        ));
    }
    let cone = treated
        .cells()
        .iter()
        .position(|c| c.is_entry_cone())
        .unwrap();
    assert!(matches!(
        treated.exact_cell_parameters(cone, 16_000_000, 0, None),
        Err(Error::Unresolved)
    ));
}
#[test]
fn an_exact_corner_below_one_is_not_snapped_when_its_enclosure_reaches_one() {
    let original = thread(false, 7);
    let source = AuthoredThread::new(
        Spec {
            length_mm: 1.9 * 7.,
            ..original.spec()
        },
        SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-7,
        },
    )
    .unwrap();
    let body = source
        .boundary_with_interrupt(BoundaryOptions::default(), None)
        .unwrap();
    let mut found = false;
    for (index, cell) in body.cells().iter().enumerate() {
        let data = body
            .exact_cell_parameters(index, 16_000_000, 0, None)
            .unwrap();
        for corner in 0..cell.count {
            let v = cell.parameters[corner][1];
            if v.lo < 1. && v.hi == 1. {
                let exact = number(data.coordinate(corner, 1).unwrap());
                assert!(exact < BigRational::from_integer(1.into()));
                found = true;
            }
        }
    }
    assert!(
        found,
        "actual7-start source must exercise a rational rounding up to one"
    );
}
