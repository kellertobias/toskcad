//! Exact coordinate-unit cost experiment, not a file/material admission.
use super::*;
use crate::authored_hanger_bolt::{BoundaryOptions, Hanger};
use num_bigint::BigInt;
use num_traits::{One, Zero};
fn dyadic(q: &N, w: &mut Work<'_>) -> Result<usize, Error> {
    let d = q.0.denom();
    w.charge(2 + (d.bits() as usize).div_ceil(64))?;
    assert!((d & (d - BigInt::one())).is_zero());
    Ok(d.bits() as usize - 1)
}
fn factor(field: &EncodedField<'_>, w: &mut Work<'_>) -> Result<N, Error> {
    let mut r = 0;
    let mut q = 0;
    for row in &field.meridians {
        for value in &row.radius {
            r = r.max(dyadic(&value.number(w)?, w)?);
        }
        for value in &row.phase {
            q = q.max(dyadic(&value.number(w)?, w)?);
        }
    }
    let mut angular = 0;
    for piece in field.angular.coefficients() {
        for row in piece {
            for value in row {
                angular = angular.max(dyadic(&N(value.clone()), w)?);
            }
        }
    }
    let mut knot = 0;
    for &value in field.angular.knots() {
        knot = knot.max(dyadic(&binary(value, w)?, w)?);
    }
    let source = field.source.source_boundary().source();
    let mut pitch = 0;
    for value in [
        source.metric().spec().pitch_mm,
        source.wood().spec().pitch_mm,
    ] {
        pitch = pitch.max(dyadic(&binary(value, w)?, w)?);
    }
    r = r.max(dyadic(&binary(source.spec().plain_radius_mm, w)?, w)?);
    for value in [
        source.metric().full_length_mm(),
        source.wood().spec().length_mm,
        source.spec().plain_length_mm,
    ] {
        q = q.max(dyadic(&binary(value, w)?, w)?);
    }
    let exponent = (r + angular).max(q).max(pitch + knot + 1);
    w.charge(2 + (exponent.div_ceil(64) + 1).pow(2))?;
    eprintln!("coordinate unit denominator=5*2^{exponent}");
    Ok(N(BigRational::from_integer(
        BigInt::from(5_u32) << exponent,
    )))
}
fn batch(field: &EncodedField<'_>, m: &mut Memory, w: &mut Work<'_>) -> Result<usize, Error> {
    let scale = factor(field, w)?;
    let mut count = 0;
    for meridian in 0..field.meridians.len() {
        let cell = field
            .cells
            .iter()
            .position(|c| c.meridian == meridian)
            .unwrap();
        let row = &field.meridians[meridian];
        let mut requests = vec![
            (None, 0, 0),
            (None, field.angular.coefficients().len() - 1, 10),
        ];
        for piece in 0..field.angular.coefficients().len() {
            requests.push((Some(0), piece, 0));
            requests.push((Some(row.degree()), piece, 0));
        }
        for (fixed, piece, column) in requests {
            let n = if fixed.is_some() {
                11
            } else {
                row.degree() + 1
            };
            for i in 0..n {
                let (u, v) = if let Some(u) = fixed {
                    (u, i)
                } else {
                    (i, column)
                };
                let (_, h) = field.homogeneous_control(cell, u, piece, v, w)?;
                for coordinate in &h[..3] {
                    let mm = coordinate.div(&h[3], w)?;
                    let integer = mm.mul(&scale, w)?;
                    assert!(integer.0.denom().is_one());
                    assert_eq!(integer.div(&scale, w)?.0, mm.0);
                    let before = m.used;
                    let text = super::super::decimal::decimal(&integer, m, w)?;
                    assert!(!text.as_str().contains('E'));
                    drop(text);
                    m.used = before;
                }
            }
            count += 1;
        }
    }
    w.charge(1)?;
    Ok(count)
}
#[test]
fn exact_scaled_coordinate_support_batch_uses_one_actual_ledger() {
    for (left, nonbinary) in [(false, false), (true, true)] {
        let source = Hanger::new(
            crate::authored_hanger_bolt::tests::fixture(left, nonbinary),
            None,
        )
        .unwrap();
        let boundary = source.boundary(BoundaryOptions::default(), None).unwrap();
        let field = EncodedField::prepare(&boundary, CurveOptions::default(), None).unwrap();
        let mut w = Work::new(16_000_000, None).unwrap();
        w.charge(field.work_used()).unwrap();
        let mut m = Memory::new(field.retained_bytes(), 64 * 1024 * 1024).unwrap();
        m.add(super::super::SCRATCH).unwrap();
        let result = batch(&field, &mut m, &mut w);
        eprintln!(
            "exact scaled support batch left={left} nonbinary={nonbinary} work={} result={result:?}",
            w.used
        );
        assert_eq!(result, Ok(180));
    }
}

#[test]
fn exact_coordinate_unit_batch_retains_real_work_byte_and_cancel_refusals() {
    let source = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, false),
        None,
    )
    .unwrap();
    let boundary = source.boundary(BoundaryOptions::default(), None).unwrap();
    let field = EncodedField::prepare(&boundary, CurveOptions::default(), None).unwrap();
    let bytes = field.retained_bytes() + super::super::SCRATCH;
    for maximum_work in [1, 64, 8192] {
        let mut w = Work::new(maximum_work, None).unwrap();
        let mut m = Memory::new(bytes, 64 * 1024 * 1024).unwrap();
        assert_eq!(batch(&field, &mut m, &mut w), Err(Error::Limit));
    }
    let mut w = Work::new(16_000_000, None).unwrap();
    let mut m = Memory::new(bytes, bytes).unwrap();
    assert_eq!(batch(&field, &mut m, &mut w), Err(Error::Limit));
    let cancel = || true;
    let mut w = Work::new(16_000_000, Some(&cancel)).unwrap();
    let mut m = Memory::new(bytes, 64 * 1024 * 1024).unwrap();
    assert_eq!(batch(&field, &mut m, &mut w), Err(Error::Cancelled));
}

fn cached_mm_batch(
    field: &EncodedField<'_>,
    m: &mut Memory,
    w: &mut Work<'_>,
) -> Result<usize, Error> {
    let mut cache = super::super::decimal::DecimalCache::new(m, w)?;
    let mut count = 0;
    for meridian in 0..field.meridians.len() {
        let cell = field
            .cells
            .iter()
            .position(|c| c.meridian == meridian)
            .unwrap();
        let row = &field.meridians[meridian];
        let mut requests = vec![
            (None, 0, 0),
            (None, field.angular.coefficients().len() - 1, 10),
        ];
        for piece in 0..field.angular.coefficients().len() {
            requests.push((Some(0), piece, 0));
            requests.push((Some(row.degree()), piece, 0));
        }
        for (fixed, piece, column) in requests {
            let n = if fixed.is_some() {
                11
            } else {
                row.degree() + 1
            };
            for i in 0..n {
                let (u, v) = if let Some(u) = fixed {
                    (u, i)
                } else {
                    (i, column)
                };
                let (_, h) = field.homogeneous_control(cell, u, piece, v, w)?;
                for coordinate in &h[..3] {
                    let mm = coordinate.div(&h[3], w)?;
                    let before = m.used;
                    let text = cache.format(&mm, m, w)?;
                    let (mantissa, exp) = text
                        .as_str()
                        .split_once("E-")
                        .unwrap_or((text.as_str(), "0"));
                    let integer =
                        BigInt::parse_bytes(mantissa.strip_suffix('.').unwrap().as_bytes(), 10)
                            .unwrap();
                    assert_eq!(
                        BigRational::new(integer, BigInt::from(10_u32).pow(exp.parse().unwrap())),
                        mm.0
                    );
                    let output_bytes = text.retained_bytes();
                    drop(text);
                    m.used -= output_bytes;
                    assert!(m.used >= before);
                }
            }
            count += 1;
        }
    }
    cache.release(m);
    w.charge(1)?;
    Ok(count)
}

#[test]
fn exact_cached_physical_mm_support_batch_uses_one_actual_ledger() {
    for (left, nb) in [(false, false), (true, true)] {
        let source =
            Hanger::new(crate::authored_hanger_bolt::tests::fixture(left, nb), None).unwrap();
        let boundary = source.boundary(BoundaryOptions::default(), None).unwrap();
        let field = EncodedField::prepare(&boundary, CurveOptions::default(), None).unwrap();
        let mut w = Work::new(16_000_000, None).unwrap();
        w.charge(field.work_used()).unwrap();
        let mut m = Memory::new(field.retained_bytes(), 64 * 1024 * 1024).unwrap();
        m.add(super::super::SCRATCH).unwrap();
        let result = cached_mm_batch(&field, &mut m, &mut w);
        eprintln!(
            "physical mm cached batch left={left} nb={nb} work={} result={result:?}",
            w.used
        );
        assert_eq!(result, Ok(180));
    }
}
