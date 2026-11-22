//! Exact side-curve support experiment; finite trim/endpoint/file admission
//! remains separate. Every input is an actual retained source tensor slice.
use super::*;
use crate::authored_hanger_bolt::{BoundaryOptions, Hanger};
fn support(
    field: &EncodedField<'_>,
    index: usize,
    piece: usize,
    fixed_u: Option<usize>,
    column: usize,
    m: &mut Memory,
    w: &mut Work<'_>,
) -> Result<Vec<([EncodedScalar; 3], EncodedScalar)>, Error> {
    let row = &field.meridians[field.cells[index].meridian];
    let count = if fixed_u.is_some() {
        11
    } else {
        row.degree() + 1
    };
    let mut output = m.allocate(count, w)?;
    for i in 0..count {
        let (u, v) = if let Some(u) = fixed_u {
            (u, i)
        } else {
            (i, column)
        };
        let control = field.control(
            index,
            u,
            piece,
            v,
            DomainOptions {
                maximum_work: w.remaining(),
                ..DomainOptions::default()
            },
            None,
        )?;
        w.charge(control.work_used)?;
        let h = control
            .homogeneous
            .each_ref()
            .map(|x| x.number(w))
            .into_iter()
            .collect::<Result<Vec<_>, Error>>()?;
        assert!(h[3].cmp(&N::integer(0), w)? == Ordering::Greater);
        let point = [
            EncodedScalar::retain(&h[0].div(&h[3], w)?, m, w)?,
            EncodedScalar::retain(&h[1].div(&h[3], w)?, m, w)?,
            EncodedScalar::retain(&h[2].div(&h[3], w)?, m, w)?,
        ];
        let weight = EncodedScalar::retain(&h[3], m, w)?;
        // No pole quantization or independent Z rounding is necessary here.
        // Rebuilding homogeneous controls yields EXACT original coefficients.
        for a in 0..3 {
            assert_eq!(point[a].number(w)?.mul(&weight.number(w)?, w)?.0, h[a].0);
        }
        for scalar in point.iter().chain([&weight]) {
            let value = scalar.number(w)?;
            let before = m.used;
            let decimal = super::super::decimal::decimal(&value, m, w)?;
            assert!(!decimal.as_str().is_empty());
            drop(decimal);
            m.used = before;
        }
        output.push((point, weight));
    }
    Ok(output)
}
#[test]
fn actual_rational_meridian_and_helical_rail_supports_have_exact_finite_decimal_poles() {
    for (left, nb) in [(false, false), (true, true)] {
        let source =
            Hanger::new(crate::authored_hanger_bolt::tests::fixture(left, nb), None).unwrap();
        let boundary = source.boundary(BoundaryOptions::default(), None).unwrap();
        let field = EncodedField::prepare(&boundary, CurveOptions::default(), None).unwrap();
        let mut total_work = field.work_used();
        let mut maximum_query_work = 0;
        let mut m = Memory::new(field.retained_bytes(), 64 * 1024 * 1024).unwrap();
        m.add(SCRATCH).unwrap();
        let mut count = 0;
        for meridian in 0..field.meridians.len() {
            let index = field
                .cells
                .iter()
                .position(|c| c.meridian == meridian)
                .unwrap();
            for u in [0, field.meridians[meridian].degree()] {
                for piece in 0..field.angular.coefficients().len() {
                    let before = m.used;
                    // A support query borrows the already prepared source.
                    // This is not a claim that all decimal encodings fit in
                    // one file's shared sixteen-million-work ledger.
                    let mut w = Work::new(16_000_000, None).unwrap();
                    let curve = support(&field, index, piece, Some(u), 0, &mut m, &mut w).unwrap();
                    total_work += w.used;
                    maximum_query_work = maximum_query_work.max(w.used);
                    drop(curve);
                    m.used = before;
                    count += 1;
                }
            }
            for (piece, column) in [(0, 0), (field.angular.coefficients().len() - 1, 10)] {
                let before = m.used;
                let mut w = Work::new(16_000_000, None).unwrap();
                let curve = support(&field, index, piece, None, column, &mut m, &mut w).unwrap();
                total_work += w.used;
                maximum_query_work = maximum_query_work.max(w.used);
                drop(curve);
                m.used = before;
                count += 1;
            }
        }
        eprintln!(
            "exact side supports left={left} nb={nb} count={count} aggregate_work={total_work} maximum_query_work={maximum_query_work}"
        );
    }
}

#[test]
fn exact_side_supports_respect_real_shared_batch_and_query_limits() {
    let source = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, false),
        None,
    )
    .unwrap();
    let boundary = source.boundary(BoundaryOptions::default(), None).unwrap();
    let field = EncodedField::prepare(&boundary, CurveOptions::default(), None).unwrap();
    let base = field.retained_bytes() + SCRATCH;
    for (work, bytes, cancel, expected) in [
        (1, 64 * 1024 * 1024, false, Error::Limit),
        (16, 64 * 1024 * 1024, false, Error::Limit),
        (16_000_000, base, false, Error::Limit),
        (16_000_000, 64 * 1024 * 1024, true, Error::Cancelled),
    ] {
        let interrupt = || cancel;
        let mut w = Work::new(work, Some(&interrupt)).unwrap();
        let mut m = Memory::new(base, bytes).unwrap();
        assert_eq!(
            support(&field, 0, 0, Some(0), 0, &mut m, &mut w).err(),
            Some(expected)
        );
    }
    // The complete representative batch genuinely exceeds the original shared
    // budget. Keep that refusal; the separate positive queries above do not
    // establish an affordable whole-file pipeline or reset a caller's ledger.
    let mut queries = Vec::new();
    for meridian in 0..field.meridians.len() {
        let index = field
            .cells
            .iter()
            .position(|c| c.meridian == meridian)
            .unwrap();
        for u in [0, field.meridians[meridian].degree()] {
            for piece in 0..field.angular.coefficients().len() {
                queries.push((index, piece, Some(u), 0));
            }
        }
        for (piece, column) in [(0, 0), (field.angular.coefficients().len() - 1, 10)] {
            queries.push((index, piece, None, column));
        }
    }
    // Even without charging the already prepared field again, these numeric
    // queries alone exceed the shared budget. Earlier successes consume work;
    // the actual later refusal must not reset it.
    let mut w = Work::new(16_000_000, None).unwrap();
    let mut m = Memory::new(base, 64 * 1024 * 1024).unwrap();
    let result = queries
        .into_iter()
        .try_for_each(|(index, piece, fixed_u, column)| {
            let before = m.used;
            let curve = support(&field, index, piece, fixed_u, column, &mut m, &mut w)?;
            drop(curve);
            m.used = before;
            Ok(())
        });
    assert_eq!(result, Err(Error::Limit));
}
