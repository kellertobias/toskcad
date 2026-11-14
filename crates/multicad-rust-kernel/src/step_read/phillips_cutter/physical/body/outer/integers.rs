//! Degree and multiplicities are STEP INTEGERs, not nearby REALs rounded to an
//! integer by binary64. Keep the original exact lexical category and value.
use super::*;
fn entity(shell: &Shell<'_, '_>, id: u32, work: &mut Work<'_>) -> Result<()> {
    let e = shell.simple(id)?;
    let (degrees, lists) = match e.name {
        "B_SPLINE_SURFACE_WITH_KNOTS" => (2, &[8, 9][..]),
        "B_SPLINE_CURVE_WITH_KNOTS" => (1, &[6][..]),
        _ => return Ok(()),
    };
    let text = e.source.ok_or(ReadError::Unsupported)?;
    work.step(text.len())?;
    let tokens = crate::step_read::part21::numeric_spellings(e).ok_or(ReadError::Unsupported)?;
    let mut expected = vec![crate::step::helical::DEGREE as u64; degrees];
    for &index in lists {
        let Value::List(values) = e.args.get(index).ok_or(ReadError::Invalid)? else {
            return Err(ReadError::Invalid);
        };
        require(values.len() >= 2)?;
        expected.extend((0..values.len()).map(|i| {
            (crate::step::helical::DEGREE + usize::from(i == 0 || i + 1 == values.len())) as u64
        }));
    }
    require(tokens.len() >= expected.len())?;
    for (token, integer) in tokens.into_iter().zip(expected) {
        require(token.parse::<u64>().ok() == Some(integer))?;
    }
    Ok(())
}
pub(super) fn validate(body: &Body<'_>, shell: &Shell<'_, '_>, work: &mut Work<'_>) -> Result<()> {
    let mut checked = BTreeSet::new();
    for &id in body.tensor_face_entities() {
        let (surface, _, bounds) = helical::face(shell, id)?;
        entity(shell, surface, work)?;
        for bound in bounds {
            for (edge, _) in helical::ring(shell, bound, true)? {
                let args = shell.args(edge, "EDGE_CURVE")?;
                let curve = reference(args.get(3).ok_or(ReadError::Invalid)?)?;
                if checked.insert(curve) {
                    entity(shell, curve, work)?;
                }
            }
        }
    }
    Ok(())
}
