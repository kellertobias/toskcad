//! A finite polynomial isoline, proved over every original Bernstein cell.
//! This yields its point set and positive traversal, never a surface/body.
use super::bezier::{self, Coordinates};
use super::*;
use crate::step_read::brep::check;

pub(super) fn segment(
    source: &Source<'_>,
    definition: &super::super::uv::Definition,
    budget: &mut Budget<'_>,
) -> Result<Option<AxisSegment>> {
    let Curve::Nurbs {
        curve,
        closed: false,
    } = &definition.geometry
    else {
        return Ok(None);
    };
    if curve.degree() < 2 {
        return Ok(None);
    }
    budget.charge(4096)?;
    let field = fields(source, definition.curve_id, 1, budget)?;
    let controls = list(field.controls)?;
    let first = reference(controls.first().ok_or(ReadError::Unsupported)?)?;
    let last = reference(controls.last().ok_or(ReadError::Unsupported)?)?;
    let a = point(source, first, 2, budget)?;
    let b = point(source, last, 2, budget)?;
    let [_, values] = source.args(first, "CARTESIAN_POINT")? else {
        return Err(ReadError::Invalid);
    };
    let values = list(values)?;
    check(values.len() == 2, ReadError::Unsupported)?;
    for moving in 0..2 {
        let fixed = 1 - moving;
        let rate = b[moving] - a[moving];
        if !rate.valid() || !(rate.lo > 0. || rate.hi < 0.) {
            continue;
        }
        let positive = rate.lo > 0.;
        let mut constant = true;
        for p in controls {
            budget.charge(128)?;
            let id = reference(p)?;
            let [_, v] = source.args(id, "CARTESIAN_POINT")? else {
                return Err(ReadError::Invalid);
            };
            let v = list(v)?;
            check(v.len() == 2, ReadError::Unsupported)?;
            if !exact_same(source, (first, &values[fixed]), (id, &v[fixed]), budget)? {
                constant = false;
                break;
            }
        }
        if !constant {
            continue;
        }
        let data = bezier::original_curve(
            source,
            definition.curve_id,
            curve,
            Coordinates::Raw(2),
            budget,
        )?;
        check(
            data.matches(source, definition.curve_id)
                && matches!(data.coordinates(), Coordinates::Raw(2))
                && data.degree() == curve.degree()
                && !data.cells().is_empty(),
            ReadError::Unsupported,
        )?;
        check(
            data.peak_bytes()
                .checked_add(4096)
                .is_some_and(|peak| peak <= 8 * 1024 * 1024),
            ReadError::Limit,
        )?;
        for cell in data.cells() {
            budget.charge(128)?;
            if cell
                .controls()
                .iter()
                .any(|p| p[3].lo != 1. || p[3].hi != 1.)
            {
                return Ok(None);
            }
        }
        // Unit weights earn a polynomial. Each derivative Bernstein control
        // is p*(P[i+1]-P[i]); all have the same strict sign on every local
        // [0,1] cell. Original clamped, ordered knots give positive global
        // cell maps and endpoint interpolation, hence precisely this segment.
        let mut monotone = true;
        for cell in data.cells() {
            for pair in cell.controls().windows(2) {
                budget.charge(64)?;
                let difference = pair[1][moving] - pair[0][moving];
                if !difference.valid()
                    || !(if positive {
                        difference.lo > 0.
                    } else {
                        difference.hi < 0.
                    })
                {
                    monotone = false;
                    break;
                }
            }
            if !monotone {
                break;
            }
        }
        if monotone {
            budget.charge(64)?;
            return Ok(Some(AxisSegment {
                moving,
                endpoints: [a, b],
            }));
        }
    }
    Ok(None)
}

#[cfg(test)]
mod tests;
