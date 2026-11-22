//! Original-REAL IP2 for a clamped bilinear boundary against an actual LINE.
//! This is ONLY a sense certificate. The caller separately earns whole finite
//! geometric coverage and raw-coordinate transfer before publishing admission.
use super::*;

fn dot(a: V, b: V) -> I {
    (0..3).fold(I::point(0.), |sum, i| sum + a[i] * b[i])
}

pub(super) fn certify(
    source: &Source<'_>,
    curve_id: u32,
    definition: &super::super::uv::Definition,
    surface: &Surface,
    edge: &Edge,
    budget: &mut Budget<'_>,
) -> Result<()> {
    budget.charge(4096)?;
    let (Surface::Nurbs { surface, closed }, Curve::Line { .. }) = (surface, &edge.curve) else {
        return Err(ReadError::Unsupported);
    };
    if *closed != [false, false] || surface.degrees() != [1, 1] || surface.counts() != [2, 2] {
        return Err(ReadError::Unsupported);
    }
    let iso = isoline(source, definition, budget)?;
    let fixed = 1 - iso.moving;
    let field = fields(source, definition.surface, 2, budget)?;
    for axis in 0..2 {
        integer(source, definition.surface, field.degrees[axis], 1, budget)?;
        let Some((multiplicities, knots)) = field.knots[axis] else {
            // Bezier defaults also have the exact natural range [0,1].
            continue;
        };
        let (multiplicities, knots) = (list(multiplicities)?, list(knots)?);
        if multiplicities.len() != 2 || knots.len() != 2 {
            return Err(ReadError::Unsupported);
        }
        for m in multiplicities {
            integer(source, definition.surface, m, 2, budget)?;
        }
        let a = scalar(source, definition.surface, &knots[0], budget)?;
        let b = scalar(source, definition.surface, &knots[1], budget)?;
        let gap = b - a;
        if !gap.valid() || gap.lo <= 0. {
            return Err(ReadError::Unsupported);
        }
    }
    let range = surface.range(fixed);
    let fixed_index = if iso.origin[fixed].lo <= range[0] && iso.origin[fixed].hi >= range[0] {
        // Interval containment is only a proposal; exact original equality
        // with this boundary is independently required immediately below.
        0
    } else if iso.origin[fixed].lo <= range[1] && iso.origin[fixed].hi >= range[1] {
        1
    } else {
        return Err(ReadError::Unsupported);
    };
    boundary(
        source,
        definition.surface,
        &field,
        fixed,
        fixed_index,
        &iso,
        budget,
    )?;
    let rows = list(field.controls)?;
    if rows.len() != 2
        || rows
            .iter()
            .any(|row| !list(row).is_ok_and(|r| r.len() == 2))
    {
        return Err(ReadError::Unsupported);
    }
    let mut controls = [[I::point(0.); 3]; 2];
    for (i, control) in controls.iter_mut().enumerate() {
        budget.charge(256)?;
        let (u, v) = if iso.moving == 0 {
            (i, fixed_index)
        } else {
            (fixed_index, i)
        };
        *control = point(source, reference(&list(&rows[u])?[v])?, 3, budget)?;
        if let Some(weights) = field.weights {
            let rows = list(weights)?;
            if rows.len() != 2 {
                return Err(ReadError::Unsupported);
            }
            let row = list(&rows[u])?;
            if row.len() != 2 {
                return Err(ReadError::Unsupported);
            }
            let weight = scalar(source, definition.surface, &row[v], budget)?;
            if !weight.valid() || weight.lo <= 0. {
                return Err(ReadError::Unsupported);
            }
        }
    }
    let id = unwrapped(source, curve_id, budget)?;
    let [_, _, direction] = source.args(id, "LINE")? else {
        return Err(ReadError::Unsupported);
    };
    let target = vector(source, reference(direction)?, 3, budget)?;
    // On an actual clamped degree-one boundary, including positive rational
    // weights, dP/ds = (w0*w1/W(s)^2)*(P1-P0)/(k1-k0).
    // Every omitted factor is STRICTLY POSITIVE over the whole natural span.
    // The UV LINE rate, or two-pole positive-rational UV endpoint direction,
    // binds orientation to its original increasing parameter. Thus this one
    // directed dot proves positive target parameter flow on EVERY finite
    // used subwindow, not merely the ordering of sampled or cached points.
    budget.charge(128)?;
    let flow = dot(difference(controls[1], controls[0]), target) * iso.rate;
    if !flow.valid() || flow.lo <= 0. {
        return Err(ReadError::Unsupported);
    }
    Ok(())
}

#[cfg(test)]
mod tests;
