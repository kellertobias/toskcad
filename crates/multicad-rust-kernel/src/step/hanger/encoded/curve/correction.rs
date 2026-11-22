//! Exact Cartesian linear endpoint correction on an UNRESTRICTED rational
//! basis. Finite trims/geometry proof remain with the owning caller.
use super::*;
pub(super) struct LinearCorrection {
    pub(super) controls: Vec<[N; 4]>,
    pub(super) position_error: I,
    pub(super) derivative_error: I,
}
/// A(u)/W(u) -> A(u)/W(u)+d0*(b-u)/(b-a)+d1*(u-a)/(b-a).
/// Every coordinate and the weight share degree scale n+1; no independent
/// pole/weight rounding and no fitted support substitution occurs here.
pub(super) fn correct(
    h: &[[N; 4]],
    delta: &[[N; 3]; 2],
    a: &N,
    b: &N,
    memory: &mut Memory,
    w: &mut Work<'_>,
) -> Result<LinearCorrection, Error> {
    if h.len() < 2 || h.len() > 23 || a.cmp(b, w)? != Ordering::Less {
        return Err(Error::Invalid);
    }
    let span = b.sub(a, w)?;
    let width = span.enclosure(w)?;
    if width.lo <= 0. {
        return Err(Error::Unresolved);
    }
    let mut zero = true;
    for d in delta.iter().flatten() {
        if d.cmp(&N::integer(0), w)? != Ordering::Equal {
            zero = false;
            break;
        }
    }
    if zero {
        let mut controls = memory.allocate(h.len(), w)?;
        for row in h {
            w.charge(
                1 + row
                    .iter()
                    .map(|v| ((v.0.numer().bits() + v.0.denom().bits()) as usize).div_ceil(64))
                    .sum::<usize>(),
            )?;
            controls.push(row.clone());
        }
        return Ok(LinearCorrection {
            controls,
            position_error: p(0.),
            derivative_error: p(0.),
        });
    }
    let n = h.len();
    let mut alpha: [N; 3] = std::array::from_fn(|_| N::integer(0));
    let mut beta: [N; 3] = std::array::from_fn(|_| N::integer(0));
    let mut norms = [p(0.); 2];
    let mut derivative = p(0.);
    for axis in 0..3 {
        alpha[axis] = b
            .mul(&delta[0][axis], w)?
            .sub(&a.mul(&delta[1][axis], w)?, w)?
            .div(&span, w)?;
        beta[axis] = delta[1][axis].sub(&delta[0][axis], w)?.div(&span, w)?;
        for e in 0..2 {
            norms[e] = norms[e] + p(delta[e][axis].enclosure(w)?.magnitude());
        }
        derivative = derivative + p(beta[axis].enclosure(w)?.magnitude());
    }
    let mut constant = true;
    for row in h {
        if row[3].cmp(&h[0][3], w)? != Ordering::Equal {
            constant = false;
            break;
        }
    }
    if constant {
        let mut controls = memory.allocate(n, w)?;
        for (i, input) in h.iter().enumerate() {
            w.charge(4 * 512)?;
            let mut row = input.clone();
            let u = N::integer(i as i64).div(&N::integer((n - 1) as i64), w)?;
            for axis in 0..3 {
                row[axis] = row[axis].add(
                    &h[0][3].mul(&alpha[axis].add(&beta[axis].mul(&u, w)?, w)?, w)?,
                    w,
                )?;
            }
            controls.push(row);
        }
        return Ok(LinearCorrection {
            controls,
            position_error: p(norms[0].hi.max(norms[1].hi)),
            derivative_error: derivative,
        });
    }
    let mut controls = memory.allocate(n + 1, w)?;
    for i in 0..=n {
        let mut row: [N; 4] = std::array::from_fn(|_| N::integer(0));
        let mut linear_weight = N::integer(0);
        if i > 0 {
            let scale = N::integer(i as i64);
            for (slot, x) in row.iter_mut().zip(&h[i - 1]) {
                *slot = x.mul(&scale, w)?;
            }
            linear_weight = row[3].clone();
            w.charge(512)?;
        }
        if i < n {
            let scale = N::integer((n - i) as i64);
            for (slot, x) in row.iter_mut().zip(&h[i]) {
                *slot = slot.add(&x.mul(&scale, w)?, w)?;
            }
        }
        for axis in 0..3 {
            row[axis] = row[axis]
                .add(&alpha[axis].mul(&row[3], w)?, w)?
                .add(&beta[axis].mul(&linear_weight, w)?, w)?;
        }
        controls.push(row);
    }
    w.charge(64)?;
    Ok(LinearCorrection {
        controls,
        position_error: p(norms[0].hi.max(norms[1].hi)),
        derivative_error: derivative,
    })
}
#[cfg(test)]
mod tests;
