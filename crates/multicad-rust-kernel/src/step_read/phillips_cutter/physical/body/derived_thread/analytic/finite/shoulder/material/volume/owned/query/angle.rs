//! Rigorous local angle, with no libm atan/heuristic ulp enclosure. On a
//! positive-X halfplane, atan2 has no branch cut; its rectangle extrema are
//! at corners because each boundary restriction has a fixed derivative sign.
//! Quarter-turn rotation selects such a halfplane before applying that fact.
use super::*;
fn small(r: I) -> I {
    let square = r.square();
    let mut power = r;
    let mut sum = r;
    for k in 1..32 {
        power = power * square;
        let term = power / I::point((2 * k + 1) as f64);
        sum = if k % 2 == 0 { sum + term } else { sum - term };
    }
    let tail = (power * square / I::point(65.)).magnitude();
    sum + I::new(-tail, tail)
}
fn ratio(a: f64, b: f64, pi: I) -> Option<I> {
    if a == 0. {
        return Some(I::point(0.));
    }
    // a/b<=1 is proved by the caller's exact scalar ordering. Half-angle
    // reduction gives r<=sqrt(2)-1<1/2, and the omitted alternating tail is
    // bounded by |r|^65/65, including every directed arithmetic error.
    let q = I::point(a) / I::point(b);
    if !q.valid() {
        return None;
    }
    let q = I::new(q.lo.max(0.), q.hi.min(1.));
    let r = q / (I::point(1.) + (I::point(1.) + q.square()).sqrt());
    let r = I::new(r.lo.max(0.), r.hi.min(0.5));
    let v = I::point(2.) * small(r);
    Some(I::new(v.lo.max(0.), v.hi.min((pi / I::point(4.)).hi)))
}
fn point(x: f64, y: f64, pi: I) -> Option<I> {
    let yabs = y.abs();
    let a = if yabs <= x {
        ratio(yabs, x, pi)?
    } else {
        pi / I::point(2.) - ratio(x, yabs, pi)?
    };
    Some(if y < 0. { -a } else { a })
}
pub(super) fn enclose(p: [I; 2], work: &mut Work<'_>) -> FResult<Option<(I, I)>> {
    work.step(256)?;
    if p.into_iter().any(|v| !v.valid()) {
        return Err(PhillipsReadError::Invalid);
    }
    let (local, turn) = if p[0].lo > 0. {
        (p, 0)
    } else if p[1].lo > 0. {
        ([p[1], -p[0]], 1)
    } else if p[0].hi < 0. {
        ([-p[0], -p[1]], 2)
    } else if p[1].hi < 0. {
        ([-p[1], p[0]], 3)
    } else {
        return Ok(None);
    };
    work.step(4096)?;
    let rotation = rotation::Rotation::new(work)?;
    let pi = rotation.pi(work)?;
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for x in [local[0].lo, local[0].hi] {
        for y in [local[1].lo, local[1].hi] {
            let Some(a) = point(x, y, pi) else {
                return Ok(None);
            };
            lo = lo.min(a.lo);
            hi = hi.max(a.hi);
        }
    }
    let angle = I::new(lo, hi) + I::point(turn as f64 / 2.) * pi;
    Ok(angle.valid().then_some((angle, pi)))
}
#[cfg(test)]
mod tests;
