//! Whole periodic radial-law enclosure, also usable independently of a spatial
//! query. Large/coarse phases return the earned global range, never a guessed
//! modulo representative. Precise domains retain every intersecting real band.
use super::*;
fn unit(q: I, a: I, b: I) -> I {
    let n = q - a;
    let d = b - a;
    // Direct directed division avoids reciprocal overflow for narrow positive
    // denominators. Caller proved overlap and d.lo>0. Thus min quotient<=1 and
    // max quotient>=0, even if other endpoint quotients overflow before clipping.
    let quotients = [n.lo / d.lo, n.lo / d.hi, n.hi / d.lo, n.hi / d.hi];
    I::new(
        quotients
            .into_iter()
            .fold(f64::INFINITY, f64::min)
            .next_down()
            .max(0.),
        quotients
            .into_iter()
            .fold(f64::NEG_INFINITY, f64::max)
            .next_up()
            .min(1.),
    )
}
pub(super) fn bounds(source: &PeriodicThreadSource, phase: I, work: &mut Work<'_>) -> FResult<I> {
    work.step(64)?;
    if !phase.valid() {
        return Err(PhillipsReadError::Invalid);
    }
    let radii = source.radius_nodes_mm();
    let global = I::new(radii[1], radii[0]);
    let pitch = I::point(source.pitch_mm());
    if phase.width() >= source.pitch_mm() {
        return Ok(global);
    }
    let nodes = source.phases();
    let period = (phase - nodes[0]) / pitch;
    if !period.valid() || period.magnitude() > 1e12 {
        return Ok(global);
    }
    let first = period.lo.floor() as i64;
    let last = period.hi.floor() as i64;
    // Input width<pitch and |period|<=1e12 leave <1/64 extra quotient width
    // after these directed binary64 operations. Floor indices differ by at most
    // two; reserve one additional neighbor. Casts remain exact below2^53.
    assert!(
        last - first <= 3,
        "bounded periodic quotient has a fixed-size integer scan"
    );
    let mut lower = global.hi;
    let mut upper = global.lo;
    let mut found = false;
    for k in first..=last {
        let q = phase - I::point(k as f64) * pitch;
        for band in 0..4 {
            work.step(128)?;
            if q.hi < nodes[band].lo || q.lo > nodes[band + 1].hi {
                continue;
            }
            let u = unit(q, nodes[band], nodes[band + 1]);
            assert!(
                u.valid(),
                "overlapping ordered phase boxes have a nonempty unit intersection"
            );
            let r = source.radii(band, u);
            lower = lower.min(r.lo.max(global.lo));
            upper = upper.max(r.hi.min(global.hi));
            found = true;
        }
    }
    // Exact nodes0..4 cover one complete real period by identity. The directed
    // quotient includes that period for every phase in the nonempty input, and
    // every intersecting band above retains its full unit intersection. Radii
    // are finite affine combinations of positive bounded source nodes.
    assert!(
        found && I::new(lower, upper).valid(),
        "complete real phase partition supplies a finite radial enclosure"
    );
    Ok(I::new(lower, upper))
}
#[cfg(test)]
mod tests;
