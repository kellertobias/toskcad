//! Actual literal z-only translation, with measured per-vertex rounding reserve.
use super::*;
#[cfg(test)]
thread_local! { pub(super) static FORCE_EXPLICIT:std::cell::Cell<bool> = const {std::cell::Cell::new(false)}; }
pub(super) fn checked(
    old: [[f64; 3]; 3],
    new: [[f64; 3]; 3],
    delta: f64,
    w: &mut Work<'_>,
) -> R<f64> {
    #[cfg(test)]
    if FORCE_EXPLICIT.get() {
        w.charge(64)?;
        return explicit(old, new);
    }
    // XY equality is the exact premise of this specialized transfer, not a
    // near-axis or near-rigid assumption about a generic occurrence matrix.
    if old
        .iter()
        .zip(new)
        .any(|(a, b)| a[0] != b[0] || a[1] != b[1])
    {
        return Err(E::Unresolved);
    }
    w.charge(32)?;
    z_only(old, delta)
}
pub(super) fn z_only(points: [[f64; 3]; 3], delta: f64) -> R<f64> {
    let p = points.map(|v| v.map(I::point));
    let a = std::array::from_fn(|i| p[1][i] - p[0][i]);
    let b = std::array::from_fn(|i| p[2][i] - p[0][i]);
    let n = cross(a, b);
    // ||N||>=||N||infinity. This lower bound avoids a redundant normalization
    // and is conservative even when source normals are almost axis aligned.
    let lower = n
        .into_iter()
        .map(|x| {
            if x.lo > 0. {
                x.lo
            } else if x.hi < 0. {
                -x.hi
            } else {
                0.
            }
        })
        .fold(0., f64::max);
    let edges = [a[0], a[1], b[0], b[1]]
        .into_iter()
        .fold(I::point(0.), |s, x| s + I::point(x.magnitude()));
    // da,db are parallel to z, so da×db=0 exactly. Each edge perturbation has
    // magnitude<=2delta; only the unchanged XY edge components enter ΔN.
    let epsilon = I::point(2.) * I::point(delta) * edges;
    let ratio = epsilon / I::point(lower);
    if lower <= 0. || !ratio.valid() || ratio.lo < 0. || ratio.hi >= 1. {
        return Err(E::Unresolved);
    }
    // The ball B(N,epsilon) excludes zero, preserves orientation, and subtends
    // at most asin(epsilon/||N||). atan form retains directed interval support.
    let angle = (ratio / (I::point(1.) - ratio.square()).sqrt()).atan();
    if !angle.valid() {
        return Err(E::Unresolved);
    }
    Ok(angle.hi)
}
#[cfg(test)]
pub(super) fn explicit(old: [[f64; 3]; 3], new: [[f64; 3]; 3]) -> R<f64> {
    let old = normal(old);
    let new = normal(new);
    let projection = dot(old, new);
    if !projection.valid() || projection.lo <= 0. {
        return Err(E::Unresolved);
    }
    let c = cross(old, new);
    let magnitude = dot(c, c).sqrt();
    let angle = (I::new(0., magnitude.hi) / I::point(projection.lo)).atan();
    if !angle.valid() {
        return Err(E::Unresolved);
    }
    Ok(angle.hi)
}
