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
    upper(angle)
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
    upper(angle)
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn real_literal_translation_cone_matches_explicit_normals_and_refuses_uncertainty() {
        let old = [[1., 0., 0.], [1., 0.1, 0.02], [0.9, 0.05, 0.2]];
        let new = old.map(|mut p| {
            p[2] = (p[2] - 12.7) - 2.15;
            p
        });
        let delta = old.into_iter().zip(new).fold(0_f64, |e, (a, b)| {
            e.max((I::point(a[2]) - I::point(12.7) - I::point(2.15) - I::point(b[2])).magnitude())
        });
        let mut w = Work {
            remaining: 10000,
            used: 0,
            interrupt: None,
        };
        let cone = checked(old, new, delta, &mut w).unwrap();
        let explicit_bound = explicit(old, new).unwrap();
        assert!(cone < 1e-10 && explicit_bound < 1e-10);
        assert!(z_only(old, 1.).is_err());
        assert!(z_only([[0.; 3]; 3], delta).is_err());
        assert!(z_only(old, f64::NAN).is_err());
        let mut changed = new;
        changed[1][0] += 0.01;
        assert!(checked(old, changed, delta, &mut w).is_err());
        assert!(explicit(old, [new[0], new[2], new[1]]).is_err());
        assert!(explicit([[0.; 3]; 3], new).is_err());
    }
}
