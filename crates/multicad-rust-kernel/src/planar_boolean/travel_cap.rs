//! Exact cap halfspaces of plane-terminated linear sweeps.
//!
//! A sweep `p = (x, y, 0) + d (a + b x + c y)` ends on the plane obtained by
//! eliminating the generator: `N = (-dz b, -dz c, 1 + dx b + dy c)` with
//! `N.p = dz a`. Every binary64 input is admitted exactly, so the halfspace is
//! the authored travel law rather than its rounded cap vertices, and the
//! comparison of two caps never overflows.
use num_rational::BigRational as Q;
use num_traits::{Signed, Zero};
use std::cmp::Ordering;

type Point = [f64; 3];

/// World halfspace `normal . p <= offset` bounding a sweep's material.
#[derive(Clone, Debug, PartialEq)]
pub struct TravelCap {
    normal: [Q; 3],
    offset: Q,
}

impl TravelCap {
    /// The start (`end == false`) or end cap of the sweep with local
    /// `direction` and travel `law = [a, b, c]`, placed by `axes = [u, v, n]`
    /// at `origin`. Only signed axis permutations place a sweep exactly, and
    /// the direction must advance along the local normal.
    pub fn new(
        direction: Point,
        law: Point,
        axes: [Point; 3],
        origin: Point,
        end: bool,
    ) -> Option<Self> {
        let exact = |v: f64| Q::from_float(v);
        if axes.iter().any(|axis| {
            axis.iter().filter(|x| x.abs() == 1.).count() != 1
                || axis.iter().any(|&x| x != 0. && x.abs() != 1.)
        }) {
            return None;
        }
        let [dx, dy, dz] = direction.map(exact);
        let [dx, dy, dz] = [dx?, dy?, dz?];
        if !dz.is_positive() {
            return None;
        }
        let [a, b, c] = law.map(exact);
        let [a, b, c] = [a?, b?, c?];
        let one = Q::from_integer(1.into());
        let local = [-(&dz * &b), -(&dz * &c), one + &dx * &b + &dy * &c];
        let sign = Q::from_integer(if end { 1 } else { -1 }.into());
        let mut normal = [Q::zero(), Q::zero(), Q::zero()];
        for (n, axis) in local.iter().zip(axes) {
            for (world, component) in normal.iter_mut().zip(axis) {
                if component != 0. {
                    *world += n * &sign * exact(component)?;
                }
            }
        }
        let mut offset = sign * dz * a;
        for (n, o) in normal.iter().zip(origin) {
            offset += n * exact(o)?;
        }
        Some(Self { normal, offset })
    }

    /// For caps with `other.normal = -k normal` (k > 0), the sign of the
    /// width of their common band `-other.offset / k <= normal . p <= offset`:
    /// `Less` is a gap, `Equal` exact contact and `Greater` a positive band.
    /// `None` when the caps are not antiparallel.
    pub fn opposed(&self, other: &Self) -> Option<Ordering> {
        let pivot = (0..3).find(|&i| !self.normal[i].is_zero())?;
        let k = -(&other.normal[pivot] / &self.normal[pivot]);
        if !k.is_positive() {
            return None;
        }
        if (0..3).any(|i| other.normal[i] != -(&k * &self.normal[i])) {
            return None;
        }
        Some((&self.offset * &k + &other.offset).cmp(&Q::zero()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const AXES: [Point; 3] = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];

    #[test]
    fn opposed_caps_classify_gap_contact_and_positive_band_exactly() {
        let up = [0., 0., 1.];
        let cap = |a, end| TravelCap::new(up, [a, 0., 0.], AXES, [0.; 3], end).unwrap();
        // z <= 10 against -z <= -10, -10.5 and -9.5.
        let top = cap(10., true);
        assert_eq!(top.opposed(&cap(10., false)), Some(Ordering::Equal));
        assert_eq!(top.opposed(&cap(10.5, false)), Some(Ordering::Less));
        assert_eq!(top.opposed(&cap(9.5, false)), Some(Ordering::Greater));
        // Parallel caps and non-antiparallel slopes are not opposed.
        assert_eq!(top.opposed(&cap(3., true)), None);
        let tilted = TravelCap::new(up, [10., 0.25, 0.], AXES, [0.; 3], false).unwrap();
        assert_eq!(top.opposed(&tilted), None);
    }

    #[test]
    fn one_ulp_of_travel_is_a_real_band_and_placement_is_exact() {
        let up = [0., 0., 1.];
        let a = 0.1_f64;
        let top = TravelCap::new(up, [a, 0., 0.], AXES, [0.; 3], true).unwrap();
        let low = TravelCap::new(up, [a.next_down(), 0., 0.], AXES, [0.; 3], false).unwrap();
        assert_eq!(top.opposed(&low), Some(Ordering::Greater));
        // A permuted, translated placement keeps exact contact: the local
        // normal is world +X here, and both caps lie on x = 0.75.
        let axes = [[0., 1., 0.], [0., 0., 1.], [1., 0., 0.]];
        let start = TravelCap::new(up, [0.5, 0., 0.], axes, [0.25, 0., 0.], true).unwrap();
        let next = TravelCap::new(up, [0., 0., 0.], axes, [0.75, 0., 0.], false).unwrap();
        assert_eq!(start.opposed(&next), Some(Ordering::Equal));
        let next = TravelCap::new(
            up,
            [0., 0., 0.],
            axes,
            [0.75_f64.next_down(), 0., 0.],
            false,
        );
        assert_eq!(start.opposed(&next.unwrap()), Some(Ordering::Greater));
        assert!(
            TravelCap::new(
                up,
                [a, 0., 0.],
                [[0.6, 0.8, 0.], AXES[1], AXES[2]],
                [0.; 3],
                true
            )
            .is_none()
        );
        assert!(TravelCap::new([0., 0., -1.], [a, 0., 0.], AXES, [0.; 3], true).is_none());
    }
}
