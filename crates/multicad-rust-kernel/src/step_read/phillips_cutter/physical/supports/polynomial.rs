//! Directed Bernstein algebra through degree12, including every mixed term.
use super::*;
pub(super) type Poly = Vec<I>;
const fn choices() -> [[u32; 13]; 13] {
    let mut c = [[0; 13]; 13];
    let mut n = 0;
    while n < 13 {
        c[n][0] = 1;
        c[n][n] = 1;
        let mut k = 1;
        while k < n {
            c[n][k] = c[n - 1][k - 1] + c[n - 1][k];
            k += 1;
        }
        n += 1;
    }
    c
}
const CHOOSE: [[u32; 13]; 13] = choices();
pub(super) fn product(a: &[I], b: &[I]) -> Poly {
    let n = a.len() - 1;
    let m = b.len() - 1;
    (0..=n + m)
        .map(|k| {
            (k.saturating_sub(m)..=n.min(k)).fold(I::point(0.), |sum, i| {
                let weight = I::point(CHOOSE[n][i] as f64) * I::point(CHOOSE[m][k - i] as f64)
                    / I::point(CHOOSE[n + m][k] as f64);
                sum + a[i] * b[k - i] * weight
            })
        })
        .collect()
}
pub(super) fn square(a: &[I]) -> Poly {
    product(a, a)
}
pub(super) fn add(a: &[I], b: &[I]) -> Poly {
    a.iter().zip(b).map(|(&a, &b)| a + b).collect()
}
pub(super) fn sub(a: &[I], b: &[I]) -> Poly {
    a.iter().zip(b).map(|(&a, &b)| a - b).collect()
}
pub(super) fn hull(a: &[I]) -> I {
    if a.iter().any(|x| !x.valid()) {
        return I::new(f64::NAN, f64::NAN);
    }

    I::new(
        a.iter().map(|x| x.lo).fold(f64::INFINITY, f64::min),
        a.iter().map(|x| x.hi).fold(f64::NEG_INFINITY, f64::max),
    )
}
pub(super) fn split(p: [[I; 3]; 4]) -> [[[I; 3]; 4]; 2] {
    let mut q = p;
    let mut a = p;
    let mut b = p;
    for level in 1..4 {
        for i in 0..4 - level {
            q[i] = std::array::from_fn(|j| (q[i][j] + q[i + 1][j]) / I::point(2.));
        }
        a[level] = q[0];
        b[3 - level] = q[3 - level];
    }
    [a, b]
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn degree_twelve_mixed_products_and_split_preserve_exact_dyadic_oracles() {
        assert_eq!(choices()[12][6], 924);
        let a = [I::point(0.), I::point(0.), I::point(0.), I::point(1.)];
        let b = [I::point(1.), I::point(0.), I::point(0.), I::point(0.)];
        let p = product(&a, &b);
        assert!(p[3].lo <= 1. / 20. && p[3].hi >= 1. / 20.);
        let q = square(&p);
        assert!(q[6].lo <= 1. / 924. && q[6].hi >= 1. / 924.);
        let controls = std::array::from_fn(|i| [I::point(i as f64); 3]);
        let [a, b] = split(controls);
        assert!(a[3][0].lo <= 1.5 && a[3][0].hi >= 1.5);
        assert!(b[0][0].lo <= 1.5 && b[0][0].hi >= 1.5);
        assert!(!hull(&[I::point(0.), I::new(f64::NAN, f64::NAN)]).valid());
    }
}
