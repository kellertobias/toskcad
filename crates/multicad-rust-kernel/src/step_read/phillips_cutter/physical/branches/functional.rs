//! One immutable section functional for the complete finite cubic edge.
use super::*;
#[derive(Clone, Copy)]
pub(in super::super) struct Functional {
    pub direction: [f64; 3],
    pub pivot: usize,
    pub free: [usize; 2],
}
fn dot(a: [f64; 3], b: [I; 3]) -> I {
    (0..3).fold(I::point(0.), |s, k| s + I::point(a[k]) * b[k])
}
impl Functional {
    pub fn new(spans: &[CubicTrimSpan]) -> BranchResult<Self> {
        let first = spans.first().ok_or(PhillipsReadError::Invalid)?.poles_mm[0];
        let last = spans.last().ok_or(PhillipsReadError::Invalid)?.poles_mm[3];
        let chord: [f64; 3] = std::array::from_fn(|k| last[k] - first[k]);
        if chord.iter().any(|x| !x.is_finite()) {
            return Err(PhillipsReadError::Unsupported);
        }
        let pivot = (0..3)
            .max_by(|&i, &j| chord[i].abs().total_cmp(&chord[j].abs()))
            .unwrap();
        if chord[pivot] == 0. {
            return Err(PhillipsReadError::Unsupported);
        }
        let mut direction = chord.map(|x| x / chord[pivot]);
        // This is a sampling-independent exact chosen functional, not a
        // fitted surface frame. Its pivot is exactly one; |other terms|<=1.
        direction[pivot] = 1.;
        Ok(Self {
            direction,
            pivot,
            free: [(pivot + 1) % 3, (pivot + 2) % 3],
        })
    }
    pub fn basis(&self) -> [[I; 3]; 2] {
        self.free.map(|j| {
            std::array::from_fn(|k| {
                I::point(if k == j {
                    1.
                } else if k == self.pivot {
                    -self.direction[j]
                } else {
                    0.
                })
            })
        })
    }
    pub fn derivative(&self, poles: [[I; 3]; 4]) -> [I; 3] {
        std::array::from_fn(|j| {
            dot(
                self.direction,
                std::array::from_fn(|k| I::point(3.) * (poles[j + 1][k] - poles[j][k])),
            )
        })
    }
    // The parsed clamped B-spline table supplies exact shared endpoint
    // CONTROL identities at adjacent spans. This checks their derivative
    // orientation; it never joins disconnected curves by numeric proximity.
    pub fn monotone(&self, spans: &[CubicTrimSpan]) -> BranchResult<bool> {
        let mut sign = None;
        for span in spans {
            let d = self.derivative(span.poles_mm.map(|p| p.map(support::adjacent)));
            let current = if d.iter().all(|v| v.valid() && v.lo > 0.) {
                true
            } else if d.iter().all(|v| v.valid() && v.hi < 0.) {
                false
            } else {
                return Err(PhillipsReadError::Unsupported);
            };
            if sign.is_some_and(|s| s != current) {
                return Err(PhillipsReadError::Unsupported);
            }
            sign = Some(current);
        }
        sign.ok_or(PhillipsReadError::Invalid)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    fn line(a: [f64; 3], b: [f64; 3]) -> CubicTrimSpan {
        CubicTrimSpan {
            parameters: [0., 1.],
            poles_mm: std::array::from_fn(|i| {
                std::array::from_fn(|k| a[k] + (b[k] - a[k]) * (i as f64 / 3.))
            }),
        }
    }
    #[test]
    fn actual_cubic_section_functional_preserves_both_free_directions_and_orientation() {
        for (a, b) in [([0.; 3], [2., 1., -0.5]), ([4., 2., 1.], [-2., 4., 0.])] {
            let span = line(a, b);
            let f = Functional::new(&[span]).unwrap();
            for basis in f.basis() {
                let d = dot(f.direction, basis);
                assert!(d.lo <= 0. && d.hi >= 0.);
            }
            assert_eq!(f.direction[f.pivot], 1.);
            assert_eq!(f.monotone(&[span]).unwrap(), b[f.pivot] > a[f.pivot]);
        }
    }
    #[test]
    fn closed_stationary_and_reversed_piecewise_carriers_are_not_a_single_section_chart() {
        assert!(Functional::new(&[]).is_err());
        assert!(Functional::new(&[line([0.; 3], [0.; 3])]).is_err());
        assert!(Functional::new(&[line([0.; 3], [f64::INFINITY, 0., 0.])]).is_err());
        let a = line([0.; 3], [1., 0., 0.]);
        let b = line([1., 0., 0.], [0.5, 0., 0.]);
        let f = Functional::new(&[a, b]).unwrap();
        assert!(f.monotone(&[a, b]).is_err());
        let f = Functional::new(&[a]).unwrap();
        assert!(f.monotone(&[]).is_err());
        let mut stopped = a;
        stopped.poles_mm[1] = stopped.poles_mm[0];
        assert!(f.monotone(&[stopped]).is_err());
    }
}
