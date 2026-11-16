//! Original P11 point lattice and complete-root Hessian bounds. No fitted
//! proposal replaces geometry: all strips are exact directed Bernstein sums.
use super::*;
type B = [I; 3];
pub(super) const ROOT_WORK: usize = 16_384;
pub(super) const STRIP_WORK: usize = 4096;
pub(super) const SPLIT_WORK: usize = 1024;
const ENDPOINT_WORK: usize = 64;
pub(super) const C0_WORK: usize = 512;
pub(super) struct Root {
    controls: [[P; 12]; 12],
    hessian: [f64; 3],
}
pub(super) struct Lattice {
    // SAME private immutable root allocation, never an external point cache.
    owner: Option<Arc<Root>>,
    points: Vec<B>,
    divisions: [usize; 2],
}
impl Root {
    pub(super) fn new(
        source: &crate::tapered_thread::RetainedTensor11,
        index: [usize; 2],
    ) -> R<Self> {
        if index.into_iter().zip(source.spans()).any(|(a, n)| a >= n) {
            return Err(PhillipsReadError::Invalid);
        }
        let height = source.dimensions()[1];
        let controls = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                source.poles()[(11 * index[0] + i) * height + 11 * index[1] + j]
            })
        });
        let mut hessian = [0f64; 3];
        for (direction, bound) in hessian.iter_mut().enumerate() {
            let count = if direction == 1 {
                [11, 11]
            } else if direction == 0 {
                [10, 12]
            } else {
                [12, 10]
            };
            for i in 0..count[0] {
                for j in 0..count[1] {
                    let norm = (0..3).fold(I::point(0.), |sum, k| {
                        let value = match direction {
                            0 => {
                                I::point(110.)
                                    * (I::point(controls[i + 2][j][k])
                                        - I::point(2.) * I::point(controls[i + 1][j][k])
                                        + I::point(controls[i][j][k]))
                            }
                            1 => {
                                I::point(121.)
                                    * (I::point(controls[i + 1][j + 1][k])
                                        - I::point(controls[i + 1][j][k])
                                        - I::point(controls[i][j + 1][k])
                                        + I::point(controls[i][j][k]))
                            }
                            _ => {
                                I::point(110.)
                                    * (I::point(controls[i][j + 2][k])
                                        - I::point(2.) * I::point(controls[i][j + 1][k])
                                        + I::point(controls[i][j][k]))
                            }
                        };
                        sum + I::point(value.magnitude())
                    });
                    if !norm.valid() {
                        return Err(PhillipsReadError::Unsupported);
                    }
                    *bound = bound.max(norm.hi);
                }
            }
        }
        Ok(Self { controls, hessian })
    }
    pub(super) fn bytes(&self) -> usize {
        std::mem::size_of::<Self>() + 64
    }
    pub(super) fn interpolation_error(&self, domain: [[f64; 2]; 2]) -> f64 {
        let [u, v] = domain.map(|[a, b]| I::point(b) - I::point(a));
        (I::point(self.hessian[0]) * u.square() / I::point(8.)
            + I::point(self.hessian[1]) * u * v / I::point(4.)
            + I::point(self.hessian[2]) * v.square() / I::point(8.))
        .hi
    }
    fn strip_at(&self, parameter: f64) -> [B; 12] {
        if parameter == 0. {
            return self.controls[0].map(|p| p.map(I::point));
        }
        if parameter == 1. {
            return self.controls[11].map(|p| p.map(I::point));
        }
        let u = I::point(parameter);
        let complement = I::point(1.) - u;
        let mut basis = [I::point(0.); 12];
        basis[0] = I::point(1.);
        for degree in 1..12 {
            for i in (1..=degree).rev() {
                basis[i] = u * basis[i - 1] + complement * basis[i];
            }
            basis[0] = basis[0] * complement;
        }
        // Exact partition unity allows the anchor once. This avoids a
        // translation-dependent interval floor without normalizing weights.
        std::array::from_fn(|j| {
            std::array::from_fn(|k| {
                let anchor = I::point(self.controls[0][j][k]);
                (1..12).fold(anchor, |sum, i| {
                    sum + basis[i] * (I::point(self.controls[i][j][k]) - anchor)
                })
            })
        })
    }
    #[cfg(test)]
    pub(super) fn lattice(&self, divisions: [usize; 2], work: &mut Work<'_>) -> R<Lattice> {
        self.lattice_impl(divisions, None, work)
    }
    pub(super) fn lattice_cached(
        self: &Arc<Self>,
        divisions: [usize; 2],
        old: Option<&Lattice>,
        work: &mut Work<'_>,
    ) -> R<Lattice> {
        work.charge(32)?;
        let old = old.filter(|old| {
            old.owner
                .as_ref()
                .is_some_and(|owner| Arc::ptr_eq(owner, self))
                && old.divisions[1] == divisions[1]
                && old.divisions[0] <= divisions[0]
        });
        let mut result = self.lattice_impl(divisions, old, work)?;
        result.owner = Some(Arc::clone(self));
        Ok(result)
    }
    fn lattice_impl(
        &self,
        divisions: [usize; 2],
        old: Option<&Lattice>,
        work: &mut Work<'_>,
    ) -> R<Lattice> {
        if divisions.iter().any(|&n| !n.is_power_of_two()) {
            return Err(PhillipsReadError::Invalid);
        }
        let count = (divisions[0] + 1)
            .checked_mul(divisions[1] + 1)
            .ok_or(PhillipsReadError::Limit)?;
        let bytes = count
            .checked_mul(std::mem::size_of::<B>())
            .ok_or(PhillipsReadError::Limit)?;
        work.storage(bytes + SCRATCH)?;
        work.charge(count)?;
        let mut points = Vec::new();
        points
            .try_reserve_exact(count)
            .map_err(|_| PhillipsReadError::Limit)?;
        work.storage(points.capacity() * std::mem::size_of::<B>() + SCRATCH)?;
        points.resize(count, [I::point(0.); 3]);
        for u in 0..=divisions[0] {
            work.charge(16)?;
            if let Some(old) = old {
                let ratio = divisions[0] / old.divisions[0];
                if u % ratio == 0 {
                    work.charge(16 * (divisions[1] + 1))?;
                    let previous = (u / ratio) * (divisions[1] + 1);
                    let next = u * (divisions[1] + 1);
                    points[next..next + divisions[1] + 1]
                        .copy_from_slice(&old.points[previous..previous + divisions[1] + 1]);
                    continue;
                }
            }
            work.charge(STRIP_WORK)?;
            let strip = self.strip_at(u as f64 / divisions[0] as f64);
            emit(
                strip,
                divisions[1],
                &mut points[u * (divisions[1] + 1)..(u + 1) * (divisions[1] + 1)],
                work,
            )?;
        }
        Ok(Lattice {
            owner: None,
            points,
            divisions,
        })
    }
    pub(super) fn priority_points(
        &self,
        domain: [[f64; 2]; 2],
        work: &mut Work<'_>,
    ) -> R<[[P; 3]; 3]> {
        work.charge(3 * STRIP_WORK + 9 * 1024)?;
        let values = domain.map(|[a, b]| [a, a + (b - a) * 0.5, b]);
        Ok(std::array::from_fn(|i| {
            let strip = self.strip_at(values[0][i]);
            std::array::from_fn(|j| {
                let mut controls = strip;
                let t = I::point(values[1][j]);
                let complement = I::point(1.) - t;
                for degree in (1..12).rev() {
                    for k in 0..degree {
                        controls[k] = std::array::from_fn(|axis| {
                            complement * controls[k][axis] + t * controls[k + 1][axis]
                        });
                    }
                }
                controls[0].map(I::midpoint)
            })
        }))
    }
    pub(super) fn materialize(
        &self,
        domain: [[f64; 2]; 2],
        work: &mut Work<'_>,
    ) -> R<TensorCell11> {
        work.charge(4096)?;
        let poles = std::array::from_fn(|n| self.controls[n / 12][n % 12].map(I::point));
        let mut cell = crate::tapered_thread::calculation_cell11(poles).map_err(tensor_error)?;
        for (axis, target) in domain.into_iter().enumerate() {
            loop {
                let bounds = cell.domain()[axis];
                if bounds == target {
                    break;
                }
                let mid = bounds[0] + (bounds[1] - bounds[0]) * 0.5;
                let child = usize::from(target[0] >= mid);
                if target[0] < bounds[0]
                    || target[1] > bounds[1]
                    || target[1] - target[0] >= bounds[1] - bounds[0]
                {
                    return Err(PhillipsReadError::Unsupported);
                }
                work.charge(super::SPLIT_WORK)?;
                let [left, right] = cell.subdivide(axis).map_err(tensor_error)?;
                cell = if child == 0 { left } else { right };
            }
        }
        Ok(cell)
    }
}
fn split(mut values: [B; 12]) -> [[B; 12]; 2] {
    let mut result = [[[I::point(0.); 3]; 12]; 2];
    result[0][0] = values[0];
    result[1][11] = values[11];
    for depth in 1..12 {
        for i in 0..12 - depth {
            values[i] = std::array::from_fn(|k| (values[i][k] + values[i + 1][k]) * I::point(0.5));
        }
        result[0][depth] = values[0];
        result[1][11 - depth] = values[11 - depth];
    }
    result
}
fn emit(values: [B; 12], divisions: usize, out: &mut [B], work: &mut Work<'_>) -> R<()> {
    if divisions == 1 {
        work.charge(ENDPOINT_WORK)?;
        out[0] = values[0];
        out[1] = values[11];
        return Ok(());
    }
    work.charge(SPLIT_WORK)?;
    let [left, right] = split(values);
    let half = divisions / 2;
    emit(left, half, &mut out[..=half], work)?;
    emit(right, half, &mut out[half..], work)
}
impl Lattice {
    pub(super) fn bytes(&self) -> usize {
        self.points.capacity() * std::mem::size_of::<B>()
    }
    pub(super) fn corners(&self, domain: [[f64; 2]; 2]) -> R<[B; 4]> {
        let u = domain[0].map(|t| (t * self.divisions[0] as f64) as usize);
        let v = domain[1].map(|t| (t * self.divisions[1] as f64) as usize);
        let at = |u, v| {
            self.points
                .get(u * (self.divisions[1] + 1) + v)
                .copied()
                .ok_or(PhillipsReadError::Unsupported)
        };
        Ok([
            at(u[0], v[0])?,
            at(u[1], v[0])?,
            at(u[1], v[1])?,
            at(u[0], v[1])?,
        ])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn source() -> crate::tapered_thread::RetainedTensor11 {
        crate::tapered_thread::RetainedTensor11::new(
            [1, 1],
            (0..12)
                .flat_map(|i| {
                    (0..12)
                        .map(move |j| [(i * (i - 1)) as f64, (i * j) as f64, (j * (j - 1)) as f64])
                })
                .collect(),
        )
        .unwrap()
    }
    fn work() -> Work<'static> {
        Work {
            remaining: MESH_WORK,
            held: 0,
            peak: 0,
            maximum_bytes: MESH_BYTES,
            maximum_triangles: 1_000_000,
            delegated: Default::default(),
            interrupt: None,
        }
    }
    fn point(u: f64, v: f64) -> P {
        [110. * u * u, 121. * u * v, 110. * v * v]
    }
    #[test]
    fn original_strips_and_midpoint_rows_enclose_independent_tensor_polynomial() {
        let root = Root::new(&source(), [0, 0]).unwrap();
        let mut budget = work();
        let grid = root.lattice([4, 8], &mut budget).unwrap();
        for u in 0..=4 {
            for v in 0..=8 {
                let exact = point(u as f64 / 4., v as f64 / 8.);
                let actual = grid.points[u * 9 + v];
                for k in 0..3 {
                    assert!(actual[k].lo <= exact[k] && exact[k] <= actual[k].hi);
                    assert!(actual[k].width() < 1e-8);
                }
            }
        }
        for (actual, exact) in root.hessian.into_iter().zip([220., 121., 220.]) {
            assert!(actual >= exact && actual < exact * (1. + 1e-12));
        }
        let domain = [[0.25, 0.5], [0.5, 0.625]];
        let expected =
            220. * 0.25f64.powi(2) / 8. + 121. * 0.25 * 0.125 / 4. + 220. * 0.125f64.powi(2) / 8.;
        let error = root.interpolation_error(domain);
        assert!(error >= expected && error < expected * (1. + 1e-12));
        let corners = grid.corners(domain).unwrap();
        for (corner, (u, v)) in
            corners
                .into_iter()
                .zip([(0.25, 0.5), (0.5, 0.5), (0.5, 0.625), (0.25, 0.625)])
        {
            for (actual, exact) in corner.into_iter().zip(point(u, v)) {
                assert!(actual.lo <= exact && exact <= actual.hi);
            }
        }
        let cell = root.materialize(domain, &mut budget).unwrap();
        assert_eq!(cell.domain(), domain);
        let c = domain.map(|d| d[0]);
        let expected = point(c[0], c[1]);
        for (actual, exact) in cell.corner_bounds()[0].into_iter().zip(expected) {
            assert!(actual[0] <= exact && exact <= actual[1]);
        }
        let corners = [
            point(0.25, 0.5),
            point(0.5, 0.5),
            point(0.5, 0.625),
            point(0.25, 0.625),
        ];
        for i in 0..=8 {
            for j in 0..=8 {
                let a = i as f64 / 8.;
                let b = j as f64 / 8.;
                let actual = point(0.25 + 0.25 * a, 0.5 + 0.125 * b);
                let weights = if a >= b {
                    [1. - a, a - b, b, 0.]
                } else {
                    [1. - b, 0., a, b - a]
                };
                let distance = (0..3)
                    .map(|k| {
                        (actual[k] - (0..4).map(|n| weights[n] * corners[n][k]).sum::<f64>()).abs()
                    })
                    .sum::<f64>();
                assert!(distance <= error);
            }
        }
    }
    #[test]
    fn earned_numeric_rows_reuse_only_same_root_and_identical_dyadic_coordinates() {
        let root = Arc::new(Root::new(&source(), [0, 0]).unwrap());
        let old = root.lattice_cached([4, 8], None, &mut work()).unwrap();
        let mut cached_work = work();
        cached_work.held = old.bytes();
        let refined = root
            .lattice_cached([8, 8], Some(&old), &mut cached_work)
            .unwrap();
        let mut fresh_work = work();
        fresh_work.held = old.bytes();
        let fresh = root.lattice_cached([8, 8], None, &mut fresh_work).unwrap();
        assert!(cached_work.remaining > fresh_work.remaining);
        assert!(Arc::ptr_eq(refined.owner.as_ref().unwrap(), &root));
        for u in 0..=8 {
            for v in 0..=8 {
                let actual = refined.points[u * 9 + v];
                for (component, expected) in
                    actual.into_iter().zip(point(u as f64 / 8., v as f64 / 8.))
                {
                    assert!(component.lo <= expected && expected <= component.hi);
                }
                if u % 2 == 0 {
                    for (actual, old) in actual.into_iter().zip(old.points[(u / 2) * 9 + v]) {
                        assert_eq!(
                            [actual.lo.to_bits(), actual.hi.to_bits()],
                            [old.lo.to_bits(), old.hi.to_bits()]
                        );
                    }
                }
            }
        }
        let other = Arc::new(Root::new(&source(), [0, 0]).unwrap());
        let mut other_work = work();
        other_work.held = old.bytes();
        let denied = other
            .lattice_cached([8, 8], Some(&old), &mut other_work)
            .unwrap();
        assert_eq!(other_work.remaining, fresh_work.remaining);
        let bits = |points: &[B]| {
            points
                .iter()
                .flat_map(|p| p.iter().map(|v| [v.lo.to_bits(), v.hi.to_bits()]))
                .collect::<Vec<_>>()
        };
        assert_eq!(bits(&denied.points), bits(&fresh.points));
        let mut poor = work();
        poor.remaining = 1;
        assert!(matches!(
            root.lattice_cached([8, 8], Some(&old), &mut poor),
            Err(PhillipsReadError::Limit)
        ));
        let mut poor = work();
        poor.maximum_bytes = old.bytes();
        poor.held = old.bytes();
        assert!(matches!(
            root.lattice_cached([8, 8], Some(&old), &mut poor),
            Err(PhillipsReadError::Limit)
        ));
        let mut cancelled = work();
        cancelled.interrupt = Some(&|| true);
        assert!(matches!(
            root.lattice_cached([8, 8], Some(&old), &mut cancelled),
            Err(PhillipsReadError::Cancelled)
        ));
        let retry = root
            .lattice_cached([8, 8], Some(&old), &mut work())
            .unwrap();
        assert_eq!(bits(&retry.points), bits(&refined.points));
        assert!(
            root.lattice_cached([8, 16], Some(&old), &mut work())
                .is_ok()
        );
        assert!(root.lattice_cached([2, 8], Some(&old), &mut work()).is_ok());
    }
    #[test]
    fn original_lattice_limits_cancellation_and_retry_keep_source_unchanged() {
        let src = source();
        let root = Root::new(&src, [0, 0]).unwrap();
        assert!(matches!(
            Root::new(&src, [1, 0]),
            Err(PhillipsReadError::Invalid)
        ));
        assert!(matches!(
            root.lattice([3, 2], &mut work()),
            Err(PhillipsReadError::Invalid)
        ));
        assert!(matches!(
            root.lattice(
                [1usize << (usize::BITS - 1), 1usize << (usize::BITS - 1)],
                &mut work()
            ),
            Err(PhillipsReadError::Limit)
        ));
        let mut budget = work();
        budget.maximum_bytes = 1;
        assert!(matches!(
            root.lattice([4, 8], &mut budget),
            Err(PhillipsReadError::Limit)
        ));
        let mut budget = work();
        budget.remaining = 1;
        assert!(matches!(
            root.lattice([4, 8], &mut budget),
            Err(PhillipsReadError::Limit)
        ));
        let mut budget = work();
        budget.interrupt = Some(&|| true);
        assert!(matches!(
            root.lattice([4, 8], &mut budget),
            Err(PhillipsReadError::Cancelled)
        ));
        assert!(matches!(
            root.materialize([[0.25, 0.75], [0., 1.]], &mut work()),
            Err(PhillipsReadError::Unsupported)
        ));
        let before = src.poles().to_vec();
        let grid = root.lattice([4, 8], &mut work()).unwrap();
        assert_eq!(grid.points.len(), 45);
        assert_eq!(src.poles(), before);
    }
}
