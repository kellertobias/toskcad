//! Exact U-restricted auxiliary polynomial, with the ORIGINAL global residual.
//! Mean-value bounds concern Q7, never a guessed derivative of N21-Q7.
use super::*;
pub(crate) struct MeanValueStrip11 {
    source: Arc<ReducedNormal11>,
    center: [[I; 3]; 8],
    derivatives: [[I; 3]; 2],
}
impl ReducedNormal11 {
    // On an already-earned U-only restricted cell: derivative construction672
    // +center contraction1344 interval operations,672 hull comparisons,
    // fixed copies/finite checks/loop/Arc metadata, conservatively <=5120.
    pub(crate) const MEAN_STRIP_WORK: usize = 5120;
    pub(crate) const MEAN_STRIP_BYTES: usize = 32 * 1024;
    // One V point252 interval operations, MVT/residual arithmetic, validation,
    // fixed copies and callback/control margin. Source holders are resident.
    pub(crate) const MEAN_POINT_WORK: usize = 512;
}
impl ReducedNormalCell11 {
    pub(crate) fn mean_value_strip(&self) -> Result<MeanValueStrip11, Error> {
        // V is the ORIGINAL root coordinate. U derivatives below use the
        // already-restricted cell coordinate, hence already contain h_u.
        if self.domain[1] != [0., 1.] {
            return Err(Error::Invalid);
        }
        let grid = &self.controls;
        let mut derivatives = [[I::new(f64::INFINITY, f64::NEG_INFINITY); 3]; 2];
        for (axis, hull) in derivatives.iter_mut().enumerate() {
            for i in 0..if axis == 0 { 7 } else { 8 } {
                for j in 0..if axis == 1 { 7 } else { 8 } {
                    for k in 0..3 {
                        let next = if axis == 0 {
                            grid[i + 1][j][k]
                        } else {
                            grid[i][j + 1][k]
                        };
                        let value = (I::point(7.) * (next - grid[i][j][k])).checked()?;
                        hull[k].lo = hull[k].lo.min(value.lo);
                        hull[k].hi = hull[k].hi.max(value.hi);
                    }
                }
            }
        }
        let mut rows = self.controls;
        for depth in 1..8 {
            for i in 0..8 - depth {
                let (left, right) = rows.split_at_mut(i + 1);
                for (value, next) in left[i].iter_mut().zip(&right[0]) {
                    *value = std::array::from_fn(|k| (value[k] + next[k]) * I::point(0.5));
                }
            }
        }
        if rows[0].iter().flatten().any(|v| !v.valid()) {
            return Err(Error::Unresolved);
        }
        Ok(MeanValueStrip11 {
            source: Arc::clone(&self.source),
            center: rows[0],
            derivatives,
        })
    }
}
impl MeanValueStrip11 {
    pub(crate) fn bounds(&self, v: [f64; 2]) -> Result<[[f64; 2]; 3], Error> {
        if !NormalPolynomial11::valid([[0., 1.], v]) {
            return Err(Error::Invalid);
        }
        let middle = v[0] + (v[1] - v[0]) * 0.5;
        let t = I::point(middle);
        let s = I::point(1.) - t;
        let mut row = self.center;
        for depth in 1..8 {
            for j in 0..8 - depth {
                row[j] = std::array::from_fn(|k| s * row[j][k] + t * row[j + 1][k]);
            }
        }
        // U derivative belongs to the normalized RESTRICTED row. Its domain
        // is [0,1], so do not multiply the original U width a second time.
        let offsets = [I::new(-0.5, 0.5), I::new(v[0], v[1]) - t];
        let result = std::array::from_fn::<_, 3, _>(|k| {
            row[0][k]
                + self.derivatives[0][k] * offsets[0]
                + self.derivatives[1][k] * offsets[1]
                + self.source.residual[k]
        });
        for value in result {
            value.checked()?;
        }
        Ok(result.map(|v| [v.lo, v.hi]))
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn restricted_mean_value_strip_encloses_independent_bilinear_normals_and_keeps_residual() {
        let tensor = RetainedTensor11::new(
            [1, 1],
            (0..12)
                .flat_map(|i| (0..12).map(move |j| [i as f64, j as f64, (i * j) as f64]))
                .collect(),
        )
        .unwrap();
        let normal = Arc::new(NormalPolynomial11::new(&tensor.cell([0, 0]).unwrap()).unwrap());
        let field = Arc::new(ReducedNormal11::new(normal).unwrap());
        for ud in [[0., 1.], [0.25, 0.375], [0.875, 1.]] {
            let mut row = field.root_cell();
            while row.domain()[0] != ud {
                let [a, b] = row.domain()[0];
                let child = usize::from(ud[0] >= a + (b - a) * 0.5);
                row = row.subdivide(0).unwrap()[child].clone();
            }
            let strip = row.mean_value_strip().unwrap();
            assert!(Arc::ptr_eq(&strip.source, &field));
            assert_eq!(
                strip.source.coefficient_residual_bounds(),
                field.coefficient_residual_bounds()
            );
            for vd in [[0., 1.], [0.125, 0.25], [0.75, 0.78125]] {
                let bounds = strip.bounds(vd).unwrap();
                for i in 0..=8 {
                    for j in 0..=8 {
                        let u = ud[0] + (ud[1] - ud[0]) * i as f64 / 8.;
                        let v = vd[0] + (vd[1] - vd[0]) * j as f64 / 8.;
                        for (actual, expected) in
                            bounds.into_iter().zip([-1331. * v, -1331. * u, 121.])
                        {
                            assert!(
                                actual[0] <= expected && expected <= actual[1],
                                "{actual:?} {expected}"
                            );
                        }
                    }
                }
            }
            assert_eq!(strip.bounds([0., f64::NAN]), Err(Error::Invalid));
            assert_eq!(strip.bounds([1., 0.]), Err(Error::Invalid));
            assert!(strip.bounds([0., 1.]).is_ok());
        }
        assert!(matches!(
            field.root_cell().subdivide(1).unwrap()[0].mean_value_strip(),
            Err(Error::Invalid)
        ));
        assert!(field.root_cell().mean_value_strip().is_ok());
    }
}
