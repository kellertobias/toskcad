use super::*;
fn p(x: f64) -> I {
    I::point(x)
}
fn power(x: I, n: usize) -> I {
    (0..n).fold(p(1.), |a, _| a * x)
}
fn norm(x: [I; 3]) -> I {
    x.into_iter().fold(p(0.), |a, b| a + p(b.magnitude()))
}
pub(super) struct Candidate {
    pub ideal: [[[I; 3]; 12]; 12],
    domain: [[f64; 2]; 2],
    remainder: [f64; 3],
}
pub(super) type ScaledBasis = [[[I; 12]; 12]; 2];
pub(super) fn scaled_basis(
    spans: [usize; 2],
    basis: &curves::Basis,
    work: &mut Work<'_>,
) -> Result<ScaledBasis, Error> {
    // The sealed source-derived uniform grid has positive dyadic spans.
    //1872 interval operations plus fixed coefficient publication fit2560.
    work.charge(2560)?;
    Ok(std::array::from_fn(|axis| {
        std::array::from_fn(|i| {
            std::array::from_fn(|k| basis[i][k] * power(p(1. / spans[axis] as f64), k))
        })
    }))
}
pub(super) fn cell(
    law: &StockClipPatchLaw<'_>,
    domain: [[f64; 2]; 2],
    scaled: &ScaledBasis,
    work: &mut Work<'_>,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Candidate, Error> {
    let center = domain.map(|d| d[0].midpoint(d[1]));
    let local = law
        .taylor(
            [center[0], center[0]],
            [center[1], center[1]],
            source_options(work)?,
            interrupt,
        )
        .map_err(source_error)?;
    work.charge(local.work_used())?;
    let whole = law
        .taylor(domain[0], domain[1], source_options(work)?, interrupt)
        .map_err(source_error)?;
    work.charge(whole.work_used())?;
    // Two separable passes use20736 interval operations; fixed arrays and
    // coefficient publication fit1792 further units. The uniform scaled
    // bases are prepared ONCE per patch, not recomputed for every cell.
    work.charge(22_528)?;
    let h = domain.map(|d| d[1] - d[0]);
    let mut first = [[[p(0.); 3]; 12]; 12];
    for (i, row) in first.iter_mut().enumerate() {
        for (n, out) in row.iter_mut().enumerate() {
            for (m, scale) in scaled[0][i].iter().enumerate() {
                for (axis, value) in out.iter_mut().enumerate() {
                    let b = local.coefficients_mm()[m][n][axis];
                    *value = *value + I::new(b[0], b[1]) * *scale;
                }
            }
        }
    }
    let mut ideal = [[[p(0.); 3]; 12]; 12];
    for (i, row) in ideal.iter_mut().enumerate() {
        for (j, out) in row.iter_mut().enumerate() {
            for n in 0..12 {
                for (axis, value) in out.iter_mut().enumerate() {
                    *value = *value + first[i][n][axis] * scaled[1][j][n];
                }
            }
        }
    }
    let coefficient = |i: usize, j: usize| p(whole.coefficient_l1_bounds_mm()[i][j]);
    let a = p(h[0] / 2.);
    let b = p(h[1] / 2.);
    // Iterated one-variable Taylor theorem. The second sum retains ALL
    // mixed derivatives; each whole rectangle query encloses every required
    // intermediate center. For d/du use degree10 in u; d/dv is symmetric.
    let position = coefficient(12, 0) * power(a, 12)
        + (0..12).fold(p(0.), |s, i| {
            s + coefficient(i, 12) * power(a, i) * power(b, 12)
        });
    let du = p(12.) * coefficient(12, 0) * power(a, 11)
        + (0..11).fold(p(0.), |s, i| {
            s + p((i + 1) as f64) * coefficient(i + 1, 12) * power(a, i) * power(b, 12)
        });
    let dv = p(12.) * coefficient(0, 12) * power(b, 11)
        + (0..11).fold(p(0.), |s, j| {
            s + p((j + 1) as f64) * coefficient(12, j + 1) * power(a, 12) * power(b, j)
        });
    Ok(Candidate {
        ideal,
        domain,
        remainder: [position.hi, du.hi, dv.hi],
    })
}
impl Candidate {
    pub fn certify(
        self,
        poles: [[[f64; 3]; 12]; 12],
        work: &mut Work<'_>,
    ) -> Result<SurfaceCellCertificate, Error> {
        work.charge(6144)?;
        let errors: [[[I; 3]; 12]; 12] = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                std::array::from_fn(|axis| {
                    I::new(poles[i][j][axis].next_down(), poles[i][j][axis].next_up())
                        - self.ideal[i][j][axis]
                })
            })
        });
        let position = errors
            .into_iter()
            .flatten()
            .map(norm)
            .fold(0f64, |a, b| a.max(b.hi));
        let mut derivative = [0f64; 2];
        for i in 0..12 {
            for j in 0..12 {
                if i < 11 {
                    derivative[0] = derivative[0].max(
                        norm(std::array::from_fn(|axis| {
                            errors[i + 1][j][axis] - errors[i][j][axis]
                        }))
                        .hi,
                    );
                }
                if j < 11 {
                    derivative[1] = derivative[1].max(
                        norm(std::array::from_fn(|axis| {
                            errors[i][j + 1][axis] - errors[i][j][axis]
                        }))
                        .hi,
                    );
                }
            }
        }
        let derivative_error = std::array::from_fn(|axis| {
            (p(11.) * p(derivative[axis]) / p(self.domain[axis][1] - self.domain[axis][0])
                + p(self.remainder[axis + 1]))
            .hi
        });
        Ok(SurfaceCellCertificate {
            domain: self.domain,
            position_error: (p(position) + p(self.remainder[0])).hi,
            derivative_error,
        })
    }
}
