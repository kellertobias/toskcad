//! Actual source spline parameterizations; no inverse or sampled fitting.
use super::*;
#[derive(Clone)]
pub(super) struct Spline {
    pub controls: Vec<[I; 3]>,
    pub knots: Vec<f64>,
}
impl Spline {
    pub fn new(controls: Vec<[I; 3]>) -> Self {
        let spans = (controls.len() - 1) / 9;
        let knots = (0..=spans).map(|i| i as f64 / spans as f64).collect();
        Self { controls, knots }
    }
    pub fn axial_error(&self, height: f64, work: &mut Work<'_>) -> Result<f64, Error> {
        let mut error: f64 = 0.;
        for span in 0..self.knots.len() - 1 {
            work.charge(10)?;
            for j in 0..=9 {
                let v = I::point(self.knots[span])
                    + (I::point(self.knots[span + 1]) - I::point(self.knots[span]))
                        * I::point(j as f64)
                        / I::point(9.);
                error = error.max(
                    (self.controls[9 * span + j][2] - I::point(height) * v)
                        .checked()?
                        .magnitude(),
                );
            }
        }
        Ok(error)
    }
    pub fn chord_error(
        &self,
        range: [f64; 2],
        points: [Point; 2],
        work: &mut Work<'_>,
    ) -> Result<f64, Error> {
        let affine = certificate::Affine::line(range, points)?;
        let first = self
            .knots
            .partition_point(|&k| k < range[0])
            .saturating_sub(1);
        let last = self
            .knots
            .partition_point(|&k| k <= range[1])
            .min(self.knots.len() - 1);
        let mut error: f64 = 0.;
        for span in first..last {
            let a = range[0].max(self.knots[span]);
            let b = range[1].min(self.knots[span + 1]);
            if a >= b {
                continue;
            }
            let grid = self.controls[9 * span..=9 * span + 9]
                .iter()
                .map(|&p| vec![p; 2])
                .collect();
            error = error.max(certificate::chord_error(
                &grid,
                [[self.knots[span], self.knots[span + 1]], [0., 1.]],
                [[a, b], [0., 0.]],
                &affine,
                work,
            )?);
        }
        Ok(error)
    }
    pub fn bounds(&self, range: [f64; 2], work: &mut Work<'_>) -> Result<[[I; 3]; 2], Error> {
        let values = self.derivatives(range, work)?;
        Ok([values[0], values[1]])
    }
    /// One-sided regular derivatives; C0 first-derivative jumps are
    /// retained independently by angular_jumps, never assumed absent.
    pub fn derivatives(&self, range: [f64; 2], work: &mut Work<'_>) -> Result<[[I; 3]; 3], Error> {
        let first = self
            .knots
            .partition_point(|&k| k < range[0])
            .saturating_sub(1);
        let last = self
            .knots
            .partition_point(|&k| k <= range[1])
            .min(self.knots.len() - 1);
        let mut bounds = [[I::new(f64::INFINITY, f64::NEG_INFINITY); 3]; 3];
        for span in first..last {
            work.charge(32)?;
            let width = (I::point(self.knots[span + 1]) - I::point(self.knots[span])).positive()?;
            let a = ((I::point(range[0]) - I::point(self.knots[span])) / width)
                .lo
                .max(0.);
            let b = ((I::point(range[1]) - I::point(self.knots[span])) / width)
                .hi
                .min(1.);
            if a > b {
                continue;
            }
            let points = &self.controls[9 * span..=9 * span + 9];
            let derivative = points
                .windows(2)
                .map(|p| std::array::from_fn(|k| (p[1][k] - p[0][k]) * I::point(9.) / width))
                .collect::<Vec<_>>();
            let second = derivative
                .windows(2)
                .map(|p| std::array::from_fn(|k| (p[1][k] - p[0][k]) * I::point(8.) / width))
                .collect::<Vec<_>>();
            for (order, mut controls) in [points.to_vec(), derivative, second]
                .into_iter()
                .enumerate()
            {
                if b < 1. {
                    controls = certificate::split(&controls, I::point(b), work)?.0;
                }
                if a > 0. {
                    controls = certificate::split(
                        &controls,
                        (I::point(a) / I::point(b)).checked()?,
                        work,
                    )?
                    .1;
                }
                for point in controls {
                    for (k, v) in point.into_iter().enumerate() {
                        bounds[order][k].lo = bounds[order][k].lo.min(v.lo);
                        bounds[order][k].hi = bounds[order][k].hi.max(v.hi);
                    }
                }
            }
        }
        for v in bounds.into_iter().flatten() {
            v.checked()?;
        }
        Ok(bounds)
    }
    pub fn angular_jumps(&self, range: [f64; 2], work: &mut Work<'_>) -> Result<I, Error> {
        let first = self.knots.partition_point(|&k| k <= range[0]).max(1);
        let last = self
            .knots
            .partition_point(|&k| k < range[1])
            .min(self.knots.len() - 1);
        let mut total = I::point(0.);
        for knot in first..last {
            work.charge(64)?;
            let p = self.controls[9 * knot];
            let left = (I::point(self.knots[knot]) - I::point(self.knots[knot - 1])).positive()?;
            let right = (I::point(self.knots[knot + 1]) - I::point(self.knots[knot])).positive()?;
            let jump = std::array::from_fn::<_, 3, _>(|k| {
                I::point(9.)
                    * ((self.controls[9 * knot + 1][k] - p[k]) / right
                        - (p[k] - self.controls[9 * knot - 1][k]) / left)
            });
            let radius = (p[0].square() + p[1].square()).positive()?;
            let angular = ((p[0] * jump[1] - p[1] * jump[0]) / radius).checked()?;
            total = total + I::point(angular.magnitude());
        }
        total.checked()
    }
}
