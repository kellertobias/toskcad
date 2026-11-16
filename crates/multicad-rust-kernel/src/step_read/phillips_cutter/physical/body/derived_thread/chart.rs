//! Inverse of the actual four-column proposed chart, without rigid projection.
use super::*;
fn cross(a: [I; 3], b: [I; 3]) -> [I; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dot(a: [I; 3], b: [I; 3]) -> I {
    a.into_iter()
        .zip(b)
        .fold(I::point(0.), |s, (a, b)| s + a * b)
}
pub(super) struct Chart {
    chart: CommonMouthChart,
    inverse: [[I; 3]; 3],
}
impl Chart {
    pub fn new(head: &ClosedRepairedHead<'_, '_, '_, '_, '_>) -> Self {
        let chart = head.chart();
        // These are exactly the four immutable columns retained by
        // branches/derived/closed.rs, including its independently stored N.
        // The inherited Gram1e-10+rounding proof earns det>1-1e-9 and bounded
        // entries; the identical directed determinant below is positive, and
        // every cofactor/division is finite. No normalization is performed.
        let [_, u, v, n] = chart.columns.map(|p| p.map(I::point));
        let det = dot(u, cross(v, n));
        let inverse = [cross(v, n), cross(n, u), cross(u, v)].map(|r| r.map(|x| x / det));
        Self { chart, inverse }
    }
    pub fn local(&self, world: Point) -> Result<[I; 3]> {
        let d = std::array::from_fn(|i| {
            I::new(world[i].next_down(), world[i].next_up()) - I::point(self.chart.columns[0][i])
        });
        let mut p = self.inverse.map(|r| dot(r, d));
        p[2] = p[2] + I::point(self.chart.outer_top_z);
        require(
            p.iter()
                .all(|x| x.valid() && x.lo.is_finite() && x.hi.is_finite() && x.magnitude() < 1e7),
        )?;
        Ok(p)
    }
    pub fn world(&self, p: [I; 3]) -> [I; 3] {
        std::array::from_fn(|i| {
            I::point(self.chart.columns[0][i])
                + I::point(self.chart.columns[1][i]) * p[0]
                + I::point(self.chart.columns[2][i]) * p[1]
                + I::point(self.chart.columns[3][i]) * (p[2] - I::point(self.chart.outer_top_z))
        })
    }
}
