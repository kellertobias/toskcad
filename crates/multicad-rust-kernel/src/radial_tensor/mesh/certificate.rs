//! Whole-patch source-to-facet certificates, including emitted vertex rounding.
//! Mesh UV is an auxiliary piecewise-affine chart over the unchanged nets.
use super::*;
type Points = Vec<[I; 3]>;
type Grid = Vec<Points>;
#[derive(Clone)]
pub(super) struct Affine {
    origin: [I; 3],
    u: [I; 3],
    v: [I; 3],
}
impl Affine {
    pub fn new(uv: [[f64; 2]; 3], vertices: [Point; 3]) -> Result<Self, Error> {
        let u = [1, 2].map(|i| I::point(uv[i][0]) - I::point(uv[0][0]));
        let v = [1, 2].map(|i| I::point(uv[i][1]) - I::point(uv[0][1]));
        let determinant = (u[0] * v[1] - u[1] * v[0]).checked()?;
        if determinant.lo <= 0. && determinant.hi >= 0. {
            return Err(Error::Unresolved);
        }
        let edges = [1, 2].map(|i| {
            std::array::from_fn::<_, 3, _>(|k| I::point(vertices[i][k]) - I::point(vertices[0][k]))
        });
        let du = std::array::from_fn(|k| (edges[0][k] * v[1] - edges[1][k] * v[0]) / determinant);
        let dv = std::array::from_fn(|k| (u[0] * edges[1][k] - u[1] * edges[0][k]) / determinant);
        let origin = std::array::from_fn(|k| {
            I::point(vertices[0][k]) - du[k] * I::point(uv[0][0]) - dv[k] * I::point(uv[0][1])
        });
        for value in origin.into_iter().chain(du).chain(dv) {
            value.checked()?;
        }
        Ok(Self {
            origin,
            u: du,
            v: dv,
        })
    }
    pub fn line(range: [f64; 2], points: [Point; 2]) -> Result<Self, Error> {
        let width = (I::point(range[1]) - I::point(range[0])).positive()?;
        let u = std::array::from_fn(|k| (I::point(points[1][k]) - I::point(points[0][k])) / width);
        let origin = std::array::from_fn(|k| I::point(points[0][k]) - u[k] * I::point(range[0]));
        for value in origin.into_iter().chain(u) {
            value.checked()?;
        }
        Ok(Self {
            origin,
            u,
            v: [I::point(0.); 3],
        })
    }
    fn at(&self, u: I, v: I) -> [I; 3] {
        std::array::from_fn(|k| self.origin[k] + self.u[k] * u + self.v[k] * v)
    }
}
/// The coefficients of P-L enclose the difference to the actual emitted
/// triangle's affine interpolant. Restriction is stronger than necessary:
/// it covers the entire parameter bounding rectangle, including the triangle.
pub(super) fn chord_error(
    poles: &Grid,
    domain: [[f64; 2]; 2],
    range: [[f64; 2]; 2],
    affine: &Affine,
    work: &mut Work<'_>,
) -> Result<f64, Error> {
    work.charge(poles.len() * poles[0].len())?;
    let degree = [poles.len() - 1, poles[0].len() - 1];
    let mut error = poles.clone();
    for (i, row) in error.iter_mut().enumerate() {
        let u = I::point(domain[0][0])
            + (I::point(domain[0][1]) - I::point(domain[0][0])) * I::point(i as f64)
                / I::point(degree[0] as f64);
        for (j, p) in row.iter_mut().enumerate() {
            let v = I::point(domain[1][0])
                + (I::point(domain[1][1]) - I::point(domain[1][0])) * I::point(j as f64)
                    / I::point(degree[1] as f64);
            let linear = affine.at(u, v);
            *p = std::array::from_fn(|k| p[k] - linear[k]);
        }
    }
    let restricted = restrict(error, local_range(domain, range)?, work)?;
    let mut upper: f64 = 0.;
    for point in restricted.iter().flatten() {
        let l1 = point
            .iter()
            .fold(I::point(0.), |s, v| s + I::point(v.magnitude()))
            .checked()?;
        upper = upper.max(l1.hi);
    }
    Ok(upper)
}
/// One restriction of the immutable source polynomial. Every actual facet
/// still subtracts its own affine polynomial in the retained global chart.
pub(super) struct SourceCell {
    poles: Grid,
    domain: [[I; 2]; 2],
}
/// A restriction in U only, owned by one immutable source patch.
/// Its V coordinate still uses the original patch chart.
pub(super) struct SourceStrip {
    poles: Grid,
    v_domain: [f64; 2],
    u: [I; 2],
    next_v: f64,
}
fn transfer(domain: [f64; 2], range: [f64; 2]) -> [I; 2] {
    range.map(|t| I::point(domain[0]) + (I::point(domain[1]) - I::point(domain[0])) * I::point(t))
}
impl SourceStrip {
    pub fn new(
        poles: &Grid,
        domain: [[f64; 2]; 2],
        u: [f64; 2],
        work: &mut Work<'_>,
    ) -> Result<Self, Error> {
        work.charge(poles.len() * poles[0].len())?;
        let local = local_range(domain, [u, domain[1]])?[0];
        let poles = restrict_axis(poles.clone(), 0, local, work)?;
        Ok(Self {
            poles,
            v_domain: domain[1],
            u: transfer(domain[0], local),
            next_v: domain[1][0],
        })
    }
    pub fn consecutive(&self, v: [f64; 2]) -> bool {
        v[0].max(self.v_domain[0]).to_bits() == self.next_v.to_bits()
    }
    /// The retained tail is the ORIGINAL U-restricted position polynomial
    /// reparameterized at an exact nominal global V boundary. Directed ratios
    /// enclose that exact change of variable; no sampled source is substituted.
    pub fn continued_cell(
        &mut self,
        v: [f64; 2],
        work: &mut Work<'_>,
    ) -> Result<SourceCell, Error> {
        // Tensor position nets are exactly10×10. The bounded continuation
        // workspace covers old/current/two split nets, row Vec capacities and
        // de Casteljau temporaries. Other degrees keep the stateless path.
        if self.poles.len() > 10 || self.poles.iter().any(|row| row.len() > 10) {
            return Err(Error::Limit);
        }
        let start = v[0].max(self.v_domain[0]);
        let end = v[1].min(self.v_domain[1]);
        if !start.is_finite() || !end.is_finite() || start < self.next_v || end <= start {
            return Err(Error::Unresolved);
        }
        work.charge(self.poles.len() * self.poles[0].len())?;
        let mut poles = self.poles.clone();
        if start > self.next_v {
            let t = ((I::point(start) - I::point(self.next_v))
                / (I::point(self.v_domain[1]) - I::point(self.next_v)).positive()?)
            .checked()?;
            poles = split_v(poles, t, work)?.1;
        }
        let cell = if end == self.v_domain[1] {
            // The complete remaining tail already has this exact nominal
            // chart. An endpoint split would be the identity, so do no split.
            poles
        } else {
            let t = ((I::point(end) - I::point(start))
                / (I::point(self.v_domain[1]) - I::point(start)).positive()?)
            .checked()?;
            let (cell, tail) = split_v(poles, t, work)?;
            self.poles = tail;
            cell
        };
        self.next_v = end;
        Ok(SourceCell {
            poles: cell,
            domain: [self.u, [I::point(start), I::point(end)]],
        })
    }
    pub fn normal_coefficients(&self, v: [f64; 2], work: &mut Work<'_>) -> Result<Grid, Error> {
        Ok(self.cell(v, work)?.poles)
    }
    pub fn cell(&self, v: [f64; 2], work: &mut Work<'_>) -> Result<SourceCell, Error> {
        // A continued position tail must never masquerade as the original
        // full-V net, including through normal/derivative helper misuse.
        if self.next_v.to_bits() != self.v_domain[0].to_bits() {
            return Err(Error::Invalid);
        }
        work.charge(self.poles.len() * self.poles[0].len())?;
        let local = local_range([[0., 1.], self.v_domain], [[0., 1.], v])?[1];
        let poles = restrict_axis(self.poles.clone(), 1, local, work)?;
        Ok(SourceCell {
            poles,
            domain: [self.u, transfer(self.v_domain, local)],
        })
    }
}
fn split_v(poles: Grid, t: I, work: &mut Work<'_>) -> Result<(Grid, Grid), Error> {
    let mut left = Vec::with_capacity(poles.len());
    let mut right = Vec::with_capacity(poles.len());
    for row in poles {
        let (a, b) = split(&row, t, work)?;
        left.push(a);
        right.push(b);
    }
    Ok((left, right))
}
impl SourceCell {
    pub fn new(
        poles: &Grid,
        domain: [[f64; 2]; 2],
        range: [[f64; 2]; 2],
        work: &mut Work<'_>,
    ) -> Result<Self, Error> {
        work.charge(poles.len() * poles[0].len())?;
        let local = local_range(domain, range)?;
        let source = restrict(poles.clone(), local, work)?;
        // local_range rounds outwards. Use those actual retained restriction
        // endpoints, including the original chart's interval transfer, rather
        // than replacing them with the requested rectangle's rounded bounds.
        let domain = std::array::from_fn(|axis| transfer(domain[axis], local[axis]));
        Ok(Self {
            poles: source,
            domain,
        })
    }
    pub fn chord_error(&self, affine: &Affine, work: &mut Work<'_>) -> Result<f64, Error> {
        work.charge(self.poles.len() * self.poles[0].len())?;
        let degree = [self.poles.len() - 1, self.poles[0].len() - 1];
        let mut upper = 0_f64;
        for (i, row) in self.poles.iter().enumerate() {
            let u = self.domain[0][0]
                + (self.domain[0][1] - self.domain[0][0]) * I::point(i as f64)
                    / I::point(degree[0] as f64);
            for (j, point) in row.iter().enumerate() {
                let v = self.domain[1][0]
                    + (self.domain[1][1] - self.domain[1][0]) * I::point(j as f64)
                        / I::point(degree[1] as f64);
                let line = affine.at(u, v);
                let l1 = (0..3)
                    .fold(I::point(0.), |s, k| {
                        s + I::point((point[k] - line[k]).magnitude())
                    })
                    .checked()?;
                upper = upper.max(l1.hi);
            }
        }
        Ok(upper)
    }
}
/// A convex positive normal cone proves the angle to the actual rounded
/// facet normal. No secant-normal theorem or sampled normal is assumed.
/// Restricting the original normal polynomial suffices for directions:
/// affine UV Jacobians are strictly positive scalars and do not change
/// this cone. These unscaled coefficients are never used as area values.
pub(super) fn normal_matches(
    normals: Grid,
    domain: [[f64; 2]; 2],
    range: [[f64; 2]; 2],
    facet: [I; 3],
    angle: f64,
    work: &mut Work<'_>,
) -> Result<bool, Error> {
    let normals = restrict(normals, local_range(domain, range)?, work)?;
    normal_coefficients_match(&normals, facet, angle, work)
}
pub(super) fn normal_coefficients_match(
    normals: &Grid,
    facet: [I; 3],
    angle: f64,
    work: &mut Work<'_>,
) -> Result<bool, Error> {
    // For positive projection the actual angle is atan(cross/dot),
    // which is <=cross/dot. This sufficient bound preserves arbitrarily
    // small requests without a separate absolute-libm preparation error.
    let tangent = if angle >= std::f64::consts::FRAC_PI_2 {
        None
    } else {
        Some(I::point(angle))
    };
    for &normal in normals.iter().flatten() {
        work.step()?;
        let projection = dot(facet, normal).checked()?;
        if projection.lo <= 0. {
            return Ok(false);
        }
        if let Some(tangent) = tangent {
            let error = cross(facet, normal)
                .into_iter()
                .fold(I::point(0.), |s, v| s + v.square())
                .sqrt()
                .checked()?;
            if error.hi > (tangent * projection).lo {
                return Ok(false);
            }
        }
    }
    Ok(true)
}
pub(super) fn facet_normal(points: [Point; 3]) -> [I; 3] {
    let edge =
        [1, 2].map(|i| std::array::from_fn(|k| I::point(points[i][k]) - I::point(points[0][k])));
    cross(edge[0], edge[1])
}
fn cross(a: [I; 3], b: [I; 3]) -> [I; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn local_range(domain: [[f64; 2]; 2], range: [[f64; 2]; 2]) -> Result<[[f64; 2]; 2], Error> {
    let mut result = [[0.; 2]; 2];
    for axis in 0..2 {
        let width = (I::point(domain[axis][1]) - I::point(domain[axis][0])).positive()?;
        result[axis] = [
            ((I::point(range[axis][0]) - I::point(domain[axis][0])) / width)
                .checked()?
                .lo
                .max(0.),
            ((I::point(range[axis][1]) - I::point(domain[axis][0])) / width)
                .checked()?
                .hi
                .min(1.),
        ];
        if result[axis][0] > result[axis][1] {
            return Err(Error::Invalid);
        }
    }
    Ok(result)
}
fn restrict(mut grid: Grid, range: [[f64; 2]; 2], work: &mut Work<'_>) -> Result<Grid, Error> {
    for (axis, &range) in range.iter().enumerate() {
        grid = restrict_axis(grid, axis, range, work)?;
    }
    Ok(grid)
}
fn restrict_axis(
    mut grid: Grid,
    axis: usize,
    [a, b]: [f64; 2],
    work: &mut Work<'_>,
) -> Result<Grid, Error> {
    if axis == 0 {
        grid = transpose(grid);
    }
    for row in &mut grid {
        if b < 1. {
            *row = split(row, I::point(b), work)?.0;
        }
        if a > 0. {
            *row = split(row, (I::point(a) / I::point(b)).checked()?, work)?.1;
        }
    }
    if axis == 0 {
        grid = transpose(grid);
    }
    Ok(grid)
}
fn transpose(grid: Grid) -> Grid {
    (0..grid[0].len())
        .map(|j| grid.iter().map(|row| row[j]).collect())
        .collect()
}
pub(super) fn split(
    points: &[[I; 3]],
    t: I,
    work: &mut Work<'_>,
) -> Result<(Points, Points), Error> {
    work.charge(points.len() * (points.len() - 1) / 2)?;
    let mut row = points.to_vec();
    let mut left = vec![row[0]];
    let mut right = vec![row[row.len() - 1]];
    let one_minus = I::point(1.) - t;
    while row.len() > 1 {
        row = row
            .windows(2)
            .map(|p| std::array::from_fn(|k| one_minus * p[0][k] + t * p[1][k]))
            .collect();
        left.push(row[0]);
        right.push(row[row.len() - 1]);
    }
    right.reverse();
    Ok((left, right))
}

#[cfg(test)]
mod continuation_tests {
    use super::*;
    fn work() -> Work<'static> {
        Work {
            remaining: 10_000_000,
            interrupt: None,
        }
    }
    fn contains(value: I, expected: f64) {
        assert!(
            value.lo <= expected && expected <= value.hi,
            "{value:?} misses {expected}"
        );
    }
    fn polynomial() -> Grid {
        (0..3)
            .map(|i| {
                (0..3)
                    .map(|j| {
                        let u = i as f64 / 2.;
                        let v = j as f64 / 2.;
                        [
                            I::point(u),
                            I::point(v),
                            I::point(u * v + if j == 2 { 1. } else { 0. }),
                        ]
                    })
                    .collect()
            })
            .collect()
    }
    fn evaluate(grid: &Grid, t: [I; 2]) -> [I; 3] {
        let row_at = |points: &[[I; 3]], t: I| {
            let mut p = points.to_vec();
            while p.len() > 1 {
                p = p
                    .windows(2)
                    .map(|p| std::array::from_fn(|k| (I::point(1.) - t) * p[0][k] + t * p[1][k]))
                    .collect();
            }
            p[0]
        };
        let rows: Vec<_> = grid.iter().map(|r| row_at(r, t[1])).collect();
        row_at(&rows, t[0])
    }
    fn derivative(grid: &Grid, axis: usize, width: I) -> Grid {
        let degree = [grid.len() - 1, grid[0].len() - 1];
        (0..grid.len() - usize::from(axis == 0))
            .map(|i| {
                (0..grid[0].len() - usize::from(axis == 1))
                    .map(|j| {
                        std::array::from_fn(|k| {
                            (grid[i + usize::from(axis == 0)][j + usize::from(axis == 1)][k]
                                - grid[i][j][k])
                                * I::point(degree[axis] as f64)
                                / width
                        })
                    })
                    .collect()
            })
            .collect()
    }
    #[test]
    fn continued_original_polynomial_coefficients_and_global_derivatives_have_independent_oracles()
    {
        let mut w = work();
        let mut source = SourceStrip::new(&polynomial(), [[0., 1.]; 2], [0., 1.], &mut w).unwrap();
        for [a, b] in [
            [0., 0.125],
            [0.125, 0.375],
            [0.375, 0.5],
            [0.5, 0.875],
            [0.875, 1.],
        ] {
            assert!(source.consecutive([a, b]));
            let cell = source.continued_cell([a, b], &mut w).unwrap();
            assert_eq!(cell.domain[1].map(|x| [x.lo, x.hi]), [[a, a], [b, b]]);
            // Exact dyadic Bernstein coefficients of ORIGINAL P=(u,v,uv+v²),
            // independently expanded in global coordinates for the entire cell.
            for i in 0..3 {
                for j in 0..3 {
                    let u = i as f64 / 2.;
                    let v = a + (b - a) * j as f64 / 2.;
                    let vv = [a * a, a * b, b * b][j];
                    for (bound, expected) in cell.poles[i][j].into_iter().zip([u, v, u * v + vv]) {
                        contains(bound, expected);
                    }
                }
            }
            let widths = cell.domain.map(|d| d[1] - d[0]);
            let du = derivative(&cell.poles, 0, widths[0]);
            let dv = derivative(&cell.poles, 1, widths[1]);
            let mixed = derivative(&du, 1, widths[1]);
            for u in [0., 0.25, 0.5, 1.] {
                for v in [a, (a + b) / 2., b] {
                    let t = std::array::from_fn(|k| {
                        (I::point([u, v][k]) - cell.domain[k][0]) / widths[k]
                    });
                    for (grid, expected) in [
                        (&cell.poles, [u, v, u * v + v * v]),
                        (&du, [1., 0., v]),
                        (&dv, [0., 1., u + 2. * v]),
                        (&mixed, [0., 0., 1.]),
                    ] {
                        for (bound, expected) in evaluate(grid, t).into_iter().zip(expected) {
                            contains(bound, expected);
                        }
                    }
                }
            }
        }
        assert_eq!(source.cell([0., 0.125], &mut w).err(), Some(Error::Invalid));
        assert_eq!(
            source.normal_coefficients([0., 0.125], &mut w).err(),
            Some(Error::Invalid)
        );
        assert!(!source.consecutive([0., 0.125]));
        assert_eq!(
            source.continued_cell([0., 0.125], &mut w).err(),
            Some(Error::Unresolved)
        );
        assert_eq!(
            source.continued_cell([1., 1.], &mut w).err(),
            Some(Error::Unresolved)
        );
    }
}
