//! Actual four-column construction map followed by the occurrence map.
//! Neither map is normalized or substituted by a rigid nominal pose.
use super::*;
mod inverse;
type Point = [f64; 3];
type BoxPoint = [I; 3];
fn dot(a: BoxPoint, b: BoxPoint) -> I {
    (0..3).fold(I::ZERO, |s, i| s + a[i] * b[i])
}
fn cross(a: BoxPoint, b: BoxPoint) -> BoxPoint {
    std::array::from_fn(|i| a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3])
}
fn difference(a: Point, b: Point) -> BoxPoint {
    std::array::from_fn(|i| I::point(a[i]) - I::point(b[i]))
}
#[derive(Clone, Copy)]
struct Affine {
    columns: [Point; 4],
    determinant: I,
    sigma: [f64; 2],
}
impl Affine {
    fn new(columns: [Point; 4]) -> Result<Self, GeometryError> {
        if !columns.iter().flatten().all(|x| x.is_finite()) {
            return Err(unsupported("Actual affine columns are not finite"));
        }
        let axes = [columns[1], columns[2], columns[3]].map(|a| a.map(I::point));
        let mut gram = I::ZERO;
        for (i, a) in axes.iter().enumerate() {
            for (j, b) in axes.iter().enumerate() {
                let delta = dot(*a, *b) - I::point(if i == j { 1. } else { 0. });
                gram = gram + I::point(delta.magnitude());
            }
        }
        let determinant = dot(cross(axes[0], axes[1]), axes[2]);
        if !gram.valid() || gram.hi > 1e-6 || !determinant.valid() || determinant.lo <= 0. {
            return Err(unsupported(
                "Actual affine metric/orientation cannot be certified",
            ));
        }
        // Entrywise L1 bounds spectral Gram error, hence the complete
        // singular-value interval without an orthogonality assumption.
        let sigma = [
            (I::ONE - I::point(gram.hi)).sqrt().lo,
            (I::ONE + I::point(gram.hi)).sqrt().hi,
        ];
        Ok(Self {
            columns,
            determinant,
            sigma,
        })
    }
    fn linear(&self, p: BoxPoint) -> BoxPoint {
        std::array::from_fn(|i| {
            (0..3).fold(I::ZERO, |s, j| s + I::point(self.columns[j + 1][i]) * p[j])
        })
    }
    fn point(&self, p: BoxPoint) -> BoxPoint {
        let q = self.linear(p);
        std::array::from_fn(|i| I::point(self.columns[0][i]) + q[i])
    }
    fn literal(&self, p: Point) -> Point {
        std::array::from_fn(|i| {
            self.columns[0][i]
                + self.columns[1][i] * p[0]
                + self.columns[2][i] * p[1]
                + self.columns[3][i] * p[2]
        })
    }
}
#[derive(Clone, Copy)]
pub(super) struct Placement {
    construction: Affine,
    occurrence: Affine,
    top_z: f64,
    pub(super) determinant: I,
    pub(super) sigma_max: f64,
    pub(super) condition: f64,
}
impl Placement {
    pub(super) fn new(
        chart: multicad_rust_kernel::step_read::PhillipsCommonMouthChart,
        occurrence: PlaneFrame,
    ) -> Result<Self, GeometryError> {
        if !chart.outer_top_z.is_finite() {
            return Err(unsupported("Actual common station is not finite"));
        }
        let construction = Affine::new(chart.columns)?;
        let occurrence = Affine::new([
            occurrence.origin,
            occurrence.u,
            occurrence.v,
            occurrence.normal(),
        ])?;
        let high = I::point(construction.sigma[1]) * I::point(occurrence.sigma[1]);
        let low = I::point(construction.sigma[0]) * I::point(occurrence.sigma[0]);
        let condition = high / low;
        let determinant = construction.determinant * occurrence.determinant;
        // Each sealed Affine earned Gram error <=1e-6 and a finite positive
        // determinant. Its singular values lie near one, so these two-factor
        // products and their ratio are finite and low is strictly positive.
        Ok(Self {
            construction,
            occurrence,
            top_z: chart.outer_top_z,
            determinant,
            sigma_max: high.hi,
            condition: condition.hi,
        })
    }
    fn shifted(&self, p: BoxPoint) -> BoxPoint {
        [p[0], p[1], p[2] - I::point(self.top_z)]
    }
    fn exact(&self, p: BoxPoint) -> BoxPoint {
        self.occurrence
            .point(self.construction.point(self.shifted(p)))
    }
    pub(super) fn is_identity(&self) -> bool {
        let identity = [[0.; 3], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
        // This is the FULL material map, including the symbolic station
        // translation. The equal literal Z values cancel in exact arithmetic;
        // no rounded collapsed origin or near-identity tolerance is used.
        self.construction.columns[1..] == identity[1..]
            && self.construction.columns[0] == [0., 0., self.top_z]
            && self.occurrence.columns == identity
    }
    fn linear(&self, p: BoxPoint) -> BoxPoint {
        self.occurrence.linear(self.construction.linear(p))
    }
    pub(super) fn point(&self, p: Point) -> Result<(Point, f64), GeometryError> {
        let shifted = [p[0], p[1], p[2] - self.top_z];
        let literal = self.occurrence.literal(self.construction.literal(shifted));
        let exact = self.exact(p.map(I::point));
        let error = (0..3).fold(I::ZERO, |s, i| {
            s + I::point((I::point(literal[i]) - exact[i]).magnitude())
        });
        if !literal.into_iter().all(f64::is_finite) || !error.valid() {
            return Err(unsupported("Literal affine point transfer is unresolved"));
        }
        Ok((literal, error.hi))
    }
    /// Enclose literal transport, then require its measured outward bound to
    /// fit the supplied nonnegative distance reserve (zero may refuse).
    pub(super) fn point_within(
        &self,
        p: Point,
        reserve: f64,
    ) -> Result<(Point, f64), GeometryError> {
        let mapped = self.point(p)?;
        if !reserve.is_finite() || mapped.1 > reserve {
            return Err(unsupported(
                "Literal world points exceed the reserved source distance",
            ));
        }
        Ok(mapped)
    }
    pub(super) fn volume_interval(&self, local: [f64; 2]) -> Result<[f64; 2], GeometryError> {
        let volume = I::new(local[0], local[1]) * self.determinant;
        if !volume.valid() || volume.lo <= 0. {
            return Err(unsupported("Placed new-source volume cannot be enclosed"));
        }
        Ok([volume.lo, volume.hi])
    }
    /// Actual emitted normals must fit the same-parameter affine normal
    /// transport; an opposite or degenerate literal facet refuses.
    pub(super) fn require_normal(
        &self,
        local: [Point; 3],
        world: [Point; 3],
        reserve: f64,
    ) -> Result<(), GeometryError> {
        if !self.normal_transfer(local, world, reserve) {
            return Err(unsupported(
                "Literal world normals exceed the reserved angular transfer",
            ));
        }
        Ok(())
    }
    pub(super) fn normal_transfer(
        &self,
        local: [Point; 3],
        world: [Point; 3],
        reserve: f64,
    ) -> bool {
        // Translation is absent from exact edge transport. Compare the true
        // composed affine facet normal with the ACTUAL emitted world normal.
        let a = self.linear(difference(local[1], local[0]));
        let b = self.linear(difference(local[2], local[0]));
        let exact = cross(a, b);
        let actual = cross(
            difference(world[1], world[0]),
            difference(world[2], world[0]),
        );
        let projection = dot(exact, actual);
        let perpendicular = cross(exact, actual)
            .into_iter()
            .fold(I::ZERO, |s, v| s + I::point(v.magnitude()));
        let ratio = perpendicular / projection;
        // atan(ratio)<=ratio, avoiding any transcendental rounding assumption.
        projection.valid() && projection.lo > 0. && ratio.valid() && ratio.hi <= reserve
    }
    pub(super) fn constructed_bounds(&self, bounds: Bounds) -> Result<Bounds, GeometryError> {
        let local = std::array::from_fn(|i| I::new(bounds.min[i], bounds.max[i]));
        Self::enclosure_bounds(self.construction.point(self.shifted(local)))
    }
    fn enclosure_bounds(world: BoxPoint) -> Result<Bounds, GeometryError> {
        if !world.iter().all(|b| b.valid()) {
            return Err(unsupported("Actual affine source bounds are unresolved"));
        }
        Ok(Bounds {
            min: world.map(|b| b.lo),
            max: world.map(|b| b.hi),
        })
    }
    pub(super) fn bounds(&self, bounds: Bounds) -> Result<Bounds, GeometryError> {
        let local = std::array::from_fn(|i| I::new(bounds.min[i], bounds.max[i]));
        let world = self.exact(local);
        if !world.iter().all(|b| b.valid()) {
            return Err(unsupported("Actual affine source bounds are unresolved"));
        }
        Ok(Bounds {
            min: world.map(|b| b.lo),
            max: world.map(|b| b.hi),
        })
    }
}
