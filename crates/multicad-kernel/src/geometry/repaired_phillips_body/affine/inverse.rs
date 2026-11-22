//! Actual interval cofactor inverse; forward transport retains its own proof.
use super::*;
impl Affine {
    fn inverse(&self, world: BoxPoint) -> BoxPoint {
        let axes = [self.columns[1], self.columns[2], self.columns[3]].map(|p| p.map(I::point));
        let rows = [
            cross(axes[1], axes[2]),
            cross(axes[2], axes[0]),
            cross(axes[0], axes[1]),
        ];
        let delta = std::array::from_fn(|i| world[i] - I::point(self.columns[0][i]));
        rows.map(|row| dot(row, delta) / self.determinant)
    }
}
impl Placement {
    /// Invert occurrence, then the actual four construction columns. Add the
    /// common station only AFTER the construction inverse, not in world space.
    pub(in crate::geometry::repaired_phillips_body) fn inverse_bounds(
        &self,
        world: [[f64; 2]; 3],
    ) -> Result<[[f64; 2]; 3], GeometryError> {
        let local = self
            .construction
            .inverse(self.occurrence.inverse(world.map(|b| I::new(b[0], b[1]))));
        let local = [local[0], local[1], local[2] + I::point(self.top_z)];
        if local.iter().any(|b| !b.valid()) {
            return Err(unsupported("Actual inverse affine enclosure is unresolved"));
        }
        Ok(local.map(|b| [b.lo, b.hi]))
    }
}
