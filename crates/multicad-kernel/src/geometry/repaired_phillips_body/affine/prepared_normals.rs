//! One actual composed cofactor map, shared by all literal facet checks.
use super::*;

pub(in crate::geometry::repaired_phillips_body) struct PreparedNormals {
    columns: [BoxPoint; 3],
}
impl Placement {
    pub(in crate::geometry::repaired_phillips_body) fn prepare_normals(&self) -> PreparedNormals {
        PreparedNormals::new(self.construction, self.occurrence)
    }
}
impl PreparedNormals {
    fn new(construction: Affine, occurrence: Affine) -> Self {
        let columns = [
            construction.columns[1],
            construction.columns[2],
            construction.columns[3],
        ]
        .map(|column| occurrence.linear(column.map(I::point)));
        // For every linear map A, cross(Au,Av)=cofactor(A)cross(u,v).
        // Directed arithmetic encloses the complete actual O*C map. No
        // orthogonality, normalized columns or nominal pose is assumed.
        Self {
            columns: [
                cross(columns[1], columns[2]),
                cross(columns[2], columns[0]),
                cross(columns[0], columns[1]),
            ],
        }
    }
    pub(in crate::geometry::repaired_phillips_body) fn require_normal(
        &self,
        local: [Point; 3],
        world: [Point; 3],
        reserve: f64,
    ) -> Result<(), GeometryError> {
        let normal = cross(
            difference(local[1], local[0]),
            difference(local[2], local[0]),
        );
        let exact = std::array::from_fn(|i| {
            (0..3).fold(I::ZERO, |sum, j| sum + self.columns[j][i] * normal[j])
        });
        let actual = cross(
            difference(world[1], world[0]),
            difference(world[2], world[0]),
        );
        let projection = dot(exact, actual);
        let perpendicular = cross(exact, actual)
            .into_iter()
            .fold(I::ZERO, |sum, v| sum + I::point(v.magnitude()));
        let ratio = perpendicular / projection;
        if projection.valid() && projection.lo > 0. && ratio.valid() && ratio.hi <= reserve {
            Ok(())
        } else {
            Err(unsupported(
                "Literal world normals exceed the reserved angular transfer",
            ))
        }
    }
}

#[cfg(test)]
mod tests;
