//! Pairwise union of one cylinder with a block stacked along its axis.
//!
//! The axial profile stack certifies the union of a circular cylinder and a
//! world-placed rectangular block whose faces are parallel or normal to the
//! cylinder axis, with the circle inscribed against at least two sides of
//! the rectangle where they overlap.
//! The pairwise Boolean names the result after the cylinder's source so a
//! re-evaluation produces the same identities.
use super::*;

impl Solid {
    pub(super) fn stacked_cylinder_union(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if operation != BooleanOperation::Union || policy != self.policy || policy != other.policy {
            return None;
        }
        let (cylinder, _) = self
            .axial_cylinder_frame()
            .or_else(|| other.axial_cylinder_frame())?;
        let feature = EntityId::from_uuid_bits(cylinder.feature());
        Self::union_round_rectangular_stack(&[self, other], feature, policy)
            .map(|result| result.map(Some))
    }
}

#[cfg(test)]
mod tests;
