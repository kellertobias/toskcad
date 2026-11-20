//! Intersection with a compound of separate solids, member by member.
//!
//! Intersection distributes over union, so a solid meets a compound exactly
//! where it meets each member. Every member pair goes through the ordinary
//! exact Boolean; a member refusal refuses the whole pair, empty members drop
//! out, and several common pieces form one compound again.
use super::*;

impl Solid {
    pub(super) fn compound_intersection(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if operation != BooleanOperation::Intersect
            || policy != self.policy
            || policy != other.policy
            || !self.uses_rust_geometry()
            || !other.uses_rust_geometry()
        {
            return None;
        }
        let (compound, solid, left) = match (self.rust_compound(), other.rust_compound()) {
            (Some(compound), _) => (compound, other, true),
            (None, Some(compound)) => (compound, self, false),
            (None, None) => return None,
        };
        let mut pieces = Vec::new();
        for member in compound.members() {
            let common = if left {
                member.boolean(solid, operation, policy)
            } else {
                solid.boolean(member, operation, policy)
            };
            match common {
                Ok(Some(piece)) => pieces.push(piece),
                Ok(None) => {}
                Err(error) => return Some(Err(error)),
            }
        }
        Some(match pieces.len() {
            0 => Ok(None),
            1 => Ok(pieces.pop()),
            _ => {
                let pieces: Vec<&Self> = pieces.iter().collect();
                Self::make_rust_compound(&pieces, compound.feature, policy, None).map(Some)
            }
        })
    }
}

#[cfg(test)]
mod tests;
