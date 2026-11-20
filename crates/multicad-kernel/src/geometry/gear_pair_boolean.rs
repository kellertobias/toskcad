//! Booleans of two placed Rust gear solids settled by their pair certificate.
//!
//! A `Clear` certificate proves the actual solids share no point, so an
//! intersection is empty and a difference keeps the unchanged left operand.
//! A proved common interior point means real interference; its exact volume
//! needs a curved surface-surface intersection, which no Rust provider builds
//! yet, so that case is a typed refusal rather than a sampled solid.
use super::*;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Pair {
    Clear,
    Interferes,
    Unresolved,
}

impl From<multicad_rust_kernel::spur_gear::clearance::ResultKind> for Pair {
    fn from(kind: multicad_rust_kernel::spur_gear::clearance::ResultKind) -> Self {
        use multicad_rust_kernel::spur_gear::clearance::ResultKind::*;
        match kind {
            Clear => Self::Clear,
            CommonInterior => Self::Interferes,
            Unresolved => Self::Unresolved,
        }
    }
}

impl From<multicad_rust_kernel::bevel_gear::clearance::ResultKind> for Pair {
    fn from(kind: multicad_rust_kernel::bevel_gear::clearance::ResultKind) -> Self {
        use multicad_rust_kernel::bevel_gear::clearance::ResultKind::*;
        match kind {
            Clear => Self::Clear,
            CommonInterior => Self::Interferes,
            Unresolved => Self::Unresolved,
        }
    }
}

fn is_gear(solid: &Solid) -> bool {
    matches!(
        solid.body,
        SolidBody::SpurGear(..)
            | SolidBody::InternalRing(..)
            | SolidBody::Herringbone(..)
            | SolidBody::Bevel(..)
            | SolidBody::Worm(..)
    )
}

impl Solid {
    fn certify_gear_pair(
        &self,
        other: &Self,
        policy: GeometryPolicy,
    ) -> Option<Result<Pair, GeometryError>> {
        use SolidBody::*;
        Some(match (&self.body, &other.body) {
            (SpurGear(..), SpurGear(..)) => self
                .certify_rust_spur_pair(other, policy, None)
                .map(Pair::from),
            (Herringbone(..), Herringbone(..)) => self
                .certify_rust_herringbone_pair(other, policy, None)
                .map(Pair::from),
            (InternalRing(..), SpurGear(..)) => self
                .certify_rust_internal_pair(other, policy, None)
                .map(Pair::from),
            (SpurGear(..), InternalRing(..)) => other
                .certify_rust_internal_pair(self, policy, None)
                .map(Pair::from),
            (Bevel(..), Bevel(..)) => self
                .certify_rust_bevel_pair(other, policy, 0., None)
                .map(|(kind, _)| Pair::from(kind)),
            (Worm(..), Worm(..)) => self
                .certify_rust_worm_pair(other, policy, 0., None)
                .map(|(kind, _)| Pair::from(kind)),
            // Other sweeps get an Unsupported refusal, left to later providers.
            (SpurGear(..), NormalSweep(..)) | (NormalSweep(..), SpurGear(..)) => self
                .certify_rust_rack_pair(other, policy, None)
                .map(Pair::from),
            _ => return None,
        })
    }

    /// `None` leaves non-gear operands, unions, fuzzy policies and
    /// unresolved or unsupported poses to the remaining providers.
    pub(in crate::geometry) fn gear_pair_boolean(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if operation == BooleanOperation::Union
            || policy.boolean_fuzzy_mm != 0.
            || !(is_gear(self) || is_gear(other))
        {
            return None;
        }
        match self.certify_gear_pair(other, policy)? {
            Ok(Pair::Clear) => Some(match operation {
                BooleanOperation::Subtract => self.share().map(Some),
                _ => Ok(None),
            }),
            Ok(Pair::Interferes) => Some(Err(error(
                "Boolean",
                GeometryErrorKind::Unsupported,
                format!(
                    "The placed {} and {} gears interfere; Rust has no exact provider for their curved interference volume",
                    self.body_kind_name(),
                    other.body_kind_name()
                ),
            ))),
            Ok(Pair::Unresolved) => None,
            Err(error) if error.kind == GeometryErrorKind::Unsupported => None,
            Err(error) => Some(Err(error)),
        }
    }
}
