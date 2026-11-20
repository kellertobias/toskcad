//! Boolean subtraction of cap-normal machining tools.
//!
//! A Boolean whose tool is an exact circular cylinder or a stadium sweep and
//! whose target is polygonal or rectilinear stock is the same exact cut as
//! the dedicated circular-pocket or rounded-slot provider. This adapter only
//! recognizes the operands; the providers certify entry, floors, clearance
//! and the result. A recognized pair that the providers refuse is left to the
//! later Boolean providers, so the final refusal stays the typed one.
use super::*;

impl Solid {
    /// Exact drilling stock for `self`: the solid itself when a machining
    /// provider accepts it, or the polygon reading of a plain normal sweep.
    /// An earlier Boolean's authored metadata describes the uncut body; the
    /// Boolean derives the result's metadata from the operands again.
    fn machining_stock(&self) -> Option<Self> {
        if self.orthogonal_body().is_some()
            || self.rust_prism().is_some()
            || self.linear_body().is_some()
            || self.scaled_orthogonal_cell().is_some()
            || self.polygon_prism().is_some()
        {
            let mut stock = self.share().ok()?;
            stock.provenance.clear();
            return Some(stock);
        }
        self.normal_sweep_polygon_stock()
    }

    fn cylinder_pocket_tool(&self) -> Option<CircularPocketTool> {
        let (cylinder, frame) = self.axial_cylinder_frame()?;
        let origin = self.cylinder_sketch().map_or_else(
            || Origin {
                source: EntityId::from_uuid_bits(cylinder.feature()),
                role: "side.cylinder".into(),
            },
            |sketch| sketch.origin.clone(),
        );
        Some(CircularPocketTool {
            frame,
            center_uv: [0.; 2],
            radius_mm: cylinder.circle().radius().as_f64(),
            displacement_mm: frame.normal().map(|n| n * cylinder.height().as_f64()),
            origin,
        })
    }

    /// `Some` only when a machining provider produced the exact cut, or
    /// when it reached a budget limit that no later provider can lift.
    pub(super) fn subtract_machining_tool(
        &self,
        tool: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if operation != BooleanOperation::Subtract
            || policy.boolean_fuzzy_mm != 0.
            || policy != self.policy
            || policy != tool.policy
        {
            return None;
        }
        let stock = self.machining_stock()?;
        let accept = |result: Result<Self, GeometryError>| match result {
            Ok(solid) => Some(Ok(Some(solid))),
            Err(error) if error.kind == GeometryErrorKind::Limit => Some(Err(error)),
            Err(_) => None,
        };
        if let Some(pocket) = tool.cylinder_pocket_tool() {
            let tools = std::slice::from_ref(&pocket);
            return accept(stock.rust_drill_circular_pockets(tools, policy))
                .or_else(|| accept(stock.rust_drill_circular_machining(tools, &[true], policy)));
        }
        let slot = tool.normal_sweep_slot_tool()?;
        accept(stock.rust_drill_rounded_slots(std::slice::from_ref(&slot), policy))
    }
}

#[cfg(test)]
mod tests;
