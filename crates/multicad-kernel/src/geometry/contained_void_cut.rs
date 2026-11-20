//! Subtract a tool lying strictly inside box stock as an enclosed void.
//!
//! A hole-free box is convex. When the tool's bounds sit inside it by more
//! than the construction tolerance on every side, every tool point is
//! interior stock material, so the difference is the stock boundary with the
//! tool boundary added, reversed, as one void shell. Both shells are the
//! operands' own exact analytic exchange boundaries, read back by the
//! retained-boundary reader that also validates the void's nesting. Nothing
//! is approximated and no cut surface is computed.
use super::*;
use crate::geometry::helical_body::contained_boolean::chart_bounds;
use crate::geometry::internal_thread_treatment_body::interval::Interval as B;
use multicad_rust_kernel::brep;

const OPERATION: &str = "contained void cut";

/// One operand's retained boundary with exactly one closed shell.
fn single_shell(solid: &Solid, policy: GeometryPolicy) -> Option<brep::Boundary> {
    let text = solid.rust_step_text(OPERATION).ok()?;
    let imported = Solid::import_rust_brep(&text, EntityId::new(), policy).ok()?;
    let boundary = imported.imported_brep_body()?.solid.boundary().clone();
    (boundary.shells.len() == 1).then_some(boundary)
}

/// Append `void` to `stock` as one reversed interior shell.
fn with_void(mut stock: brep::Boundary, void: brep::Boundary) -> brep::Boundary {
    let [vertices, edges, faces] = [stock.vertices.len(), stock.edges.len(), stock.faces.len()];
    stock.vertices.extend(void.vertices);
    stock.edges.extend(void.edges.into_iter().map(|mut edge| {
        edge.vertices = edge.vertices.map(|v| v + vertices);
        edge
    }));
    stock.faces.extend(void.faces.into_iter().map(|mut face| {
        for ring in &mut face.loops {
            for (edge, _) in &mut ring.uses {
                *edge += edges;
            }
            if let Some(vertex) = &mut ring.vertex {
                *vertex += vertices;
            }
        }
        face
    }));
    stock
        .shells
        .extend(void.shells.into_iter().map(|shell| brep::Shell {
            faces: shell.faces.into_iter().map(|f| f + faces).collect(),
            reversed: !shell.reversed,
        }));
    stock
}

impl Solid {
    pub(super) fn contained_void_cut(
        &self,
        tool: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if operation != BooleanOperation::Subtract
            || policy.boolean_fuzzy_mm != 0.
            || policy != self.policy
            || policy != tool.policy
            || !tool.uses_rust_geometry()
        {
            return None;
        }
        let (prism, frame) = self.rust_prism_frame()?;
        if prism.rounded_end().is_some() || prism.inner_loops().next().is_some() {
            return None;
        }
        let size = prism.dimensions().map(Rational::as_f64);
        let bounds = tool.bounds().ok()?;
        let world: [B; 3] = std::array::from_fn(|i| B::new(bounds.min[i], bounds.max[i]));
        let local = chart_bounds(world, frame.origin, [frame.u, frame.v, frame.normal()]);
        let magnitude = size
            .into_iter()
            .chain(frame.origin)
            .chain(bounds.min)
            .chain(bounds.max)
            .map(f64::abs)
            .fold(1., f64::max);
        let margin = 2. * policy.linear_tolerance_mm + 64. * f64::EPSILON * magnitude;
        let inside = (0..3).all(|i| {
            let [lo, hi] = local[i];
            lo.is_finite() && hi.is_finite() && lo > margin && hi < size[i] - margin
        });
        if !inside {
            return None;
        }
        let feature = self.rust_feature()?;
        let stock = single_shell(self, policy)?;
        let void = single_shell(tool, policy)?;
        let result = brep::Solid::new(with_void(stock, void), policy.linear_tolerance_mm, None)
            .map_err(|cause| {
                error(
                    OPERATION,
                    GeometryErrorKind::InvalidResult,
                    format!("The enclosed void boundary did not validate: {cause:?}"),
                )
            })
            .and_then(|solid| Self::from_assembled_brep(solid, feature, policy));
        Some(result.map(Some))
    }
}

#[cfg(test)]
mod tests;
