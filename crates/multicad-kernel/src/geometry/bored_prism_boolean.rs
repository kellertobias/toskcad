//! Subtraction from a box prism that already carries through bores.
//!
//! The through bores of a rectangular Rust prism are full-height cylinders.
//! When a tool's bounds are separated from every bore cylinder by more than
//! the linear tolerance, removing the tool and removing the bores commute:
//! `(box - bores) - tool = (box - tool) - bores`. The plain box cut is left
//! to the exact rectilinear providers and the same bores, with their source
//! identities, are drilled through the cut stock again.
use super::*;

/// The world bounds of the through cylinder `circle` in `prism`'s chart.
fn bore_bounds(frame: PlaneFrame, circle: Circle, height: f64) -> Bounds {
    let [x, y] = circle.center().map(Rational::as_f64);
    let r = circle.radius().as_f64();
    let normal = frame.normal();
    let corners = [-r, r].into_iter().flat_map(|du| {
        [-r, r].into_iter().flat_map(move |dv| {
            [0., height].map(move |z| -> [f64; 3] {
                std::array::from_fn(|i| {
                    frame.origin[i] + (x + du) * frame.u[i] + (y + dv) * frame.v[i] + z * normal[i]
                })
            })
        })
    });
    let mut bounds = Bounds {
        min: [f64::INFINITY; 3],
        max: [f64::NEG_INFINITY; 3],
    };
    for corner in corners {
        for i in 0..3 {
            bounds.min[i] = bounds.min[i].min(corner[i]);
            bounds.max[i] = bounds.max[i].max(corner[i]);
        }
    }
    bounds
}

fn separated(a: Bounds, b: Bounds, tolerance: f64) -> bool {
    (0..3).any(|i| a.max[i] + tolerance < b.min[i] || b.max[i] + tolerance < a.min[i])
}

impl Solid {
    pub(super) fn bored_prism_subtract(
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
        let (prism, frame) = self.rust_prism_frame()?;
        if prism.rounded_end().is_some() || prism.inner_loops().len() == 0 {
            return None;
        }
        let height = prism.dimensions()[2];
        let reach = tool.bounds().ok()?;
        let tolerance = policy.linear_tolerance_mm;
        let bores: Vec<_> = prism.inner_loops().collect();
        if !bores.iter().all(|&(_, circle)| {
            separated(
                bore_bounds(frame, circle, height.as_f64()),
                reach,
                tolerance,
            )
        }) {
            return None;
        }
        let plain = Self {
            body: SolidBody::Prism(
                Rc::new(Prism::rectangular(prism.feature(), prism.dimensions()).ok()?),
                frame,
            ),
            policy,
            provenance: vec![],
            rust_origins: vec![],
            retained_placement: None,
        };
        let cut = match plain.boolean_geometry(tool, operation, policy) {
            Ok(Some(cut)) => cut,
            Ok(None) => return Some(Ok(None)),
            Err(error) if error.kind == GeometryErrorKind::Unsupported => return None,
            Err(error) => return Some(Err(error)),
        };
        // Drill each bore from one unit below the bottom cap to one unit
        // above the top cap, along the prism normal.
        let normal = frame.normal();
        let length = height.as_f64() + 2.;
        let start = PlaneFrame {
            origin: std::array::from_fn(|i| frame.origin[i] - normal[i]),
            ..frame
        };
        let drills: Vec<_> = bores
            .iter()
            .map(|&(feature, circle)| CircularPocketTool {
                frame: start,
                center_uv: circle.center().map(Rational::as_f64),
                radius_mm: circle.radius().as_f64(),
                displacement_mm: normal.map(|n| n * length),
                origin: Origin {
                    source: EntityId::from_uuid_bits(feature),
                    role: "side.cylinder".into(),
                },
            })
            .collect();
        Some(
            cut.rust_drill_circular_machining(&drills, &vec![true; drills.len()], policy)
                .map(Some),
        )
    }
}

#[cfg(test)]
mod tests;
