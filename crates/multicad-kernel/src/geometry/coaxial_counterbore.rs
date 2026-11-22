//! Boolean cuts that seat or add through bores in polygonal bore stock.
//!
//! A blind cylinder coaxial with an existing circular through bore, wider
//! than it and entering from outside one cap, is a counterbore: the exact
//! seated-bore provider certifies its shoulder, wall and clearance. A
//! through cylinder clear of everything in stock that already has
//! counterbores is a new hole; the stock's counterbores are reapplied to the
//! untreated section with the hole added, so `(stock - seats) - hole` equals
//! `(stock - hole) - seats` exactly. Only counterbore treatments are
//! rebuilt; chamfers and countersinks are left to other providers.
use super::*;
use multicad_rust_kernel::chamfered_bores::Counterbore;

/// The untreated polygon stock, its counterbore cuts with their origins.
type Seated = (Solid, Vec<(Counterbore, Origin)>);

impl Solid {
    fn seated_stock(&self) -> Option<Seated> {
        if let Some((prism, frame)) = self.rust_prism_frame() {
            if prism.rounded_end().is_some() {
                return None;
            }
            let feature = EntityId::from_uuid_bits(prism.feature());
            let [l, w, h] = prism.dimensions();
            let zero = Rational::new(0, 1).expect("zero is exact");
            let ring = vec![[zero, zero], [l, zero], [l, w], [zero, w]];
            let (circles, holes): (Vec<_>, Vec<_>) = prism
                .inner_loops()
                .map(|(source, circle)| {
                    (
                        circle,
                        Origin {
                            source: EntityId::from_uuid_bits(source),
                            role: "side.cylinder".into(),
                        },
                    )
                })
                .unzip();
            let origins: Vec<_> = (0..4)
                .map(|i| Origin {
                    source: feature,
                    role: format!("profile.edge.{i}"),
                })
                .chain(holes)
                .collect();
            let polygon = PolygonPrism::with_circular_holes(vec![ring], circles, h).ok()?;
            let stock = Self::from_polygon_prism_sweep(
                polygon,
                frame,
                false,
                feature,
                &origins,
                &vec![false; origins.len()],
                self.policy,
            )
            .ok()?;
            return Some((stock, vec![]));
        }
        let placed = self.polygon_prism()?;
        if !placed.pockets.is_empty() {
            return None;
        }
        let Some(treatments) = &placed.bore_treatments else {
            let mut stock = self.share().ok()?;
            stock.provenance.clear();
            return Some((stock, vec![]));
        };
        let geometry = &treatments.geometry;
        if !geometry.chamfers().is_empty()
            || !geometry.countersinks().is_empty()
            || treatments.face_origins.is_some()
        {
            return None;
        }
        let cuts = geometry
            .counterbores()
            .iter()
            .map(|&cut| Some((cut, treatments.origins.get(&(cut.hole, cut.end))?.clone())))
            .collect::<Option<Vec<_>>>()?;
        let mut untreated = placed.clone();
        untreated.bore_treatments = None;
        Some((
            Self {
                body: SolidBody::PolygonPrism(Rc::new(untreated)),
                policy: self.policy,
                provenance: vec![],
                rust_origins: self.rust_origins.clone(),
                retained_placement: None,
            },
            cuts,
        ))
    }

    pub(super) fn coaxial_counterbore_cut(
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
        let (cylinder, placement) = tool.axial_cylinder_frame()?;
        let (stock, mut cuts) = self.seated_stock()?;
        let placed = stock.polygon_prism()?;
        let frame = placed.frame;
        let normal = frame.normal();
        let axis = placement.normal();
        let side = cross(axis, normal);
        if dot(side, side).sqrt() > 64. * f64::EPSILON {
            return None;
        }
        let tolerance = policy.linear_tolerance_mm;
        let height = placed.prism.height().as_f64();
        let sign = if placed.reverse_sweep { -1. } else { 1. };
        let local = |p: [f64; 3]| {
            let d: [f64; 3] = std::array::from_fn(|i| p[i] - frame.origin[i]);
            [dot(d, frame.u), dot(d, frame.v), sign * dot(d, normal)]
        };
        let start = local(placement.origin);
        let finish = start[2] + sign * dot(axis, normal) * cylinder.height().as_f64();
        let center = [start[0], start[1]];
        let radius = cylinder.circle().radius().as_f64();
        let origin = tool.cylinder_sketch().map_or_else(
            || Origin {
                source: EntityId::from_uuid_bits(cylinder.feature()),
                role: "side.cylinder".into(),
            },
            |sketch| sketch.origin.clone(),
        );
        let [low, high] = [start[2].min(finish), start[2].max(finish)];
        let coaxial = placed
            .prism
            .circular_holes()
            .iter()
            .position(|hole| hole.center().map(Rational::as_f64) == center);
        let result = if let Some(hole) = coaxial {
            let bore = placed.prism.circular_holes()[hole].radius().as_f64();
            // Enter from outside one cap and stop inside the stock.
            let (end, depth) = if low < -tolerance && high > tolerance && high < height - tolerance
            {
                (false, high)
            } else if high > height + tolerance && low < height - tolerance && low > tolerance {
                (true, height - low)
            } else {
                return None;
            };
            if radius <= bore + tolerance
                || cuts.iter().any(|(c, _)| (c.hole, c.end) == (hole, end))
            {
                return None;
            }
            cuts.push((
                Counterbore {
                    hole,
                    end,
                    radius_mm: radius,
                    depth_mm: depth,
                },
                origin,
            ));
            stock
        } else {
            // A new through hole clear of every boundary, below counterbores.
            if cuts.is_empty() || low >= -tolerance || high <= height + tolerance {
                return None;
            }
            let [x, y] =
                center.map(|value| rational_coordinate(value, 0.).filter(|r| r.as_f64() == value));
            let r = rational_coordinate(radius, 0.).filter(|r| r.as_f64() == radius)?;
            let circle = Circle::new([x?, y?], r).ok()?;
            match Self::with_polygon_through_holes(placed, [(circle, origin)], policy) {
                Ok(holed) => holed,
                Err(_) => return None,
            }
        };
        let (cuts, origins): (Vec<_>, Vec<_>) = cuts.into_iter().unzip();
        match result.rust_counterbore_polygon(&cuts, &origins, policy) {
            Ok(seated) => Some(Ok(Some(seated))),
            Err(error) if error.kind == GeometryErrorKind::Limit => Some(Err(error)),
            Err(_) => None,
        }
    }
}

#[cfg(test)]
mod tests;
