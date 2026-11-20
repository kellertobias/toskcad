//! Pairwise Booleans that stay in the warped-union chart.
//!
//! A quarter-round fillet sweep, with or without a base and back wall, is
//! exactly a planar union in the chart whose height coordinate straightens
//! its arc. Planar cutters whose walls through the rounded band are vertical
//! keep their planes in that chart, so a pairwise subtraction runs through
//! the same certified planar Boolean the compiler uses for a whole tool list.
//! A retained union joins another planar operand once its positive cells are
//! recovered from the certified cut, and planar stock or a retained union
//! takes a cylinder as an attached boss whose attachment is proved exactly.
//! Every refusal of these proofs leaves the pair to the caller's refusal.
use super::*;
use multicad_rust_kernel::planar_boolean::CutMaterial;

impl Chart {
    /// Replace stock-and-cutter pieces with the certified positive cells of
    /// their difference, so later operands join the cut material.
    fn into_material(mut self, tolerance: f64) -> Result<Self, GeometryError> {
        if self.pieces.len() == self.stock_count {
            return Ok(self);
        }
        let (stock, tools) = self.pieces.split_at(self.stock_count);
        let boundary = Union::subtract(stock, tools, tolerance)
            .map_err(source_error)?
            .ok_or_else(|| input("The retained cut has no material left to join"))?;
        let material =
            CutMaterial::recover(&boundary, stock, tools, tolerance).map_err(source_error)?;
        let operand = self.sources.first().map_or(0, |part| part.operand);
        self.pieces = material.pieces().to_vec();
        self.sources = vec![
            SourcePart {
                operand,
                faces: Vec::new(),
            };
            self.pieces.len()
        ];
        self.stock_count = self.pieces.len();
        Ok(self)
    }
}

fn removed_all(error: &GeometryError) -> bool {
    error.kind == GeometryErrorKind::InvalidResult && error.message.contains("removed all material")
}

impl Solid {
    pub(in crate::geometry) fn pairwise_warped_boolean(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if policy.boolean_fuzzy_mm != 0. || self.policy != policy || other.policy != policy {
            return None;
        }
        let result = match operation {
            BooleanOperation::Subtract => self.warped_subtract(other, policy)?,
            BooleanOperation::Union => self
                .warped_join(other, policy)
                .or_else(|| other.warped_join(self, policy))?,
            BooleanOperation::Intersect => return None,
        };
        match result {
            Ok(solid) => Some(Ok(Some(solid))),
            Err(error) if removed_all(&error) => Some(Ok(None)),
            Err(error) if error.kind == GeometryErrorKind::Limit => Some(Err(error)),
            Err(_) => None,
        }
    }

    /// One quarter-round fillet sweep as a retained chart union.
    fn quarter_round_chart(&self, policy: GeometryPolicy) -> Option<Result<Self, GeometryError>> {
        let tolerance = policy.linear_tolerance_mm;
        let fillet = QuarterFillet::recognize_open(self, tolerance)?;
        let feature = self.rust_feature()?;
        let mut chart = Chart::new(Some(Band {
            center_z_mm: fillet.center[2],
            radius_mm: fillet.radius,
        }));
        if let Err(cause) = chart.append(self, Some(fillet), policy)? {
            return Some(Err(cause));
        }
        chart.stock_count = chart.pieces.len();
        Some(
            Union::new(&chart.pieces, tolerance)
                .map_err(source_error)
                .and_then(|union| chart.finish(union, feature, policy)),
        )
    }

    fn warped_subtract(
        &self,
        tool: &Self,
        policy: GeometryPolicy,
    ) -> Option<Result<Self, GeometryError>> {
        if self.warped_union_body().is_some() {
            return self.subtract_warped_tools(&[tool], self.rust_feature()?, policy);
        }
        let stock = match self.quarter_round_chart(policy)? {
            Ok(stock) => stock,
            Err(cause) => return Some(Err(cause)),
        };
        stock.subtract_warped_tools(&[tool], self.rust_feature()?, policy)
    }

    /// Join `addition` to this retained union or planar stock: a cylinder
    /// becomes an attached boss, anything else must be planar chart material
    /// clear of existing bosses. Axial cuts stay last in the retained order,
    /// so a retained union with axial cuts takes nothing more.
    fn warped_join(
        &self,
        addition: &Self,
        policy: GeometryPolicy,
    ) -> Option<Result<Self, GeometryError>> {
        let tolerance = policy.linear_tolerance_mm;
        let boss = addition.axial_cylinder_frame().is_some();
        let (chart, feature, frame) = if let Some((body, frame)) = self.warped_union_body() {
            let chart = body.chart.clone()?;
            if !chart.axial.is_empty() || !boss && !chart.bosses.is_empty() {
                return None;
            }
            (chart, body.feature, frame)
        } else {
            // Plain stock only takes a boss here; planar pairs and fillets
            // have their own providers.
            if !boss || self.axial_cylinder_frame().is_some() || self.normal_sweep_body().is_some()
            {
                return None;
            }
            let mut chart = Chart::new(None);
            if let Err(cause) = chart.append(self, None, policy)? {
                return Some(Err(cause));
            }
            chart.stock_count = chart.pieces.len();
            (chart, addition.rust_feature()?, PlaneFrame::default())
        };
        let axes = [frame.u, frame.v, frame.normal()];
        let local = if frame == PlaneFrame::default() {
            addition.share()
        } else {
            addition.transformed(RigidTransform {
                rotation: axes,
                translation_mm: axes.map(|axis| -dot(axis, frame.origin)),
            })
        };
        let result = local.and_then(|local| {
            let mut chart = chart.into_material(tolerance)?;
            let appended = if boss {
                chart.append_boss(&local)
            } else {
                chart.append(&local, None, policy)
            };
            appended.ok_or_else(|| input("The joined operand is not chart material"))??;
            chart.stock_count = chart.pieces.len();
            let union = Union::new(&chart.pieces, tolerance).map_err(source_error)?;
            let joined = chart.finish(union, feature, policy)?;
            if frame == PlaneFrame::default() {
                Ok(joined)
            } else {
                joined.transformed(RigidTransform {
                    rotation: std::array::from_fn(|i| axes.map(|axis| axis[i])),
                    translation_mm: frame.origin,
                })
            }
        });
        Some(result)
    }
}

#[cfg(test)]
mod tests;
