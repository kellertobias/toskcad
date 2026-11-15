//! Whole actual screw versus retained stock, tap and crown relief.
use super::*;
use multicad_rust_kernel::wood_thread_contact as source;
fn frame(matrix: [[f64; 3]; 3], origin_mm: Point) -> source::Frame {
    source::Frame {
        origin_mm,
        axes: std::array::from_fn(|j| std::array::from_fn(|i| matrix[i][j])),
    }
}
impl Solid {
    /// Proves the actual profiled screw against the original stock and both
    /// finite receiving domains. The stock translation, bore offset and outer
    /// placement are retained as separate exact symbolic expressions.
    pub fn rust_line8_nut_mating<'a>(
        &'a self,
        nut: &'a Self,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<source::ClockedLine8Mating<'a>, GeometryError> {
        self.policy.validate()?;
        if self.policy != nut.policy || self.policy.boolean_fuzzy_mm != 0. {
            return Err(error(
                "Line8 mating",
                GeometryErrorKind::Input,
                "Mating requires equal source policies without Boolean fuzz",
            ));
        }
        let unsupported = || {
            error(
                "Line8 mating",
                GeometryErrorKind::Unsupported,
                "Mating requires the actual profiled hex-socket screw and retained Line8 stock",
            )
        };
        let (screw, placed) = self
            .rust_profiled_screw_geometry()
            .ok_or_else(unsupported)?;
        let SolidBody::Line8Nut(body) = &nut.body else {
            return Err(unsupported());
        };
        let screw_frame = source::Frame {
            origin_mm: placed.origin,
            axes: [placed.u, placed.v, placed.normal()],
        };
        source::certify_clocked_line8_mating(
            source::ClockedLine8Placement {
                source: &body.source,
                stock_frame: frame(body.stock.matrix(), body.stock.translation_mm()),
                outer_frame: frame(body.outer.rotation, body.outer.translation_mm),
                minimum_wall_mm: body.requested_wall_mm,
            },
            source::MatingExternal {
                source: source::MatingScrew::Profiled(screw),
                frame: screw_frame,
            },
            source::ClockedMatingOptions {
                requirements: source::MatingRequirements {
                    minimum_engagement_mm: 2. * screw.thread().spec().pitch_mm,
                    minimum_clearance_mm: 4. * self.policy.linear_tolerance_mm,
                },
                maximum_cells: 4096,
                linear_tolerance_mm: self.policy.linear_tolerance_mm,
            },
            &|| interrupt.is_some_and(|f| f()),
        )
        .map_err(|cause| {
            error(
                "Line8 mating",
                match cause {
                    source::Error::Invalid => GeometryErrorKind::Input,
                    source::Error::Limit => GeometryErrorKind::Limit,
                    source::Error::Cancelled => GeometryErrorKind::Cancelled,
                    _ => GeometryErrorKind::Unsupported,
                },
                format!("Original stock/tap/relief mating proof failed: {cause:?}"),
            )
        })
    }
}
