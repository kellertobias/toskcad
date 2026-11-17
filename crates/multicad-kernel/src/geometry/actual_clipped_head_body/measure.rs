//! Selected actual reference properties. Normal plans are transient; only
//! completed scalar enclosures are cached, without retaining their big nets.
use super::*;
impl ActualClippedHeadBody {
    fn local_face(&self, index: usize) -> Result<LocalFace, GeometryError> {
        if index >= 17 {
            return Err(failure(read::Error::Invalid));
        }
        if let Some(value) = self.local_faces.borrow()[index] {
            return Ok(value);
        }
        let scale = self.source.actual_frame().scale_bounds();
        let sigma = self
            .occurrence
            .as_ref()
            .map_or(1., |a| a.linear_l2_norm_upper());
        let stretch = I::point(scale[1]) * I::point(sigma);
        let area_width = (I::point(1e-6) / stretch.square()).lo;
        let centroid_width = (I::point(self.policy.linear_tolerance_mm) / stretch).lo;
        if !area_width.is_finite()
            || area_width <= 0.
            || !centroid_width.is_finite()
            || centroid_width <= 0.
        {
            return Err(failure(read::Error::Unsupported));
        }
        let options = HeadMeasureOptions {
            measure: read::FaceMeasureOptions {
                area_width_mm2: area_width,
                centroid_width_mm: centroid_width,
                ..Default::default()
            },
            additional_resident_bytes: self.external_bytes()
                + self
                    .occurrence
                    .as_ref()
                    .map_or(0, |a| a.retained_budget_bytes()),
        };
        let value = if index < 8 {
            let plan = HeadGeometryPlan::new_mapped(
                Arc::clone(&self.source),
                index,
                HeadGeometryOptions {
                    additional_resident_bytes: options.additional_resident_bytes,
                    ..Default::default()
                },
                None,
            )
            .map_err(failure)?;
            let value = plan.face_measure(index, options, None).map_err(failure)?;
            LocalFace {
                area: value.area_mm2(),
                center: value.centroid_mm(),
            }
        } else {
            // This helper owns only the local physical arena; declare the
            // mapped file/basis still resident as well, once.
            let local = self.source.local_geometry();
            let extra = self
                .source
                .receipt()
                .retained_bytes
                .checked_sub(local.receipt().retained_bytes)
                .and_then(|n| n.checked_add(options.additional_resident_bytes))
                .ok_or_else(|| failure(read::Error::Limit))?;
            let value = local
                .planar_face_measure(
                    index,
                    HeadMeasureOptions {
                        additional_resident_bytes: extra,
                        ..options
                    },
                    None,
                )
                .map_err(failure)?;
            LocalFace {
                area: value.area_mm2(),
                center: value.centroid_mm(),
            }
        };
        // Validate current world scalar widths before publishing a new local
        // cache entry. Failed requests preserve prior successful descriptors.
        self.world_face(value)?;
        self.local_faces.borrow_mut()[index] = Some(value);
        self.face_builds.set(self.face_builds.get() + 1);
        Ok(value)
    }
    fn world_face(&self, local: LocalFace) -> Result<LocalFace, GeometryError> {
        let frame = self.source.actual_frame();
        let scale = frame.scale_bounds();
        let mut area = I::new(local.area[0], local.area[1]) * I::new(scale[0], scale[1]).square();
        let mut center = frame.map_box_mm(local.center).map_err(failure)?;
        if let Some(a) = &self.occurrence {
            let b = self.source.world_sewn_enclosure_mm();
            let mut radius = I::ZERO;
            for (b, c) in b.into_iter().zip(ranges(center)) {
                radius = radius + I::point((I::new(b[0], b[1]) - c).magnitude());
            }
            let delta = I::point(a.metric_error());
            let shift = I::point(a.linear_l1_norm()) * radius * delta / (I::ONE - delta);
            center = a.map_box_mm(center).map_err(failure)?;
            if !shift.valid() {
                return Err(failure(read::Error::Unsupported));
            }
            center = output(ranges(center).map(|v| v + I::new(-shift.hi, shift.hi)));
            let factor = a.area_scale_bounds();
            area = area * I::new(factor[0], factor[1]);
        }
        if !area.valid()
            || area.lo <= 0.
            || area.width() > 1e-6
            || ranges(center)
                .iter()
                .any(|c| !c.valid() || c.width() > self.policy.linear_tolerance_mm)
        {
            return Err(failure(read::Error::Unsupported));
        }
        Ok(LocalFace {
            area: [area.lo, area.hi],
            center,
        })
    }
    pub(in crate::geometry) fn face(&self, index: u32) -> Result<Face, GeometryError> {
        let world = self.world_face(self.local_face(index as usize)?)?;
        let center_mm = world.center.map(|v| v[0].midpoint(v[1]));
        self.policy.point(center_mm)?;
        Ok(Face {
            index,
            plane: self.plane(index)?.map(|f| (f.origin_mm, f.normal())),
            center_mm,
            area_mm2: world.area[0].midpoint(world.area[1]),
        })
    }
    pub(in crate::geometry) fn faces(&self) -> Result<Vec<Face>, GeometryError> {
        (0..17).map(|i| self.face(i)).collect()
    }
}
