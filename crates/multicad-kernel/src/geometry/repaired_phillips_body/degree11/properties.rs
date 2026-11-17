//! Selected actual-source measures, independent of unrelated curved faces.
use super::*;
use multicad_rust_kernel::step_read::PhillipsDegree11FacePropertyOptions;
impl RepairedPhillipsBody {
    pub(crate) fn degree11_face_plane(
        &self,
        frame: PlaneFrame,
        index: u32,
    ) -> Result<Option<(Point, Point)>, GeometryError> {
        let source = self
            .source
            .degree11()
            .ok_or_else(|| unsupported("No degree-eleven material"))?;
        let plane =
            super::planes::plane(source, index as usize, self.placement(frame)?, self.policy)?;
        Ok(plane.map(|p| {
            let normal = std::array::from_fn(|k| {
                p.u[(k + 1) % 3] * p.v[(k + 2) % 3] - p.u[(k + 2) % 3] * p.v[(k + 1) % 3]
            });
            (p.origin_mm, normal)
        }))
    }
    pub(crate) fn degree11_face(
        &self,
        frame: PlaneFrame,
        index: u32,
    ) -> Result<Face, GeometryError> {
        let source = self
            .source
            .degree11()
            .ok_or_else(|| unsupported("No degree-eleven material"))?;
        let placement = self.placement(frame)?;
        let plane = super::planes::plane(source, index as usize, placement, self.policy)?;
        let extent = source
            .local_enclosure_bounds_mm()
            .into_iter()
            .flatten()
            .map(f64::abs)
            .fold(1_f64, f64::max);
        // Actual Native placement/plane/returned-value stack reservation, separate
        // from the borrowed immutable source and its delegated integration.
        let outside = self
            .resident_bytes()
            .checked_sub(source.retained_budget_bytes())
            .and_then(|n| n.checked_add(16 * 1024))
            .ok_or_else(query::exhausted)?;
        if source
            .retained_budget_bytes()
            .checked_add(outside)
            .is_none_or(|n| n > BYTES)
        {
            return Err(query::exhausted());
        }
        let cached = self
            .degree11_properties
            .borrow()
            .get(index as usize)
            .copied()
            .flatten();
        let mut pending = None;
        let mut outside = outside;
        let properties = if let Some(properties) = cached {
            properties
        } else {
            if self.degree11_properties.borrow().is_empty() {
                let count = source.face_count();
                let bytes = count
                    .checked_mul(std::mem::size_of::<Option<Degree11LocalFace>>())
                    .ok_or_else(query::exhausted)?;
                if self
                    .resident_bytes()
                    .checked_add(bytes)
                    .and_then(|n| n.checked_add(16 * 1024))
                    .is_none_or(|n| n > BYTES)
                {
                    return Err(query::exhausted());
                }
                let cache = vec![None; count];
                let allocated = cache
                    .capacity()
                    .checked_mul(std::mem::size_of::<Option<Degree11LocalFace>>())
                    .ok_or_else(query::exhausted)?;
                outside = outside
                    .checked_add(allocated)
                    .ok_or_else(query::exhausted)?;
                if source
                    .retained_budget_bytes()
                    .checked_add(outside)
                    .is_none_or(|n| n > BYTES)
                {
                    return Err(query::exhausted());
                }
                pending = Some(cache);
            }
            let properties = source
                .face_properties(
                    index as usize,
                    PhillipsDegree11FacePropertyOptions {
                        area_width_mm2: self.policy.linear_tolerance_mm * extent / 16.,
                        first_moment_width_mm3: self.policy.linear_tolerance_mm * extent * extent
                            / 16.,
                        centroid_width_mm: self.policy.linear_tolerance_mm / 16.,
                        maximum_work: WORK - 16_384,
                        maximum_bytes: BYTES.checked_sub(outside).ok_or_else(query::exhausted)?,
                        ..Default::default()
                    },
                    None,
                )
                .map_err(source_failure)?;
            Degree11LocalFace {
                area_mm2: properties.area_mm2,
                centroid_mm: properties.centroid_mm,
            }
        };
        let scale = placement.area_scale_bounds();
        let centroid = properties.centroid_mm.map(|v| I::new(v[0], v[1]));
        let midpoint = scale.lo.midpoint(scale.hi);
        let variation = if plane.is_some() {
            I::ZERO
        } else {
            I::point((scale - I::point(midpoint)).magnitude()) / I::point(scale.lo)
        };
        // The actual affine cofactor can reweight curved area elements. Bound
        // the centroid displacement using its complete certified variation.
        let domain = source
            .local_enclosure_bounds_mm()
            .map(|v| I::new(v[0], v[1]));
        let weighted = std::array::from_fn(|k| {
            let error = variation * I::point((domain[k] - centroid[k]).magnitude());
            let p = centroid[k] + I::new(-error.hi, error.hi);
            I::new(p.lo.max(domain[k].lo), p.hi.min(domain[k].hi))
        });
        let center = placement.bounded_point(weighted, self.policy.linear_tolerance_mm / 8.)?;
        let area = I::new(properties.area_mm2[0], properties.area_mm2[1]) * scale;
        if !area.valid()
            || area.lo <= 0.
            || area.width() > self.policy.linear_tolerance_mm * extent / 8.
        {
            return Err(unsupported(
                "Placed degree-eleven face area reaches the precision floor",
            ));
        }
        let plane = plane.map(|p| {
            let n = std::array::from_fn(|k| {
                p.u[(k + 1) % 3] * p.v[(k + 2) % 3] - p.u[(k + 2) % 3] * p.v[(k + 1) % 3]
            });
            (p.origin_mm, n)
        });
        // Publish only after successful source and occurrence validation.
        // The body owns the immutable source and policy, so transforms share
        // valid local integrals while recomputing their actual world transport.
        if cached.is_none() {
            if let Some(mut cache) = pending {
                cache[index as usize] = Some(properties);
                *self.degree11_properties.borrow_mut() = cache;
            } else {
                self.degree11_properties.borrow_mut()[index as usize] = Some(properties);
            }
        }
        Ok(Face {
            index,
            plane,
            center_mm: center,
            area_mm2: area.lo.midpoint(area.hi),
        })
    }
    pub(crate) fn degree11_faces(&self, frame: PlaneFrame) -> Result<Vec<Face>, GeometryError> {
        let source = self
            .source
            .degree11()
            .ok_or_else(|| unsupported("No degree-eleven material"))?;
        (0..source.face_count())
            .map(|i| self.degree11_face(frame, i as u32))
            .collect()
    }
}
