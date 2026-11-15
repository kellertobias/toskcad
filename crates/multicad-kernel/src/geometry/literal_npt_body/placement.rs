//! Additional exact Native occurrences retain their original STEP source.
use super::*;
use crate::geometry::internal_thread_treatment_body::interval::Interval as I;
fn ranges(value: [[f64; 2]; 3]) -> [I; 3] {
    value.map(|b| I::new(b[0], b[1]))
}
fn output(value: [I; 3]) -> [[f64; 2]; 3] {
    value.map(|i| [i.lo, i.hi])
}
impl LiteralNptBody {
    pub(in crate::geometry) fn transformed(
        &self,
        transform: RigidTransform,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        transform.validate(policy)?;
        let options = read::Options {
            maximum_bytes: OVERHEAD,
            ..Default::default()
        };
        let occurrence = match &self.occurrence {
            Some(old) => {
                old.transformed(transform.rotation, transform.translation_mm, options, None)
            }
            None => read::LiteralOccurrencePlacement::new(
                transform.rotation,
                transform.translation_mm,
                options,
                None,
            ),
        }
        .map_err(failure)?;
        let source_box = self.source.world_enclosure_mm();
        let world = occurrence.map_box_mm(source_box).map_err(failure)?;
        policy.point(world.map(|b| b[0]))?;
        policy.point(world.map(|b| b[1]))?;
        let span = source_box
            .into_iter()
            .fold(I::ZERO, |s, b| s + (I::point(b[1]) - I::point(b[0])));
        // A near-isometric literal matrix also changes area weighting. This
        // condition reserves a real centroid budget for that effect, instead
        // of pretending the additional Native matrix is exactly orthogonal.
        let delta = I::point(occurrence.metric_error());
        let deviation = I::point(occurrence.linear_l1_norm()) * span * delta / (I::ONE - delta);
        if !deviation.valid() || deviation.hi > policy.linear_tolerance_mm / 8. {
            return Err(failure(read::Error::Unsupported));
        }
        let mut vertices = Vec::with_capacity(self.source.vertices().len());
        let available = self
            .available_query_bytes()?
            .checked_sub(occurrence.retained_budget_bytes())
            .ok_or_else(|| failure(read::Error::Limit))?;
        for index in 0..self.source.vertices().len() {
            let sample = self
                .source
                .sample_vertex(
                    index,
                    read::Options {
                        maximum_bytes: available,
                        ..Default::default()
                    },
                    None,
                )
                .map_err(failure)?;
            let original = sample.position_bounds_mm();
            let exact = occurrence.map_box_mm(original).map_err(failure)?;
            let literal = exact.map(|b| b[0].midpoint(b[1]));
            let error = exact.into_iter().zip(literal).fold(I::ZERO, |s, (b, p)| {
                s + I::point((p - b[0]).abs().max((b[1] - p).abs()).next_up())
            });
            if !error.valid() || error.hi > policy.linear_tolerance_mm / 8. {
                return Err(failure(read::Error::Unsupported));
            }
            policy.point(literal)?;
            vertices.push(literal);
        }
        Ok(Self {
            feature: self.feature,
            source: self.source.clone(),
            plan: Rc::clone(&self.plan),
            occurrence: Some(occurrence),
            vertices,
            policy,
            faces: Rc::clone(&self.faces),
            display: RefCell::new(None),
            last_mesh: Cell::new(None),
            last_actual: Cell::new(None),
            mesh_builds: Cell::new(0),
            step_ready: Cell::new(None),
            last_step: Cell::new(None),
        })
    }
    pub(super) fn world_bounds(&self) -> [[f64; 2]; 3] {
        let b = self.source.world_enclosure_mm();
        self.occurrence.as_ref().map_or(b, |p| {
            p.map_box_mm(b)
                .expect("qualified immutable occurrence bounds")
        })
    }
    pub(super) fn world_volume(&self) -> [f64; 2] {
        let v = self.source.world_volume_mm3();
        let v = self.occurrence.as_ref().map_or(I::new(v[0], v[1]), |p| {
            let d = p.determinant_bounds();
            I::new(v[0], v[1]) * I::new(d[0], d[1])
        });
        [v.lo, v.hi]
    }
    pub(super) fn face_properties(
        &self,
        original: OriginalFace,
    ) -> Result<OriginalFace, GeometryError> {
        let Some(p) = &self.occurrence else {
            return Ok(original);
        };
        let center = p.map_box_mm(original.center).map_err(failure)?;
        let mut radius = I::ZERO;
        for (b, c) in self
            .source
            .world_enclosure_mm()
            .into_iter()
            .zip(ranges(original.center))
        {
            let d = I::new(b[0], b[1]) - c;
            radius = radius + I::point(d.magnitude());
        }
        let delta = I::point(p.metric_error());
        let error = I::point(p.linear_l1_norm()) * radius * delta / (I::ONE - delta);
        let center = ranges(center).map(|c| c + I::new(-error.hi, error.hi));
        if center
            .iter()
            .any(|c| !c.valid() || c.width() > self.policy.linear_tolerance_mm)
        {
            return Err(failure(read::Error::Unsupported));
        }
        let scale = p.area_scale_bounds();
        let area = I::new(original.area[0], original.area[1]) * I::new(scale[0], scale[1]);
        if !area.valid() || area.lo <= 0. {
            return Err(failure(read::Error::Unsupported));
        }
        Ok(OriginalFace {
            area: [area.lo, area.hi],
            center: output(center),
        })
    }
    pub(super) fn place_polyline(
        &self,
        source: &read::LiteralEdgePolyline,
        q: Tessellation,
    ) -> Result<Vec<[f64; 3]>, MeshError> {
        let Some(p) = &self.occurrence else {
            return Ok(source.points_mm().to_vec());
        };
        let mut out = Vec::with_capacity(source.points_mm().len());
        let mut rounding = 0f64;
        for point in source.points_mm() {
            let bound = p
                .map_box_mm(point.map(|v| [v, v]))
                .map_err(|_| MeshError::UnresolvedBoundary)?;
            let point = bound.map(|b| b[0].midpoint(b[1]));
            let error = bound.into_iter().zip(point).fold(I::ZERO, |s, (b, x)| {
                s + I::point((x - b[0]).abs().max((b[1] - x).abs()).next_up())
            });
            rounding = rounding.max(error.hi);
            out.push(point);
        }
        let error = I::point(p.linear_l1_norm()) * I::point(source.position_error_mm())
            + I::point(rounding);
        if !error.valid() || error.hi > q.deflection {
            return Err(MeshError::UnresolvedBoundary);
        }
        Ok(out)
    }
}
