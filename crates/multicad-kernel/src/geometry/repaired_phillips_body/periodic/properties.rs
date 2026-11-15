//! Lazy actual-source face measures; display triangles never supply measures.
use super::*;
use multicad_rust_kernel::{
    phillips_drive::FacePropertyOptions, step_read::PhillipsPeriodicFaceProperties,
};
impl RepairedPhillipsBody {
    pub(crate) fn periodic_faces(&self, frame: PlaneFrame) -> Result<Vec<Face>, GeometryError> {
        let source = self.source.periodic().ok_or_else(|| {
            unsupported("Retained polynomial face properties have not been admitted")
        })?;
        let placement = self.placement(frame)?;
        let extent = source
            .local_enclosure_bounds()
            .into_iter()
            .flatten()
            .map(f64::abs)
            .fold(1_f64, f64::max);
        if self.periodic_properties.borrow().is_none() {
            // The delegated source API retains its original cumulative-cell
            // contract. It does not claim a scratch-byte receipt. This cache's
            // actual capacity is charged to every subsequent mesh query.
            let approximation = source
                .mesh_approximation(16_000_000, None)
                .map_err(source_failure)?;
            let properties = approximation
                .source_face_properties(
                    FacePropertyOptions {
                        area_tolerance_mm2: self.policy.linear_tolerance_mm * extent / 16.,
                        first_moment_tolerance_mm3: self.policy.linear_tolerance_mm
                            * extent
                            * extent
                            / 16.,
                        centroid_tolerance_mm: self.policy.linear_tolerance_mm / 16.,
                        maximum_cells: 65536,
                    },
                    16_000_000 - approximation.work_used(),
                    None,
                )
                .map_err(source_failure)?;
            let added =
                properties.capacity() * std::mem::size_of::<PhillipsPeriodicFaceProperties>();
            if self.resident_bytes().saturating_add(added) > BYTES {
                return Err(failure(
                    GeometryErrorKind::Limit,
                    "Source face cache exceeds retained bytes",
                ));
            }
            *self.periodic_properties.borrow_mut() = Some(properties);
        }
        let properties = self.periodic_properties.borrow();
        let properties = properties
            .as_ref()
            .expect("successful source cache publication");
        let mut result = Vec::with_capacity(properties.len());
        for p in properties {
            let source_face = source
                .face_source(p.source_face)
                .expect("immutable property owner");
            let plane = super::topology_planes::plane(source, source_face, placement, self.policy)?;
            let scale = placement.area_scale_bounds();
            let centroid = p.centroid_mm.map(|v| I::new(v[0], v[1]));
            let midpoint = scale.lo.midpoint(scale.hi);
            let variation = if plane.is_some() {
                I::ZERO
            } else {
                I::point((scale - I::point(midpoint)).magnitude()) / I::point(scale.lo)
            };
            // A nonrigid map weights curved surface elements differently.
            // E[w*(X-C)]/E[w] is bounded by the actual Jacobian variation
            // times each source-domain radius. Plane Jacobians are constant.
            let domain = source.local_enclosure_bounds().map(|v| I::new(v[0], v[1]));
            let weighted = std::array::from_fn(|i| {
                let error = variation * I::point((domain[i] - centroid[i]).magnitude());
                let value = centroid[i] + I::new(-error.hi, error.hi);
                I::new(value.lo.max(domain[i].lo), value.hi.min(domain[i].hi))
            });
            let center = placement.bounded_point(weighted, self.policy.linear_tolerance_mm / 8.)?;
            let area = I::new(p.area_mm2[0], p.area_mm2[1]) * placement.area_scale_bounds();
            if !area.valid()
                || area.lo <= 0.
                || area.width() > self.policy.linear_tolerance_mm * extent / 8.
            {
                return Err(unsupported(
                    "Placed source face area reaches the precision floor",
                ));
            }
            let plane = plane.map(|p| {
                let n = std::array::from_fn(|i| {
                    p.u[(i + 1) % 3] * p.v[(i + 2) % 3] - p.u[(i + 2) % 3] * p.v[(i + 1) % 3]
                });
                (p.origin_mm, n)
            });
            result.push(Face {
                index: p.source_face as u32,
                plane,
                center_mm: center,
                area_mm2: area.lo.midpoint(area.hi),
            });
        }
        Ok(result)
    }
}
