//! Literal affine pose of the unchanged local NPT source.
use super::*;
use crate::geometry::internal_thread_treatment_body::interval::Interval as I;
use multicad_rust_kernel::axial_pipe_thread::authored::{StockClipPlacement, StockClipPolyline};
fn point(pose: RigidTransform, p: [f64; 3]) -> [f64; 3] {
    let v = pose.vector(p);
    std::array::from_fn(|i| v[i] + pose.translation_mm[i])
}
fn reference(pose: RigidTransform, f: ReferenceFrame) -> ReferenceFrame {
    ReferenceFrame {
        origin_mm: point(pose, f.origin_mm),
        u: pose.vector(f.u),
        v: pose.vector(f.v),
    }
}
impl ClippedNptBody {
    pub(super) fn calibration(&self) -> StockClipPlacement {
        StockClipPlacement::new(self.pose.rotation, self.pose.translation_mm)
            .expect("qualified immutable source pose")
    }
    pub(in crate::geometry) fn transformed(
        &self,
        outer: RigidTransform,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        let mut body = self.clone();
        body.pose = RigidTransform {
            rotation: std::array::from_fn(|i| {
                std::array::from_fn(|j| {
                    (0..3)
                        .map(|k| outer.rotation[i][k] * self.pose.rotation[k][j])
                        .sum()
                })
            }),
            translation_mm: point(outer, self.pose.translation_mm),
        };
        body.pose.validate(policy)?;
        let calibration = StockClipPlacement::new(body.pose.rotation, body.pose.translation_mm)
            .map_err(failure)?;
        let extent = body
            .source
            .bounds_mm()
            .into_iter()
            .flatten()
            .fold(0_f64, |m, x| m.max(x.abs()));
        if (I::point(calibration.metric_error()) * I::point(extent)).hi
            > policy.linear_tolerance_mm / 4.
        {
            return Err(failure(Cause::Unresolved));
        }
        let bounds = calibration
            .point_bounds(body.source.bounds_mm())
            .map_err(failure)?;
        policy.point(bounds.map(|b| b[0]))?;
        policy.point(bounds.map(|b| b[1]))?;
        Ok(body)
    }
    pub(super) fn placed_volume(&self) -> [f64; 2] {
        if self.pose == RigidTransform::default() {
            return self.volume;
        }
        let d = self.calibration().determinant_bounds();
        let value = I::new(self.volume[0], self.volume[1]) * I::new(d[0], d[1]);
        [value.lo, value.hi]
    }
    pub(super) fn placed_faces(&self) -> Vec<Face> {
        if self.pose == RigidTransform::default() {
            return self.faces.clone();
        }
        let scale = self.calibration().area_scale_bounds();
        self.faces
            .iter()
            .zip(&self.measures)
            .map(|(face, measure)| {
                let mut face = face.clone();
                face.center_mm = point(self.pose, face.center_mm);
                let a = measure.area_mm2();
                let area = I::new(a[0], a[1]) * I::new(scale[0], scale[1]);
                face.area_mm2 = area.lo.midpoint(area.hi);
                face.plane = face.plane.map(|(p, n)| {
                    // A plane's outward normal follows the actual cofactor, not R*n.
                    let matrix = self.pose.rotation;
                    let pairs = [(1, 2), (2, 0), (0, 1)];
                    let normal = std::array::from_fn(|i| {
                        let (a, b) = pairs[i];
                        (0..3)
                            .map(|j| {
                                let (c, d) = pairs[j];
                                (matrix[a][c] * matrix[b][d] - matrix[a][d] * matrix[b][c]) * n[j]
                            })
                            .sum()
                    });
                    (point(self.pose, p), normal)
                });
                face
            })
            .collect()
    }
    pub(super) fn place_topology(&self, snapshot: &mut TopologySnapshot) {
        if self.pose == RigidTransform::default() {
            return;
        }
        for element in &mut snapshot.elements {
            for frame in &mut element.frames {
                frame.frame = reference(self.pose, frame.frame);
            }
            if let Some(geometry) = &mut element.geometry {
                match geometry {
                    AnalyticGeometry::Point { position_mm } => {
                        *position_mm = point(self.pose, *position_mm)
                    }
                    AnalyticGeometry::Plane { frame } => *frame = reference(self.pose, *frame),
                    _ => unreachable!("clipped source topology contains only points and planes"),
                }
            }
        }
    }
    pub(super) fn place_polyline(
        &self,
        source: &StockClipPolyline,
        quality: Tessellation,
        calibration: StockClipPlacement,
    ) -> Result<Vec<[f64; 3]>, MeshError> {
        if self.pose == RigidTransform::default() {
            return Ok(source.points_mm().to_vec());
        }
        let mut rounding = 0_f64;
        let mut points = Vec::with_capacity(source.points_mm().len());
        for &p in source.points_mm() {
            let bounds = calibration
                .point_bounds(p.map(|x| [x, x]))
                .map_err(|_| MeshError::UnresolvedBoundary)?;
            let literal = bounds.map(|b| b[0].midpoint(b[1]));
            let squared = (0..3).fold(I::point(0.), |e, i| {
                e + (I::point(literal[i]) - I::new(bounds[i][0], bounds[i][1])).square()
            });
            // The sum of squares is nonnegative independently of directed
            // addition's possible subnormal lower-end rounding.
            let error = I::new(squared.lo.max(0.), squared.hi).sqrt();
            rounding = rounding.max(error.hi);
            points.push(literal);
        }
        let scale = (I::point(1.) + I::point(calibration.metric_error())).sqrt();
        if (scale * I::point(source.source_chord_error_mm()) + I::point(rounding)).hi
            > quality.deflection
        {
            return Err(MeshError::UnresolvedBoundary);
        }
        Ok(points)
    }
}
