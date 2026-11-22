//! Recognition earns every face and oriented shared edge from the literal file.
use super::*;
use crate::axial_profile_stack::RectangularBlindCounterbore;
impl Shell<'_, '_> {
    pub(in crate::step_read::import) fn rectangular_blind_counterbore(
        &self,
        faces: &[Face],
        readback: &Readback,
    ) -> Option<ImportSolid> {
        if self.points.len() != 24 || self.edges.len() != 44 || faces.len() != 22 {
            return None;
        }
        let arc = readback.arcs.first()?;
        for sign in [1., -1.] {
            let axis = arc.normal.map(|v| v * sign);
            let side = faces.iter().find_map(|f| match f {
                Face::Planar { normal, .. } if dot(*normal, axis).abs() <= FRAME_EPS => {
                    Some(*normal)
                }
                _ => None,
            })?;
            let [u, v] = frame_axes(axis, side);
            let low = self
                .points
                .values()
                .map(|&p| dot(sub(p, arc.center_mm), axis))
                .min_by(f64::total_cmp)
                .expect("the candidate has exactly 24 vertices");
            let frame = Frame {
                origin: std::array::from_fn(|i| arc.center_mm[i] + low * axis[i]),
                axes: [u, v, axis],
            };
            if let Some(stack) = self.match_blind_counterbore(faces, readback, &frame) {
                return Some(ImportSolid {
                    origin_mm: frame.origin.map(zero),
                    u,
                    v,
                    height_mm: stack.bounds_mm()[1][2],
                    profile: ImportProfile::AxialProfileStack { stack },
                });
            }
        }
        None
    }
    fn match_blind_counterbore(
        &self,
        faces: &[Face],
        readback: &Readback,
        frame: &Frame,
    ) -> Option<ProfileStack> {
        let mut levels = Levels::from_points(&self.points, frame, self.tolerance);
        if levels.values.len() != 4 {
            return None;
        }
        levels.values[0] = 0.;
        let mut min = [f64::INFINITY; 2];
        let mut max = [f64::NEG_INFINITY; 2];
        for point in self.points.values() {
            let p = frame.project(*point);
            for i in 0..2 {
                min[i] = min[i].min(p[i]);
                max[i] = max[i].max(p[i]);
            }
        }
        let mut radii: Vec<_> = readback.arcs.iter().map(|a| a.radius_mm).collect();
        radii.sort_by(f64::total_cmp);
        let bore = *radii.first()?;
        let seat = radii[radii.len() - 1];
        let minimum = self.tolerance.next_up();
        let stack = ProfileStack::blind_counterbore(
            RectangularBlindCounterbore {
                min_mm: min,
                max_mm: max,
                height_mm: levels.values[3],
                bore_radius_mm: bore,
                bore_depth_mm: levels.values[2],
                seat_radius_mm: seat,
                seat_depth_mm: levels.values[1],
                minimum_wall_mm: minimum,
                minimum_floor_mm: minimum,
            },
            self.tolerance,
            minimum,
        )
        .ok()?;
        self.match_stack(
            &stack,
            frame,
            &levels,
            faces,
            &self.indexed_edges(faces.len()),
        )?;
        stack_metrics_match(&stack, frame, readback, self.tolerance).then_some(stack)
    }
}

#[cfg(test)]
mod tests;
