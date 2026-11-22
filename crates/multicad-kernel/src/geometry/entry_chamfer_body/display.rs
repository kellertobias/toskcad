//! Interval-bounded source curves; these polylines never become material.
use super::*;
impl EntryChamferBody {
    pub(in crate::geometry) fn display_edges(
        &self,
        frame: PlaneFrame,
        deflection: f64,
    ) -> Result<Vec<DisplayEdge>, GeometryError> {
        const MAX_POINTS: usize = 500_000;
        // Full placement is separately metric checked. Reserve a conservative
        // finite forward-evaluation allowance before sampling the source.
        let scale = frame.origin.into_iter().map(f64::abs).fold(
            self.geometry
                .full_length_mm()
                .max(self.geometry.source_thread().spec().major_radius_mm),
            f64::max,
        );
        let arithmetic = 128. * f64::EPSILON * scale;
        let budget = deflection / 4. - arithmetic - 2. * self.topology.maximum_source_error_mm();
        if budget <= 0. {
            return Err(failure(source::Error::Unresolved));
        }
        let mut used = 0;
        let mut output = Vec::new();
        for (logical, reference) in self.references.edges.iter().enumerate() {
            let source = self.topology.edges()[reference.pieces[0]].source();
            if matches!(source, source::SourceCurve::NeckSeam { .. }) {
                continue;
            }
            let mut points = Vec::new();
            for &edge in &reference.pieces {
                let sample = |t| {
                    if t == 0. || t == 1. {
                        Ok(self.topology.vertices()
                            [self.topology.edges()[edge].vertices()[usize::from(t == 1.)]]
                        .nominal_mm)
                    } else {
                        self.topology
                            .curve_sample(edge, t)
                            .map(|p| p.nominal_mm)
                            .map_err(failure)
                    }
                };
                if points.is_empty() {
                    points.push(position(frame, sample(0.)?));
                }
                let mut pending = vec![[0., 1.]];
                while let Some(range) = pending.pop() {
                    if used + pending.len() + points.len() >= MAX_POINTS {
                        return Err(failure(source::Error::Limit));
                    }
                    let endpoints = [sample(range[0])?, sample(range[1])?];
                    let error = self
                        .topology
                        .curve_chord_error(edge, range, endpoints)
                        .map_err(failure)?;
                    if error <= budget {
                        points.push(position(frame, endpoints[1]));
                    } else {
                        let mid = range[0].midpoint(range[1]);
                        if mid <= range[0] || mid >= range[1] {
                            return Err(failure(source::Error::Unresolved));
                        }
                        pending.push([mid, range[1]]);
                        pending.push([range[0], mid]);
                    }
                }
            }
            used += points.len();
            output.push(DisplayEdge {
                index: logical as u32,
                points_mm: points,
            });
        }
        Ok(output)
    }
}
