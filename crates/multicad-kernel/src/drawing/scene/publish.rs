//! Publish exact scene curves as drawing vectors within the output budget.
use super::*;
use multicad_rust_kernel::section_hatch;
use std::{collections::HashMap, f64::consts::FRAC_PI_2};

pub(in crate::drawing) struct Output {
    options: DrawingProjection,
    vectors: Vec<ProjectedVector>,
    curves: usize,
    points: usize,
    hatch_budget: section_hatch::Budget,
}

fn limit() -> GeometryError {
    err(
        GeometryErrorKind::Limit,
        "Rust scene drawing exceeded its curve or point budget",
    )
}

impl Output {
    pub(in crate::drawing) fn new(options: DrawingProjection) -> Self {
        Self {
            options,
            vectors: Vec::new(),
            curves: options.maximum_curves as usize,
            points: options.maximum_points as usize,
            hatch_budget: section_hatch::Budget::default(),
        }
    }

    fn push(
        &mut self,
        source_index: usize,
        role: VectorRole,
        curve: VectorCurve,
        preview_mm: Vec<P2>,
    ) -> Result<(), GeometryError> {
        self.curves = self.curves.checked_sub(1).ok_or_else(limit)?;
        self.points = self
            .points
            .checked_sub(preview_mm.len())
            .ok_or_else(limit)?;
        self.vectors.push(ProjectedVector {
            source_index,
            role,
            curve,
            preview_mm,
        });
        Ok(())
    }

    /// Parameters of the path's certified chords within a range.
    fn stops(&self, path: &dyn Path, [from, to]: [f64; 2], tolerance: f64) -> Vec<f64> {
        let [lo, hi] = [from.min(to), from.max(to)];
        let mut stops = vec![lo];
        stops.extend(
            path.flatten(tolerance)
                .into_iter()
                .filter(|t| *t > lo && *t < hi),
        );
        stops.push(hi);
        stops
    }

    /// Chords of a range of a path in the paper, for hatch regions.
    pub(in crate::drawing) fn chords(
        &self,
        path: &dyn Path,
        range: [f64; 2],
        tolerance: f64,
    ) -> Vec<[P2; 2]> {
        let points: Vec<P2> = self
            .stops(path, range, tolerance)
            .into_iter()
            .map(|t| {
                let p = path.point(t);
                [p[0], p[1]]
            })
            .collect();
        points.windows(2).map(|pair| [pair[0], pair[1]]).collect()
    }

    /// Publish a range of a path: an exact line or conic arc where the paper
    /// curve is one, otherwise certified chords of the exact path.
    pub(in crate::drawing) fn curve(
        &mut self,
        source: usize,
        role: VectorRole,
        path: &dyn Path,
        range: [f64; 2],
    ) -> Result<(), GeometryError> {
        match path.paper(range) {
            Some(Curve::Line([a, b])) => {
                self.push(source, role, VectorCurve::Line { a, b }, vec![a, b])
            }
            Some(Curve::Ellipse(ellipse)) => {
                let [start, end] = ellipse.trim;
                let step = 4.
                    * (self.options.chord_tolerance_mm / (2. * ellipse.major))
                        .sqrt()
                        .min(1.)
                        .asin();
                let segments = (((end - start) / step.min(FRAC_PI_2)).ceil() as usize).max(1);
                if segments >= self.points {
                    return Err(limit());
                }
                let preview = (0..=segments)
                    .map(|i| {
                        ellipse.point((end - start).mul_add(i as f64 / segments as f64, start))
                    })
                    .collect();
                self.push(source, role, curved::ellipse_vector(ellipse), preview)
            }
            None => {
                let stops = self.stops(path, range, self.options.chord_tolerance_mm / 2.);
                if stops.len() > self.curves + 1 {
                    return Err(limit());
                }
                let points: Vec<P2> = stops
                    .into_iter()
                    .map(|t| {
                        let p = path.point(t);
                        [p[0], p[1]]
                    })
                    .collect();
                for pair in points.windows(2) {
                    let [a, b] = [pair[0], pair[1]];
                    self.push(source, role, VectorCurve::Line { a, b }, vec![a, b])?;
                }
                Ok(())
            }
        }
    }

    pub(in crate::drawing) fn hatch(
        &mut self,
        source: usize,
        hatch: DrawingHatch,
        boundary: &[[P2; 2]],
        tolerance: f64,
        weld: f64,
    ) -> Result<(), GeometryError> {
        // Region boundaries meet at shared exact events that were evaluated
        // on different paths; give every such event one representative.
        let mut known: HashMap<[i64; 2], Vec<P2>> = HashMap::new();
        let mut snap = |point: P2| {
            let [x, y] = point.map(|v| (v / weld).floor() as i64);
            for dx in -1..=1 {
                for dy in -1..=1 {
                    if let Some(found) = known.get(&[x + dx, y + dy]).and_then(|cell| {
                        cell.iter()
                            .find(|other| (other[0] - point[0]).hypot(other[1] - point[1]) <= weld)
                    }) {
                        return *found;
                    }
                }
            }
            known.entry([x, y]).or_default().push(point);
            point
        };
        let welded: Vec<[P2; 2]> = boundary
            .iter()
            .map(|&[a, b]| [snap(a), snap(b)])
            .filter(|[a, b]| a != b)
            .collect();
        let boundary = &welded;
        let strokes = section::publish_hatch(
            source,
            hatch,
            boundary,
            tolerance,
            &mut self.curves,
            &mut self.points,
            &mut self.hatch_budget,
        )?;
        self.vectors.extend(strokes);
        Ok(())
    }

    /// Stable order: by source, then visible, hidden, section and hatch.
    pub(in crate::drawing) fn finish(mut self) -> Vec<ProjectedVector> {
        self.vectors.sort_by_key(|vector| {
            (
                vector.source_index,
                match vector.role {
                    VectorRole::Visible => 0,
                    VectorRole::Hidden => 1,
                    VectorRole::Section => 2,
                    VectorRole::Hatch => 3,
                },
            )
        });
        self.vectors
    }
}
