//! Publish exact scene curves as drawing vectors within the output budget.
use super::arrange::flattened;
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
    fn stops(
        &self,
        path: &dyn Path,
        [from, to]: [f64; 2],
        tolerance: f64,
    ) -> Result<Vec<f64>, GeometryError> {
        let [lo, hi] = [from.min(to), from.max(to)];
        let mut stops = vec![lo];
        stops.extend(
            flattened(path, tolerance)?
                .into_iter()
                .filter(|t| *t > lo && *t < hi),
        );
        stops.push(hi);
        Ok(stops)
    }

    /// Chords of a range of a path in the paper, for hatch regions.
    pub(in crate::drawing) fn chords(
        &self,
        path: &dyn Path,
        range: [f64; 2],
        tolerance: f64,
    ) -> Result<Vec<[P2; 2]>, GeometryError> {
        let points: Vec<P2> = self
            .stops(path, range, tolerance)?
            .into_iter()
            .map(|t| {
                let p = path.point(t);
                [p[0], p[1]]
            })
            .collect();
        Ok(points.windows(2).map(|pair| [pair[0], pair[1]]).collect())
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
                let stops = self.stops(path, range, self.options.chord_tolerance_mm / 2.)?;
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

    /// Join collinear straight vectors of one source that touch or overlap:
    /// the many rim pieces of an edge-on cap publish as its covered
    /// stretches. Hidden stretches lose whatever a visible one covers.
    fn merge_lines(&mut self, eps: f64) {
        struct Group {
            source: usize,
            origin: P2,
            direction: P2,
            spans: [Vec<[f64; 2]>; 3],
        }
        let rank = |role: VectorRole| match role {
            VectorRole::Visible => Some(0),
            VectorRole::Hidden => Some(1),
            VectorRole::Section => Some(2),
            VectorRole::Hatch => None,
        };
        let mut lines: Vec<(usize, usize, [P2; 2])> = Vec::new();
        let mut rest = Vec::new();
        for vector in std::mem::take(&mut self.vectors) {
            match (vector.curve.clone(), rank(vector.role)) {
                (VectorCurve::Line { a, b }, Some(role)) => {
                    lines.push((vector.source_index, role, [a, b]));
                }
                _ => rest.push(vector),
            }
        }
        // Longest first, so that each group is anchored on a well-conditioned line.
        let length = |[a, b]: &[P2; 2]| (b[0] - a[0]).hypot(b[1] - a[1]);
        lines.sort_by(|x, y| length(&y.2).total_cmp(&length(&x.2)));
        let mut groups: Vec<Group> = Vec::new();
        for (source, role, [a, b]) in lines {
            let station = |group: &Group, p: P2| {
                let offset = [p[0] - group.origin[0], p[1] - group.origin[1]];
                [
                    offset[0] * group.direction[0] + offset[1] * group.direction[1],
                    offset[0] * group.direction[1] - offset[1] * group.direction[0],
                ]
            };
            let found = groups.iter().position(|group| {
                group.source == source && [a, b].iter().all(|p| station(group, *p)[1].abs() <= eps)
            });
            let index = found.unwrap_or_else(|| {
                let d = length(&[a, b]);
                groups.push(Group {
                    source,
                    origin: a,
                    direction: [(b[0] - a[0]) / d, (b[1] - a[1]) / d],
                    spans: Default::default(),
                });
                groups.len() - 1
            });
            let [s, t] = [station(&groups[index], a)[0], station(&groups[index], b)[0]];
            groups[index].spans[role].push([s.min(t), s.max(t)]);
        }
        let join = |mut spans: Vec<[f64; 2]>| {
            spans.sort_by(|a, b| a[0].total_cmp(&b[0]));
            let mut merged: Vec<[f64; 2]> = Vec::new();
            for [a, b] in spans {
                match merged.last_mut() {
                    Some(last) if a <= last[1] + eps => last[1] = last[1].max(b),
                    _ => merged.push([a, b]),
                }
            }
            merged
        };
        for group in groups {
            let [visible, hidden, section] = group.spans.map(join);
            let mut uncovered = Vec::new();
            for [mut a, b] in hidden {
                for [c, d] in &visible {
                    if *d <= a + eps || *c >= b - eps {
                        continue;
                    }
                    if *c > a + eps {
                        uncovered.push([a, *c]);
                    }
                    a = a.max(*d);
                }
                if b > a + eps {
                    uncovered.push([a, b]);
                }
            }
            let roles = [VectorRole::Visible, VectorRole::Hidden, VectorRole::Section];
            for (role, spans) in roles.into_iter().zip([visible, uncovered, section]) {
                for span in spans {
                    let [a, b] = span.map(|s| {
                        [
                            group.origin[0] + s * group.direction[0],
                            group.origin[1] + s * group.direction[1],
                        ]
                    });
                    rest.push(ProjectedVector {
                        source_index: group.source,
                        role,
                        curve: VectorCurve::Line { a, b },
                        preview_mm: vec![a, b],
                    });
                }
            }
        }
        self.vectors = rest;
    }

    /// Stable order: by source, then visible, hidden, section and hatch.
    pub(in crate::drawing) fn finish(mut self, eps: f64) -> Vec<ProjectedVector> {
        self.merge_lines(eps);
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
