//! A cylinder in certified holes and clear of the stock meets no material.
//!
//! A through bore removes the whole cylindrical column between its two
//! original stock faces, and a blind pocket the column from its entry face
//! to its floor. Along a parallel cylinder whose disk lies in a hole's disk,
//! that hole's axial interval is void. Every remaining interval of the
//! cylinder must lie in a box that is separated from, or only touches,
//! every stock cell; cuts only remove stock, so the cells enclose the
//! material. Then the cylinder shares no interior with the stock: the
//! intersection is empty and cutting it removes nothing. Ends may meet
//! within coordinate rounding, which is far below the stock's minimum
//! feature.
//! A rectangular prism with circular through holes is read as the identical
//! orthogonal stock with those bores.
use super::*;

impl Solid {
    pub(in crate::geometry) fn through_bore_clearance(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if operation == BooleanOperation::Union
            || policy.boolean_fuzzy_mm != 0.
            || self.policy != policy
            || other.policy != policy
        {
            return None;
        }
        let bored = |solid: &Self| {
            solid.orthogonal_body().is_some() || solid.drawing_bored_equivalent().is_some()
        };
        let (stock, tool) = if bored(self) && other.axial_cylinder_frame().is_some() {
            (self, other)
        } else if operation == BooleanOperation::Intersect
            && bored(other)
            && self.axial_cylinder_frame().is_some()
        {
            (other, self)
        } else {
            return None;
        };
        // A rectangular prism with circular through holes is the identical
        // point set as orthogonal stock with those through bores.
        let equivalent;
        let stock = if stock.orthogonal_body().is_some() {
            stock
        } else {
            equivalent = stock.drawing_bored_equivalent()?.ok()?;
            &equivalent
        };
        let (body, frame) = stock.orthogonal_body()?;
        let (cylinder, tool_frame) = tool.axial_cylinder_frame()?;
        let axes = [frame.u, frame.v, frame.normal()];
        let local = |point| axes.map(|axis| dot(subtract(point, frame.origin), axis));
        let start = local(tool_frame.origin);
        let end = local(transform_point(
            tool_frame,
            [0., 0., cylinder.height().as_f64()],
        ));
        let radius = cylinder.circle().radius().as_f64();
        let bounds = body.geometry.bounds();
        let numeric = 64.
            * f64::EPSILON
            * bounds
                .min
                .into_iter()
                .chain(bounds.max)
                .chain(start)
                .chain(end)
                .map(f64::abs)
                .fold(1., f64::max);
        let delta = subtract(end, start);
        let axis = (0..3).max_by(|&a, &b| delta[a].abs().total_cmp(&delta[b].abs()))?;
        if (0..3).any(|i| i != axis && delta[i].abs() > numeric) {
            return None;
        }
        let [i, j] = [(axis + 1) % 3, (axis + 2) % 3];
        let inside = |center: Point, hole_radius: f64| {
            (start[i] - center[i]).portable_hypot(start[j] - center[j]) + radius + numeric
                <= hole_radius
        };
        // Axial intervals of the tool's own column that a hole leaves void.
        let mut voids: Vec<[f64; 2]> = Vec::new();
        for bore in body.through_bores.iter() {
            let [a, b] = bore.centers_mm;
            if (0..3).all(|k| k == axis || a[k] == b[k])
                && a[axis] != b[axis]
                && inside(a, bore.radius_mm)
            {
                voids.push([a[axis].min(b[axis]), a[axis].max(b[axis])]);
            }
        }
        for pocket in body.pockets.iter() {
            let normal = body.boundary.faces.get(pocket.face)?.normal;
            if normal[axis].abs() == 1. && inside(pocket.center_mm, pocket.radius_mm) {
                let floor = pocket.floor_center_mm(normal)[axis];
                let entry = pocket.center_mm[axis];
                voids.push([entry.min(floor), entry.max(floor)]);
            }
        }
        // What remains of the tool's interval outside every void must clear
        // every stock cell, contact allowed.
        voids.sort_by(|a, b| a[0].total_cmp(&b[0]));
        let mut remaining = Vec::new();
        let mut low = start[axis].min(end[axis]);
        let high = start[axis].max(end[axis]);
        for [a, b] in voids {
            if b < low - numeric || a > high + numeric {
                continue;
            }
            if a > low + numeric {
                remaining.push([low, a.min(high)]);
            }
            low = low.max(b);
        }
        if low < high - numeric {
            remaining.push([low, high]);
        }
        let inside = remaining.iter().all(|&[a, b]| {
            let mut min = [0.; 3];
            let mut max = [0.; 3];
            min[axis] = a;
            max[axis] = b;
            for k in [i, j] {
                min[k] = start[k] - radius;
                max[k] = start[k] + radius;
            }
            body.geometry.cells().iter().all(|cell| {
                (0..3).any(|k| cell.min[k] >= max[k] - numeric || cell.max[k] <= min[k] + numeric)
            })
        });
        inside.then(|| match operation {
            BooleanOperation::Intersect => Ok(None),
            _ => self.share().map(Some),
        })
    }
}

/// A circular void column: axis index, center, radius and axial interval.
struct Column {
    axis: usize,
    center: Point,
    radius: f64,
    span: [f64; 2],
}

/// How a box relates to one void column.
enum Relation {
    /// The interiors are disjoint.
    Apart,
    /// The box's cross-section lies in the column's disk.
    Within,
    Partial,
}

fn relate(column: &Column, min: Point, max: Point, numeric: f64) -> Relation {
    let k = column.axis;
    if max[k] <= column.span[0] + numeric || min[k] >= column.span[1] - numeric {
        return Relation::Apart;
    }
    let [i, j] = [(k + 1) % 3, (k + 2) % 3];
    let c = column.center;
    // Nearest and farthest points of the cross-section rectangle.
    let near = |n: usize| c[n].clamp(min[n], max[n]) - c[n];
    let far = |n: usize| (min[n] - c[n]).abs().max((max[n] - c[n]).abs());
    if near(i).portable_hypot(near(j)) >= column.radius + numeric {
        Relation::Apart
    } else if far(i).portable_hypot(far(j)) + numeric <= column.radius {
        Relation::Within
    } else {
        Relation::Partial
    }
}

impl Solid {
    /// Intersect drilled orthogonal stock with box cells. Cuts only remove
    /// stock, so the material is the stock cells less the void columns of
    /// its through bores and blind pockets. A box apart from every column
    /// meets exactly the uncut cells, whose exact orthogonal intersection is
    /// the result. A box whose cross-section lies in a column's disk is void
    /// over that column's interval; when what remains of it is separated
    /// from, or only touches, every cell, it holds no material. Any other
    /// relation is left to the caller's refusal.
    pub(in crate::geometry) fn drilled_box_intersection(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if operation != BooleanOperation::Intersect
            || policy.boolean_fuzzy_mm != 0.
            || self.policy != policy
            || other.policy != policy
        {
            return None;
        }
        let drilled = |solid: &Self| {
            solid
                .orthogonal_body()
                .is_some_and(|(body, _)| !body.through_bores.is_empty() || !body.pockets.is_empty())
        };
        let (stock, probe) = if drilled(self) {
            (self, other)
        } else if drilled(other) {
            (other, self)
        } else {
            return None;
        };
        let (body, frame) = stock.orthogonal_body()?;
        if !body.slots.is_empty() {
            return None;
        }
        let axes = [frame.u, frame.v, frame.normal()];
        if !axes.iter().all(|axis| {
            axis.iter().all(|c| c.abs() == 1. || *c == 0.)
                && axis.iter().filter(|c| c.abs() == 1.).count() == 1
        }) {
            return None;
        }
        let cells = probe.orthogonal_cells().ok()??;
        let local = |point: Point| axes.map(|axis| dot(subtract(point, frame.origin), axis));
        let boxes: Vec<(Point, Point)> = cells
            .iter()
            .map(|cell| {
                let [a, b] = [local(cell.min), local(cell.max)];
                (
                    std::array::from_fn(|i| a[i].min(b[i])),
                    std::array::from_fn(|i| a[i].max(b[i])),
                )
            })
            .collect();
        let mut columns = Vec::new();
        for bore in body.through_bores.iter() {
            let [a, b] = bore.centers_mm;
            let axis = (0..3).find(|&k| a[k] != b[k])?;
            if (0..3).any(|k| k != axis && a[k] != b[k]) {
                return None;
            }
            columns.push(Column {
                axis,
                center: a,
                radius: bore.radius_mm,
                span: [a[axis].min(b[axis]), a[axis].max(b[axis])],
            });
        }
        for pocket in body.pockets.iter() {
            let normal = body.boundary.faces.get(pocket.face)?.normal;
            let axis = (0..3).find(|&k| normal[k].abs() == 1.)?;
            let [entry, floor] = [pocket.center_mm[axis], pocket.floor_center_mm(normal)[axis]];
            columns.push(Column {
                axis,
                center: pocket.center_mm,
                radius: pocket.radius_mm,
                span: [entry.min(floor), entry.max(floor)],
            });
        }
        let bounds = body.geometry.bounds();
        let numeric = 64.
            * f64::EPSILON
            * bounds
                .min
                .into_iter()
                .chain(bounds.max)
                .chain(boxes.iter().flat_map(|(a, b)| a.iter().chain(b).copied()))
                .map(f64::abs)
                .fold(1., f64::max);
        let clear = |min: Point, max: Point| {
            body.geometry.cells().iter().all(|cell| {
                (0..3).any(|k| cell.min[k] >= max[k] - numeric || cell.max[k] <= min[k] + numeric)
            })
        };
        let mut solid_boxes = 0;
        for &(min, max) in &boxes {
            let relations: Vec<_> = columns
                .iter()
                .map(|column| relate(column, min, max, numeric))
                .collect();
            if relations.iter().any(|r| matches!(r, Relation::Partial)) {
                return None;
            }
            let within: Vec<&Column> = columns
                .iter()
                .zip(&relations)
                .filter(|(_, r)| matches!(r, Relation::Within))
                .map(|(column, _)| column)
                .collect();
            if within.is_empty() {
                solid_boxes += 1;
                continue;
            }
            // The box is void over each containing column; the rest of it
            // along that axis must clear every cell.
            let axis = within[0].axis;
            if within.iter().any(|column| column.axis != axis) {
                return None;
            }
            let mut spans: Vec<_> = within.iter().map(|column| column.span).collect();
            spans.sort_by(|a, b| a[0].total_cmp(&b[0]));
            let mut low = min[axis];
            for [a, b] in spans {
                if a > low + numeric {
                    let mut top = max;
                    top[axis] = a.min(max[axis]);
                    let mut bottom = min;
                    bottom[axis] = low;
                    if !clear(bottom, top) {
                        return None;
                    }
                }
                low = low.max(b);
            }
            if low < max[axis] - numeric {
                let mut bottom = min;
                bottom[axis] = low;
                if !clear(bottom, max) {
                    return None;
                }
            }
        }
        if solid_boxes == 0 {
            return Some(Ok(None));
        }
        if solid_boxes != boxes.len() {
            return None;
        }
        // Every box is apart from every void: the material it meets is the
        // uncut stock's.
        let uncut = match OrthogonalBody::build(body.feature, (*body.geometry).clone(), None) {
            Ok(uncut) => uncut,
            Err(cause) => return Some(Err(cause)),
        };
        let uncut = Self {
            body: SolidBody::Orthogonal(Rc::new(uncut), frame),
            policy,
            provenance: vec![],
            rust_origins: vec![],
            retained_placement: None,
        };
        Some(uncut.boolean_geometry(probe, BooleanOperation::Intersect, policy))
    }
}
