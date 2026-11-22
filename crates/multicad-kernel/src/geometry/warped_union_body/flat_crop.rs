//! Intersect a retained warped shell with box cells where it is locally a
//! straight prism, chart or not.
//!
//! A face whose support or extent clears a cell (by the same conservative
//! tests the whole-probe relation uses) cannot meet it. When every face that
//! may meet the cell is a plane normal to one cell axis, a plane containing
//! that axis, or a cylinder along it, all edges inside the cell are parallel
//! to the axis or lie in those normal planes, because two such surfaces meet
//! only in lines along the axis. Between consecutive normal-plane heights the
//! boundary inside the cell is therefore made of segments along the axis, and
//! the material there is the plane section at mid-height, clipped to the
//! cell's footprint, swept over the slab. The exact section and its clip keep
//! every line and circular trim. A cell whose material spans more than one
//! slab, a non-circular section trim, or a footprint split into several
//! regions is left to the caller's refusal.
use super::probe::{
    SECTION_SEGMENTS, cylinder_axial_gap, cylinder_clearance, face_box_gap, loop_clearance,
    plane_clearance,
};
use super::*;
use crate::geometry::internal_thread_treatment_body::interval::Interval as B;
use crate::{CurvedEdge2d, CurvedProfile};
use multicad_rust_kernel::warped_union::{SectionClipPlane, SectionEdge, clip_section_loops};

fn unit_axis(k: usize) -> Point {
    std::array::from_fn(|i| f64::from(i == k))
}

/// Whether `a` is parallel (`true`) or perpendicular (`false`) to unit `axis`
/// up to rounding, or neither.
fn relation(a: Point, axis: Point) -> Option<bool> {
    let length = dot(a, a).sqrt();
    if !length.is_finite() || length == 0. {
        return None;
    }
    let along = dot(a, axis) / length;
    let slack = 64. * f64::EPSILON;
    if (along.abs() - 1.).abs() <= slack {
        Some(true)
    } else if along.abs() <= slack {
        Some(false)
    } else {
        None
    }
}

/// Faces of `shell` that may meet the closed box `local`, widened by `margin`.
fn touching(
    shell: &Shell,
    local: [B; 3],
    margin: f64,
) -> impl Iterator<Item = &warped_union::Face> {
    shell.faces().iter().filter(move |face| {
        let gap = match face.support {
            Support::Plane { normal } => plane_clearance(normal, face.point_mm, local, margin)
                .or_else(|| face_box_gap(shell, face, local, margin))
                .or_else(|| loop_clearance(shell, face, normal, local, margin)),
            support => cylinder_clearance(support, local, margin)
                .or_else(|| cylinder_axial_gap(shell, face, local, margin)),
        };
        gap.is_none_or(|gap| gap <= margin)
    })
}

/// Heights along axis `k` of the faces normal to it, when every face that
/// may meet the cell is a prism face along `k`.
fn prism_heights(faces: &[&warped_union::Face], k: usize) -> Option<Vec<f64>> {
    let axis = unit_axis(k);
    let mut heights = Vec::new();
    for face in faces {
        match face.support {
            Support::Plane { normal } => {
                if relation(normal, axis)? {
                    heights.push(face.point_mm[k]);
                }
            }
            Support::Cylinder { axis_direction, .. } => {
                if !relation(axis_direction, axis)? {
                    return None;
                }
            }
        }
    }
    Some(heights)
}

/// A section edge as a profile edge in the axes `[i, j]`, with the turning
/// sense of a full circle about axis `k` (`+1` counter-clockwise).
fn profile_edge(
    edge: &SectionEdge,
    [i, j]: [usize; 2],
    k: usize,
    scale: f64,
) -> Option<(CurvedEdge2d, f64)> {
    let flat = |p: Point| [p[i], p[j]];
    let [a, b] = edge.ends_mm.map(flat);
    match edge.curve {
        warped_union::Curve::Line => Some((CurvedEdge2d::Line { a, b }, 0.)),
        warped_union::Curve::Ellipse(arc) => {
            let ellipse = arc.ellipse();
            let [major, minor] = ellipse.semiaxes_mm();
            if (major - minor).abs() > 256. * f64::EPSILON * scale
                || !relation(ellipse.normal(), unit_axis(k))?
            {
                return None;
            }
            let center = flat(ellipse.center_mm());
            if a == b {
                let sense = arc.sweep_rad().signum() * ellipse.normal()[k].signum();
                return Some((
                    CurvedEdge2d::Circle {
                        center,
                        radius_mm: major,
                    },
                    sense,
                ));
            }
            let mid = flat(ellipse.point_at(arc.start_rad() + arc.sweep_rad() / 2.));
            let chord = [b[0] - a[0], b[1] - a[1]];
            let side = chord[0] * (mid[1] - a[1]) - chord[1] * (mid[0] - a[0]);
            if side == 0. {
                return None;
            }
            // A counter-clockwise arc lies to the right of its chord.
            Some((
                CurvedEdge2d::Arc {
                    center,
                    a,
                    b,
                    clockwise: side > 0.,
                },
                0.,
            ))
        }
    }
}

fn signed_area(ring: &[CurvedEdge2d]) -> f64 {
    // Shoelace over the chords plus each arc's circular segment.
    ring.iter()
        .map(|edge| match *edge {
            CurvedEdge2d::Line { a, b } => (a[0] * b[1] - b[0] * a[1]) / 2.,
            CurvedEdge2d::Arc {
                center,
                a,
                b,
                clockwise,
            } => {
                let r = (a[0] - center[0]).portable_hypot(a[1] - center[1]);
                let start = (a[1] - center[1]).portable_atan2(a[0] - center[0]);
                let end = (b[1] - center[1]).portable_atan2(b[0] - center[0]);
                let mut sweep = (end - start).rem_euclid(std::f64::consts::TAU);
                if clockwise {
                    sweep -= std::f64::consts::TAU;
                }
                let chord = (a[0] * b[1] - b[0] * a[1]) / 2.;
                chord + r * r * (sweep - sweep.portable_sin()) / 2.
            }
            CurvedEdge2d::Circle { .. } => f64::NAN,
            CurvedEdge2d::BSpline(_) => f64::NAN,
        })
        .sum()
}

impl Solid {
    pub(in crate::geometry) fn warped_prism_intersection(
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
        let (probe, (body, frame)) = match self.warped_union_body() {
            Some(found) => (other, found),
            None => (self, other.warped_union_body()?),
        };
        let axes = [frame.u, frame.v, frame.normal()];
        if !axes.iter().all(|axis| {
            axis.iter().filter(|c| c.abs() == 1.).count() == 1
                && axis.iter().all(|c| c.abs() == 1. || *c == 0.)
        }) {
            return None;
        }
        let cells = probe.orthogonal_cells().ok()??;
        let shell = &body.geometry;
        let mut pieces = Vec::new();
        for cell in &cells {
            let corners = [cell.min, cell.max].map(|corner| {
                let d: Point = std::array::from_fn(|i| corner[i] - frame.origin[i]);
                axes.map(|axis| dot(axis, d))
            });
            let low: Point = std::array::from_fn(|i| corners[0][i].min(corners[1][i]));
            let high: Point = std::array::from_fn(|i| corners[0][i].max(corners[1][i]));
            let local: [B; 3] = std::array::from_fn(|i| B::new(low[i], high[i]));
            if let Some(piece) =
                self.prism_cell(shell, frame, low, high, local, body.feature, policy)?
            {
                pieces.push(piece);
            }
        }
        match pieces.len() {
            0 => Some(Ok(None)),
            1 => Some(Ok(pieces.pop())),
            _ => None,
        }
    }

    /// `None` when the cell is not a prism zone, `Some(None)` when it holds
    /// no material, and the swept clipped section otherwise.
    #[allow(clippy::too_many_arguments)]
    fn prism_cell(
        &self,
        shell: &Shell,
        frame: PlaneFrame,
        low: Point,
        high: Point,
        local: [B; 3],
        feature: EntityId,
        policy: GeometryPolicy,
    ) -> Option<Option<Self>> {
        let margin = policy.linear_tolerance_mm;
        let faces: Vec<_> = touching(shell, local, margin).collect();
        let (k, mut heights) = (0..3)
            .rev()
            .find_map(|k| Some((k, prism_heights(&faces, k)?)))?;
        let scale = low
            .into_iter()
            .chain(high)
            .chain(frame.origin)
            .map(f64::abs)
            .fold(1., f64::max);
        let rounding = 256. * f64::EPSILON * scale;
        heights.sort_by(f64::total_cmp);
        heights.dedup_by(|a, b| (*a - *b).abs() <= rounding);
        let mut levels = vec![low[k]];
        for &height in &heights {
            // A normal plane on or outside the cell's own faces bounds no
            // slab inside it.
            if height <= low[k] + rounding || height >= high[k] - rounding {
                continue;
            }
            // One barely inside would leave a sliver below the tolerance.
            if height - levels.last()? <= margin || high[k] - height <= margin {
                return None;
            }
            levels.push(height);
        }
        levels.push(high[k]);
        let [i, j] = match k {
            0 => [1, 2],
            1 => [2, 0],
            _ => [0, 1],
        };
        let mut found = None;
        for slab in levels.windows(2) {
            let mut center: Point = std::array::from_fn(|n| low[n].midpoint(high[n]));
            center[k] = slab[0].midpoint(slab[1]);
            let loops = shell
                .section(center, unit_axis(k), margin.max(1e-7), SECTION_SEGMENTS)
                .ok()?;
            let planes = [i, j].into_iter().flat_map(|n| {
                [(low[n], 1.), (high[n], -1.)].map(|(level, sign)| {
                    let mut origin_mm = center;
                    origin_mm[n] = level;
                    SectionClipPlane {
                        origin_mm,
                        inward_normal: unit_axis(n).map(|c| c * sign),
                    }
                })
            });
            let planes: Vec<_> = planes.collect();
            let loops = clip_section_loops(
                loops,
                unit_axis(k),
                &planes,
                margin.max(1e-7),
                SECTION_SEGMENTS,
            )
            .ok()?;
            if loops.is_empty() {
                continue;
            }
            if found.is_some() {
                return None;
            }
            found = Some((slab[0], slab[1], loops));
        }
        let Some((bottom, top, loops)) = found else {
            return Some(None);
        };
        let mut outer = Vec::new();
        let mut holes = Vec::new();
        for ring in &loops {
            let edges = ring
                .iter()
                .map(|edge| profile_edge(edge, [i, j], k, scale))
                .collect::<Option<Vec<_>>>()?;
            let area = match edges.as_slice() {
                [(CurvedEdge2d::Circle { .. }, sense)] => *sense,
                _ => signed_area(
                    &edges
                        .iter()
                        .map(|(edge, _)| edge.clone())
                        .collect::<Vec<_>>(),
                ),
            };
            let edges: Vec<_> = edges.into_iter().map(|(edge, _)| edge).collect();
            if area > 0. {
                outer.push(edges);
            } else if area < 0. {
                holes.push(edges);
            } else {
                return None;
            }
        }
        let [outer] = <[_; 1]>::try_from(outer).ok()?;
        // The profile chart: `u`, `v` along local axes `i`, `j`, normal `k`.
        let axes = [frame.u, frame.v, frame.normal()];
        let origin = std::array::from_fn(|n| frame.origin[n] + bottom * axes[k][n]);
        let profile = CurvedProfile {
            frame: PlaneFrame {
                origin,
                u: axes[i],
                v: axes[j],
            },
            outer,
            holes,
        };
        let count = profile.outer.len() + profile.holes.iter().map(Vec::len).sum::<usize>();
        let origins: Vec<_> = (0..count)
            .map(|n| Origin {
                source: feature,
                role: format!("crop.side/{n}"),
            })
            .collect();
        Some(Some(
            Self::rust_extrude_normal_profile_named(
                &profile,
                top - bottom,
                feature,
                &origins,
                policy,
            )
            .ok()?,
        ))
    }
}
