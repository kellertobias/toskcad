//! Write box stock with one oblique cylindrical cut on its retained supports.
//!
//! Every face keeps its support: a stock or cap plane, or the tool cylinder.
//! Every edge keeps its curve: a line, a circle or an ellipse in principal
//! form. A wall face bounded by two loops that both wrap the tool axis is an
//! annulus; STEP bounds it by one loop that crosses a straight seam
//! generator twice, as for every other periodic wall in this writer.
use super::*;
use crate::oblique_bore::{Arc, ObliqueBore, Support};
use std::f64::consts::TAU;

/// Whether the cut fits STEP's single-shell profile: every planar face has
/// one outer loop and every wall face is a patch or an annulus.
pub fn supports_oblique_bore(bore: &ObliqueBore) -> bool {
    Plan::new(bore).is_ok()
}

/// Write a rigidly placed oblique cut with analytic line, circle and
/// ellipse trims.
pub fn write_oblique_bore(
    bore: &ObliqueBore,
    origin: Point,
    u: Point,
    v: Point,
) -> Result<String, StepError> {
    write_one(|writer| writer.append_oblique_bore(bore, origin, u, v))
}

#[derive(Clone, Copy)]
struct PlanEdge {
    ends: [usize; 2],
    curve: Option<Arc>,
}

enum Bounds {
    /// Outer loop first, then holes.
    Planar(Vec<Vec<(usize, bool)>>),
    Patch(Vec<(usize, bool)>),
    /// One loop through the seam edge, which is traversed in both senses.
    Annulus {
        ring: Vec<(usize, bool)>,
        seam: f64,
    },
}

/// The retained boundary with the seam vertices and edges an annulus needs.
struct Plan {
    vertices: Vec<Point>,
    edges: Vec<PlanEdge>,
    faces: Vec<(Support, Bounds)>,
}

fn sub(a: Point, b: Point) -> Point {
    std::array::from_fn(|i| a[i] - b[i])
}

fn norm(a: Point) -> f64 {
    dot(a, a).sqrt()
}

fn start(edges: &[PlanEdge], (edge, forward): (usize, bool)) -> usize {
    edges[edge].ends[usize::from(!forward)]
}

impl Plan {
    fn new(bore: &ObliqueBore) -> Result<Self, StepError> {
        let tool = bore.tool();
        let [e1, e2, _] = tool.axes;
        let mut plan = Self {
            vertices: bore.vertices().to_vec(),
            edges: bore
                .edges()
                .iter()
                .map(|edge| PlanEdge {
                    ends: edge.ends,
                    curve: edge.curve,
                })
                .collect(),
            faces: Vec::new(),
        };
        let mut loops: Vec<Vec<Vec<(usize, bool)>>> =
            bore.faces().iter().map(|face| face.loops.clone()).collect();
        if loops.iter().flatten().any(|ring| ring.is_empty()) {
            return Err(StepError::Unsupported);
        }
        let angle = |p: Point| {
            let d = sub(p, tool.base_mm);
            dot(d, e2).portable_atan2(dot(d, e1))
        };
        let scale = bore
            .vertices()
            .iter()
            .flatten()
            .chain(&tool.base_mm)
            .fold(tool.radius_mm, |m, x| m.max(x.abs()));
        let close = 64. * f64::EPSILON * scale;
        let mut seams = Vec::new();
        for (index, face) in bore.faces().iter().enumerate() {
            if face.support != Support::ToolSide {
                continue;
            }
            // Signed number of turns about the tool axis of each loop.
            let turns: Vec<f64> = loops[index]
                .iter()
                .map(|ring| {
                    ring.iter()
                        .map(|&(edge, forward)| {
                            let sweep = bore.edges()[edge].curve.map_or(0., |arc| arc.sweep_rad);
                            if forward { sweep } else { -sweep }
                        })
                        .sum::<f64>()
                        / TAU
                })
                .collect();
            match turns[..] {
                [t] if t.abs() < 0.5 => {}
                [a, b] if (a.abs() - 1.).abs() < 1e-9 && (a + b).abs() < 1e-9 => {
                    seams.push(index);
                }
                _ => return Err(StepError::Unsupported),
            }
        }
        for &index in &seams {
            // The seam leaves the first vertex of the first loop along the
            // generator at its tool angle and meets the second loop there.
            let from = start(&plan.edges, loops[index][0][0]);
            let theta = angle(plan.vertices[from]);
            let target = plan.wall_point(&loops[index][1], theta, angle);
            let to = match target {
                Some(vertex) => vertex,
                None => plan.split(&mut loops, index, theta, close)?,
            };
            if norm(sub(plan.vertices[from], plan.vertices[to])) <= close {
                return Err(StepError::InvalidGeometry);
            }
            let seam = plan.edges.len();
            plan.edges.push(PlanEdge {
                ends: [from, to],
                curve: None,
            });
            let rotate = |ring: &[(usize, bool)], vertex: usize| {
                let first = ring
                    .iter()
                    .position(|&use_| start(&plan.edges, use_) == vertex)
                    .ok_or(StepError::InvalidGeometry)?;
                Ok::<_, StepError>(
                    ring[first..]
                        .iter()
                        .chain(&ring[..first])
                        .copied()
                        .collect::<Vec<_>>(),
                )
            };
            let mut ring = rotate(&loops[index][0], from)?;
            ring.push((seam, true));
            ring.extend(rotate(&loops[index][1], to)?);
            ring.push((seam, false));
            loops[index] = vec![ring];
            plan.faces.push((
                Support::ToolSide,
                Bounds::Annulus {
                    ring: loops[index][0].clone(),
                    seam: theta,
                },
            ));
        }
        let annuli = std::mem::take(&mut plan.faces);
        let mut annuli = annuli.into_iter();
        for (index, face) in bore.faces().iter().enumerate() {
            let bounds = if face.support == Support::ToolSide {
                if seams.contains(&index) {
                    annuli.next().ok_or(StepError::InvalidGeometry)?.1
                } else {
                    Bounds::Patch(loops[index][0].clone())
                }
            } else {
                let (_, normal) = bore.plane(face.support).ok_or(StepError::InvalidGeometry)?;
                let areas: Vec<f64> = loops[index]
                    .iter()
                    .map(|ring| plan.signed_area(ring, normal))
                    .collect();
                let outer: Vec<usize> = (0..areas.len()).filter(|&i| areas[i] > 0.).collect();
                let [outer] = outer[..] else {
                    return Err(StepError::Unsupported);
                };
                let mut rings = vec![loops[index][outer].clone()];
                rings.extend(
                    (0..areas.len())
                        .filter(|&i| i != outer)
                        .map(|i| loops[index][i].clone()),
                );
                Bounds::Planar(rings)
            };
            plan.faces.push((face.support, bounds));
        }
        let total: usize = plan
            .faces
            .iter()
            .map(|(_, bounds)| match bounds {
                Bounds::Planar(rings) => rings.iter().map(Vec::len).sum(),
                Bounds::Patch(ring) | Bounds::Annulus { ring, .. } => ring.len(),
            })
            .sum();
        if total > 64 * MAX_VERTICES {
            return Err(StepError::Limit);
        }
        Ok(plan)
    }

    /// The vertex of `ring` on the generator at tool angle `theta`, if any.
    fn wall_point(
        &self,
        ring: &[(usize, bool)],
        theta: f64,
        angle: impl Fn(Point) -> f64,
    ) -> Option<usize> {
        ring.iter()
            .map(|&use_| start(&self.edges, use_))
            .find(|&vertex| {
                let gap = (angle(self.vertices[vertex]) - theta).rem_euclid(TAU);
                gap.min(TAU - gap) <= 64. * f64::EPSILON
            })
    }

    /// Split the arc of the second loop of face `index` that spans `theta`,
    /// in every loop that uses it, and return the new vertex.
    fn split(
        &mut self,
        loops: &mut [Vec<Vec<(usize, bool)>>],
        index: usize,
        theta: f64,
        close: f64,
    ) -> Result<usize, StepError> {
        let (edge, arc) = loops[index][1]
            .iter()
            .find_map(|&(edge, _)| {
                let arc = self.edges[edge].curve?;
                let offset = (theta - arc.start_rad).rem_euclid(TAU);
                (offset > 0. && offset < arc.sweep_rad).then_some((edge, arc))
            })
            .ok_or(StepError::InvalidGeometry)?;
        let offset = (theta - arc.start_rad).rem_euclid(TAU);
        let point = arc.point(arc.start_rad + offset);
        let [a, b] = self.edges[edge].ends;
        if [a, b]
            .iter()
            .any(|&v| norm(sub(self.vertices[v], point)) <= close)
        {
            return Err(StepError::InvalidGeometry);
        }
        let vertex = self.vertices.len();
        self.vertices.push(point);
        let first = PlanEdge {
            ends: [a, vertex],
            curve: Some(Arc {
                sweep_rad: offset,
                ..arc
            }),
        };
        let second = PlanEdge {
            ends: [vertex, b],
            curve: Some(Arc {
                start_rad: arc.start_rad + offset,
                sweep_rad: arc.sweep_rad - offset,
                ..arc
            }),
        };
        self.edges[edge] = first;
        let tail = self.edges.len();
        self.edges.push(second);
        for ring in loops.iter_mut().flatten() {
            if let Some(at) = ring.iter().position(|&(e, _)| e == edge) {
                let forward = ring[at].1;
                let pieces = if forward {
                    [(edge, true), (tail, true)]
                } else {
                    [(tail, false), (edge, false)]
                };
                ring.splice(at..=at, pieces);
            }
        }
        Ok(vertex)
    }

    /// The middle of the widest tool-angle range that no arc of a patch
    /// loop covers.
    fn free_angle(&self, ring: &[(usize, bool)]) -> Option<f64> {
        let mut spans: Vec<[f64; 2]> = ring
            .iter()
            .filter_map(|&(edge, _)| self.edges[edge].curve)
            .map(|arc| {
                let start = arc.start_rad.rem_euclid(TAU);
                [start, start + arc.sweep_rad]
            })
            .collect();
        spans.sort_by(|a, b| a[0].total_cmp(&b[0]));
        let first = spans.first()?[0];
        let mut reach = first;
        let (mut best, mut at) = (0., None);
        for span in &spans {
            if span[0] - reach > best {
                (best, at) = (span[0] - reach, Some((reach + span[0]) / 2.));
            }
            reach = reach.max(span[1]);
        }
        if first + TAU - reach > best {
            (best, at) = (first + TAU - reach, Some((reach + first + TAU) / 2.));
        }
        (best > 0.).then_some(at?)
    }

    /// Signed area about `normal` of a planar loop with its arc segments.
    fn signed_area(&self, ring: &[(usize, bool)], normal: Point) -> f64 {
        let reference = self.vertices[start(&self.edges, ring[0])];
        ring.iter()
            .map(|&(edge, forward)| {
                let PlanEdge { ends, curve } = self.edges[edge];
                let [p, q] = if forward { ends } else { [ends[1], ends[0]] }
                    .map(|v| sub(self.vertices[v], reference));
                let chord = dot(normal, cross(p, q)) / 2.;
                let lens = curve.map_or(0., |arc| {
                    let sweep = arc.sweep_rad;
                    let sign = if forward { 1. } else { -1. };
                    sign * (sweep - sweep.portable_sin()) / 2.
                        * dot(normal, cross(arc.cos_mm, arc.sin_mm))
                });
                chord + lens
            })
            .sum()
    }
}

fn place(placement: &Placement, [x, y, z]: Point) -> Point {
    placement.point_f64([x, y], z)
}

fn arc_edge(
    file: &mut StepFile,
    [from, to]: [Vertex; 2],
    arc: Arc,
    placement: &Placement,
) -> Result<Entity, StepError> {
    let principal = arc.principal();
    let [major, minor] = principal.semiaxes_mm;
    if !(major.is_finite() && minor.is_finite() && minor > 0. && arc.sweep_rad > 0.)
        || arc.sweep_rad > TAU
    {
        return Err(StepError::InvalidGeometry);
    }
    let complete = arc.sweep_rad == TAU;
    if (from.vertex.0 == to.vertex.0) != complete {
        return Err(StepError::InvalidGeometry);
    }
    let center = file.point(place(placement, principal.center_mm));
    let normal = placement.vector(direction([0.; 3], cross(principal.u, principal.v))?);
    let circular = (major - minor).abs() <= 64. * f64::EPSILON * major;
    let curve = if circular {
        // A cap trim: its own reference keeps the vertex at angle zero.
        let reference = placement.vector(direction([0.; 3], arc.cos_mm)?);
        let axis = file.axis(center, normal, reference);
        let r = norm(arc.cos_mm);
        file.entity(format_args!("CIRCLE('',{axis},{r:.16E})"))
    } else {
        let axis = file.axis(center, normal, placement.vector(principal.u));
        file.entity(format_args!("ELLIPSE('',{axis},{major:.16E},{minor:.16E})"))
    };
    let (start, end) = (from.vertex, to.vertex);
    Ok(file.entity(format_args!("EDGE_CURVE('',{start},{end},{curve},.T.)")))
}

fn plane_reference(normal: Point) -> Result<Point, StepError> {
    let coordinate = (0..3)
        .min_by(|&a, &b| normal[a].abs().total_cmp(&normal[b].abs()))
        .expect("three coordinate axes");
    let axis = std::array::from_fn(|i| f64::from(i == coordinate));
    direction([0.; 3], cross(axis, normal))
}

pub(super) fn encode(
    file: &mut StepFile,
    placement: Placement,
    bore: &ObliqueBore,
) -> Result<Entity, StepError> {
    let plan = Plan::new(bore)?;
    let tool = bore.tool();
    let [e1, e2, w] = tool.axes;
    let vertices: Vec<_> = plan
        .vertices
        .iter()
        .map(|&point| file.vertex(place(&placement, point)))
        .collect();
    let edges = plan
        .edges
        .iter()
        .map(|edge| {
            let ends = edge.ends.map(|id| vertices[id]);
            match edge.curve {
                None => line_edge(file, ends[0], ends[1]),
                Some(arc) => arc_edge(file, ends, arc, &placement),
            }
        })
        .collect::<Result<Vec<_>, _>>()?;
    let bound = |file: &mut StepFile, ring: &[(usize, bool)], outer: bool| {
        let oriented: Vec<_> = ring
            .iter()
            .map(|&(edge, forward)| file.oriented_edge(edges[edge], forward))
            .collect();
        file.bound(&oriented, outer)
    };
    let mut faces = Vec::with_capacity(plan.faces.len());
    for (support, bounds) in &plan.faces {
        let face = match bounds {
            Bounds::Planar(rings) => {
                let (_, normal) = bore.plane(*support).ok_or(StepError::InvalidGeometry)?;
                let entities: Vec<_> = rings
                    .iter()
                    .enumerate()
                    .map(|(i, ring)| bound(file, ring, i == 0))
                    .collect();
                let origin = vertices[start(&plan.edges, rings[0][0])].cartesian;
                let normal = placement.vector(normal);
                file.face(&entities, origin, normal, plane_reference(normal)?)
            }
            Bounds::Patch(ring) | Bounds::Annulus { ring, .. } => {
                let theta = match bounds {
                    Bounds::Annulus { seam, .. } => *seam,
                    // Place the periodic seam inside the widest angle range
                    // the patch does not reach.
                    _ => plan.free_angle(ring).ok_or(StepError::Unsupported)?,
                };
                let (sin, cos) = theta.portable_sin_cos();
                let reference: Point = std::array::from_fn(|i| cos * e1[i] + sin * e2[i]);
                let entity = bound(file, ring, true);
                let center = file.point(place(&placement, tool.base_mm));
                let axis = file.axis(
                    center,
                    placement.vector(w),
                    placement.vector(direction([0.; 3], reference)?),
                );
                let radius = tool.radius_mm;
                let cylinder =
                    file.entity(format_args!("CYLINDRICAL_SURFACE('',{axis},{radius:.16E})"));
                // Material lies outside the removed tool.
                file.face_on(&[entity], cylinder, false)
            }
        };
        faces.push(face);
    }
    Ok(file.solid(&faces))
}
