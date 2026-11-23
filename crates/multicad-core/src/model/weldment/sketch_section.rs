//! Custom profile-stock sections authored as a constraint sketch (CAD-385).
//!
//! The stock owns the sketch as its canonical editable source. Its stored
//! point positions are the solved geometry: every authoring path (UI modal,
//! MCP) solves the sketch before it is written, and section dimensions are
//! plain numbers, so no later parameter change can leave the stored geometry
//! stale. The section is the one closed region the sketch encloses, with its
//! holes, found by the same exact region extraction that profiles use for
//! extrusion. Open, crossing, ambiguous or unsupported sketches are refused
//! with the extraction's own diagnostic; nothing is guessed or repaired.
use super::{SectionEdge, SectionLoop};
use crate::{
    constraint_sketch::{CurveKind, Dimension, Sketch},
    document::DocumentError,
    sketch_regions::{self, Boundary, ExactCurve, Problem},
};
use multicad_precision::PortableMath as _;
use std::f64::consts::{FRAC_PI_2, TAU};

/// Upper bound on authored section elements, matching an interactive sketch.
pub const MAX_SECTION_CURVES: usize = 256;

/// The resolved section in section millimeters: the bounding box starts at
/// (0, 0) and spans `size_mm`, as for the closed-form profile shapes.
#[derive(Clone, Debug, PartialEq)]
pub struct SketchSection {
    /// Counterclockwise outside boundary.
    pub outer: SectionLoop,
    /// Counterclockwise inner boundaries of a hollow section.
    pub holes: Vec<SectionLoop>,
    /// Net material area: outside minus holes.
    pub area_mm2: f64,
    pub size_mm: [f64; 2],
    /// Sketch coordinates of section (0, 0); adding it maps back to the sketch.
    pub offset_mm: [f64; 2],
}

fn invalid(message: impl Into<String>) -> DocumentError {
    DocumentError::Invalid(message.into())
}

/// Actionable refusal text for an extraction problem.
fn problem_message(problem: Problem, detail: &str) -> String {
    let lead = match problem {
        Problem::OpenBoundary => "The section outline is open; close every loop",
        Problem::Intersection | Problem::Overlap => {
            "The section outline crosses or overlaps itself; separate the loops"
        }
        Problem::AmbiguousJunction => {
            "More than two section edges meet at one point; keep each loop simple"
        }
        Problem::DegenerateCurve | Problem::DegenerateBoundary => {
            "The section contains a zero-size edge or loop; remove it"
        }
        Problem::InvalidSource | Problem::QueryFailure => "The section sketch cannot be read",
    };
    if detail.is_empty() {
        lead.into()
    } else {
        format!("{lead} ({detail})")
    }
}

/// Resolve and validate a section sketch.
pub fn resolve_sketch_section(sketch: &Sketch) -> Result<SketchSection, DocumentError> {
    if !sketch.texts.is_empty() {
        return Err(invalid(
            "Text cannot form a profile section; remove the text element",
        ));
    }
    if sketch.curves.len() > MAX_SECTION_CURVES {
        return Err(invalid(format!(
            "A profile section may contain at most {MAX_SECTION_CURVES} sketch curves"
        )));
    }
    if sketch
        .constraints
        .iter()
        .any(|c| matches!(c.kind.dimension(), Some(Dimension::Expression { .. })))
    {
        return Err(invalid(
            "Section dimensions must be numbers; project-parameter expressions are not supported in shared stock",
        ));
    }
    let solid: Vec<_> = sketch
        .curves
        .iter()
        .filter(|c| !c.construction)
        .cloned()
        .collect();
    if solid.is_empty() {
        return Err(invalid(
            "Draw a closed section outline; construction geometry alone has no material",
        ));
    }
    if solid
        .iter()
        .any(|c| matches!(c.kind, CurveKind::BSpline { .. }))
    {
        return Err(invalid(
            "Splines are not supported in a profile section; use lines, arcs and circles",
        ));
    }
    let extraction = sketch_regions::extract_geometry(&sketch.points, &solid);
    if let Some(d) = extraction.diagnostics.first() {
        return Err(invalid(problem_message(d.problem, &d.message)));
    }
    // A region whose outside is another region's hole is that hole's disk;
    // the section is the one region not enclosed that way.
    let holes: Vec<_> = extraction
        .regions
        .iter()
        .flat_map(|r| r.holes.iter().map(|h| h.intent.clone()))
        .collect();
    let top: Vec<_> = extraction
        .regions
        .iter()
        .filter(|r| !holes.contains(&r.outer.intent))
        .collect();
    let region = match top.as_slice() {
        [] => return Err(invalid("Draw a closed section outline")),
        [region] => *region,
        many => {
            return Err(invalid(format!(
                "A profile section must be one closed outline; this sketch has {} separate outlines",
                many.len()
            )));
        }
    };
    if extraction.regions.len() != 1 + region.holes.len() {
        return Err(invalid(
            "Islands inside a hole are not supported; a section has one outline and its holes",
        ));
    }
    let outer = boundary_loop(&region.outer)?;
    let holes = region
        .holes
        .iter()
        .map(boundary_loop)
        .collect::<Result<Vec<_>, _>>()?;
    let mut min = [f64::INFINITY; 2];
    let mut max = [f64::NEG_INFINITY; 2];
    for (lo, hi) in [bounds(&outer)] {
        for i in 0..2 {
            min[i] = min[i].min(lo[i]);
            max[i] = max[i].max(hi[i]);
        }
    }
    let size_mm = [max[0] - min[0], max[1] - min[1]];
    if !size_mm.iter().all(|v| v.is_finite() && *v > 0.) {
        return Err(invalid("The section outline has no area"));
    }
    let area_mm2 = region.outer.signed_area_mm2.abs()
        - region
            .holes
            .iter()
            .map(|h| h.signed_area_mm2.abs())
            .sum::<f64>();
    if area_mm2.is_nan() || area_mm2 <= 0. {
        return Err(invalid(
            "The section has no material between its outline and holes",
        ));
    }
    let shift = |p: [f64; 2]| [p[0] - min[0], p[1] - min[1]];
    Ok(SketchSection {
        outer: translated(&outer, shift),
        holes: holes.iter().map(|h| translated(h, shift)).collect(),
        area_mm2,
        size_mm,
        offset_mm: min,
    })
}

/// One extracted boundary as a counterclockwise section loop.
fn boundary_loop(boundary: &Boundary) -> Result<SectionLoop, DocumentError> {
    if let [edge] = boundary.edges.as_slice()
        && let ExactCurve::Circle { center, radius_mm } = edge.curve
    {
        return Ok(SectionLoop::Circle { center, radius_mm });
    }
    let mut edges = boundary
        .edges
        .iter()
        .map(|edge| match edge.curve {
            ExactCurve::Line { a, b } => Ok(SectionEdge::Line { a, b }),
            ExactCurve::Arc {
                center,
                a,
                b,
                clockwise,
            } => Ok(SectionEdge::Arc {
                center,
                a,
                b,
                clockwise,
            }),
            ExactCurve::Circle { .. } | ExactCurve::BSpline { .. } => Err(invalid(
                "A full circle must be its own loop and splines are not supported in a profile section",
            )),
        })
        .collect::<Result<Vec<_>, _>>()?;
    if boundary.signed_area_mm2 < 0. {
        edges.reverse();
        for edge in &mut edges {
            *edge = match *edge {
                SectionEdge::Line { a, b } => SectionEdge::Line { a: b, b: a },
                SectionEdge::Arc {
                    center,
                    a,
                    b,
                    clockwise,
                } => SectionEdge::Arc {
                    center,
                    a: b,
                    b: a,
                    clockwise: !clockwise,
                },
            };
        }
    }
    Ok(SectionLoop::Curved(edges))
}

fn translated(section_loop: &SectionLoop, f: impl Fn([f64; 2]) -> [f64; 2]) -> SectionLoop {
    match section_loop {
        SectionLoop::Polygon(points) => {
            SectionLoop::Polygon(points.iter().map(|p| f(*p)).collect())
        }
        SectionLoop::Circle { center, radius_mm } => SectionLoop::Circle {
            center: f(*center),
            radius_mm: *radius_mm,
        },
        SectionLoop::Curved(edges) => SectionLoop::Curved(
            edges
                .iter()
                .map(|e| match *e {
                    SectionEdge::Line { a, b } => SectionEdge::Line { a: f(a), b: f(b) },
                    SectionEdge::Arc {
                        center,
                        a,
                        b,
                        clockwise,
                    } => SectionEdge::Arc {
                        center: f(center),
                        a: f(a),
                        b: f(b),
                        clockwise,
                    },
                })
                .collect(),
        ),
    }
}

/// Angle of `p` about `center`.
fn angle(center: [f64; 2], p: [f64; 2]) -> f64 {
    (p[1] - center[1]).portable_atan2(p[0] - center[0])
}

/// Whether direction `theta` lies on the arc swept from `a` to `b`.
pub(crate) fn arc_contains(
    center: [f64; 2],
    a: [f64; 2],
    b: [f64; 2],
    clockwise: bool,
    theta: f64,
) -> bool {
    let (start, end) = (angle(center, a), angle(center, b));
    let (from, to) = if clockwise {
        (end, start)
    } else {
        (start, end)
    };
    let sweep = (to - from).rem_euclid(TAU);
    let sweep = if sweep == 0. { TAU } else { sweep };
    (theta - from).rem_euclid(TAU) <= sweep
}

/// Exact axis-aligned bounds of a loop, including arc extremes.
pub(crate) fn bounds(section_loop: &SectionLoop) -> ([f64; 2], [f64; 2]) {
    let mut min = [f64::INFINITY; 2];
    let mut max = [f64::NEG_INFINITY; 2];
    let mut add = |p: [f64; 2]| {
        for i in 0..2 {
            min[i] = min[i].min(p[i]);
            max[i] = max[i].max(p[i]);
        }
    };
    match section_loop {
        SectionLoop::Polygon(points) => points.iter().for_each(|p| add(*p)),
        SectionLoop::Circle { center, radius_mm } => {
            add([center[0] - radius_mm, center[1] - radius_mm]);
            add([center[0] + radius_mm, center[1] + radius_mm]);
        }
        SectionLoop::Curved(edges) => {
            for edge in edges {
                match *edge {
                    SectionEdge::Line { a, b } => {
                        add(a);
                        add(b);
                    }
                    SectionEdge::Arc {
                        center,
                        a,
                        b,
                        clockwise,
                    } => {
                        add(a);
                        add(b);
                        let r = (a[0] - center[0]).portable_hypot(a[1] - center[1]);
                        for k in 0..4 {
                            let theta = k as f64 * FRAC_PI_2;
                            if arc_contains(center, a, b, clockwise, theta) {
                                add([
                                    center[0] + r * theta.portable_cos(),
                                    center[1] + r * theta.portable_sin(),
                                ]);
                            }
                        }
                    }
                }
            }
        }
    }
    (min, max)
}

/// Exact largest value of `dot(p, direction)` over a loop: the support of
/// the section in that direction.
pub(crate) fn support(section_loop: &SectionLoop, direction: [f64; 2]) -> f64 {
    let d = |p: [f64; 2]| p[0] * direction[0] + p[1] * direction[1];
    let theta = direction[1].portable_atan2(direction[0]);
    let norm = direction[0].portable_hypot(direction[1]);
    match section_loop {
        SectionLoop::Polygon(points) => points
            .iter()
            .map(|p| d(*p))
            .fold(f64::NEG_INFINITY, f64::max),
        SectionLoop::Circle { center, radius_mm } => d(*center) + radius_mm * norm,
        SectionLoop::Curved(edges) => edges
            .iter()
            .map(|edge| match *edge {
                SectionEdge::Line { a, b } => d(a).max(d(b)),
                SectionEdge::Arc {
                    center,
                    a,
                    b,
                    clockwise,
                } => {
                    let r = (a[0] - center[0]).portable_hypot(a[1] - center[1]);
                    let ends = d(a).max(d(b));
                    if norm > 0. && arc_contains(center, a, b, clockwise, theta) {
                        ends.max(d(center) + r * norm)
                    } else {
                        ends
                    }
                }
            })
            .fold(f64::NEG_INFINITY, f64::max),
    }
}

/// Boundary points for comparing placed sections and measuring cut
/// perimeters: every vertex, plus points at most 1° apart along arcs and
/// circles. Never used as geometry.
pub(crate) fn sampled(section_loop: &SectionLoop) -> Vec<[f64; 2]> {
    let arc_points = |center: [f64; 2], from: f64, sweep: f64, r: f64| {
        let steps = ((sweep.abs().to_degrees()).ceil() as usize).max(1);
        (0..steps)
            .map(|i| {
                let t = from + sweep * i as f64 / steps as f64;
                [
                    center[0] + r * t.portable_cos(),
                    center[1] + r * t.portable_sin(),
                ]
            })
            .collect::<Vec<_>>()
    };
    match section_loop {
        SectionLoop::Polygon(points) => points.clone(),
        SectionLoop::Circle { center, radius_mm } => arc_points(*center, 0., TAU, *radius_mm),
        SectionLoop::Curved(edges) => edges
            .iter()
            .flat_map(|edge| match *edge {
                SectionEdge::Line { a, .. } => vec![a],
                SectionEdge::Arc {
                    center,
                    a,
                    b,
                    clockwise,
                } => {
                    let r = (a[0] - center[0]).portable_hypot(a[1] - center[1]);
                    let start = angle(center, a);
                    let mut sweep = (angle(center, b) - start).rem_euclid(TAU);
                    if sweep == 0. {
                        sweep = TAU;
                    }
                    if clockwise {
                        sweep -= TAU;
                    }
                    arc_points(center, start, sweep, r)
                }
            })
            .collect(),
    }
}

#[cfg(test)]
mod tests;
