//! Weldment member schematics (CAD-395): the chosen stock's actual
//! cross-section, rotated about the path anchor it is placed by, above a
//! side view of the member with its end extensions. Stock sizes are notes;
//! only member inputs carry callouts.
use super::super::diagram::{Diagram, Key, Line, Notation};
use multicad_core::model::{AxisAnchor, ProfileSection, SectionEdge, SectionLoop};
use multicad_precision::PortableMath as _;

pub(in crate::modeling) const ANCHOR: Key = "weld.anchor";
pub(in crate::modeling) const ROTATION: Key = "weld.rotation";
pub(in crate::modeling) const EXTEND_START: Key = "weld.extend_start";
pub(in crate::modeling) const EXTEND_END: Key = "weld.extend_end";
pub(in crate::modeling) const START: [Key; 3] = ["weld.start.x", "weld.start.y", "weld.start.z"];
pub(in crate::modeling) const END: [Key; 3] = ["weld.end.x", "weld.end.y", "weld.end.z"];
pub(in crate::modeling) const EXTENSIONS: [Key; 2] = [EXTEND_START, EXTEND_END];

#[derive(Clone, Debug, Default)]
pub(in crate::modeling) struct Values {
    /// The selected stock's section, when the stock resolves.
    pub section: Option<ProfileSection>,
    pub anchor: [AxisAnchor; 2],
    pub rotation_deg: Option<f64>,
    pub extend: [Option<f64>; 2],
    /// Point-pair path length, when both points resolve.
    pub length: Option<f64>,
    /// The member is placed between typed points (start/end callouts).
    pub points: bool,
}

fn sample(boundary: &SectionLoop) -> Vec<[f64; 2]> {
    match boundary {
        SectionLoop::Polygon(points) => points.clone(),
        SectionLoop::Circle { center, radius_mm } => (0..48)
            .map(|i| {
                let a = f64::from(i) / 48. * std::f64::consts::TAU;
                let (s, c) = a.portable_sin_cos();
                [center[0] + radius_mm * c, center[1] + radius_mm * s]
            })
            .collect(),
        SectionLoop::Curved(edges) => {
            let mut points = vec![];
            for edge in edges {
                match *edge {
                    SectionEdge::Line { a, .. } => points.push(a),
                    SectionEdge::Arc {
                        center: [cx, cy],
                        a,
                        b,
                        clockwise,
                    } => {
                        let start = (a[1] - cy).portable_atan2(a[0] - cx);
                        let end = (b[1] - cy).portable_atan2(b[0] - cx);
                        let ccw = (end - start).rem_euclid(std::f64::consts::TAU);
                        let sweep = if clockwise {
                            ccw - std::f64::consts::TAU
                        } else {
                            ccw
                        };
                        let radius = (a[0] - cx).portable_hypot(a[1] - cy);
                        points.extend((0..12).map(|i| {
                            let (s, c) = (start + sweep * f64::from(i) / 12.).portable_sin_cos();
                            [cx + radius * c, cy + radius * s]
                        }));
                    }
                }
            }
            points
        }
    }
}

fn anchor_at(anchor: AxisAnchor, size: f64) -> f64 {
    match anchor {
        AxisAnchor::Min => 0.,
        AxisAnchor::Center => size / 2.,
        AxisAnchor::Max => size,
    }
}

pub(in crate::modeling) fn anchor_label(anchor: AxisAnchor) -> &'static str {
    match anchor {
        AxisAnchor::Min => "Min",
        AxisAnchor::Center => "Center",
        AxisAnchor::Max => "Max",
    }
}

/// The section loops in diagram coordinates: rotated about the path anchor,
/// which sits at the origin.
pub(in crate::modeling) fn placed_loops(
    section: &ProfileSection,
    anchor: [AxisAnchor; 2],
    rotation_deg: f64,
) -> Vec<Vec<[f32; 2]>> {
    let pivot = [
        anchor_at(anchor[0], section.width_mm),
        anchor_at(anchor[1], section.height_mm),
    ];
    let (s, c) = rotation_deg.to_radians().portable_sin_cos();
    std::iter::once(section.outer(0.))
        .chain(section.holes())
        .map(|boundary| {
            sample(&boundary)
                .into_iter()
                .map(|[x, y]| {
                    let (dx, dy) = (x - pivot[0], y - pivot[1]);
                    [(c * dx - s * dy) as f32, (s * dx + c * dy) as f32]
                })
                .collect()
        })
        .collect()
}

/// Stock section (top) rotated about the path anchor at the origin, then
/// the member side view below. Labels sit in rows: anchor and rotation
/// above the section, path points above the bar, extensions below it.
pub(in crate::modeling) fn diagram(v: &Values) -> Diagram {
    let mut d = Diagram::new();
    let (w, h) = v
        .section
        .as_ref()
        .map_or((40., 40.), |s| (s.width_mm as f32, s.height_mm as f32));
    let s = w.max(h).max(1e-3);
    match &v.section {
        Some(section) => {
            for points in placed_loops(section, v.anchor, v.rotation_deg.unwrap_or(0.)) {
                d = d.path(&points, true, Line::Outline);
            }
            d = d.note([0., 1.65 * s], section.label());
        }
        None => {
            d = d
                .outline_rect(-w / 2., -h / 2., w / 2., h / 2., Line::Hidden)
                .note([0., 0.], "Choose profile stock");
        }
    }
    // The member path passes through the anchor at the origin.
    let r = 0.09 * s;
    d = d
        .line([-r * 2., 0.], [r * 2., 0.], Line::Center)
        .line([0., -r * 2.], [0., r * 2.], Line::Center)
        .circle([0., 0.], r * 0.5, Line::Outline, true)
        .leader(
            ANCHOR,
            [0., 0.],
            [-0.75 * s, 1.25 * s],
            Notation::Length,
            "Anchor",
        )
        .leader(
            ROTATION,
            [r * 1.6, r * 1.6],
            [1.05 * s, 0.6 * s],
            Notation::Angle,
            "Rotation",
        );
    // Side view of the member along its path, extensions dashed.
    let length = 2.4 * s;
    let real = v.length.unwrap_or(1000.).max(1e-3) as f32;
    let ext = |e: Option<f64>| {
        let e = e.unwrap_or(0.) as f32;
        if e == 0. {
            0.
        } else {
            e.signum() * (e.abs() / real * length).clamp(0.06 * length, 0.25 * length)
        }
    };
    let (es, ee) = (ext(v.extend[0]), ext(v.extend[1]));
    let (x0, x1) = (-length / 2., length / 2.);
    let top = -1.55 * s;
    let bar = 0.3 * s;
    let mid = top - bar / 2.;
    d = d
        .rect(x0 - es.min(0.), top - bar, x1 + ee.min(0.), top)
        .line([x0 - 0.3 * s, mid], [x1 + 0.3 * s, mid], Line::Center);
    if es > 0. {
        d = d.outline_rect(x0 - es, top - bar, x0, top, Line::Hidden);
    }
    if ee > 0. {
        d = d.outline_rect(x1, top - bar, x1 + ee, top, Line::Hidden);
    }
    let below = top - bar - 0.5 * s;
    d = d
        .leader(
            EXTEND_START,
            [x0 - es, mid],
            [x0 + 0.2 * s, below],
            Notation::Length,
            "Extend start",
        )
        .leader(
            EXTEND_END,
            [x1 + ee, mid],
            [x1 - 0.2 * s, below - 0.35 * s],
            Notation::Length,
            "Extend end",
        );
    if v.points {
        let above = top + 0.45 * s;
        d = d
            .leader(
                START[0],
                [x0, top],
                [x0 + 0.2 * s, above],
                Notation::Length,
                "Start",
            )
            .leader(
                END[0],
                [x1, top],
                [x1 - 0.2 * s, above + 0.35 * s],
                Notation::Length,
                "End",
            );
    }
    d.note(
        [0., below - 0.85 * s],
        match v.length {
            Some(l) => format!("Side view · path {l:.1} mm"),
            None => "Side view · follows its sketch line".into(),
        },
    )
}
