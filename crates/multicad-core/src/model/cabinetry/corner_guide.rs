//! Why two panels do not form a corner joint, and what would fix it.
//!
//! A butt, miter or finger joint needs a corner block: the region inside
//! both finished panels that is one full thickness of each. That means
//! panel B must cover all of panel A's thickness, panel A must cover all of
//! panel B's thickness, and both must share some length along the corner.
//! Panels that only touch, butted face to edge, do not share that block.
//!
//! The diagnosis works in panel A's frame (x: A's thickness, y: B's
//! thickness, e: along the corner). It never changes a panel: it names the
//! exact length or width each panel lacks and where that material would
//! go, so the user can resize or move the panel with the normal controls.
use super::*;
use crate::{model::RigidPose, parameters::EvaluatedParameters};

const TOL: f64 = 1e-6;

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// A box in panel A's frame, placed in world space: `origin` is one corner
/// and `edges` are its three edge vectors.
#[derive(Clone, Debug, PartialEq)]
pub struct GuideBox {
    pub origin: [f64; 3],
    pub edges: [[f64; 3]; 3],
}
impl GuideBox {
    /// The eight corners, bit k of the index selecting edge k.
    pub fn corners(&self) -> [[f64; 3]; 8] {
        std::array::from_fn(|n| {
            std::array::from_fn(|c| {
                self.origin[c]
                    + (0..3)
                        .filter(|k| n >> k & 1 == 1)
                        .map(|k| self.edges[k][c])
                        .sum::<f64>()
            })
        })
    }
    /// The twelve edges as corner index pairs.
    pub const EDGES: [(usize, usize); 12] = [
        (0, 1),
        (2, 3),
        (4, 5),
        (6, 7),
        (0, 2),
        (1, 3),
        (4, 6),
        (5, 7),
        (0, 4),
        (1, 5),
        (2, 6),
        (3, 7),
    ];
}

/// Which input of a joint pair a panel fills.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PairRole {
    A,
    B,
}
impl PairRole {
    pub fn label(self) -> &'static str {
        match self {
            Self::A => "A",
            Self::B => "B",
        }
    }
}

/// Material one panel lacks to reach the outside corner.
#[derive(Clone, Debug, PartialEq)]
pub struct PanelExtension {
    pub panel: EntityId,
    pub role: PairRole,
    /// 0: the panel's length, 1: its width.
    pub axis: usize,
    pub current_mm: f64,
    pub amount_mm: f64,
    /// The material is added at the panel's origin end, so its origin moves.
    pub at_start: bool,
    /// Where the origin moves to when `at_start`; otherwise unchanged.
    pub origin_mm: [f64; 3],
    pub new_origin_mm: [f64; 3],
    /// The panel is a plain rectangle whose size and origin the user edits;
    /// otherwise a cabinet or box generates it.
    pub editable: bool,
    /// The missing material.
    pub region: GuideBox,
}
impl PanelExtension {
    pub fn dimension(&self) -> &'static str {
        if self.axis == 0 { "length" } else { "width" }
    }
    pub fn required_mm(&self) -> f64 {
        self.current_mm + self.amount_mm
    }
}

/// What prevents the corner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CornerProblem {
    /// The panels never share a length along the corner; `gap_mm` apart.
    NoSharedLength { gap_mm: f64 },
    /// The panels touch or stand apart; no corner block exists yet.
    NoOverlap,
    /// They overlap, but not by the full thickness of both panels.
    PartialThickness,
}

/// The corner cross-section in panel A's frame, x across A's thickness and
/// y across B's, for a diagram: each rectangle is [[x0, y0], [x1, y1]].
#[derive(Clone, Debug, PartialEq)]
pub struct CornerSection {
    pub a: [[f64; 2]; 2],
    pub b: [[f64; 2]; 2],
    /// The corner block both panels must fill.
    pub block: [[f64; 2]; 2],
    /// Missing material of A, then of B.
    pub add: Vec<(PairRole, [[f64; 2]; 2])>,
    /// World directions of the section's x and y.
    pub x_dir: [f64; 3],
    pub y_dir: [f64; 3],
}

/// The diagnosis of a refused corner, with the exact fix.
#[derive(Clone, Debug, PartialEq)]
pub struct CornerGuide {
    pub problem: CornerProblem,
    /// The corner block both panels must fill, when they share a length.
    pub block: Option<GuideBox>,
    pub extensions: Vec<PanelExtension>,
    pub section: CornerSection,
    /// One explanation naming both panels and every change.
    pub message: String,
}

/// Panel B in panel A's frame: which of A's axes B's length, width and
/// thickness run along, and B's extents along each of A's axes.
pub(super) struct BInA {
    pub axis: [usize; 3],
    pub lo: [f64; 3],
    pub hi: [f64; 3],
}
pub(super) fn b_in_a(a: &PanelBox, b: &PanelBox) -> Result<BInA, DocumentError> {
    let axes = [a.u, a.v, a.normal()];
    let d = [
        b.origin[0] - a.origin[0],
        b.origin[1] - a.origin[1],
        b.origin[2] - a.origin[2],
    ];
    let start = [dot(d, axes[0]), dot(d, axes[1]), dot(d, axes[2])];
    let mut b_axis = [usize::MAX; 3];
    let mut lo = [0.; 3];
    let mut hi = [0.; 3];
    let sizes = [b.length_mm, b.width_mm, b.thickness_mm];
    for (k, dir) in [b.u, b.v, b.normal()].into_iter().enumerate() {
        let i = (0..3)
            .find(|&i| dot(dir, axes[i]).abs() > 1. - 1e-9)
            .ok_or_else(|| invalid("Joint panels must be square to each other"))?;
        b_axis[k] = i;
        let extent = sizes[k] * dot(dir, axes[i]).signum();
        (lo[i], hi[i]) = if extent >= 0. {
            (start[i], start[i] + extent)
        } else {
            (start[i] + extent, start[i])
        };
    }
    Ok(BInA {
        axis: b_axis,
        lo,
        hi,
    })
}

fn fmt(v: f64) -> String {
    let s = format!("{v:.2}");
    let s = s.trim_end_matches('0').trim_end_matches('.');
    if s == "-0" { "0".into() } else { s.into() }
}
fn axis_name(dir: [f64; 3]) -> Option<&'static str> {
    (0..3)
        .find(|&i| dir[i].abs() > 1. - 1e-9)
        .map(|i| ["X", "Y", "Z"][i])
}

struct Side<'a> {
    part: EntityId,
    role: PairRole,
    editable: bool,
    panel: &'a PanelBox,
    /// From where the panel now is back to its panel feature's own frame,
    /// where its origin is edited (identity unless it was moved).
    back: RigidPose,
    /// The panel's extents in A's frame.
    lo: [f64; 3],
    hi: [f64; 3],
}

/// Diagnose one corner pair. `Ok(None)` means the panels share a full
/// corner block; errors are panels that cannot form a corner at all
/// (not panels, not square to each other, or parallel).
pub fn diagnose_panel_corner(
    document: &Document,
    parameters: &EvaluatedParameters,
    pair: &PanelJointPair,
) -> Result<Option<CornerGuide>, DocumentError> {
    let a = super::joint::panel_box(document, parameters, pair.a)?;
    let b = super::joint::panel_box(document, parameters, pair.b)?;
    let info = |part| {
        let name = document
            .get(part)
            .map_or_else(|| "Panel".to_owned(), |e| e.name.clone());
        let editable = panel_feature(document, part)
            .and_then(|f| read_panel(document, f).ok())
            .is_some_and(|p| matches!(p.source, PanelSource::Rectangle { .. }));
        (name, editable)
    };
    let back = |part| -> Result<RigidPose, DocumentError> {
        let feature = panel_feature(document, part).expect("checked panel");
        super::joint::panel_pose_after(document, parameters, part, feature)?.inverse()
    };
    let (a_name, a_editable) = info(pair.a);
    let (b_name, b_editable) = info(pair.b);
    guide([
        (pair.a, &a_name, a_editable, &a, back(pair.a)?),
        (pair.b, &b_name, b_editable, &b, back(pair.b)?),
    ])
}

fn guide(
    [
        (a_part, a_name, a_editable, a, a_back),
        (b_part, b_name, b_editable, b, b_back),
    ]: [(EntityId, &str, bool, &PanelBox, RigidPose); 2],
) -> Result<Option<CornerGuide>, DocumentError> {
    let BInA {
        axis: b_axis,
        lo: b_lo,
        hi: b_hi,
    } = b_in_a(a, b)?;
    let (x, y) = (2, b_axis[2]);
    if y == x {
        return Err(invalid(
            "Joint panels are parallel; a corner joint needs panels at a right angle",
        ));
    }
    let e = 3 - x - y;
    let axes = [a.u, a.v, a.normal()];
    let world = |l: [f64; 3]| -> [f64; 3] {
        std::array::from_fn(|c| a.origin[c] + (0..3).map(|i| l[i] * axes[i][c]).sum::<f64>())
    };
    let to_box = |lo: [f64; 3], hi: [f64; 3]| GuideBox {
        origin: world(lo),
        edges: std::array::from_fn(|i| axes[i].map(|c| c * (hi[i] - lo[i]))),
    };
    let a_hi = [a.length_mm, a.width_mm, a.thickness_mm];
    let sa = Side {
        part: a_part,
        role: PairRole::A,
        editable: a_editable,
        panel: a,
        back: a_back,
        lo: [0.; 3],
        hi: a_hi,
    };
    let sb = Side {
        part: b_part,
        role: PairRole::B,
        editable: b_editable,
        panel: b,
        back: b_back,
        lo: b_lo,
        hi: b_hi,
    };
    let overlap = |i: usize| b_hi[i].min(a_hi[i]) - b_lo[i].max(0.);
    let section_rect = |s: &Side| [[s.lo[x], s.lo[y]], [s.hi[x], s.hi[y]]];
    let block_lo_hi = |i: usize| match i {
        _ if i == x => (0., a.thickness_mm),
        _ if i == y => (b_lo[y], b_hi[y]),
        _ => (b_lo[e].max(0.), b_hi[e].min(a_hi[e])),
    };
    let mut section = CornerSection {
        a: section_rect(&sa),
        b: section_rect(&sb),
        block: [[0., b_lo[y]], [a.thickness_mm, b_hi[y]]],
        add: vec![],
        x_dir: axes[x],
        y_dir: axes[y],
    };
    if overlap(e) <= TOL {
        let gap = (b_lo[e] - a_hi[e]).max(-b_hi[e]).max(0.);
        return Ok(Some(CornerGuide {
            problem: CornerProblem::NoSharedLength { gap_mm: gap },
            block: None,
            extensions: vec![],
            section,
            message: format!(
                "{a_name} (A) and {b_name} (B) do not overlap at a corner: they share no length along the corner{}. Move one panel along the corner until they overlap there; a corner joint is only cut where both panels run side by side.",
                if gap > TOL {
                    format!(" and are {} mm apart along it", fmt(gap))
                } else {
                    String::new()
                }
            ),
        }));
    }
    // B must span A's full thickness (x); A must span B's full thickness (y).
    let mut wants = vec![];
    for (side, i, need) in [(&sb, x, (0., a.thickness_mm)), (&sa, y, (b_lo[y], b_hi[y]))] {
        if side.lo[i] - need.0 > TOL {
            wants.push((side, i, -1., side.lo[i] - need.0));
        }
        if need.1 - side.hi[i] > TOL {
            wants.push((side, i, 1., need.1 - side.hi[i]));
        }
    }
    if wants.is_empty() {
        return Ok(None);
    }
    wants.sort_by_key(|w| w.0.role == PairRole::B);
    let problem = if overlap(x) > TOL && overlap(y) > TOL {
        CornerProblem::PartialThickness
    } else {
        CornerProblem::NoOverlap
    };
    let mut extensions = vec![];
    for (side, i, sign, amount) in wants {
        let p = side.panel;
        let k = local_axis(side.role, &b_axis, i);
        let dir = [p.u, p.v][k];
        let local = sign * dot(dir, axes[i]).signum();
        let at_start = local < 0.;
        let delta = if at_start {
            dir.map(|c| -c * amount)
        } else {
            [0.; 3]
        };
        let mut lo = side.lo;
        let mut hi = side.hi;
        if sign < 0. {
            (lo[i], hi[i]) = (side.lo[i] - amount, side.lo[i]);
        } else {
            (lo[i], hi[i]) = (side.hi[i], side.hi[i] + amount);
        }
        section
            .add
            .push((side.role, [[lo[x], lo[y]], [hi[x], hi[y]]]));
        extensions.push(PanelExtension {
            panel: side.part,
            role: side.role,
            axis: k,
            current_mm: [p.length_mm, p.width_mm][k],
            amount_mm: amount,
            at_start,
            origin_mm: side.back.point(p.origin),
            new_origin_mm: side
                .back
                .point(std::array::from_fn(|c| p.origin[c] + delta[c])),
            editable: side.editable,
            region: to_box(lo, hi),
        });
    }
    let (blo, bhi): ([f64; 3], [f64; 3]) = {
        let r: [(f64, f64); 3] = std::array::from_fn(block_lo_hi);
        (r.map(|v| v.0), r.map(|v| v.1))
    };
    let headline = match problem {
        CornerProblem::PartialThickness => format!(
            "The corner block must be the full thickness of both panels: {a_name} (A) and {b_name} (B) overlap only partly"
        ),
        _ => format!(
            "{a_name} (A) and {b_name} (B) do not overlap at a corner: they only touch or stand apart, so there is no shared corner block to divide"
        ),
    };
    let steps: Vec<String> = extensions
        .iter()
        .map(|x| {
            let name = if x.role == PairRole::A {
                a_name
            } else {
                b_name
            };
            let other = if x.role == PairRole::A {
                b_name
            } else {
                a_name
            };
            let verb = if x.axis == 0 { "lengthen" } else { "widen" };
            let mut s = format!(
                "{verb} {name} ({}) by {} mm at its {} ({} {} → {} mm)",
                x.role.label(),
                fmt(x.amount_mm),
                if x.at_start { "start" } else { "end" },
                x.dimension(),
                fmt(x.current_mm),
                fmt(x.required_mm()),
            );
            if x.at_start {
                let moved: Vec<String> = (0..3)
                    .filter(|&c| (x.new_origin_mm[c] - x.origin_mm[c]).abs() > TOL)
                    .map(|c| {
                        format!(
                            "{} {} → {} mm",
                            ["X", "Y", "Z"][c],
                            fmt(x.origin_mm[c]),
                            fmt(x.new_origin_mm[c])
                        )
                    })
                    .collect();
                s += &format!(" and move its origin {}", moved.join(", "));
            }
            s += &format!(" so it reaches the outside face of {other}");
            if !x.editable {
                s += " (its cabinet or box generates it: change the size there)";
            }
            s
        })
        .collect();
    let along = axis_name(axes[e]).map_or(String::new(), |n| format!(" along {n}"));
    let message = format!(
        "{headline}. To form the {} × {} mm corner block{along}, {}. Nothing is changed for you; the joint is refused until the panels overlap.",
        fmt(a.thickness_mm),
        fmt(b.thickness_mm),
        steps.join("; and ")
    );
    Ok(Some(CornerGuide {
        problem,
        block: Some(to_box(blo, bhi)),
        extensions,
        section,
        message,
    }))
}

/// The panel-local axis (0 length, 1 width) of a pair panel that runs
/// along A's frame axis `i`; never the thickness, which runs along x or y.
fn local_axis(role: PairRole, b_axis: &[usize; 3], i: usize) -> usize {
    match role {
        PairRole::A => i,
        PairRole::B => b_axis.iter().position(|&j| j == i).expect("an axis of B"),
    }
}
