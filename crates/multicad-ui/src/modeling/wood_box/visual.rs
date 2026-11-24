//! Wood box form schematic (CAD-395): an oblique view of the outside
//! envelope with each present face shaded and each omitted face dashed,
//! dimensioned by the staged width, depth and height. The panel layout
//! itself stays in the native resolution and the viewport guide.
use super::super::diagram::{Diagram, Key, Line, Notation, legible};
use multicad_core::model::WoodBoxFace;

pub(in crate::modeling) const SIZE: [Key; 3] = ["box.width", "box.depth", "box.height"];
pub(in crate::modeling) const ORIGIN: [Key; 3] = ["box.origin.x", "box.origin.y", "box.origin.z"];
pub(in crate::modeling) const ROTATION: Key = "box.rotation";
/// Fields of the named placement section.
pub(in crate::modeling) const PLACEMENT: [Key; 4] = [ORIGIN[0], ORIGIN[1], ORIGIN[2], ROTATION];

#[derive(Clone, Debug, Default)]
pub(in crate::modeling) struct Values {
    pub size: [Option<f64>; 3],
    pub faces: Vec<WoodBoxFace>,
    /// Panel stock thickness, when stock is chosen.
    pub thickness: Option<f64>,
}

/// Oblique projection: X right, Z up, depth receding up-right at half scale.
fn at(p: [f32; 3]) -> [f32; 2] {
    [p[0] + 0.5 * p[1], p[2] + 0.35 * p[1]]
}

pub(in crate::modeling) fn diagram(v: &Values) -> Diagram {
    let w = legible(v.size[0], 600., 1e-3, 1e7);
    let dep = legible(v.size[1], 400., 0.15 * w, 2. * w);
    let h = legible(v.size[2], 300., 0.15 * w, 2. * w);
    let corner = |x: bool, y: bool, z: bool| {
        at([
            if x { w } else { 0. },
            if y { dep } else { 0. },
            if z { h } else { 0. },
        ])
    };
    // Each face: its four corners, whether it faces the viewer, a label spot.
    let face = |f: WoodBoxFace| -> ([[f32; 2]; 4], bool) {
        let c = &corner;
        match f {
            WoodBoxFace::Front => (
                [
                    c(false, false, false),
                    c(true, false, false),
                    c(true, false, true),
                    c(false, false, true),
                ],
                true,
            ),
            WoodBoxFace::Back => (
                [
                    c(false, true, false),
                    c(true, true, false),
                    c(true, true, true),
                    c(false, true, true),
                ],
                false,
            ),
            WoodBoxFace::Left => (
                [
                    c(false, false, false),
                    c(false, true, false),
                    c(false, true, true),
                    c(false, false, true),
                ],
                false,
            ),
            WoodBoxFace::Right => (
                [
                    c(true, false, false),
                    c(true, true, false),
                    c(true, true, true),
                    c(true, false, true),
                ],
                true,
            ),
            WoodBoxFace::Bottom => (
                [
                    c(false, false, false),
                    c(true, false, false),
                    c(true, true, false),
                    c(false, true, false),
                ],
                false,
            ),
            WoodBoxFace::Top => (
                [
                    c(false, false, true),
                    c(true, false, true),
                    c(true, true, true),
                    c(false, true, true),
                ],
                true,
            ),
        }
    };
    let mut d = Diagram::new();
    // Hidden faces first, then the visible ones over them.
    for visible in [false, true] {
        for f in WoodBoxFace::ALL {
            let (points, faces_viewer) = face(f);
            if faces_viewer != visible {
                continue;
            }
            let present = v.faces.contains(&f);
            d = if present && visible {
                d.solid(&points)
            } else {
                d.path(
                    &points,
                    true,
                    if present { Line::Outline } else { Line::Hidden },
                )
            };
            let centre = [
                points.iter().map(|p| p[0]).sum::<f32>() / 4.,
                points.iter().map(|p| p[1]).sum::<f32>() / 4.,
            ];
            if !present {
                d = d.note(centre, format!("{} open", f.label()));
            }
        }
    }
    let present = v.faces.len();
    d.dim(
        SIZE[0],
        corner(false, false, false),
        corner(true, false, false),
        -16.,
        Notation::Length,
        "Width",
    )
    .dim(
        SIZE[2],
        corner(false, false, false),
        corner(false, false, true),
        16.,
        Notation::Length,
        "Height",
    )
    .dim(
        SIZE[1],
        corner(true, false, false),
        corner(true, true, false),
        -16.,
        Notation::Length,
        "Depth",
    )
    .leader(
        ORIGIN[0],
        corner(false, false, false),
        [0.05 * w, -0.32 * w],
        Notation::Length,
        "Origin",
    )
    .note(
        [0.5 * w + 0.25 * dep, h + 0.35 * dep + 0.16 * h],
        match v.thickness {
            Some(t) => format!("Outside envelope · {present} of 6 faces · {t} mm panels"),
            None => format!("Outside envelope · {present} of 6 faces"),
        },
    )
}
