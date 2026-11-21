//! Cabinet schematics (CAD-389): a front elevation that shows which panels
//! run full length for the chosen construction, the toe-kick and every shelf
//! at its height, beside a side section showing the depth, the back
//! arrangement and the shelf setback. Construction and back choices get their
//! own small drawings so each option shows what it produces.
use super::super::diagram::{Diagram, Key, Line, Notation, legible};
use super::BackKind;
use multicad_core::model::Construction;
use std::sync::OnceLock;

pub(super) const WIDTH: Key = "cabinet.width";
pub(super) const HEIGHT: Key = "cabinet.height";
pub(super) const DEPTH: Key = "cabinet.depth";
pub(super) const TOE: Key = "cabinet.toe";
pub(super) const RECESS: Key = "cabinet.recess";
pub(super) const SETBACK: Key = "cabinet.setback";
pub(super) const CLEARANCE: Key = "cabinet.clearance";
pub(super) const BAND: [Key; 2] = ["cabinet.band", "cabinet.band_allowance"];
pub(super) const PLACEMENT: [Key; 4] = ["cabinet.x", "cabinet.y", "cabinet.z", "cabinet.rotation"];

/// Stable keys for each retained shelf slot: (underside, depth).
pub(super) fn shelf_keys(slot: u32) -> (Key, Key) {
    static KEYS: OnceLock<Vec<(String, String)>> = OnceLock::new();
    let keys = KEYS.get_or_init(|| {
        (0..multicad_core::model::MAX_SHELVES)
            .map(|i| {
                (
                    format!("cabinet.shelf.{i}.underside"),
                    format!("cabinet.shelf.{i}.depth"),
                )
            })
            .collect()
    });
    let (underside, depth) = &keys[(slot as usize).min(keys.len() - 1)];
    (underside.as_str(), depth.as_str())
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct Shelf {
    pub slot: u32,
    /// Underside above the local floor datum, where it resolves.
    pub underside: Option<f64>,
    pub depth: Option<f64>,
    pub adjustable: bool,
}

#[derive(Clone, Debug, Default)]
pub(super) struct Values {
    pub width: Option<f64>,
    pub height: Option<f64>,
    pub depth: Option<f64>,
    pub toe: Option<f64>,
    pub body: Option<f64>,
    pub back: Option<f64>,
    pub recess: Option<f64>,
    pub setback: Option<f64>,
    pub shelves: Vec<Shelf>,
    /// Legacy carcasses have no toe-kick or explicit shelves.
    pub explicit: bool,
}

pub(super) struct Choices {
    pub construction: Construction,
    pub back: BackKind,
}

/// Front elevation panels for a construction: sides, top and bottom.
fn carcass_front(
    d: Diagram,
    construction: Construction,
    w: f32,
    h: f32,
    y0: f32,
    t: f32,
) -> Diagram {
    let top = y0 + h;
    match construction {
        Construction::SidesFullHeight => d
            .rect(0., y0, t, top)
            .rect(w - t, y0, w, top)
            .rect(t, top - t, w - t, top)
            .rect(t, y0, w - t, y0 + t),
        Construction::TopBottomFullWidth => d
            .rect(0., top - t, w, top)
            .rect(0., y0, w, y0 + t)
            .rect(0., y0 + t, t, top - t)
            .rect(w - t, y0 + t, w, top - t),
    }
}

pub(super) fn construction_art(construction: Construction) -> Diagram {
    carcass_front(Diagram::new(), construction, 60., 70., 0., 7.)
}

/// Plan section through the rear of a carcass for a back arrangement.
pub(super) fn back_art(back: BackKind) -> Diagram {
    let (w, depth, t, bt, recess) = (60., 40., 6., 4., 7.);
    let d = Diagram::new()
        .rect(0., 0., t, depth)
        .rect(w - t, 0., w, depth)
        .line([-4., -6.], [w + 4., -6.], Line::Hidden)
        .note([w / 2., -12.], "front");
    match back {
        BackKind::None => d.line([t, depth], [w - t, depth], Line::Hidden),
        BackKind::Inset => d.rect(t, depth - recess - bt, w - t, depth - recess),
        BackKind::Overlay => d.rect(0., depth, w, depth + bt),
    }
}

/// The draft carcass, proportioned from its current values.
pub(super) fn carcass(choices: &Choices, v: &Values) -> Diagram {
    let w = legible(v.width, 600., 50., 1e6);
    let h = legible(v.height, 720., 0.2 * w, 6. * w);
    let depth = legible(v.depth, 560., 0.15 * w, 3. * w);
    let toe = if v.explicit {
        legible(v.toe, 100., 0., 0.4 * h)
    } else {
        0.
    };
    let t = legible(v.body, 18., 0.02 * w, 0.08 * w);
    let bt = legible(v.back, 8., 0.01 * w, 0.06 * w);
    let recess = legible(v.recess, 10., 0., 0.3 * depth);
    let setback = legible(v.setback, 20., 0., 0.5 * depth);
    let box_h = h - toe;
    let mut d = carcass_front(Diagram::new(), choices.construction, w, box_h, toe, t);
    if toe > 0. {
        d = d
            .line([0., 0.], [w, 0.], Line::Center)
            .outline_rect(0.06 * w, 0., 0.94 * w, toe, Line::Hidden)
            .dim(TOE, [w, 0.], [w, toe], 18., Notation::Length, "Toe-kick");
    }
    // Side section to the right: front at the left edge.
    let x0 = w * 1.25 + 30.;
    let side_depth = if choices.back == BackKind::Overlay {
        depth - bt
    } else {
        depth
    };
    d = d.outline_rect(x0, toe, x0 + side_depth, h, Line::Outline);
    d = match choices.back {
        BackKind::None => d,
        BackKind::Inset => d.rect(
            x0 + depth - recess - bt,
            toe + t,
            x0 + depth - recess,
            h - t,
        ),
        BackKind::Overlay => d.rect(x0 + side_depth, toe, x0 + depth, h),
    };
    if choices.back == BackKind::Inset {
        d = d.dim(
            RECESS,
            [x0 + depth - recess, h],
            [x0 + depth, h],
            12.,
            Notation::Length,
            "Back recess",
        );
    }
    let interior = (toe + t, h - t);
    for (index, shelf) in v.shelves.iter().enumerate() {
        let fallback = interior.0
            + (interior.1 - interior.0) * (index + 1) as f32 / (v.shelves.len() + 1) as f32;
        let y = shelf
            .underside
            .map(|u| u as f32)
            .filter(|u| u.is_finite())
            .unwrap_or(fallback)
            .clamp(interior.0, interior.1 - t);
        let gap = if shelf.adjustable { 0.012 * w } else { 0. };
        d = d.rect(t + gap, y, w - t - gap, y + t);
        if shelf.adjustable {
            d = d
                .circle([t + 0.03 * w, y - 0.01 * w], 0.008 * w, Line::Outline, true)
                .circle(
                    [w - t - 0.03 * w, y - 0.01 * w],
                    0.008 * w,
                    Line::Outline,
                    true,
                );
        }
        let shelf_depth = shelf
            .depth
            .map(|v| v as f32)
            .filter(|v| *v > 0.)
            .unwrap_or(side_depth - setback)
            .min(side_depth - setback);
        d = d.rect(x0 + setback, y, x0 + setback + shelf_depth, y + t);
        let (underside, shelf_depth_key) = shelf_keys(shelf.slot);
        if !shelf.adjustable {
            d = d.leader(
                underside,
                [0.5 * w, y],
                [0.5 * w, y - 0.08 * h],
                Notation::Length,
                format!("Shelf {} underside", shelf.slot + 1),
            );
        }
        if index == 0 {
            d = d
                .dim(
                    SETBACK,
                    [x0, y + t],
                    [x0 + setback, y + t],
                    12.,
                    Notation::Length,
                    "Setback",
                )
                .dim(
                    shelf_depth_key,
                    [x0 + setback, y],
                    [x0 + setback + shelf_depth, y],
                    -12.,
                    Notation::Length,
                    "Shelf depth",
                );
        }
    }
    d.dim(WIDTH, [0., toe], [w, toe], -16., Notation::Length, "Width")
        .dim(
            HEIGHT,
            [0., 0.],
            [0., h],
            18.,
            Notation::Length,
            if v.explicit { "Floor to top" } else { "Height" },
        )
        .dim(
            DEPTH,
            [x0, toe],
            [x0 + depth, toe],
            -16.,
            Notation::Length,
            "Depth",
        )
        .note([0.5 * w, h + 0.05 * h], "Front")
        .note([x0 + 0.5 * depth, h + 0.05 * h], "Side section")
}
