//! Architecture schematics (CAD-395): a stair plan with its unrolled
//! walking-line elevation, and a floor slab or pedestal plan with its
//! section at the level floor. Diagrams explain the staged fields only;
//! the stair layout mirrors the native flight split for illustration.
use super::super::diagram::{Diagram, Key, Line, Notation, legible};
use multicad_core::model::{StairForm, Turn};

pub(in crate::modeling) const SLAB_THICKNESS: Key = "slab.thickness";
pub(in crate::modeling) const SLAB_OVERRIDE: Key = "slab.override";
pub(in crate::modeling) const PEDESTAL_HEIGHT: Key = "pedestal.height";
pub(in crate::modeling) const STAIR_WIDTH: Key = "stair.width";
pub(in crate::modeling) const STAIR_GOING: Key = "stair.going";
pub(in crate::modeling) const STAIR_RISERS: Key = "stair.risers";
pub(in crate::modeling) const STAIR_DIRECTION: Key = "stair.direction";
pub(in crate::modeling) const STAIR_START: Key = "stair.start";
/// Fields of the stair's Position section.
pub(in crate::modeling) const STAIR_POSITION: [Key; 2] = [STAIR_START, STAIR_DIRECTION];

#[derive(Clone, Copy, Debug, Default)]
pub(in crate::modeling) struct StairValues {
    pub width: Option<f64>,
    pub going: Option<f64>,
    pub risers: u32,
    /// Riser height from the resolved levels, when they resolve.
    pub riser: Option<f64>,
}

/// Plan rectangles (local frame, walking along +X, left is +Y) of each
/// flight and landing, plus the tread lines across them.
struct Layout {
    pieces: Vec<[[f32; 2]; 4]>,
    treads: Vec<[[f32; 2]; 2]>,
    walk: Vec<[f32; 2]>,
}

fn rect(x0: f32, y0: f32, x1: f32, y1: f32) -> [[f32; 2]; 4] {
    [[x0, y0], [x1, y0], [x1, y1], [x0, y1]]
}

fn layout(form: StairForm, risers: u32, g: f32, w: f32) -> Layout {
    let mut l = Layout {
        pieces: vec![],
        treads: vec![],
        walk: vec![[0., 0.]],
    };
    let across_x = |x: f32, yc: f32| [[x, yc - w / 2.], [x, yc + w / 2.]];
    match form {
        StairForm::Straight => {
            let run = risers as f32 * g;
            l.pieces.push(rect(0., -w / 2., run, w / 2.));
            l.treads
                .extend((1..risers).map(|k| across_x(k as f32 * g, 0.)));
            l.walk.push([run, 0.]);
        }
        StairForm::QuarterTurn { turn } | StairForm::HalfTurn { turn } => {
            let n1 = risers / 2;
            let n2 = risers - n1;
            let s = if turn == Turn::Left { 1. } else { -1. };
            let flight1 = n1.saturating_sub(1) as f32 * g;
            // Flight 1 ends in a landing one width deep.
            l.pieces.push(rect(0., -w / 2., flight1 + w, w / 2.));
            l.treads.extend((1..n1).map(|k| across_x(k as f32 * g, 0.)));
            l.treads.push(across_x(flight1, 0.));
            if matches!(form, StairForm::HalfTurn { .. }) {
                let yc = s * w;
                l.pieces
                    .push(rect(flight1, yc - w / 2., flight1 + w, yc + w / 2.));
                let end = flight1 - n2 as f32 * g;
                l.pieces.push(rect(end, yc - w / 2., flight1, yc + w / 2.));
                l.treads
                    .extend((0..n2).map(|k| across_x(flight1 - k as f32 * g, yc)));
                l.walk
                    .extend([[flight1 + w / 2., 0.], [flight1 + w / 2., yc], [end, yc]]);
            } else {
                let xc = flight1 + w / 2.;
                let y0 = s * w / 2.;
                let y1 = s * (w / 2. + n2 as f32 * g);
                l.pieces
                    .push(rect(xc - w / 2., y0.min(y1), xc + w / 2., y0.max(y1)));
                l.treads.extend((0..n2).map(|k| {
                    let y = y0 + s * k as f32 * g;
                    [[xc - w / 2., y], [xc + w / 2., y]]
                }));
                l.walk.extend([[xc, 0.], [xc, y1]]);
            }
        }
    }
    l
}

fn plan(mut d: Diagram, l: &Layout) -> Diagram {
    for piece in &l.pieces {
        d = d.solid(piece);
    }
    for tread in &l.treads {
        d = d.line(tread[0], tread[1], Line::Outline);
    }
    d.path(&l.walk, false, Line::Center)
}

/// Plan art for a stair form tile.
pub(in crate::modeling) fn form_art(form: StairForm) -> Diagram {
    plan(Diagram::new(), &layout(form, 10, 1., 3.))
}

pub(in crate::modeling) fn form_label(form: StairForm) -> &'static str {
    match form {
        StairForm::Straight => "Straight",
        StairForm::QuarterTurn { turn: Turn::Left } => "Quarter turn left",
        StairForm::QuarterTurn { turn: Turn::Right } => "Quarter turn right",
        StairForm::HalfTurn { turn: Turn::Left } => "Half turn left",
        StairForm::HalfTurn { turn: Turn::Right } => "Half turn right",
    }
}

pub(in crate::modeling) const FORMS: [StairForm; 5] = [
    StairForm::Straight,
    StairForm::QuarterTurn { turn: Turn::Left },
    StairForm::QuarterTurn { turn: Turn::Right },
    StairForm::HalfTurn { turn: Turn::Left },
    StairForm::HalfTurn { turn: Turn::Right },
];

/// The stair diagram: local plan (first riser at the origin, walking along
/// +X) above the walking-line elevation of every riser.
pub(in crate::modeling) fn stair(form: StairForm, v: &StairValues) -> Diagram {
    let risers = v.risers.max(2);
    let g = legible(v.going, 270., 1., 1e6);
    let w = legible(v.width, 1000., 0.3 * g, 12. * g);
    let riser = legible(v.riser, 0.65 * g, 0.05 * g, 3. * g);
    let l = layout(form, risers, g, w);
    let mut d = plan(Diagram::new(), &l);
    let min_y = l
        .pieces
        .iter()
        .flatten()
        .map(|p| p[1])
        .fold(f32::INFINITY, f32::min);
    // Walking-line elevation below the plan, unrolled over every riser.
    let base = min_y - 0.9 * w - risers as f32 * riser;
    let mut profile = vec![[0., base]];
    for k in 0..risers {
        let x = k as f32 * g;
        let z = base + (k + 1) as f32 * riser;
        profile.push([x, z]);
        profile.push([x + g, z]);
    }
    profile.push([risers as f32 * g, base]);
    d = d
        .path(&profile, true, Line::Outline)
        .line(
            [-0.4 * g, base],
            [(risers as f32 + 0.4) * g, base],
            Line::Center,
        )
        .dim(
            STAIR_WIDTH,
            [0., -w / 2.],
            [0., w / 2.],
            -16.,
            Notation::Length,
            "Width",
        )
        .dim(
            STAIR_GOING,
            [0., base + riser],
            [g, base + riser],
            14.,
            Notation::Length,
            "Going",
        )
        .leader(
            STAIR_RISERS,
            [risers as f32 * g * 0.6, base + risers as f32 * riser * 0.6],
            [
                risers as f32 * g * 0.95,
                base + risers as f32 * riser * 1.05,
            ],
            Notation::Length,
            "Risers",
        )
        .leader(
            STAIR_DIRECTION,
            [0.6 * g, 0.],
            [2.5 * g, w * 0.9],
            Notation::Angle,
            "Direction",
        )
        .leader(
            STAIR_START,
            [0., 0.],
            [-0.6 * g, -w * 0.95],
            Notation::Length,
            "First riser",
        )
        .note(
            [risers as f32 * g * 0.5, base - 0.35 * w],
            match v.riser {
                Some(r) => format!("{} · riser {r:.1} mm · walking line", form_label(form)),
                None => format!("{} · choose levels for the rise", form_label(form)),
            },
        );
    d
}

#[derive(Clone, Debug, Default)]
pub(in crate::modeling) struct SlabValues {
    /// Resolved outline corners in millimetres, when the outline resolves.
    pub outline: Option<Vec<[f64; 2]>>,
    pub pedestal: bool,
    pub override_thickness: bool,
    pub thickness: Option<f64>,
    /// The level's own slab thickness, when the level resolves.
    pub level_thickness: Option<f64>,
}

/// Floor slab or pedestal: plan footprint above a section at the level
/// floor. A slab hangs below the floor; a pedestal stands on it.
pub(in crate::modeling) fn slab(v: &SlabValues) -> Diagram {
    let points: Vec<[f32; 2]> = v
        .outline
        .as_ref()
        .filter(|p| p.len() >= 3)
        .map(|p| p.iter().map(|q| [q[0] as f32, q[1] as f32]).collect())
        .unwrap_or_default();
    let mut d = Diagram::new();
    let (x0, x1, y0) = if points.is_empty() {
        d = d
            .outline_rect(0., 0., 6000., 4000., Line::Hidden)
            .note([3000., 2000.], "Outline incomplete");
        (0., 6000., 0.)
    } else {
        let fold =
            |f: fn(f32, f32) -> f32, i: usize, init: f32| points.iter().map(|p| p[i]).fold(init, f);
        d = d.path(&points, true, Line::Outline);
        (
            fold(f32::min, 0, f32::INFINITY),
            fold(f32::max, 0, f32::NEG_INFINITY),
            fold(f32::min, 1, f32::INFINITY),
        )
    };
    let span = (x1 - x0).max(1.);
    let floor = y0 - 0.22 * span;
    let value = if v.pedestal || v.override_thickness {
        v.thickness
    } else {
        v.level_thickness
    };
    let t = legible(value, 0.04 * span, 0.025 * span, 0.2 * span);
    d = d.line(
        [x0 - 0.08 * span, floor],
        [x1 + 0.08 * span, floor],
        Line::Center,
    );
    let (bottom, top) = if v.pedestal {
        (floor, floor + t)
    } else {
        (floor - t, floor)
    };
    d = d.rect(x0, bottom, x1, top).note(
        [x0 + 0.5 * span, bottom - 0.07 * span],
        if v.pedestal {
            "Section · raised on the level floor"
        } else {
            "Section · below the level floor"
        },
    );
    let (key, name) = if v.pedestal {
        (PEDESTAL_HEIGHT, "Height")
    } else if v.override_thickness {
        (SLAB_THICKNESS, "Thickness")
    } else {
        (SLAB_OVERRIDE, "Level thickness")
    };
    d.dim(key, [x1, bottom], [x1, top], -14., Notation::Length, name)
}
