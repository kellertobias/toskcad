//! Thread schematics (CAD-394): the axial profile of the selected family on an
//! external shaft or inside a receiving bore. The tooth count is schematic;
//! dimension labels carry the actual staged values.
use super::super::diagram::{Diagram, Key, Line, Notation, legible};
use multicad_core::model::{ThreadFamily, ThreadHand, ThreadMode};
use multicad_precision::PortableMath as _;

pub(in crate::modeling) const DIAMETER: Key = "thread.diameter";
pub(in crate::modeling) const PITCH: Key = "thread.pitch";
pub(in crate::modeling) const LENGTH: Key = "thread.length";
pub(in crate::modeling) const CUSTOM: [Key; 5] = [
    "thread.custom.depth",
    "thread.custom.crest",
    "thread.custom.lower",
    "thread.custom.upper",
    "thread.custom.taper",
];
pub(in crate::modeling) const START_DEPTH: Key = "thread.start_depth";
pub(in crate::modeling) const RELIEF_LENGTH: Key = "thread.relief_length";
pub(in crate::modeling) const RELIEF_DEPTH: Key = "thread.relief_depth";
pub(in crate::modeling) const ALLOWANCE: [Key; 3] =
    ["thread.radial", "thread.axial", "thread.wall"];
pub(in crate::modeling) const ORIGIN: [Key; 3] = ["thread.x", "thread.y", "thread.z"];
pub(in crate::modeling) const ENGAGEMENT: Key = "thread.engagement";

/// Tooth form in pitch units: flank angles from the radial direction, depth
/// and crest/root flats. Taper is on radius per unit length.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(in crate::modeling) struct Form {
    pub lower: f32,
    pub upper: f32,
    pub depth: f32,
    pub rounded: bool,
    pub taper: f32,
    pub included: &'static str,
}

pub(in crate::modeling) fn form(
    family: ThreadFamily,
    custom: Option<[f64; 5]>,
    pitch: f32,
) -> Form {
    let v = |half: f32, depth: f32, rounded, taper, included| Form {
        lower: half,
        upper: half,
        depth,
        rounded,
        taper,
        included,
    };
    match family {
        ThreadFamily::IsoMetric | ThreadFamily::Unified | ThreadFamily::Unc | ThreadFamily::Unf => {
            v(30., 0.54, false, 0., "60° V")
        }
        ThreadFamily::Whitworth | ThreadFamily::Bsw | ThreadFamily::Bsf | ThreadFamily::Bspp => {
            v(27.5, 0.64, true, 0., "55° rounded V")
        }
        ThreadFamily::Bspt => v(27.5, 0.64, true, 1. / 32., "55° rounded V · taper 1:16"),
        ThreadFamily::Npt => v(30., 0.8, false, 1. / 32., "60° V · taper 1:16"),
        ThreadFamily::Trapezoidal => v(15., 0.5, false, 0., "30° trapezoidal"),
        ThreadFamily::Acme => v(14.5, 0.5, false, 0., "29° Acme"),
        ThreadFamily::Buttress => Form {
            lower: 45.,
            upper: 7.,
            depth: 0.6,
            rounded: false,
            taper: 0.,
            included: "Buttress 7°/45°",
        },
        ThreadFamily::Round => v(15., 0.5, true, 0., "30° round"),
        ThreadFamily::Custom => {
            let [depth, _crest, lower, upper, taper] = custom.unwrap_or([0.5, 0.2, 30., 30., 0.]);
            Form {
                lower: (lower as f32).clamp(0., 60.),
                upper: (upper as f32).clamp(0., 60.),
                depth: (depth as f32 / pitch.max(1e-3)).clamp(0.1, 0.9),
                rounded: false,
                taper: (taper as f32).clamp(-0.2, 0.2),
                included: "Custom flanks",
            }
        }
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub(in crate::modeling) struct Values {
    pub diameter: Option<f64>,
    pub pitch: Option<f64>,
    pub length: Option<f64>,
    pub custom: Option<[f64; 5]>,
    pub start_chamfer: bool,
    pub end_relief: bool,
}

pub(in crate::modeling) struct Variant {
    pub family: ThreadFamily,
    pub mode: ThreadMode,
    pub hand: ThreadHand,
    pub starts: u32,
}

/// One tooth row along +x at the given base radius; `sign` flips the
/// profile inward for an internal thread.
fn silhouette(form: &Form, pitch: f32, length: f32, root: f32, sign: f32) -> Vec<[f32; 2]> {
    let depth = form.depth * pitch;
    let rise = |angle: f32| depth * angle.to_radians().portable_tan();
    let (a, b) = (rise(form.lower), rise(form.upper));
    let flats = (pitch - a - b).max(0.1 * pitch);
    let (crest, bottom) = (flats * 0.4, flats * 0.6);
    let mut points = vec![];
    let mut x = 0.;
    let radius = |x: f32, h: f32| root + sign * (h + form.taper * 3. * x);
    while x < length - 1e-4 {
        let steps = [
            (x, 0.),
            (x + bottom * 0.5, 0.),
            (x + bottom * 0.5 + a, depth),
            (x + bottom * 0.5 + a + crest, depth),
            (x + bottom * 0.5 + a + crest + b, 0.),
        ];
        for (px, h) in steps {
            if px <= length {
                let h = if form.rounded && h > 0. { h * 0.92 } else { h };
                points.push([px, radius(px, h)]);
            }
        }
        x += pitch;
    }
    points.push([length, radius(length, 0.)]);
    points
}

/// The thread diagram for the selected family and mode.
pub(in crate::modeling) fn diagram(variant: &Variant, v: &Values) -> Diagram {
    let d = legible(v.diameter, 6., 0.1, 1e6);
    let true_pitch = legible(v.pitch, 1., 0.01 * d, d);
    let l = legible(v.length, 2. * d, 0.6 * d, 4. * d);
    let teeth = (l / true_pitch).round().clamp(3., 8.);
    let p = l / teeth;
    let form = form(variant.family, v.custom, true_pitch);
    let depth = (form.depth * p).min(0.22 * d);
    let form = Form {
        depth: depth / p,
        ..form
    };
    let r = d / 2.;
    let internal = variant.mode == ThreadMode::CutInternal;
    let mut diagram = Diagram::new();
    let axis = [[-0.15 * l, 0.], [1.15 * l, 0.]];
    if internal {
        // Receiving material above and below a threaded bore.
        let outer = r * 1.9;
        let tops = silhouette(&form, p, l, r, -1.);
        let mut bottom: Vec<[f32; 2]> = tops.iter().map(|q| [q[0], -q[1]]).collect();
        diagram = diagram
            .rect(-0.1 * l, r, l * 1.1, outer)
            .rect(-0.1 * l, -outer, l * 1.1, -r)
            .path(&tops, false, Line::Outline);
        bottom.reverse();
        diagram = diagram.path(&bottom, false, Line::Outline);
    } else {
        let core = r - depth;
        let top = silhouette(&form, p, l, core, 1.);
        let mut outline = top.clone();
        outline.extend(top.iter().rev().map(|q| [q[0], -q[1]]));
        diagram = diagram
            .rect(0., -core, l, core)
            .path(&outline, true, Line::Outline)
            .rect(-0.25 * l, -core * 0.98, 0., core * 0.98);
        if v.start_chamfer {
            diagram = diagram.line([0., core], [depth, r], Line::Outline).leader(
                START_DEPTH,
                [depth * 0.5, core + depth * 0.5],
                [-0.12 * l, r * 1.45],
                Notation::Length,
                "Entry chamfer",
            );
        }
        if v.end_relief {
            let x0 = l * 0.86;
            diagram = diagram
                .outline_rect(x0, -core * 0.9, l, core * 0.9, Line::Hidden)
                .dim(
                    RELIEF_LENGTH,
                    [x0, -r],
                    [l, -r],
                    -30.,
                    Notation::Length,
                    "Relief",
                )
                .leader(
                    RELIEF_DEPTH,
                    [l * 0.93, core * 0.9],
                    [l * 1.1, r * 1.5],
                    Notation::Length,
                    "Relief depth",
                );
        }
    }
    let crest = if internal { r - depth } else { r };
    diagram = diagram
        .line(axis[0], axis[1], Line::Center)
        .dim(
            DIAMETER,
            [l * 0.5, -r],
            [l * 0.5, r],
            -(0.0) + 0.,
            Notation::Diameter,
            "Major",
        )
        .dim(
            PITCH,
            [0.5 * p, crest],
            [1.5 * p, crest],
            14.,
            Notation::Length,
            "Pitch",
        )
        .dim(
            LENGTH,
            [0., -r * if internal { 1.9 } else { 1. }],
            [l, -r * if internal { 1.9 } else { 1. }],
            -16.,
            Notation::Length,
            "Length",
        )
        .note(
            [l * 0.5, r * if internal { 2.3 } else { 1.65 }],
            format!(
                "{} · {} · {} · {} start{}",
                variant.family.label(),
                form.included,
                match variant.hand {
                    ThreadHand::Right => "right hand",
                    ThreadHand::Left => "left hand",
                },
                variant.starts,
                if variant.starts == 1 { "" } else { "s" }
            ),
        );
    if variant.family == ThreadFamily::Custom {
        let x = 2.5 * p;
        diagram = diagram
            .leader(
                CUSTOM[0],
                [x, crest - depth * 0.5],
                [x + 0.6 * l, r * 1.35],
                Notation::Length,
                "Depth",
            )
            .leader(
                CUSTOM[2],
                [x - 0.15 * p, crest - depth * 0.3],
                [x - 0.9 * p, r * 1.35],
                Notation::Angle,
                "Lower flank",
            )
            .leader(
                CUSTOM[3],
                [x + 0.25 * p, crest - depth * 0.3],
                [x + 1.6 * p, r * 1.35],
                Notation::Angle,
                "Upper flank",
            );
    }
    diagram
}
