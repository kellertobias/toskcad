//! Surface relief schematic (CAD-395): one picked face with its flat border
//! and the rotated, offset relief window filled with the chosen pattern,
//! above a section of the pattern showing pitch and relief height. Each
//! picked face maps its own window; the face outline is representative.
use super::super::diagram::{Diagram, Key, Line, Notation, legible};
use multicad_precision::PortableMath as _;

pub(in crate::modeling) const PITCH: Key = "surface.pitch";
pub(in crate::modeling) const AMPLITUDE: Key = "surface.amplitude";
pub(in crate::modeling) const SIZE: [Key; 2] = ["surface.width", "surface.height"];
pub(in crate::modeling) const OFFSET: [Key; 2] = ["surface.offset_u", "surface.offset_v"];
pub(in crate::modeling) const ROTATION: Key = "surface.rotation";
pub(in crate::modeling) const BORDER: Key = "surface.border";
/// Fields of the named placement section.
pub(in crate::modeling) const PLACEMENT: [Key; 4] = [OFFSET[0], OFFSET[1], ROTATION, BORDER];

/// The pattern families the tiles offer; a height map is imported instead.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(in crate::modeling) enum Kind {
    Ribs,
    Waves,
    Stipple,
    HeightMap,
}

#[derive(Clone, Copy, Debug, Default)]
pub(in crate::modeling) struct Values {
    pub pitch: Option<f64>,
    pub amplitude: Option<f64>,
    pub size: [Option<f64>; 2],
    pub offset: [Option<f64>; 2],
    pub rotation_deg: Option<f64>,
    pub border: Option<f64>,
    pub recessed: bool,
}

fn rotate(p: [f32; 2], c: [f32; 2], (s, k): (f32, f32)) -> [f32; 2] {
    let (dx, dy) = (p[0] - c[0], p[1] - c[1]);
    [c[0] + k * dx - s * dy, c[1] + s * dx + k * dy]
}

/// One period-repeating section of the pattern, `periods` long.
fn profile(kind: Kind, p: f32, a: f32, periods: usize) -> Vec<[f32; 2]> {
    let steps = 16;
    let mut points = vec![];
    for i in 0..=periods * steps {
        let x = i as f32 / steps as f32 * p;
        let phase = (i % steps) as f32 / steps as f32;
        let h = match kind {
            Kind::Ribs => {
                if phase < 0.5 {
                    a
                } else {
                    0.
                }
            }
            Kind::Waves => a * 0.5 * (1. - (phase * std::f32::consts::TAU).portable_cos()),
            Kind::Stipple => a * (phase * std::f32::consts::PI).portable_sin().max(0.),
            // A representative grey-level profile.
            Kind::HeightMap => a * [0.2, 0.9, 0.5, 1., 0.1, 0.7][i / 3 % 6],
        };
        points.push([x, h]);
    }
    points
}

pub(in crate::modeling) fn tile_art(kind: Kind) -> Diagram {
    let points = profile(kind, 4., 1.6, 3);
    let mut d = Diagram::new().line([0., 0.], [12., 0.], Line::Center);
    d = d.path(&points, false, Line::Outline);
    if kind == Kind::Stipple {
        for i in 0..3 {
            d = d.circle([2. + 4. * i as f32, -1.8], 0.6, Line::Outline, true);
        }
    } else {
        d = d.line([0., -1.8], [12., -1.8], Line::Hidden);
    }
    d
}

pub(in crate::modeling) fn diagram(kind: Kind, v: &Values, image: Option<[u32; 2]>) -> Diagram {
    let u = legible(v.size[0], 40., 1e-3, 1e7);
    let vv = legible(v.size[1], 30., 0.2 * u, 3. * u);
    let b = legible(v.border, 0.06 * u, 0.04 * u, 0.2 * u);
    let p = legible(v.pitch, u / 8., u / 16., u / 2.);
    let a = legible(v.amplitude, 0.3 * p, 0.12 * p, 0.8 * p) * if v.recessed { -1. } else { 1. };
    let off = [
        v.offset[0]
            .map_or(0., |o| o as f32)
            .clamp(-0.15 * u, 0.15 * u),
        v.offset[1]
            .map_or(0., |o| o as f32)
            .clamp(-0.15 * vv, 0.15 * vv),
    ];
    let rot = v.rotation_deg.unwrap_or(0.) as f32;
    let sc = rot.to_radians().portable_sin_cos();
    // Representative face, centred at the origin, with its flat border.
    let (fw, fh) = (u * 0.7 + b + 0.12 * u, vv * 0.7 + b + 0.12 * u);
    let face = [[-fw, -fh], [fw, -fh], [fw, fh], [-fw, fh]];
    let inner = [
        [-fw + b, -fh + b],
        [fw - b, -fh + b],
        [fw - b, fh - b],
        [-fw + b, fh - b],
    ];
    let c = off;
    let corner = |x: f32, y: f32| rotate([c[0] + x, c[1] + y], c, sc);
    let window = [
        corner(-u / 2., -vv / 2.),
        corner(u / 2., -vv / 2.),
        corner(u / 2., vv / 2.),
        corner(-u / 2., vv / 2.),
    ];
    let mut d = Diagram::new()
        .path(&face, true, Line::Outline)
        .path(&inner, true, Line::Hidden)
        .path(&window, true, Line::Outline);
    // Pattern hint inside the window, at the schematic pitch.
    let count = ((u / p).floor() as usize).clamp(2, 14);
    let step = u / count as f32;
    for i in 1..count {
        let x = -u / 2. + i as f32 * step;
        match kind {
            Kind::Ribs => {
                d = d.line(corner(x, -vv / 2.), corner(x, vv / 2.), Line::Outline);
            }
            Kind::Waves => {
                let points: Vec<_> = (0..=12)
                    .map(|j| {
                        let y = -vv / 2. + j as f32 / 12. * vv;
                        let wobble = 0.25 * step * (j as f32 * 1.1).portable_sin();
                        corner(x + wobble, y)
                    })
                    .collect();
                d = d.path(&points, false, Line::Outline);
            }
            Kind::Stipple => {
                let rows = ((vv / step).floor() as usize).clamp(2, 14);
                for j in 1..rows {
                    let y = -vv / 2. + j as f32 * vv / rows as f32;
                    d = d.circle(corner(x, y), 0.12 * step, Line::Outline, true);
                }
            }
            Kind::HeightMap => {
                d = d.line(corner(x, -vv / 2.), corner(x, vv / 2.), Line::Hidden);
            }
        }
    }
    d = d
        .dim(
            SIZE[0],
            window[0],
            window[1],
            -14.,
            Notation::Length,
            "Width U",
        )
        .dim(
            SIZE[1],
            window[1],
            window[2],
            -14.,
            Notation::Length,
            "Height V",
        )
        .leader(
            OFFSET[0],
            c,
            [-0.55 * fw, fh + 0.22 * u],
            Notation::Length,
            "Offset",
        )
        .leader(
            ROTATION,
            corner(u / 2., vv / 4.),
            [0.75 * fw, fh + 0.22 * u],
            Notation::Angle,
            "Rotation",
        )
        .dim(
            BORDER,
            [-fw, -fh + 0.5 * b],
            [-fw + b, -fh + 0.5 * b],
            -14.,
            Notation::Length,
            "Border",
        );
    // Section of the pattern below the face.
    let periods = 3;
    let y0 = -fh - 0.42 * u - a.abs();
    let x0 = -1.5 * p;
    let section: Vec<_> = profile(kind, p, a, periods)
        .into_iter()
        .map(|q| [x0 + q[0], y0 + q[1]])
        .collect();
    d = d
        .line(
            [x0 - 0.3 * p, y0],
            [x0 + (periods as f32 + 0.3) * p, y0],
            Line::Center,
        )
        .path(&section, false, Line::Outline)
        .dim(
            AMPLITUDE,
            [x0 + 3.2 * p, y0],
            [x0 + 3.2 * p, y0 + a],
            -14.,
            Notation::Length,
            if v.recessed { "Depth" } else { "Height" },
        );
    match image {
        Some([w, h]) => d.note(
            [x0 + 1.5 * p, y0 - 0.35 * u],
            format!("Height map {w} × {h} px · black is zero"),
        ),
        None => d.dim(
            PITCH,
            [x0, y0 + a],
            [x0 + p, y0 + a],
            if v.recessed { -14. } else { 14. },
            Notation::Length,
            "Pitch",
        ),
    }
}
