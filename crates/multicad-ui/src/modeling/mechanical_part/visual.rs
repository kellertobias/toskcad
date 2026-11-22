//! Family-specific hardware schematics (CAD-386). Each diagram draws the
//! selected family with the dimensions its form actually edits; every callout
//! key names one canonical field of the hardware form.
use super::super::diagram::{Diagram, Key, Line, Notation, legible};
use super::Kind;
use super::screw_shape::HeadKind;
use multicad_core::model::ScrewLengthConvention;
use multicad_precision::PortableMath as _;

pub(super) const WASHER: [Key; 3] = ["washer.bore", "washer.outside", "washer.thickness"];
pub(super) const BEARING: [Key; 3] = ["bearing.bore", "bearing.outside", "bearing.width"];
pub(super) const NUT_FLATS: Key = "nut.flats";
pub(super) const SCREW_LENGTH: Key = "screw.length";
pub(super) const SCREW_SHANK: Key = "screw.shank";
pub(super) const THREAD_DIAMETER: Key = "thread.diameter";
pub(super) const THREAD_LENGTH: Key = "thread.length";
pub(super) const HEAD_DIAMETER: Key = "head.diameter";
pub(super) const HEAD_FLATS: Key = "head.flats";
pub(super) const HEAD_HEIGHT: Key = "head.height";
pub(super) const HEAD_ANGLE: Key = "head.angle";
pub(super) const WOOD: [Key; 8] = [
    "wood.length",
    "wood.diameter",
    "wood.head_flats",
    "wood.head_height",
    "wood.pitch",
    "wood.core",
    "wood.angle",
    "wood.extent",
];

/// Current values in millimetres (angles in degrees); `None` where a staged
/// field does not resolve, in which case a legible default proportion is used.
#[derive(Clone, Copy, Debug, Default)]
pub(super) struct Values {
    pub a: Option<f64>,
    pub b: Option<f64>,
    pub c: Option<f64>,
    pub length: Option<f64>,
    pub shank: Option<f64>,
    pub thread_diameter: Option<f64>,
    pub thread_length: Option<f64>,
    pub head_size: Option<f64>,
    pub head_height: Option<f64>,
}

/// Washer and bearing: plan view (two circles) beside a sectioned elevation.
fn ring(keys: [Key; 3], v: &Values, names: [&str; 3], defaults: [f32; 3]) -> Diagram {
    let od = legible(v.b, defaults[1], 1., 1e6);
    let bore = legible(v.a, defaults[0], 0.05 * od, 0.85 * od);
    let t = legible(v.c, defaults[2], 0.08 * od, 0.9 * od);
    let (r, b) = (od / 2., bore / 2.);
    let x = od * 0.85 + t;
    Diagram::new()
        .circle([0., 0.], r, Line::Outline, true)
        .circle([0., 0.], b, Line::Outline, false)
        .line([-r * 1.15, 0.], [r * 1.15, 0.], Line::Center)
        .rect(x, -r, x + t, -b)
        .rect(x, b, x + t, r)
        .line([x, -b], [x, b], Line::Hidden)
        .line([x + t, -b], [x + t, b], Line::Hidden)
        .line([x - t * 0.4, 0.], [x + t * 1.4, 0.], Line::Center)
        .dim(keys[0], [-b, 0.], [b, 0.], 0., Notation::Diameter, names[0])
        .dim(
            keys[1],
            [-r, -r],
            [r, -r],
            16.,
            Notation::Diameter,
            names[1],
        )
        .dim(
            keys[2],
            [x, r],
            [x + t, r],
            -16.,
            Notation::Length,
            names[2],
        )
}

fn nut(v: &Values) -> Diagram {
    let af = legible(v.a, 10., 1., 1e6);
    let major = legible(v.thread_diameter, 6., 0.1 * af, 0.9 * af);
    let h = legible(v.thread_length, 5.2, 0.2 * af, 2. * af);
    let corner = af / 3f32.sqrt();
    let hex: Vec<[f32; 2]> = (0..6)
        .map(|i| {
            let a = i as f32 * std::f32::consts::FRAC_PI_3;
            [corner * a.portable_cos(), corner * a.portable_sin()]
        })
        .collect();
    let x = corner * 1.9;
    let r = major / 2.;
    Diagram::new()
        .solid(&hex)
        .circle([0., 0.], r * 0.85, Line::Outline, false)
        .circle([0., 0.], r, Line::Hidden, false)
        .rect(x, -corner, x + h, corner)
        .line([x, -r], [x + h, -r], Line::Hidden)
        .line([x, r], [x + h, r], Line::Hidden)
        .line([x - h * 0.3, 0.], [x + h * 1.3, 0.], Line::Center)
        .dim(
            NUT_FLATS,
            [0., -af / 2.],
            [0., af / 2.],
            -(af * 0.0) - 28.,
            Notation::Length,
            "Across flats",
        )
        .leader(
            THREAD_DIAMETER,
            [r * 0.7, r * 0.7],
            [corner * 0.9, corner * 1.25],
            Notation::Diameter,
            "Thread",
        )
        .dim(
            THREAD_LENGTH,
            [x, -corner],
            [x + h, -corner],
            16.,
            Notation::Length,
            "Height",
        )
}

pub(super) struct ScrewShape {
    pub head: HeadKind,
    pub convention: ScrewLengthConvention,
    /// Wood envelopes have a smooth stock shank and a minimum thread extent.
    pub wood: bool,
}

/// Elevation with the axis vertical: head above the neck datum (y = 0), shank
/// to the tip below. Length endpoints follow the selected convention and the
/// actual head form, so an overall countersunk length includes its cone.
fn screw(shape: &ScrewShape, v: &Values) -> Diagram {
    let d = legible(v.shank, 3., 0.1, 1e6);
    let hh = legible(v.head_height, d, 0.25 * d, 1.5 * d);
    let hd = legible(v.head_size, 1.8 * d, 1.1 * d, 3. * d);
    let typed = legible(v.length, 10., 0.2 * d, 1e6);
    let underhead = match shape.convention {
        ScrewLengthConvention::Underhead => typed,
        ScrewLengthConvention::Overall => typed - hh,
    }
    .clamp(1.5 * d, 12. * d);
    let tl = legible(v.thread_length, underhead, 0.15 * underhead, underhead);
    let major = if shape.wood {
        d
    } else {
        legible(v.thread_diameter, d, 0.6 * d, 1.2 * d)
    };
    let (r, m, hr) = (d / 2., major / 2., hd / 2.);
    let tip = -underhead;
    let thread_top = tip + tl;
    let mut diagram = match shape.head {
        HeadKind::Countersunk => Diagram::new().solid(&[[-r, 0.], [r, 0.], [hr, hh], [-hr, hh]]),
        HeadKind::Pan | HeadKind::Button => {
            let crown: Vec<[f32; 2]> = (0..=12)
                .map(|i| {
                    let a = i as f32 / 12. * std::f32::consts::PI;
                    let x = hr * a.portable_cos();
                    let y = hh
                        * (0.35
                            + 0.65
                                * a.portable_sin()
                                    .portable_powf(if shape.head == HeadKind::Pan {
                                        0.5
                                    } else {
                                        0.9
                                    }));
                    [x, y]
                })
                .collect();
            let mut points = vec![[hr, 0.]];
            points.extend(crown);
            points.push([-hr, 0.]);
            Diagram::new().solid(&points)
        }
        HeadKind::ExternalHex => Diagram::new()
            .rect(-hr, 0., hr, hh)
            .line([-hr / 2., 0.], [-hr / 2., hh], Line::Outline)
            .line([hr / 2., 0.], [hr / 2., hh], Line::Outline),
        HeadKind::SocketCap => Diagram::new().rect(-hr, 0., hr, hh).outline_rect(
            -hr * 0.45,
            hh * 0.45,
            hr * 0.45,
            hh,
            Line::Hidden,
        ),
    };
    diagram = diagram
        .rect(-r, thread_top.min(0.), r, 0.)
        .rect(-m, tip, m, thread_top)
        .line([0., tip - d * 0.4], [0., hh + d * 0.4], Line::Center);
    if !shape.wood {
        // Conventional thin lines at the thread root.
        let root = m * 0.82;
        diagram = diagram
            .line([-root, tip], [-root, thread_top], Line::Hidden)
            .line([root, tip], [root, thread_top], Line::Hidden);
    }
    let (length_top, length_name) = match shape.convention {
        ScrewLengthConvention::Underhead => (0., "Length under head"),
        ScrewLengthConvention::Overall => (hh, "Overall length"),
    };
    let (length_key, diameter_key, head_size_key, head_height_key, thread_key) = if shape.wood {
        (WOOD[0], WOOD[1], WOOD[2], WOOD[3], WOOD[7])
    } else {
        (
            SCREW_LENGTH,
            SCREW_SHANK,
            if shape.head == HeadKind::ExternalHex {
                HEAD_FLATS
            } else {
                HEAD_DIAMETER
            },
            HEAD_HEIGHT,
            THREAD_LENGTH,
        )
    };
    let x = hr.max(m) * 1.05;
    diagram = diagram
        .dim(
            length_key,
            [x, length_top],
            [x, tip],
            20.,
            Notation::Length,
            length_name,
        )
        .dim(
            thread_key,
            [-m, tip],
            [-m, thread_top],
            20.,
            Notation::Length,
            if shape.wood {
                "Min. thread extent"
            } else {
                "Thread"
            },
        )
        .dim(
            head_size_key,
            [-hr, hh],
            [hr, hh],
            14.,
            if head_size_key == HEAD_FLATS || shape.wood {
                Notation::Length
            } else {
                Notation::Diameter
            },
            if head_size_key == HEAD_FLATS || shape.wood {
                "Head flats"
            } else {
                "Head"
            },
        )
        .dim(
            head_height_key,
            [hr, 0.],
            [hr, hh],
            -14.,
            Notation::Length,
            "Head height",
        );
    if shape.wood {
        diagram = diagram.dim(
            diameter_key,
            [-r, tip * 0.15],
            [r, tip * 0.15],
            0.,
            Notation::Diameter,
            "Stock",
        );
    } else {
        let y = thread_top.min(0.) / 2.;
        if thread_top < -0.05 * underhead {
            diagram = diagram.dim(
                diameter_key,
                [-r, y],
                [r, y],
                0.,
                Notation::Diameter,
                "Shank",
            );
        } else {
            // Fully threaded: the shank diameter is still the plain neck size.
            diagram = diagram.leader(
                diameter_key,
                [r, -0.02 * underhead],
                [x + d, -0.3 * underhead],
                Notation::Diameter,
                "Shank",
            );
        }
        diagram = diagram.dim(
            THREAD_DIAMETER,
            [-m, tip],
            [m, tip],
            -14.,
            Notation::Diameter,
            "Thread",
        );
        if shape.head == HeadKind::Countersunk {
            diagram = diagram.leader(
                HEAD_ANGLE,
                [hr * 0.75, hh * 0.75],
                [-hr * 1.6, hh * 1.4],
                Notation::Angle,
                "Head angle",
            );
        }
    }
    diagram
}

/// The dimension diagram of the selected family with the current values.
pub(super) fn family(kind: Kind, shape: &ScrewShape, v: &Values) -> Diagram {
    match kind {
        Kind::Washer => ring(WASHER, v, ["Bore", "Outside", "Thickness"], [3.2, 7., 0.5]),
        Kind::BearingEnvelope => ring(BEARING, v, ["Bore", "Outside", "Width"], [10., 30., 9.]),
        Kind::HexNut => nut(v),
        Kind::MachineScrew | Kind::WoodScrew => screw(shape, v),
    }
}

/// Selector artwork: the family drawn without dimensions at typical size.
pub(super) fn art(kind: Kind) -> Diagram {
    let shape = ScrewShape {
        head: if kind == Kind::WoodScrew {
            HeadKind::ExternalHex
        } else {
            HeadKind::SocketCap
        },
        convention: ScrewLengthConvention::Underhead,
        wood: kind == Kind::WoodScrew,
    };
    let values = match kind {
        Kind::Washer => Values {
            c: Some(1.4),
            ..Default::default()
        },
        _ => Values::default(),
    };
    let mut diagram = family(kind, &shape, &values);
    diagram.callouts.clear();
    diagram
}
