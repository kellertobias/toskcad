//! Gear schematics (CAD-394): each family is drawn as itself — external
//! cylindrical gears in front and side view with their tooth traces, the
//! internal ring with inward teeth, the rack as a straight toothed bar, the
//! bevel as its pitch cone and the worm or worm wheel. Callout keys name the
//! gear form's canonical fields.
use super::super::diagram::{Diagram, Key, Line, Notation, legible};
use multicad_core::model::{GearFamily, GearHand, GearPitchInput, WormMember};
use std::f32::consts::{PI, TAU};

pub(super) const FIELDS: [Key; 11] = [
    "gear.pitch",
    "gear.teeth",
    "gear.pressure",
    "gear.face",
    "gear.helix",
    "gear.gap",
    "gear.bore",
    "gear.backlash",
    "gear.clearance",
    "gear.fillet",
    "gear.wall",
];
pub(super) const OUTSIDE: Key = "gear.outside";
pub(super) const BACKING: Key = "gear.backing";
pub(super) const SHAFT: Key = "gear.shaft";
pub(super) const BEVEL_MATE: Key = "gear.bevel_mate";
pub(super) const WORM_MATE: Key = "gear.worm_mate";
pub(super) const QUOTIENT: Key = "gear.quotient";
pub(super) const ORIGIN: [Key; 3] = ["gear.x", "gear.y", "gear.z"];
/// Fields grouped in the tolerance section.
pub(super) const TOLERANCES: [Key; 4] = [FIELDS[7], FIELDS[8], FIELDS[9], FIELDS[10]];

#[derive(Clone, Copy, Debug)]
pub(super) struct Variant {
    pub family: GearFamily,
    pub hand: GearHand,
    pub pitch_input: GearPitchInput,
    pub worm_role: WormMember,
}

#[derive(Clone, Copy, Debug, Default)]
pub(super) struct Values {
    /// Module in millimetres (diametral pitch already converted).
    pub module: Option<f64>,
    pub teeth: Option<f64>,
    pub face: Option<f64>,
    pub helix: Option<f64>,
    pub gap: Option<f64>,
    pub bore: Option<f64>,
    pub outside: Option<f64>,
    pub backing: Option<f64>,
    pub mate: Option<f64>,
}

/// Closed tooth outline around the axis; `inward` points teeth to the center.
fn tooth_ring(z: usize, pitch_radius: f32, module: f32, inward: bool) -> Vec<[f32; 2]> {
    let (tip, root) = if inward {
        (pitch_radius - module, pitch_radius + 1.25 * module)
    } else {
        (pitch_radius + module, pitch_radius - 1.25 * module)
    };
    let step = TAU / z as f32;
    let mut points = Vec::with_capacity(z * 4);
    for i in 0..z {
        let c = i as f32 * step + PI / 2.;
        for (r, a) in [
            (root, c - 0.29 * step),
            (tip, c - 0.14 * step),
            (tip, c + 0.14 * step),
            (root, c + 0.29 * step),
        ] {
            points.push([r * a.cos(), r * a.sin()]);
        }
    }
    points.push(points[0]);
    points
}

fn pitch_label(variant: &Variant) -> &'static str {
    match (variant.pitch_input, variant.family) {
        (GearPitchInput::DiametralPitch, _) => "Diametral pitch",
        (_, GearFamily::Worm) => "Axial module",
        (_, GearFamily::Bevel) => "Outer module",
        _ => "Module",
    }
}

/// Front view (tooth ring, pitch circle, bore) beside a side view whose
/// tooth traces show the family: straight, helical, chevron or split bands.
fn cylindrical(variant: &Variant, v: &Values) -> Diagram {
    let z = legible(v.teeth, 24., 6., 400.).round();
    let m = legible(v.module, 2., 0.01, 1e5);
    let r = m * z / 2.;
    let drawn = (z as usize).min(96);
    let dm = 2. * r / drawn as f32;
    let tip = r + m;
    let internal = variant.family == GearFamily::InternalRing;
    let mut d = Diagram::new();
    if internal {
        let outer = legible(v.outside, 2. * (r + 4. * m), 2. * (r + 1.6 * m), 4. * r) / 2.;
        d = d
            .circle([0., 0.], outer, Line::Outline, true)
            .path(&tooth_ring(drawn, r, dm, true), true, Line::Outline)
            .dim(
                OUTSIDE,
                [-outer, -outer],
                [outer, -outer],
                16.,
                Notation::Diameter,
                "Outside",
            );
    } else {
        let ring = tooth_ring(drawn, r, dm, false);
        d = d
            .circle([0., 0.], r - 1.25 * dm, Line::Outline, true)
            .path(&ring, true, Line::Outline);
        let bore = legible(v.bore, 0., 0., 1.6 * (r - 1.25 * dm));
        if bore > 0. {
            d = d.circle([0., 0.], bore / 2., Line::Outline, false).dim(
                FIELDS[6],
                [-bore / 2., 0.],
                [bore / 2., 0.],
                0.,
                Notation::Diameter,
                "Bore",
            );
        }
    }
    d = d
        .circle([0., 0.], r, Line::Center, false)
        .leader(
            FIELDS[1],
            [0., if internal { r - dm } else { tip }],
            [-0.75 * tip, 1.25 * tip],
            Notation::Length,
            "Teeth",
        )
        .leader(
            FIELDS[0],
            [0.12 * r, 0.99 * r],
            [0.85 * tip, 1.25 * tip],
            Notation::Length,
            pitch_label(variant),
        )
        .leader(
            FIELDS[2],
            [-0.35 * r, 0.93 * r],
            [-1.05 * tip, 0.55 * tip],
            Notation::Angle,
            "Pressure",
        );
    // Side view: axis horizontal, face width along x.
    let b = legible(v.face, 0.5 * r, 0.15 * tip, 1.5 * tip);
    let x0 = tip * 1.45;
    let h = if internal {
        legible(v.outside, 2. * (r + 4. * m), 2. * tip, 4. * r) / 2.
    } else {
        tip
    };
    d = d
        .rect(x0, -h, x0 + b, h)
        .line([x0 - 0.15 * b, 0.], [x0 + 1.15 * b, 0.], Line::Center)
        .dim(
            FIELDS[3],
            [x0, -h],
            [x0 + b, -h],
            -16.,
            Notation::Length,
            "Face width",
        );
    let beta = legible(v.helix, 20., 5., 45.).to_radians();
    let sign = if variant.hand == GearHand::Right {
        1.
    } else {
        -1.
    };
    let traces = 5;
    let band = |d: Diagram, xa: f32, xb: f32, s: f32| {
        (0..traces).fold(d, |d, i| {
            let y = -0.8 * h + 1.6 * h * i as f32 / (traces - 1) as f32;
            let dy = s * (xb - xa) * beta.tan();
            d.line([xa, y - dy / 2.], [xb, y + dy / 2.], Line::Hidden)
        })
    };
    d = match variant.family {
        GearFamily::Helical => {
            let d = band(d, x0, x0 + b, sign);
            d.leader(
                FIELDS[4],
                [x0 + 0.5 * b, 0.4 * h],
                [x0 + 1.2 * b, 1.25 * tip],
                Notation::Angle,
                "Helix",
            )
        }
        GearFamily::Herringbone => {
            let d = band(d, x0, x0 + b / 2., sign);
            let d = band(d, x0 + b / 2., x0 + b, -sign);
            d.leader(
                FIELDS[4],
                [x0 + 0.25 * b, 0.4 * h],
                [x0 + 1.2 * b, 1.25 * tip],
                Notation::Angle,
                "Helix",
            )
        }
        GearFamily::DoubleHelical => {
            let gap = legible(v.gap, 0.2 * b, 0.05 * b, 0.6 * b);
            let band_w = (b - gap) / 2.;
            let d = band(d, x0, x0 + band_w, sign);
            let d = band(d, x0 + band_w + gap, x0 + b, -sign);
            d.outline_rect(
                x0 + band_w,
                -0.8 * h,
                x0 + band_w + gap,
                0.8 * h,
                Line::Hidden,
            )
            .dim(
                FIELDS[5],
                [x0 + band_w, h],
                [x0 + band_w + gap, h],
                14.,
                Notation::Length,
                "Center gap",
            )
            .leader(
                FIELDS[4],
                [x0 + 0.25 * b, 0.4 * h],
                [x0 + 1.2 * b, 1.4 * tip],
                Notation::Angle,
                "Helix",
            )
        }
        _ => band(d, x0, x0 + b, 0.),
    };
    d
}

fn rack(variant: &Variant, v: &Values) -> Diagram {
    let cells = legible(v.teeth, 12., 4., 64.).round();
    let m = legible(v.module, 2., 0.01, 1e5);
    let drawn = cells.min(12.);
    let p = PI * m * cells / drawn;
    let dm = p / PI;
    let length = p * drawn;
    let backing = legible(v.backing, 4., 0.3 * dm, 6. * dm);
    let root = -1.25 * dm;
    let mut teeth = vec![[0., root]];
    for i in 0..drawn as usize {
        let x = i as f32 * p;
        teeth.extend([
            [x + 0.1 * p, root],
            [x + 0.3 * p, dm],
            [x + 0.5 * p, dm],
            [x + 0.7 * p, root],
        ]);
    }
    teeth.push([length, root]);
    Diagram::new()
        .rect(0., root - backing, length, root)
        .path(&teeth, false, Line::Outline)
        .line([-0.05 * length, 0.], [1.05 * length, 0.], Line::Center)
        .dim(
            FIELDS[0],
            [0.4 * p, 0.],
            [1.4 * p, 0.],
            0.,
            Notation::Length,
            format!("{} (pitch = π·m)", pitch_label(variant)),
        )
        .dim(
            BACKING,
            [length, root - backing],
            [length, root],
            16.,
            Notation::Length,
            "Backing",
        )
        .dim(
            FIELDS[1],
            [0., root - backing],
            [length, root - backing],
            -16.,
            Notation::Length,
            "Pitch cells",
        )
        .leader(
            FIELDS[2],
            [2.62 * p, 0.6 * dm],
            [3.6 * p, 2.6 * dm],
            Notation::Angle,
            "Pressure",
        )
        .leader(
            FIELDS[3],
            [0.6 * length, root - 0.5 * backing],
            [0.8 * length, -3.5 * dm - backing],
            Notation::Length,
            "Face width (along Z)",
        )
}

/// Axial section through the bevel: pitch cone from the apex, the tooth band
/// between heel and toe along the generator, and the shaft axis.
fn bevel(variant: &Variant, v: &Values) -> Diagram {
    let z = legible(v.teeth, 30., 6., 400.);
    let mate = legible(v.mate, z, 6., 400.);
    let m = legible(v.module, 2., 0.01, 1e5);
    let delta = (z / mate).atan();
    let cone = m * z / 2. / delta.sin();
    let b = legible(v.face, cone / 4., 0.1 * cone, 0.45 * cone);
    let dir = [delta.sin(), delta.cos()];
    let at = |s: f32, off: f32| [dir[0] * s + dir[1] * off, dir[1] * s - dir[0] * off];
    let (heel, toe) = (cone, cone - b);
    // Half tooth height either side of the pitch cone (addendum ≈ module).
    let t = (1.1 * m).max(0.04 * cone);
    let band = [at(toe, -t), at(heel, -t), at(heel, t), at(toe, t)];
    let mirror = band.map(|p| [-p[0], p[1]]);
    Diagram::new()
        .solid(&band)
        .solid(&mirror)
        .line([0., 0.], at(heel * 1.08, 0.), Line::Center)
        .line(
            [0., 0.],
            [-dir[0] * heel * 1.08, dir[1] * heel * 1.08],
            Line::Center,
        )
        .line([0., -0.1 * cone], [0., cone * dir[1] * 1.15], Line::Center)
        .dim(
            FIELDS[3],
            at(toe, t),
            at(heel, t),
            16.,
            Notation::Length,
            "Face width",
        )
        .leader(
            FIELDS[0],
            at(heel, 0.),
            [dir[0] * heel * 1.25, dir[1] * heel * 0.75],
            Notation::Length,
            pitch_label(variant),
        )
        .leader(
            FIELDS[1],
            at(heel - b / 2., -t),
            [dir[0] * heel * 1.25, dir[1] * heel * 0.35],
            Notation::Length,
            "Teeth",
        )
        .leader(
            BEVEL_MATE,
            [-dir[0] * heel * 0.7, dir[1] * heel * 0.7],
            [-dir[0] * heel * 1.1, dir[1] * heel * 1.15],
            Notation::Length,
            "Mate teeth",
        )
        .leader(
            SHAFT,
            [0., 0.],
            [-dir[0] * heel * 0.9, -0.08 * cone],
            Notation::Angle,
            "Shaft angle",
        )
        .note([0., -0.16 * cone], "Pitch apex")
}

fn worm(variant: &Variant, v: &Values) -> Diagram {
    let m = legible(v.module, 1., 0.01, 1e5);
    if variant.worm_role == WormMember::Wheel {
        let mut wheel = cylindrical(
            &Variant {
                family: GearFamily::Spur,
                ..*variant
            },
            v,
        );
        wheel.callouts.retain(|c| c.key != FIELDS[6]);
        return wheel.leader(
            WORM_MATE,
            [0., 0.],
            [0., -0.7 * m * legible(v.teeth, 32., 6., 400.)],
            Notation::Length,
            "Worm starts",
        );
    }
    let d = 10. * m;
    let l = legible(v.face, 24. * m, 8. * m, 40. * m);
    let p = PI * m;
    let mut diagram = Diagram::new().rect(0., -d / 2. - m, l, d / 2. + m).line(
        [-0.1 * l, 0.],
        [1.1 * l, 0.],
        Line::Center,
    );
    let mut x = 0.;
    while x + 0.6 * p < l {
        diagram = diagram.line([x, -d / 2. - m], [x + 0.6 * p, d / 2. + m], Line::Hidden);
        x += p;
    }
    diagram
        .dim(
            FIELDS[3],
            [0., -d / 2. - m],
            [l, -d / 2. - m],
            -16.,
            Notation::Length,
            "Threaded length",
        )
        .dim(
            FIELDS[0],
            [0., d / 2. + m],
            [p, d / 2. + m],
            14.,
            Notation::Length,
            format!("{} (axial pitch = π·m)", pitch_label(variant)),
        )
        .leader(
            FIELDS[1],
            [0.6 * l, d / 4.],
            [0.75 * l, d],
            Notation::Length,
            "Starts",
        )
        .leader(
            QUOTIENT,
            [0.3 * l, -d / 4.],
            [0.2 * l, -d],
            Notation::Length,
            "Diameter quotient (pitch Ø = q·m)",
        )
        .leader(
            WORM_MATE,
            [0.9 * l, 0.],
            [1.05 * l, d],
            Notation::Length,
            "Wheel teeth",
        )
}

/// The diagram of the selected gear family with its staged values.
pub(super) fn diagram(variant: &Variant, v: &Values) -> Diagram {
    match variant.family {
        GearFamily::Rack => rack(variant, v),
        GearFamily::Bevel => bevel(variant, v),
        GearFamily::Worm => worm(variant, v),
        _ => cylindrical(variant, v),
    }
}

/// Family selector artwork without dimensions.
pub(super) fn art(family: GearFamily) -> Diagram {
    let mut d = diagram(
        &Variant {
            family,
            hand: GearHand::Right,
            pitch_input: GearPitchInput::Module,
            worm_role: WormMember::Worm,
        },
        &Values {
            teeth: Some(if family == GearFamily::Rack { 6. } else { 16. }),
            ..Default::default()
        },
    );
    d.callouts.clear();
    d.marks
        .retain(|m| !matches!(m, super::super::diagram::Mark::Note { .. }));
    d
}
