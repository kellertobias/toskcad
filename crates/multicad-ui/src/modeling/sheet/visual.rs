//! Base-sheet schematic (CAD-395): a section along the base's u axis with
//! the stock thickness, ending in a representative 90° bend that shows the
//! inside radius and the neutral axis at K × thickness, above the base plan.
//! Thickness belongs to the shared stock and is a note, not a callout.
use super::super::diagram::{Diagram, Key, Line, Notation, legible};
use multicad_precision::PortableMath as _;

pub(in crate::modeling) const LENGTH: Key = "sheet.length";
pub(in crate::modeling) const WIDTH: Key = "sheet.width";
pub(in crate::modeling) const RADIUS: Key = "sheet.radius";
pub(in crate::modeling) const RADIUS_OVERRIDE: Key = "sheet.radius_override";
pub(in crate::modeling) const K_FACTOR: Key = "sheet.k_factor";
pub(in crate::modeling) const K_OVERRIDE: Key = "sheet.k_override";
pub(in crate::modeling) const ORIGIN: [Key; 3] =
    ["sheet.origin.x", "sheet.origin.y", "sheet.origin.z"];
pub(in crate::modeling) const BEND: [Key; 4] = [RADIUS_OVERRIDE, RADIUS, K_OVERRIDE, K_FACTOR];

#[derive(Clone, Copy, Debug, Default)]
pub(in crate::modeling) struct Values {
    /// Stock thickness, radius and K-factor, when stock is chosen.
    pub stock: Option<(f64, f64, f64)>,
    pub length: Option<f64>,
    pub width: Option<f64>,
    pub radius_override: bool,
    pub radius: Option<f64>,
    pub k_override: bool,
    pub k_factor: Option<f64>,
    pub reverse: bool,
    /// The base is a kept sketch region or converted solid, not a rectangle.
    pub kept_base: bool,
}

fn arc(center: [f32; 2], radius: f32, from_deg: f32, to_deg: f32) -> Vec<[f32; 2]> {
    (0..=12)
        .map(|i| {
            let a = (from_deg + (to_deg - from_deg) * i as f32 / 12.).to_radians();
            let (s, c) = a.portable_sin_cos();
            [center[0] + radius * c, center[1] + radius * s]
        })
        .collect()
}

pub(in crate::modeling) fn diagram(v: &Values) -> Diagram {
    let (t_mm, stock_r, stock_k) = v.stock.unwrap_or((1., 1., 0.4));
    let t = 1.;
    let r_mm = if v.radius_override {
        v.radius.unwrap_or(stock_r)
    } else {
        stock_r
    };
    let k = if v.k_override {
        v.k_factor.unwrap_or(stock_k)
    } else {
        stock_k
    }
    .clamp(0., 1.) as f32;
    let r = legible(Some(r_mm / t_mm.max(1e-6)), 1., 0.3, 5.);
    // Material above the base face, or below it when reversed; the
    // representative bend turns up from the material's top surface.
    let b = if v.reverse { -t } else { 0. };
    let flat = 9.;
    let leg = 3.5;
    let c = [flat, b + t + r];
    let inner = arc(c, r, -90., 0.);
    let outer = arc(c, r + t, -90., 0.);
    let neutral = arc(c, r + k * t, -90., 0.);
    let top = b + t + r + leg;
    let mut d = Diagram::new()
        .line([-1., 0.], [flat + 1.5, 0.], Line::Center)
        .rect(0., b, flat, b + t)
        .path(&inner, false, Line::Outline)
        .path(&outer, false, Line::Outline)
        .line([flat + r, b + t + r], [flat + r, top], Line::Outline)
        .line(
            [flat + r + t, b + t + r],
            [flat + r + t, top],
            Line::Outline,
        )
        .line([flat + r, top], [flat + r + t, top], Line::Outline)
        .path(&neutral, false, Line::Center)
        .note(
            [flat * 0.45, b + t + 1.1],
            match v.stock {
                Some((t, ..)) => format!("t {t} mm · stock"),
                None => "Choose sheet stock".into(),
            },
        );
    let (r_key, r_name) = if v.radius_override {
        (RADIUS, "Inside radius")
    } else {
        (RADIUS_OVERRIDE, "Stock radius")
    };
    let (k_key, k_name) = if v.k_override {
        (K_FACTOR, "K-factor")
    } else {
        (K_OVERRIDE, "Stock K")
    };
    d = d
        .leader(
            r_key,
            inner[6],
            [flat - 2.5, top - 0.3],
            Notation::Radius,
            r_name,
        )
        .leader(
            k_key,
            neutral[4],
            [flat + r + t + 2., b - 0.8],
            Notation::Length,
            k_name,
        );
    if v.kept_base {
        return d.note([flat * 0.5, b - 1.6], "Base kept from its sketch region");
    }
    d = d.dim(
        LENGTH,
        [0., b],
        [flat, b],
        -14.,
        Notation::Length,
        "Length (u)",
    );
    // Base plan below the section, at its own scale.
    let l = legible(v.length, 200., 1e-3, 1e7);
    let w = legible(v.width, 100., 0.15 * l, 2. * l);
    let s = flat / l;
    let (pl, pw) = (flat, w * s);
    let y0 = -t - 4.5 - pw;
    d.rect(0., y0, pl, y0 + pw)
        .dim(
            WIDTH,
            [pl, y0],
            [pl, y0 + pw],
            -14.,
            Notation::Length,
            "Width (v)",
        )
        .note([pl * 0.5, y0 + pw * 0.5], "Base plan")
}
