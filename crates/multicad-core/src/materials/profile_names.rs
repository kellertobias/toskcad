//! Plain-language names for profile stock (CAD-434).
//!
//! A stock's designation is whatever the user or a supplier catalogue calls
//! it ("SHS 40×40×2", "item 0.0.026.03"). The description says what the
//! section actually is, from its canonical shape and size, so a picker never
//! shows an unexplained catalogue ID on its own.
use super::{ProfileShape, StockKind};

/// Millimetre values as people write them: no trailing zeros, at most three
/// decimals.
pub fn format_mm(value: f64) -> String {
    let text = format!("{value:.3}");
    let text = text.trim_end_matches('0').trim_end_matches('.');
    if text == "-0" {
        "0".into()
    } else {
        text.into()
    }
}

/// The family a profile section belongs to, as a user names it.
pub fn profile_family(shape: ProfileShape, wall_mm: Option<f64>, square: bool) -> &'static str {
    match (shape, wall_mm.is_some(), square) {
        (ProfileShape::Rectangular, true, true) => "Square tube",
        (ProfileShape::Rectangular, true, false) => "Rectangular tube",
        (ProfileShape::Rectangular, false, true) => "Square bar",
        (ProfileShape::Rectangular, false, false) => "Flat bar",
        (ProfileShape::Round, true, _) => "Round tube",
        (ProfileShape::Round, false, _) => "Round bar",
        (ProfileShape::Angle, _, true) => "Equal angle",
        (ProfileShape::Angle, _, false) => "Unequal angle",
        (ProfileShape::ItemLine8_40x40, ..) => "Aluminium slot profile",
        (ProfileShape::Sketch, ..) => "Custom section",
        (ProfileShape::Channel, ..) => "Channel",
        (ProfileShape::Tee, ..) => "T section",
        (ProfileShape::Beam, ..) => "I/H beam",
    }
}

/// What a profile stock is, in plain words, for example "Square tube
/// 40×40×2 mm, corner radius R4 mm" or "Aluminium slot profile 40×40, 8 mm
/// slot (item Line 8 system)". None for stock that is not a profile.
pub fn profile_description(specification: &StockKind) -> Option<String> {
    let StockKind::Profile {
        width_mm,
        height_mm,
        wall_mm,
        shape,
        corner_radius_mm,
        web_mm,
        ..
    } = specification
    else {
        return None;
    };
    let (w, h) = (format_mm(*width_mm), format_mm(*height_mm));
    let square = width_mm == height_mm;
    let family = profile_family(*shape, *wall_mm, square);
    let wall = wall_mm.map(format_mm);
    let mut text = match shape {
        ProfileShape::Rectangular | ProfileShape::Angle => match &wall {
            Some(t) => format!("{family} {w}×{h}×{t} mm"),
            None => format!("{family} {w}×{h} mm"),
        },
        ProfileShape::Round => match &wall {
            Some(t) => format!("{family} Ø{w}×{t} mm"),
            None => format!("{family} Ø{w} mm"),
        },
        ProfileShape::ItemLine8_40x40 => {
            return Some(format!("{family} 40×40, 8 mm slot (item Line 8 system)"));
        }
        ProfileShape::Sketch => format!("{family} {w}×{h} mm, drawn as a sketch"),
        ProfileShape::Channel | ProfileShape::Tee | ProfileShape::Beam => {
            let mut text = format!("{family} {w}×{h} mm");
            if let Some(t) = &wall {
                text.push_str(&format!(", flange {t} mm"));
            }
            if let Some(t) = web_mm {
                text.push_str(&format!(", web {} mm", format_mm(*t)));
            }
            text
        }
    };
    if let Some(r) = corner_radius_mm {
        text.push_str(&format!(", corner radius R{} mm", format_mm(*r)));
    }
    Some(text)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn profile(
        shape: ProfileShape,
        w: f64,
        h: f64,
        wall: Option<f64>,
        r: Option<f64>,
    ) -> StockKind {
        StockKind::Profile {
            designation: "catalogue id".into(),
            length_mm: 6000.,
            width_mm: w,
            height_mm: h,
            wall_mm: wall,
            shape,
            section: None,
            corner_radius_mm: r,
            web_mm: None,
        }
    }
    fn open(shape: ProfileShape, w: f64, h: f64, flange: f64, web: f64) -> StockKind {
        let mut spec = profile(shape, w, h, Some(flange), None);
        if let StockKind::Profile { web_mm, .. } = &mut spec {
            *web_mm = Some(web);
        }
        spec
    }
    #[test]
    fn every_profile_family_gets_a_plain_description_without_its_designation() {
        use ProfileShape::*;
        for (spec, expected) in [
            (
                profile(Rectangular, 40., 40., Some(2.), None),
                "Square tube 40×40×2 mm",
            ),
            (
                profile(Rectangular, 60., 40., Some(3.), Some(6.)),
                "Rectangular tube 60×40×3 mm, corner radius R6 mm",
            ),
            (
                profile(Rectangular, 20., 20., None, None),
                "Square bar 20×20 mm",
            ),
            (
                profile(Rectangular, 40., 5., None, None),
                "Flat bar 40×5 mm",
            ),
            (
                profile(Round, 33.7, 33.7, Some(2.6), None),
                "Round tube Ø33.7×2.6 mm",
            ),
            (profile(Round, 20., 20., None, None), "Round bar Ø20 mm"),
            (
                profile(Angle, 40., 40., Some(4.), None),
                "Equal angle 40×40×4 mm",
            ),
            (
                profile(Angle, 60., 40., Some(5.), None),
                "Unequal angle 60×40×5 mm",
            ),
            (
                profile(ItemLine8_40x40, 40., 40., None, None),
                "Aluminium slot profile 40×40, 8 mm slot (item Line 8 system)",
            ),
            (
                open(Channel, 55., 100., 7.5, 4.),
                "Channel 55×100 mm, flange 7.5 mm, web 4 mm",
            ),
            (
                open(Tee, 50., 50., 6., 6.),
                "T section 50×50 mm, flange 6 mm, web 6 mm",
            ),
            (
                open(Beam, 100., 200., 8.5, 5.6),
                "I/H beam 100×200 mm, flange 8.5 mm, web 5.6 mm",
            ),
        ] {
            let text = profile_description(&spec).unwrap();
            assert_eq!(text, expected);
            assert!(!text.contains("catalogue id"));
        }
        assert_eq!(
            profile_description(&StockKind::Print {
                process: "FDM".into(),
                nozzle_mm: None,
                layer_height_mm: 0.2
            }),
            None
        );
    }
    #[test]
    fn millimetres_are_written_without_trailing_zeros() {
        assert_eq!(format_mm(40.), "40");
        assert_eq!(format_mm(2.50), "2.5");
        assert_eq!(format_mm(0.1234), "0.123");
        assert_eq!(format_mm(-0.0001), "0");
    }
}
