//! CAD-434 (Q11): parametric channel, T and I/H profile stock. Flange and
//! web thicknesses are validated with plain messages, files using them need
//! their own capability, and the resolved section is the exact sharp-cornered
//! outline with its exact area.
use super::*;
use crate::model::{ProfileSection, SectionLoop};

fn open(
    shape: ProfileShape,
    w: f64,
    h: f64,
    flange: Option<f64>,
    web: Option<f64>,
) -> LibraryDefinition {
    LibraryDefinition::Stock(Stock {
        material: EntityId::new(),
        specification: StockKind::Profile {
            designation: "UPE 100".into(),
            length_mm: 6000.,
            width_mm: w,
            height_mm: h,
            wall_mm: flange,
            shape,
            section: None,
            corner_radius_mm: None,
            web_mm: web,
        },
    })
}
fn section(definition: &LibraryDefinition) -> ProfileSection {
    let LibraryDefinition::Stock(stock) = definition else {
        unreachable!()
    };
    ProfileSection::from_stock(stock).unwrap()
}
/// Shoelace area of a counterclockwise polygon: positive when the winding
/// is counterclockwise.
fn shoelace(points: &[[f64; 2]]) -> f64 {
    points
        .iter()
        .zip(points.iter().cycle().skip(1))
        .map(|(a, b)| a[0] * b[1] - b[0] * a[1])
        .sum::<f64>()
        / 2.
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-9, "{a} != {b}");
}

#[test]
fn channel_tee_and_beam_validate_flange_and_web_in_plain_words() {
    use ProfileShape::*;
    for shape in [Channel, Tee, Beam] {
        open(shape, 50., 100., Some(8.), Some(5.))
            .validate()
            .unwrap();
    }
    for (definition, needle) in [
        (
            open(Channel, 50., 100., Some(8.), None),
            "Give a flange and a web thickness for a channel",
        ),
        (
            open(Beam, 50., 100., None, Some(5.)),
            "Give a flange and a web thickness for an I/H beam",
        ),
        (open(Tee, 50., 100., Some(8.), Some(0.)), "above 0 mm"),
        (
            open(Tee, 50., 100., Some(8.), Some(50.)),
            "The web must be thinner than the width (50 mm)",
        ),
        (
            open(Tee, 60., 20., Some(20.), Some(5.)),
            "The flange must be thinner than the height (20 mm)",
        ),
        (
            open(Beam, 60., 30., Some(15.), Some(5.)),
            "Both flanges must fit the height: flange thickness under 15 mm",
        ),
        (
            open(Rectangular, 40., 40., Some(2.), Some(2.)),
            "Only channel, T and I/H sections have a web thickness",
        ),
    ] {
        let error = definition.validate().unwrap_err().to_string();
        assert!(error.contains(needle), "{error}");
    }
}

#[test]
fn open_sections_are_exact_sharp_outlines_with_exact_areas() {
    use ProfileShape::*;
    let (w, h, tf, tw) = (50., 100., 8., 5.);
    for (shape, area, vertices) in [
        (Channel, w * h - (w - tw) * (h - 2. * tf), 8),
        (Tee, w * tf + tw * (h - tf), 8),
        (Beam, 2. * w * tf + tw * (h - 2. * tf), 12),
    ] {
        let s = section(&open(shape, w, h, Some(tf), Some(tw)));
        near(s.area_mm2(), area);
        assert!(!s.is_hollow() && s.holes().is_empty() && s.is_open_section());
        let SectionLoop::Polygon(points) = s.outer(0.) else {
            panic!("{shape:?} is straight-sided");
        };
        assert_eq!(points.len(), vertices, "{shape:?}");
        // The drawn outline encloses exactly the stated area, counterclockwise.
        near(shoelace(&points), area);
        // It fills the stock's bounding box.
        for (axis, size) in [(0, w), (1, h)] {
            let lo = points.iter().map(|p| p[axis]).fold(f64::INFINITY, f64::min);
            let hi = points
                .iter()
                .map(|p| p[axis])
                .fold(f64::NEG_INFINITY, f64::max);
            near(lo, 0.);
            near(hi, size);
        }
        assert_eq!(s.outline_points(), points);
    }
    // The web of a channel is on the −x side; T and I/H webs are centred.
    let channel = section(&open(Channel, w, h, Some(tf), Some(tw)));
    near(channel.support([-1., 0.]), 0.);
    let SectionLoop::Polygon(tee) = section(&open(Tee, w, h, Some(tf), Some(tw))).outer(0.) else {
        unreachable!()
    };
    assert!(tee.contains(&[(w - tw) / 2., 0.]) && tee.contains(&[(w + tw) / 2., 0.]));
    // Labels name the family; the web is not kept on other shapes.
    assert!(channel.label().contains("channel"));
    let tube = section(&open(Rectangular, 40., 40., Some(2.), None));
    assert_eq!(tube.web_mm, None);
}

#[test]
fn open_sections_store_their_web_and_need_their_own_capability() {
    let mut document = Document::new("Open sections").unwrap();
    let plain = open(ProfileShape::Rectangular, 40., 40., Some(2.), None)
        .entity(EntityId::new(), "SHS")
        .unwrap();
    document.insert(plain).unwrap();
    assert!(!requires_open_section_profile(&document));
    let LibraryDefinition::Stock(stock) =
        open(ProfileShape::Beam, 100., 200., Some(8.5), Some(5.6))
    else {
        unreachable!()
    };
    let json = serde_json::to_value(&stock.specification).unwrap();
    assert_eq!(json["shape"], "beam");
    assert_eq!(json["web_mm"], 5.6);
    let back: StockKind = serde_json::from_value(json).unwrap();
    assert_eq!(back, stock.specification);
    // Other shapes write no web field at all.
    let LibraryDefinition::Stock(tube) = open(ProfileShape::Rectangular, 40., 40., Some(2.), None)
    else {
        unreachable!()
    };
    assert!(
        serde_json::to_value(&tube.specification)
            .unwrap()
            .get("web_mm")
            .is_none()
    );
    let beam = LibraryDefinition::Stock(stock)
        .entity(EntityId::new(), "IPE 200")
        .unwrap();
    document.insert(beam).unwrap();
    assert!(requires_open_section_profile(&document));
}
