//! CAD-434: the outside corner radius of rectangular profile stock is
//! validated against the section and wall, stored only when set, requires
//! its own file capability and shapes the resolved section exactly.
use super::*;
use crate::model::{ProfileSection, SectionEdge, SectionLoop};
use multicad_precision::PortableMath as _;

fn tube(
    w: f64,
    h: f64,
    wall: Option<f64>,
    r: Option<f64>,
    shape: ProfileShape,
) -> LibraryDefinition {
    LibraryDefinition::Stock(Stock {
        material: EntityId::new(),
        specification: StockKind::Profile {
            designation: "RHS".into(),
            length_mm: 6000.,
            width_mm: w,
            height_mm: h,
            wall_mm: wall,
            shape,
            section: None,
            corner_radius_mm: r,
        },
    })
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-9, "{a} != {b}");
}

#[test]
fn corner_radius_is_validated_against_the_section_and_wall() {
    use ProfileShape::*;
    for r in [0.5, 2., 4., 20.] {
        tube(40., 40., Some(2.), Some(r), Rectangular)
            .validate()
            .unwrap();
    }
    tube(60., 40., None, Some(20.), Rectangular)
        .validate()
        .unwrap();
    for (definition, needle) in [
        (
            tube(40., 40., Some(2.), Some(20.5), Rectangular),
            "at most half the smaller side (20 mm)",
        ),
        (
            tube(60., 30., Some(2.), Some(16.), Rectangular),
            "at most half the smaller side (15 mm)",
        ),
        (
            tube(40., 40., Some(2.), Some(0.), Rectangular),
            "above 0 mm",
        ),
        (
            tube(40., 40., Some(2.), Some(-1.), Rectangular),
            "above 0 mm",
        ),
        (
            tube(40., 40., Some(2.), Some(f64::NAN), Rectangular),
            "above 0 mm",
        ),
        (
            tube(40., 40., Some(2.), Some(3.), Round),
            "Only rectangular",
        ),
        (
            tube(40., 40., Some(4.), Some(3.), Angle),
            "Only rectangular",
        ),
    ] {
        let error = definition.validate().unwrap_err().to_string();
        assert!(error.contains(needle), "{error}");
    }
    // The derived inside radius: the outside radius less the wall.
    assert_eq!(inner_corner_radius(Some(5.), Some(2.)), Some(3.));
    assert_eq!(
        inner_corner_radius(Some(2.), Some(2.)),
        None,
        "sharp inside"
    );
    assert_eq!(
        inner_corner_radius(Some(1.), Some(2.)),
        None,
        "sharp inside"
    );
    assert_eq!(inner_corner_radius(None, Some(2.)), None);
    assert_eq!(
        inner_corner_radius(Some(5.), None),
        None,
        "a bar has no bore"
    );
}

#[test]
fn corner_radius_is_stored_only_when_set_and_requires_its_capability() {
    let sharp = tube(40., 40., Some(2.), None, ProfileShape::Rectangular);
    let text = serde_json::to_string(&sharp).unwrap();
    assert!(!text.contains("corner_radius"), "{text}");
    let rounded = tube(40., 40., Some(2.), Some(4.), ProfileShape::Rectangular);
    let json = serde_json::to_value(&rounded).unwrap();
    assert_eq!(
        json.pointer("/value/specification/corner_radius_mm"),
        Some(&4.0.into())
    );
    assert_eq!(
        serde_json::from_value::<LibraryDefinition>(json).unwrap(),
        rounded
    );
    let mut document = Document::new("Rounded tube").unwrap();
    document
        .insert(sharp.entity(EntityId::new(), "SHS").unwrap())
        .unwrap();
    assert!(!requires_profile_corner_radius(&document));
    let project = crate::project_file::ProjectFile::new(document.clone()).unwrap();
    assert!(
        !project
            .capabilities()
            .contains(CORNER_RADIUS_PROFILE_CAPABILITY)
    );
    document
        .insert(rounded.entity(EntityId::new(), "SHS R4").unwrap())
        .unwrap();
    assert!(requires_profile_corner_radius(&document));
    let project = crate::project_file::ProjectFile::new(document).unwrap();
    assert!(
        project
            .capabilities()
            .contains(CORNER_RADIUS_PROFILE_CAPABILITY)
    );
    assert!(
        crate::project_file::supported_capabilities().contains(CORNER_RADIUS_PROFILE_CAPABILITY)
    );
}

fn section(definition: &LibraryDefinition) -> ProfileSection {
    let LibraryDefinition::Stock(stock) = definition else {
        unreachable!()
    };
    ProfileSection::from_stock(stock).unwrap()
}

#[test]
fn rounded_sections_have_exact_arcs_and_areas() {
    let pi = std::f64::consts::PI;
    // RHS 60×40×3 with R5 outside: inside R2.
    let s = section(&tube(
        60.,
        40.,
        Some(3.),
        Some(5.),
        ProfileShape::Rectangular,
    ));
    near(
        s.area_mm2(),
        (60. * 40. - (4. - pi) * 25.) - (54. * 34. - (4. - pi) * 4.),
    );
    let SectionLoop::Curved(outer) = s.outer(0.) else {
        panic!("rounded outside")
    };
    let arcs: Vec<_> = outer
        .iter()
        .filter_map(|e| match *e {
            SectionEdge::Arc {
                center,
                a,
                clockwise,
                ..
            } => {
                assert!(!clockwise);
                Some(((a[0] - center[0]).portable_hypot(a[1] - center[1]), center))
            }
            SectionEdge::Line { .. } => None,
        })
        .collect();
    assert_eq!(arcs.len(), 4);
    assert!(arcs.iter().all(|(r, _)| (r - 5.).abs() < 1e-12));
    assert!(arcs.iter().any(|(_, c)| *c == [55., 35.]));
    let holes = s.holes();
    let [SectionLoop::Curved(inner)] = holes.as_slice() else {
        panic!("rounded bore")
    };
    assert!(inner.iter().any(|e| matches!(*e,
        SectionEdge::Arc { center, a, .. } if center == [55., 35.] && a == [57., 35.])));
    // The cope envelope grows the radius with the clearance.
    let SectionLoop::Curved(grown) = s.outer(1.) else {
        panic!("rounded envelope")
    };
    assert!(grown.iter().any(|e| matches!(*e,
        SectionEdge::Arc { center, a, .. } if center == [55., 5.] && a == [55., -1.])));
    // Outside radius at most the wall: the bore keeps sharp corners.
    let s = section(&tube(
        40.,
        40.,
        Some(3.),
        Some(3.),
        ProfileShape::Rectangular,
    ));
    assert!(matches!(s.holes().as_slice(), [SectionLoop::Polygon(p)] if p.len() == 4));
    near(s.area_mm2(), (1600. - (4. - pi) * 9.) - 34. * 34.);
    // A solid bar with half-height radius has round ends and no flat sides
    // to spare; a square one becomes a circle.
    let bar = section(&tube(60., 20., None, Some(10.), ProfileShape::Rectangular));
    near(bar.area_mm2(), 60. * 20. - (4. - pi) * 100.);
    let SectionLoop::Curved(edges) = bar.outer(0.) else {
        panic!("slot outline")
    };
    assert_eq!(edges.len(), 6, "two flats and four quarter arcs");
    let round = section(&tube(20., 20., None, Some(10.), ProfileShape::Rectangular));
    assert!(matches!(round.outer(0.), SectionLoop::Circle { radius_mm, .. } if radius_mm == 10.));
    near(round.area_mm2(), pi * 100.);
    // Sampled outline points carry the arcs (perimeters, extents).
    let s = section(&tube(
        40.,
        40.,
        Some(2.),
        Some(4.),
        ProfileShape::Rectangular,
    ));
    let points = s.outline_points();
    let perimeter: f64 = (0..points.len())
        .map(|i| {
            let (a, b) = (points[i], points[(i + 1) % points.len()]);
            (b[0] - a[0]).portable_hypot(b[1] - a[1])
        })
        .sum();
    let exact = 4. * (40. - 8.) + 2. * pi * 4.;
    assert!((perimeter - exact).abs() < 1e-3, "{perimeter} vs {exact}");
    assert!(s.label().contains("rounded corners"));
}
