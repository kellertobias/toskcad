use super::*;
use crate::{
    constraint_sketch::{Constraint, ConstraintKind, Curve, Point, Sketch},
    document::{Document, EntityId},
    materials::{
        LibraryDefinition, ProfileShape, ProfileSketch, Stock, StockKind, requires_sketch_profile,
    },
    model::ProfileSection,
    parameters::BoundExpression,
};
use std::f64::consts::PI;

#[derive(Default)]
struct Draft {
    sketch: Sketch,
}
impl Draft {
    fn point(&mut self, p: [f64; 2]) -> EntityId {
        if let Some(existing) = self.sketch.points.iter().find(|q| q.position_mm == p) {
            return existing.id;
        }
        let id = EntityId::new();
        self.sketch.points.push(Point { id, position_mm: p });
        id
    }
    fn curve(&mut self, kind: CurveKind, construction: bool) -> EntityId {
        let id = EntityId::new();
        self.sketch.curves.push(Curve {
            id,
            construction,
            kind,
        });
        id
    }
    fn line(&mut self, a: [f64; 2], b: [f64; 2]) -> EntityId {
        let (start, end) = (self.point(a), self.point(b));
        self.curve(CurveKind::Line { start, end }, false)
    }
    fn polygon(&mut self, points: &[[f64; 2]]) {
        for i in 0..points.len() {
            self.line(points[i], points[(i + 1) % points.len()]);
        }
    }
    fn rect(&mut self, x0: f64, y0: f64, x1: f64, y1: f64) {
        self.polygon(&[[x0, y0], [x1, y0], [x1, y1], [x0, y1]]);
    }
    fn circle(&mut self, c: [f64; 2], r: f64) -> EntityId {
        let center = self.point(c);
        self.curve(
            CurveKind::Circle {
                center,
                radius_mm: r,
            },
            false,
        )
    }
    fn arc(&mut self, c: [f64; 2], a: [f64; 2], b: [f64; 2], clockwise: bool) {
        let (center, start, end) = (self.point(c), self.point(a), self.point(b));
        self.curve(
            CurveKind::Arc {
                center,
                start,
                end,
                clockwise,
            },
            false,
        );
    }
}

fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-9, "{a} != {b}");
}
fn refused(sketch: &Sketch, fragment: &str) {
    let error = resolve_sketch_section(sketch).unwrap_err().to_string();
    assert!(error.contains(fragment), "{error:?} lacks {fragment:?}");
}
fn signed_area(section_loop: &SectionLoop) -> f64 {
    sampled(section_loop)
        .windows(2)
        .chain(std::iter::once(
            &[
                *sampled(section_loop).last().unwrap(),
                sampled(section_loop)[0],
            ][..],
        ))
        .map(|w| w[0][0] * w[1][1] - w[1][0] * w[0][1])
        .sum::<f64>()
        / 2.
}

#[test]
fn a_rectangle_normalizes_to_its_bounding_box_counterclockwise() {
    let mut d = Draft::default();
    // Drawn clockwise and away from the origin.
    d.polygon(&[[10., 5.], [10., 25.], [50., 25.], [50., 5.]]);
    let s = resolve_sketch_section(&d.sketch).unwrap();
    assert_eq!(s.size_mm, [40., 20.]);
    assert_eq!(s.offset_mm, [10., 5.]);
    near(s.area_mm2, 800.);
    assert!(s.holes.is_empty());
    assert!(signed_area(&s.outer) > 0., "outer loop is counterclockwise");
    let (min, max) = bounds(&s.outer);
    assert_eq!((min, max), ([0., 0.], [40., 20.]));
}

#[test]
fn a_tube_keeps_its_hole_and_net_area() {
    let mut d = Draft::default();
    d.rect(0., 0., 40., 40.);
    d.circle([20., 20.], 10.);
    let s = resolve_sketch_section(&d.sketch).unwrap();
    assert_eq!(s.holes.len(), 1);
    assert!(matches!(s.holes[0], SectionLoop::Circle { radius_mm, .. } if radius_mm == 10.));
    near(s.area_mm2, 1600. - PI * 100.);
    // A rectangular hole is a counterclockwise loop as well.
    let mut d = Draft::default();
    d.rect(0., 0., 30., 20.);
    d.polygon(&[[5., 5.], [5., 15.], [25., 15.], [25., 5.]]);
    let s = resolve_sketch_section(&d.sketch).unwrap();
    near(s.area_mm2, 600. - 200.);
    assert!(signed_area(&s.holes[0]) > 0.);
}

#[test]
fn arcs_keep_their_direction_and_exact_extents() {
    // A 20 × 10 slot: two half circles bulging beyond the straight sides.
    let mut d = Draft::default();
    d.line([0., 0.], [20., 0.]);
    d.arc([20., 5.], [20., 0.], [20., 10.], false);
    d.line([20., 10.], [0., 10.]);
    d.arc([0., 5.], [0., 10.], [0., 0.], false);
    let s = resolve_sketch_section(&d.sketch).unwrap();
    assert_eq!(s.size_mm, [30., 10.], "arc extremes count in the size");
    near(s.area_mm2, 200. + PI * 25.);
    near(support(&s.outer, [1., 0.]), 30.);
    near(support(&s.outer, [-1., 0.]), 0.);
    let diagonal = [0.6, 0.8];
    near(support(&s.outer, diagonal), 25. * 0.6 + 5. * 0.8 + 5.);
    // A concave notch is traversed clockwise inside a counterclockwise loop.
    let mut d = Draft::default();
    d.line([0., 0.], [30., 0.]);
    d.line([30., 0.], [30., 20.]);
    d.line([30., 20.], [20., 20.]);
    d.arc([15., 20.], [20., 20.], [10., 20.], true);
    d.line([10., 20.], [0., 20.]);
    d.line([0., 20.], [0., 0.]);
    let s = resolve_sketch_section(&d.sketch).unwrap();
    near(s.area_mm2, 600. - PI * 25. / 2.);
    assert_eq!(s.size_mm, [30., 20.]);
    let SectionLoop::Curved(edges) = &s.outer else {
        panic!("curved outline")
    };
    assert!(edges.iter().any(|e| matches!(
        e,
        SectionEdge::Arc {
            clockwise: true,
            ..
        }
    )));
    assert!(signed_area(&s.outer) > 0.);
}

#[test]
fn open_crossing_separate_and_unsupported_sketches_are_refused() {
    let mut open = Draft::default();
    open.line([0., 0.], [10., 0.]);
    open.line([10., 0.], [10., 10.]);
    refused(&open.sketch, "open");
    let mut bowtie = Draft::default();
    bowtie.polygon(&[[0., 0.], [10., 10.], [10., 0.], [0., 10.]]);
    refused(&bowtie.sketch, "crosses");
    let mut two = Draft::default();
    two.rect(0., 0., 10., 10.);
    two.rect(20., 0., 30., 10.);
    refused(&two.sketch, "one closed outline");
    let mut island = Draft::default();
    island.rect(0., 0., 40., 40.);
    island.circle([20., 20.], 15.);
    island.circle([20., 20.], 5.);
    refused(&island.sketch, "Islands");
    let mut construction = Draft::default();
    construction.rect(0., 0., 10., 10.);
    for c in &mut construction.sketch.curves {
        c.construction = true;
    }
    refused(&construction.sketch, "construction geometry alone");
    refused(&Sketch::default(), "closed section outline");
    // Construction geometry next to a valid outline is ignored.
    let mut guided = Draft::default();
    guided.rect(0., 0., 10., 10.);
    let (a, b) = (guided.point([-5., 5.]), guided.point([15., 5.]));
    guided.curve(CurveKind::Line { start: a, end: b }, true);
    assert_eq!(
        resolve_sketch_section(&guided.sketch).unwrap().size_mm,
        [10., 10.]
    );
    let mut spline = Draft::default();
    spline.rect(0., 0., 10., 10.);
    let poles = vec![
        spline.point([2., 2.]),
        spline.point([4., 6.]),
        spline.point([6., 2.]),
    ];
    spline.curve(
        CurveKind::BSpline {
            poles,
            basis: crate::constraint_sketch::SplineBasis::open_uniform(2, 3).unwrap(),
        },
        false,
    );
    refused(&spline.sketch, "Splines");
}

#[test]
fn parameter_expressions_are_refused_and_literal_dimensions_kept() {
    let document = Document::new("Section expressions").unwrap();
    let mut d = Draft::default();
    d.rect(0., 0., 10., 10.);
    let line = d.sketch.curves[0].id;
    d.sketch.constraints.push(Constraint {
        id: EntityId::new(),
        reference: false,
        kind: ConstraintKind::Length {
            line,
            dimension: Dimension::Literal {
                value: 10.,
                unit: "mm".into(),
            },
        },
    });
    assert!(resolve_sketch_section(&d.sketch).is_ok());
    d.sketch.constraints[0].kind = ConstraintKind::Length {
        line,
        dimension: Dimension::Expression {
            expression: BoundExpression::bind(&document, "=10", "mm", None).unwrap(),
        },
    };
    refused(&d.sketch, "must be numbers");
}

fn stock(section: Option<Sketch>, shape: ProfileShape, size: [f64; 2]) -> LibraryDefinition {
    LibraryDefinition::Stock(Stock {
        material: EntityId::new(),
        specification: StockKind::Profile {
            designation: "Custom T".into(),
            length_mm: 3000.,
            width_mm: size[0],
            height_mm: size[1],
            wall_mm: None,
            shape,
            section: section.map(|sketch| Box::new(ProfileSketch { sketch })),
            corner_radius_mm: None,
            web_mm: None,
        },
    })
}

fn tee() -> Sketch {
    let mut d = Draft::default();
    d.polygon(&[
        [0., 0.],
        [40., 0.],
        [40., 5.],
        [22.5, 5.],
        [22.5, 30.],
        [17.5, 30.],
        [17.5, 5.],
        [0., 5.],
    ]);
    d.sketch
}

#[test]
fn stock_validation_keeps_shape_section_and_size_consistent() {
    assert!(
        stock(Some(tee()), ProfileShape::Sketch, [40., 30.])
            .validate()
            .is_ok()
    );
    let mismatch = stock(Some(tee()), ProfileShape::Sketch, [40., 31.])
        .validate()
        .unwrap_err()
        .to_string();
    assert!(
        mismatch.contains("does not match its section sketch"),
        "{mismatch}"
    );
    assert!(
        stock(None, ProfileShape::Sketch, [40., 30.])
            .validate()
            .is_err()
    );
    assert!(
        stock(Some(tee()), ProfileShape::Rectangular, [40., 30.])
            .validate()
            .is_err()
    );
    let mut open = tee();
    open.curves.pop();
    let error = stock(Some(open), ProfileShape::Sketch, [40., 30.])
        .validate()
        .unwrap_err()
        .to_string();
    assert!(error.contains("open"), "{error}");
}

#[test]
fn sketch_stock_round_trips_requires_its_capability_and_resolves_as_a_member_section() {
    let definition = stock(Some(tee()), ProfileShape::Sketch, [40., 30.]);
    let entity = definition.entity(EntityId::new(), "Custom T").unwrap();
    let json = serde_json::to_value(&definition).unwrap();
    assert_eq!(
        json.pointer("/value/specification/shape"),
        Some(&"sketch".into())
    );
    let back: LibraryDefinition = serde_json::from_value(json).unwrap();
    assert_eq!(back, definition);
    // Existing profile stock serializes exactly as before: no section key.
    let plain = stock(None, ProfileShape::Rectangular, [40., 30.]);
    let text = serde_json::to_string(&plain).unwrap();
    assert!(
        !text.contains("section") && !text.contains("shape"),
        "{text}"
    );
    let mut document = Document::new("Sketch stock").unwrap();
    assert!(!requires_sketch_profile(&document));
    document.insert(entity).unwrap();
    assert!(requires_sketch_profile(&document));
    let project = crate::project_file::ProjectFile::new(document).unwrap();
    assert!(
        project
            .capabilities()
            .contains(crate::materials::SKETCH_PROFILE_CAPABILITY)
    );
    assert!(
        crate::project_file::supported_capabilities()
            .contains(crate::materials::SKETCH_PROFILE_CAPABILITY)
    );
    let LibraryDefinition::Stock(stock) = definition else {
        unreachable!()
    };
    let section = ProfileSection::from_stock(&stock).unwrap();
    assert!(section.is_sketch() && !section.is_hollow());
    near(section.area_mm2(), 40. * 5. + 5. * 25.);
    assert!(section.label().contains("custom sketch section"));
    assert!(section.holes().is_empty());
    near(section.support([0., 1.]), 30.);
    assert_eq!(section.outline_points().len(), 8);
}
