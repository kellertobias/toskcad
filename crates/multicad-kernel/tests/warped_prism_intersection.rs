//! A box through the flat web of a formed sheet meets exactly the web
//! material in its footprint, opening and all; a box reaching a bend is not
//! claimed by the flat-zone proof.
use multicad_core::{document::EntityId, references::Origin};
use multicad_kernel::{
    BooleanOperation, GeometryPolicy, PlaneFrame, RigidTransform, Solid,
    rust_kernel::{
        normal_sweep::Edge as ProfileEdge,
        warped_union::{SheetFlange, SheetWebEdge, formed_sheet},
    },
};
use std::{collections::BTreeMap, f64::consts::PI};

fn policy() -> GeometryPolicy {
    GeometryPolicy::default()
}

/// A 10 × 8 web, 1 mm thick, with a hole of radius 1 at (5, 4) and two
/// flanges bent down from its long edges between x = 2 and 8.
fn formed() -> Solid {
    let corners = [
        [0., 0.],
        [2., 0.],
        [8., 0.],
        [10., 0.],
        [10., 8.],
        [8., 8.],
        [2., 8.],
        [0., 8.],
    ];
    let web: Vec<_> = (0..corners.len())
        .map(|i| SheetWebEdge {
            ends: [corners[i], corners[(i + 1) % corners.len()]],
            source: i,
            flange: match i {
                1 => Some(0),
                5 => Some(1),
                _ => None,
            },
        })
        .collect();
    let line = |a, b| ProfileEdge::Line { a, b };
    let arc = |a, b, clockwise| ProfileEdge::Arc {
        center: [0., -1.],
        a,
        b,
        clockwise,
    };
    let section: Vec<_> = [
        arc([0., 0.], [1., -1.], true),
        line([1., -1.], [1., -3.]),
        line([1., -3.], [2., -3.]),
        line([2., -3.], [2., -1.]),
        arc([2., -1.], [0., 1.], false),
        line([0., 1.], [0., 0.]),
    ]
    .into_iter()
    .enumerate()
    .map(|(i, edge)| (edge, i + 10))
    .collect();
    let flanges = [[20, 21], [22, 23]].map(|ends| SheetFlange {
        section: &section,
        ends,
    });
    let hole = ProfileEdge::Circle {
        center: [5., 4.],
        radius_mm: 1.,
        clockwise: true,
    };
    let shell = formed_sheet(
        &web,
        &[vec![(hole, 50)]],
        [30, 31],
        &flanges,
        1.,
        1e-7,
        1e-5,
    )
    .unwrap();
    let feature = EntityId::new();
    let origins: Vec<_> = (0..64)
        .map(|i| Origin {
            source: feature,
            role: format!("face/{i}"),
        })
        .collect();
    Solid::rust_place_formed_sheet(
        shell,
        PlaneFrame::default(),
        feature,
        &origins,
        [30, 31],
        &BTreeMap::new(),
        policy(),
    )
    .unwrap()
}

fn block(min: [f64; 3], size: [f64; 3]) -> Solid {
    Solid::box_solid(size, policy())
        .unwrap()
        .transformed(RigidTransform {
            translation_mm: min,
            ..Default::default()
        })
        .unwrap()
}

fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-9 * b.abs().max(1.), "{a} != {b}");
}

#[test]
fn a_box_through_the_web_keeps_its_opening_and_web_thickness() {
    let sheet = formed();
    for (probe, area) in [
        // Around the whole hole.
        // The bends' cylinders keep a clear margin from y = 2.5 .. 5.5.
        (block([3., 2.5, -1.], [4., 3., 3.]), 12. - PI),
        // Through part of it, 0.5 mm off its centre: the clip keeps the
        // circular trim, less the circular segment beyond x = 5.5.
        (
            block([5.5, 2.5, -1.], [1.5, 3., 3.]),
            4.5 - (PI / 3. - 3_f64.sqrt() / 4.),
        ),
        // Plain web beside it.
        (block([7.5, 2.5, -1.], [1., 3., 3.]), 3.),
    ] {
        for (a, b) in [(&sheet, &probe), (&probe, &sheet)] {
            let common = a
                .boolean(b, BooleanOperation::Intersect, policy())
                .unwrap()
                .unwrap();
            near(common.volume_mm3().unwrap(), area);
            let bounds = common.bounds().unwrap();
            near(bounds.min[2], 0.);
            near(bounds.max[2], 1.);
        }
    }
    // A box beside the web holds none of it.
    let apart = block([3., 2.5, 1.5], [4., 3., 3.]);
    assert!(!matches!(
        sheet.boolean(&apart, BooleanOperation::Intersect, policy()),
        Ok(Some(_))
    ));
}

#[test]
fn a_box_over_a_bend_is_not_a_flat_zone() {
    let sheet = formed();
    // Spans the web edge at y = 0 and the start of the bend below it.
    let probe = block([3., -1.5, -1.5], [4., 3., 3.]);
    if let Ok(Some(common)) = sheet.boolean(&probe, BooleanOperation::Intersect, policy()) {
        // Another exact provider may answer, never with the flat-web volume.
        assert!((common.volume_mm3().unwrap() - 4. * 1.5).abs() > 1e-3);
    }
}
