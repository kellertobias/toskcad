use super::*;
use crate::authored_hanger_bolt::{BoundaryOptions, Hanger};
#[test]
fn plain_wall_reuses_actual_closed_angular_polynomial_and_cardinal_vertices() {
    let source = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, false),
        None,
    )
    .unwrap();
    let boundary = source.boundary(BoundaryOptions::default(), None).unwrap();
    let field = EncodedField::prepare(&boundary, CurveOptions::default(), None).unwrap();
    let domains = FiniteDomains::prepare(&field, DomainOptions::default(), None).unwrap();
    let closed = ClosedField::prepare(&domains, DomainOptions::default(), None).unwrap();
    let surfaces = CompactSurfaces::prepare(&field, DomainOptions::default(), None).unwrap();
    let compact =
        CompactBoundary::prepare(&closed, &surfaces, DomainOptions::default(), None).unwrap();
    let vertices = closed
        .literal_vertices(CurveOptions::default(), None)
        .unwrap();
    assert!(matches!(
        compact.plain_wall(usize::MAX, DomainOptions::default(), None),
        Err(Error::Invalid)
    ));
    let mut walls = 0;
    for face in 0..compact.faces().len() {
        let Some(wall) = compact
            .plain_wall(face, DomainOptions::default(), None)
            .unwrap()
        else {
            continue;
        };
        walls += 1;
        assert!(std::ptr::eq(wall.source(), &compact));
        assert_eq!(wall.face(), face);
        assert_eq!(wall.degrees(), [1, 10]);
        let ChartFaceKind::PlainWall { quarter } = compact.faces()[face].kind else {
            unreachable!()
        };
        assert_eq!(wall.knots()[0], f64::from(quarter) / 4.);
        assert_eq!(*wall.knots().last().unwrap(), f64::from(quarter + 1) / 4.);
        let width = wall.dimensions()[1];
        for end in 0..2 {
            for side in 0..2 {
                let q = (usize::from(quarter) + side) % 4;
                let v = closed
                    .vertices()
                    .iter()
                    .position(|v| {
                        *v == ClosedVertex::Plain {
                            end: if end == 0 { End::Start } else { End::End },
                            quarter: q as u8,
                        }
                    })
                    .unwrap();
                assert_eq!(
                    &wall.points()[end * width + side * (width - 1)],
                    vertices.vertices()[v].coordinates()
                );
            }
        }
        let mut w = Work::new(16_000_000, None).unwrap();
        let r = binary(source.spec().plain_radius_mm, &mut w).unwrap();
        let a = field
            .angular
            .knots()
            .iter()
            .position(|&k| k == wall.knots()[0])
            .unwrap();
        for (i, p) in wall.points()[..width].iter().enumerate() {
            let piece = a + (i / 10).min(wall.knots().len() - 2);
            let j = if i == width - 1 { 10 } else { i % 10 };
            for (axis, coordinate) in p.iter().enumerate().take(2) {
                assert_eq!(
                    coordinate.number(&mut w).unwrap().0,
                    r.mul(
                        &N(field.angular.coefficients()[piece][j][axis].clone()),
                        &mut w
                    )
                    .unwrap()
                    .0
                );
            }
        }
        assert!(wall.source_bounds().iter().all(|&e| e > 0. && e < 1e-6));
        assert!(wall.retained_bytes() > 0);
        assert!(matches!(
            compact.plain_wall(
                face,
                DomainOptions {
                    maximum_work: wall.work_used() - 1,
                    ..DomainOptions::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
        assert!(matches!(
            compact.plain_wall(face, DomainOptions::default(), Some(&|| true)),
            Err(Error::Cancelled)
        ));
        assert!(
            compact
                .plain_wall(
                    face,
                    DomainOptions {
                        maximum_bytes: 1,
                        ..DomainOptions::default()
                    },
                    None
                )
                .is_err()
        );
    }
    assert_eq!(walls, 4);
}
