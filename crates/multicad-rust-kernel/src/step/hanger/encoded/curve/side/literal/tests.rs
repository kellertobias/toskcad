use super::*;
use crate::authored_hanger_bolt::{BoundaryOptions, Hanger};
use crate::step::hanger::encoded::cap::tests::bezier;
use num_traits::Signed;
#[test]
fn every_real_side_support_retains_spatial_z_and_whole_curve_error_receipt() {
    for (left, nb) in [(false, false), (true, true)] {
        let source =
            Hanger::new(crate::authored_hanger_bolt::tests::fixture(left, nb), None).unwrap();
        let boundary = source.boundary(BoundaryOptions::default(), None).unwrap();
        let field = EncodedField::prepare(&boundary, CurveOptions::default(), None).unwrap();
        let mut requests = Vec::new();
        for meridian in 0..field.meridians.len() {
            let cell = field
                .cells
                .iter()
                .position(|c| c.meridian == meridian)
                .unwrap();
            for turn_end in [End::Start, End::End] {
                requests.push((cell, SideSupportKind::Meridian { turn_end }));
            }
            for angular_piece in 0..field.angular.coefficients().len() {
                for phase_end in [End::Start, End::End] {
                    requests.push((
                        cell,
                        SideSupportKind::Helix {
                            phase_end,
                            angular_piece,
                        },
                    ));
                }
            }
        }
        let batch = field
            .side_support_batch(&requests, DomainOptions::default(), None)
            .unwrap();
        assert_eq!(batch.supports().len(), 180);
        for (index, support) in batch.supports().iter().enumerate() {
            let literal = support.literal(CurveOptions::default(), None).unwrap();
            assert!(std::ptr::eq(literal.source(), support));
            assert!(literal.retained_bytes() > 0);
            assert_eq!(literal.points().len(), literal.weights().len());
            let [c0, d1] = literal.source_bounds();
            assert!(c0 <= 1e-7 && d1 <= 1e-6);
            let mut w = Work::new(16_000_000, None).unwrap();
            let z0 = literal.points()[0][2].number(&mut w).unwrap();
            let z1 = literal.points().last().unwrap()[2].number(&mut w).unwrap();
            assert_ne!(z0.0, z1.0, "real varying-Z source must not flatten");
            let original = support
                .poles()
                .iter()
                .map(|pole| {
                    let weight = pole.weight().number(&mut w).unwrap().0;
                    std::array::from_fn::<_, 4, _>(|a| {
                        if a == 3 {
                            weight.clone()
                        } else {
                            pole.point()[a].number(&mut w).unwrap().0 * &weight
                        }
                    })
                })
                .collect::<Vec<_>>();
            let actual = literal
                .points()
                .iter()
                .zip(literal.weights())
                .map(|(p, weight)| {
                    let weight = weight.number(&mut w).unwrap().0;
                    std::array::from_fn::<_, 4, _>(|a| {
                        if a == 3 {
                            weight.clone()
                        } else {
                            p[a].number(&mut w).unwrap().0 * &weight
                        }
                    })
                })
                .collect::<Vec<_>>();
            for fraction in [0, 1, 2, 3, 4] {
                let t = BigRational::new(fraction.into(), 4.into());
                let o = bezier(original.clone(), &t);
                let a = bezier(actual.clone(), &t);
                let error = (0..3)
                    .map(|axis| (&o[axis] / &o[3] - &a[axis] / &a[3]).abs())
                    .fold(BigRational::from_integer(0.into()), |a, b| a + b);
                assert!(N(error).enclosure(&mut w).unwrap().hi <= literal.field_bounds()[0]);
            }
            if index == 0 {
                let used = literal.work_used();
                assert!(
                    support
                        .literal(
                            CurveOptions {
                                work: used,
                                ..Default::default()
                            },
                            None
                        )
                        .is_ok()
                );
                assert!(matches!(
                    support.literal(
                        CurveOptions {
                            work: used - 1,
                            ..Default::default()
                        },
                        None
                    ),
                    Err(Error::Limit)
                ));
                assert!(matches!(
                    support.literal(
                        CurveOptions {
                            position_error_mm: 1e-30,
                            ..Default::default()
                        },
                        None
                    ),
                    Err(Error::Unresolved)
                ));
                assert!(matches!(
                    support.literal(CurveOptions::default(), Some(&|| true)),
                    Err(Error::Cancelled)
                ));
                assert!(matches!(
                    support.literal(
                        CurveOptions {
                            bytes: field.retained_bytes() + support.retained_bytes(),
                            ..Default::default()
                        },
                        None
                    ),
                    Err(Error::Limit)
                ));
            }
        }
    }
}
#[test]
fn actual_closed_side_edges_have_original_parameter_boundary_owners() {
    for (left, nb) in [(false, false), (true, true)] {
        let source =
            Hanger::new(crate::authored_hanger_bolt::tests::fixture(left, nb), None).unwrap();
        let boundary = source.boundary(BoundaryOptions::default(), None).unwrap();
        let field = EncodedField::prepare(&boundary, CurveOptions::default(), None).unwrap();
        let domains = FiniteDomains::prepare(&field, DomainOptions::default(), None).unwrap();
        let closed = ClosedField::prepare(&domains, DomainOptions::default(), None).unwrap();
        let mut cap = vec![false; domains.edges().len()];
        for face in closed.faces() {
            if matches!(face.kind, ClosedFaceKind::Plane(_)) {
                for use_ in face.loops().iter().flatten() {
                    if let ClosedCurve::Thread(edge) = closed.edges()[use_.edge].curve {
                        cap[edge] = true;
                    }
                }
            }
        }
        let (mut helices, mut meridians, mut root_meridians) = (0, 0, 0);
        let mut w = Work::new(16_000_000, None).unwrap();
        for edge in closed.edges() {
            let ClosedCurve::Thread(id) = edge.curve else {
                continue;
            };
            if cap[id] {
                continue;
            }
            let domain_edge = &domains.edges()[id];
            assert!(domain_edge.faces[1].is_some());
            let (face, _) = domain_edge.faces[0].unwrap();
            let cell = &domains.cells()[face];
            let points = domain_edge
                .vertices
                .map(|vertex| cell.corners().iter().find(|c| c.vertex == vertex).unwrap());
            let q = points.map(|p| p.phase.number(&mut w).unwrap());
            let v = points.map(|p| p.turn.number(&mut w).unwrap());
            let row = &field.meridians[field.cells[cell.source_cell].meridian];
            if q[0].0 == q[1].0 {
                helices += 1;
                assert!(
                    q[0].0
                        == row.phase[0]
                            .number(&mut w)
                            .unwrap()
                            .div(&row.weights[0].number(&mut w).unwrap(), &mut w)
                            .unwrap()
                            .0
                        || q[0].0
                            == row.phase[row.degree()]
                                .number(&mut w)
                                .unwrap()
                                .div(&row.weights[row.degree()].number(&mut w).unwrap(), &mut w)
                                .unwrap()
                                .0
                );
                assert_ne!(v[0].0, v[1].0)
            } else {
                assert_eq!(v[0].0, v[1].0);
                assert!(v[0].0 == N::integer(0).0 || v[0].0 == N::integer(1).0);
                meridians += 1;
                if row.degree() == 2 {
                    root_meridians += 1
                }
            }
        }
        eprintln!(
            "actual side ownership left={left} nb={nb} helix_edges={helices} meridian_edges={meridians} root_meridians={root_meridians}"
        );
        assert!(helices > 0 && meridians > 0 && root_meridians > 0);
    }
}
