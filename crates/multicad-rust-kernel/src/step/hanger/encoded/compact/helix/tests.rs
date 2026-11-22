use super::*;
use crate::authored_hanger_bolt::{BoundaryOptions, Hanger};
use crate::step::hanger::encoded::cap::tests::bezier;
#[test]
fn long_helix_geometry_is_exactly_both_surface_boundary_rows_with_shared_vertices() {
    for (left, nb, metric_starts, wood_starts) in [
        (false, false, 1, 1),
        (true, true, 1, 1),
        (false, true, 2, 3),
    ] {
        let mut spec = crate::authored_hanger_bolt::tests::fixture(left, nb);
        spec.metric = crate::authored_metric_thread::AuthoredMetricThread::new(
            crate::authored_metric_thread::Spec {
                starts: metric_starts,
                ..spec.metric.spec()
            },
        )
        .unwrap();
        spec.wood = crate::authored_flat_thread::AuthoredThread::new(
            crate::helical_thread::Spec {
                starts: wood_starts,
                ..spec.wood.spec()
            },
            crate::authored_flat_thread::SourceOptions {
                position_limit_mm: 1e-7,
                normal_angle_limit_rad: 1e-6,
            },
        )
        .unwrap();
        let source = Hanger::new(spec, None).unwrap();
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
        let mut count = 0;
        assert!(matches!(
            compact.helix(usize::MAX, &vertices, DomainOptions::default(), None),
            Err(Error::Invalid)
        ));
        for index in 0..compact.edges().len() {
            let Some(curve) = compact
                .helix(index, &vertices, DomainOptions::default(), None)
                .unwrap()
            else {
                continue;
            };
            count += 1;
            assert!(std::ptr::eq(curve.source(), &compact));
            assert_eq!(curve.edge(), index);
            assert!(curve.retained_bytes() > 0);
            let mut w = Work::new(16_000_000, None).unwrap();
            assert!(curve.weight().number(&mut w).unwrap().0 > BigRational::from_integer(0.into()));
            let range = curve
                .range()
                .each_ref()
                .map(|x| x.number(&mut w).unwrap().0);
            assert_eq!(curve.same_sense(), range[0] < range[1]);
            assert!(curve.vertex_residual_mm().iter().all(|&e| e <= 1e-7));
            for k in 0..=4 {
                let t = &range[0] + (&range[1] - &range[0]) * BigRational::new(k.into(), 4.into());
                let actual = eval_points(curve.points(), curve.knots(), &t, &mut w);
                for pc in curve.pcurves() {
                    let ChartFaceKind::Thread(surface) = compact.faces()[pc.face()].kind else {
                        panic!("side pcurve must ownthreadsurface")
                    };
                    let surface = &compact.surfaces().surfaces()[surface];
                    let row = if pc.u() == 0. {
                        0
                    } else {
                        surface.dimensions()[0] - 1
                    };
                    let width = surface.dimensions()[1];
                    let pcurve_t = &t + pc.v_shift().number(&mut w).unwrap().0;
                    let peer = eval_points(
                        &surface.points()[row * width..(row + 1) * width],
                        surface.knots(),
                        &pcurve_t,
                        &mut w,
                    );
                    assert_eq!(actual, peer);
                }
            }
            if count == 1 {
                let used = curve.work_used();
                assert!(
                    compact
                        .helix(
                            index,
                            &vertices,
                            DomainOptions {
                                maximum_work: used,
                                ..Default::default()
                            },
                            None
                        )
                        .unwrap()
                        .is_some()
                );
                assert!(matches!(
                    compact.helix(
                        index,
                        &vertices,
                        DomainOptions {
                            maximum_work: used - 1,
                            ..Default::default()
                        },
                        None
                    ),
                    Err(Error::Limit)
                ));
                assert!(matches!(
                    compact.helix(index, &vertices, DomainOptions::default(), Some(&|| true)),
                    Err(Error::Cancelled)
                ));
            }
        }
        assert_eq!(count, surfaces.surfaces().len());
    }
}
fn eval_points(
    points: &[[EncodedScalar; 3]],
    knots: &[f64],
    t: &BigRational,
    w: &mut Work<'_>,
) -> [BigRational; 3] {
    let span = knots
        .windows(2)
        .position(|p| {
            BigRational::from_float(p[0]).unwrap() <= *t
                && *t <= BigRational::from_float(p[1]).unwrap()
        })
        .unwrap();
    let a = BigRational::from_float(knots[span]).unwrap();
    let b = BigRational::from_float(knots[span + 1]).unwrap();
    let normalized = (t - a.clone()) / (b - a);
    let c = points[span * 10..span * 10 + 11]
        .iter()
        .map(|p| p.each_ref().map(|x| x.number(w).unwrap().0))
        .collect();
    bezier(c, &normalized)
}
