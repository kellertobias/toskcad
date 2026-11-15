use super::*;
use crate::authored_hanger_bolt::{BoundaryOptions, Hanger};
use crate::step::hanger::encoded::cap::tests::bezier;
use num_traits::Signed;
#[test]
fn published_root_cap_geometry_binds_original_events_with_both_endpoint_algorithms() {
    for (left, nb) in [(false, false), (true, true)] {
        let source =
            Hanger::new(crate::authored_hanger_bolt::tests::fixture(left, nb), None).unwrap();
        let boundary = source.boundary(BoundaryOptions::default(), None).unwrap();
        let field = EncodedField::prepare(&boundary, CurveOptions::default(), None).unwrap();
        let domains = FiniteDomains::prepare(&field, DomainOptions::default(), None).unwrap();
        let closed = ClosedField::prepare(&domains, DomainOptions::default(), None).unwrap();
        let vertices = closed
            .literal_vertices(CurveOptions::default(), None)
            .unwrap();
        let events = vertices.cap_events(DomainOptions::default(), None).unwrap();
        assert!(matches!(
            events.bound_root_cap(usize::MAX, CurveOptions::default(), None),
            Err(Error::Invalid)
        ));
        let mut positive = 0;
        let mut affine = 0;
        for id in 0..events.segments().len() {
            let Some(cap) = events
                .bound_root_cap(id, CurveOptions::default(), None)
                .unwrap()
            else {
                affine += 1;
                continue;
            };
            positive += 1;
            assert!(std::ptr::eq(cap.source(), &events));
            assert_eq!(cap.segment_id(), id);
            assert!(cap.retained_bytes() > 0);
            let [c0, d1, original] = cap.source_bounds();
            assert!(c0 <= 1e-7 && d1 <= 1e-6 && original <= 1e-6);
            assert!(cap.field_bounds()[0] <= 1e-7 && cap.field_bounds()[1] <= 1e-6);
            assert!(cap.certificate_work_used() <= cap.work_used());
            assert!(cap.stage_work().iter().sum::<usize>() <= cap.certificate_work_used());
            assert_eq!(cap.points().len(), cap.weights().len());
            assert_eq!(cap.surface_points().len(), cap.surface_weights().len());
            check_endpoints(&cap, &events, id);
            if positive == 1 {
                let used = cap.work_used();
                assert!(
                    events
                        .bound_root_cap(
                            id,
                            CurveOptions {
                                work: used,
                                ..Default::default()
                            },
                            None
                        )
                        .unwrap()
                        .is_some()
                );
                assert!(matches!(
                    events.bound_root_cap(
                        id,
                        CurveOptions {
                            work: used - 1,
                            ..Default::default()
                        },
                        None
                    ),
                    Err(Error::Limit)
                ));
                assert!(matches!(
                    events.bound_root_cap(
                        id,
                        CurveOptions {
                            position_error_mm: 1e-30,
                            ..Default::default()
                        },
                        None
                    ),
                    Err(Error::Unresolved)
                ));
                assert!(matches!(
                    events.bound_root_cap(
                        id,
                        CurveOptions {
                            bytes: events.retained_bytes() + 1024 * 1024,
                            ..Default::default()
                        },
                        None
                    ),
                    Err(Error::Limit)
                ));
                assert!(matches!(
                    events.bound_root_cap(id, CurveOptions::default(), Some(&|| true)),
                    Err(Error::Cancelled)
                ));
                let exact = events
                    .bound_root_cap_with_arithmetic(
                        id,
                        CurveOptions::default(),
                        EndpointArithmetic::Exact,
                        None,
                    )
                    .unwrap()
                    .unwrap();
                check_endpoints(&exact, &events, id);
            }
        }
        assert_eq!(positive, if nb { 5 } else { 4 });
        assert_eq!(affine, if nb { 43 } else { 40 });
    }
}
fn check_endpoints(
    cap: &BoundRootCap<'_, '_, '_, '_, '_, '_>,
    events: &CapEvents<'_, '_, '_, '_, '_>,
    id: usize,
) {
    let mut w = Work::new(16_000_000, None).unwrap();
    let weights = cap
        .weights()
        .iter()
        .map(|x| x.number(&mut w).unwrap().0)
        .collect::<Vec<_>>();
    let ids = events.segments()[id].event_ids();
    for (e, &event_id) in ids.iter().enumerate() {
        let t = cap.parameter_range()[e].number(&mut w).unwrap().0;
        let denominator = bezier(weights.iter().map(|x| [x.clone()]).collect(), &t)[0].clone();
        let target = events.coordinates(event_id).unwrap();
        let mut error = N::integer(0);
        for axis in 0..3 {
            let numerators = cap
                .points()
                .iter()
                .zip(&weights)
                .map(|(p, weight)| [p[axis].number(&mut w).unwrap().0 * weight])
                .collect::<Vec<_>>();
            let actual = bezier(numerators, &t)[0].clone() / &denominator;
            let delta = N((actual - target[axis].number(&mut w).unwrap().0).abs());
            error = error.add(&delta, &mut w).unwrap();
        }
        assert!(error.enclosure(&mut w).unwrap().hi <= cap.vertex_incidence_error_mm());
    }
}
