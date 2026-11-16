use super::*;
use crate::authored_hanger_bolt::{BoundaryOptions, Hanger};
use crate::step::hanger::encoded::cap::tests::bezier;
use num_traits::Signed;
#[test]
fn actual_published_affine_cap_buffers_bind_original_event_points_and_keep_source_owner() {
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
            events.bound_affine_cap(usize::MAX, CurveOptions::default(), None),
            Err(Error::Invalid)
        ));
        let mut positive = 0;
        let mut roots = 0;
        for id in 0..events.segments().len() {
            let Some(cap) = events
                .bound_affine_cap(id, CurveOptions::default(), None)
                .unwrap()
            else {
                roots += 1;
                continue;
            };
            positive += 1;
            assert!(std::ptr::eq(cap.source(), &events));
            assert_eq!(cap.segment_id(), id);
            assert!(cap.original_phase_reversed());
            assert!(cap.retained_bytes() > 0);
            let [c0, d1, original] = cap.source_bounds();
            assert!(c0 <= 1e-7 && d1 <= 1e-6 && original <= 1e-6);
            assert_eq!(cap.points().len(), cap.weights().len());
            assert_eq!(cap.surface_points().len(), cap.surface_weights().len());
            let mut w = Work::new(16_000_000, None).unwrap();
            let weights = cap
                .weights()
                .iter()
                .map(|x| x.number(&mut w).unwrap().0)
                .collect::<Vec<_>>();
            let ids = events.segments()[id].event_ids();
            for e in 0..2 {
                let t = cap.parameter_range()[e].number(&mut w).unwrap().0;
                let denominator =
                    bezier(weights.iter().map(|x| [x.clone()]).collect(), &t)[0].clone();
                let target = events.coordinates(ids[1 - e]).unwrap();
                let mut error = N::integer(0);
                for axis in 0..3 {
                    let numerators = cap
                        .points()
                        .iter()
                        .zip(&weights)
                        .map(|(p, weight)| p[axis].number(&mut w).unwrap().0 * weight)
                        .collect::<Vec<_>>();
                    let actual = bezier(numerators.into_iter().map(|x| [x]).collect(), &t)[0]
                        .clone()
                        / &denominator;
                    let delta = N((actual - target[axis].number(&mut w).unwrap().0).abs());
                    error = error.add(&delta, &mut w).unwrap();
                }
                assert!(error.enclosure(&mut w).unwrap().hi <= cap.vertex_incidence_error_mm());
            }
            if positive == 1 && !left {
                use std::cell::Cell;
                let calls = Cell::new(0_usize);
                let count = || {
                    calls.set(calls.get() + 1);
                    false
                };
                events
                    .bound_affine_cap(id, CurveOptions::default(), Some(&count))
                    .unwrap()
                    .unwrap();
                for checkpoint in 1..=calls.get().min(4096) {
                    let current = Cell::new(0_usize);
                    let cancel = || {
                        current.set(current.get() + 1);
                        current.get() == checkpoint
                    };
                    assert!(
                        matches!(
                            events.bound_affine_cap(id, CurveOptions::default(), Some(&cancel)),
                            Err(Error::Cancelled)
                        ),
                        "checkpoint={checkpoint}"
                    );
                }
            }
            if positive == 1 {
                let used = cap.work_used();
                assert!(
                    events
                        .bound_affine_cap(
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
                    events.bound_affine_cap(
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
                    events.bound_affine_cap(
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
                    events.bound_affine_cap(
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
                    events.bound_affine_cap(id, CurveOptions::default(), Some(&|| true)),
                    Err(Error::Cancelled)
                ));
            }
        }
        assert_eq!(positive, if nb { 43 } else { 40 });
        assert_eq!(roots, if nb { 5 } else { 4 });
    }
}
