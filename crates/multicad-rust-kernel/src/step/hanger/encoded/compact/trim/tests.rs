use super::*;
use crate::authored_hanger_bolt::{BoundaryOptions, Hanger};
use std::collections::BTreeSet;
#[test]
fn finite_cap_subdivision_preserves_shared_owner_closed_oriented_shell() {
    for (left, nb) in [(false, false), (true, true)] {
        let source =
            Hanger::new(crate::authored_hanger_bolt::tests::fixture(left, nb), None).unwrap();
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
        let events = vertices.cap_events(DomainOptions::default(), None).unwrap();
        let other_closed = ClosedField::prepare(&domains, DomainOptions::default(), None).unwrap();
        let other_vertices = other_closed
            .literal_vertices(CurveOptions::default(), None)
            .unwrap();
        let other_events = other_vertices
            .cap_events(DomainOptions::default(), None)
            .unwrap();
        assert!(matches!(
            compact.trim_caps(&other_events, DomainOptions::default(), None),
            Err(Error::Invalid)
        ));
        let trimmed = compact
            .trim_caps(&events, DomainOptions::default(), None)
            .unwrap();
        assert!(std::ptr::eq(trimmed.source(), &compact));
        assert!(std::ptr::eq(trimmed.events(), &events));
        assert_eq!(trimmed.faces().len(), compact.faces().len());
        assert_eq!(
            trimmed
                .vertices()
                .iter()
                .filter(|v| matches!(v, TrimVertex::CapEvent(_)))
                .count(),
            events
                .events()
                .iter()
                .filter(|e| e.graph_vertex().is_none())
                .count()
        );
        assert!(trimmed.coordinates(usize::MAX).is_none());
        for (v, owner) in trimmed.vertices().iter().enumerate() {
            let expected = match owner {
                TrimVertex::Graph(g) => vertices.vertices()[*g].coordinates(),
                TrimVertex::CapEvent(e) => events.coordinates(*e).unwrap(),
            };
            assert_eq!(trimmed.coordinates(v).unwrap(), expected);
        }
        let mut uses = vec![Vec::new(); trimmed.edges().len()];
        let mut used_vertices = BTreeSet::new();
        let mut contribution = 0_isize;
        for (f, face) in trimmed.faces().iter().enumerate() {
            contribution += 2 - face.loops().len() as isize;
            for ring in face.loops() {
                assert!(ring.len() >= 3);
                for (i, c) in ring.iter().enumerate() {
                    uses[c.edge].push((f, c.sense));
                    let edge = &trimmed.edges()[c.edge];
                    let ends = if c.sense {
                        edge.vertices
                    } else {
                        [edge.vertices[1], edge.vertices[0]]
                    };
                    used_vertices.extend(ends);
                    let next = &ring[(i + 1) % ring.len()];
                    let next_edge = &trimmed.edges()[next.edge];
                    assert_eq!(ends[1], next_edge.vertices[usize::from(!next.sense)]);
                }
            }
        }
        assert_eq!(
            used_vertices.len() as isize - trimmed.edges().len() as isize + contribution,
            2
        );
        for owners in &uses {
            assert_eq!(owners.len(), 2);
            assert_ne!(owners[0].0, owners[1].0);
            assert_ne!(owners[0].1, owners[1].1);
        }
        let mut segments = BTreeSet::new();
        for edge in trimmed.edges() {
            if let TrimEdgeKind::CapSegment(segment) = edge.kind {
                assert!(segments.insert(segment));
                assert_eq!(
                    edge.vertices,
                    events.segments()[segment].event_ids().map(|e| {
                        match events.events()[e].graph_vertex() {
                            Some(v) => v,
                            None => trimmed
                                .vertices()
                                .iter()
                                .position(|v| *v == TrimVertex::CapEvent(e))
                                .unwrap(),
                        }
                    })
                );
            }
        }
        assert_eq!(segments.len(), events.segments().len());
        let options = DomainOptions {
            maximum_work: trimmed.work_used() - 1,
            ..DomainOptions::default()
        };
        assert!(matches!(
            compact.trim_caps(&events, options, None),
            Err(Error::Limit)
        ));
        assert!(matches!(
            compact.trim_caps(&events, DomainOptions::default(), Some(&|| true)),
            Err(Error::Cancelled)
        ));
        let options = DomainOptions {
            maximum_bytes: 1,
            ..DomainOptions::default()
        };
        assert!(compact.trim_caps(&events, options, None).is_err());
        eprintln!(
            "trimmed left={left} faces={} edges={} vertices={} work={} bytes={}",
            trimmed.faces().len(),
            trimmed.edges().len(),
            used_vertices.len(),
            trimmed.work_used(),
            trimmed.retained_bytes()
        );
    }
}
