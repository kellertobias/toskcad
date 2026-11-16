use super::*;
use crate::authored_hanger_bolt::{BoundaryOptions, Hanger};
#[test]
fn complete_connected_field_has_real_annular_shoulders_and_closed_coedges() {
    for (left, nb) in [(false, false), (true, true)] {
        let s = Hanger::new(crate::authored_hanger_bolt::tests::fixture(left, nb), None).unwrap();
        let b = s.boundary(BoundaryOptions::default(), None).unwrap();
        let field = EncodedField::prepare(&b, CurveOptions::default(), None).unwrap();
        let domains = FiniteDomains::prepare(&field, DomainOptions::default(), None).unwrap();
        let closed = ClosedField::prepare(&domains, DomainOptions::default(), None).unwrap();
        let holes = closed
            .faces()
            .iter()
            .map(|f| f.loops().len() - 1)
            .sum::<usize>();
        assert_eq!(holes, 2);
        assert_eq!(
            closed.vertices().len() + closed.faces().len(),
            closed.edges().len() + holes + 2
        );
        eprintln!(
            "closed left={left} nonbinary={nb} V={} E={} F={} holes={holes} work={} bytes={}",
            closed.vertices().len(),
            closed.edges().len(),
            closed.faces().len(),
            closed.work_used(),
            closed.retained_bytes()
        );
    }
}
fn reject(candidate: &ClosedField<'_, '_, '_>) {
    let mut w = Work::new(16_000_000, None).unwrap();
    let mut m = Memory::new(candidate.bytes + SCRATCH + 1024 * 1024, 64 * 1024 * 1024).unwrap();
    assert!(matches!(
        candidate.certify(&mut m, &mut w),
        Err(Error::Unresolved)
    ));
}
#[test]
fn actual_generated_graph_tampering_cannot_claim_owned_closed_material() {
    let source = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, false),
        None,
    )
    .unwrap();
    let boundary = source.boundary(BoundaryOptions::default(), None).unwrap();
    let field = EncodedField::prepare(&boundary, CurveOptions::default(), None).unwrap();
    let domains = FiniteDomains::prepare(&field, DomainOptions::default(), None).unwrap();
    let mut c = ClosedField::prepare(&domains, DomainOptions::default(), None).unwrap();
    let saved = c.vertices.pop().unwrap();
    reject(&c);
    c.vertices.push(saved);
    let original = c.vertices[0];
    c.vertices[0] = ClosedVertex::Thread(usize::MAX);
    reject(&c);
    c.vertices[0] = original;
    let original = c.vertices[1];
    c.vertices[1] = c.vertices[0];
    reject(&c);
    c.vertices[1] = original;
    let plain = domains.vertices().len();
    let original = c.vertices[plain];
    c.vertices[plain] = ClosedVertex::Plain {
        end: End::Start,
        quarter: 4,
    };
    reject(&c);
    c.vertices[plain] = original;
    let original = c.edges[1];
    c.edges[1].curve = c.edges[0].curve;
    reject(&c);
    c.edges[1] = original;
    let original = c.edges[0];
    c.edges[0].vertices[0] = usize::MAX;
    reject(&c);
    c.edges[0] = original;
    c.edges[0].curve = ClosedCurve::Thread(usize::MAX);
    reject(&c);
    c.edges[0] = original;
    c.edges[0].curve = ClosedCurve::PlainRail { quarter: 4 };
    reject(&c);
    c.edges[0] = original;
    let original = c.faces[0].kind;
    c.faces[0].kind = ClosedFaceKind::ThreadCell(usize::MAX);
    reject(&c);
    c.faces[0].kind = original;
    let original = c.faces[1].kind;
    c.faces[1].kind = c.faces[0].kind;
    reject(&c);
    c.faces[1].kind = original;
    c.faces[0].same_sense = !c.faces[0].same_sense;
    reject(&c);
    c.faces[0].same_sense = !c.faces[0].same_sense;
    let loops = std::mem::take(&mut c.faces[0].loops);
    reject(&c);
    c.faces[0].loops = loops;
    c.faces[0].loops[0].swap(0, 1);
    reject(&c);
    c.faces[0].loops[0].swap(0, 1);
    let plane = c
        .faces
        .iter()
        .position(|f| f.kind == ClosedFaceKind::Plane(FaceKey::MetricTip))
        .unwrap();
    let original = c.faces[plane].kind;
    c.faces[plane].kind = ClosedFaceKind::Plane(FaceKey::PlainWall);
    reject(&c);
    c.faces[plane].kind = original;
    c.faces[plane].same_sense = !c.faces[plane].same_sense;
    reject(&c);
    c.faces[plane].same_sense = !c.faces[plane].same_sense;
    let wire = std::mem::take(&mut c.faces[plane].loops[0]);
    c.faces[plane].loops[0].extend_from_slice(&wire[..2]);
    reject(&c);
    c.faces[plane].loops[0] = wire;
    let shoulder = c
        .faces
        .iter()
        .position(|f| f.kind == ClosedFaceKind::Plane(FaceKey::MetricShoulder))
        .unwrap();
    let inner = c.faces[shoulder].loops.pop().unwrap();
    reject(&c);
    c.faces[shoulder].loops.push(inner);
    let wall = c
        .faces
        .iter()
        .position(|f| f.kind == ClosedFaceKind::PlainWall { quarter: 0 })
        .unwrap();
    let original = c.faces[wall].kind;
    c.faces[wall].kind = ClosedFaceKind::PlainWall { quarter: 4 };
    reject(&c);
    c.faces[wall].kind = original;
    c.faces[wall].same_sense = true;
    reject(&c);
    c.faces[wall].same_sense = false;
    // Separate metric torus plus plain/wood sphere: closed opposite coedges
    // and global Euler2 are insufficient without connected face adjacency.
    let moved = c.faces[shoulder].loops.pop().unwrap();
    c.faces[plane].loops.push(moved);
    reject(&c);
    let moved = c.faces[plane].loops.pop().unwrap();
    c.faces[shoulder].loops.push(moved);
    let mut w = Work::new(16_000_000, None).unwrap();
    let mut m = Memory::new(c.bytes + SCRATCH, 64 * 1024 * 1024).unwrap();
    c.certify(&mut m, &mut w).unwrap();
}
#[test]
fn link_and_boundary_cycle_certificates_reject_branches_contacts_and_open_chains() {
    let mut w = Work::new(16_000_000, None).unwrap();
    let mut m = Memory::new(SCRATCH, 64 * 1024 * 1024).unwrap();
    for bad in [
        &[][..],
        &[[1, 2], [2, 3]][..],
        &[[1, 2], [2, 1], [3, 4], [4, 3]][..],
        &[[1, 2], [2, 1], [1, 3], [3, 1]][..],
    ] {
        assert!(matches!(
            check_link(bad, &mut m, &mut w),
            Err(Error::Unresolved)
        ));
    }
    check_link(&[[1, 2], [2, 3], [3, 1]], &mut m, &mut w).unwrap();
    let edges = [
        ClosedEdge {
            vertices: [0, 1],
            curve: ClosedCurve::PlainRail { quarter: 0 },
        },
        ClosedEdge {
            vertices: [1, 2],
            curve: ClosedCurve::PlainRail { quarter: 1 },
        },
        ClosedEdge {
            vertices: [2, 0],
            curve: ClosedCurve::PlainRail { quarter: 2 },
        },
    ];
    let mut query = |values: Vec<ClosedCoedge>| {
        let mut memory = Memory::new(SCRATCH, 64 * 1024 * 1024).unwrap();
        memory
            .add(values.capacity() * size_of::<ClosedCoedge>())
            .unwrap();
        ordered(values, &edges, &mut memory, &mut w)
    };
    assert!(matches!(query(vec![]), Err(Error::Unresolved)));
    assert!(matches!(
        query(vec![ClosedCoedge {
            edge: 0,
            sense: true
        }]),
        Err(Error::Unresolved)
    ));
    assert!(matches!(
        query(vec![
            ClosedCoedge {
                edge: 0,
                sense: true
            },
            ClosedCoedge {
                edge: 2,
                sense: true
            }
        ]),
        Err(Error::Unresolved)
    ));
    assert!(matches!(
        query(vec![
            ClosedCoedge {
                edge: 0,
                sense: true
            },
            ClosedCoedge {
                edge: 1,
                sense: true
            },
            ClosedCoedge {
                edge: 2,
                sense: true
            },
            ClosedCoedge {
                edge: 0,
                sense: true
            }
        ]),
        Err(Error::Unresolved)
    ));
}
#[test]
fn support_bindings_and_measured_cap_allocation_frontiers_are_not_bypassed() {
    let source = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, false),
        None,
    )
    .unwrap();
    let boundary_ = source.boundary(BoundaryOptions::default(), None).unwrap();
    let field = EncodedField::prepare(&boundary_, CurveOptions::default(), None).unwrap();
    let domains = FiniteDomains::prepare(&field, DomainOptions::default(), None).unwrap();
    let mut c = ClosedField::prepare(&domains, DomainOptions::default(), None).unwrap();
    assert!(std::ptr::eq(c.domains(), &domains));
    let mut w = Work::new(16_000_000, None).unwrap();
    assert!(matches!(
        boundary(
            &domains,
            domains
                .edges()
                .iter()
                .find(|e| e.faces[1].is_some())
                .unwrap(),
            &mut w
        ),
        Err(Error::Unresolved)
    ));
    let plane_bytes = c
        .faces
        .iter()
        .filter(|f| matches!(f.kind, ClosedFaceKind::Plane(_)))
        .map(|f| {
            f.loops.capacity() * size_of::<Vec<ClosedCoedge>>()
                + f.loops
                    .iter()
                    .map(|l| l.capacity() * size_of::<ClosedCoedge>())
                    .sum::<usize>()
        })
        .sum::<usize>();
    let before_caps = c.bytes + SCRATCH - plane_bytes;
    let o = DomainOptions::default();
    for bytes in [
        before_caps + size_of::<Vec<ClosedCoedge>>() - 1,
        before_caps
            + size_of::<Vec<ClosedCoedge>>()
            + domains.edges().len() * size_of::<ClosedCoedge>()
            - 1,
    ] {
        assert!(matches!(
            ClosedField::prepare(
                &domains,
                DomainOptions {
                    maximum_bytes: bytes,
                    ..o
                },
                None
            ),
            Err(Error::Limit)
        ));
    }
    let mut bind_reject = |candidate: &ClosedField<'_, '_, '_>| {
        assert!(matches!(
            candidate.bind_faces(&mut w),
            Err(Error::Unresolved)
        ));
    };
    let saved = c.faces[0].loops[0][0];
    let own = c.faces[0].loops[0]
        .iter()
        .map(|x| x.edge)
        .collect::<Vec<_>>();
    let foreign = c
        .faces
        .iter()
        .skip(1)
        .flat_map(|f| f.loops.iter())
        .flat_map(|l| l.iter())
        .copied()
        .find(|x| !own.contains(&x.edge))
        .unwrap();
    c.faces[0].loops[0][0] = foreign;
    bind_reject(&c);
    c.faces[0].loops[0][0] = saved;
    let wall = c
        .faces
        .iter()
        .position(|f| f.kind == ClosedFaceKind::PlainWall { quarter: 0 })
        .unwrap();
    let saved = c.faces[wall].loops[0][0];
    c.faces[wall].loops[0][0] = foreign;
    bind_reject(&c);
    c.faces[wall].loops[0][0] = saved;
    let plane = c
        .faces
        .iter()
        .position(|f| f.kind == ClosedFaceKind::Plane(FaceKey::MetricTip))
        .unwrap();
    let saved = c.faces[plane].loops[0][0];
    let internal = domains
        .edges()
        .iter()
        .position(|e| e.faces[1].is_some())
        .unwrap();
    c.faces[plane].loops[0][0].edge = internal;
    bind_reject(&c);
    c.faces[plane].loops[0][0].edge = domains.edges().len() + 8;
    bind_reject(&c);
    c.faces[plane].loops[0][0] = saved;
    let shoulder = c
        .faces
        .iter()
        .position(|f| f.kind == ClosedFaceKind::Plane(FaceKey::MetricShoulder))
        .unwrap();
    let saved = c.faces[shoulder].loops[0][0];
    c.faces[shoulder].loops[0][0] = c.faces[plane].loops[0][0];
    bind_reject(&c);
    c.faces[shoulder].loops[0][0] = saved;
    // A second traversal of an otherwise real cap cycle is a double cover.
    let extra = c.faces[plane].loops[0].clone();
    c.faces[plane].loops[0].extend_from_slice(&extra);
    reject(&c);
    c.faces[plane].loops[0].truncate(extra.len());
    let mut spec = source.spec();
    spec.plain_radius_mm = spec.metric.maximum_radius_mm().next_up();
    let near = Hanger::new(spec, None).unwrap();
    let b = near.boundary(BoundaryOptions::default(), None).unwrap();
    let f = EncodedField::prepare(&b, CurveOptions::default(), None).unwrap();
    let d = FiniteDomains::prepare(&f, o, None).unwrap();
    assert!(matches!(
        ClosedField::prepare(&d, o, None),
        Err(Error::Unresolved)
    ));
}
#[test]
fn initial_closed_field_ledger_validates_and_polls_before_any_publication() {
    let source = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, false),
        None,
    )
    .unwrap();
    let b = source.boundary(BoundaryOptions::default(), None).unwrap();
    let f = EncodedField::prepare(&b, CurveOptions::default(), None).unwrap();
    let d = FiniteDomains::prepare(&f, DomainOptions::default(), None).unwrap();
    let o = DomainOptions::default();
    assert!(matches!(
        ClosedField::prepare(
            &d,
            DomainOptions {
                maximum_bytes: 1,
                ..o
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        ClosedField::prepare(&d, o, Some(&|| true)),
        Err(Error::Cancelled)
    ));
}
