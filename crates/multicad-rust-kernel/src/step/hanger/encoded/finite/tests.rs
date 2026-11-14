use super::*;
use crate::authored_hanger_bolt::{BoundaryOptions, Hanger};
#[test]
fn actual_exact_finite_domains_cover_both_component_cylinders_and_share_all_seams() {
    for (left, nb) in [(false, false), (true, true)] {
        let source =
            Hanger::new(crate::authored_hanger_bolt::tests::fixture(left, nb), None).unwrap();
        let boundary = source.boundary(BoundaryOptions::default(), None).unwrap();
        let field = EncodedField::prepare(&boundary, CurveOptions::default(), None).unwrap();
        let domains = FiniteDomains::prepare(&field, DomainOptions::default(), None).unwrap();
        eprintln!(
            "finite left={left} nonbinary={nb} faces={} vertices={} edges={} boundary={} work={} bytes={}",
            domains.cells().len(),
            domains.vertices().len(),
            domains.edges().len(),
            domains
                .edges()
                .iter()
                .filter(|e| e.faces[1].is_none())
                .count(),
            domains.work_used(),
            domains.retained_bytes()
        );
        for component in [Component::Metric, Component::Wood] {
            let v = domains
                .vertices()
                .iter()
                .filter(|v| v.component == component)
                .count();
            let e = domains
                .edges()
                .iter()
                .filter(|e| domains.vertices()[e.vertices[0]].component == component)
                .count();
            let f = domains
                .cells()
                .iter()
                .filter(|c| {
                    field.meridians[field.cells[c.source_cell].meridian].component == component
                })
                .count();
            assert_eq!(v + f, e, "actual sewn side annulus Euler0");
        }
    }
}
#[test]
fn exact_candidate_orientation_incidence_and_domain_loss_are_rejected() {
    let mut w = Work::new(16_000_000, None).unwrap();
    let rectangle = || {
        vec![
            Point {
                q: N::integer(0),
                v: N::integer(0),
            },
            Point {
                q: N::integer(3),
                v: N::integer(0),
            },
            Point {
                q: N::integer(3),
                v: N::integer(1),
            },
            Point {
                q: N::integer(0),
                v: N::integer(1),
            },
        ]
    };
    let mut reversed = rectangle();
    reversed.reverse();
    assert!(matches!(
        positive_polygon(&reversed, &mut w),
        Err(Error::Unresolved)
    ));
    assert!(positive_polygon(&[], &mut w).unwrap().is_none());
    let tangent = clip(
        rectangle(),
        &N::integer(0),
        &N::integer(1),
        &N::integer(0),
        true,
        &mut w,
    )
    .unwrap();
    assert!(positive_polygon(&tangent, &mut w).unwrap().is_none());
    let too_many = vec![(0, 0), (1, 0), (2, 0), (3, 0), (3, 1), (2, 1), (0, 1)]
        .into_iter()
        .map(|(q, v)| Point {
            q: N::integer(q),
            v: N::integer(v),
        })
        .collect::<Vec<_>>();
    assert!(matches!(
        positive_polygon(&too_many, &mut w),
        Err(Error::Unresolved)
    ));
    let mut edges = Vec::with_capacity(4);
    assert!(matches!(
        add_edge(&mut edges, [0, 0], 0, &mut w),
        Err(Error::Unresolved)
    ));
    add_edge(&mut edges, [0, 1], 0, &mut w).unwrap();
    assert!(matches!(
        add_edge(&mut edges, [0, 1], 1, &mut w),
        Err(Error::Unresolved)
    ));
    add_edge(&mut edges, [1, 0], 1, &mut w).unwrap();
    assert!(matches!(
        add_edge(&mut edges, [1, 0], 2, &mut w),
        Err(Error::Unresolved)
    ));
    let source = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, false),
        None,
    )
    .unwrap();
    let boundary = source.boundary(BoundaryOptions::default(), None).unwrap();
    let field = EncodedField::prepare(&boundary, CurveOptions::default(), None).unwrap();
    assert!(matches!(
        check_totals(&field, &[N::integer(1), N::integer(1)], &mut w),
        Err(Error::Unresolved)
    ));
    assert!(matches!(
        FiniteDomains::prepare(
            &field,
            DomainOptions {
                maximum_work: 1,
                ..DomainOptions::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        FiniteDomains::prepare(&field, DomainOptions::default(), Some(&|| true)),
        Err(Error::Cancelled)
    ));
}
#[test]
fn exact_clipper_normalizes_only_identical_vertices_and_all_first_cell_checkpoints_cancel() {
    let mut w = Work::new(16_000_000, None).unwrap();
    let point = || Point {
        q: N::integer(0),
        v: N::integer(0),
    };
    let mut out = Vec::with_capacity(6);
    append(&mut out, point(), &mut w).unwrap();
    append(&mut out, point(), &mut w).unwrap();
    assert_eq!(out.len(), 1);
    assert!(
        clip(
            vec![],
            &N::integer(0),
            &N::integer(1),
            &N::integer(0),
            false,
            &mut w
        )
        .unwrap()
        .is_empty()
    );
    let closed = vec![
        point(),
        Point {
            q: N::integer(3),
            v: N::integer(0),
        },
        Point {
            q: N::integer(3),
            v: N::integer(1),
        },
        Point {
            q: N::integer(0),
            v: N::integer(1),
        },
        point(),
    ];
    assert_eq!(
        clip(
            closed,
            &N::integer(0),
            &N::integer(1),
            &N::integer(-1),
            false,
            &mut w
        )
        .unwrap()
        .len(),
        4
    );
    let source = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, true),
        None,
    )
    .unwrap();
    let boundary = source.boundary(BoundaryOptions::default(), None).unwrap();
    let field = EncodedField::prepare(&boundary, CurveOptions::default(), None).unwrap();
    for prefix in 0..1200 {
        let n = std::cell::Cell::new(0);
        assert!(matches!(
            FiniteDomains::prepare(
                &field,
                DomainOptions::default(),
                Some(&|| {
                    n.set(n.get() + 1);
                    n.get() > prefix
                })
            ),
            Err(Error::Cancelled)
        ));
    }
}
#[test]
fn crossing_append_cancellation_propagates_from_exact_comparison() {
    let rectangle = || {
        vec![
            Point {
                q: N::integer(0),
                v: N::integer(0),
            },
            Point {
                q: N::integer(3),
                v: N::integer(0),
            },
            Point {
                q: N::integer(3),
                v: N::integer(1),
            },
            Point {
                q: N::integer(0),
                v: N::integer(1),
            },
        ]
    };
    let count = std::cell::Cell::new(0);
    let callback = || {
        count.set(count.get() + 1);
        false
    };
    let mut w = Work::new(16_000_000, Some(&callback)).unwrap();
    clip(
        rectangle(),
        &N::integer(0),
        &N::integer(1),
        &N(num_rational::BigRational::new(1.into(), 2.into())),
        false,
        &mut w,
    )
    .unwrap();
    let maximum = count.get();
    for prefix in 0..maximum {
        let count = std::cell::Cell::new(0);
        let callback = || {
            count.set(count.get() + 1);
            count.get() > prefix
        };
        let mut w = Work::new(16_000_000, Some(&callback)).unwrap();
        let limit = N(num_rational::BigRational::new(1.into(), 2.into()));
        assert!(matches!(
            clip(
                rectangle(),
                &N::integer(0),
                &N::integer(1),
                &limit,
                false,
                &mut w
            ),
            Err(Error::Cancelled)
        ));
    }
}
