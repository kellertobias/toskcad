use super::*;
use crate::authored_hanger_bolt::{BoundaryOptions, Hanger};
#[test]
fn actual_closed_graph_vertices_share_all_original_owner_tokens_and_cap_planes() {
    for (left, nb) in [(false, false), (true, true)] {
        let s = Hanger::new(crate::authored_hanger_bolt::tests::fixture(left, nb), None).unwrap();
        let b = s.boundary(BoundaryOptions::default(), None).unwrap();
        let f = EncodedField::prepare(&b, CurveOptions::default(), None).unwrap();
        let d = FiniteDomains::prepare(&f, DomainOptions::default(), None).unwrap();
        let c = ClosedField::prepare(&d, DomainOptions::default(), None).unwrap();
        let v = c.literal_vertices(CurveOptions::default(), None).unwrap();
        assert!(std::ptr::eq(v.source(), &c));
        assert_eq!(v.vertices().len(), c.vertices().len());
        let mut worst = 0_f64;
        for (actual, owner) in v.vertices().iter().zip(c.vertices()) {
            assert_eq!(actual.owner(), *owner);
            assert!(actual.field_incidence_error_mm() < 1e-7);
            assert!(actual.source_incidence_error_mm() < 1e-7);
            worst = worst.max(actual.source_incidence_error_mm());
            for x in actual.coordinates() {
                // Actual finite decimal emission, not a merely dyadic label.
                x.exact_step_real(DomainOptions::default(), None).unwrap();
            }
        }
        // Independent EVERY-peer interval oracle: the optimized source
        // triangle receipt must bound the actual different tensor owner.
        let mut oracle = Work::new(16_000_000, None).unwrap();
        let mut direct_max = 0_f64;
        let field_max = v
            .vertices()
            .iter()
            .map(|p| p.field_incidence_error_mm())
            .fold(0_f64, f64::max);
        for cell in d.cells() {
            let source_cell = f.cells[cell.source_cell];
            let row = &f.meridians[source_cell.meridian];
            let (pitch, starts) = if row.component == Component::Metric {
                let s = s.metric().spec();
                (s.pitch_mm, s.starts)
            } else {
                let s = s.wood().spec();
                (s.pitch_mm, s.starts)
            };
            let pitch = binary(pitch, &mut oracle).unwrap();
            for corner in cell.corners() {
                // Independently reconstruct BOTH original periodic carries;
                // exact tokens, never approximate intervals, decide identity.
                let original = &d.vertices()[corner.vertex];
                let mut q = corner.phase.number(&mut oracle).unwrap();
                let mut turn = corner.turn.number(&mut oracle).unwrap();
                let mut period = source_cell.source_key.total_period;
                if q.cmp(&pitch, &mut oracle).unwrap() == Ordering::Equal {
                    q = N::integer(0);
                    period += 1;
                }
                if turn.cmp(&N::integer(1), &mut oracle).unwrap() == Ordering::Equal {
                    turn = N::integer(0);
                    period += starts as i32;
                }
                assert_eq!(original.component, row.component);
                assert_eq!(original.period, period);
                assert_eq!(original.phase.number(&mut oracle).unwrap().0, q.0);
                assert_eq!(original.turn.number(&mut oracle).unwrap().0, turn.0);
                let (xy, z, _) = footpoint(&f, cell, corner, &mut oracle).unwrap();
                let receipt = &v.vertices()[corner.vertex];
                let direct = bind(receipt.coordinates(), xy, &z, &mut oracle).unwrap();
                direct_max = direct_max.max(direct);
                assert!(
                    direct <= receipt.field_incidence_error_mm(),
                    "direct={direct} bound={}",
                    receipt.field_incidence_error_mm()
                );
            }
        }
        assert!(v.work_used() <= 16_000_000 && v.retained_bytes() <= 64 * 1024 * 1024);
        eprintln!(
            "literal vertices left={left} nb={nb} count={} work={} bytes={} source_incidence={worst} field_receipt={field_max} direct_every_peer_oracle={direct_max}",
            v.vertices().len(),
            v.work_used(),
            v.retained_bytes()
        );
    }
}

#[test]
fn exact_source_cap_z_and_shared_counterparts_are_not_replaced_by_rounded_planes() {
    let s = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(true, true),
        None,
    )
    .unwrap();
    let b = s.boundary(BoundaryOptions::default(), None).unwrap();
    let f = EncodedField::prepare(&b, CurveOptions::default(), None).unwrap();
    let d = FiniteDomains::prepare(&f, DomainOptions::default(), None).unwrap();
    let c = ClosedField::prepare(&d, DomainOptions::default(), None).unwrap();
    let v = c.literal_vertices(CurveOptions::default(), None).unwrap();
    let mut w = Work::new(16_000_000, None).unwrap();
    for face in c.faces() {
        let ClosedFaceKind::Plane(role) = face.kind else {
            continue;
        };
        let z = match role {
            FaceKey::MetricTip => N::integer(0)
                .sub(
                    &binary(s.metric().full_length_mm(), &mut w).unwrap(),
                    &mut w,
                )
                .unwrap()
                .sub(&binary(s.spec().plain_length_mm, &mut w).unwrap(), &mut w)
                .unwrap(),
            FaceKey::MetricShoulder => N::integer(0)
                .sub(&binary(s.spec().plain_length_mm, &mut w).unwrap(), &mut w)
                .unwrap(),
            FaceKey::WoodShoulder => N::integer(0),
            FaceKey::WoodTip => binary(s.wood().full_length_mm(), &mut w).unwrap(),
            _ => panic!("only actual source plane roles may own a planar fill"),
        };
        for coedge in face.loops().iter().flatten() {
            for vertex in c.edges()[coedge.edge].vertices {
                assert_eq!(
                    v.vertices()[vertex].coordinates()[2]
                        .number(&mut w)
                        .unwrap()
                        .0,
                    z.0
                );
            }
        }
    }
    let point = v.vertices()[0].coordinates();
    let xy = [
        point[0].number(&mut w).unwrap().enclosure(&mut w).unwrap(),
        point[1].number(&mut w).unwrap().enclosure(&mut w).unwrap(),
    ];
    let z = point[2]
        .number(&mut w)
        .unwrap()
        .add(&N::integer(1), &mut w)
        .unwrap();
    assert!(matches!(
        bind(point, xy, &z, &mut w),
        Err(Error::Unresolved)
    ));
    let z = point[2].number(&mut w).unwrap();
    assert!(matches!(
        bind(point, [I::new(f64::NAN, f64::NAN); 2], &z, &mut w),
        Err(Error::Unresolved)
    ));
}

#[test]
fn original_source_precision_and_resource_frontiers_cancel_atomically() {
    use std::cell::Cell;
    let s = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, false),
        None,
    )
    .unwrap();
    let b = s.boundary(BoundaryOptions::default(), None).unwrap();
    let f = EncodedField::prepare(&b, CurveOptions::default(), None).unwrap();
    let d = FiniteDomains::prepare(&f, DomainOptions::default(), None).unwrap();
    let c = ClosedField::prepare(&d, DomainOptions::default(), None).unwrap();
    for options in [
        CurveOptions {
            work: 1,
            ..CurveOptions::default()
        },
        CurveOptions {
            bytes: 1,
            ..CurveOptions::default()
        },
    ] {
        assert!(matches!(
            c.literal_vertices(options, None),
            Err(Error::Limit)
        ));
    }
    assert!(matches!(
        c.literal_vertices(
            CurveOptions {
                work: 0,
                ..CurveOptions::default()
            },
            None
        ),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        c.literal_vertices(
            CurveOptions {
                position_error_mm: 1e-12,
                ..CurveOptions::default()
            },
            None
        ),
        Err(Error::Unresolved)
    ));
    let calls = Cell::new(0);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    let good = c
        .literal_vertices(CurveOptions::default(), Some(&count))
        .unwrap();
    // Actual allocation/reference-map prefix, then the complete bounded
    // first-point arithmetic frontier (including axial sums and quantizer).
    let first_point = 4 + d.cells().iter().map(|c| c.corners().len()).sum::<usize>();
    for stop in [1, 2, 3, calls.get() / 2, calls.get()]
        .into_iter()
        .chain(first_point..first_point + 512)
    {
        let n = Cell::new(0);
        let cancel = || {
            n.set(n.get() + 1);
            n.get() >= stop
        };
        assert!(
            matches!(
                c.literal_vertices(CurveOptions::default(), Some(&cancel)),
                Err(Error::Cancelled)
            ),
            "stop={stop}"
        );
    }
    assert_eq!(s.spec(), b.source().spec());
    assert_eq!(good.vertices().len(), c.vertices().len());
    eprintln!(
        "registry inherited field={} domains={} closed={} result={} callbackcount={}",
        f.work_used(),
        d.work_used(),
        c.work_used(),
        good.work_used(),
        calls.get()
    );
}

#[test]
fn numeric_proposals_retain_invalid_enclosure_and_exact_quantizer_failures() {
    let mut w = Work::new(16_000_000, None).unwrap();
    assert!(matches!(
        proposal(I::new(f64::NAN, f64::NAN), &mut w),
        Err(Error::Unresolved)
    ));
    for x in [0., 3., -3., f64::MIN_POSITIVE] {
        let n = proposal(p(x), &mut w).unwrap();
        assert!((n.enclosure(&mut w).unwrap() - p(x)).magnitude() <= 2_f64.powi(-64));
    }
}
