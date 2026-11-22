use super::*;
#[test]
fn all_actual_authored_sources_have_complete_embedded_domain_certificates() {
    for index in 0..6 {
        let (s, d) = crate::phillips_drive::test_fixture(index);
        let source = Cutter::new(s, d).unwrap();
        let certificate = source
            .closed_boundary_certificate(Default::default(), None)
            .unwrap();
        assert!(std::ptr::eq(certificate.source(), &source));
        assert_eq!(certificate.topology().euler_characteristic, 2);
        assert_eq!(certificate.chart_topology().euler_characteristic, 2);
        assert_eq!(
            (
                certificate.topology().vertices,
                certificate.topology().edges,
                certificate.topology().faces
            ),
            (49, 76, 29)
        );
        let atlas = source.surface_patches().unwrap();
        assert_eq!(certificate.patch_domains().len(), atlas.patches().len());
        let mut covered = std::collections::BTreeSet::new();
        let mut edge_uses = std::collections::BTreeMap::<usize, Vec<bool>>::new();
        let mut collapsed = std::collections::BTreeMap::<usize, usize>::new();
        for (i, domain) in certificate.patch_domains().iter().enumerate() {
            let p = &atlas.patches()[i];
            assert_eq!(domain.source_patch, p.source_identity());
            assert_eq!(
                domain.source_face,
                atlas.topology().faces()[p.face_index()].source_identity()
            );
            assert_eq!(domain.winding, atlas.patch_winding(i).unwrap());
            assert_eq!(domain.normalized_domain, [[0., 1.]; 2]);
            covered.insert(p.face_index());
            for (side, positive) in [(2, true), (1, true), (3, false), (0, false)] {
                match domain.sides[side] {
                    SourceChartSide::Curve(id) => edge_uses
                        .entry(id)
                        .or_default()
                        .push(positive == (domain.winding > 0)),
                    SourceChartSide::Collapsed(id) => *collapsed.entry(id).or_default() += 1,
                }
            }
        }
        assert_eq!(covered.len(), 29);
        assert_eq!(edge_uses.len(), certificate.chart_topology().edges);
        assert!(edge_uses.values().all(|v| v.len() == 2 && v[0] != v[1]));
        assert_eq!(collapsed.len(), 10);
        let mut multiplicities = collapsed.values().copied().collect::<Vec<_>>();
        multiplicities.sort_unstable();
        assert_eq!(multiplicities, [1, 1, 1, 1, 1, 1, 1, 1, 4, 36]);
        assert_eq!(
            certificate
                .patch_domains()
                .iter()
                .filter(|p| p.degenerate_contact_parameter)
                .count(),
            4
        );
        for trim in certificate.trim_domains() {
            for b in [
                trim.parameter_width,
                trim.y_decrease_mm,
                trim.depth_increase_mm,
            ] {
                assert!(b[0] > 0. && b[0] <= b[1]);
            }
        }
        let r = certificate.reserves();
        for b in [
            r.front_intercept_mm,
            r.front_slope_below_one,
            r.adjacent_corner_diagonal_mm,
            r.finite_far_mm,
            r.finite_start_mm,
            r.finite_end_mm,
        ]
        .into_iter()
        .chain(r.depth_transversality)
        .chain(r.cone_footpoints_mm)
        .chain(r.rounded_pole_mm2)
        {
            assert!(b[0] > 0. && b[0] <= b[1]);
        }
        let witness = certificate.outward_orientation();
        assert_eq!(witness.normalized_point, [0.5; 2]);
        assert_eq!(
            certificate.patch_domains()[witness.patch_index].source_patch,
            PatchKey::FloorFan { quarter: 0 }
        );
        let t = s.floor_angle_deg.to_radians().tan();
        let independent = std::f64::consts::PI / 64. * (t + 1. / t);
        assert!(witness.positive_gradient_dot[0] > 0.);
        assert!(
            witness.positive_gradient_dot[0] <= independent
                && independent <= witness.positive_gradient_dot[1]
        );
        assert!(certificate.work_used() < MAX_WORK);
        assert!(certificate.retained_budget_bytes() < MAX_BYTES);
        eprintln!(
            "sourceH{index} closed {}V/{}E/{}F domains{} collapsed{} work{}",
            certificate.chart_topology().vertices,
            certificate.chart_topology().edges,
            certificate.chart_topology().faces,
            certificate.patch_domains().len(),
            collapsed.len(),
            certificate.work_used()
        );
    }
}

#[test]
fn whole_domain_obstructions_refuse_without_modifying_source() {
    let (s, d) = crate::phillips_drive::test_fixture(0);
    for (spec, domain) in [
        (
            s,
            Domain {
                corner_start_depth_mm: 0.,
                ..d
            },
        ),
        (
            s,
            Domain {
                corner_end_depth_mm: 0.01,
                ..d
            },
        ),
        (
            s,
            Domain {
                corner_reach_mm: 0.001,
                ..d
            },
        ),
        (
            Spec {
                tilt_angle_deg: 30.,
                ..s
            },
            d,
        ),
        (
            Spec {
                inner: Inner::ShallowV {
                    angle_deg: 90.1,
                    e_width_mm: 0.435,
                },
                ..s
            },
            d,
        ),
    ] {
        let source = Cutter::new(spec, domain).unwrap();
        let old = source.clone();
        assert!(matches!(
            source.closed_boundary_certificate(Default::default(), None),
            Err(Error::Unresolved)
        ));
        assert_eq!(source, old);
    }
    // These valid source domains enter the inner diagonal formula, but their
    // entry intersects an inner corner: current finite atlas event order refuses.
    for inner in [
        Inner::Rounded { radius_mm: 0.6 },
        Inner::ShallowV {
            angle_deg: 138.,
            e_width_mm: 0.8,
        },
    ] {
        let spec = Spec {
            mouth_diameter_mm: 0.94,
            mouth_radius_mm: 0.01,
            b_width_mm: 0.85,
            g_diameter_mm: 0.90,
            inner,
            ..s
        };
        let source = Cutter::new(spec, d).unwrap();
        assert!(matches!(
            source.closed_boundary_certificate(Default::default(), None),
            Err(Error::Unresolved)
        ));
    }
    // No branch is guessed when the exact diagonal enclosure straddles its join.
    let seed = Cutter::new(
        Spec {
            mouth_diameter_mm: 0.94,
            mouth_radius_mm: 0.01,
            b_width_mm: 0.85,
            g_diameter_mm: 0.90,
            inner: Inner::ShallowV {
                angle_deg: 138.,
                e_width_mm: 0.8,
            },
            ..s
        },
        d,
    )
    .unwrap();
    let e = seed.envelope.entry.midpoint() * 2_f64.sqrt();
    let source = Cutter::new(
        Spec {
            inner: Inner::ShallowV {
                angle_deg: 138.,
                e_width_mm: e,
            },
            ..seed.spec()
        },
        d,
    )
    .unwrap();
    assert!(matches!(
        source.closed_boundary_certificate(Default::default(), None),
        Err(Error::Unresolved)
    ));
}
#[test]
fn resource_boundaries_and_every_source_cancellation_stage_are_atomic() {
    let (s, d) = crate::phillips_drive::test_fixture(0);
    let source = Cutter::new(s, d).unwrap();
    let before = source.clone();
    for options in [
        ClosedBoundaryOptions {
            maximum_work: 0,
            ..Default::default()
        },
        ClosedBoundaryOptions {
            maximum_work: MAX_WORK + 1,
            ..Default::default()
        },
        ClosedBoundaryOptions {
            maximum_bytes: 0,
            ..Default::default()
        },
        ClosedBoundaryOptions {
            maximum_bytes: MAX_BYTES + 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            source.closed_boundary_certificate(options, None),
            Err(Error::Invalid)
        ));
    }
    let peak = std::mem::size_of::<ClosedBoundaryCertificate<'_>>()
        + std::mem::size_of::<Cutter>()
        + SCRATCH;
    assert!(matches!(
        source.closed_boundary_certificate(
            ClosedBoundaryOptions {
                maximum_bytes: peak - 1,
                ..Default::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    let c = source
        .closed_boundary_certificate(
            ClosedBoundaryOptions {
                maximum_bytes: peak,
                ..Default::default()
            },
            None,
        )
        .unwrap();
    let cost = c.work_used();
    assert!(matches!(
        source.closed_boundary_certificate(
            ClosedBoundaryOptions {
                maximum_work: cost - 1,
                ..Default::default()
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert_eq!(
        source
            .closed_boundary_certificate(
                ClosedBoundaryOptions {
                    maximum_work: cost,
                    ..Default::default()
                },
                None
            )
            .unwrap()
            .work_used(),
        cost
    );
    let calls = std::cell::Cell::new(0);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    source
        .closed_boundary_certificate(Default::default(), Some(&count))
        .unwrap();
    for prefix in 0..calls.get() {
        let n = std::cell::Cell::new(0);
        let cancel = || {
            let i = n.get();
            n.set(i + 1);
            i >= prefix
        };
        assert!(
            matches!(
                source.closed_boundary_certificate(Default::default(), Some(&cancel)),
                Err(Error::Cancelled)
            ),
            "prefix{prefix}/{}",
            calls.get()
        );
    }
    assert_eq!(source, before);
    assert_eq!(
        source
            .closed_boundary_certificate(Default::default(), None)
            .unwrap()
            .reserves(),
        c.reserves()
    );
}
#[test]
fn legitimate_quotient_and_nonregular_parameter_boundaries_remain_explicit() {
    let (s, d) = crate::phillips_drive::test_fixture(0);
    let source = Cutter::new(s, d).unwrap();
    let c = source
        .closed_boundary_certificate(Default::default(), None)
        .unwrap();
    let atlas = source.surface_patches().unwrap();
    for (i, domain) in c.patch_domains().iter().enumerate() {
        if matches!(
            domain.source_patch,
            PatchKey::Envelope {
                envelope: EnvelopeRole::EntryRoundover,
                ..
            }
        ) {
            let p = atlas.patch_bounds(i, [[0.; 2], [0.5; 2]]).unwrap();
            assert!(p.du_mm[2][0] <= 0. && p.du_mm[2][1] >= 0.);
            assert!(p.normal[2][0] > 0.);
        }
        if domain.degenerate_contact_parameter {
            let p = atlas.patch_bounds(i, [[1.; 2], [0.5; 2]]).unwrap();
            assert!(p.du_mm[2][0] <= 0. && p.du_mm[2][1] >= 0.);
            assert!(matches!(domain.sides[1], SourceChartSide::Curve(_)));
        }
        if matches!(domain.source_patch, PatchKey::FloorFan { .. }) {
            assert!(matches!(domain.sides[1], SourceChartSide::Collapsed(_)));
        }
    }
    let tip = source.envelope.tip.midpoint();
    for z in [0.1 * tip, 0.5 * tip, 0.9 * tip] {
        assert_eq!(
            source.classify_point_mm([0., 0., z]).unwrap(),
            Classification::Inside
        );
        assert_eq!(
            source
                .classify_point_mm([source.envelope.entry.hi * 1.1, 0., z])
                .unwrap(),
            Classification::Outside
        );
    }
}
