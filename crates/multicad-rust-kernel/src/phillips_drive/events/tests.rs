use super::*;
use multicad_precision::PortableMath as _;
use std::cell::Cell;
fn zero(v: I) {
    assert!(v.valid() && v.lo <= 0. && v.hi >= 0., "{v:?}")
}
fn envelope(c: &Cutter, p: [I; 3], role: EnvelopeRole) {
    let [x, y, d] = p;
    let rho = (x.square() + y.square()).sqrt();
    let e = &c.envelope;
    match role {
        EnvelopeRole::EntryRoundover => zero(
            (rho - e.entry).square() + (d - I::point(c.spec.mouth_radius_mm)).square()
                - I::point(c.spec.mouth_radius_mm).square(),
        ),
        EnvelopeRole::OuterFlank => {
            zero(rho + d * e.sin_a / e.cos_a - I::point(c.spec.mouth_diameter_mm) / I::point(2.))
        }
        EnvelopeRole::Floor => zero(rho + (d - e.tip) * e.cos_c / e.sin_c),
    }
}
fn wall(c: &Cutter, p: [I; 3], kind: WallKind) {
    let k = &c.corners[0];
    let [x, y, d] = p;
    let s = k.cos_beta * (x - k.b_half) + k.sin_beta * (d - k.b_depth);
    match kind {
        WallKind::Outer => zero((s - k.join_s) * k.sin_theta - (y - k.join_y) * k.cos_theta),
        WallKind::Inner => match k.inner {
            SourceInner::ShallowV { sine, cosine } => zero(s * sine - y * cosine),
            SourceInner::Rounded { radius_mm } => {
                zero((s - I::point(radius_mm)).square() + y.square() - I::point(radius_mm).square())
            }
        },
    }
}
#[test]
fn metric_event_graph_owns_true_shared_join_points_and_dual_support_trims() {
    for index in 0..=5 {
        let (s, d) = super::super::tests::fixture(index);
        let c = Cutter::new(s, d).unwrap();
        let atlas = c.boundary_events().unwrap();
        assert_eq!(atlas.source(), &c);
        assert_eq!(atlas.mouth().depth_mm(), 0.);
        let mouth_vertex = atlas
            .mouth()
            .vertices()
            .iter()
            .find(|v| v.source_identity() == (0, VertexRole::UpperEnvelope))
            .unwrap();
        assert_eq!(
            atlas.event(EventRole::MouthEnvelope).point_mm_bounds()[..2],
            mouth_vertex.point_mm_bounds()
        );
        for role in [
            EventRole::MouthEnvelope,
            EventRole::EntryJoin,
            EventRole::FloorJoin,
            EventRole::WallJoin,
            EventRole::WallContact,
            EventRole::FloorApex,
        ] {
            let event = atlas.event(role);
            let p = event.point_mm_bounds().map(b);
            envelope(&c, p, event.envelope_role());
            for q in 0..4 {
                for lower in [false, true] {
                    let point = atlas.event_orbit_point_mm_bounds(role, q, lower).unwrap();
                    let [x, y] = rotate([p[0], if lower { -p[1] } else { p[1] }], q);
                    assert_eq!(point, [bounds(x), bounds(y), bounds(p[2])]);
                }
            }
        }
        envelope(
            &c,
            atlas.event(EventRole::EntryJoin).point,
            EnvelopeRole::OuterFlank,
        );
        envelope(
            &c,
            atlas.event(EventRole::FloorJoin).point,
            EnvelopeRole::Floor,
        );
        wall(&c, atlas.event(EventRole::WallJoin).point, WallKind::Inner);
        wall(&c, atlas.event(EventRole::WallJoin).point, WallKind::Outer);
        assert_eq!(
            atlas.event(EventRole::WallJoin).envelope_role(),
            if index == 0 {
                EnvelopeRole::OuterFlank
            } else {
                EnvelopeRole::Floor
            }
        );
        for (n, trim) in atlas.trims().iter().enumerate() {
            let [start, end] = trim.end_events();
            if n > 0 {
                assert_eq!(atlas.trims()[n - 1].end_events()[1], start)
            }
            let params = trim.parameter_bounds();
            if n == 0 {
                assert!(params[0][1] < params[1][0])
            } else {
                assert!(params[0][0] > params[1][1])
            }
            assert_eq!(
                atlas.trim_point_mm_bounds(n, 0.).unwrap(),
                atlas.event(start).point_mm_bounds()
            );
            assert_eq!(
                atlas.trim_point_mm_bounds(n, 1.).unwrap(),
                atlas.event(end).point_mm_bounds()
            );
            let mut last_depth = -1.;
            for j in 0..=40 {
                let t = j as f64 / 40.;
                let p = atlas.trim_point_mm_bounds(n, t).unwrap().map(b);
                envelope(&c, p, trim.envelope_role());
                wall(&c, p, trim.wall_kind());
                let midpoint = p.map(I::midpoint);
                assert!(midpoint[2] >= last_depth);
                last_depth = midpoint[2];
                if midpoint[2] > 0. {
                    let rho = midpoint[0].portable_hypot(midpoint[1]);
                    let epsilon = 1e-5;
                    let inside = [
                        midpoint[0] * (1. - epsilon / rho),
                        midpoint[1] * (1. - epsilon / rho),
                        midpoint[2],
                    ];
                    let outside = [
                        midpoint[0] * (1. + epsilon / rho),
                        midpoint[1] * (1. + epsilon / rho),
                        midpoint[2],
                    ];
                    assert!(super::super::tests::oracle(s, d, inside));
                    assert!(!super::super::tests::oracle(s, d, outside));
                }
            }
        }
    }
}
#[test]
fn source_graph_limits_unknown_event_order_and_cancel_are_explicit() {
    let (s, d) = super::super::tests::fixture(1);
    let c = Cutter::new(s, d).unwrap();
    let before = c.clone();
    for stop in 1..=4 {
        let calls = Cell::new(0);
        assert!(matches!(
            c.boundary_events_with_cancel(&|| {
                calls.set(calls.get() + 1);
                calls.get() == stop
            }),
            Err(Error::Cancelled)
        ));
    }
    let atlas = c.boundary_events().unwrap();
    assert!(matches!(
        atlas.event_orbit_point_mm_bounds(EventRole::MouthEnvelope, 4, false),
        Err(Error::Invalid)
    ));
    for (index, t) in [(4, 0.5), (0, f64::NAN), (0, -0.1), (0, 1.1)] {
        assert!(matches!(
            atlas.trim_point_mm_bounds(index, t),
            Err(Error::Invalid)
        ))
    }
    let domain = Domain {
        corner_start_depth_mm: -100.,
        corner_end_depth_mm: 1000.,
        corner_reach_mm: 1000.,
    };
    for a in [1., 8.] {
        let other = Cutter::new(
            Spec {
                outer_flank_angle_deg: a,
                ..s
            },
            domain,
        )
        .unwrap();
        assert!(matches!(other.boundary_events(), Err(Error::Unresolved)));
    }
    let tiny_angle = Cutter::new(
        Spec {
            tilt_angle_deg: 1e-300,
            ..s
        },
        d,
    )
    .unwrap();
    assert!(matches!(
        tiny_angle.boundary_events(),
        Err(Error::Unresolved)
    ));
    let broad = Cutter::new(
        Spec {
            mouth_radius_mm: 0.001,
            b_width_mm: 2.9,
            g_diameter_mm: 2.95,
            inner: Inner::ShallowV {
                angle_deg: 170.,
                e_width_mm: 2.85,
            },
            ..s
        },
        d,
    )
    .unwrap();
    assert!(matches!(broad.boundary_events(), Err(Error::Unresolved)));
    let tiny_join = Cutter::new(
        Spec {
            inner: Inner::ShallowV {
                angle_deg: 138.,
                e_width_mm: 1e-300,
            },
            ..s
        },
        d,
    )
    .unwrap();
    assert!(matches!(
        tiny_join.boundary_events(),
        Err(Error::Unresolved)
    ));
    assert_eq!(c, before);
    assert!(c.boundary_events().is_ok());
}

#[test]
fn real_coincident_floor_wall_join_and_cone_chart_conditioning_fail_closed() {
    let (s, d) = super::super::tests::fixture(1);
    let c = Cutter::new(s, d).unwrap();
    let k = &c.corners[0];
    let radius = I::point(s.g_diameter_mm) / I::point(2.);
    let beta = k.sin_beta / k.cos_beta;
    let base = k.b_half + k.b_depth * beta;
    let rb = base - c.envelope.flank * beta;
    let SourceInner::ShallowV { sine, cosine } = k.inner else {
        unreachable!()
    };
    let y =
        volume::line_intersection(j(radius), j(rb), j(cosine / sine / k.cos_beta)).unwrap()[1].v;
    let coincident = Cutter::new(
        Spec {
            inner: Inner::ShallowV {
                angle_deg: 138.,
                e_width_mm: 2. * y.midpoint(),
            },
            ..s
        },
        d,
    )
    .unwrap();
    assert!(matches!(
        coincident.boundary_events(),
        Err(Error::Unresolved)
    ));
    let p = cone_point(
        I::point(1.),
        I::point(0.1),
        I::point(1.),
        I::point(3.),
        I::point(0.2),
        I::point(0.2),
    )
    .unwrap();
    zero(p[0] + I::point(0.1) * p[2] - I::point(1.2));
    zero((p[0].square() + p[1].square()).sqrt() + p[2] - I::point(3.));
    assert!(matches!(
        cone_point(
            I::point(0.),
            I::point(0.5),
            I::point(1.),
            I::point(3.),
            I::point(0.),
            I::point(0.)
        ),
        Err(Error::Unresolved)
    ));
    assert!(matches!(
        cone_point(
            I::point(1e308),
            I::point(0.1),
            I::point(1.),
            I::point(1.),
            I::point(1e308),
            I::point(0.)
        ),
        Err(Error::Unresolved)
    ));
}

#[test]
fn inner_linear_contact_and_round_chart_domains_are_not_repaired() {
    let (s, d) = super::super::tests::fixture(1);
    let c = Cutter::new(s, d).unwrap();
    assert!(matches!(
        xy(&c, I::point(1.), I::point(1.), WallKind::Inner),
        Err(Error::Unresolved)
    ));
    let (s, d) = super::super::tests::fixture(0);
    let c = Cutter::new(s, d).unwrap();
    let k = &c.corners[0];
    let SourceInner::Rounded { radius_mm } = k.inner else {
        unreachable!()
    };
    assert!(matches!(
        evaluate::wall_offset(k, WallKind::Inner, I::new(0., radius_mm)),
        Err(Error::Unresolved)
    ));
    let offset = evaluate::wall_offset(k, WallKind::Inner, k.join_y).unwrap();
    zero(offset - k.join_s / k.cos_beta);
}
fn scaled(s: Spec, f: f64) -> Spec {
    Spec {
        mouth_diameter_mm: s.mouth_diameter_mm * f,
        mouth_radius_mm: s.mouth_radius_mm * f,
        b_width_mm: s.b_width_mm * f,
        g_diameter_mm: s.g_diameter_mm * f,
        inner: match s.inner {
            Inner::ShallowV {
                angle_deg,
                e_width_mm,
            } => Inner::ShallowV {
                angle_deg,
                e_width_mm: e_width_mm * f,
            },
            Inner::Rounded { radius_mm } => Inner::Rounded {
                radius_mm: radius_mm * f,
            },
        },
        ..s
    }
}
#[test]
fn authored_positive_features_refuse_an_underresolved_cone_join() {
    let (s, d) = super::super::tests::fixture(1);
    let beta = s.tilt_angle_deg.to_radians();
    let q = 0.9999;
    let flank = (beta.portable_tan() / q).portable_atan();
    let qf = beta.portable_tan() * s.floor_angle_deg.to_radians().portable_tan();
    let width = (s.g_diameter_mm * (q - qf) + 2e-10) / (1. - qf);
    let candidate = scaled(
        Spec {
            mouth_radius_mm: 1e-9,
            b_width_mm: width,
            outer_flank_angle_deg: flank.to_degrees(),
            inner: Inner::ShallowV {
                angle_deg: 138.,
                e_width_mm: 4e-4,
            },
            ..s
        },
        1e-158,
    );
    let c = Cutter::new(candidate, d).unwrap();
    assert!(volume::Sections::new(&c).is_ok());
    let k = &c.corners[0];
    let bt = k.sin_beta / k.cos_beta;
    let base = k.b_half + k.b_depth * bt;
    let slope = c.envelope.sin_a / c.envelope.cos_a;
    assert!(matches!(
        cone_point(
            base,
            bt,
            slope,
            I::point(c.spec.mouth_diameter_mm) / I::point(2.),
            k.join_y,
            k.join_s / k.cos_beta
        ),
        Err(Error::Unresolved)
    ));
    assert!(matches!(c.boundary_events(), Err(Error::Unresolved)));
}
