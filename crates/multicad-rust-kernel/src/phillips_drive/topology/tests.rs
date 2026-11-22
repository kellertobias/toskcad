use super::*;
use multicad_precision::PortableMath as _;
use std::cell::Cell;
fn zero(x: I) {
    assert!(x.valid() && x.lo <= 0. && x.hi >= 0., "{x:?}")
}
fn residual(surface: Surface, p: [I; 3]) {
    match surface {
        Surface::Plane {
            normal_bounds,
            offset_mm_bounds,
        } => zero(
            normal_bounds
                .map(i)
                .into_iter()
                .zip(p)
                .fold(I::point(0.), |s, (n, x)| s + n * x)
                - i(offset_mm_bounds),
        ),
        Surface::Cylinder {
            origin_mm_bounds,
            axis_bounds,
            radius_mm,
            outward_sign,
        } => {
            assert_eq!(outward_sign, -1);
            let origin = origin_mm_bounds.map(i);
            let axis = axis_bounds.map(i);
            zero(axis.into_iter().fold(I::point(0.), |s, v| s + v.square()) - I::point(1.));
            let delta = std::array::from_fn::<_, 3, _>(|j| p[j] - origin[j]);
            let projection = delta
                .into_iter()
                .zip(axis)
                .fold(I::point(0.), |s, (d, a)| s + d * a);
            zero(
                delta.into_iter().fold(I::point(0.), |s, v| s + v.square())
                    - projection.square()
                    - I::point(radius_mm).square(),
            );
        }
        Surface::Torus {
            major_radius_mm_bounds,
            minor_radius_mm,
            center_depth_mm,
            outward_sign,
        } => {
            assert_eq!(outward_sign, -1);
            let rho = (p[0].square() + p[1].square()).sqrt();
            zero(
                (rho - i(major_radius_mm_bounds)).square()
                    + (p[2] - I::point(center_depth_mm)).square()
                    - I::point(minor_radius_mm).square(),
            );
        }
        Surface::Cone {
            radius_at_zero_mm_bounds,
            radius_decrease_per_depth_bounds,
        } => zero(
            (p[0].square() + p[1].square()).sqrt() + p[2] * i(radius_decrease_per_depth_bounds)
                - i(radius_at_zero_mm_bounds),
        ),
    }
}
#[test]
fn metric_closed_topologies_own_real_shared_cross_surface_curves_and_manifold_apex() {
    for fixture in 0..=5 {
        let (s, d) = super::super::tests::fixture(fixture);
        let cutter = Cutter::new(s, d).unwrap();
        let body = cutter.boundary_topology().unwrap();
        let cert = body.certificate();
        assert_eq!(
            cert,
            TopologyCertificate {
                vertices: 49,
                edges: 76,
                faces: 29,
                euler_characteristic: 2
            }
        );
        assert_eq!(body.events().source(), &cutter);
        for j in 0..20 {
            assert_eq!(
                body.vertices()[j].point_mm_bounds()[..2],
                body.events().mouth().vertices()[j].point_mm_bounds()
            );
        }
        let mut incidents = vec![Vec::new(); body.edges().len()];
        for (face_id, face) in body.faces().iter().enumerate() {
            let _key = face.source_identity();
            let surface = body.face_surface(face_id).unwrap();
            for u in face.edge_uses() {
                incidents[u.edge].push((face_id, u.forward));
                for t in [0., 0.25, 0.5, 0.75, 1.] {
                    residual(
                        surface,
                        body.curve_point_mm_bounds(u.edge, t).unwrap().map(i),
                    );
                }
            }
        }
        for (index, edge) in body.edges().iter().enumerate() {
            let _key = edge.source_identity();
            if let Some(line) = body.curve_line_support(index).unwrap() {
                let direction = line.direction_bounds.map(i);
                zero(
                    direction
                        .into_iter()
                        .fold(I::point(0.), |s, v| s + v.square())
                        - I::point(1.),
                );
                assert_eq!(
                    line.origin_mm_bounds,
                    body.vertices()[edge.vertices()[0]].point_mm_bounds()
                );
                let a = body.curve_point_mm_bounds(index, 0.25).unwrap().map(i);
                let b = body.curve_point_mm_bounds(index, 0.75).unwrap().map(i);
                let delta = std::array::from_fn::<_, 3, _>(|j| b[j] - a[j]);
                for j in 0..3 {
                    zero(
                        delta[(j + 1) % 3] * direction[(j + 2) % 3]
                            - delta[(j + 2) % 3] * direction[(j + 1) % 3],
                    )
                }
                assert!(
                    delta
                        .into_iter()
                        .zip(direction)
                        .fold(I::point(0.), |s, (d, n)| s + d * n)
                        .lo
                        > 0.
                );
            }
            let [a, b] = edge.vertices();
            assert_eq!(
                body.curve_point_mm_bounds(index, 0.).unwrap(),
                body.vertices()[a].point_mm_bounds()
            );
            assert_eq!(
                body.curve_point_mm_bounds(index, 1.).unwrap(),
                body.vertices()[b].point_mm_bounds()
            );
            assert_eq!(incidents[index].len(), 2);
            assert_ne!(incidents[index][0].1, incidents[index][1].1);
        }
        let apices: Vec<_> = body
            .vertices()
            .iter()
            .enumerate()
            .filter(|(_, v)| {
                matches!(
                    v.source_identity(),
                    VertexKey::Event {
                        role: EventRole::FloorApex,
                        ..
                    }
                )
            })
            .collect();
        assert_eq!(apices.len(), 1);
        let apex = apices[0].0;
        let generators: Vec<_> = body
            .edges()
            .iter()
            .filter(|e| matches!(e.source_identity(), CurveKey::FloorGenerator { .. }))
            .collect();
        assert_eq!(generators.len(), 4);
        for edge in generators {
            assert_eq!(edge.vertices()[1], apex);
            assert_ne!(edge.vertices()[0], apex)
        }
        for quarter in 0..4 {
            for t in [0., 0.25, 0.75, 1.] {
                let n = body.floor_apex_normal_limit(quarter, t).unwrap().map(i);
                zero(n.into_iter().fold(I::point(0.), |s, v| s + v.square()) - I::point(1.));
                let angle = (quarter as f64 + t) * PI / 2.;
                let tangent = [
                    I::point(angle.portable_cos()),
                    I::point(angle.portable_sin()),
                    -cutter.envelope.sin_c / cutter.envelope.cos_c,
                ];
                zero(
                    n.into_iter()
                        .zip(tangent)
                        .fold(I::point(0.), |s, (n, v)| s + n * v),
                );
            }
        }
        assert_ne!(
            body.floor_apex_normal_limit(0, 0.).unwrap(),
            body.floor_apex_normal_limit(0, 1.).unwrap()
        );
    }
}
#[test]
fn source_topology_queries_limits_and_cancellation_are_atomic() {
    let (s, d) = super::super::tests::fixture(1);
    let source = Cutter::new(s, d).unwrap();
    let before = source.clone();
    let calls = Cell::new(0);
    let body = source
        .boundary_topology_with_cancel(&|| {
            calls.set(calls.get() + 1);
            false
        })
        .unwrap();
    let total = calls.get();
    for stop in 1..=total {
        let n = Cell::new(0);
        assert!(matches!(
            source.boundary_topology_with_cancel(&|| {
                n.set(n.get() + 1);
                n.get() == stop
            }),
            Err(Error::Cancelled)
        ))
    }
    for (index, t) in [(76, 0.5), (0, -0.1), (0, 1.1), (0, f64::NAN)] {
        assert!(matches!(
            body.curve_point_mm_bounds(index, t),
            Err(Error::Invalid)
        ))
    }
    assert!(matches!(body.face_surface(29), Err(Error::Invalid)));
    assert!(matches!(body.curve_line_support(76), Err(Error::Invalid)));
    for (q, t) in [(4, 0.5), (0, -0.1), (0, 1.1), (0, f64::NAN)] {
        assert!(matches!(
            body.floor_apex_normal_limit(q, t),
            Err(Error::Invalid)
        ))
    }
    assert_eq!(source, before);
    assert!(source.boundary_topology().is_ok());
    // Lookup failure is a source identity error, never a default to another edge.
    assert!(
        vertex_id(
            body.vertices(),
            VertexKey::Mouth {
                quarter: 4,
                role: VertexRole::Axis
            }
        )
        .is_err()
    );
    assert!(edge_id(body.edges(), CurveKey::Mouth { edge: 99 }).is_err());
    assert!(mouth_edge(body.events(), 4, EdgeRole::Envelope).is_err());
}

#[test]
fn edit_across_event_order_preserves_owner_keys_and_never_reassigns_absent_trim() {
    let (s, d) = super::super::tests::fixture(1);
    let old = Cutter::new(s, d).unwrap();
    let old = old.boundary_topology().unwrap();
    let changed = Cutter::new(
        Spec {
            inner: Inner::ShallowV {
                angle_deg: 138.,
                e_width_mm: 0.55,
            },
            ..s
        },
        d,
    )
    .unwrap();
    let changed = changed.boundary_topology().unwrap();
    let common = CurveKey::Trim {
        envelope: EnvelopeRole::OuterFlank,
        wall: WallKind::Outer,
        quarter: 0,
        lower: false,
    };
    assert!(old.edge_index(common).is_some());
    assert!(changed.edge_index(common).is_some());
    let removed = CurveKey::Trim {
        envelope: EnvelopeRole::Floor,
        wall: WallKind::Outer,
        quarter: 0,
        lower: false,
    };
    assert!(old.edge_index(removed).is_some());
    assert_eq!(changed.edge_index(removed), None);
    let added = CurveKey::Trim {
        envelope: EnvelopeRole::OuterFlank,
        wall: WallKind::Inner,
        quarter: 0,
        lower: false,
    };
    assert_eq!(old.edge_index(added), None);
    assert!(changed.edge_index(added).is_some());
}

#[test]
fn schema_mapping_and_publication_reject_missing_roles_and_incomplete_incidence() {
    let (s, d) = super::super::tests::fixture(1);
    let c = Cutter::new(s, d).unwrap();
    let body = c.boundary_topology().unwrap();
    let circle = c.finite_section(c.envelope.tip.lo * 0.95).unwrap();
    assert!(MouthNodes::read(circle.vertices()).is_err());
    let mut duplicate = body.events().mouth().vertices().to_vec();
    duplicate.push(duplicate[0]);
    assert!(MouthNodes::read(&duplicate).is_err());
    let mut faces = body.faces.clone();
    faces.remove(0);
    let assembly = Assembly {
        events: body.events.clone(),
        vertices: body.vertices.clone(),
        edges: body.edges.clone(),
        faces,
    };
    assert!(assembly.finish().is_err());
    // The valid source geometry and previously published body remain unchanged.
    assert_eq!(body.certificate().euler_characteristic, 2);
    assert!(c.boundary_topology().is_ok());
}
