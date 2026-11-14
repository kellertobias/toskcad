use super::*;
use num_rational::BigRational as Q;
use std::cell::Cell;
fn q(x: f64) -> Q {
    Q::from_float(x).unwrap()
}
fn verify(body: &RadialBoundary, traced: &TracedMesh<'_>) {
    assert!(std::ptr::eq(body, traced.source()));
    let rims = traced.cap_rims();
    assert_eq!(rims[0].len(), rims[1].len());
    let mut spline = 0;
    let mut circle = 0;
    let mut vertices = 0;
    for (end, rim) in rims.into_iter().enumerate() {
        let mut indices = BTreeSet::new();
        for point in rim {
            assert!(indices.insert(point.vertex));
            assert_eq!(
                point.cap_role,
                if end == 0 {
                    CapRole::AxialStart
                } else {
                    CapRole::AxialEnd
                }
            );
            assert_eq!(point.uv[1], end as f64);
            assert_eq!(point.side_face, body.side_faces()[point.side]);
            assert_eq!(
                point.side_surface,
                body.input().faces[point.side_face].surface
            );
            assert_eq!(
                point.cap_surface,
                body.input().faces[point.cap_face].surface
            );
            assert_eq!(point.cap_face, body.ordered.caps.axial.caps[end]);
            let literal = super::super::point(&traced.mesh, point.vertex as usize);
            for (k, bound) in point.source_position_mm.iter().enumerate() {
                assert_eq!(literal[k], I::new(bound[0], bound[1]).midpoint());
            }
            if let Some(vertex) = point.source_vertex {
                vertices += 1;
                assert_eq!(
                    point.source_position_mm,
                    body.input().vertices[vertex].map(|x| [x, x])
                );
                continue;
            }
            let group = &body.ordered.caps.rims[end][point.side];
            let t = if group.len() == 1 {
                q(point.uv[0])
            } else if point.uv[0] >= 0.5 {
                q(2.) * q(point.uv[0]) - q(1.)
            } else {
                q(2.) * q(point.uv[0])
            };
            assert!(q(point.edge_parameter[0]) <= t && q(point.edge_parameter[1]) >= t);
            assert!(group.contains(&(point.edge, point.edge_forward)));
            match &body.input().curves[body.input().edges[point.edge].curve] {
                Curve::Spline { spans, poles } => {
                    spline += 1;
                    let mut poles = poles.clone();
                    if body.input().edges[point.edge].forward != point.edge_forward {
                        poles.reverse();
                    }
                    let knots = (0..=*spans)
                        .map(|i| i as f64 / *spans as f64)
                        .collect::<Vec<_>>();
                    let span = (0..*spans)
                        .find(|&i| q(knots[i]) <= t && t <= q(knots[i + 1]))
                        .unwrap();
                    let local = (&t - q(knots[span])) / (q(knots[span + 1]) - q(knots[span]));
                    for k in 0..3 {
                        let mut p = poles[9 * span..=9 * span + 9]
                            .iter()
                            .map(|p| q(p[k]))
                            .collect::<Vec<_>>();
                        for n in (1..10).rev() {
                            for i in 0..n {
                                p[i] = (&q(1.) - &local) * &p[i] + &local * &p[i + 1];
                            }
                        }
                        let b = point.source_position_mm[k];
                        assert!(q(b[0]) <= p[0] && q(b[1]) >= p[0]);
                    }
                }
                Curve::Circle { axis, radius_mm } => {
                    circle += 1;
                    let p = literal;
                    // Independent actual support check, with only the source
                    // enclosure propagated; no caller deflection as residual.
                    let delta = std::array::from_fn::<_, 3, _>(|k| {
                        I::new(
                            point.source_position_mm[k][0],
                            point.source_position_mm[k][1],
                        ) - I::point(axis.origin[k])
                    });
                    let basis = checked_axis(*axis).unwrap();
                    let h = dot(delta, basis[2]);
                    assert!(h.lo <= 0. && h.hi >= 0.);
                    let r2 = dot(delta, basis[0]).square() + dot(delta, basis[1]).square();
                    assert!(r2.lo <= radius_mm * radius_mm && r2.hi >= radius_mm * radius_mm);
                    assert!(p.into_iter().all(f64::is_finite));
                }
            }
        }
    }
    assert!(vertices > 0);
    // Cam fixture has genuine two-arc sides, thread fixture has only splines.
    if body
        .ordered
        .caps
        .kinds
        .iter()
        .any(|k| matches!(k, caps::SideKind::Cone))
    {
        assert!(circle > 0);
    } else {
        assert!(spline > 0);
    }
}
#[test]
fn original_ring_trace_keeps_mesh_bytes_and_actual_cap_source_ownership() {
    let body = super::super::tests::cam();
    for frame in [
        None,
        Some([[10., 20., 30.], [0.8, 0.6, 0.], [-0.36, 0.48, 0.8]]),
    ] {
        let body = if let Some(frame) = frame {
            body.transformed_with_interrupt(frame, 1e-6, 1e-5, None)
                .unwrap()
        } else {
            body.clone()
        };
        let options = super::super::tests::options();
        let old = body.tessellate_with_interrupt(options, None).unwrap();
        let traced = body
            .tessellate_with_cap_trace(options, Default::default(), None)
            .unwrap();
        assert_eq!(old, *traced.mesh());
        verify(&body, &traced);
        assert!(traced.retained_metadata_bytes() < BYTES);
    }
}
#[test]
fn handed_nonbinary_threads_keep_original_trace_parameters_and_mesh() {
    for (starts, left_hand, pitch) in [(1, false, 1.25), (2, true, 1.3)] {
        let source = crate::tapered_thread::Thread::new(
            crate::helical_thread::Spec {
                major_radius_mm: 4.,
                pitch_mm: pitch,
                length_mm: pitch,
                depth_mm: 0.5,
                crest_width_mm: 0.25,
                lower_flank_deg: 12.,
                upper_flank_deg: 18.,
                starts,
                left_hand,
                root_radius_mm: None,
            },
            0.03125,
            1e-7,
            1e-6,
        )
        .unwrap();
        let body = source
            .retained_boundary_with_interrupt(
                [[0.; 3], [1., 0., 0.], [0., 1., 0.]],
                1e-7,
                1e-6,
                None,
            )
            .unwrap()
            .boundary;
        let options = Tessellation {
            deflection: 0.1,
            angle: 0.4,
            triangle_limit: 100_000,
        };
        let old = body.tessellate_with_interrupt(options, None).unwrap();
        let traced = body
            .tessellate_with_cap_trace(options, Default::default(), None)
            .unwrap();
        assert_eq!(old, *traced.mesh());
        verify(&body, &traced);
    }
}
#[test]
fn trace_metadata_work_and_cancel_refusals_leave_source_and_old_path_reusable() {
    let body = super::super::tests::cam();
    let options = super::super::tests::options();
    let calls = Cell::new(0);
    let poll = || {
        calls.set(calls.get() + 1);
        false
    };
    let old = body
        .tessellate_with_interrupt(options, Some(&poll))
        .unwrap();
    let base = calls.get();
    let trace = body
        .tessellate_with_cap_trace(options, Default::default(), None)
        .unwrap();
    let exact = trace.retained_metadata_bytes();
    let work = trace.work_used();
    assert!(
        body.tessellate_with_cap_trace_budget(options, Default::default(), work, None)
            .is_ok()
    );
    assert!(matches!(
        body.tessellate_with_cap_trace_budget(options, Default::default(), work - 1, None),
        Err(Error::Limit)
    ));
    assert!(matches!(
        body.tessellate_with_cap_trace_budget(
            options,
            Default::default(),
            DEFAULT_LIMITS.work + 1,
            None
        ),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        body.tessellate_with_cap_trace_budget(options, Default::default(), work, Some(&|| true)),
        Err(Error::Cancelled)
    ));

    let traced_calls = Cell::new(0);
    let traced_poll = || {
        traced_calls.set(traced_calls.get() + 1);
        false
    };
    assert!(
        body.tessellate_with_cap_trace(options, Default::default(), Some(&traced_poll))
            .is_ok()
    );
    let final_poll = traced_calls.get() - 1;

    assert_eq!(*trace.mesh(), old);
    assert!(
        body.tessellate_with_cap_trace(
            options,
            CapTraceOptions {
                maximum_bytes: exact
            },
            None
        )
        .is_ok()
    );
    for bytes in [0, exact - 1] {
        assert!(matches!(
            body.tessellate_with_cap_trace(
                options,
                CapTraceOptions {
                    maximum_bytes: bytes
                },
                None
            ),
            Err(Error::Limit)
        ));
    }
    assert!(matches!(
        body.tessellate_with_cap_trace(
            options,
            CapTraceOptions {
                maximum_bytes: BYTES + 1
            },
            None
        ),
        Err(Error::Invalid)
    ));
    for cutoff in [0, base, base + 1, base + 10, final_poll] {
        let calls = Cell::new(0);
        let cancel = || {
            let n = calls.get();
            calls.set(n + 1);
            n >= cutoff
        };
        assert!(matches!(
            body.tessellate_with_cap_trace(options, Default::default(), Some(&cancel)),
            Err(Error::Cancelled)
        ));
    }
    let base_work = Cell::new(0);
    build_with_output(
        &body,
        options,
        None,
        DEFAULT_LIMITS,
        |_, _, _, mesh, work| {
            base_work.set(DEFAULT_LIMITS.work - work.remaining);
            Ok(mesh)
        },
    )
    .unwrap();
    let limits = Limits {
        work: base_work.get(),
        ..DEFAULT_LIMITS
    };
    assert!(matches!(
        build_with_output(
            &body,
            options,
            None,
            limits,
            |prepared, grid, ring, mesh, work| finish(
                &body,
                prepared,
                grid,
                ring,
                mesh,
                Default::default(),
                (limits.work, work)
            )
        ),
        Err(Error::Limit)
    ));
    let calls = Cell::new(0);
    let poll = || {
        calls.set(calls.get() + 1);
        false
    };
    assert_eq!(
        body.tessellate_with_interrupt(options, Some(&poll))
            .unwrap(),
        old
    );
    assert_eq!(calls.get(), base);
}
