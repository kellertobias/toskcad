use super::*;
fn spec(quarter: bool) -> multicad_rust_kernel::axial_pipe_thread::authored::AuthoredSpec {
    use multicad_rust_kernel::axial_pipe_thread::authored::{AuthoredSection, AuthoredSpec};
    let (tpi, e1, l1, l2, od, max) = if quarter {
        (18., 0.49163, 0.2278, 0.4018, 0.540, 0.088)
    } else {
        (27., 0.37360, 0.1615, 0.2639, 0.405, 0.096)
    };
    let p = 25.4 / tpi;
    let k = 1. / 32.;
    let m = 3f64.sqrt();
    let s = 1. / (m - k) + 1. / (m + k);
    let d = p / s - (0.033 + max) * p;
    AuthoredSpec {
        major_radius_mm: (e1 * 25.4 - 2. * k * l1 * 25.4 + d) / 2.,
        pitch_mm: p,
        length_mm: l2 * 25.4,
        radial_taper: k,
        section: AuthoredSection::NptFlat {
            depth_mm: d,
            flat_width_mm: (p - d * s) / 2.,
        },
        crest_radius_limit_mm: Some(od * 25.4 / 2.),
    }
}
#[test]
fn retained_exact_stock_cut_has_one_actual_manufacturing_shell() {
    for quarter in [false, true] {
        let spec = spec(quarter);
        let policy = GeometryPolicy::default();
        let stock = crate::geometry::construction::cylinder(
            PlaneFrame::default(),
            if quarter { 6.858 } else { 5.1435 },
            spec.length_mm + 2.,
            EntityId::new(),
            policy,
        )
        .unwrap();
        let before = *stock.axial_cylinder().unwrap();
        let solid = stock
            .rust_cut_external_npt_stock(spec, EntityId::new(), PlaneFrame::default(), policy, None)
            .unwrap();
        assert_eq!(
            *solid.rust_external_npt_stock_source().unwrap().stock(),
            before
        );
        assert_eq!(*stock.axial_cylinder().unwrap(), before);
        assert!(solid.rust_clipped_npt_source().is_none());
        assert!(!solid.supports_step_export());
        assert_eq!(solid.faces().unwrap().len(), 9);
        let topology = solid.topology_snapshot(EntityId::new(), 1).unwrap();
        assert_eq!(topology.elements.len(), 51);
        assert!(
            solid
                .rust_external_npt_stock_contains([0., 0., spec.length_mm + 1.])
                .unwrap()
        );
        assert!(
            !solid
                .rust_external_npt_stock_contains([20., 0., spec.length_mm + 1.])
                .unwrap()
        );
        assert!(solid.volume_mm3().unwrap() < stock.volume_mm3().unwrap());
        let receipt = solid
            .manufacturing_mesh_with_receipt(MeshOptions {
                linear_deflection_mm: 0.02,
                angular_deflection_rad: 0.2,
                maximum_triangles: 200000,
            })
            .unwrap();
        assert_eq!(receipt.mesh.face_indices.iter().copied().max(), Some(8));
        assert_eq!(receipt.maximum_displacement_mm, 0.);
        assert!(
            solid
                .rust_external_npt_stock_transfer()
                .unwrap()
                .position_error_mm
                < policy.linear_tolerance_mm
        );
    }
}
#[test]
fn exact_stock_cut_rejects_actual_bad_stock_and_cancel() {
    let s = spec(false);
    let p = GeometryPolicy::default();
    let stock = crate::geometry::construction::cylinder(
        PlaneFrame::default(),
        5.1435,
        s.length_mm + 2.,
        EntityId::new(),
        p,
    )
    .unwrap();
    assert!(
        stock
            .rust_cut_external_npt_stock(
                s,
                EntityId::new(),
                PlaneFrame::default(),
                p,
                Some(&|| true)
            )
            .is_err()
    );
    let mut frame = PlaneFrame::default();
    frame.origin[0] = 1.;
    assert!(
        stock
            .rust_cut_external_npt_stock(s, EntityId::new(), frame, p, None)
            .is_err()
    );
    let short = crate::geometry::construction::cylinder(
        PlaneFrame::default(),
        5.1435,
        1.,
        EntityId::new(),
        p,
    )
    .unwrap();
    assert!(
        short
            .rust_cut_external_npt_stock(s, EntityId::new(), PlaneFrame::default(), p, None)
            .is_err()
    );
}
#[test]
fn actual_source_edges_and_placed_literal_manufacturing_are_qualified() {
    let s = spec(false);
    let p = GeometryPolicy::default();
    let stock = crate::geometry::construction::cylinder(
        PlaneFrame::default(),
        5.1435,
        s.length_mm + 2.,
        EntityId::new(),
        p,
    )
    .unwrap();
    let solid = stock
        .rust_cut_external_npt_stock(s, EntityId::new(), PlaneFrame::default(), p, None)
        .unwrap();
    let q = multicad_rust_kernel::display::Tessellation {
        deflection: 0.02,
        angle: 0.2,
        triangle_limit: 200000,
    };
    for i in 0..24 {
        let points = solid
            .sample_non_analytic_edge(i, q)
            .unwrap()
            .unwrap_or_else(|e| panic!("edge {i}: {e:?}"));
        assert!(points.len() >= 2);
    }
    let placed = solid
        .transformed(RigidTransform {
            rotation: [[0., 0., 1.], [0., 1., 0.], [-1., 0., 0.]],
            translation_mm: [3., -2., 7.],
        })
        .unwrap();
    let receipt = placed
        .manufacturing_mesh_with_receipt(MeshOptions {
            linear_deflection_mm: 0.02,
            angular_deflection_rad: 0.2,
            maximum_triangles: 200000,
        })
        .unwrap();
    assert!(receipt.vertex_links > 0);
    assert!(
        placed
            .rust_external_npt_stock_contains([3. + s.length_mm + 1., -2., 7.])
            .unwrap()
    );
}

#[test]
fn original_world_exact_and_local_posed_cylinder_charts_keep_their_actual_datums() {
    let s = spec(false);
    let p = GeometryPolicy::default();
    let frames = [
        PlaneFrame {
            origin: [1., -2., 0.5],
            ..Default::default()
        },
        PlaneFrame {
            origin: [1., -2., 0.5],
            u: [0., 0., -1.],
            v: [0., 1., 0.],
        },
    ];
    for frame in frames {
        let stock = crate::geometry::construction::cylinder(
            frame,
            5.1435,
            s.length_mm + 2.,
            EntityId::new(),
            p,
        )
        .unwrap();
        let exact = *stock.axial_cylinder_frame().unwrap().0;
        let solid = stock
            .rust_cut_external_npt_stock(s, EntityId::new(), frame, p, None)
            .unwrap();
        assert_eq!(
            *solid.rust_external_npt_stock_source().unwrap().stock(),
            exact
        );
        let n = frame.normal();
        let point = std::array::from_fn(|i| frame.origin[i] + n[i] * (s.length_mm + 1.));
        assert!(solid.rust_external_npt_stock_contains(point).unwrap());
        solid
            .manufacturing_mesh(MeshOptions {
                linear_deflection_mm: 0.02,
                angular_deflection_rad: 0.2,
                maximum_triangles: 200000,
            })
            .unwrap();
    }
}
#[test]
fn actual_face_measure_cancellation_keeps_its_public_kind_and_context() {
    use std::cell::Cell;
    let s = spec(false);
    let p = GeometryPolicy::default();
    let stock = crate::geometry::construction::cylinder(
        PlaneFrame::default(),
        5.1435,
        s.length_mm + 2.,
        EntityId::new(),
        p,
    )
    .unwrap();
    let count = Cell::new(0usize);
    stock
        .rust_cut_external_npt_stock(
            s,
            EntityId::new(),
            PlaneFrame::default(),
            p,
            Some(&|| {
                count.set(count.get() + 1);
                false
            }),
        )
        .unwrap();
    let total = count.get();
    let mut found = false;
    // Actual interruption frontiers, across the real production call graph.
    // The bounded stride finds the comparatively long source measure stage.
    for cutoff in (1..total).step_by((total / 64).max(1)) {
        let calls = Cell::new(0usize);
        let error = stock
            .rust_cut_external_npt_stock(
                s,
                EntityId::new(),
                PlaneFrame::default(),
                p,
                Some(&|| {
                    calls.set(calls.get() + 1);
                    calls.get() >= cutoff
                }),
            )
            .err()
            .unwrap();
        assert_eq!(error.kind, GeometryErrorKind::Cancelled);
        if error.message.starts_with("Face ") {
            found = true;
            break;
        }
    }
    assert!(
        found,
        "no actual face measure interruption among {total} checkpoints"
    );
}
#[test]
fn actual_public_topology_and_input_budget_frontiers() {
    let s = spec(false);
    let p = GeometryPolicy::default();
    let feature = EntityId::new();
    let stock = crate::geometry::construction::cylinder(
        PlaneFrame::default(),
        5.1435,
        s.length_mm + 2.,
        feature,
        p,
    )
    .unwrap();
    assert!(
        stock
            .rust_cut_external_npt_stock(s, feature, PlaneFrame::default(), p, None)
            .is_err()
    );
    let solid = stock
        .rust_cut_external_npt_stock(s, EntityId::new(), PlaneFrame::default(), p, None)
        .unwrap();
    let body = solid.external_npt_stock_body().unwrap();
    for i in 0..9 {
        let rings = body.planar_loops(i);
        assert_eq!(
            rings.is_empty(),
            !matches!(
                body.boundary.faces()[i as usize].face(),
                SourceFace::Thread(StockClipFace::StartCap)
                    | SourceFace::Shoulder
                    | SourceFace::TailCap
            )
        );
    }
    assert!(body.planar_loops(99).is_empty());
    let q = multicad_rust_kernel::display::Tessellation {
        deflection: 0.02,
        angle: 0.2,
        triangle_limit: 200000,
    };
    assert!(body.sample_edge(99, q).is_none());
    let mut bad = q;
    bad.deflection = 0.;
    assert!(body.sample_edge(0, bad).unwrap().is_err());
    bad = q;
    bad.triangle_limit = 1;
    assert!(body.sample_edge(0, bad).unwrap().is_err());
    assert!(body.mesh(bad).is_err());
    assert!(stock.rust_external_npt_stock_contains([0.; 3]).is_err());
    let wide = crate::geometry::construction::cylinder(
        PlaneFrame::default(),
        5.15,
        s.length_mm + 2.,
        EntityId::new(),
        p,
    )
    .unwrap();
    assert!(
        wide.rust_cut_external_npt_stock(s, EntityId::new(), PlaneFrame::default(), p, None)
            .is_err()
    );
    assert!(
        solid
            .share()
            .unwrap()
            .rust_external_npt_stock_source()
            .is_some()
    );
}
#[test]
fn genuine_unsupported_chart_boundary_and_curve_frontiers() {
    let s = spec(false);
    let p = GeometryPolicy::default();
    let feature = EntityId::new();
    let frame = PlaneFrame {
        u: [0., 1., 0.],
        v: [-1., 0., 0.],
        ..Default::default()
    };
    let rotated_phase =
        crate::geometry::construction::cylinder(frame, 5.1435, s.length_mm + 2., feature, p)
            .unwrap();
    assert!(
        rotated_phase
            .rust_cut_external_npt_stock(s, EntityId::new(), frame, p, None)
            .is_err()
    );
    let box_stock = Solid::box_solid([5.; 3], p).unwrap();
    assert!(
        box_stock
            .rust_cut_external_npt_stock(s, EntityId::new(), PlaneFrame::default(), p, None)
            .is_err()
    );
    let stock = crate::geometry::construction::cylinder(
        PlaneFrame::default(),
        5.1435,
        s.length_mm + 2.,
        feature,
        p,
    )
    .unwrap();
    let solid = stock
        .rust_cut_external_npt_stock(s, EntityId::new(), PlaneFrame::default(), p, None)
        .unwrap();
    assert!(
        solid
            .rust_external_npt_stock_contains([5.1435, 0., s.length_mm + 1.])
            .is_err()
    );
    let body = solid.external_npt_stock_body().unwrap();
    let q = multicad_rust_kernel::display::Tessellation {
        deflection: 1e-30,
        angle: 0.2,
        triangle_limit: 100,
    };
    assert!(body.sample_edge(0, q).unwrap().is_err());
    let q = multicad_rust_kernel::display::Tessellation {
        deflection: 0.02,
        angle: 0.2,
        triangle_limit: 1,
    };
    assert!(body.sample_edge(17, q).unwrap().is_err());
}
#[test]
fn rounded_affine_composition_accumulates_real_correspondence_and_eventually_refuses() {
    let s = spec(false);
    let p = GeometryPolicy::default();
    let stock = crate::geometry::construction::cylinder(
        PlaneFrame::default(),
        5.1435,
        s.length_mm + 2.,
        EntityId::new(),
        p,
    )
    .unwrap();
    let source = stock
        .rust_cut_external_npt_stock(s, EntityId::new(), PlaneFrame::default(), p, None)
        .unwrap();
    let shift = |x| RigidTransform {
        translation_mm: [x, 0., 0.],
        ..Default::default()
    };
    let first = source.transformed(shift(0.1)).unwrap();
    let second = first.transformed(shift(0.2)).unwrap();
    let before = first.rust_external_npt_stock_transfer().unwrap();
    let after = second.rust_external_npt_stock_transfer().unwrap();
    // Exact dyadic 0.1+0.2 differs from the retained rounded sum by2^-55.
    assert!(after.position_error_mm >= before.position_error_mm + 2f64.powi(-55));
    assert!(after.derivative_error_mm > before.derivative_error_mm);
    assert!(after.axis_alignment_mm > before.axis_alignment_mm);
    let mut far = source.transformed(shift(1_000_000.)).unwrap();
    let mut refused = false;
    for _ in 0..1024 {
        match far.transformed(shift(0.1)) {
            Ok(next) => {
                assert!(
                    next.rust_external_npt_stock_transfer()
                        .unwrap()
                        .position_error_mm
                        <= p.linear_tolerance_mm
                );
                far = next
            }
            Err(e) => {
                assert_eq!(e.kind, GeometryErrorKind::Unsupported);
                refused = true;
                break;
            }
        }
    }
    assert!(
        refused,
        "repeated real far-world rounding must not receive an understated transfer"
    );
}
