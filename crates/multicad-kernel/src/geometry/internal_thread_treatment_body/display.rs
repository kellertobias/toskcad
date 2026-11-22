//! Viewport curves follow the same immutable source edges as native references.
use super::interval::Interval as B;
use super::*;
const POINT_LIMIT: usize = 500_000;
impl InternalThreadTreatmentBody {
    pub(in crate::geometry) fn display_edges(
        &self,
        frame: PlaneFrame,
        deflection: f64,
        policy: GeometryPolicy,
    ) -> Result<Vec<DisplayEdge>, GeometryError> {
        self.display_edges_limited(frame, deflection, POINT_LIMIT, policy)
    }
    fn display_edges_limited(
        &self,
        frame: PlaneFrame,
        deflection: f64,
        point_limit: usize,
        policy: GeometryPolicy,
    ) -> Result<Vec<DisplayEdge>, GeometryError> {
        if !(2..=POINT_LIMIT).contains(&point_limit) {
            return Err(input("Invalid aggregate internal curve point budget"));
        }
        if !deflection.is_finite() || deflection <= 0. {
            return Err(input("Invalid internal curve display quality"));
        }
        self.validate_frame(frame, policy)?;
        let boundary = self.qualified_boundary()?;
        let (_, lower) = metric(frame)?;
        // lower² <=1-e, with e bounding ||AᵀA-I||. Thus sigma_max² <=2-lower².
        let upper = (B::point(2.) - B::point(lower).square()).sqrt().hi;
        let quality = (B::point(deflection) / (B::point(2.) * B::point(upper))).lo;
        if !quality.is_finite() || quality <= 0. {
            return Err(unsupported(
                "Placed curve quality is below numerical resolution",
            ));
        }
        // The source proof caps the table at8192 descriptors. Point payloads
        // have the separate aggregate500k budget, before each source query.
        let mut output = Vec::with_capacity(boundary.topology_counts().1);
        let axes = [frame.u, frame.v, frame.normal()];
        let mut remaining = point_limit;
        for (ordinal, (id, _)) in boundary.curves().enumerate() {
            if remaining < 2 {
                return Err(failure(source::Error::Limit));
            }
            let mut line = boundary
                .curve_polyline(id, quality, remaining, None)
                .map_err(failure)?;
            remaining -= line.points_mm.len();
            let mut rounding = 0_f64;
            for point in &mut line.points_mm {
                let mut error = 0_f64;
                let placed = std::array::from_fn(|i| {
                    let (v, e) = dot_error(axes.map(|axis| axis[i]), *point);
                    let (v, r) = sum_error(frame.origin[i], v);
                    error = up(error, up(e, r));
                    v
                });
                policy.point(placed)?;
                rounding = rounding.max(error);
                *point = placed;
            }
            if !rounding.is_finite()
                || up(product(upper, line.maximum_error_mm), rounding) > deflection
            {
                return Err(unsupported(
                    "Literal placed curve points exceed the requested chord error",
                ));
            }
            output.push(DisplayEdge {
                index: ordinal as u32,
                points_mm: line.points_mm,
            });
        }
        Ok(output)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use multicad_rust_kernel::Rational;
    #[test]
    fn actual_source_edges_and_vertices_match_native_topology_at_both_hands_and_poses() {
        for (left, starts) in [(false, 1), (true, 2)] {
            let policy = GeometryPolicy::default();
            let frame = PlaneFrame::default();
            let stock = super::super::tests::fixture_stock(Rational::new(4, 1).unwrap(), frame);
            let original = stock
                .rust_cut_authored_internal_thread_treatments_with_interrupt(
                    super::super::tests::cut(frame, left, starts),
                    EntityId::from_uuid_bits(8),
                    policy,
                    None,
                )
                .unwrap();
            for solid in [
                original.share().unwrap(),
                original
                    .transformed(RigidTransform {
                        rotation: [[0.8, -0.36, 0.48], [0.6, 0.48, -0.64], [0., 0.8, 0.6]],
                        translation_mm: [100., -200., 300.],
                    })
                    .unwrap(),
            ] {
                let display = solid
                    .topology_with_display(EntityId::from_uuid_bits(9), 0, Some(0.01))
                    .unwrap();
                let (body, placed) = solid.internal_thread_treatment_body().unwrap();
                let boundary = body.qualified_boundary().unwrap();
                assert_eq!(display.edges.len(), boundary.topology_counts().1);
                assert_eq!(display.vertices.len(), boundary.topology_counts().0);
                assert_eq!(solid.display_edges(0.01).unwrap(), display.edges);
                assert_eq!(solid.display_vertices().unwrap(), display.vertices);
                let n = placed.normal();
                for (edge, (_, curve)) in display.edges.iter().zip(boundary.curves()) {
                    assert!(edge.points_mm.len() >= 2);
                    let first = boundary.points()[curve.endpoints()[0] as usize].position_bounds();
                    let local = first.map(|b| b.lower.midpoint(b.upper));
                    for (i, p) in edge.points_mm[0].into_iter().enumerate() {
                        let expected = placed.origin[i]
                            + placed.u[i] * local[0]
                            + placed.v[i] * local[1]
                            + n[i] * local[2];
                        assert!((p - expected).abs() < policy.linear_tolerance_mm);
                    }
                }
                let first_count = display.edges[0].points_mm.len();
                assert!(
                    body.display_edges_limited(placed, 0.01, first_count, policy)
                        .is_err()
                );
                assert!(body.display_edges_limited(placed, 0.01, 1, policy).is_err());
                assert!(body.display_edges(placed, f64::NAN, policy).is_err());
            }
        }
    }
    #[test]
    fn actual_curve_queries_refuse_unsupported_domains_frames_and_numeric_quality() {
        let policy = GeometryPolicy::default();
        let frame = PlaneFrame::default();
        let stock = super::super::tests::fixture_stock(Rational::new(4, 1).unwrap(), frame);
        for depth in [0.5, 0.1] {
            let mut request = super::super::tests::cut(frame, false, 1);
            request.intent.entry_depth_mm = Some(depth);
            let solid = stock
                .rust_cut_authored_internal_thread_treatments_with_interrupt(
                    request,
                    EntityId::from_uuid_bits(8),
                    policy,
                    None,
                )
                .unwrap();
            let (body, _) = solid.internal_thread_treatment_body().unwrap();
            if depth == 0.1 {
                assert!(body.display_edges(frame, 0.01, policy).is_err());
                continue;
            }
            assert!(
                body.display_edges(frame, f64::from_bits(1), policy)
                    .is_err()
            );
            assert!(
                body.display_edges(
                    PlaneFrame {
                        u: [0.; 3],
                        ..frame
                    },
                    0.01,
                    policy
                )
                .is_err()
            );
            assert!(
                body.display_edges(
                    PlaneFrame {
                        origin: [policy.maximum_coordinate_mm, 0., 0.],
                        ..frame
                    },
                    0.01,
                    policy
                )
                .is_err()
            );
            let far = PlaneFrame {
                origin: [1_000_000.1; 3],
                ..frame
            };
            let error = body.display_edges(far, 1e-10, policy).unwrap_err();
            assert!(error.message.contains("Literal placed curve points"));
            assert!(body.display_edges(far, 1e-9, policy).is_err());
        }
    }
}
