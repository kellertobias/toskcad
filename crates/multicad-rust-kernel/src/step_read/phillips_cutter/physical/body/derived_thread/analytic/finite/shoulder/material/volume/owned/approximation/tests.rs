use super::*;
#[test]
fn actual_authored_m2_proposal_is_only_a_bounded_approximation_of_the_owned_material() {
    for left in [false, true] {
        let owner = super::super::query::tests::owned(left, 0.03);
        let receipt = owner.authored_mesh_approximation(WORK, None).unwrap();
        assert!(std::ptr::eq(receipt.source(), &owner));
        assert!(
            receipt.maximum_position_error_mm() > 0. && receipt.maximum_position_error_mm() < 1e-8
        );
        assert!(
            receipt.maximum_normal_angle_radians() > 0.
                && receipt.maximum_normal_angle_radians() < 1e-7
        );
        assert!(
            receipt
                .maximum_thread_derivative_error_mm()
                .into_iter()
                .all(|v| v > 0. && v < 1e-8)
        );
        assert_eq!(receipt.proposal().thread().spec().left_hand, left);
        assert!(receipt.work_used() > 2_000_000 && receipt.work_used() < 2_100_000);
        assert!(
            owner
                .authored_mesh_approximation(receipt.work_used(), None)
                .is_ok()
        );
        assert!(matches!(
            owner.authored_mesh_approximation(receipt.work_used() - 1, None),
            Err(PhillipsReadError::Limit)
        ));
        let coefficients = receipt.periodic_outward_normal_coefficients();
        assert!(
            coefficients
                .into_iter()
                .flatten()
                .flatten()
                .all(f64::is_finite)
        );
        assert!(coefficients[0][1][0] > 0. && coefficients[2][1][1] < 0.);
        let head_limit = 1 + owner.head().view().segment_count() * 4096;
        assert!(matches!(
            owner.authored_mesh_approximation(head_limit, None),
            Err(PhillipsReadError::Limit)
        ));
        assert!(matches!(
            owner.authored_mesh_approximation(WORK + 1, None),
            Err(PhillipsReadError::Invalid)
        ));
        assert!(matches!(
            owner.authored_mesh_approximation(0, None),
            Err(PhillipsReadError::Limit)
        ));
        assert!(matches!(
            owner.authored_mesh_approximation(WORK, Some(&|| true)),
            Err(PhillipsReadError::Cancelled)
        ));
        let mesh = crate::helical_thread_mesh::tessellate_phillips_screw(
            receipt.proposal(),
            crate::display::Tessellation {
                deflection: 0.05,
                angle: 0.3,
                triangle_limit: 1_000_000,
            },
            &|| false,
        )
        .unwrap();
        let literal = receipt
            .certify_literal(&mesh, 0.1, 512_000_000, 64 * 1024 * 1024, None)
            .unwrap();
        assert!(std::ptr::eq(literal.source(), &owner) && std::ptr::eq(literal.candidate(), &mesh));
        assert!(
            literal.maximum_surface_error_mm() > 0. && literal.maximum_surface_error_mm() <= 0.1
        );
        assert!(
            literal.maximum_outer_vertex_error_mm() > 0.
                && literal.maximum_outer_vertex_error_mm() < 1e-7
        );
        assert!(literal.work_used() > 0 && literal.accounted_peak_bytes() <= 64 * 1024 * 1024);
        eprintln!(
            "DIRECTED literal point error{} whole{} work{} peak{}",
            literal.maximum_outer_vertex_error_mm(),
            literal.maximum_surface_error_mm(),
            literal.work_used(),
            literal.accounted_peak_bytes()
        );
        if !left {
            assert!(
                receipt
                    .certify_literal(
                        &mesh,
                        0.1,
                        literal.work_used(),
                        literal.accounted_peak_bytes(),
                        None
                    )
                    .is_ok()
            );
            assert!(matches!(
                receipt.certify_literal(
                    &mesh,
                    0.1,
                    literal.work_used() - 1,
                    64 * 1024 * 1024,
                    None
                ),
                Err(PhillipsReadError::Limit)
            ));
            assert!(matches!(
                receipt.certify_literal(
                    &mesh,
                    0.1,
                    512_000_000,
                    literal.accounted_peak_bytes() - 1,
                    None
                ),
                Err(PhillipsReadError::Limit)
            ));
            for tolerance in [0., f64::NAN] {
                assert!(matches!(
                    receipt.certify_literal(&mesh, tolerance, 512_000_000, 64 * 1024 * 1024, None),
                    Err(PhillipsReadError::Invalid)
                ));
            }
            assert!(matches!(
                receipt.certify_literal(&mesh, 0.001, 512_000_000, 64 * 1024 * 1024, None),
                Err(PhillipsReadError::Unsupported)
            ));
            assert!(matches!(
                receipt.certify_literal(&mesh, 0.1, 512_000_001, 64 * 1024 * 1024, None),
                Err(PhillipsReadError::Invalid)
            ));
            assert!(matches!(
                receipt.certify_literal(&mesh, 0.1, 512_000_000, 64 * 1024 * 1024 + 1, None),
                Err(PhillipsReadError::Invalid)
            ));
            for at in [0, 1, 20, 100, 1000] {
                let calls = Cell::new(0);
                let cancel = || {
                    let n = calls.get();
                    calls.set(n + 1);
                    n == at
                };
                assert!(matches!(
                    receipt.certify_literal(
                        &mesh,
                        0.1,
                        512_000_000,
                        64 * 1024 * 1024,
                        Some(&cancel)
                    ),
                    Err(PhillipsReadError::Cancelled)
                ));
            }
        }
        eprintln!(
            "Periodic approximation hand{left}: C0{} C1{:?} normal{} facets{} oldmass{:?} newmass{:?}",
            receipt.maximum_position_error_mm(),
            receipt.maximum_thread_derivative_error_mm(),
            receipt.maximum_normal_angle_radians(),
            mesh.mesh().indices.len() / 3,
            mesh.source_volume().volume_mm3(),
            owner.local_volume_mm3()
        );
        assert!(mesh.connectivity().signed_volume_mm3 > 0.);
    }
}
#[test]
fn approximation_component_cancellation_uses_one_public_balance() {
    let owner = super::super::query::tests::owned(false, 0.03);
    let calls = Cell::new(0);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    owner
        .authored_mesh_approximation(WORK, Some(&count))
        .unwrap();
    let total = calls.get();
    for at in (0..32).chain([total / 2, total - 2, total - 1]) {
        calls.set(0);
        let cancel = || {
            let n = calls.get();
            calls.set(n + 1);
            n == at
        };
        assert!(matches!(
            owner.authored_mesh_approximation(WORK, Some(&cancel)),
            Err(PhillipsReadError::Cancelled)
        ));
    }
}
#[test]
fn actual_decoded_owner_does_not_fabricate_authored_provenance_for_a_mesh_proposal() {
    degree11::tests::with_head(&mut |head| {
        let field = head
            .common_chart_periodic_thread_repair(Default::default(), None)
            .unwrap();
        let finite = head
            .periodic_thread_finite_curves(&field, Default::default(), None)
            .unwrap();
        let join = head
            .periodic_shared_crest_shoulder(&finite, Default::default(), None)
            .unwrap();
        let material = head
            .periodic_closed_material(&join, Default::default(), None)
            .unwrap();
        let mass = material.volume_bounds(Default::default(), None).unwrap();
        let owner = mass.owned_source(Default::default(), None).unwrap();
        assert!(matches!(
            owner.authored_mesh_approximation(WORK, None),
            Err(PhillipsReadError::Unsupported)
        ));
    });
}
