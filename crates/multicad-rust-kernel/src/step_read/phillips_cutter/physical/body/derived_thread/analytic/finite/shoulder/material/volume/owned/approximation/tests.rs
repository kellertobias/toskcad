use super::*;
#[test]
fn actual_authored_m2_proposal_is_only_a_bounded_approximation_of_the_owned_material() {
    for left in [false, true] {
        let owner = super::super::query::tests::owned(left, 0.03);
        let receipt = owner.authored_mesh_approximation(WORK, None).unwrap();
        assert!(std::ptr::eq(receipt.source(), &owner));
        assert_eq!(
            receipt.proposal_origin(),
            PeriodicProposalOrigin::AuthoredSeed
        );
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
        let receipt = owner.mesh_approximation(WORK, None).unwrap();
        assert_eq!(
            receipt.proposal_origin(),
            PeriodicProposalOrigin::DecodedSeed
        );
        assert!(std::ptr::eq(receipt.source(), &owner));
        assert!(matches!(
            owner.provenance(),
            PeriodicMaterialProvenance::Decoded { .. }
        ));
        assert!(receipt.maximum_position_error_mm() < 1e-7);
        assert!(owner.mesh_approximation(receipt.work_used(), None).is_ok());
        assert!(matches!(
            owner.mesh_approximation(receipt.work_used() - 1, None),
            Err(PhillipsReadError::Limit)
        ));
        assert!(matches!(
            owner.mesh_approximation(0, None),
            Err(PhillipsReadError::Limit)
        ));
        assert!(matches!(
            owner.mesh_approximation(WORK + 1, None),
            Err(PhillipsReadError::Invalid)
        ));
        assert!(matches!(
            owner.mesh_approximation(WORK, Some(&|| true)),
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
        assert!(std::ptr::eq(literal.source(), &owner));
        assert!(literal.maximum_outer_vertex_error_mm() < 1e-7);
        eprintln!(
            "decoded proposal literal C0 {} peak {} work {}",
            literal.maximum_surface_error_mm(),
            literal.accounted_peak_bytes(),
            literal.work_used()
        );
    });
}

#[test]
fn actual_periodic_m2_face_measures_transfer_whole_parameter_domains() {
    for left in [false, true] {
        let owner = super::super::query::tests::owned(left, 0.03);
        let approximation = owner.authored_mesh_approximation(WORK, None).unwrap();
        let now = std::time::Instant::now();
        let properties = approximation
            .source_face_properties(
                crate::phillips_drive::FacePropertyOptions {
                    area_tolerance_mm2: 1e-5,
                    first_moment_tolerance_mm3: 1e-5,
                    centroid_tolerance_mm: 1e-5,
                    maximum_cells: 65536,
                },
                WORK,
                None,
            )
            .unwrap();
        eprintln!(
            "actual periodic measures hand{left} elapsed {:?}",
            now.elapsed()
        );
        assert_eq!(properties.len(), owner.face_count());
        let mut work = Work {
            remaining: WORK,
            interrupt: None,
            cancelled: Cell::new(false),
        };
        let rotation = rotation::Rotation::new(&mut work).unwrap();
        let pi = rotation.pi(&mut work).unwrap();
        let phases = owner.thread().phases();
        let radii = owner.thread().radius_nodes_mm();
        let pitch = I::point(owner.thread().pitch_mm());
        let caps = owner.cap_expressions_mm().map(expression);
        let height = caps[1] - caps[0];
        let mut disk = I::point(0.);
        for band in 0..4 {
            let a = I::point(radii[band]);
            let b = I::point(radii[band + 1]);
            disk = disk
                + pi * (phases[band + 1] - phases[band]) * (a.square() + a * b + b.square())
                    / (I::point(3.) * pitch);
        }
        let crest =
            I::point(2.) * pi * I::point(radii[3]) * (phases[4] - phases[3]) * height / pitch;
        let h = owner.head().view().segment_count();
        for (index, oracle) in [
            (h + 3, crest),
            (h + 4, disk),
            (h + 5, pi * I::point(radii[0]).square() - disk),
        ] {
            let enclosure = properties[index].area_mm2;
            assert!(
                enclosure[0] <= oracle.lo && enclosure[1] >= oracle.hi,
                "source area{index} {:?} vs {:?}",
                enclosure,
                bounds_for_test(oracle)
            );
        }
        assert!(
            properties[h + 4].centroid_mm[2][0] <= caps[1].lo
                && properties[h + 4].centroid_mm[2][1] >= caps[1].hi
        );

        for (i, p) in properties.iter().enumerate() {
            assert_eq!(p.source_face, i);
            assert!(p.area_mm2[0] > 0. && p.area_mm2[1] >= p.area_mm2[0]);
            assert!(p.area_transfer_error_mm2.is_finite());
            assert!(
                p.centroid_mm
                    .iter()
                    .all(|b| b[0].is_finite() && b[1] >= b[0] && b[1] - b[0] <= 1e-5)
            );
        }
    }
}

fn bounds_for_test(v: I) -> [f64; 2] {
    [v.lo, v.hi]
}
#[test]
fn actual_source_face_measure_frontiers_refuse_without_partial_publication() {
    let owner = super::super::query::tests::owned(false, 0.03);
    let approximation = owner.authored_mesh_approximation(WORK, None).unwrap();
    let options = crate::phillips_drive::FacePropertyOptions::default();
    for (maximum, expected) in [
        (0, PhillipsReadError::Limit),
        (1, PhillipsReadError::Limit),
        (WORK + 1, PhillipsReadError::Invalid),
    ] {
        assert!(
            matches!(approximation.source_face_properties(options,maximum,None),Err(e) if e==expected)
        );
    }
    assert!(matches!(
        approximation.source_face_properties(
            crate::phillips_drive::FacePropertyOptions {
                area_tolerance_mm2: 0.,
                ..options
            },
            WORK,
            None
        ),
        Err(PhillipsReadError::Invalid)
    ));
    assert!(matches!(
        approximation.source_face_properties(options, WORK, Some(&|| true)),
        Err(PhillipsReadError::Cancelled)
    ));
    let calls = Cell::new(0);
    let cancel = || {
        let n = calls.get() + 1;
        calls.set(n);
        n > 10
    };
    assert!(matches!(
        approximation.source_face_properties(options, WORK, Some(&cancel)),
        Err(PhillipsReadError::Cancelled)
    ));
    assert!(calls.get() > 10);
}
