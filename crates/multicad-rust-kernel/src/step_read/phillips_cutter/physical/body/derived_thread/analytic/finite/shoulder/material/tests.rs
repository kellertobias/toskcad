use super::*;
#[test]
fn actual_m2_distinct_periodic_material_retains_shaped_cavity_and_complete_boundary() {
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
        let start = std::time::Instant::now();
        let material = head
            .periodic_closed_material(&join, Default::default(), None)
            .unwrap();
        assert!(std::ptr::eq(material.shoulder(), &join));
        assert!(std::ptr::eq(material.containment().head(), head.source()));
        assert_eq!(material.chart(), head.chart());
        assert!(material.determinant_bounds()[0] > 0.);
        assert_eq!(material.inherited_work_used(), join.total_work_used());
        assert_eq!(
            material.total_work_used(),
            material.inherited_work_used() + material.new_work_used()
        );
        assert!(material.retained_budget_bytes() < BYTES);
        let faces = material.boundary_faces().collect::<Vec<_>>();
        assert_eq!(
            faces
                .iter()
                .filter(|f| matches!(f, PeriodicMaterialFace::ThreadBand(_)))
                .count(),
            4
        );
        assert_eq!(
            faces
                .iter()
                .filter(|f| matches!(f, PeriodicMaterialFace::ShoulderBand(_)))
                .count(),
            3
        );
        assert!(
            faces.contains(&PeriodicMaterialFace::TopAnnulus)
                && faces.contains(&PeriodicMaterialFace::TipCap)
        );
        assert!(!faces.contains(&PeriodicMaterialFace::ReversedCavity(CavityFace::Mouth)));
        for domain in head
            .common_mouth()
            .closed_cavity()
            .local_source_certificate()
            .patch_domains()
        {
            if domain.source_face != CavityFace::Mouth {
                assert!(faces.contains(&PeriodicMaterialFace::ReversedCavity(domain.source_face)));
            }
        }
        eprintln!(
            "PERIODIC MATERIAL {:?} new{} total{} retained{} faces{} radial{} floor{}",
            start.elapsed(),
            material.new_work_used(),
            material.total_work_used(),
            material.retained_budget_bytes(),
            faces.len(),
            material.containment().minimum_radial_clearance_mm(),
            material.containment().axial_floor_clearance_mm()
        );
    });
}
#[test]
fn actual_periodic_material_public_frontiers_preserve_stage_sum_and_own_head_cavity() {
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
        let o = PeriodicMaterialOptions::default();
        let good = head.periodic_closed_material(&join, o, None).unwrap();
        let bytes = good.retained_budget_bytes()
            + RepairedHeadContainment::scratch_budget_bytes(o.maximum_cells);
        for bad in [
            PeriodicMaterialOptions {
                minimum_clearance_mm: 0.,
                ..o
            },
            PeriodicMaterialOptions {
                minimum_clearance_mm: f64::NAN,
                ..o
            },
            PeriodicMaterialOptions {
                maximum_cells: 0,
                ..o
            },
            PeriodicMaterialOptions {
                maximum_cells: CELLS + 1,
                ..o
            },
            PeriodicMaterialOptions {
                maximum_work: WORK + 1,
                ..o
            },
            PeriodicMaterialOptions {
                maximum_total_work: TOTAL + 1,
                ..o
            },
            PeriodicMaterialOptions {
                maximum_bytes: BYTES + 1,
                ..o
            },
        ] {
            assert!(matches!(
                head.periodic_closed_material(&join, bad, None),
                Err(PhillipsReadError::Invalid)
            ));
        }
        for bad in [
            PeriodicMaterialOptions {
                maximum_work: 0,
                ..o
            },
            PeriodicMaterialOptions {
                maximum_work: 1,
                ..o
            },
            PeriodicMaterialOptions {
                maximum_work: METADATA,
                ..o
            },
            PeriodicMaterialOptions {
                maximum_work: METADATA + 1,
                ..o
            },
            PeriodicMaterialOptions {
                maximum_work: good.new_work_used() - 1,
                ..o
            },
            PeriodicMaterialOptions {
                maximum_total_work: join.total_work_used() - 1,
                ..o
            },
            PeriodicMaterialOptions {
                maximum_total_work: join.total_work_used(),
                ..o
            },
            PeriodicMaterialOptions {
                maximum_total_work: good.total_work_used() - 1,
                ..o
            },
            PeriodicMaterialOptions {
                maximum_bytes: bytes - 1,
                ..o
            },
            PeriodicMaterialOptions {
                maximum_cells: 1,
                ..o
            },
        ] {
            assert!(matches!(
                head.periodic_closed_material(&join, bad, None),
                Err(PhillipsReadError::Limit)
            ));
        }
        assert!(matches!(
            head.periodic_closed_material(
                &join,
                PeriodicMaterialOptions {
                    minimum_clearance_mm: 1.,
                    ..o
                },
                None
            ),
            Err(PhillipsReadError::Unsupported)
        ));
        assert!(
            head.periodic_closed_material(
                &join,
                PeriodicMaterialOptions {
                    maximum_work: good.new_work_used(),
                    maximum_total_work: good.total_work_used(),
                    maximum_bytes: bytes,
                    ..o
                },
                None
            )
            .is_ok()
        );
        degree11::tests::with_head(&mut |other| {
            assert!(matches!(
                other.periodic_closed_material(&join, o, None),
                Err(PhillipsReadError::Invalid)
            ));
        });
        let calls = Cell::new(0);
        let count = || {
            calls.set(calls.get() + 1);
            false
        };
        head.periodic_closed_material(&join, o, Some(&count))
            .unwrap();
        let complete = calls.get();
        assert!(complete > 10);
        for stop in 1..=complete {
            calls.set(0);
            let cancel = || {
                calls.set(calls.get() + 1);
                calls.get() == stop
            };
            assert!(
                matches!(
                    head.periodic_closed_material(&join, o, Some(&cancel)),
                    Err(PhillipsReadError::Cancelled)
                ),
                "callback{stop}/{complete}"
            );
        }
        assert!(head.periodic_closed_material(&join, o, None).is_ok());
    });
}
