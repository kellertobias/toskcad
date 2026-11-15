use super::*;
#[test]
fn actual_retained_npt_native_occurrences_earn_whole_source_quality_and_literal_embedding() {
    for (quarter, originally_mapped) in [(false, false), (true, false), (false, true)] {
        let source = super::super::tests::source(quarter, originally_mapped);
        let plan = source.geometry_plan(Default::default(), None).unwrap();
        let original = source.world_volume_mm3();
        let p = LiteralOccurrencePlacement::new(
            [[0.8, -0.36, 0.48], [0.6, 0.48, -0.64], [0., 0.8, 0.6]],
            [17., -8., 13.],
            Default::default(),
            None,
        )
        .unwrap();
        let p = p
            .transformed(
                [[0., 0., 1.], [0., 1., 0.], [-1., 0., 0.]],
                [2., 3., 4.],
                Default::default(),
                None,
            )
            .unwrap();
        let started = std::time::Instant::now();
        let result = plan
            .placed_mesh(p.clone(), Default::default(), None)
            .unwrap();
        eprintln!(
            "literal occurrence quarter={quarter} original_mapped={originally_mapped} wall={:?} vertices={} facets={} C0euclidean={} normalchord={} work={} actual={} peak={}",
            started.elapsed(),
            result.mesh().vertices.len() / 3,
            result.mesh().face_indices.len(),
            result.position_error_mm(),
            result.normal_chord(),
            result.work_used(),
            result.query_work_used(),
            result.accounted_peak_bytes()
        );
        assert!(std::sync::Arc::ptr_eq(&result.source().0, &source.0));
        assert_eq!(result.placement().transforms().len(), 2);
        assert!(
            result.position_error_mm() <= 0.1 && result.normal_chord() <= 2. * (0.3f64 / 2.).sin()
        );
        assert_eq!(
            crate::manufacturing::certify(result.mesh()).unwrap(),
            result.embedding().connectivity()
        );
        assert!(result.embedding().minimum_radial_determinant() > 0.);
        assert!(result.query_work_used() < result.work_used());
        assert_eq!(result.work_used() + result.remaining_work(), WORK);
        assert!(result.accounted_peak_bytes() <= BYTES);
        let det = p.determinant_bounds();
        let volume = interval(original) * interval(det);
        assert!(volume.valid() && volume.lo > 0.);
        assert!(!result.into_mesh().indices.is_empty());
    }
}
