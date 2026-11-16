use super::*;
use num_rational::BigRational as Q;
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
        assert_eq!(
            plan.placed_mesh(
                p.clone(),
                LiteralMeshOptions {
                    maximum_work: 1,
                    ..Default::default()
                },
                None,
            )
            .err(),
            Some(Error::Limit)
        );
        assert_eq!(
            plan.placed_mesh(p.clone(), Default::default(), Some(&|| true))
                .err(),
            Some(Error::Cancelled)
        );
        let source_bounds = source.world_enclosure_mm();
        let mut corners = Vec::new();
        for mask in 0..8 {
            let mut point = std::array::from_fn::<_, 3, _>(|axis| {
                Q::from_float(source_bounds[axis][(mask >> axis) & 1]).unwrap()
            });
            for transform in p.transforms() {
                let matrix = transform.matrix();
                let translation = transform.translation_mm();
                point = std::array::from_fn(|i| {
                    (0..3)
                        .map(|j| Q::from_float(matrix[i][j]).unwrap() * point[j].clone())
                        .sum::<Q>()
                        + Q::from_float(translation[i]).unwrap()
                });
            }
            corners.push(point);
        }
        let bounds = std::array::from_fn(|axis| {
            [
                corners.iter().map(|p| p[axis].clone()).min().unwrap(),
                corners.iter().map(|p| p[axis].clone()).max().unwrap(),
            ]
        });
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
        assert_eq!(
            result.preparation_work_used(),
            plan.preparation_work_used() + p.work_used()
        );
        assert_eq!(
            result.query_work_used() + result.preparation_work_used(),
            result.work_used()
        );
        assert_eq!(result.work_used() + result.remaining_work(), WORK);
        assert!(result.accounted_peak_bytes() <= BYTES);
        let det = p.determinant_bounds();
        let volume = interval(original) * interval(det);
        assert!(volume.valid() && volume.lo > 0.);
        let point = source.point(0).unwrap();
        let entity = point.entity();
        let record = source.original_component(entity, 0).unwrap().to_owned();
        let proof = result.embedding().connectivity();
        let distance = result.position_error_mm();
        drop(plan);
        drop(source);
        assert_eq!(
            result.source().original_component(entity, 0),
            Some(record.as_str())
        );
        let mesh = result.into_mesh();
        let independent =
            super::super::tests::detached_geometry(&mesh, bounds, [volume.lo, volume.hi], distance);
        assert_eq!(independent, proof);
    }
}
