use super::*;
use crate::{authored_flat_thread::SourceOptions, helical_thread::Spec};
fn source(left_hand: bool, starts: u32) -> AuthoredThread {
    AuthoredThread::new(
        Spec {
            major_radius_mm: 3.1,
            pitch_mm: 0.9,
            length_mm: 4.7,
            depth_mm: 0.45,
            crest_width_mm: 0.17,
            lower_flank_deg: 24.,
            upper_flank_deg: 32.,
            starts,
            left_hand,
            root_radius_mm: None,
        },
        SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-7,
        },
    )
    .unwrap()
}
#[test]
fn bare_actual_source_has_closed_oriented_owned_geometry_without_a_treatment() {
    for left in [false, true] {
        for starts in [1, 2, 3] {
            let source = source(left, starts);
            let graph = source
                .boundary_with_interrupt(BoundaryOptions::default(), None)
                .unwrap();
            assert_eq!(graph.source(), &source);
            assert!(graph.material().is_none());
            assert_eq!(graph.face_count(), 2 + 4 * starts);
            assert_eq!(graph.cell_surface_same_sense(), left);
            assert!(graph.cells().iter().all(|c| !c.is_entry_cone()));
            let (v, e) = graph.topology_counts();
            assert_eq!(v as isize - e as isize + graph.face_count() as isize, 2);
            let mut uses = std::collections::BTreeMap::<usize, Vec<bool>>::new();
            for face in 0..graph.face_count() {
                let loops = graph.face_loops(face, None).unwrap();
                assert_eq!(loops.len(), 1);
                for chain in loops {
                    for (i, edge) in chain.iter().enumerate() {
                        let next = chain[(i + 1) % chain.len()];
                        let a = &graph.graph.curves[edge.curve];
                        let b = &graph.graph.curves[next.curve];
                        assert_eq!(
                            a.endpoints[usize::from(edge.forward)],
                            b.endpoints[usize::from(!next.forward)]
                        );
                        uses.entry(edge.curve).or_default().push(edge.forward);
                    }
                }
            }
            assert_eq!(uses.len(), e);
            for sides in uses.values() {
                assert_eq!(sides.len(), 2);
                assert_ne!(sides[0], sides[1]);
            }
            for (id, curve) in graph.curves() {
                let polyline = graph.curve_polyline(id, 0.05, 10000, None).unwrap();
                assert!(polyline.maximum_error_mm <= 0.05);
                assert_eq!(
                    polyline.points_mm.first().copied(),
                    Some(graph.curve_point(id, 0., None).unwrap().point_mm)
                );
                assert_eq!(
                    polyline.points_mm.last().copied(),
                    Some(graph.curve_point(id, 1., None).unwrap().point_mm)
                );
                for (end, t) in [0., 1.].into_iter().enumerate() {
                    let sample = graph.curve_point(id, t, None).unwrap();
                    let bounds = graph.points()[curve.endpoints()[end] as usize].position_bounds();
                    for (p, b) in sample.point_mm.into_iter().zip(bounds) {
                        assert!(p >= b.lower && p <= b.upper);
                    }
                }
            }
            assert_eq!(graph.source(), &source);
        }
    }
}
#[test]
fn bare_source_cancellation_and_resource_failures_do_not_substitute_an_owner() {
    let source = source(false, 1);
    assert!(matches!(
        source.boundary_with_interrupt(BoundaryOptions::default(), Some(&|| true)),
        Err(Error::Cancelled)
    ));
    for options in [
        BoundaryOptions {
            cells: 0,
            ..BoundaryOptions::default()
        },
        BoundaryOptions {
            work: 0,
            ..BoundaryOptions::default()
        },
    ] {
        assert!(matches!(
            source.boundary_with_interrupt(options, None),
            Err(Error::Invalid)
        ));
    }
    for options in [
        BoundaryOptions {
            work: 1,
            ..BoundaryOptions::default()
        },
        BoundaryOptions {
            cells: 1,
            ..BoundaryOptions::default()
        },
        BoundaryOptions {
            bytes: EXACT_SCRATCH - 1,
            ..BoundaryOptions::default()
        },
    ] {
        assert!(matches!(
            source.boundary_with_interrupt(options, None),
            Err(Error::Limit)
        ));
    }
    let graph = source
        .boundary_with_interrupt(BoundaryOptions::default(), None)
        .unwrap();
    assert_eq!(graph.source(), &source);
}
