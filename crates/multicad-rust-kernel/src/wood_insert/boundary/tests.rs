use super::*;
use crate::{authored_flat_thread::SourceOptions, helical_thread::Spec};
fn source(outer_left: bool, inner_left: bool, outer_starts: u32, inner_starts: u32) -> Insert {
    let options = SourceOptions {
        position_limit_mm: 1e-7,
        normal_angle_limit_rad: 1e-7,
    };
    let spec = Spec {
        major_radius_mm: 5.,
        pitch_mm: 1.9,
        length_mm: 4.7,
        depth_mm: 1.,
        crest_width_mm: 0.4,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: outer_starts,
        left_hand: outer_left,
        root_radius_mm: None,
    };
    let outer = AuthoredThread::new(spec, options).unwrap();
    let inner = AuthoredThread::new(
        Spec {
            major_radius_mm: 1.5,
            pitch_mm: 0.5,
            depth_mm: 0.3,
            crest_width_mm: 0.1,
            starts: inner_starts,
            left_hand: inner_left,
            ..spec
        },
        options,
    )
    .unwrap();
    Insert::new(outer, inner, 0.5, None).unwrap()
}
#[test]
fn real_paired_sources_have_two_annuli_and_one_closed_oriented_torus() {
    for a in [false, true] {
        for b in [false, true] {
            for (outer, inner) in [(1, 1), (2, 3)] {
                let source = source(a, b, outer, inner);
                let graph = source
                    .boundary_with_interrupt(BoundaryOptions::default(), None)
                    .unwrap();
                assert_eq!(graph.source(), &source);
                assert_eq!(graph.source_boundary(Side::Outer).source(), source.outer());
                assert_eq!(graph.source_boundary(Side::Inner).source(), source.inner());
                assert_eq!(graph.face_count(), source.face_count());
                for cap in [FaceKey::StartCap, FaceKey::EndCap] {
                    assert_eq!(graph.face_loops(cap).unwrap().len(), 2);
                }
                let (v, e) = graph.topology_counts();
                assert_eq!(v as isize - e as isize + graph.face_count() as isize - 2, 0);
                assert_eq!(graph.vertices().count(), v);
                assert_eq!(graph.curves().count(), e);
                assert!(graph.work_used() <= BoundaryOptions::default().work);
                assert!(graph.retained_bytes() <= BoundaryOptions::default().bytes);
                for side in [Side::Outer, Side::Inner] {
                    assert_eq!(
                        graph.cell_surface_same_sense(side),
                        graph.source_boundary(side).cell_surface_same_sense()
                            != (side == Side::Inner)
                    );
                    for (i, _) in graph.cells(side).iter().enumerate() {
                        assert!(graph.cell_face(side, i).is_some());
                        assert!(graph.cell_coedges(side, i).unwrap().count() >= 3);
                    }
                }
            }
        }
    }
}
