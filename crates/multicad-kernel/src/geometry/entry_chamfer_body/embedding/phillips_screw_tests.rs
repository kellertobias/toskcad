//! Actual whole-screw literal gate; no native Phillips factory admission.
use super::*;
use multicad_rust_kernel::{
    display::Tessellation,
    head_profile::{Head, Profile},
    helical_thread::{self, Thread},
    helical_thread_mesh::tessellate_phillips_screw,
    phillips_drive::{Cavity, Cutter, PhillipsScrew},
    profiled_screw::SourceSpec,
};
fn source(index: usize) -> PhillipsScrew {
    let (spec, domain) = super::phillips_tests::fixture(index);
    let cavity = Cavity::new(Cutter::new(spec, domain).unwrap(), 1e-5, 100_000).unwrap();
    let (rs, pitch, r, h, crown, land, neck) = [
        (1., 0.4, 2., 1.6, 3.2, 1.25, 0.1),
        (1.5, 0.5, 2.8, 2.4, 3., 2.1, 0.1),
        (1.5, 0.5, 2.75, 1.65, 0., 0., 0.1),
        (2., 0.7, 4., 3.1, 6.5, 3., 0.2),
        (3., 1., 6., 4.6, 10., 4.2, 0.25),
        (4., 1.25, 8., 6., 13., 5.4, 0.4),
    ][index];
    let head = if index == 2 {
        Head::Countersunk {
            radius_mm: r,
            height_mm: h,
            included_angle_deg: 90.,
            edge_chamfer_mm: 0.05,
            neck_radius_mm: neck,
        }
    } else {
        Head::Crowned {
            radius_mm: r,
            height_mm: h,
            crown_radius_mm: crown,
            top_land_radius_mm: land,
            neck_radius_mm: neck,
        }
    };
    let thread = Thread::new(helical_thread::Spec {
        major_radius_mm: rs - 0.03,
        pitch_mm: pitch,
        length_mm: 8.,
        depth_mm: 0.6 * pitch,
        crest_width_mm: 0.125 * pitch,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: false,
        root_radius_mm: None,
    })
    .unwrap();
    PhillipsScrew::new(
        thread,
        Profile::new_for_thread_union(head, rs, 0.).unwrap(),
        SourceSpec {
            head: Some(head),
            shank_length_mm: 0.,
            minimum_wall_mm: 1e-5,
        },
        cavity,
        100_000,
        &|| false,
    )
    .unwrap()
}
#[test]
fn whole_metric_phillips_screw_candidates_pass_actual_literal_embedding() {
    for index in 0..6 {
        let source = source(index);
        for (deflection, angle) in [(0.1, 0.3), (0.02, 0.5)] {
            let start = std::time::Instant::now();
            let candidate = tessellate_phillips_screw(
                &source,
                Tessellation {
                    deflection,
                    angle,
                    triangle_limit: 1_000_000,
                },
                &|| false,
            )
            .unwrap();
            let native = Mesh {
                vertices: candidate.mesh().vertices.clone(),
                indices: candidate.mesh().indices.clone(),
                face_indices: candidate.mesh().face_indices.clone(),
            };
            let pairs = certify(&native).unwrap_or_else(|e| {
                panic!("actualPhillipsScrew {index}/{deflection}/{angle}: {e:?}")
            });
            eprintln!(
                "actualPhillipsScrew {index}/{deflection}/{angle}: {}facets/{pairs}pairs/{:?}",
                native.indices.len() / 3,
                start.elapsed()
            );
            assert!(pairs > 0);
            assert!(candidate.connectivity().signed_volume_mm3 > 0.);
            if index == 0 && deflection == 0.02 {
                // Real top-annulus triangle vertices are pulled into the
                // opposite thread wall. A prior source proof cannot admit it.
                let face = candidate
                    .faces()
                    .iter()
                    .position(|f| *f == multicad_rust_kernel::phillips_drive::ScrewFace::HeadTop)
                    .unwrap();
                let tri = native
                    .face_indices
                    .iter()
                    .position(|&f| f == face as u32)
                    .unwrap();
                let mut changed = native.clone();
                let id = changed.indices[3 * tri] as usize;
                changed.vertices[3 * id] = 0.;
                changed.vertices[3 * id + 1] = 0.;
                changed.vertices[3 * id + 2] = 6.;
                assert!(certify(&changed).is_err());
            }
        }
    }
}
