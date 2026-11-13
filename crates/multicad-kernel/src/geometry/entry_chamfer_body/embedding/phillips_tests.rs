//! Independent whole-literal embedding checks; no native Phillips admission.
use super::*;
use multicad_rust_kernel::{
    display::Tessellation,
    phillips_drive::{Cutter, Domain, Inner, Spec},
};
fn fixture(index: usize) -> (Spec, Domain) {
    let (m, r, b, g, inner, beta, head_r, head_h) = match index {
        0 => (
            1.9,
            0.3,
            0.61,
            0.81,
            Inner::Rounded { radius_mm: 0.25 },
            7.,
            2.,
            1.6,
        ),
        1 => (
            3.,
            0.5,
            0.955,
            1.27,
            Inner::ShallowV {
                angle_deg: 138.,
                e_width_mm: 0.435,
            },
            7.,
            2.8,
            2.4,
        ),
        2 => (
            3.2,
            0.5,
            0.955,
            1.27,
            Inner::ShallowV {
                angle_deg: 138.,
                e_width_mm: 0.435,
            },
            7.,
            2.75,
            1.65,
        ),
        3 => (
            4.4,
            0.6,
            1.455,
            2.29,
            Inner::ShallowV {
                angle_deg: 140.,
                e_width_mm: 0.815,
            },
            5.75,
            4.,
            3.1,
        ),
        4 => (
            6.9,
            0.8,
            2.395,
            3.81,
            Inner::ShallowV {
                angle_deg: 146.,
                e_width_mm: 2.005,
            },
            5.75,
            6.,
            4.6,
        ),
        _ => (
            9.,
            1.,
            3.465,
            5.08,
            Inner::ShallowV {
                angle_deg: 153.,
                e_width_mm: 2.415,
            },
            7.,
            8.,
            6.,
        ),
    };
    let spec = Spec {
        mouth_diameter_mm: m,
        mouth_radius_mm: r,
        b_width_mm: b,
        g_diameter_mm: g,
        tilt_angle_deg: beta,
        outer_flank_angle_deg: 26.5,
        floor_angle_deg: 28.,
        outer_corner_angle_deg: 92.,
        inner,
    };
    let tip = (m - g) / (2. * 26.5_f64.to_radians().tan()) + g / 2. * 28_f64.to_radians().tan();
    (
        spec,
        Domain {
            corner_start_depth_mm: -2. * m,
            corner_end_depth_mm: tip + m,
            corner_reach_mm: 4. * (head_r + m + (tip - head_h).abs() + 1.),
        },
    )
}
#[test]
fn actual_metric_phillips_candidates_pass_literal_embedding_at_both_qualities() {
    for metric in 0..=5 {
        let (spec, domain) = fixture(metric);
        let source = Cutter::new(spec, domain).unwrap();
        let atlas = source.surface_patches().unwrap();
        let volume = source.volume_bounds_mm3(0.001, 100_000).unwrap();
        for (deflection, angle) in [(0.1, 0.3), (0.02, 0.5)] {
            let start = std::time::Instant::now();
            let candidate = atlas
                .boundary_mesh(Tessellation {
                    deflection,
                    angle,
                    triangle_limit: 1_000_000,
                })
                .unwrap();
            let native = Mesh {
                vertices: candidate.mesh().vertices.clone(),
                indices: candidate.mesh().indices.clone(),
                face_indices: candidate.mesh().face_indices.clone(),
            };
            let pairs = certify(&native).unwrap_or_else(|error| {
                panic!("metric {metric}, delta {deflection}, angle {angle}: {error:?}")
            });
            assert!(candidate.connectivity().signed_volume_mm3 > 0.);
            assert!(
                (candidate.connectivity().signed_volume_mm3 - volume.volume_mm3[0]).abs()
                    < volume.volume_mm3[1] * 0.03
            );
            eprintln!(
                "literal Phillips metric {metric}, delta {deflection}, angle {angle}: {} facets, {pairs} checked pairs, {:?}, source volume {:?}, mesh volume {}",
                native.indices.len() / 3,
                start.elapsed(),
                volume.volume_mm3,
                candidate.connectivity().signed_volume_mm3
            );
            if metric == 0 && deflection == 0.02 {
                let mut changed = native.clone();
                let center = changed
                    .vertices
                    .as_chunks::<3>()
                    .0
                    .iter()
                    .position(|p| *p == [0.; 3])
                    .unwrap();
                let below = changed
                    .vertices
                    .as_chunks::<3>()
                    .0
                    .iter()
                    .map(|p| p[2])
                    .fold(0., f64::max)
                    + 1.;
                changed.vertices[3 * center + 2] = below;
                assert!(
                    certify(&changed).is_err(),
                    "real intersecting cap must not inherit candidate source proof"
                );
            }
        }
    }
}
