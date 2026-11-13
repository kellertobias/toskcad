use super::*;
use crate::authored_flat_thread::{AuthoredThread, SourceOptions};
use crate::helical_thread::Spec;
pub(super) fn source(left_hand: bool, starts: u32, entry: Option<f64>, tail: bool) -> Cavity {
    let pitch_mid = (2.675 + 2.775) / 2.;
    let minor_mid = (2.459 + 2.599) / 2.;
    let thread = AuthoredThread::new(
        Spec {
            major_radius_mm: 1.5,
            pitch_mm: 0.5,
            length_mm: 3.3,
            depth_mm: (3. - minor_mid) / 2.,
            crest_width_mm: 0.25 - (3. - pitch_mid) / 3_f64.sqrt(),
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            starts,
            left_hand,
            root_radius_mm: None,
        },
        SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-7,
        },
    )
    .unwrap();
    Cavity::new(
        thread,
        Intent {
            entry_depth_mm: entry,
            terminal: tail.then_some(Terminal {
                width_mm: 0.3,
                radial_depth_mm: 0.05,
            }),
        },
    )
    .unwrap()
}
pub(super) fn options() -> Tessellation {
    Tessellation {
        deflection: 0.1,
        angle: 0.3,
        triangle_limit: 1_000_000,
    }
}
#[test]
fn actual_combined_source_candidate_has_whole_facet_receipts_and_closed_indexed_connectivity() {
    for (left, starts) in [(false, 1), (true, 2)] {
        let source = source(left, starts, Some(0.4), true);
        let timer = std::time::Instant::now();
        let candidate = tessellate_candidate(&source, options(), None).unwrap();
        eprintln!(
            "authored cavity candidate {left}/{starts}: {} triangles in {:?}, {} work",
            candidate.mesh.face_indices.len(),
            timer.elapsed(),
            candidate.work_used
        );
        assert!(candidate.source_chord_bound_mm <= options().deflection);
        assert!(candidate.source_normal_angle_bound_rad <= options().angle);
        let proof = crate::manufacturing::certify(&candidate.mesh).unwrap();
        assert!(proof.edge_pairs > 0 && proof.vertex_links > 0 && proof.signed_volume_mm3 > 0.);
        // Connectivity is intentionally not a global embedding proof. Native
        // publication still requires the literal triangle intersection gate.
    }
}
#[test]
fn unsupported_threshold_cancel_and_resource_failure_leave_source_unchanged() {
    let body = source(false, 1, Some(0.4), true);
    let before = body;
    assert_eq!(
        tessellate_candidate(
            &body,
            Tessellation {
                deflection: 0.,
                ..options()
            },
            None
        )
        .unwrap_err(),
        Error::Invalid
    );
    assert_eq!(
        tessellate_candidate(
            &body,
            Tessellation {
                triangle_limit: 1,
                ..options()
            },
            None
        )
        .unwrap_err(),
        Error::Limit
    );
    assert_eq!(
        tessellate_candidate(&body, options(), Some(&|| true)).unwrap_err(),
        Error::Cancelled
    );
    let small = source(false, 1, Some(0.01), true);
    assert_eq!(
        tessellate_candidate(&small, options(), None).unwrap_err(),
        Error::Unresolved
    );
    assert_eq!(body, before);
}
#[test]
fn actual_terminal_and_entry_only_sources_close_for_both_hands_and_multiple_starts() {
    for (left, starts, entry, tail) in [
        (false, 1, None, true),
        (true, 2, None, true),
        (false, 1, Some(0.4), false),
        (true, 8, Some(0.4), true),
    ] {
        let body = source(left, starts, entry, tail);
        let result = tessellate_candidate(&body, options(), None).unwrap();
        let topology = crate::manufacturing::certify(&result.mesh).unwrap();
        assert!(topology.signed_volume_mm3 > 0.);
        assert!(result.source_chord_bound_mm <= options().deflection);
        assert!(result.source_normal_angle_bound_rad <= options().angle);
        assert!(result.work_used <= 16_000_000);
    }
}
#[test]
fn interrupted_exact_clipping_never_publishes_partial_source_and_retry_is_identical() {
    let body = source(true, 2, Some(0.4), true);
    let original = body;
    let before = tessellate_candidate(&body, options(), None).unwrap();
    for stop in [10, 10_000, 100_000] {
        let calls = std::cell::Cell::new(0);
        let check = || {
            calls.set(calls.get() + 1);
            calls.get() > stop
        };
        assert_eq!(
            tessellate_candidate(&body, options(), Some(&check)).unwrap_err(),
            Error::Cancelled
        );
        assert_eq!(body, original);
    }
    let after = tessellate_candidate(&body, options(), None).unwrap();
    assert_eq!(before.mesh, after.mesh);
    assert_eq!(before.work_used, after.work_used);
}
#[test]
fn exact_source_mesh_refuses_unrepresentable_quality_and_limits_during_shell_and_end_closure() {
    let body = source(false, 1, Some(0.4), true);
    let unchanged = body;
    for limit in [1000, 4160, 4230, 4232, 4432] {
        assert_eq!(
            tessellate_candidate(
                &body,
                Tessellation {
                    triangle_limit: limit,
                    ..options()
                },
                None
            )
            .unwrap_err(),
            Error::Limit
        );
    }
    assert_eq!(
        tessellate_candidate(
            &body,
            Tessellation {
                deflection: 1e-300,
                ..options()
            },
            None
        )
        .unwrap_err(),
        Error::Unresolved
    );
    assert_eq!(
        tessellate_candidate(
            &body,
            Tessellation {
                angle: 0.01,
                ..options()
            },
            None
        )
        .unwrap_err(),
        Error::Limit
    );
    assert_eq!(body, unchanged);
}
#[test]
fn exact_seam_and_caps_on_grid_corners_keep_closed_source_boundary_coverage() {
    use regular::Patch;
    use std::cmp::Ordering;
    for (left_hand, starts) in [(false, 1), (true, 1), (false, 2), (true, 2)] {
        let thread = AuthoredThread::new(
            Spec {
                major_radius_mm: 1.5,
                pitch_mm: 0.5,
                length_mm: 4.,
                depth_mm: 0.25,
                crest_width_mm: 0.125,
                lower_flank_deg: 30.,
                upper_flank_deg: 30.,
                starts,
                left_hand,
                root_radius_mm: None,
            },
            SourceOptions {
                position_limit_mm: 1e-7,
                normal_angle_limit_rad: 1e-7,
            },
        )
        .unwrap();
        let body = Cavity::new(
            thread,
            Intent {
                entry_depth_mm: Some(0.5),
                terminal: Some(Terminal {
                    width_mm: 0.5,
                    radial_depth_mm: 0.125,
                }),
            },
        )
        .unwrap();
        let w = &mut Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        let patch = Patch::new(&body, 0, 0, 1, w).unwrap();
        let seam_row = 32 / starts;
        let cone = patch.cell(0, seam_row - 1, 64, Kind::Cone, w).unwrap();
        let wall = patch.cell(0, seam_row, 64, Kind::Wall, w).unwrap();
        let a = cone
            .points
            .iter()
            .find(|p| p.active & 0b101001 == 0b101001)
            .unwrap();
        let b = wall
            .points
            .iter()
            .find(|p| p.active & 0b010101 == 0b010101)
            .unwrap();
        let (ka, ca, _) = patch.owned(a, &cone, starts, w).unwrap();
        let (kb, cb, _) = patch.owned(b, &wall, starts, w).unwrap();
        assert_eq!(ka, kb);
        assert_eq!(ca, cb);
        assert_eq!(
            a.q.cmp(&exact::Q::integer(0, w).unwrap(), w).unwrap(),
            Ordering::Equal
        );
        // At topZ3.5, period7/q0/v0 has only a zero-area cell.
        // The adjacent positive cell supplies that SAME source boundary point.
        let zero = Patch::new(&body, 7, 0, 1, w)
            .unwrap()
            .cell(0, 0, 64, Kind::Wall, w)
            .unwrap();
        assert!(zero.points.is_empty());
        let adjacent = Patch::new(&body, 6, 4, 1, w)
            .unwrap()
            .cell(0, 0, 64, Kind::Wall, w)
            .unwrap();
        assert!(
            adjacent
                .points
                .iter()
                .any(|p| p.active & 0b100110 == 0b100110)
        );
        let result = tessellate_candidate(&body, options(), None).unwrap();
        let proof = crate::manufacturing::certify(&result.mesh).unwrap();
        assert!(proof.signed_volume_mm3 > 0.);
        assert!(result.mesh.face_indices.contains(&0) && result.mesh.face_indices.contains(&1));
        let theta = I::point(0.)
            * I::new(
                std::f64::consts::TAU.next_down(),
                std::f64::consts::TAU.next_up(),
            )
            * I::point(if left_hand { -1. } else { 1. });
        for (r, z) in [(1.75, 0.), (1.625, 4.)] {
            let actual = [
                I::point(r) * theta.cos(),
                I::point(r) * theta.sin(),
                I::point(z),
            ];
            assert!(
                result
                    .mesh
                    .vertices
                    .as_chunks::<3>()
                    .0
                    .iter()
                    .any(|p| p.iter().zip(actual).all(|(p, a)| a.lo <= *p && *p <= a.hi))
            );
        }
    }
}
