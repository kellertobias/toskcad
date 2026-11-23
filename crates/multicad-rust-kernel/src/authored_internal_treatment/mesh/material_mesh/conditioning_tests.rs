use super::*;
use crate::{
    Circle, Rational,
    authored_flat_thread::{AuthoredThread, SourceOptions},
    cylinder::AxialCylinder,
    helical_thread::Spec,
};
fn source(left: bool, width: f64) -> Material {
    let thread = AuthoredThread::new(
        Spec {
            major_radius_mm: 1.5,
            pitch_mm: 0.5,
            length_mm: 3.3000000000000003,
            depth_mm: 0.27063293868263705,
            crest_width_mm: 0.0625,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            starts: 1,
            left_hand: left,
            root_radius_mm: None,
        },
        SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-7,
        },
    )
    .unwrap();
    let end = Rational::new(33_000_000_000_000_003, 10_000_000_000_000_000).unwrap();
    let body = Cavity::new_with_exact_receiving_extent(
        thread,
        Intent {
            entry_depth_mm: None,
            terminal: Some(Terminal {
                width_mm: width,
                radial_depth_mm: 0.075,
            }),
        },
        end,
    )
    .unwrap();
    let zero = Rational::new(0, 1).unwrap();
    Material::new(
        AxialCylinder::new(
            7,
            Circle::new([zero; 2], Rational::new(7, 2).unwrap()).unwrap(),
            zero,
            end,
        )
        .unwrap(),
        body,
        0.5,
    )
    .unwrap()
}
#[test]
fn actual_unshifted_tail_chart_refusal_keeps_failed_work_and_source() {
    let options = Tessellation {
        deflection: 0.1,
        angle: 0.3,
        triangle_limit: 1_000_000,
    };
    let material = source(false, 0.35);
    let mut work = Work {
        remaining: 16_000_000,
        used: 0,
        interrupt: None,
    };
    // General refusal of the first (unshifted, 64-row) chart at whatever stage
    // refuses it. Under portable maths band 1's knot is one ulp off row 24, so
    // the chart repartitions before closure; the assembler's stock-trace
    // refusal is covered directly by
    // constructed_coarse_stock_trace_refuses_the_outer_wall_without_publishing_stock.
    let first =
        match chart::traced_grid_attempt(material.cavity(), options, &Grid::new(), &mut work)
            .unwrap()
        {
            Attempt::Complete((shell, trace)) => {
                match close::prepare(shell, material.cavity(), options, &mut work).unwrap() {
                    Some(prepared) => {
                        assemble(&material, prepared, trace, options, &mut work).unwrap()
                    }
                    None => None,
                }
            }
            Attempt::Refine | Attempt::Repartition => None,
        };
    assert!(first.is_none());
    let failed_work = work.used;
    assert_eq!(failed_work, 582_905);
    let before = material;
    let calls = std::cell::Cell::new(0usize);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    let candidate = material
        .tessellate_candidate(options, Some(&count))
        .unwrap();
    assert!(candidate.work_used > failed_work && candidate.work_used < 16_000_000);
    let total_calls = calls.get();
    for stop in [1, total_calls / 2, total_calls - 1] {
        calls.set(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() >= stop
        };
        assert!(matches!(
            material.tessellate_candidate(options, Some(&cancel)),
            Err(Error::Cancelled)
        ));
        assert_eq!(material, before);
    }
    assert_eq!(
        material.tessellate_candidate(options, None).unwrap().mesh,
        candidate.mesh
    );
    let mut tiny = options;
    tiny.triangle_limit = 16;
    assert!(matches!(
        material.tessellate_candidate(tiny, None),
        Err(Error::Limit)
    ));
}
#[test]
fn actual_nearby_terminal_stations_keep_original_domains_and_all_facet_receipts() {
    let center = 0.35_f64;
    let widths = [
        0.3,
        0.4,
        f64::from_bits(center.to_bits() - 2),
        center,
        f64::from_bits(center.to_bits() + 2),
    ];
    for left in [false, true] {
        for width in widths {
            let material = source(left, width);
            let before = material;
            for deflection in [0.1, 0.02] {
                let options = Tessellation {
                    deflection,
                    angle: 0.3,
                    triangle_limit: 500_000,
                };
                let mesh = material.tessellate_candidate(options, None).unwrap();
                assert!(mesh.work_used < 16_000_000);
                assert!(mesh.source_chord_bound_mm <= deflection);
                assert!(mesh.source_normal_angle_bound_rad <= 0.3);
                assert!(
                    crate::manufacturing::certify(&mesh.mesh)
                        .unwrap()
                        .signed_volume_mm3
                        > 0.
                );
                assert_eq!(material, before);
            }
        }
    }
}

#[test]
fn actual_canonical_nominal_entry_source_qualifies_the_posed_native_local_quality() {
    let options = Tessellation {
        deflection: 0.01,
        angle: 0.3 * 0.75,
        triangle_limit: 500_000,
    };
    for left in [false, true] {
        for end in [
            Rational::new(33, 10).unwrap(),
            Rational::new(4, 1).unwrap(),
            Rational::new(33_000_000_000_000_001, 10_000_000_000_000_000).unwrap(),
            Rational::new(33_000_000_000_000_003, 10_000_000_000_000_000).unwrap(),
            Rational::new(33_000_000_000_000_005, 10_000_000_000_000_000).unwrap(),
        ] {
            let thread = AuthoredThread::new_with_exact_axial_domain(
                Spec {
                    major_radius_mm: 1.5,
                    pitch_mm: 0.5,
                    length_mm: 3.3000000000000003,
                    depth_mm: 0.27063293868263705,
                    crest_width_mm: 0.0625,
                    lower_flank_deg: 30.,
                    upper_flank_deg: 30.,
                    starts: 1,
                    left_hand: left,
                    root_radius_mm: None,
                },
                end,
                SourceOptions {
                    position_limit_mm: 1e-7,
                    normal_angle_limit_rad: 1e-7,
                },
                None,
            )
            .unwrap();
            let end = thread.exact_axial_extent().unwrap();
            let cavity = Cavity::new_with_exact_receiving_extent(
                thread,
                Intent {
                    entry_depth_mm: Some(0.4),
                    terminal: None,
                },
                end,
            )
            .unwrap();
            let zero = Rational::new(0, 1).unwrap();
            let material = Material::new(
                AxialCylinder::new(
                    7,
                    Circle::new([zero; 2], Rational::new(7, 2).unwrap()).unwrap(),
                    zero,
                    end,
                )
                .unwrap(),
                cavity,
                0.5,
            )
            .unwrap();
            let calls = std::cell::Cell::new(0usize);
            let count = || {
                calls.set(calls.get() + 1);
                false
            };
            let before = material;
            let mesh = material
                .tessellate_candidate(options, Some(&count))
                .unwrap();
            eprintln!(
                "canonicalEntryStrict left={left} end={:?} work={} facets={}",
                end.ratio(),
                mesh.work_used,
                mesh.mesh.face_indices.len()
            );
            assert!(mesh.work_used < 16_000_000);
            if !left
                && end == Rational::new(33_000_000_000_000_003, 10_000_000_000_000_000).unwrap()
            {
                let total = calls.get();
                for stop in [1, total / 2, total - 1] {
                    calls.set(0);
                    let cancel = || {
                        calls.set(calls.get() + 1);
                        calls.get() >= stop
                    };
                    assert!(matches!(
                        material.tessellate_candidate(options, Some(&cancel)),
                        Err(Error::Cancelled)
                    ));
                    assert_eq!(material, before);
                }
                assert_eq!(
                    material.tessellate_candidate(options, None).unwrap().mesh,
                    mesh.mesh
                );
                assert!(matches!(
                    material.tessellate_candidate(
                        Tessellation {
                            triangle_limit: 16,
                            ..options
                        },
                        None
                    ),
                    Err(Error::Limit)
                ));
            }
            assert!(
                mesh.source_chord_bound_mm <= options.deflection
                    && mesh.source_normal_angle_bound_rad <= options.angle
            );
            assert!(
                crate::manufacturing::certify(&mesh.mesh)
                    .unwrap()
                    .signed_volume_mm3
                    > 0.
            );
        }
    }
}
