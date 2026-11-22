use super::super::super::common_mouth::tests::derived_thread_probe::with_head;
use super::*;
fn text(scale: f64, posed: bool) -> String {
    super::super::super::super::branches::representative::scaled_text(scale, posed)
}
#[test]
fn actual_polynomial_flux_matches_independent_signed_integrals_and_large_z_translation() {
    let w = tensor::weights();
    for (kind, expected) in [(0, -243.), (1, -1640.25), (2, -93312.)] {
        for shift in [0., 1e6] {
            let p: [[f64; 3]; 100] = std::array::from_fn(|i| {
                let [u, v] = [i / 10, i % 10];
                match kind {
                    0 => [0.5 + u as f64, 3., shift + 2. * v as f64],
                    1 => [u as f64, v as f64, shift + 1. + (u * v) as f64],
                    _ => [
                        u as f64,
                        v as f64,
                        shift + (u * u.saturating_sub(1) * v * v.saturating_sub(1)) as f64,
                    ],
                }
            });
            let t = crate::tapered_thread::RetainedTensor::new([1, 1], p.to_vec()).unwrap();
            let value = tensor::cell(std::array::from_fn(|i| t.poles()[i]), &w);
            assert!(
                value.lo <= expected && value.hi >= expected,
                "{kind} shift{shift} {value:?}"
            );
        }
    }
}
#[test]
fn actual_decoded_closed_material_has_new_component_and_affine_volume_intervals() {
    for (scale, posed) in [(1., false), (1.137, true)] {
        let input = text(scale, posed);
        with_head(&input, &mut |head| {
            let thread = head
                .derived_thread_carrier(Default::default(), None)
                .unwrap()
                .certify_occupied(1e-7, 1e-5, None)
                .unwrap();
            let material = head
                .closed_repaired_material(&thread, Default::default(), None)
                .unwrap();
            let value = material.volume_bounds(Default::default(), None).unwrap();
            assert!(std::ptr::eq(value.source(), &material));
            let local = interval(value.local_volume_mm3());
            let world = interval(value.world_volume_mm3());
            assert!(local.lo > 0. && world.lo > 0. && world.width() <= 1e-6);
            // Independent exact-dyadic determinant of the ACTUAL four-column
            // chosen source frame, not an assumed rigid/normalized pose.
            use num_rational::BigRational as Q;
            use num_traits::ToPrimitive;
            let q = |x: f64| Q::from_float(x).unwrap();
            let [_, u, v, n] = material.chart().columns.map(|p| p.map(q));
            let d = &u[0] * (&v[1] * &n[2] - &v[2] * &n[1])
                - &u[1] * (&v[0] * &n[2] - &v[2] * &n[0])
                + &u[2] * (&v[0] * &n[1] - &v[1] * &n[0]);
            let near = d.to_f64().unwrap();
            let rounded = q(near);
            let lo = if rounded <= d { near } else { near.next_down() };
            let hi = if rounded >= d { near } else { near.next_up() };
            assert!(q(lo) <= d && q(hi) >= d);
            if posed {
                assert_ne!(d, Q::from_integer(1.into()));
            }
            let transported = local * I::new(lo, hi);
            assert!(world.lo <= transported.lo && world.hi >= transported.hi);

            let components = value.component_volumes_mm3().map(interval);
            assert!(components.into_iter().all(|x| x.valid() && x.lo > 0.));
            let reference = material.thread().volume_with_interrupt(1e-5, None).unwrap();
            assert!(
                components[1].lo >= reference.volume_mm3[0]
                    && components[1].hi <= reference.volume_mm3[1]
            );
            assert!(local.lo >= components[1].lo);
            assert_eq!(value.thread_work(), value.thread_cells() * CELL_WORK);
            assert!(
                value.head_cells() > 0 && value.cavity_evaluations() > 0 && value.head_work() > 0
            );
            assert!(value.work_used() > material.work_used() && value.work_used() < WORK_LIMIT);
            eprintln!(
                "MATERIALVOLUME scale{scale} pose{posed} local{:?} world{:?} components{:?} threadcells{} headcells{} cavevals{} work{}",
                value.local_volume_mm3(),
                value.world_volume_mm3(),
                value.component_volumes_mm3(),
                value.thread_cells(),
                value.head_cells(),
                value.cavity_evaluations(),
                value.work_used()
            );
            assert_eq!(material.body().literal_input(), input);
        });
    }
}
#[test]
fn real_volume_queries_enforce_work_memory_precision_and_publication_cancellation() {
    let input = text(1., false);
    with_head(&input, &mut |head| {
        let thread = head
            .derived_thread_carrier(Default::default(), None)
            .unwrap()
            .certify_occupied(1e-7, 1e-5, None)
            .unwrap();
        let source = head
            .closed_repaired_material(&thread, Default::default(), None)
            .unwrap();
        let receipt = source.volume_bounds(Default::default(), None).unwrap();
        let prefix = source.retained_budget_bytes()
            + std::mem::size_of::<MaterialVolume<'_, '_, '_>>()
            + BASE_SCRATCH;
        let scratch = STREAM_SCRATCH
            .max(RepairedHeadVolume::scratch_budget_bytes(
                receipt.head_cells(),
            ))
            .max(Cutter::volume_scratch_budget_bytes(
                receipt.cavity_evaluations(),
            ));
        let options = MaterialVolumeOptions {
            maximum_head_cells: receipt.head_cells(),
            maximum_cavity_evaluations: receipt.cavity_evaluations(),
            maximum_work: receipt.work_used(),
            maximum_bytes: prefix + scratch,
            ..Default::default()
        };
        assert!(source.volume_bounds(options, None).is_ok());
        for options in [
            MaterialVolumeOptions {
                maximum_work: 0,
                ..options
            },
            MaterialVolumeOptions {
                maximum_work: 1 + source.work_used() + META_WORK - 1,
                ..options
            },
            MaterialVolumeOptions {
                maximum_work: 1 + source.work_used() + META_WORK + CELL_WORK - 1,
                ..options
            },
            MaterialVolumeOptions {
                maximum_work: options.maximum_work - 1,
                ..options
            },
            MaterialVolumeOptions {
                maximum_work: 1
                    + source.work_used()
                    + META_WORK
                    + receipt.thread_work()
                    + receipt.head_work(),
                ..options
            },
            MaterialVolumeOptions {
                maximum_head_cells: 1,
                ..options
            },
            MaterialVolumeOptions {
                maximum_cavity_evaluations: 1,
                ..options
            },
            MaterialVolumeOptions {
                maximum_bytes: 0,
                ..options
            },
            MaterialVolumeOptions {
                maximum_bytes: prefix + STREAM_SCRATCH - 1,
                ..options
            },
            MaterialVolumeOptions {
                maximum_bytes: options.maximum_bytes - 1,
                ..options
            },
        ] {
            assert!(matches!(
                source.volume_bounds(options, None),
                Err(PhillipsReadError::Limit)
            ));
        }
        for options in [
            MaterialVolumeOptions {
                maximum_width_mm3: f64::NAN,
                ..options
            },
            MaterialVolumeOptions {
                maximum_width_mm3: 0.,
                ..options
            },
            MaterialVolumeOptions {
                maximum_head_cells: 0,
                ..options
            },
            MaterialVolumeOptions {
                maximum_head_cells: 65537,
                ..options
            },
            MaterialVolumeOptions {
                maximum_cavity_evaluations: 0,
                ..options
            },
            MaterialVolumeOptions {
                maximum_cavity_evaluations: 1_000_001,
                ..options
            },
            MaterialVolumeOptions {
                maximum_work: WORK_LIMIT + 1,
                ..options
            },
            MaterialVolumeOptions {
                maximum_bytes: BYTES + 1,
                ..options
            },
        ] {
            assert!(matches!(
                source.volume_bounds(options, None),
                Err(PhillipsReadError::Invalid)
            ));
        }
        assert!(matches!(
            source.volume_bounds(
                MaterialVolumeOptions {
                    maximum_width_mm3: 1e-300,
                    ..options
                },
                None
            ),
            Err(PhillipsReadError::Unsupported)
        ));
        let count = std::cell::Cell::new(0);
        source
            .volume_bounds(
                options,
                Some(&|| {
                    count.set(count.get() + 1);
                    false
                }),
            )
            .unwrap();
        let after_thread = 3 + 2 * receipt.thread_cells();
        // Observe the published work-cell counts to reach the real adaptive
        // head publication/debit boundary without production stage hooks.
        let lines = source
            .head()
            .profile_proposal()
            .segments()
            .filter(|s| matches!(s.curve, crate::head_profile::Curve::Line))
            .count();
        let head_polls = source.head().segment_count() + 3 + 2 * (receipt.head_cells() - lines);
        let after_head = after_thread + head_polls;
        for at in [
            0,
            1,
            2,
            3,
            4,
            after_thread - 2,
            after_thread - 1,
            after_thread,
            after_thread + 1,
            after_thread + 2,
            after_head - 1,
            after_head,
            count.get() / 2,
            count.get() - 3,
            count.get() - 2,
            count.get() - 1,
        ] {
            let calls = std::cell::Cell::new(0);
            assert!(
                matches!(
                    source.volume_bounds(
                        options,
                        Some(&|| {
                            let n = calls.get();
                            calls.set(n + 1);
                            n == at
                        })
                    ),
                    Err(PhillipsReadError::Cancelled)
                ),
                "poll{at}/{}",
                count.get()
            );
        }
        let loose = source
            .volume_bounds(
                MaterialVolumeOptions {
                    maximum_width_mm3: 100.,
                    ..Default::default()
                },
                None,
            )
            .unwrap();
        assert!(
            loose.world_volume_mm3()[0] > 0. && interval(loose.world_volume_mm3()).width() <= 100.
        );
        assert!(source.volume_bounds(options, None).is_ok());
        assert_eq!(source.body().literal_input(), input);
    });
}
fn huge_head_text() -> String {
    use crate::{
        head_profile::{Head, Profile},
        helical_thread::{Spec, Thread},
        phillips_drive::{Cavity, PhillipsScrew},
        profiled_screw::SourceSpec,
    };
    let (s, d) = crate::phillips_drive::test_fixture(0);
    let cavity = Cavity::new(Cutter::new(s, d).unwrap(), 1e-5, 100_000).unwrap();
    let head = Head::Cylindrical {
        radius_mm: 1e5,
        height_mm: 1e5,
    };
    let thread = Thread::new(Spec {
        major_radius_mm: 0.97,
        pitch_mm: 0.4,
        length_mm: 2.,
        depth_mm: 0.24,
        crest_width_mm: 0.05,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: false,
        root_radius_mm: None,
    })
    .unwrap();
    let screw = PhillipsScrew::new(
        thread,
        Profile::new_for_thread_union(head, 1., 0.).unwrap(),
        SourceSpec {
            head: Some(head),
            shank_length_mm: 0.,
            minimum_wall_mm: 1e-5,
        },
        cavity,
        100_000,
        &|| false,
    )
    .unwrap();
    crate::step::write_phillips_screw(&screw, [0.; 3], [1., 0., 0.], [0., 1., 0.], 1e-7).unwrap()
}
#[test]
fn actual_large_material_refuses_unrepresentable_absolute_volume_accuracy() {
    let input = huge_head_text();
    with_head(&input, &mut |head| {
        let thread = head
            .derived_thread_carrier(Default::default(), None)
            .unwrap()
            .certify_occupied(1e-7, 1e-5, None)
            .unwrap();
        let material = head
            .closed_repaired_material(&thread, Default::default(), None)
            .unwrap();
        // Actual volume is ~pi*1e15. Its data-derived affine/rounding reserve
        // cannot certify an absolute1e-6 interval and must remain visible.
        assert!(matches!(
            material.volume_bounds(Default::default(), None),
            Err(PhillipsReadError::Unsupported)
        ));
        assert!(
            material
                .volume_bounds(
                    MaterialVolumeOptions {
                        maximum_width_mm3: 1e4,
                        ..Default::default()
                    },
                    None
                )
                .is_ok()
        );
    });
}
