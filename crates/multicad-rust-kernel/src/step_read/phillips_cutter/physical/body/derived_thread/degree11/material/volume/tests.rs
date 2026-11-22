use super::super::super::tests::actual_m3;
use super::*;
use crate::step_read::phillips_cutter::physical::body::common_mouth::tests::derived_thread_probe::with_head;

fn choose(n: usize, k: usize) -> u64 {
    (0..k.min(n - k)).fold(1u64, |v, i| v * (n - i) as u64 / (i + 1) as u64)
}
fn bernstein_power(i: usize, k: usize) -> I {
    if i < k {
        I::point(0.)
    } else {
        I::point(choose(i, k) as f64) / I::point(choose(11, k) as f64)
    }
}
// x=1, y=-1+2u+a*u²(1-u)v(1-v),
// z=2v+b*u(1-u)v²(1-v). All four boundary curves stay unchanged.
// Nonzero z_u makes the full cross term necessary. The mixed correction
// integrals cancel exactly: integral F_u G_v = integral F_v G_u =1/3600.
fn polynomial() -> [[I; 3]; 144] {
    std::array::from_fn(|n| {
        let i = n / 12;
        let j = n % 12;
        let u = |k| bernstein_power(i, k);
        let v = |k| bernstein_power(j, k);
        [
            I::point(1.),
            -I::point(1.) + I::point(2.) * u(1) + I::point(0.25) * (u(2) - u(3)) * (v(1) - v(2)),
            I::point(2.) * v(1) + I::point(0.25) * (u(1) - u(2)) * (v(2) - v(3)),
        ]
    })
}
fn subdivide(p: [[I; 3]; 144], axis: usize) -> [[[I; 3]; 144]; 2] {
    let mut out = [[[I::point(0.); 3]; 144]; 2];
    for r in 0..12 {
        let index = |s| if axis == 0 { 12 * s + r } else { 12 * r + s };
        let mut row = std::array::from_fn::<_, 12, _>(|s| p[index(s)]);
        out[0][index(0)] = row[0];
        out[1][index(11)] = row[11];
        for depth in 1..12 {
            for j in 0..12 - depth {
                row[j] = std::array::from_fn(|k| (row[j][k] + row[j + 1][k]) * I::point(0.5));
            }
            out[0][index(depth)] = row[0];
            out[1][index(11 - depth)] = row[11 - depth];
        }
    }
    out
}
#[test]
fn full_d11_flux_matches_independent_polynomial_and_both_axis_subdivision_invariance() {
    let q = tensor::Quadrature::new();
    let p = polynomial();
    let whole = tensor::cell(p, &q);
    let cell = crate::tapered_thread::calculation_cell11(p).unwrap();
    let normal = NormalPolynomial11::new(&cell).unwrap();
    let cached = tensor::cached_cell(p, &q, &normal);
    assert_eq!(
        (cached.lo.to_bits(), cached.hi.to_bits()),
        (whole.lo.to_bits(), whole.hi.to_bits())
    );
    assert!(normal.retained_budget_bytes() <= NormalPolynomial11::retained_capacity_bound_bytes());
    assert!(whole.lo <= 2. && whole.hi >= 2.);
    assert!(whole.hi - whole.lo < 1e-10);
    let mut sum = I::point(0.);
    for u in subdivide(p, 0) {
        for uv in subdivide(u, 1) {
            sum = sum + tensor::cell(uv, &q);
        }
    }
    assert!(sum.lo <= 2. && sum.hi >= 2.);
    assert!(sum.hi - sum.lo < 1e-10);
    // Four rotations of this wall close the square prism. Its independently
    // known area4 times height2 gives volume8; horizontal caps contribute0.
    let mut prism = I::point(0.);
    for turn in 0..4 {
        let wall = p.map(|[x, y, z]| match turn {
            0 => [x, y, z],
            1 => [-y, x, z],
            2 => [-x, -y, z],
            _ => [y, -x, z],
        });
        prism = prism + tensor::cell(wall, &q);
    }
    assert!(prism.lo <= 8. && prism.hi >= 8. && prism.hi - prism.lo < 1e-9);
    let reversed = p;
    let reversed = std::array::from_fn(|n| reversed[12 * (11 - n / 12) + n % 12]);
    let negative = tensor::cell(reversed, &q);
    assert!(negative.lo <= -2. && negative.hi >= -2.);
}
#[test]
fn actual_original_flat_and_fitted_m3_volume_both_hands_and_pose() {
    for fitted in [false, true] {
        for left in [false, true] {
            let text = actual_m3(fitted, left);
            with_head(&text, &mut |head| {
                let occupied = head
                    .degree11_thread_carrier(Default::default(), None)
                    .unwrap()
                    .certify_occupied(Default::default(), None)
                    .unwrap();
                let material = head
                    .closed_degree11_material(&occupied, Default::default(), None)
                    .unwrap();
                let volume = material.volume_bounds(Default::default(), None).unwrap();
                eprintln!(
                    "D11_VOLUME fitted={fitted} left={left} local={:?} world={:?} components={:?} total={} fresh={} threadwork={} headwork={} cells={} headcells={} cavityevals={} peak={}",
                    volume.local_volume_mm3(),
                    volume.world_volume_mm3(),
                    volume.component_volumes_mm3(),
                    volume.work_used(),
                    volume.new_work_used(),
                    volume.thread_work(),
                    volume.head_work(),
                    volume.thread_cells(),
                    volume.head_cells(),
                    volume.cavity_evaluations(),
                    volume.accounted_peak_bytes()
                );
                assert!(std::ptr::eq(volume.source(), &material));
                assert!(volume.world_volume_mm3()[0] > 0.);
                assert!(volume.world_volume_mm3()[1] - volume.world_volume_mm3()[0] <= 1e-6);
                assert!(
                    volume.work_used() <= WORK_LIMIT && volume.accounted_peak_bytes() <= MAX_BYTES
                );
                assert_eq!(
                    volume.work_used(),
                    material.work_used() + volume.new_work_used()
                );
                assert_eq!(volume.thread_cells(), if fitted { 315 } else { 270 });
                assert_eq!(volume.normals.len(), volume.thread_cells());
                let original = material.thread().carrier().faces()[0]
                    .source()
                    .cell([0, 0])
                    .unwrap();
                let independent = NormalPolynomial11::new(&original).unwrap();
                let bits = |field: &NormalPolynomial11| {
                    field
                        .coefficients()
                        .iter()
                        .flatten()
                        .map(|p| p.map(|x| (x.lo.to_bits(), x.hi.to_bits())))
                        .collect::<Vec<_>>()
                };
                assert_eq!(bits(&volume.normals[0]), bits(&independent));
                assert!(volume.retained_budget_bytes() < volume.accounted_peak_bytes());
            });
        }
    }
}

#[test]
fn directed_volume_numeric_contract_and_workspace_fitting() {
    assert_eq!(
        certified(I::new(1., 1.000001), 0.01).unwrap(),
        [1., 1.000001]
    );
    for value in [
        I::new(-1., 1.),
        I::new(0., 0.),
        I::new(1., 2.),
        I::new(1., f64::INFINITY),
    ] {
        assert!(certified(value, 0.01).is_err());
    }
    assert_eq!(fitting_limit(16, 45, |n| n * 10), 4);
    assert_eq!(fitting_limit(16, 160, |n| n * 10), 16);
    let width = allocate_component_width(0.01, I::new(1., 1.), I::new(99., 100.)).unwrap();
    assert!(width > 0. && width <= 0.00125);
    assert!(allocate_component_width(0.01, I::new(1., 2.), I::new(99., 100.)).is_err());
    assert!(allocate_component_width(1., I::new(1., f64::MAX), I::new(1., f64::MAX)).is_err());
}
#[test]
fn actual_volume_invalid_resource_stage_cancel_and_exact_retry() {
    use std::cell::Cell;
    let text = actual_m3(true, false);
    with_head(&text, &mut |head| {
        let occupied = head
            .degree11_thread_carrier(Default::default(), None)
            .unwrap()
            .certify_occupied(Default::default(), None)
            .unwrap();
        let material = head
            .closed_degree11_material(&occupied, Default::default(), None)
            .unwrap();
        let default = Degree11MaterialVolumeOptions::default();
        for bad in [
            Degree11MaterialVolumeOptions {
                maximum_width_mm3: f64::NAN,
                ..default
            },
            Degree11MaterialVolumeOptions {
                maximum_width_mm3: 0.,
                ..default
            },
            Degree11MaterialVolumeOptions {
                maximum_head_cells: 0,
                ..default
            },
            Degree11MaterialVolumeOptions {
                maximum_head_cells: 65537,
                ..default
            },
            Degree11MaterialVolumeOptions {
                maximum_cavity_evaluations: 0,
                ..default
            },
            Degree11MaterialVolumeOptions {
                maximum_cavity_evaluations: 1_000_001,
                ..default
            },
            Degree11MaterialVolumeOptions {
                maximum_work: WORK_LIMIT + 1,
                ..default
            },
            Degree11MaterialVolumeOptions {
                maximum_bytes: MAX_BYTES + 1,
                ..default
            },
        ] {
            assert!(matches!(
                material.volume_bounds(bad, None),
                Err(PhillipsReadError::Invalid)
            ));
        }
        assert!(matches!(
            material.volume_bounds(
                Degree11MaterialVolumeOptions {
                    maximum_work: material.work_used(),
                    ..default
                },
                None
            ),
            Err(PhillipsReadError::Limit)
        ));
        let minimum_bytes = material.retained_budget_bytes()
            + std::mem::size_of::<Degree11MaterialVolume<'_, '_, '_, '_>>()
            + BASE_SCRATCH
            + STREAM_SCRATCH;
        assert!(matches!(
            material.volume_bounds(
                Degree11MaterialVolumeOptions {
                    maximum_bytes: minimum_bytes - 1,
                    ..default
                },
                None
            ),
            Err(PhillipsReadError::Limit)
        ));
        let count = Cell::new(0usize);
        let volume = material
            .volume_bounds(
                default,
                Some(&|| {
                    count.set(count.get() + 1);
                    false
                }),
            )
            .unwrap();
        let expected = volume.world_volume_mm3();
        let exact_work = volume.work_used();
        let after_head =
            material.work_used() + 1 + META_WORK + volume.thread_work() + volume.head_work();
        assert!(matches!(
            material.volume_bounds(
                Degree11MaterialVolumeOptions {
                    maximum_work: after_head,
                    ..default
                },
                None
            ),
            Err(PhillipsReadError::Limit)
        ));

        // Independently run the actual same-width public head integral. This
        // checks fresh-prefix accounting and identifies its real return/debit
        // boundary without injecting receipts or altering production geometry.
        let determinant = interval(material.determinant_bounds());
        let mut radius = 0f64;
        for i in 0..material.head().segment_count() {
            radius = radius.max(
                material
                    .head()
                    .certified_interval_bounds(i, [0., 1.])
                    .position_mm[0][1],
            );
        }
        let top = material.head().station(0).unwrap()[1];
        let end = material
            .head()
            .station(material.head().segment_count())
            .unwrap()[1];
        let pi = I::new(
            std::f64::consts::PI.next_down(),
            std::f64::consts::PI.next_up(),
        );
        let outer = pi * I::point(radius).square() * (I::point(end) - I::point(top))
            + I::point(volume.component_volumes_mm3()[1][1]);
        let component_width =
            allocate_component_width(default.maximum_width_mm3, determinant, outer).unwrap();
        let head_polls = Cell::new(0usize);
        let separate_head = material
            .head()
            .volume_bounds(
                RepairedHeadVolumeOptions {
                    maximum_width_mm3: component_width,
                    maximum_cells: default.maximum_head_cells,
                    maximum_work: 16_000_000,
                    maximum_bytes: MAX_BYTES,
                },
                Some(&|| {
                    head_polls.set(head_polls.get() + 1);
                    false
                }),
            )
            .unwrap();
        assert_eq!(
            separate_head.volume_mm3(),
            volume.component_volumes_mm3()[0]
        );
        assert_eq!(
            separate_head.work_used() - material.head().work_used(),
            volume.head_work()
        );
        // Four setup Work polls and two polls per original cell precede the
        // separately counted head query; its next actual poll is the debit.
        let head_debit = 4 + 2 * volume.thread_cells() + head_polls.get() + 1;
        // Actual production callbacks: immediate, original-cell flux, then
        // after full flux in head/cavity integration, and final publication.
        for stop in [1, 10, 640, head_debit, count.get() - 1, count.get()] {
            let calls = Cell::new(0usize);
            assert!(
                matches!(
                    material.volume_bounds(
                        default,
                        Some(&|| {
                            calls.set(calls.get() + 1);
                            calls.get() >= stop
                        })
                    ),
                    Err(PhillipsReadError::Cancelled)
                ),
                "checkpoint {stop}"
            );
            assert_eq!(calls.get(), stop);
        }
        assert!(matches!(
            material.volume_bounds(
                Degree11MaterialVolumeOptions {
                    maximum_work: exact_work - 1,
                    ..default
                },
                None
            ),
            Err(PhillipsReadError::Limit)
        ));
        let retry = material
            .volume_bounds(
                Degree11MaterialVolumeOptions {
                    maximum_work: exact_work,
                    ..default
                },
                None,
            )
            .unwrap();
        assert_eq!(retry.world_volume_mm3(), expected);
        assert_eq!(retry.work_used(), exact_work);
        assert_eq!(retry.head_work(), volume.head_work());
        // Explicitly insufficient actual component integration budgets refuse.
        assert!(matches!(
            material.volume_bounds(
                Degree11MaterialVolumeOptions {
                    maximum_head_cells: 1,
                    ..default
                },
                None
            ),
            Err(PhillipsReadError::Limit)
        ));
        assert!(matches!(
            material.volume_bounds(
                Degree11MaterialVolumeOptions {
                    maximum_cavity_evaluations: 1,
                    ..default
                },
                None
            ),
            Err(PhillipsReadError::Limit)
        ));
        assert!(matches!(
            material.volume_bounds(
                Degree11MaterialVolumeOptions {
                    maximum_width_mm3: 1e-14,
                    ..default
                },
                None
            ),
            Err(PhillipsReadError::Unsupported)
        ));
    });
}

#[test]
fn actual_two_start_volume_and_owned_normal_cache_preserve_original_roots_within_caps() {
    let (single, frame) = super::super::super::tests::actual_m3_parts(true, false);
    let mut spec = single.thread().spec();
    spec.starts = 2;
    let multi = crate::phillips_drive::PhillipsScrew::new(
        crate::helical_thread::Thread::new(spec).unwrap(),
        *single.profile(),
        single.spec(),
        single.cavity().clone(),
        100_000,
        &|| false,
    )
    .unwrap();
    let text =
        crate::step::write_phillips_screw(&multi, frame[0], frame[1], frame[2], 1e-7).unwrap();
    with_head(&text, &mut |head| {
        let occupied = head
            .degree11_thread_carrier(Default::default(), None)
            .unwrap()
            .certify_occupied(Default::default(), None)
            .unwrap();
        let material = head
            .closed_degree11_material(&occupied, Default::default(), None)
            .unwrap();
        let volume = material.volume_bounds(Default::default(), None).unwrap();
        let owned = volume.owned_snapshot(Default::default(), None).unwrap();
        assert_eq!(owned.thread_faces().len(), 8);
        assert_eq!(volume.normals.len(), volume.thread_cells());
        assert!(volume.work_used() <= WORK_LIMIT && owned.work_used() <= WORK_LIMIT);
        assert!(
            volume.accounted_peak_bytes() <= MAX_BYTES && owned.accounted_peak_bytes() <= MAX_BYTES
        );
        let mut index = 0;
        for (face, f) in owned.thread_faces().iter().enumerate() {
            let [nu, nv] = f.source().spans();
            for u in 0..nu {
                for v in 0..nv {
                    assert!(Arc::ptr_eq(
                        owned.thread_normal(face, [u, v]).unwrap(),
                        &volume.normals[index]
                    ));
                    index += 1;
                }
            }
        }
        assert_eq!(index, volume.thread_cells());
        eprintln!(
            "D11_CACHE_2START cells={index} volumework={} ownedwork={} cachebytes={} resident={} volumepeak={} ownedpeak={}",
            volume.work_used(),
            owned.work_used(),
            cached_bytes(&volume.normals, volume.normals.capacity()),
            owned.retained_budget_bytes(),
            volume.accounted_peak_bytes(),
            owned.accounted_peak_bytes()
        );
    });
}
