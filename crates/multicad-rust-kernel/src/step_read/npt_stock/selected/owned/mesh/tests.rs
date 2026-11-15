use super::*;
use crate::step_read::owned_insert::DocumentOptions;
pub(super) fn source(quarter: bool, mapped: bool) -> OwnedMaterial {
    let text = super::super::super::super::tests::fixture(quarter);
    let text = if mapped {
        super::super::super::tests::mapped(&text, "(-3.,0.,0.)", "1.")
    } else {
        text
    };
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    Selected::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap()
    .retain(Default::default(), None)
    .unwrap()
}
#[test]
fn actual_literal_mesh_probe() {
    for (quarter, mapped) in [(false, false), (true, false), (false, true)] {
        let source = source(quarter, mapped);
        let started = std::time::Instant::now();
        let result = source.mesh_candidate(Default::default(), None);
        eprintln!("NPT LITERAL MESH WALL {:?}", started.elapsed());
        eprintln!(
            "NPT LITERAL MESH quarter={quarter} mapped={mapped} {:?}",
            result.as_ref().map(|r| (
                r.mesh().vertices.len() / 3,
                r.mesh().face_indices.len(),
                r.position_error_mm(),
                r.normal_error_radians(),
                r.work_used(),
                r.accounted_peak_bytes(),
                r.original_root_triangle_counts()
            ))
        );
        let receipt = result.unwrap();
        assert!(std::sync::Arc::ptr_eq(&receipt.source().0, &source.0));
        assert!(receipt.position_error_mm() <= 0.1 && receipt.normal_error_radians() <= 0.3);
        assert_eq!(
            crate::manufacturing::certify(receipt.mesh()).unwrap(),
            receipt.connectivity()
        );
        assert!(
            receipt
                .mesh()
                .face_indices
                .iter()
                .all(|f| (*f as usize) < source.faces().len())
        );
        let proof_started = std::time::Instant::now();
        let qualified = receipt.certify_embedding(None).unwrap();
        eprintln!(
            "NPT LITERAL EMBEDDING quarter={quarter} mapped={mapped} time={:?} work={} peak={} radial={} height={}",
            proof_started.elapsed(),
            qualified.candidate().work_used(),
            qualified.candidate().accounted_peak_bytes(),
            qualified.embedding().minimum_radial_determinant(),
            qualified.embedding().minimum_axis_height_clearance_mm()
        );
        assert!(qualified.embedding().minimum_radial_determinant() > 0.);
        assert!(qualified.embedding().minimum_axis_height_clearance_mm() > 0.);
        assert_eq!(
            qualified.embedding().connectivity(),
            qualified.candidate().connectivity()
        );
        assert!(std::sync::Arc::ptr_eq(&source.0, &qualified.source().0));
    }
}
#[test]
fn actual_constant_angle_generators_have_whole_original_field_triangle_proofs() {
    let source = source(false, false);
    let face = &source.faces()[1];
    let Support::Tensor(t) = face.support else {
        panic!("fixture")
    };
    let tensor = &source.tensors()[t];
    let geometry = source.cell(t, [0, 0]).unwrap();
    let local = |uv: [f64; 2]| {
        geometry
            .point_bounds(uv.map(|x| [x, x]))
            .unwrap()
            .map(|b| b[0].midpoint(b[1]))
    };
    let a = local([0., 0.]);
    let delta = 0.1f64;
    let rays = [
        [a[0], a[1]],
        [
            delta.cos() * a[0] - delta.sin() * a[1],
            delta.sin() * a[0] + delta.cos() * a[1],
        ],
    ];
    let solve = |u: f64, ray: [f64; 2]| {
        let f = |v: f64| {
            let p = local([u, v]);
            (ray[0] * p[1] - ray[1] * p[0], ray[0] * p[0] + ray[1] * p[1])
        };
        let mut interval = None;
        let mut previous = 0.;
        let mut old = f(0.);
        if old.0.abs() < 1e-12 && old.1 > 0. {
            return 0.;
        }
        for j in 1..=128 {
            let v = j as f64 / 128.;
            let current = f(v);
            if old.0 * current.0 <= 0. && old.1 > 0. && current.1 > 0. {
                interval = Some((previous, v));
                break;
            }
            previous = v;
            old = current;
        }
        let (mut lo, mut hi) = interval.expect("actual local angular crossing");
        for _ in 0..48 {
            let mid = lo.midpoint(hi);
            if f(lo).0 * f(mid).0 <= 0. {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        lo.midpoint(hi)
    };
    let bytes = source.retained_budget_bytes() + SCRATCH;
    let mut budget = Budget {
        used: 0,
        remaining: WORK,
        bytes,
        peak: bytes,
        limit: BYTES,
        interrupt: None,
    };
    let root = normal::Normal::new(&source, tensor, 0, 0, &mut budget).unwrap();
    let mut qualified = false;
    for segments in [1usize, 2, 4, 8] {
        let mut max_position = 0f64;
        let mut max_angle = 0f64;
        for index in 0..segments {
            let u = [
                index as f64 / segments as f64,
                (index + 1) as f64 / segments as f64,
            ];
            let uv = [
                [u[0], solve(u[0], rays[0])],
                [u[1], solve(u[1], rays[0])],
                [u[1], solve(u[1], rays[1])],
                [u[0], solve(u[0], rays[1])],
            ];
            for indices in [[0usize, 3, 1], [1, 3, 2]] {
                let at = indices.map(|i| uv[i]);
                let intervals = at.map(|uv| {
                    source
                        .map_box_mm(geometry.point_bounds(uv.map(|x| [x, x])).unwrap())
                        .unwrap()
                        .map(interval)
                });
                let points = intervals.map(|p| p.map(I::midpoint));
                let rounding = intervals
                    .into_iter()
                    .zip(points)
                    .map(|(p, q)| cells::l1(std::array::from_fn(|k| p[k] - I::point(q[k]))))
                    .fold(0f64, f64::max);
                let quality = sheared::qualify(&root, at, points, rounding, &mut budget).unwrap();
                max_position = max_position.max(quality.position);
                max_angle = max_angle.max(quality.angle);
            }
        }
        eprintln!(
            "NPT SHEARED TRIANGLES root=[0,0] face={} generator_segments={segments} position={max_position} angle={max_angle} work={}",
            face.entity(),
            budget.used
        );
        if max_position <= 0.1 && max_angle <= 0.3 {
            qualified = true;
            break;
        }
    }
    assert!(
        qualified,
        "actual source must earn full-domain C0/normal correspondence"
    );
}

#[test]
fn actual_geometry_plan_reuse_and_mesh_refusals_preserve_the_literal_source() {
    let source = source(false, false);
    for options in [
        LiteralGeometryOptions {
            maximum_work: 0,
            ..Default::default()
        },
        LiteralGeometryOptions {
            maximum_work: 512_000_001,
            ..Default::default()
        },
        LiteralGeometryOptions {
            maximum_bytes: 0,
            ..Default::default()
        },
        LiteralGeometryOptions {
            maximum_bytes: BYTES + 1,
            ..Default::default()
        },
    ] {
        assert_eq!(
            source.geometry_plan(options, None).err(),
            Some(Error::Invalid)
        );
    }
    assert_eq!(
        source
            .geometry_plan(
                LiteralGeometryOptions {
                    maximum_work: 1,
                    ..Default::default()
                },
                None
            )
            .err(),
        Some(Error::Limit)
    );
    assert_eq!(
        source
            .geometry_plan(
                LiteralGeometryOptions {
                    maximum_bytes: 1,
                    ..Default::default()
                },
                None
            )
            .err(),
        Some(Error::Limit)
    );
    assert_eq!(
        source
            .geometry_plan(Default::default(), Some(&|| true))
            .err(),
        Some(Error::Cancelled)
    );
    let plan = source.geometry_plan(Default::default(), None).unwrap();
    assert!(std::sync::Arc::ptr_eq(&source.0, &plan.source().0));
    assert!(
        plan.preparation_work_used() > 0
            && plan.retained_budget_bytes() > source.retained_budget_bytes()
    );
    assert!(plan.construction_peak_bytes() >= plan.retained_budget_bytes());
    assert_eq!(
        plan.original_normal_coefficients(8, [0, 0]).err(),
        Some(Error::Invalid)
    );
    assert_eq!(
        plan.original_normal_coefficients(0, [2, 0]).err(),
        Some(Error::Invalid)
    );
    assert_eq!(
        plan.original_normal_coefficients(0, [0, 32]).err(),
        Some(Error::Invalid)
    );
    assert_eq!(
        plan.original_normal_coefficients(0, [0, 0]).unwrap().len(),
        484
    );
    for options in [
        LiteralMeshOptions {
            maximum_work: 0,
            ..Default::default()
        },
        LiteralMeshOptions {
            maximum_work: WORK + 1,
            ..Default::default()
        },
        LiteralMeshOptions {
            maximum_bytes: 0,
            ..Default::default()
        },
        LiteralMeshOptions {
            maximum_bytes: BYTES + 1,
            ..Default::default()
        },
        LiteralMeshOptions {
            tessellation: Tessellation {
                angle: f64::NAN,
                ..LiteralMeshOptions::default().tessellation
            },
            ..Default::default()
        },
    ] {
        assert_eq!(
            source.mesh_candidate(options, None).err(),
            Some(Error::Invalid)
        );
        assert_eq!(
            plan.mesh_candidate(options, None).err(),
            Some(Error::Invalid)
        );
    }
    for options in [
        LiteralMeshOptions {
            maximum_work: 1,
            ..Default::default()
        },
        LiteralMeshOptions {
            maximum_work: 512_000_000,
            ..Default::default()
        },
        LiteralMeshOptions {
            maximum_bytes: plan.retained_budget_bytes() - 1,
            ..Default::default()
        },
        LiteralMeshOptions {
            tessellation: Tessellation {
                triangle_limit: 1,
                ..LiteralMeshOptions::default().tessellation
            },
            ..Default::default()
        },
    ] {
        assert_eq!(plan.mesh_candidate(options, None).err(), Some(Error::Limit));
    }
    assert_eq!(
        plan.mesh_candidate(Default::default(), Some(&|| true))
            .err(),
        Some(Error::Cancelled)
    );
    let calls = std::cell::Cell::new(0usize);
    let stop = || {
        calls.set(calls.get() + 1);
        calls.get() > 20
    };
    assert_eq!(
        plan.mesh_candidate(Default::default(), Some(&stop)).err(),
        Some(Error::Cancelled)
    );
    let first = plan.mesh_candidate(Default::default(), None).unwrap();
    let second = plan.mesh_candidate(Default::default(), None).unwrap();
    assert_eq!(first.mesh(), second.mesh());
    assert_eq!(first.work_used(), second.work_used());
    assert_eq!(first.work_used() + first.remaining_work(), WORK);
    assert_eq!(first.preparation_work_used(), plan.preparation_work_used());
    assert_eq!(
        first.preparation_work_used() + first.query_work_used(),
        first.work_used()
    );
    assert!(first.accounted_peak_bytes() >= plan.retained_budget_bytes());
    assert!(std::sync::Arc::ptr_eq(&source.0, &first.source().0));
}
#[test]
fn actual_tilted_literal_mesh() {
    let text = super::super::super::super::tests::fixture(false);
    let text = super::super::super::tests::mapped(&text, "(-3.,0.,0.)", "1.")
        .replace("DIRECTION('',(0.,0.,5.))", "DIRECTION('',(0.,-4.,3.))");
    let document = CompoundDocument::parse(&text, DocumentOptions::default(), None).unwrap();
    let source = Selected::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap()
    .retain(Default::default(), None)
    .unwrap();
    let receipt = source.mesh_candidate(Default::default(), None).unwrap();
    receipt.certify_embedding(None).unwrap();
}
#[test]
fn actual_mesh_wrapper_cannot_bypass_preparation_resources_or_cancellation() {
    let source = source(false, false);
    for options in [
        LiteralMeshOptions {
            maximum_work: 1,
            ..Default::default()
        },
        LiteralMeshOptions {
            maximum_bytes: 1,
            ..Default::default()
        },
    ] {
        assert_eq!(
            source.mesh_candidate(options, None).err(),
            Some(Error::Limit)
        );
    }
    assert_eq!(
        source
            .mesh_candidate(Default::default(), Some(&|| true))
            .err(),
        Some(Error::Cancelled)
    );
    assert!(source.sample_vertex(0, Default::default(), None).is_ok());
}
