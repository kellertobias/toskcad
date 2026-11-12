use super::*;
fn spec(k: f64) -> AuthoredSpec {
    let p = 25.4 / 27.;
    let m = 3_f64.sqrt();
    let depth = p / (1. / (m - k) + 1. / (m + k)) - 2. * 0.0645 * p;
    AuthoredSpec {
        major_radius_mm: (0.37360 * 25.4 - if k > 0. { 0.1615 * 25.4 / 16. } else { 0. } + depth)
            / 2.,
        pitch_mm: p,
        length_mm: 0.2639 * 25.4,
        radial_taper: k,
        section: AuthoredSection::NptFlat {
            depth_mm: depth,
            flat_width_mm: (p - depth * (1. / (m - k) + 1. / (m + k))) / 2.,
        },
        crest_radius_limit_mm: None,
    }
}
pub(super) fn source(k: f64) -> AuthoredThread {
    AuthoredThread::new(spec(k), 1e-7, 1e-6).unwrap()
}
#[test]
fn authored_npt_endpoints_share_expressions_without_changing_literal_sections() {
    for k in [-1. / 32., 1. / 32.] {
        let t = source(k);
        assert_eq!(t.band_count(), 3);
        for i in 0..2 {
            for z in [0., t.spec.length_mm] {
                assert_eq!(
                    t.surface_point_bounds(i, [1.; 2], [z; 2]).unwrap(),
                    t.surface_point_bounds(i + 1, [0.; 2], [z; 2]).unwrap()
                );
            }
        }
        let root = t.endpoint(0, false).unwrap();
        let end = t.endpoint(2, true).unwrap();
        assert_eq!(bounds(root.1), bounds(end.1));
        assert!(t.radial_bounds([0.; 2], [1.; 2]).unwrap()[0] <= t.spec.major_radius_mm + k);
        assert_eq!(t.derived_rounding_radius_bounds_mm().unwrap(), None);
    }
}
#[test]
fn actual_authored_npt_raw_material_mass_and_default_mesh_qualify() {
    let t = source(1. / 32.);
    let started = std::time::Instant::now();
    let raw = t.raw_boundary_ownership(1e-7, None).unwrap();
    println!(
        "raw source bound={}, time={:?}, faces={}, controls={}",
        raw.source_error_mm,
        started.elapsed(),
        raw.ownership.input().faces.len(),
        raw.ownership
            .input()
            .surfaces
            .iter()
            .map(|s| match s {
                crate::radial_tensor::Surface::Tensor(t) => t.poles().len(),
                _ => 0,
            })
            .sum::<usize>()
    );
    assert_eq!(raw.ownership.input().faces.len(), 6);
    let proof_started = std::time::Instant::now();
    let bound = raw.ownership.bind_surfaces(1e-7, None).unwrap();
    println!("bound surfaces");
    let axial = bound.certify_axial(1e-7, 1e-6, None).unwrap();
    println!("axial");
    let caps = axial.certify_caps(None).unwrap();
    println!("caps");
    let ordered = caps.certify_ordering(None).unwrap();
    println!("ordered");
    let material = ordered.certify_material(1e-7, None).unwrap();
    println!("material proof {:?}", proof_started.elapsed());
    let mass_started = std::time::Instant::now();
    let properties = material.properties_with_interrupt(1e-5, None).unwrap();
    let s = t.spec();
    let d = t.depth();
    let m = 3_f64.sqrt();
    let AuthoredSection::NptFlat {
        flat_width_mm: w, ..
    } = s.section
    else {
        panic!()
    };
    let slopes = 1. / (m - s.radial_taper) + 1. / (m + s.radial_taper);
    let root_width = s.pitch_mm - w - d * slopes;
    let area = |z: f64| {
        let r = s.major_radius_mm + s.radial_taper * z;
        let root = r - d;
        std::f64::consts::PI / s.pitch_mm
            * (w * r * r + root_width * root * root + slopes * (r.powi(3) - root.powi(3)) / 3.)
    };
    let half = s.length_mm / 2.;
    let expected = half * (area(half - half / 3_f64.sqrt()) + area(half + half / 3_f64.sqrt()));
    assert!(properties.volume_mm3[0] <= expected && expected <= properties.volume_mm3[1]);
    println!(
        "mass {:?}, independent {expected}, elapsed {:?}",
        properties.volume_mm3,
        mass_started.elapsed()
    );
    let mesh_started = std::time::Instant::now();
    let checks = std::cell::Cell::new(0_u64);
    let count = || {
        checks.set(checks.get() + 1);
        false
    };
    let triangle_limit = 1_000_000;
    println!("mesh begin actual app quality .1/.3, triangle limit={triangle_limit}");
    let result = material.tessellate_with_interrupt(
        crate::display::Tessellation {
            deflection: 0.1,
            angle: 0.3,
            triangle_limit,
        },
        Some(&count),
    );
    let mesh = result.unwrap_or_else(|error| {
        panic!(
            "actual default mesh {error:?}, elapsed {:?}, checks={}",
            mesh_started.elapsed(),
            checks.get()
        )
    });
    let proof = crate::manufacturing::certify(&mesh).unwrap();
    assert!(proof.signed_volume_mm3 > 0.);
    println!(
        "actual default mesh facets={}, volume={}, elapsed {:?}, checks={}",
        mesh.indices.len() / 3,
        proof.signed_volume_mm3,
        mesh_started.elapsed(),
        checks.get()
    );
}

#[test]
fn authored_negative_taper_uses_actual_reverse_cone_without_family_alias() {
    let t = source(-1. / 32.);
    let raw = t.raw_boundary_ownership(1e-7, None).unwrap();
    let material = raw
        .ownership
        .bind_surfaces(1e-7, None)
        .unwrap()
        .certify_axial(1e-7, 1e-6, None)
        .unwrap()
        .certify_caps(None)
        .unwrap()
        .certify_ordering(None)
        .unwrap()
        .certify_material(1e-7, None)
        .unwrap();
    assert!(material.boundary_error_mm() > 0.);
    assert_eq!(t.spec(), spec(-1. / 32.));
}

#[test]
fn authored_bspt_has_derived_circle_supports_and_shared_actual_joins() {
    for k in [-1. / 32., 1. / 32.] {
        let p = 25.4 / 28.;
        let alpha = 27.5_f64.to_radians();
        let depth = p / (3. * alpha.tan());
        let s = AuthoredSpec {
            major_radius_mm: 9.728 / 2. - k * if k > 0. { 4. } else { p / 2. },
            pitch_mm: p,
            length_mm: 6.5,
            radial_taper: k,
            section: AuthoredSection::BsptRounded { depth_mm: depth },
            crest_radius_limit_mm: None,
        };
        let t = AuthoredThread::new(s, 1e-7, 1e-6).unwrap();
        let r = t.derived_rounding_radius_bounds_mm().unwrap().unwrap();
        let m = 1. / alpha.tan();
        let expected = (p * (m * m - k * k) / (2. * m) - depth)
            / (2. * (1. / alpha.sin() - (1. + k * k).sqrt()));
        assert!(r[0] <= expected && expected <= r[1]);
        assert_eq!(t.tolerance_mm(), 1e-7);
        assert_eq!(t.minimum_feature_mm(), 1e-6);
        assert!(t.endpoint(5, false).is_err());
        // Query every actual finite arc/flank, including both sides of each
        // tangent. The independent scalar formulas use the canonical pitch
        // and depth; a rounded radius is never fed back as a physical input.
        let h = expected * alpha.cos();
        let valley = p * (m * m - k * k) / (2. * m * (m + k));
        let center = s.major_radius_mm - expected * (1. + k * k).sqrt();
        let tangent = center + expected * alpha.sin();
        for n in 0..=200 {
            let q = p * n as f64 / 200.;
            let a = if q <= h {
                center + (expected * expected - q * q).sqrt() - k * q
            } else if q < valley - h {
                tangent - k * h - (m + k) * (q - h)
            } else if q <= valley + h {
                let x = q - valley;
                s.major_radius_mm - depth + expected * (1. + k * k).sqrt()
                    - (expected * expected - x * x).sqrt()
                    - k * x
            } else if q < p - h {
                tangent + k * h + (m - k) * (q - (p - h))
            } else {
                let x = q - p;
                center + (expected * expected - x * x).sqrt() - k * x
            } + k * 2.;
            let bound = t.radial_bounds([q; 2], [2.; 2]).unwrap();
            assert!(bound[0] <= a && a <= bound[1], "q={q}: {bound:?} vs {a}");
        }
        let whole = t.radial_bounds([0., p], [0., s.length_mm]).unwrap();
        assert!(whole[0] > 0. && whole[1] < 6.);
        for band in 0..4 {
            assert_eq!(
                t.surface_point_bounds(band, [1.; 2], [2.; 2]).unwrap(),
                t.surface_point_bounds(band + 1, [0.; 2], [2.; 2]).unwrap()
            );
        }
        let a = t.surface_point_bounds(0, [0.; 2], [2.; 2]).unwrap();
        let b = t.surface_point_bounds(4, [1.; 2], [2.; 2]).unwrap();
        for axis in 0..3 {
            assert!(a[axis][0] <= b[axis][1] && b[axis][0] <= a[axis][1]);
        }
        assert!(matches!(
            t.raw_boundary_ownership(1e-7, None),
            Err(Error::Unsupported)
        ));
    }
}

#[test]
fn retained_long_rails_are_identical_and_interior_fold_is_not_authority() {
    use crate::{
        radial_tensor::{BoundaryOwnership, Surface},
        tapered_thread::RetainedTensor,
    };
    let t = source(1. / 32.);
    let raw = t.raw_boundary_ownership(1e-7, None).unwrap();
    let original = raw.ownership.input().clone();
    for pair in raw.band_faces.windows(2) {
        let a = &original.surfaces[original.faces[pair[0]].surface];
        let b = &original.surfaces[original.faces[pair[1]].surface];
        let (Surface::Tensor(a), Surface::Tensor(b)) = (a, b) else {
            panic!()
        };
        let [wa, ha] = a.spans().map(|s| 9 * s + 1);
        let hb = 9 * b.spans()[1] + 1;
        assert_eq!(ha, hb);
        assert_eq!(&a.poles()[(wa - 1) * ha..], &b.poles()[..hb]);
    }
    let mut folded = original.clone();
    let Surface::Tensor(tensor) = &mut folded.surfaces[0] else {
        panic!()
    };
    let spans = tensor.spans();
    let [w, h] = spans.map(|s| 9 * s + 1);
    let mut poles = tensor.poles().to_vec();
    poles[(w / 2) * h + h / 2][0] *= -20.;
    poles[(w / 2) * h + h / 2][1] *= -20.;
    *tensor = RetainedTensor::new(spans, poles).unwrap();
    let owned = BoundaryOwnership::new(folded, 1e-7, None).unwrap();
    let result = owned
        .bind_surfaces(1e-7, None)
        .and_then(|b| b.certify_axial(1e-7, 1e-6, None));
    assert!(
        result.is_err(),
        "Actual finite interior fold cannot gain authority from matching counts/caps"
    );
    assert_eq!(raw.ownership.input(), &original);
    assert_eq!(t.spec(), spec(1. / 32.));
}

#[test]
fn source_refusals_bounds_clipping_and_cancelled_raw_builds_are_atomic() {
    let t = source(1. / 32.);
    let s = t.spec();
    for bad in [
        AuthoredSpec {
            major_radius_mm: f64::NAN,
            ..s
        },
        AuthoredSpec {
            radial_taper: 0.,
            ..s
        },
        AuthoredSpec {
            pitch_mm: 0.00001,
            ..s
        },
        AuthoredSpec {
            length_mm: 0.1,
            ..s
        },
        AuthoredSpec {
            length_mm: 1e8,
            ..s
        },
        AuthoredSpec {
            section: AuthoredSection::NptFlat {
                depth_mm: 0.5,
                flat_width_mm: 2.,
            },
            ..s
        },
        AuthoredSpec {
            section: AuthoredSection::NptFlat {
                depth_mm: 0.5,
                flat_width_mm: f64::NAN,
            },
            ..s
        },
        AuthoredSpec {
            crest_radius_limit_mm: Some(f64::NAN),
            ..s
        },
        AuthoredSpec {
            crest_radius_limit_mm: Some(1e8),
            ..s
        },
        AuthoredSpec {
            crest_radius_limit_mm: Some(1.),
            ..s
        },
    ] {
        assert!(AuthoredThread::new(bad, 1e-7, 1e-6).is_err());
    }
    assert!(AuthoredThread::new(s, 1e-30, 1e-6).is_err());
    assert!(AuthoredThread::new(s, 1e-7, 10.).is_err());
    let asymmetric = AuthoredSpec {
        major_radius_mm: 5.,
        pitch_mm: 1.,
        length_mm: 5.,
        radial_taper: -1. / 32.,
        section: AuthoredSection::BsptRounded {
            depth_mm: 0.303428437595858,
        },
        crest_radius_limit_mm: None,
    };
    assert_eq!(
        AuthoredThread::new(asymmetric, 1e-7, 1e-6),
        Err(Error::Unresolved)
    );
    for (band, u, z) in [
        (3, [0.; 2], [0.; 2]),
        (0, [-1.; 2], [0.; 2]),
        (0, [2.; 2], [0.; 2]),
        (0, [f64::NAN; 2], [0.; 2]),
        (0, [0.; 2], [-1.; 2]),
        (0, [0.; 2], [100.; 2]),
        (0, [0.; 2], [f64::NAN; 2]),
    ] {
        assert!(t.surface_point_bounds(band, u, z).is_err());
    }
    for (q, z) in [
        ([-1.; 2], [0.; 2]),
        ([2.; 2], [0.; 2]),
        ([f64::NAN; 2], [0.; 2]),
        ([0.; 2], [-1.; 2]),
        ([0.; 2], [100.; 2]),
        ([0.; 2], [f64::NAN; 2]),
    ] {
        assert!(t.radial_bounds(q, z).is_err());
    }
    let clipped = AuthoredThread::new(
        AuthoredSpec {
            crest_radius_limit_mm: Some(s.major_radius_mm + 0.1),
            ..s
        },
        1e-7,
        1e-6,
    )
    .unwrap();
    assert!(matches!(
        clipped.raw_boundary_ownership(1e-7, None),
        Err(Error::Unsupported)
    ));
    assert!(
        clipped.radial_bounds([0.; 2], [s.length_mm; 2]).unwrap()[1] <= s.major_radius_mm + 0.1
    );
    assert!(matches!(
        t.raw_boundary_ownership(0., None),
        Err(Error::Invalid)
    ));
    assert!(t.raw_boundary_ownership(1e-30, None).is_err());
    for stop in [0, 2, 10, 100, 300] {
        let calls = std::cell::Cell::new(0);
        let cancel = || {
            let n = calls.get();
            calls.set(n + 1);
            n >= stop
        };
        assert!(matches!(
            t.raw_boundary_ownership(1e-7, Some(&cancel)),
            Err(Error::Cancelled)
        ));
    }
    assert_eq!(t.spec(), s);
    let retry = t.raw_boundary_ownership(1e-7, None).unwrap();
    assert_eq!(retry.band_faces, vec![2, 3, 4]);
}

#[test]
fn source_raw_extreme_scale_and_long_turn_budgets_refuse_without_mutation() {
    // These are valid authored source laws; a bounded raw representation is
    // allowed to refuse them, and cannot silently coarsen its source receipt.
    for (radius, turns) in [(1e6, 1.), (5., 100.), (1e4, 32.)] {
        let s = AuthoredSpec {
            major_radius_mm: radius,
            length_mm: spec(1. / 32.).pitch_mm * turns,
            ..spec(1. / 32.)
        };
        let t = AuthoredThread::new(s, 0.001, 1e-6).unwrap();
        let before = t;
        let result = t.raw_boundary_ownership(1e-7, None);
        assert!(matches!(result, Err(Error::Limit | Error::Unresolved)));
        assert_eq!(t, before);
        let point = t.surface_point_bounds(1, [0.5; 2], [0.; 2]).unwrap();
        assert!(point.into_iter().flatten().all(f64::is_finite));
    }
}
