use super::*;
fn file(left: bool) -> String {
    file_length(left, 4.)
}
fn file_length(left: bool, length: f64) -> String {
    use crate::{
        authored_flat_thread::{AuthoredThread, SourceOptions},
        helical_thread::Spec,
        step::insert::{
            CurveOptions, DomainOptions, DomainPlan, SurfacePlan, TopologyPlan,
            write_local_candidate,
        },
        wood_insert::{BoundaryOptions, Insert},
    };
    let spec = Spec {
        major_radius_mm: 5.,
        pitch_mm: 2.,
        length_mm: length,
        depth_mm: 1.,
        crest_width_mm: 0.4,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: left,
        root_radius_mm: None,
    };
    let options = SourceOptions {
        position_limit_mm: 1e-7,
        normal_angle_limit_rad: 1e-7,
    };
    let outer = AuthoredThread::new(spec, options).unwrap();
    let inner = AuthoredThread::new(
        Spec {
            major_radius_mm: 1.5,
            pitch_mm: 0.5,
            depth_mm: (3. - (2.459 + 2.599) / 2.) / 2.,
            crest_width_mm: 0.25 - (3. - (2.675 + 2.775) / 2.) / 3_f64.sqrt(),
            ..spec
        },
        options,
    )
    .unwrap();
    let body = Insert::new(outer, inner, 0.5, None)
        .unwrap()
        .boundary_with_interrupt(BoundaryOptions::default(), None)
        .unwrap();
    let d = DomainPlan::prepare(&body, DomainOptions::default(), None).unwrap();
    let s = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
    let t = TopologyPlan::prepare(&s, 16_000_000, 64 * 1024 * 1024, None).unwrap();
    write_local_candidate(&t, 16_000_000, 64 * 1024 * 1024, None).unwrap()
}
#[test]
fn actual_literal_insert_flux_matches_independent_exact_coefficient_integrals() {
    for left in [false, true] {
        let text = file(left);
        let parsed = parse(&text, Options::default(), None).unwrap();
        let bindings = parsed.support_bindings(Options::default(), None).unwrap();
        let flux = bindings
            .oriented_flux(FluxOptions::default(), None)
            .unwrap();
        let [lo, hi] = flux.oriented_mm3();
        // Independent exact-rational integration of the actual emitted decimal
        // tensor and cap coefficients; no nominal or authoring-source volume.
        let independent = 231.37317864323353;
        assert!(
            lo <= independent && independent <= hi,
            "actual file flux [{lo},{hi}] vs{independent}"
        );
        assert!(std::ptr::eq(flux.bindings(), &bindings));
        assert!(flux.work_used() <= FluxOptions::default().maximum_work);
        assert!(flux.query_budget_bytes() <= FluxOptions::default().maximum_bytes);
        eprintln!(
            "literal flux [{lo},{hi}], work{}, bytes{}",
            flux.work_used(),
            flux.query_budget_bytes()
        );
    }
}

fn rewrite_points(text: &str, scale: f64, z_shift: f64) -> String {
    use std::collections::BTreeSet;
    let p = parse(text, Options::default(), None).unwrap();
    let b = p.support_bindings(Options::default(), None).unwrap();
    let mut cap_uv = BTreeSet::<u32>::new();
    for index in b.cap_face_indices() {
        for bound in p.faces()[index].loops() {
            for coedge in bound.coedges() {
                let edge = p.edge(coedge.edge).unwrap();
                for uv in edge.pcurves() {
                    if uv.support_entity == p.faces()[index].surface_entity() {
                        cap_uv.extend(
                            p.curve(uv.curve_entity)
                                .unwrap()
                                .pole_entities()
                                .iter()
                                .copied(),
                        );
                    }
                }
            }
        }
    }
    let mut result = text.to_owned();
    for point in p.points() {
        let coords = point.coordinates();
        if coords.len() != 3 && !cap_uv.contains(&point.entity()) {
            continue;
        }
        let prefix = format!("#{}=CARTESIAN_POINT(", point.entity());
        let start = result.find(&prefix).unwrap();
        let end = start + result[start..].find(';').unwrap() + 1;
        let values = coords
            .iter()
            .enumerate()
            .map(|(i, x)| {
                format!(
                    "{:.16E}",
                    x.value() * scale + if i == 2 { z_shift } else { 0. }
                )
            })
            .collect::<Vec<_>>()
            .join(",");
        result.replace_range(
            start..end,
            &format!("#{}=CARTESIAN_POINT('',({values}));", point.entity()),
        );
    }
    result
}
fn read_flux(text: &str) -> [f64; 2] {
    let p = parse(text, Options::default(), None).unwrap();
    let b = p.support_bindings(Options::default(), None).unwrap();
    b.oriented_flux(FluxOptions::default(), None)
        .unwrap()
        .oriented_mm3()
}
#[test]
fn actual_file_translation_and_scaling_use_literal_geometry_without_authoring_headers() {
    let text = file(false);
    let original = read_flux(&text);
    let shifted = read_flux(&rewrite_points(&text, 1., -10.));
    assert!(shifted[0] <= original[1] && original[0] <= shifted[1]);
    let doubled = read_flux(&rewrite_points(&text, 2., 0.));
    assert!(doubled[0] <= 8. * original[1] && 8. * original[0] <= doubled[1]);
    eprintln!("actual shifted flux{shifted:?}, doubled{doubled:?}");
}
#[test]
fn exact_query_budgets_numeric_floor_and_cancelled_publication_keep_the_actual_packet() {
    use std::cell::Cell;
    let text = file(false);
    let p = parse(&text, Options::default(), None).unwrap();
    let b = p.support_bindings(Options::default(), None).unwrap();
    let actual = b.oriented_flux(FluxOptions::default(), None).unwrap();
    let expected = actual.oriented_mm3();
    for options in [
        FluxOptions {
            maximum_width_mm3: 0.,
            ..Default::default()
        },
        FluxOptions {
            maximum_width_mm3: f64::NAN,
            ..Default::default()
        },
        FluxOptions {
            maximum_work: 0,
            ..Default::default()
        },
        FluxOptions {
            maximum_work: WORK + 1,
            ..Default::default()
        },
        FluxOptions {
            maximum_bytes: 0,
            ..Default::default()
        },
        FluxOptions {
            maximum_bytes: BYTES + 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            b.oriented_flux(options, None),
            Err(Error::Invalid)
        ));
    }
    let exact = FluxOptions {
        maximum_work: actual.work_used(),
        maximum_bytes: actual.query_budget_bytes(),
        ..Default::default()
    };
    assert_eq!(
        b.oriented_flux(exact, None).unwrap().oriented_mm3(),
        expected
    );
    assert!(matches!(
        b.oriented_flux(
            FluxOptions {
                maximum_work: exact.maximum_work - 1,
                ..exact
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        b.oriented_flux(
            FluxOptions {
                maximum_work: b.work_used() - 1,
                ..exact
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        b.oriented_flux(
            FluxOptions {
                maximum_bytes: exact.maximum_bytes - 1,
                ..exact
            },
            None
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        b.oriented_flux(
            FluxOptions {
                maximum_width_mm3: f64::from_bits(1),
                ..exact
            },
            None
        ),
        Err(Error::Unsupported)
    ));
    let calls = Cell::new(0usize);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    b.oriented_flux(exact, Some(&count)).unwrap();
    let last = calls.get();
    for at in [1, last / 2, last] {
        calls.set(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() == at
        };
        assert!(matches!(
            b.oriented_flux(exact, Some(&cancel)),
            Err(Error::Cancelled)
        ));
    }
    assert_eq!(actual.oriented_mm3(), expected);
    assert_eq!(p.literal_input(), text);
}
#[test]
fn literal_coordinate_overflow_is_refused_without_a_nominal_volume_substitute() {
    let huge = rewrite_points(&file(false), 1e110, 0.);
    let p = parse(&huge, Options::default(), None).unwrap();
    let b = p.support_bindings(Options::default(), None).unwrap();
    assert!(matches!(
        b.oriented_flux(FluxOptions::default(), None),
        Err(Error::Unsupported)
    ));
}
#[test]
fn checked_literal_coordinate_and_polynomial_integral_primitives_refuse_invalid_domains() {
    let text = file(false);
    let parsed = parse(&text, Options::default(), None).unwrap();
    let uv = parsed
        .points()
        .iter()
        .find(|p| p.coordinates().len() == 2)
        .unwrap();
    let mut work = Ledger {
        used: 0,
        remaining: WORK,
        interrupt: None,
    };
    assert!(matches!(
        point(&parsed, uv.entity(), &mut work),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        point(&parsed, u32::MAX, &mut work),
        Err(Error::Invalid)
    ));
    // Finite coefficients can still overflow their outward accumulated
    // integral. Refuse at the numerical primitive rather than publish NaN/Inf.
    let mut polynomial = Poly::zero(0, 1);
    polynomial.c[0] = I::point(f64::MAX);
    polynomial.c[1] = I::point(f64::MAX);
    assert!(matches!(
        polynomial.integral(&mut work),
        Err(Error::Unsupported)
    ));
}
#[test]
fn actual_occupied_volume_requires_the_same_global_material_and_accumulates_query_resources() {
    use std::cell::Cell;
    let text = file(false);
    let p = parse(&text, Options::default(), None).unwrap();
    let b = p.support_bindings(Options::default(), None).unwrap();
    let material = b.material_bounds(Options::default(), None).unwrap();
    let result = material
        .occupied_volume(FluxOptions::default(), None)
        .unwrap();
    assert!(std::ptr::eq(result.material(), &material));
    let bounds = result.volume_mm3();
    assert!(bounds[0] > 0. && bounds[0] <= 231.37317864323353 && bounds[1] >= 231.37317864323353);
    let exact = FluxOptions {
        maximum_work: result.work_used(),
        maximum_bytes: result.query_budget_bytes(),
        ..Default::default()
    };
    assert_eq!(
        material.occupied_volume(exact, None).unwrap().volume_mm3(),
        bounds
    );
    for options in [
        FluxOptions {
            maximum_work: 1,
            ..exact
        },
        FluxOptions {
            maximum_bytes: 1,
            ..exact
        },
        FluxOptions {
            maximum_work: exact.maximum_work - 1,
            ..exact
        },
        FluxOptions {
            maximum_bytes: exact.maximum_bytes - 1,
            ..exact
        },
    ] {
        assert!(matches!(
            material.occupied_volume(options, None),
            Err(Error::Limit)
        ));
    }
    assert!(matches!(
        material.occupied_volume(
            FluxOptions {
                maximum_work: 0,
                ..exact
            },
            None
        ),
        Err(Error::Invalid)
    ));
    let calls = Cell::new(0);
    let counter = || {
        calls.set(calls.get() + 1);
        false
    };
    material.occupied_volume(exact, Some(&counter)).unwrap();
    let last = calls.get();
    for at in [1, last] {
        calls.set(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() == at
        };
        assert!(matches!(
            material.occupied_volume(exact, Some(&cancel)),
            Err(Error::Cancelled)
        ));
    }
    assert_eq!(result.volume_mm3(), bounds);
}
#[test]
fn actual_far_axial_packet_refuses_uncertifiable_requested_mass_precision() {
    let text = rewrite_points(&file(false), 1., 1e12);
    let p = parse(&text, Options::default(), None).unwrap();
    let b = p.support_bindings(Options::default(), None).unwrap();
    let material = b.material_bounds(Options::default(), None).unwrap();
    assert!(matches!(
        material.occupied_volume(FluxOptions::default(), None),
        Err(Error::Unsupported)
    ));
    let broad = material
        .occupied_volume(
            FluxOptions {
                maximum_width_mm3: 10.,
                ..Default::default()
            },
            None,
        )
        .unwrap();
    let [lo, hi] = broad.volume_mm3();
    assert!(lo > 0. && hi - lo > 1. && hi - lo <= 10.);
    assert_eq!(p.literal_input(), text);
}
