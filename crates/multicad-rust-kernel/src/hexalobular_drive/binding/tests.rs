use super::*;
use crate::screw_recess::Recess;

fn parts() -> (Drive, FiniteVertices, [Vec<CubicIntersectionSpan>; 12]) {
    let authored = Drive::from_authored_tangent(
        TangentSpec {
            major_diameter_mm: 2.8,
            minor_diameter_mm: 2.05,
            outer_lobe_radius_mm: 0.28,
        },
        Cone {
            tip_depth_mm: 1.1 + 1.4 / 3_f64.sqrt(),
            depth_per_radius: 1. / 3_f64.sqrt(),
        },
    )
    .unwrap();
    let source =
        Drive::from_retained(*authored.arcs(), authored.cone(), TransportError::ZERO).unwrap();
    let carriers = std::array::from_fn(|i| source.intersection_carrier(i, 1e-8, 10_000).unwrap());
    let mouth_mm = source.arcs().map(|a| [a.start_mm[0], a.start_mm[1], 0.]);
    let floor_mm = std::array::from_fn(|i| carriers[i][0].poles_mm[0]);
    let vertices = FiniteVertices {
        mouth_mm,
        floor_mm,
        apex_mm: [0., 0., source.cone().tip_depth_mm],
    };
    (source, vertices, carriers)
}
fn binding() -> FiniteBinding {
    let (source, vertices, carriers) = parts();
    FiniteBinding::new(source, vertices, carriers, 1e-7).unwrap()
}

#[test]
fn finite_binding_preserves_literal_controls_and_separate_error_volume_budgets() {
    let (source, vertices, carriers) = parts();
    let value = FiniteBinding::new(source.clone(), vertices, carriers.clone(), 1e-7).unwrap();
    assert_eq!(value.source(), &source);
    assert_eq!(value.vertices(), &vertices);
    assert_eq!(value.carriers(), &carriers);
    assert_eq!(
        value.source().boundary_source(),
        BoundarySource::RetainedFamily
    );
    assert!(value.nominal_binding_error_mm() > 0.);
    assert!(
        value.spatial_source_error_bound_mm()
            >= value.nominal_binding_error_mm() + source.spatial_source_error_bound_mm()
    );
    let trim = value.additional_volume_error_bound_mm3().unwrap();
    assert!(trim > 0.);
    let recess = Recess::conical_hexalobular_from_binding(value.clone(), 1e-5, 100_000).unwrap();
    let Recess::ConicalHexalobular(actual) = &recess else {
        unreachable!()
    };
    assert_eq!(actual.finite_binding(), Some(&value));
    assert_eq!(
        recess.source_error_mm(),
        value.spatial_source_error_bound_mm()
    );
    let base = source.volume_bounds_mm3(1e-7, 100_000).unwrap();
    let volume = recess.volume_bounds_mm3();
    assert!(volume[0] <= base[0] && base[1] <= volume[1]);
    assert!(volume[1] - volume[0] <= 1e-5);
    assert!(volume[1] - volume[0] >= 2. * trim);
    assert!(recess.nominal_volume_mm3() >= volume[0] && recess.nominal_volume_mm3() <= volume[1]);
    let cloned = recess.clone();
    assert_eq!(cloned, recess);
    assert_eq!(actual.source(), value.source());
}

#[test]
fn finite_binding_rejects_positive_gaps_unbound_vertices_and_changed_control_geometry() {
    let (source, vertices, carriers) = parts();
    for mode in 0..9 {
        let mut v = vertices;
        let mut c = carriers.clone();
        match mode {
            0 => v.mouth_mm[0][0] += 1e-6,
            1 => v.floor_mm[0] = v.mouth_mm[0],
            2 => v.apex_mm = v.floor_mm[0],
            3 => c[0].clear(),
            4 => c[0][0].parameters[0] = 0.1,
            5 => c[0].last_mut().unwrap().parameters[1] = 0.9,
            6 => c[0][0].poles_mm[0][0] += 1e-6,
            7 => c[0][0].parameters[1] += 1e-6,
            8 => c[0][0].poles_mm[3][0] += 1e-6,
            _ => unreachable!(),
        }
        assert_eq!(
            FiniteBinding::new(source.clone(), v, c, 1e-7),
            Err(Error::Invalid)
        );
    }
    let mut moved = vertices;
    moved.mouth_mm[0][2] = 1e-3;
    assert_eq!(
        FiniteBinding::new(source.clone(), moved, carriers.clone(), 1e-7),
        Err(Error::Unresolved)
    );
    let mut changed = carriers.clone();
    changed[0][0].poles_mm[1][2] += 1e-3;
    assert_eq!(
        FiniteBinding::new(source.clone(), vertices, changed, 1e-7),
        Err(Error::Unresolved)
    );
    let mut invalid = carriers.clone();
    invalid[0][0].poles_mm[1][0] = f64::NAN;
    assert_eq!(
        FiniteBinding::new(source.clone(), vertices, invalid, 1e-7),
        Err(Error::Invalid)
    );
    let mut huge = carriers.clone();
    huge[0][0].poles_mm[1][0] = 1e10;
    assert_eq!(
        FiniteBinding::new(source.clone(), vertices, huge, 1e-7),
        Err(Error::Limit)
    );
    assert!(FiniteBinding::new(source, vertices, carriers, 1e-7).is_ok());
}

#[test]
fn finite_binding_precision_work_cancellation_and_volume_exhaustion_fail_without_mutating_retry() {
    use std::cell::Cell;
    let (source, vertices, carriers) = parts();
    let snapshot = source.clone();
    for width in [0., -1., f64::NAN, f64::INFINITY] {
        assert_eq!(
            FiniteBinding::new(source.clone(), vertices, carriers.clone(), width),
            Err(Error::Invalid)
        );
    }
    let mut invalid = vertices;
    invalid.apex_mm[0] = f64::NAN;
    assert_eq!(
        FiniteBinding::new(source.clone(), invalid, carriers.clone(), 1e-7),
        Err(Error::Invalid)
    );
    let mut huge = vertices;
    huge.apex_mm[0] = 1e10;
    assert_eq!(
        FiniteBinding::new(source.clone(), huge, carriers.clone(), 1e-7),
        Err(Error::Limit)
    );
    assert_eq!(
        FiniteBinding::new_with_cancel(source.clone(), vertices, carriers.clone(), 1e-7, &|| true),
        Err(Error::Cancelled)
    );
    let calls = Cell::new(0);
    assert_eq!(
        FiniteBinding::new_with_cancel(source.clone(), vertices, carriers.clone(), 1e-7, &|| {
            calls.set(calls.get() + 1);
            calls.get() > 3
        }),
        Err(Error::Cancelled)
    );
    let bound = binding();
    for width in [0., -1., f64::NAN, f64::INFINITY] {
        assert_eq!(
            Recess::conical_hexalobular_from_binding(bound.clone(), width, 100_000),
            Err(Error::Invalid)
        );
    }
    assert_eq!(
        Recess::conical_hexalobular_from_binding(bound.clone(), f64::MIN_POSITIVE, 100_000),
        Err(Error::Unresolved)
    );
    assert_eq!(
        Recess::conical_hexalobular_from_binding(bound.clone(), 1e-5, 0),
        Err(Error::Limit)
    );
    assert_eq!(
        Recess::conical_hexalobular_from_binding_with_cancel(bound.clone(), 1e-5, 100_000, &|| {
            true
        }),
        Err(Error::Cancelled)
    );
    assert_eq!(source, snapshot);
    assert!(Recess::conical_hexalobular_from_binding(bound, 1e-5, 100_000).is_ok());
}

#[test]
fn finite_graph_rejects_retracing_and_multiple_winding_even_inside_a_large_position_budget() {
    let (source, vertices, mut outside) = parts();
    outside[0][0].poles_mm[1][0] *= -1.;
    outside[0][0].poles_mm[1][1] *= -1.;
    assert_eq!(
        FiniteBinding::new(source, vertices, outside, 1e6),
        Err(Error::Unresolved)
    );
    let (source, vertices, mut carriers) = parts();
    let a = carriers[0][0].poles_mm[0];
    let b = carriers[0][0].poles_mm[3];
    carriers[0][0].poles_mm[1] = std::array::from_fn(|i| {
        if i < 2 {
            100. * (0.01 * a[i] + 0.99 * b[i])
        } else {
            a[i]
        }
    });
    carriers[0][0].poles_mm[2] = std::array::from_fn(|i| {
        if i < 2 {
            100. * (0.99 * a[i] + 0.01 * b[i])
        } else {
            b[i]
        }
    });
    assert_eq!(
        FiniteBinding::new(source.clone(), vertices, carriers, 1e6),
        Err(Error::Unresolved)
    );
    let mut twice = vertices;
    twice.floor_mm = std::array::from_fn(|i| {
        let t = i as f64 * TAU / 6.;
        [t.cos(), t.sin(), 1.]
    });
    let curves = std::array::from_fn(|i| {
        let a = twice.floor_mm[i];
        let b = twice.floor_mm[(i + 1) % 12];
        vec![CubicIntersectionSpan {
            parameters: [0., 1.],
            poles_mm: [
                a,
                std::array::from_fn(|k| a[k] + (b[k] - a[k]) / 3.),
                std::array::from_fn(|k| a[k] + 2. * (b[k] - a[k]) / 3.),
                b,
            ],
        }]
    });
    assert_eq!(
        FiniteBinding::new(source, twice, curves, 1e6),
        Err(Error::Invalid)
    );
}

#[test]
fn a_large_retained_boundary_error_cannot_mint_a_positive_volume_certificate() {
    let (source, vertices, mut carriers) = parts();
    carriers[0][0].poles_mm[1][2] += 1.;
    let value = FiniteBinding::new(source, vertices, carriers, 2.).unwrap();
    assert_eq!(
        Recess::conical_hexalobular_from_binding(value, 1e6, 100_000),
        Err(Error::Unresolved)
    );
}

#[test]
fn retained_binding_recess_keeps_original_fitted_thread_and_a_closed_metric_mesh() {
    use crate::{
        display::Tessellation,
        head_profile::Head,
        helical_thread::{Spec, Thread},
        helical_thread_mesh::{profiled_screw_volume_error_bound, tessellate_profiled_screw},
        manufacturing::certify,
        profiled_screw::Screw,
    };
    let bound = binding();
    let recess = Recess::conical_hexalobular_from_binding(bound.clone(), 1e-5, 100_000).unwrap();
    let thread = Thread::new(Spec {
        major_radius_mm: 1.4635,
        pitch_mm: 0.5,
        length_mm: 3.,
        depth_mm: 0.29908756729740676,
        crest_width_mm: 0.07131009168581062,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: true,
        root_radius_mm: Some(0.07216878364870323),
    })
    .unwrap();
    let screw = Screw::new_with_recess(
        thread,
        1.5,
        Head::Crowned {
            radius_mm: 2.8,
            height_mm: 2.4,
            crown_radius_mm: 5.,
            top_land_radius_mm: 2.05,
            neck_radius_mm: 0.125,
        },
        0.,
        0.1,
        recess,
    )
    .unwrap();
    assert_eq!(*screw.thread(), thread);
    assert!(screw.has_head_thread_overlap());
    let settings = Tessellation {
        deflection: 0.02,
        angle: 0.25,
        triangle_limit: 500_000,
    };
    let mesh = tessellate_profiled_screw(&screw, settings).unwrap();
    let proof = certify(&mesh).unwrap();
    assert!(
        (proof.signed_volume_mm3 - screw.volume_mm3()).abs()
            <= profiled_screw_volume_error_bound(&screw, settings).unwrap()
    );
    let Recess::ConicalHexalobular(recess) = screw.recess() else {
        unreachable!()
    };
    assert_eq!(recess.finite_binding(), Some(&bound));
    assert!(!crate::step::supports_profiled_screw(&screw, [0.; 3], 1e-7));
}
