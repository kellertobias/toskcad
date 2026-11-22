use super::*;
use multicad_rust_kernel::{
    Circle, Rational,
    authored_flat_thread::{AuthoredThread, SourceOptions},
    cylinder::AxialCylinder,
    helical_thread::Spec,
};
fn material() -> Material {
    let t = AuthoredThread::new(
        Spec {
            major_radius_mm: 1.5,
            pitch_mm: 0.5,
            length_mm: 4.,
            depth_mm: 0.25,
            crest_width_mm: 0.125,
            lower_flank_deg: 18.,
            upper_flank_deg: 12.,
            starts: 1,
            left_hand: false,
            root_radius_mm: None,
        },
        SourceOptions {
            position_limit_mm: 1e-7,
            normal_angle_limit_rad: 1e-7,
        },
    )
    .unwrap();
    let cavity = Cavity::new(
        t,
        source::Intent {
            entry_depth_mm: Some(0.5),
            terminal: Some(source::Terminal {
                width_mm: 0.5,
                radial_depth_mm: 0.125,
            }),
        },
    )
    .unwrap();
    let q = |n, d| Rational::new(n, d).unwrap();
    Material::new(
        AxialCylinder::new(
            7,
            Circle::new([q(3, 7), q(-5, 11)], q(5, 2)).unwrap(),
            q(7, 13),
            q(4, 1),
        )
        .unwrap(),
        cavity,
        0.2,
    )
    .unwrap()
}
fn options() -> FacePropertyOptions {
    FacePropertyOptions::default()
}
#[test]
fn actual_identity_and_planar_strain_scale_area_without_normalizing() {
    let m = material();
    let policy = GeometryPolicy::default();
    let f = FaceKey::StartCap;
    let a = transport(
        &m,
        Placement::from(PlaneFrame::default()),
        f,
        options(),
        policy,
        None,
    )
    .unwrap();
    assert_eq!(
        a.area_jacobian,
        source::Enclosure {
            lower: 1.,
            upper: 1.
        }
    );
    assert_eq!(a.local_centroid_reweighting_error_mm, [0.; 3]);
    let scale = 1. + 4e-11;
    let frame = PlaneFrame {
        u: [scale, 0., 0.],
        ..PlaneFrame::default()
    };
    frame.checked_metric_extent([2.5, 2.5, 4.], policy).unwrap();
    let b = transport(
        &m,
        Placement::from(frame),
        f,
        FacePropertyOptions {
            area_tolerance_mm2: 1e-11,
            ..options()
        },
        policy,
        None,
    )
    .unwrap();
    assert!(b.world.area_mm2.lower > a.world.area_mm2.upper);
    assert!(b.area_jacobian.lower <= scale && scale <= b.area_jacobian.upper);
    assert_eq!(b.local_centroid_reweighting_error_mm, [0.; 3]);
}

fn loose() -> FacePropertyOptions {
    FacePropertyOptions {
        area_tolerance_mm2: 1e-6,
        first_moment_tolerance_mm3: 1e-5,
        centroid_tolerance_mm: 1e-6,
        ..options()
    }
}
fn encloses(b: source::Enclosure, x: f64) {
    assert!(b.lower <= x && x <= b.upper, "{b:?} misses {x}");
}
#[test]
fn actual_all_roles_proper_transport_and_source_area_bounds() {
    let m = material();
    let p = GeometryPolicy::default();
    let frame = PlaneFrame {
        origin: [11., -7., 3.],
        u: [0., 1., 0.],
        v: [0., 0., 1.],
    };
    for index in 0..m.face_count() {
        let face = m.face_key(index).unwrap();
        let local = transport(
            &m,
            Placement::from(PlaneFrame::default()),
            face,
            loose(),
            p,
            None,
        )
        .unwrap();
        let world = transport(&m, Placement::from(frame), face, loose(), p, None).unwrap();
        assert_eq!(world.area_jacobian, B::ONE.source());
        assert!(
            area_upper(&m, face).unwrap() >= local.local.area_mm2.lower,
            "{face:?}: upper {} local {:?}",
            area_upper(&m, face).unwrap(),
            local.local.area_mm2
        );
        for (i, j) in [(0, 2), (1, 0), (2, 1)] {
            let expected = B::point(frame.origin[i]) + B::from(local.local.centroid_mm[j]);
            assert!(
                world.world.centroid_mm[i].lower <= expected.hi
                    && world.world.centroid_mm[i].upper >= expected.lo
            );
        }
    }
    let flip = PlaneFrame {
        u: [-1., 0., 0.],
        ..PlaneFrame::default()
    };
    assert!(
        transport(
            &m,
            Placement::from(flip),
            FaceKey::StartCap,
            loose(),
            p,
            None
        )
        .is_ok()
    );
    // The original stock's rational origin is ancestry, not an extra world shift.
    let cap = transport(
        &m,
        Placement::from(PlaneFrame::default()),
        FaceKey::StartCap,
        loose(),
        p,
        None,
    )
    .unwrap();
    encloses(cap.world.centroid_mm[0], 0.);
    encloses(cap.world.centroid_mm[1], 0.);
    encloses(cap.world.centroid_mm[2], 0.);
}
#[test]
fn nonplanar_weighting_actual_affine_metrics_and_precision_refusal() {
    let m = material();
    let p = GeometryPolicy::default();
    let frame = PlaneFrame {
        u: [1. + 4e-11, 0., 0.],
        ..PlaneFrame::default()
    };
    let result = transport(
        &m,
        Placement::from(frame),
        FaceKey::MaterialOuterWall,
        loose(),
        p,
        None,
    )
    .unwrap();
    assert!(
        result
            .local_centroid_reweighting_error_mm
            .into_iter()
            .all(|x| x > 0.)
    );
    let a = frame.u[0];
    // Actual cylinder normals x and y give cofactor magnitudes a and a².
    encloses(result.area_jacobian, a);
    encloses(result.area_jacobian, a * a);
    assert!(
        transport(
            &m,
            Placement::from(frame),
            FaceKey::MaterialOuterWall,
            FacePropertyOptions {
                area_tolerance_mm2: 1e-9,
                ..loose()
            },
            p,
            None
        )
        .is_err()
    );
    for axes in [
        [[-1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        [[1.01, 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        [[0.; 3]; 3],
        [[f64::NAN, 0., 0.], [0., 1., 0.], [0., 0., 1.]],
    ] {
        assert!(
            transport(
                &m,
                Placement {
                    origin: [0.; 3],
                    axes
                },
                FaceKey::StartCap,
                loose(),
                p,
                None
            )
            .is_err()
        );
    }
    let far = PlaneFrame {
        origin: [1e6, -1e6, 1e6],
        ..PlaneFrame::default()
    };
    assert!(
        transport(
            &m,
            Placement::from(far),
            FaceKey::StartCap,
            FacePropertyOptions {
                first_moment_tolerance_mm3: 1e-3,
                ..loose()
            },
            p,
            None
        )
        .is_ok()
    );
    assert!(
        transport(
            &m,
            Placement::from(PlaneFrame {
                origin: [1e7, 0., 0.],
                ..far
            }),
            FaceKey::StartCap,
            loose(),
            p,
            None
        )
        .is_err()
    );
}
#[test]
fn invalid_options_actual_budget_and_delayed_cancellation() {
    use std::cell::Cell;
    let m = material();
    let f = FaceKey::Thread {
        start: 0,
        band: source::ThreadBandRole::Lower,
    };
    let run = |o, i| {
        transport(
            &m,
            Placement::from(PlaneFrame::default()),
            f,
            o,
            GeometryPolicy::default(),
            i,
        )
    };
    for o in [
        FacePropertyOptions {
            maximum_cells: 0,
            ..loose()
        },
        FacePropertyOptions {
            maximum_cells: 65537,
            ..loose()
        },
        FacePropertyOptions {
            centroid_tolerance_mm: f64::NAN,
            ..loose()
        },
    ] {
        assert!(run(o, None).is_err());
    }
    assert!(
        run(
            FacePropertyOptions {
                maximum_cells: 1,
                ..loose()
            },
            None
        )
        .is_err()
    );
    assert!(run(loose(), Some(&|| true)).is_err());
    let count = Cell::new(0);
    assert!(
        run(
            loose(),
            Some(&|| {
                count.set(count.get() + 1);
                count.get() >= 3
            })
        )
        .is_err()
    );
    assert!(count.get() >= 3);
    let count = Cell::new(0);
    let callback = || {
        count.set(count.get() + 1);
        false
    };
    transport(
        &m,
        Placement::from(PlaneFrame::default()),
        FaceKey::StartCap,
        loose(),
        GeometryPolicy::default(),
        Some(&callback),
    )
    .unwrap();
    let total = count.get();
    count.set(0);
    assert!(
        transport(
            &m,
            Placement::from(PlaneFrame::default()),
            FaceKey::StartCap,
            loose(),
            GeometryPolicy::default(),
            Some(&|| {
                count.set(count.get() + 1);
                count.get() == total
            })
        )
        .is_err()
    );
    assert!(
        transport(
            &m,
            Placement::from(PlaneFrame::default()),
            FaceKey::Thread {
                start: 99,
                band: source::ThreadBandRole::Root
            },
            loose(),
            GeometryPolicy::default(),
            None
        )
        .is_err()
    );
}
#[test]
fn native_typed_handoff_and_faces_use_actual_occurrence() {
    let p = GeometryPolicy::default();
    let q = |n| Rational::new(n, 1).unwrap();
    let frame = PlaneFrame {
        origin: [10., 20., 30.],
        u: [0.8, 0.6, 0.],
        v: [-0.36, 0.48, 0.8],
    };
    let stock = Solid::build_axial_cylinder(
        AxialCylinder::new(
            7,
            Circle::new([q(0); 2], Rational::new(5, 2).unwrap()).unwrap(),
            q(0),
            q(4),
        )
        .unwrap(),
        frame,
        None,
        p,
    )
    .unwrap();
    assert!(
        stock
            .rust_internal_thread_treatment_world_face_properties_with_interrupt(
                FaceKey::StartCap,
                loose(),
                None
            )
            .unwrap()
            .is_none()
    );
    let cut = InternalThreadTreatmentCut {
        profile: HelicalThread {
            major_radius_mm: 1.5,
            pitch_mm: 0.5,
            length_mm: 4.,
            depth_mm: 0.25,
            crest_width_mm: 0.125,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            radial_taper: 0.,
            starts: 1,
            left_hand: false,
        },
        intent: source::Intent {
            entry_depth_mm: Some(0.5),
            terminal: Some(source::Terminal {
                width_mm: 0.5,
                radial_depth_mm: 0.125,
            }),
        },
        minimum_wall_mm: 0.5,
        frame,
    };
    let solid = stock
        .rust_cut_authored_internal_thread_treatments_with_interrupt(
            cut,
            EntityId::from_uuid_bits(9),
            p,
            None,
        )
        .unwrap();
    let actual = solid
        .rust_internal_thread_treatment_world_face_properties_with_interrupt(
            FaceKey::StartCap,
            loose(),
            None,
        )
        .unwrap()
        .unwrap();
    for (b, x) in actual.world.centroid_mm.into_iter().zip(frame.origin) {
        encloses(b, x);
    }
    assert!(solid.faces().unwrap().iter().all(|face| face.area_mm2 > 0.));
    assert!(
        solid
            .rust_internal_thread_treatment_world_face_properties_with_interrupt(
                FaceKey::StartCap,
                loose(),
                Some(&|| true)
            )
            .is_err()
    );
    let invalid = FacePropertyOptions {
        area_tolerance_mm2: 0.,
        ..loose()
    };
    assert!(
        solid
            .rust_internal_thread_treatment_world_face_properties_with_interrupt(
                FaceKey::StartCap,
                invalid,
                None
            )
            .is_err()
    );
}
#[test]
fn scalar_validation_cardinal_and_subnormal_accuracy_guards() {
    assert!(!cardinal([[1., 1., 0.], [0., 0., 1.], [0., 1., 0.]]));
    assert!(!cardinal([[1., 0., 0.], [0.; 3], [0., 0., 1.]]));
    assert!(!cardinal([[1., 0., 0.], [1., 0., 0.], [0., 0., 1.]]));
    assert!(positive(B::ZERO).is_err());
    assert!(positive(B::new(f64::NAN, 1.)).is_err());
    assert!(
        local_options(
            FacePropertyOptions {
                area_tolerance_mm2: f64::from_bits(1),
                ..loose()
            },
            B::ONE,
            1.,
            1.,
            1.
        )
        .is_err()
    );
}
