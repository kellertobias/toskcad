use super::*;
use crate::rounded_thread::Spec;
use multicad_precision::PortableMath as _;
use num_rational::BigRational;

fn cutter(left_hand: bool) -> Thread {
    let pitch = 1.27;
    let a = 27.5_f64.to_radians();
    let rounding = pitch * a.portable_cos() / (12. * (1. - a.portable_sin()));
    Thread::new(Spec {
        major_radius_mm: 3.175,
        pitch_mm: pitch,
        length_mm: 3. * pitch,
        depth_mm: pitch / (3. * a.portable_tan()),
        crest_radius_mm: rounding,
        root_radius_mm: rounding,
        flank_angle_deg: 27.5,
        starts: 1,
        left_hand,
    })
    .unwrap()
}

#[test]
fn blind_stock_retains_original_rounded_law_and_exact_independent_floor() {
    for hand in [false, true] {
        let source = cutter(hand);
        let height = Rational::new(127, 25).unwrap();
        let radius = Rational::new(127, 20).unwrap();
        let blind = BlindBore::new(source, radius, height).unwrap();
        assert_eq!(blind.cutter(), &source);
        assert_eq!(blind.stock_height_mm(), height);
        assert_eq!(blind.stock_radius_mm(), radius);
        let exact_floor = BigRational::new(127.into(), 25.into())
            - BigRational::from_float(source.spec().length_mm).unwrap();
        let [lo, hi] = blind
            .floor_thickness_bounds_mm()
            .map(|x| BigRational::from_float(x).unwrap());
        assert!(lo > BigRational::from_integer(0.into()));
        assert!(lo <= exact_floor && exact_floor <= hi);
        // Fractional physical phase makes the floor rim genuinely noncircular.
        let radii = [0., 0.7, 1.5, 3.].map(|theta| source.radius(theta, source.spec().length_mm));
        assert!(radii.windows(2).any(|r| (r[1] - r[0]).abs() > 0.01));
        assert_eq!(blind.cutter().profile_breaks(), source.profile_breaks());
        for band in 0..5 {
            assert_eq!(
                blind.cutter().profile_circle(band),
                source.profile_circle(band)
            );
        }
    }
}

#[test]
fn whole_floor_boxes_require_strict_finite_stock_membership() {
    use BoxClassification::*;
    let b = BlindBore::new(
        cutter(false),
        Rational::new(6, 1).unwrap(),
        Rational::new(127, 25).unwrap(),
    )
    .unwrap();
    assert_eq!(
        b.classify_floor_box([[-1., 1.], [-1., 1.], [4., 5.]])
            .unwrap(),
        Inside
    );
    for bounds in [
        [[-1., 1.], [-1., 1.], [3., 4.]],
        [[-1., 1.], [-1., 1.], [0., 0.]],
        [[-1., 1.], [-1., 1.], [5.08, 5.08]],
        [[-7., 7.], [-7., 7.], [4., 5.]],
    ] {
        assert_eq!(b.classify_floor_box(bounds).unwrap(), UnresolvedBoundary);
    }
    for bounds in [
        [[7., 8.], [-1., 1.], [4., 5.]],
        [[0., 0.], [0., 0.], [6., 7.]],
        [[0., 0.], [0., 0.], [-2., -1.]],
    ] {
        assert_eq!(b.classify_floor_box(bounds).unwrap(), Outside);
    }
    assert_eq!(
        b.classify_floor_box([[f64::NAN, 0.], [0., 0.], [4., 5.]]),
        Err(Error::Invalid)
    );
}

#[test]
fn through_overrun_and_unproved_wall_or_floor_do_not_become_blind_sources() {
    let thread = cutter(false);
    let radius = Rational::new(6, 1).unwrap();
    for h in [
        Rational::new(3, 1).unwrap(),
        Rational::new(381, 100).unwrap(),
        Rational::new(-5, 1).unwrap(),
    ] {
        assert!(BlindBore::new(thread, radius, h).is_err());
    }
    assert!(
        BlindBore::new(
            thread,
            Rational::new(3, 1).unwrap(),
            Rational::new(6, 1).unwrap()
        )
        .is_err()
    );
    assert!(BlindBore::new(thread, radius, Rational::new(10_000_001, 1).unwrap()).is_err());
    let good = BlindBore::new(thread, radius, Rational::new(127, 25).unwrap()).unwrap();
    assert_eq!(good.cutter(), &thread);
    // The existing through contract remains exactly its cutter's own extent.
    let through = super::super::ThroughBore::new(thread, 6.).unwrap();
    assert_eq!(through.length_mm(), thread.spec().length_mm);
    assert_ne!(through.length_mm(), good.stock_height_mm().as_f64());
}

#[test]
fn whole_meridian_support_encloses_both_exact_affine_endpoint_laws() {
    let q = |v| BigRational::from_float(v).unwrap();
    for hand in [false, true] {
        let thread = cutter(hand);
        let blind = BlindBore::new(
            thread,
            Rational::new(6, 1).unwrap(),
            Rational::new(6, 1).unwrap(),
        )
        .unwrap();
        let [lo, hi] = blind.cutter_radial_bounds_mm().map(q);
        for band in [1, 3] {
            let (anchor, radius, slope) = if band == 1 {
                (thread.breaks[1], thread.tangent_radius[0], thread.slope)
            } else {
                (thread.breaks[3], thread.tangent_radius[1], -thread.slope)
            };
            for phase in [thread.breaks[band], thread.breaks[band + 1]] {
                let endpoint = q(radius) + q(slope) * (q(phase) - q(anchor));
                assert!(lo <= endpoint && endpoint <= hi);
            }
        }
        for band in [0, 2, 4] {
            let c = thread.profile_circle(band).unwrap();
            // The retained support includes symbolic centre uncertainty,
            // rather than treating the cached binary centre as exact source.
            assert!(c.center_error_mm >= 0.);
            for i in 0..=32 {
                let phase = thread.breaks[band]
                    + (thread.breaks[band + 1] - thread.breaks[band]) * i as f64 / 32.;
                let radius = q(thread.radius_at_phase(phase));
                assert!(lo <= radius && radius <= hi);
            }
        }
        assert!(blind.cutter_radial_bounds_mm()[1] >= thread.spec().major_radius_mm);
        assert!(blind.stock_radius_bounds_mm()[0] > blind.cutter_radial_bounds_mm()[1]);
    }
}

#[test]
fn finite_material_support_encloses_independent_exact_rational_cylinder_extrema() {
    let q = |v| BigRational::from_float(v).unwrap();
    for hand in [false, true] {
        let b = BlindBore::new(
            cutter(hand),
            Rational::new(127, 20).unwrap(),
            Rational::new(127, 25).unwrap(),
        )
        .unwrap();
        for direction in [
            [1., 0., 0.],
            [0., 0., 1.],
            [0., 0., -1.],
            [0.4, -0.2, 0.7],
            [-2., 3., -1.],
            [0.; 3],
        ] {
            let [lo, hi] = b.support_bounds(direction, 1e-10).unwrap().map(q);
            let radial_squared = BigRational::new(127.into(), 20.into()).pow(2)
                * (q(direction[0]).pow(2) + q(direction[1]).pow(2));
            let axial = BigRational::new(127.into(), 25.into()) * q(direction[2].max(0.));
            let low = lo - axial.clone();
            let high = hi - axial;
            assert!(high >= BigRational::from_integer(0.into()));
            if low > BigRational::from_integer(0.into()) {
                assert!(low.pow(2) <= radial_squared);
            }
            assert!(radial_squared <= high.pow(2));
        }
        assert_eq!(
            b.support_bounds([f64::NAN, 0., 0.], 1e-10),
            Err(Error::Invalid)
        );
        assert_eq!(b.support_bounds([1., 0., 0.], 0.), Err(Error::Invalid));
        assert_eq!(
            b.support_bounds([1e308, 1e308, 1e308], 1e-10),
            Err(Error::Unresolved)
        );
        assert!(b.support_bounds([1., 0., 0.], 1e-300).is_err());
        assert!(b.support_bounds([1., 0., 0.], 1e-10).is_ok());
    }
}
