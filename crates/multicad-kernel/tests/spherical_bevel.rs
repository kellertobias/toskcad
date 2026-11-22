use multicad_core::document::EntityId;
use multicad_kernel::{GeometryPolicy, Solid, SphericalBevelGear};
use multicad_precision::PortableMath as _;
use std::f64::consts::{FRAC_PI_2, PI};
fn specimen(teeth: u32, mate: u32) -> SphericalBevelGear {
    let pitch = (f64::from(teeth) / f64::from(mate)).portable_atan();
    let distance = f64::from(teeth) / (2. * pitch.portable_sin());
    let base = (pitch.portable_sin() * 20_f64.to_radians().portable_cos()).portable_asin();
    SphericalBevelGear {
        teeth,
        pitch_cone_angle_rad: pitch,
        base_cone_angle_rad: base,
        tip_cone_angle_rad: pitch + (1. / distance).portable_atan(),
        root_cone_angle_rad: pitch - (1.25 / distance).portable_atan(),
        outer_cone_distance_mm: distance,
        face_width_mm: (distance / 6.).min(5.),
        bore_radius_mm: 2.,
        half_base_tooth_angle_rad: (PI / 2. - 0.07) / f64::from(teeth) + involute(pitch, base),
        root_fillet_radius_mm: 0.2,
    }
}
fn involute(gamma: f64, base: f64) -> f64 {
    let sigma = (gamma.portable_cos() / base.portable_cos()).portable_acos();
    sigma / base.portable_sin()
        - (sigma.portable_sin()).portable_atan2(base.portable_sin() * sigma.portable_cos())
}
fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a} != {b}; budget {tol}");
}
#[test]
fn native_spherical_bevel_bootstrap_is_one_named_volumetric_solid() {
    let id = EntityId::new();
    let solid =
        Solid::rust_spherical_bevel_gear_named(specimen(24, 24), id, GeometryPolicy::default())
            .unwrap();
    assert_eq!(solid.solid_count().unwrap(), 1);
    assert!(solid.volume_mm3().unwrap() > 0.);
    let topology = solid.topology_snapshot(id, 0).unwrap();
    assert!(topology.elements.iter().all(|e| !e.origins.is_empty()));
    assert!(
        topology
            .elements
            .iter()
            .flat_map(|e| &e.origins)
            .any(|o| o.role.contains("bevel/heel"))
    );
    assert!(
        topology
            .elements
            .iter()
            .flat_map(|e| &e.origins)
            .any(|o| o.role.contains("bevel/toe"))
    );
}
#[test]
fn invalid_bevel_domains_and_toe_bore_wall_are_refused() {
    let p = GeometryPolicy::default();
    let good = specimen(24, 24);
    for invalid in [
        SphericalBevelGear { teeth: 65, ..good },
        SphericalBevelGear {
            tip_cone_angle_rad: FRAC_PI_2,
            ..good
        },
        SphericalBevelGear {
            face_width_mm: good.outer_cone_distance_mm,
            ..good
        },
        SphericalBevelGear {
            bore_radius_mm: good.outer_cone_distance_mm,
            ..good
        },
        SphericalBevelGear {
            root_fillet_radius_mm: f64::NAN,
            ..good
        },
    ] {
        assert!(Solid::rust_spherical_bevel_gear_named(invalid, EntityId::new(), p).is_err());
    }
    assert!(Solid::rust_spherical_bevel_gear_named(good, EntityId::nil(), p).is_err());
    close(
        involute(good.base_cone_angle_rad, good.base_cone_angle_rad),
        0.,
        1e-12,
    );
}

struct RootOracle {
    rounding: f64,
    center_gamma: f64,
    center_phi: f64,
    touch: f64,
    meridian: bool,
}
fn root_oracle(s: SphericalBevelGear) -> RootOracle {
    let rho = (s.root_fillet_radius_mm / s.outer_cone_distance_mm).portable_asin();
    let center = s.root_cone_angle_rad + rho;
    let meridian =
        center <= (s.base_cone_angle_rad.portable_cos() * rho.portable_cos()).portable_acos();
    let (touch, phi) = if meridian {
        (
            (center.portable_cos() / rho.portable_cos()).portable_acos(),
            s.half_base_tooth_angle_rad
                + (rho.portable_sin() / center.portable_sin()).portable_asin(),
        )
    } else {
        let sigma =
            (center.portable_cos() / s.base_cone_angle_rad.portable_cos()).portable_acos() - rho;
        let gamma = (s.base_cone_angle_rad.portable_cos() * sigma.portable_cos()).portable_acos();
        let half = s.half_base_tooth_angle_rad - involute(gamma, s.base_cone_angle_rad);
        (
            gamma,
            half + (rho.portable_sin() * s.base_cone_angle_rad.portable_sin()
                / (gamma.portable_sin() * center.portable_sin()))
            .portable_asin(),
        )
    };
    RootOracle {
        rounding: rho,
        center_gamma: center,
        center_phi: phi,
        touch,
        meridian,
    }
}
fn half_tooth(s: SphericalBevelGear, r: &RootOracle, gamma: f64) -> f64 {
    if gamma <= s.root_cone_angle_rad {
        return r.center_phi;
    }
    if gamma < r.touch {
        let k = (r.rounding.portable_cos() - gamma.portable_cos() * r.center_gamma.portable_cos())
            / (gamma.portable_sin() * r.center_gamma.portable_sin());
        // Exactly tangent endpoints can round one ULP outside acos's domain.
        assert!(k.abs() <= 1. + 1e-12);
        return r.center_phi - k.clamp(-1., 1.).portable_acos();
    }
    if gamma < s.base_cone_angle_rad {
        s.half_base_tooth_angle_rad
    } else {
        s.half_base_tooth_angle_rad - involute(gamma, s.base_cone_angle_rad)
    }
}
fn integral(f: &impl Fn(f64) -> f64, a: f64, b: f64, eps: f64) -> f64 {
    fn walk(
        f: &impl Fn(f64) -> f64,
        a: f64,
        b: f64,
        samples: [f64; 3],
        whole: f64,
        eps: f64,
        n: u32,
    ) -> f64 {
        let [fa, fm, fb] = samples;
        let m = (a + b) / 2.;
        let l = (a + m) / 2.;
        let r = (m + b) / 2.;
        let fl = f(l);
        let fr = f(r);
        let left = (m - a) * (fa + 4. * fl + fm) / 6.;
        let right = (b - m) * (fm + 4. * fr + fb) / 6.;
        let delta = left + right - whole;
        if delta.abs() < 15. * eps {
            return left + right + delta / 15.;
        }
        assert!(n > 0, "independent quadrature did not converge");
        walk(f, a, m, [fa, fl, fm], left, eps / 2., n - 1)
            + walk(f, m, b, [fm, fr, fb], right, eps / 2., n - 1)
    }
    if a == b {
        return 0.;
    }
    let m = (a + b) / 2.;
    let fa = f(a);
    let fm = f(m);
    let fb = f(b);
    walk(
        f,
        a,
        b,
        [fa, fm, fb],
        (b - a) * (fa + 4. * fm + fb) / 6.,
        eps,
        28,
    )
}
fn independent_volume(s: SphericalBevelGear) -> f64 {
    let r = root_oracle(s);
    let d = s.pitch_cone_angle_rad;
    let density = |g: f64| g.portable_sin() / (g - d).portable_cos().powi(3);
    let mut angular = 2. * PI * integral(&density, 0., s.root_cone_angle_rad, 1e-13);
    // Squaring the coordinate removes the square-root tangency at the root.
    let root_span = r.touch - s.root_cone_angle_rad;
    angular += 2.
        * f64::from(s.teeth)
        * integral(
            &|u| {
                if u == 0. {
                    return 0.;
                }
                let gamma = s.root_cone_angle_rad + root_span * u * u;
                half_tooth(s, &r, gamma) * density(gamma) * 2. * root_span * u
            },
            0.,
            1.,
            1e-13,
        );
    let start = if r.meridian {
        angular += 2.
            * f64::from(s.teeth)
            * s.half_base_tooth_angle_rad
            * integral(&density, r.touch, s.base_cone_angle_rad, 1e-13);
        s.base_cone_angle_rad
    } else {
        r.touch
    };
    let span = s.tip_cone_angle_rad - start;
    angular += 2.
        * f64::from(s.teeth)
        * integral(
            &|u| {
                let gamma = start + span * u * u;
                half_tooth(s, &r, gamma) * density(gamma) * 2. * span * u
            },
            0.,
            1.,
            1e-13,
        );
    let radial = (s.outer_cone_distance_mm.powi(3)
        - (s.outer_cone_distance_mm - s.face_width_mm).powi(3))
        / 3.;
    radial * angular - PI * s.bore_radius_mm.powi(2) * s.face_width_mm / d.portable_cos()
}
fn position(s: SphericalBevelGear, distance: f64, gamma: f64, phi: f64) -> [f64; 3] {
    let radial = distance / (gamma - s.pitch_cone_angle_rad).portable_cos();
    [
        radial * gamma.portable_sin() * phi.portable_cos(),
        radial * gamma.portable_sin() * phi.portable_sin(),
        radial * gamma.portable_cos(),
    ]
}
fn material(solid: &Solid, p: [f64; 3]) -> f64 {
    use multicad_kernel::{BooleanOperation, RigidTransform};
    let size: f64 = 0.001;
    let cube = Solid::box_solid([size; 3], Default::default())
        .unwrap()
        .transformed(RigidTransform {
            translation_mm: p.map(|v| v - size / 2.),
            ..Default::default()
        })
        .unwrap();
    solid
        .boolean(&cube, BooleanOperation::Intersect, Default::default())
        .unwrap()
        .map_or(0., |s| s.volume_mm3().unwrap())
        / size.powi(3)
}
fn offset(p: [f64; 3], n: [f64; 3], a: f64) -> [f64; 3] {
    std::array::from_fn(|i| p[i] + a * n[i])
}
fn radial_side_normal(gamma: f64, phi: f64, slope: f64) -> [f64; 3] {
    let eg = [
        gamma.portable_cos() * phi.portable_cos(),
        gamma.portable_cos() * phi.portable_sin(),
        -gamma.portable_sin(),
    ];
    let ep = [-phi.portable_sin(), phi.portable_cos(), 0.];
    let n: [f64; 3] = std::array::from_fn(|i| ep[i] - gamma.portable_sin() * slope * eg[i]);
    let length = n.iter().map(|v| v * v).sum::<f64>().sqrt();
    n.map(|v| v / length)
}
fn verify_material(solid: &Solid, s: SphericalBevelGear) {
    let root = root_oracle(s);
    let gamma = (root.touch.max(s.base_cone_angle_rad) + s.tip_cone_angle_rad) / 2.;
    let half = half_tooth(s, &root, gamma);
    let slope = -(s.base_cone_angle_rad.portable_cos().powi(2) - gamma.portable_cos().powi(2))
        .sqrt()
        / (s.base_cone_angle_rad.portable_sin() * gamma.portable_sin());
    let normal = radial_side_normal(gamma, half, slope);
    for fraction in [0.2, 0.5, 0.8] {
        let distance = s.outer_cone_distance_mm - s.face_width_mm * fraction;
        let p = position(s, distance, gamma, half);
        assert!(
            material(solid, offset(p, normal, -0.003)) > 0.99,
            "missing spherical involute material"
        );
        assert!(
            material(solid, offset(p, normal, 0.003)) < 1e-5,
            "extra spherical involute material"
        );
        assert!(material(solid, position(s, distance, s.pitch_cone_angle_rad, 0.)) > 0.99);
        assert!(
            material(
                solid,
                position(s, distance, s.pitch_cone_angle_rad, PI / f64::from(s.teeth))
            ) < 1e-5
        );
    }
    // Both genuine backcones, away from the pitch circle so flat-Z caps fail.
    let gamma = (s.pitch_cone_angle_rad + s.tip_cone_angle_rad) / 2.;
    let normal = [
        s.pitch_cone_angle_rad.portable_sin(),
        0.,
        s.pitch_cone_angle_rad.portable_cos(),
    ];
    for (distance, direction) in [
        (s.outer_cone_distance_mm, 1.),
        (s.outer_cone_distance_mm - s.face_width_mm, -1.),
    ] {
        let cap = position(s, distance, gamma, 0.);
        assert!(
            material(solid, offset(cap, normal, -direction * 0.003)) > 0.99,
            "missing backcone cap interior"
        );
        assert!(
            material(solid, offset(cap, normal, direction * 0.003)) < 1e-5,
            "extra material beyond backcone cap"
        );
    }
    // Three independently intersected spherical latitudes constrain root-arc
    // curvature, at three different cone depths rather than one tangent point.
    for (arc_fraction, cone_fraction) in [(0.25, 0.2), (0.5, 0.5), (0.75, 0.8)] {
        let gamma = s.root_cone_angle_rad + (root.touch - s.root_cone_angle_rad) * arc_fraction;
        let half = half_tooth(s, &root, gamma);
        let k = (root.rounding.portable_cos()
            - gamma.portable_cos() * root.center_gamma.portable_cos())
            / (gamma.portable_sin() * root.center_gamma.portable_sin());
        let kp = (root.center_gamma.portable_cos()
            - root.rounding.portable_cos() * gamma.portable_cos())
            / (gamma.portable_sin().powi(2) * root.center_gamma.portable_sin());
        let normal = radial_side_normal(gamma, half, kp / (1. - k * k).sqrt());
        let p = position(
            s,
            s.outer_cone_distance_mm - s.face_width_mm * cone_fraction,
            gamma,
            half,
        );
        assert!(
            material(solid, offset(p, normal, -0.003)) > 0.99,
            "missing root-circle material"
        );
        assert!(
            material(solid, offset(p, normal, 0.003)) < 1e-5,
            "extra root-circle material"
        );
    }
}
#[test]
fn both_root_branches_and_64_tooth_boundary_match_independent_volume_cap_flank_and_raw_step_oracles()
 {
    use multicad_core::references::TopologyKind;
    let policy = GeometryPolicy::default();
    for (teeth, mate) in [(24, 24), (64, 64), (24, 64), (64, 24)] {
        let s = specimen(teeth, mate);
        if teeth == mate {
            assert_eq!(root_oracle(s).meridian, teeth == 24);
        }
        let feature = EntityId::new();
        let solid = Solid::rust_spherical_bevel_gear_named(s, feature, policy).unwrap();
        let expected = independent_volume(s);
        let area: f64 = solid.faces().unwrap().iter().map(|f| f.area_mm2).sum();
        let budget = policy.linear_tolerance_mm * area + 1e-7 * expected;
        close(solid.volume_mm3().unwrap(), expected, budget);
        verify_material(&solid, s);
        let snapshot = solid.topology_snapshot(feature, 0).unwrap();
        let origins: std::collections::BTreeMap<_, _> = snapshot
            .elements
            .iter()
            .filter(|e| e.key.kind == TopologyKind::Edge)
            .map(|e| (e.key.ordinal, &e.origins))
            .collect();
        let mut checked = 0;
        for edge in solid.display_edges(0.001).unwrap() {
            let names = &origins[&edge.index];
            if names.iter().any(|o| o.role.contains("bevel/heel")) {
                for p in &edge.points_mm {
                    close(
                        p[0].portable_hypot(p[1]) * s.pitch_cone_angle_rad.portable_sin()
                            + p[2] * s.pitch_cone_angle_rad.portable_cos(),
                        s.outer_cone_distance_mm,
                        16. * policy.linear_tolerance_mm,
                    );
                }
            }
            if names
                .iter()
                .any(|o| o.role == "bevel/tooth/0/spherical_involute.right")
                && names.iter().any(|o| o.role.contains("bevel/heel"))
            {
                for p in &edge.points_mm {
                    let gamma = p[0].portable_hypot(p[1]).portable_atan2(p[2]);
                    if gamma > s.base_cone_angle_rad + 1e-5 && gamma < s.tip_cone_angle_rad - 1e-5 {
                        close(
                            p[1].portable_atan2(p[0]),
                            s.half_base_tooth_angle_rad - involute(gamma, s.base_cone_angle_rad),
                            16. * policy.linear_tolerance_mm / s.outer_cone_distance_mm,
                        );
                        checked += 1;
                    }
                }
            }
        }
        assert!(
            checked >= 2,
            "native cap must expose measured spherical flank points"
        );
        let directory = std::env::temp_dir().join(format!("multicad-bevel-{}", EntityId::new()));
        std::fs::create_dir(&directory).unwrap();
        let path = directory.join("gear.step");
        solid.write_step(&path).unwrap();
        let step = Solid::read_step(&path, policy).unwrap();
        assert_eq!(step.solid_count().unwrap(), 1);
        close(step.volume_mm3().unwrap(), expected, budget);
        verify_material(&step, s);
        for (a, b) in solid
            .bounds()
            .unwrap()
            .min
            .into_iter()
            .chain(solid.bounds().unwrap().max)
            .zip(
                step.bounds()
                    .unwrap()
                    .min
                    .into_iter()
                    .chain(step.bounds().unwrap().max),
            )
        {
            close(a, b, 16. * policy.linear_tolerance_mm);
        }
        std::fs::remove_dir_all(directory).unwrap();
    }
}

fn volume_and_raw_step(solid: &Solid, expected: f64) {
    let policy = GeometryPolicy::default();
    let area: f64 = solid
        .faces()
        .unwrap()
        .iter()
        .map(|face| face.area_mm2)
        .sum();
    let budget = policy.linear_tolerance_mm * area + 1e-7 * expected;
    close(solid.volume_mm3().unwrap(), expected, budget);
    let directory = std::env::temp_dir().join(format!("multicad-flux-{}", EntityId::new()));
    std::fs::create_dir(&directory).unwrap();
    let path = directory.join("shape.step");
    solid.write_step(&path).unwrap();
    let reopened = Solid::read_step(&path, policy).unwrap();
    assert_eq!(
        reopened.solid_count().unwrap(),
        solid.solid_count().unwrap()
    );
    close(reopened.volume_mm3().unwrap(), expected, budget);
    std::fs::remove_dir_all(directory).unwrap();
}
fn posed_and_scaled_volume(solid: &Solid, expected: f64) {
    use multicad_kernel::{AffineTransform, RigidTransform};
    let posed = solid
        .transformed(RigidTransform {
            rotation: [[0.36, -0.8, 0.48], [0.48, 0.6, 0.64], [-0.8, 0., 0.6]],
            translation_mm: [10000., -20000., 30000.],
        })
        .unwrap();
    volume_and_raw_step(&posed, expected);
    let scaled = solid
        .affine_transformed_named(
            AffineTransform {
                matrix: [[1.5, 0., 0.], [0., 1.5, 0.], [0., 0., 1.5]],
                translation_mm: [-20000., 30000., 10000.],
            },
            EntityId::new(),
            GeometryPolicy::default(),
        )
        .unwrap();
    volume_and_raw_step(&scaled, expected * 1.5_f64.powi(3));
    close(solid.volume_mm3().unwrap(), expected, 1e-7 * expected);
}
#[test]
fn zero_bore_asymmetric_bevels_keep_independent_volume_after_rigid_scale_and_raw_step() {
    for (teeth, mate) in [(24, 32), (24, 64), (64, 24)] {
        let mut s = specimen(teeth, mate);
        s.face_width_mm = 2.;
        s.bore_radius_mm = 0.;
        let solid =
            Solid::rust_spherical_bevel_gear_named(s, EntityId::new(), GeometryPolicy::default())
                .unwrap();
        let expected = independent_volume(s);
        volume_and_raw_step(&solid, expected);
        verify_material(&solid, s);
        posed_and_scaled_volume(&solid, expected);
    }
}
#[test]
fn directional_native_flux_preserves_known_box_cone_and_holed_volumes_and_raw_step() {
    use multicad_kernel::{
        BooleanOperation, CurvedEdge2d, CurvedProfile, PlaneFrame, RevolutionAxis,
    };
    let policy = GeometryPolicy::default();
    let box_solid = Solid::box_solid([3., 4., 5.], policy).unwrap();
    let cylinder = Solid::cylinder(PlaneFrame::default(), 3., 5., policy).unwrap();
    let bore = Solid::cylinder(PlaneFrame::default(), 1., 5., policy).unwrap();
    let holed = cylinder
        .boolean(&bore, BooleanOperation::Subtract, policy)
        .unwrap()
        .unwrap();
    let vertices = [[0., 0.], [3., 0.], [0., 5.]];
    let profile = CurvedProfile {
        frame: PlaneFrame {
            origin: [0.; 3],
            u: [1., 0., 0.],
            v: [0., 0., 1.],
        },
        outer: (0..3)
            .map(|i| CurvedEdge2d::Line {
                a: vertices[i],
                b: vertices[(i + 1) % 3],
            })
            .collect(),
        holes: vec![],
    };
    let origins: Vec<_> = (0..3)
        .map(|i| multicad_core::references::Origin {
            source: EntityId::new(),
            role: format!("profile/{i}"),
        })
        .collect();
    let cone = Solid::rust_revolve_curved_profile_named(
        &profile,
        RevolutionAxis {
            origin_mm: [0.; 3],
            direction: [0., 0., 1.],
        },
        2. * PI,
        EntityId::new(),
        &origins,
        policy,
    )
    .unwrap();
    for (solid, expected) in [
        (&box_solid, 60.),
        (&cone, PI * 3_f64.powi(2) * 5. / 3.),
        (&holed, PI * (9. - 1.) * 5.),
    ] {
        volume_and_raw_step(solid, expected);
        posed_and_scaled_volume(solid, expected);
    }
}

#[test]
fn centered_basis_flux_preserves_offcenter_void_shells_and_compound_multiplicity() {
    use multicad_kernel::{BooleanOperation, RigidTransform};
    let policy = GeometryPolicy::default();
    let outer = Solid::box_solid([20., 16., 12.], policy).unwrap();
    let inner = Solid::box_solid([3., 4., 5.], policy)
        .unwrap()
        .transformed(RigidTransform {
            translation_mm: [2., 3., 4.],
            ..Default::default()
        })
        .unwrap();
    let cavity = outer
        .boolean(&inner, BooleanOperation::Subtract, policy)
        .unwrap()
        .unwrap();
    assert_eq!(cavity.solid_count().unwrap(), 1);
    assert_eq!(cavity.faces().unwrap().len(), 12);
    volume_and_raw_step(&cavity, 20. * 16. * 12. - 3. * 4. * 5.);
    posed_and_scaled_volume(&cavity, 20. * 16. * 12. - 3. * 4. * 5.);
    let directory = std::env::temp_dir().join(format!("multicad-void-shells-{}", EntityId::new()));
    std::fs::create_dir(&directory).unwrap();
    let path = directory.join("offcenter-void.step");
    cavity.write_step(&path).unwrap();
    assert_eq!(
        // Count CLOSED_SHELL entity records regardless of writer spacing;
        // the leading `=` excludes ORIENTED_CLOSED_SHELL references.
        std::fs::read_to_string(&path)
            .unwrap()
            .replace(" = ", "=")
            .matches("=CLOSED_SHELL(")
            .count(),
        2,
        "raw STEP must retain outer and inward void shells"
    );
    std::fs::remove_dir_all(directory).unwrap();
    let small = Solid::box_solid([3., 4., 5.], policy).unwrap();
    let distant = small
        .transformed(RigidTransform {
            translation_mm: [70., -20., 30.],
            ..Default::default()
        })
        .unwrap();
    // Repeated source geometry remains two separate physical members; a third
    // displaced member exercises shell-local references rather than one center.
    let compound =
        Solid::compound_named(&[&small, &distant, &small], EntityId::new(), policy).unwrap();
    assert_eq!(compound.solid_count().unwrap(), 3);
    volume_and_raw_step(&compound, 3. * 3. * 4. * 5.);
    posed_and_scaled_volume(&compound, 3. * 3. * 4. * 5.);
}
