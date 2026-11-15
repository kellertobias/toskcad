//! Non-axial, cropped, sectioned and overlapping drawings of spherical bevel
//! gears through the exact-curve scene engine (CAD-346), against closed-form
//! geometry and, when linked, the independent native OpenCascade 7.9.3 bevel.
use multicad_core::{
    document::EntityId,
    drawings::{VectorCurve, VectorRole},
};
use multicad_kernel::{
    DrawingHatch, DrawingProjection, DrawingSolid, GeometryPolicy, PlaneFrame, ProjectedVector,
    RigidTransform, Solid, SphericalBevelGear, project_drawing, rust_kernel::Rational,
};
use std::f64::consts::PI;

fn specification(module: f64, teeth: u32, mate: u32, width: f64, bore: f64) -> SphericalBevelGear {
    let (z, mate) = (f64::from(teeth), f64::from(mate));
    let alpha = 20_f64.to_radians();
    let delta = 1_f64.atan2(mate / z);
    let distance = module * z / (2. * delta.sin());
    let beta = (delta.sin() * alpha.cos()).asin();
    let inv = |g: f64| {
        let sigma = (g.cos() / beta.cos()).acos();
        sigma / beta.sin() - sigma.sin().atan2(beta.sin() * sigma.cos())
    };
    SphericalBevelGear {
        teeth,
        pitch_cone_angle_rad: delta,
        base_cone_angle_rad: beta,
        tip_cone_angle_rad: delta + (module / distance).atan(),
        root_cone_angle_rad: delta - (1.25 * module / distance).atan(),
        outer_cone_distance_mm: distance,
        face_width_mm: width,
        bore_radius_mm: bore,
        half_base_tooth_angle_rad: (PI * module / 2. - 0.05) / (module * z) + inv(delta),
        root_fillet_radius_mm: module * 0.2,
    }
}

fn rust(spec: SphericalBevelGear) -> Solid {
    Solid::rust_spherical_bevel_gear_named(spec, EntityId::new(), GeometryPolicy::default())
        .unwrap()
}

fn project(solids: &[&Solid], options: DrawingProjection) -> Vec<ProjectedVector> {
    let inputs: Vec<_> = solids
        .iter()
        .map(|&solid| DrawingSolid { solid, hatch: None })
        .collect();
    project_drawing(&inputs, options).unwrap_or_else(|error| panic!("{error:?}"))
}

fn length(vector: &ProjectedVector) -> f64 {
    vector
        .preview_mm
        .windows(2)
        .map(|w| (w[1][0] - w[0][0]).hypot(w[1][1] - w[0][1]))
        .sum()
}

fn total(vectors: &[ProjectedVector], role: VectorRole, source: usize) -> f64 {
    vectors
        .iter()
        .filter(|vector| vector.role == role && vector.source_index == source)
        .map(length)
        .sum()
}

/// View across the axis: paper x is the gear's x, paper y its axis, and the
/// viewer looks along +y of the gear.
fn side() -> PlaneFrame {
    PlaneFrame {
        origin: [0.; 3],
        u: [1., 0., 0.],
        v: [0., 0., 1.],
    }
}

/// Meridian points of a gear: the back-cone point at polar angle γ.
fn meridian(spec: &SphericalBevelGear, toe: bool, gamma: f64) -> [f64; 2] {
    let distance = spec.outer_cone_distance_mm - if toe { spec.face_width_mm } else { 0. };
    let r = distance / (gamma - spec.pitch_cone_angle_rad).cos();
    [r * gamma.sin(), r * gamma.cos()]
}

fn bore_point(spec: &SphericalBevelGear, toe: bool) -> [f64; 2] {
    let distance = spec.outer_cone_distance_mm - if toe { spec.face_width_mm } else { 0. };
    let (sin, cos) = spec.pitch_cone_angle_rad.sin_cos();
    [
        spec.bore_radius_mm,
        (distance - spec.bore_radius_mm * sin) / cos,
    ]
}

#[test]
fn a_cut_through_the_axis_and_two_tooth_crests_is_two_exact_quadrilaterals() {
    // Twenty teeth: the plane y = 0 passes through the crests of teeth 0
    // and 10. Each half of the cut is bounded by the bore, both back cones
    // and the tip ruling.
    let spec = specification(2., 20, 30, 8., 3.);
    let gear = rust(spec);
    let options = DrawingProjection {
        frame: side(),
        section_offset_mm: Some(0.),
        section_only: true,
        ..Default::default()
    };
    let vectors = project_drawing(
        &[DrawingSolid {
            solid: &gear,
            hatch: Some(DrawingHatch {
                pattern: multicad_core::materials::HatchPattern::Lines,
                angle_rad: 0.7,
                spacing_mm: 1.5,
                line_width_mm: 0.1,
            }),
        }],
        options,
    )
    .unwrap();
    let corners = [
        bore_point(&spec, true),
        meridian(&spec, true, spec.tip_cone_angle_rad),
        meridian(&spec, false, spec.tip_cone_angle_rad),
        bore_point(&spec, false),
    ];
    let mut expected = Vec::new();
    for sign in [-1., 1.] {
        for index in 0..4 {
            let [a, b] = [corners[index], corners[(index + 1) % 4]];
            expected.push([[sign * a[0], a[1]], [sign * b[0], b[1]]]);
        }
    }
    let sections: Vec<_> = vectors
        .iter()
        .filter(|vector| vector.role == VectorRole::Section)
        .collect();
    assert_eq!(sections.len(), 8);
    let mut worst = 0_f64;
    for vector in &sections {
        let VectorCurve::Line { a, b } = vector.curve else {
            panic!("a cut along rulings is a straight line: {:?}", vector.curve);
        };
        let gap = expected
            .iter()
            .map(|[p, q]| {
                let direct = (a[0] - p[0]).hypot(a[1] - p[1]) + (b[0] - q[0]).hypot(b[1] - q[1]);
                let turned = (a[0] - q[0]).hypot(a[1] - q[1]) + (b[0] - p[0]).hypot(b[1] - p[1]);
                direct.min(turned)
            })
            .fold(f64::INFINITY, f64::min);
        worst = worst.max(gap);
    }
    println!("{{\"axis_cut_corner_deviation_mm\":{worst}}}");
    assert!(worst <= 1e-10, "section corners deviate by {worst} mm");
    assert!(
        vectors
            .iter()
            .any(|vector| vector.role == VectorRole::Hatch)
    );
    // The same cut a hair off the axis keeps its shape.
    let near = project(
        &[&gear],
        DrawingProjection {
            section_offset_mm: Some(1e-12),
            ..options
        },
    );
    let perimeter: f64 = expected
        .iter()
        .map(|[p, q]| (p[0] - q[0]).hypot(p[1] - q[1]))
        .sum();
    let drawn = total(&near, VectorRole::Section, 0);
    assert!(
        (drawn - perimeter).abs() <= 1e-9,
        "{drawn} against {perimeter}"
    );
}

#[test]
fn the_side_view_silhouette_reaches_the_exact_extents() {
    for spec in [
        specification(2., 20, 30, 8., 3.),
        specification(2., 16, 16, 6., 0.),
    ] {
        let gear = rust(spec);
        let vectors = project(
            &[&gear],
            DrawingProjection {
                frame: side(),
                hidden_lines: true,
                ..Default::default()
            },
        );
        assert!(total(&vectors, VectorRole::Hidden, 0) > 0.);
        let mut lo = [f64::INFINITY; 2];
        let mut hi = [f64::NEG_INFINITY; 2];
        for vector in vectors.iter().filter(|v| v.role == VectorRole::Visible) {
            for point in &vector.preview_mm {
                for k in 0..2 {
                    lo[k] = lo[k].min(point[k]);
                    hi[k] = hi[k].max(point[k]);
                }
            }
        }
        // The widest circle is the heel tip circle; the paper x of the view
        // is the gear's −x seen from +y, symmetric either way.
        let [radius, _] = meridian(&spec, false, spec.tip_cone_angle_rad);
        let top = if spec.bore_radius_mm > 0. {
            bore_point(&spec, false)[1]
        } else {
            spec.outer_cone_distance_mm / spec.pitch_cone_angle_rad.cos()
        };
        let bottom = meridian(&spec, true, spec.tip_cone_angle_rad)[1];
        // A tooth crest reaches the widest point only when a tooth lies in
        // the paper plane, which holds for these even tooth counts.
        println!(
            "{{\"extent_x\":[{},{}],\"radius\":{radius},\"extent_y\":[{},{}],\"expected_y\":[{bottom},{top}]}}",
            lo[0], hi[0], lo[1], hi[1]
        );
        assert!((hi[0] - radius).abs() <= 1e-9 && (lo[0] + radius).abs() <= 1e-9);
        assert!((hi[1] - top).abs() <= 1e-9 && (lo[1] - bottom).abs() <= 1e-9);
    }
}

#[test]
fn an_oblique_cut_publishes_the_bore_as_exact_ellipse_arcs() {
    let spec = specification(2., 20, 30, 8., 3.);
    let gear = rust(spec);
    let s = 0.5_f64.sqrt();
    let t = (1_f64 / 6.).sqrt();
    // The cut plane has the normal (1, 1, 1) / √3 and passes through the
    // bore above the toe.
    let vectors = project(
        &[&gear],
        DrawingProjection {
            frame: PlaneFrame {
                origin: [0., 0., 36.],
                u: [s, -s, 0.],
                v: [t, t, -2. * t],
            },
            section_offset_mm: Some(0.),
            section_only: true,
            ..Default::default()
        },
    );
    let mut bore = 0;
    let mut conics = 0;
    for vector in vectors.iter().filter(|v| v.role == VectorRole::Section) {
        if let VectorCurve::Ellipse {
            center,
            major_mm,
            minor_mm,
            ..
        } = vector.curve
        {
            conics += 1;
            // A cylinder of radius r cut at 1/√3 to its axis.
            if (minor_mm - 3.).abs() <= 1e-9 && (major_mm - 3. * 3_f64.sqrt()).abs() <= 1e-9 {
                assert!(center[0].hypot(center[1]) <= 1e-9, "{center:?}");
                bore += 1;
            }
        }
    }
    assert!(bore >= 2, "{bore} bore arcs among {conics} conics");
    // The back cones are cut in conics as well.
    assert!(conics > bore);
}

fn rational(value: i64) -> Rational {
    Rational::new(value, 1).unwrap()
}

#[test]
fn stock_in_front_hides_the_gear_and_stock_behind_is_hidden_by_it() {
    let spec = specification(2., 20, 30, 8., 3.);
    let gear = rust(spec);
    let block = |y: f64, width: i64| {
        Solid::rust_box(
            EntityId::new().uuid_bits(),
            [width, 4, 80].map(rational),
            PlaneFrame::default(),
            GeometryPolicy::default(),
        )
        .unwrap()
        .transformed(RigidTransform {
            rotation: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
            translation_mm: [40. - width as f64, y, -10.],
        })
        .unwrap()
    };
    // The viewer of the side frame sits at −y.
    let options = DrawingProjection {
        frame: side(),
        hidden_lines: true,
        ..Default::default()
    };
    let alone = project(&[&gear], options);
    let lines = total(&alone, VectorRole::Visible, 0) + total(&alone, VectorRole::Hidden, 0);
    let front = project(&[&block(-60., 80), &gear], options);
    assert_eq!(total(&front, VectorRole::Visible, 1), 0.);
    let hidden = total(&front, VectorRole::Hidden, 1);
    assert!(
        (hidden - lines).abs() <= 1e-6 * lines,
        "{hidden} against {lines}"
    );
    let behind = project(&[&block(60., 40), &gear], options);
    for role in [VectorRole::Visible, VectorRole::Hidden] {
        let [with, without] = [total(&behind, role, 1), total(&alone, role, 0)];
        assert!(
            (with - without).abs() <= 1e-6 * lines,
            "{role:?}: {with} against {without}"
        );
    }
    // The block's edge behind the gear axis is hidden where the gear covers it.
    assert!(total(&behind, VectorRole::Hidden, 0) > 0.);
    assert!(total(&behind, VectorRole::Visible, 0) > 0.);
}

/// The same drawings from the independently constructed native OpenCascade
/// bevel gear.
#[cfg(feature = "occt")]
mod native {
    use super::*;

    const PLACED: RigidTransform = RigidTransform {
        rotation: [[0.8, 0., 0.6], [0., 1., 0.], [-0.6, 0., 0.8]],
        translation_mm: [10., -20., 30.],
    };

    fn pair(spec: SphericalBevelGear, transform: RigidTransform) -> [Solid; 2] {
        let policy = GeometryPolicy::default();
        let feature = EntityId::new();
        [
            Solid::rust_spherical_bevel_gear_named(spec, feature, policy),
            Solid::spherical_bevel_gear_named(spec, feature, policy),
        ]
        .map(|solid| solid.unwrap().transformed(transform).unwrap())
    }

    /// Views along the placed gear axis, across it and oblique to it.
    fn frames(transform: RigidTransform, origin: [f64; 3]) -> [(&'static str, PlaneFrame); 3] {
        let column = |k: usize| -> [f64; 3] { std::array::from_fn(|i| transform.rotation[i][k]) };
        let [x, y, z] = [column(0), column(1), column(2)];
        let mix = |a: f64, b: f64, c: f64| -> [f64; 3] {
            std::array::from_fn(|i| a * x[i] + b * y[i] + c * z[i])
        };
        let s = 0.5_f64.sqrt();
        let t = (1_f64 / 6.).sqrt();
        [
            ("axial", [x, y]),
            ("side", [x, z]),
            ("oblique", [mix(s, -s, 0.), mix(t, t, -2. * t)]),
        ]
        .map(|(name, [u, v])| (name, PlaneFrame { origin, u, v }))
    }

    fn cases(transform: RigidTransform, origin: [f64; 3]) -> Vec<(String, DrawingProjection)> {
        let mut cases = Vec::new();
        for (name, frame) in frames(transform, origin) {
            let view = DrawingProjection {
                frame,
                ..Default::default()
            };
            let section = DrawingProjection {
                section_offset_mm: Some(0.),
                ..view
            };
            cases.extend([
                (name.to_string(), view),
                (
                    format!("{name} hidden"),
                    DrawingProjection {
                        hidden_lines: true,
                        ..view
                    },
                ),
                (
                    format!("{name} detail"),
                    DrawingProjection {
                        detail: Some(([9., 2.], 6.)),
                        hidden_lines: true,
                        ..view
                    },
                ),
                (format!("{name} section"), section),
                (
                    format!("{name} section hidden"),
                    DrawingProjection {
                        hidden_lines: true,
                        section_offset_mm: Some(0.3),
                        ..section
                    },
                ),
                (
                    format!("{name} section-only"),
                    DrawingProjection {
                        section_only: true,
                        ..section
                    },
                ),
                (
                    format!("{name} finite depth"),
                    DrawingProjection {
                        section_depth_mm: Some(2.),
                        ..section
                    },
                ),
                (
                    format!("{name} window"),
                    DrawingProjection {
                        section_u_mm: Some([-14., 12.5]),
                        section_v_mm: Some([-3.5, 9.]),
                        section_depth_mm: Some(4.5),
                        ..section
                    },
                ),
                (
                    format!("{name} section detail"),
                    DrawingProjection {
                        detail: Some(([9., 1.5], 7.)),
                        ..section
                    },
                ),
            ]);
        }
        cases
    }

    /// The largest distance from a preview point or chord midpoint of `from`
    /// to the previews of `to` with the same source and one of the roles.
    fn deviation(
        from: &[ProjectedVector],
        to: &[ProjectedVector],
        roles: &dyn Fn(VectorRole) -> Vec<VectorRole>,
    ) -> (f64, [f64; 2], VectorRole, usize) {
        let mut worst = (0., [0.; 2], VectorRole::Visible, 0);
        for vector in from.iter().filter(|v| v.role != VectorRole::Hatch) {
            let accepted = roles(vector.role);
            let spans: Vec<_> = to
                .iter()
                .filter(|other| {
                    other.source_index == vector.source_index && accepted.contains(&other.role)
                })
                .flat_map(|other| other.preview_mm.windows(2).map(|w| [w[0], w[1]]))
                .collect();
            for pair in vector.preview_mm.windows(2) {
                let mid = [
                    pair[0][0].midpoint(pair[1][0]),
                    pair[0][1].midpoint(pair[1][1]),
                ];
                for [x, y] in [pair[0], mid, pair[1]] {
                    let distance = spans
                        .iter()
                        .map(|&[[a, b], [c, d]]| {
                            let [u, v] = [c - a, d - b];
                            let along = (((x - a) * u + (y - b) * v) / (u * u + v * v).max(1e-300))
                                .clamp(0., 1.);
                            (x - a - along * u).hypot(y - b - along * v)
                        })
                        .fold(f64::INFINITY, f64::min);
                    if distance > worst.0 {
                        worst = (distance, [x, y], vector.role, vector.source_index);
                    }
                }
            }
        }
        worst
    }

    fn center(solid: &Solid) -> [f64; 3] {
        let b = solid.bounds().unwrap();
        std::array::from_fn(|i| (b.min[i] + b.max[i]) / 2.)
    }

    /// Rust curves lie on native curves and native curves on Rust curves
    /// within three chord tolerances (measured worst case 0.0122 mm).
    ///
    /// Measured OpenCascade deviations, each bounded below: it refuses
    /// oblique cuts and hidden-line lengthwise cuts of its own bevel; it
    /// bridges transverse flank cuts with chords (sagitta up to 0.174 mm);
    /// it draws the tangent-continuous flank seams; and in hidden-line views
    /// it publishes up to 15.4 mm of rear tooth edges as visible although
    /// material lies in front of them, for example 2.49 mm at
    /// (12.185, 13.898) of the bored oblique view. The member's unit test
    /// `material_in_front_of_an_edge_is_found_by_marching_the_membership`
    /// decides those cases from the exact membership predicate. No Rust
    /// visible curve lies on a native curve of another role.
    fn check(
        label: &str,
        rust: &[&Solid],
        native: &[&Solid],
        transform: RigidTransform,
        origin: [f64; 3],
    ) {
        use VectorRole::{Hidden, Section, Visible};
        let only = std::env::var("SCENE_CASE").ok();
        let mut failures = Vec::new();
        for (name, options) in cases(transform, origin) {
            let label = format!("{label} {name}");
            if only
                .as_ref()
                .is_some_and(|only| !label.contains(only.as_str()))
            {
                continue;
            }
            let ours = project(rust, options);
            let inputs: Vec<_> = native
                .iter()
                .map(|&solid| DrawingSolid { solid, hatch: None })
                .collect();
            let reference = match project_drawing(&inputs, options) {
                Ok(reference) => reference,
                Err(error) => {
                    println!(
                        "{{\"case\":\"{label}\",\"native_error\":\"{}\"}}",
                        error.message
                    );
                    continue;
                }
            };
            let bound = 3. * options.chord_tolerance_mm;
            let any = |_: VectorRole| vec![Visible, Hidden, Section];
            let complete = project(
                rust,
                DrawingProjection {
                    tangent_edges: true,
                    hidden_lines: true,
                    ..options
                },
            );
            // OpenCascade replaces stretches of a cut through a flank by one
            // straight chord between two true boundary points. Such a native
            // chord has both ends on the Rust curves and its middle off them.
            let (bridges, reference): (Vec<_>, Vec<_>) =
                reference.into_iter().partition(|vector| {
                    vector.preview_mm.len() == 2
                        && options.section_offset_mm.is_some()
                        && deviation(std::slice::from_ref(vector), &complete, &any).0 > bound
                        && vector.preview_mm.iter().all(|end| {
                            let probe = ProjectedVector {
                                preview_mm: vec![*end, *end],
                                ..vector.clone()
                            };
                            deviation(&[probe], &complete, &any).0 <= bound
                        })
                });
            let sagitta = deviation(&bridges, &complete, &any).0;
            // Geometry: every Rust curve lies on a native curve or under a
            // native bridge, and every other native curve on a Rust curve
            // (the native view draws the tangent-continuous flank seams).
            let covered: Vec<_> = reference.iter().chain(&bridges).cloned().collect();
            let exact: Vec<_> = ours
                .iter()
                .filter(|vector| {
                    deviation(std::slice::from_ref(vector), &reference, &any).0 <= bound
                        || deviation(std::slice::from_ref(vector), &bridges, &any).0 > 0.25
                })
                .cloned()
                .collect();
            let forward = deviation(&exact, &covered, &any);
            let backward = deviation(&reference, &complete, &any);
            // Roles: the length of Rust chords whose native curve carries
            // another role. OpenCascade repeats rear edges under coincident
            // visible ones, so a hidden native twin is no disagreement.
            let mut disagreement = [0.; 2];
            for vector in ours.iter().filter(|v| v.role != VectorRole::Hatch) {
                for pair in vector.preview_mm.windows(2) {
                    let span = (pair[1][0] - pair[0][0]).hypot(pair[1][1] - pair[0][1]);
                    let steps = (span / 0.05).ceil().max(1.) as usize;
                    for step in 0..steps {
                        let at = (step as f64 + 0.5) / steps as f64;
                        let point = [0, 1].map(|k| pair[0][k] + (pair[1][k] - pair[0][k]) * at);
                        let probe = ProjectedVector {
                            preview_mm: vec![point, point],
                            ..vector.clone()
                        };
                        let probe = std::slice::from_ref(&probe);
                        if deviation(probe, &reference, &|role| vec![role]).0 > bound
                            && deviation(probe, &bridges, &any).0 > 0.25
                        {
                            disagreement[usize::from(vector.role == Hidden)] += span / steps as f64;
                        }
                    }
                }
            }
            let [shown, concealed] = disagreement;
            println!(
                "{{\"case\":\"{label}\",\"rust_to_native_mm\":{},\"native_to_rust_mm\":{},\"rust_visible_native_other_mm\":{shown},\"rust_hidden_native_other_mm\":{concealed},\"native_bridges\":{},\"bridge_sagitta_mm\":{sagitta},\"rust\":{},\"native\":{}}}",
                forward.0,
                backward.0,
                bridges.len(),
                ours.len(),
                reference.len()
            );
            if std::env::var("SCENE_DUMP").is_ok() {
                for (name, set) in [("rust", &ours), ("native", &covered)] {
                    for v in set.iter() {
                        println!(
                            "DUMP|{name}|{:?}|{}|{:?}",
                            v.role, v.source_index, v.preview_mm
                        );
                    }
                }
            }
            if forward.0 > bound {
                failures.push(format!("{label}: Rust {forward:?}"));
            }
            if backward.0 > bound {
                failures.push(format!("{label}: OpenCascade {backward:?}"));
            }
            if shown > 0.5 || concealed > 16. || sagitta > 0.25 {
                failures.push(format!(
                    "{label}: roles {shown} mm visible, {concealed} mm hidden, bridge {sagitta} mm"
                ));
            }
        }
        assert!(failures.is_empty(), "{failures:#?}");
    }

    #[test]
    fn bevel_views_and_sections_match_the_native_gear() {
        for (label, spec, transform) in [
            (
                "bored",
                specification(2., 20, 30, 8., 3.),
                RigidTransform::default(),
            ),
            ("solid placed", specification(2., 16, 16, 6., 0.), PLACED),
        ] {
            let [rust, native] = pair(spec, transform);
            check(label, &[&rust], &[&native], transform, center(&rust));
        }
    }

    #[test]
    fn stock_and_bevel_scenes_match_the_native_scene() {
        let spec = specification(2., 20, 30, 8., 3.);
        let [rust, native] = pair(
            spec,
            RigidTransform {
                rotation: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
                translation_mm: [2., 0., 20.],
            },
        );
        let stock = Solid::rust_box(
            EntityId::new().uuid_bits(),
            [20, 12, 8].map(rational),
            PlaneFrame::default(),
            GeometryPolicy::default(),
        )
        .unwrap();
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("stock.step");
        stock.write_step(&path).unwrap();
        let imported = Solid::read_step(&path, stock.geometry_policy()).unwrap();
        check(
            "box+bevel",
            &[&stock, &rust],
            &[&imported, &native],
            RigidTransform::default(),
            [2., 0., 45.],
        );
    }
}
