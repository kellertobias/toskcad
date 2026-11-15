//! Independent OpenCascade 7.9.3 references for the exact-curve scene engine
//! (CAD-346): every newly drawn class of view, section and overlapping scene
//! is compared by role with native hidden-line removal of the same solids.
#![cfg(feature = "occt")]
use multicad_core::drawings::VectorRole;
use multicad_core::{document::EntityId, references::Origin};
use multicad_kernel::{
    DrawingProjection, DrawingSolid, GearHub, GearSweep, GeometryPolicy, InvoluteGear, PlaneFrame,
    ProjectedVector, RigidTransform, Solid, project_drawing,
    rust_kernel::{Circle, Rational, cylinder::AxialCylinder},
};
use std::f64::consts::PI;

fn spur(hub: Option<GearHub>) -> Solid {
    let a = 20_f64.to_radians();
    let pitch = 12.;
    Solid::rust_spur_gear_named(
        InvoluteGear {
            hub,
            teeth: 24,
            pitch_radius_mm: pitch,
            base_radius_mm: pitch * a.cos(),
            tip_radius_mm: pitch + 1.,
            root_radius_mm: pitch - 1.25,
            root_fillet_radius_mm: 0.2,
            half_base_angle_rad: (PI / 2. - 0.05) / (2. * pitch) + a.tan() - a,
            face_width_mm: 4.,
            bore_radius_mm: 3.,
            twist_rad_per_mm: 0.,
            sweep: GearSweep::Spur,
        },
        EntityId::new(),
        GeometryPolicy::default(),
    )
    .unwrap()
}

fn ring() -> Solid {
    let a = 20_f64.to_radians();
    let pitch = 24.;
    Solid::rust_internal_ring_named(
        InvoluteGear {
            hub: None,
            teeth: 48,
            pitch_radius_mm: pitch,
            base_radius_mm: pitch * a.cos(),
            tip_radius_mm: pitch - 1.,
            root_radius_mm: pitch + 1.25,
            root_fillet_radius_mm: 0.2,
            half_base_angle_rad: (PI / 2. - 0.05) / (2. * pitch) - a.tan() + a,
            face_width_mm: 2.,
            bore_radius_mm: 0.,
            twist_rad_per_mm: 0.,
            sweep: GearSweep::InternalSpur {
                outer_radius_mm: pitch + 4.,
            },
        },
        EntityId::new(),
        GeometryPolicy::default(),
    )
    .unwrap()
}

/// Turn a gear about its axis so that no face seam of the exchanged solid
/// falls on a silhouette of the axis-aligned views.
fn turned(solid: Solid, translation_mm: [f64; 3]) -> Solid {
    let (sin, cos) = 0.1_f64.sin_cos();
    solid
        .transformed(RigidTransform {
            rotation: [[cos, -sin, 0.], [sin, cos, 0.], [0., 0., 1.]],
            translation_mm,
        })
        .unwrap()
}

fn rational(value: i64) -> Rational {
    Rational::new(value, 1).unwrap()
}

fn stock() -> Solid {
    Solid::rust_box(
        EntityId::new().uuid_bits(),
        [20, 12, 8].map(rational),
        PlaneFrame::default(),
        GeometryPolicy::default(),
    )
    .unwrap()
}

fn cylinder() -> Solid {
    Solid::from_axial_cylinder(
        AxialCylinder::new(
            EntityId::new().uuid_bits(),
            Circle::new([rational(0); 2], rational(5)).unwrap(),
            rational(0),
            rational(8),
        )
        .unwrap(),
        PlaneFrame::default(),
        GeometryPolicy::default(),
    )
    .unwrap()
}

fn annulus() -> Solid {
    Solid::rust_annulus(
        3.,
        6.,
        2.5,
        EntityId::new(),
        PlaneFrame::default(),
        GeometryPolicy::default(),
    )
    .unwrap()
}

fn slot() -> Solid {
    let owner = EntityId::new();
    Solid::rust_extrude_rounded_slot(
        PlaneFrame::default(),
        [[0., 0.], [6., 0.]],
        2.,
        5.,
        owner,
        std::array::from_fn(|edge| Origin {
            source: owner,
            role: format!("profile.{edge}"),
        }),
        GeometryPolicy::default(),
    )
    .unwrap()
}

fn center(solid: &Solid) -> [f64; 3] {
    let b = solid.bounds().unwrap();
    std::array::from_fn(|i| (b.min[i] + b.max[i]) / 2.)
}

fn frames(origin: [f64; 3]) -> [(&'static str, PlaneFrame); 3] {
    let s = 0.5_f64.sqrt();
    let t = (1_f64 / 6.).sqrt();
    [
        ("axial", [[1., 0., 0.], [0., 1., 0.]]),
        ("side", [[1., 0., 0.], [0., 0., 1.]]),
        ("oblique", [[s, -s, 0.], [t, t, -2. * t]]),
    ]
    .map(|(name, [u, v])| (name, PlaneFrame { origin, u, v }))
}

/// Plain and hidden-line views, then sections through `origin`: whole,
/// section-only, finite depth, lateral window and detail.
fn cases(origin: [f64; 3]) -> Vec<(String, DrawingProjection)> {
    let mut cases = Vec::new();
    for (name, frame) in frames(origin) {
        let view = DrawingProjection {
            frame,
            ..Default::default()
        };
        let section = DrawingProjection {
            section_offset_mm: Some(0.),
            ..view
        };
        cases.extend([
            (format!("{name}"), view),
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
                    detail: Some(([0.5, 0.25], 3.)),
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
                    section_depth_mm: Some(1.),
                    ..section
                },
            ),
            (
                format!("{name} window"),
                DrawingProjection {
                    section_u_mm: Some([-4., 4.5]),
                    section_v_mm: Some([-1.5, 3.]),
                    section_depth_mm: Some(2.5),
                    ..section
                },
            ),
            (
                format!("{name} section detail"),
                DrawingProjection {
                    detail: Some(([1., 0.5], 4.)),
                    ..section
                },
            ),
        ]);
    }
    cases
}

fn project(solids: &[&Solid], options: DrawingProjection) -> Vec<ProjectedVector> {
    let inputs: Vec<_> = solids
        .iter()
        .map(|&solid| DrawingSolid { solid, hatch: None })
        .collect();
    project_drawing(&inputs, options).unwrap_or_else(|error| panic!("{error:?}"))
}

/// The largest distance from a preview point or chord midpoint of `from` to
/// the previews of `to` with the same source and, if required, role.
fn deviation(
    from: &[ProjectedVector],
    to: &[ProjectedVector],
    same_role: bool,
) -> (f64, [f64; 2], VectorRole, usize) {
    let mut worst = (0., [0.; 2], VectorRole::Visible, 0);
    for vector in from {
        let spans: Vec<_> = to
            .iter()
            .filter(|other| {
                other.source_index == vector.source_index
                    && (!same_role || other.role == vector.role)
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

/// Rust curves lie on native curves of the same role and source, and native
/// curves on Rust curves of the same source, within three chord tolerances.
/// Native rear edges that coincide with a visible edge are published twice,
/// so the second direction accepts any role. With `seams`, the native curves
/// are compared with the Rust drawing that includes tangent-continuous seams:
/// the STEP import of an involute flank does not mark its seam with the
/// radial root face as smooth, so OpenCascade draws it as a sharp edge.
fn check_with(label: &str, solids: &[&Solid], origin: [f64; 3], seams: bool, skip: &[&str]) {
    let only = std::env::var("SCENE_CASE").ok();
    let directory = tempfile::tempdir().unwrap();
    let native: Vec<Solid> = solids
        .iter()
        .enumerate()
        .map(|(index, solid)| {
            let path = directory.path().join(format!("solid-{index}.step"));
            solid.write_step(&path).unwrap();
            Solid::read_step(&path, solid.geometry_policy()).unwrap()
        })
        .collect();
    let native: Vec<&Solid> = native.iter().collect();
    let mut failures = Vec::new();
    for (name, options) in cases(origin) {
        let label = format!("{label} {name}");
        if only
            .as_ref()
            .is_some_and(|only| !label.contains(only.as_str()))
            || skip.iter().any(|skipped| name == *skipped)
        {
            continue;
        }
        let rust = project(solids, options);
        let mut reference = project(&native, options);
        let complete = project(
            solids,
            DrawingProjection {
                tangent_edges: seams,
                hidden_lines: options.hidden_lines || seams,
                ..options
            },
        );
        // OpenCascade's section of the imported flank B-spline surfaces by a
        // plane across the axis bridges short stretches of the boundary with
        // chords up to 0.104 mm off the profile (measured here and in
        // `rust_gear_section_drawing.rs`). Its axial view of the same import
        // is the section profile of a straight band, so that view is the
        // independent reference and the native gap chords are not required
        // to lie on the Rust rims.
        let bridged = seams && name.starts_with("axial") && options.section_offset_mm.is_some();
        if bridged {
            let profile = project(
                &native,
                DrawingProjection {
                    section_offset_mm: None,
                    section_only: false,
                    section_depth_mm: None,
                    section_u_mm: None,
                    section_v_mm: None,
                    detail: None,
                    hidden_lines: true,
                    ..options
                },
            );
            reference.extend(profile.into_iter().flat_map(|vector| {
                [VectorRole::Visible, VectorRole::Hidden, VectorRole::Section].map(|role| {
                    ProjectedVector {
                        role,
                        ..vector.clone()
                    }
                })
            }));
        }
        if seams {
            // Measured on the imported gears: OpenCascade publishes rear root
            // fillet outlines and flank seams as visible although the sight
            // line through them first crosses the solid root cylinder (for
            // example at u = -8.9057 mm, 12 mm of material lie in front).
            // A hidden Rust curve may therefore lie on a native curve of
            // either role; `rust_scene_gear.rs` certifies the classification
            // against the exact profile instead.
            let relabeled: Vec<_> = reference
                .iter()
                .filter(|vector| vector.role == VectorRole::Visible)
                .map(|vector| ProjectedVector {
                    role: VectorRole::Hidden,
                    ..vector.clone()
                })
                .collect();
            reference.extend(relabeled);
        }
        // OpenCascade bridges stretches of an oblique or lengthwise cut
        // through an imported flank with one straight chord between two true
        // boundary points (sagitta up to 0.104 mm, as measured for the axial
        // cuts above). A native chord is such a bridge when both of its ends
        // lie on the Rust curves but its middle leaves them. Bridges are not
        // required to lie on the Rust curves, and the Rust curve under a
        // bridge is compared with the bridge's ends instead.
        let bound = 3. * options.chord_tolerance_mm;
        let mut bridges = Vec::new();
        if seams && options.section_offset_mm.is_some() {
            reference.retain(|vector| {
                let ends = [vector.preview_mm[0], vector.preview_mm[1]];
                let bridge = vector.preview_mm.len() == 2
                    && deviation(std::slice::from_ref(vector), &complete, false).0 > bound
                    && ends.iter().all(|end| {
                        let probe = ProjectedVector {
                            preview_mm: vec![*end, *end],
                            ..vector.clone()
                        };
                        deviation(&[probe], &complete, false).0 <= bound
                    });
                if bridge {
                    bridges.push(vector.clone());
                }
                !bridge
            });
        }
        let under_bridge = |vector: &ProjectedVector| {
            bridges.iter().any(|bridge| {
                bridge.source_index == vector.source_index
                    && deviation(
                        std::slice::from_ref(vector),
                        std::slice::from_ref(bridge),
                        false,
                    )
                    .0 <= 0.11
            })
        };
        let exact: Vec<ProjectedVector> = rust
            .iter()
            .filter(|vector| {
                !under_bridge(vector)
                    || deviation(std::slice::from_ref(vector), &reference, true).0 <= bound
            })
            .cloned()
            .collect();
        assert!(
            bridges.len() <= 4 * 24,
            "{label}: {} bridges",
            bridges.len()
        );
        let forward = deviation(&exact, &reference, true);
        let backward = if bridged {
            forward
        } else {
            deviation(&reference, &complete, false)
        };
        println!(
            "{{\"case\":\"{label}\",\"rust_to_native_mm\":{},\"native_to_rust_mm\":{},\"native_bridges\":{}}}",
            forward.0,
            backward.0,
            bridges.len()
        );
        if forward.0 > bound {
            failures.push(format!("{label}: Rust {forward:?}"));
        }
        if backward.0 > bound {
            failures.push(format!("{label}: OpenCascade {backward:?}"));
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}

fn check(label: &str, solids: &[&Solid], origin: [f64; 3]) {
    check_with(label, solids, origin, false, &[]);
}

#[test]
fn prism_families_match_native_open_cascade() {
    for (label, solid) in [
        ("box", stock()),
        ("cylinder", cylinder()),
        ("annulus", annulus()),
        ("slot", slot()),
    ] {
        check(label, &[&solid], center(&solid));
    }
}

#[test]
fn spur_gear_matches_native_open_cascade() {
    let solid = turned(spur(None), [0.; 3]);
    check_with("spur", &[&solid], center(&solid), true, &[]);
}

#[test]
fn hubbed_spur_gear_matches_native_open_cascade() {
    let solid = turned(
        spur(Some(GearHub {
            radius_mm: 6.,
            axial_start_mm: -1.5,
            axial_end_mm: 6.,
        })),
        [0.; 3],
    );
    check_with("hub", &[&solid], [0.4, 0.3, 5.], true, &[]);
}

#[test]
fn internal_ring_matches_native_open_cascade() {
    let solid = turned(ring(), [0.; 3]);
    check_with("ring", &[&solid], center(&solid), true, &[]);
}

#[test]
fn overlapping_scenes_match_native_open_cascade() {
    let place = |solid: &Solid| {
        solid
            .transformed(RigidTransform {
                rotation: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
                translation_mm: [2., 0., 20.],
            })
            .unwrap()
    };
    let stock = stock();
    for (label, first, second) in [
        ("box+spur", &stock, turned(spur(None), [2., 0., 20.])),
        ("box+ring", &stock, turned(ring(), [2., 0., 20.])),
        ("spur+spur", &spur(None), turned(spur(None), [2., 0., 20.])),
        // Offset so that the cylinder rim is not tangent to the annulus bore.
        (
            "cylinder+annulus",
            &cylinder(),
            annulus()
                .transformed(RigidTransform {
                    rotation: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
                    translation_mm: [2.3, 0., 20.],
                })
                .unwrap(),
        ),
    ] {
        check_with(label, &[first, &second], [2., 0., 10.], true, &[]);
    }
}

/// Planes lying exactly in a cap or in the junction of a hub and its tooth
/// band, and a scene whose rim is tangent to an occluding bore.
#[test]
fn coincident_planes_and_tangent_scenes_match_native_open_cascade() {
    let solid = cylinder();
    for level in [0., 8.] {
        check(&format!("cylinder cap {level}"), &[&solid], [0., 0., level]);
    }
    let gear = turned(spur(None), [0.; 3]);
    for level in [0., 4.] {
        check_with(
            &format!("spur cap {level}"),
            &[&gear],
            [0., 0., level],
            true,
            &[],
        );
    }
    let hubbed = turned(
        spur(Some(GearHub {
            radius_mm: 6.,
            axial_start_mm: -1.5,
            axial_end_mm: 6.,
        })),
        [0.; 3],
    );
    check_with("hub junction", &[&hubbed], [0., 0., 4.], true, &[]);
    let ring = annulus()
        .transformed(RigidTransform {
            rotation: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
            translation_mm: [2., 0., 20.],
        })
        .unwrap();
    check("tangent cylinder+annulus", &[&solid, &ring], [2., 0., 10.]);
}
