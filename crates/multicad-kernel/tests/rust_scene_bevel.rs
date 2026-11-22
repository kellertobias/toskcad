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
use multicad_precision::PortableMath as _;
use std::f64::consts::PI;

fn specification(module: f64, teeth: u32, mate: u32, width: f64, bore: f64) -> SphericalBevelGear {
    let (z, mate) = (f64::from(teeth), f64::from(mate));
    let alpha = 20_f64.to_radians();
    let delta = 1_f64.portable_atan2(mate / z);
    let distance = module * z / (2. * delta.portable_sin());
    let beta = (delta.portable_sin() * alpha.portable_cos()).portable_asin();
    let inv = |g: f64| {
        let sigma = (g.portable_cos() / beta.portable_cos()).portable_acos();
        sigma / beta.portable_sin()
            - sigma
                .portable_sin()
                .portable_atan2(beta.portable_sin() * sigma.portable_cos())
    };
    SphericalBevelGear {
        teeth,
        pitch_cone_angle_rad: delta,
        base_cone_angle_rad: beta,
        tip_cone_angle_rad: delta + (module / distance).portable_atan(),
        root_cone_angle_rad: delta - (1.25 * module / distance).portable_atan(),
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
        .map(|w| (w[1][0] - w[0][0]).portable_hypot(w[1][1] - w[0][1]))
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
    let r = distance / (gamma - spec.pitch_cone_angle_rad).portable_cos();
    [r * gamma.portable_sin(), r * gamma.portable_cos()]
}

fn bore_point(spec: &SphericalBevelGear, toe: bool) -> [f64; 2] {
    let distance = spec.outer_cone_distance_mm - if toe { spec.face_width_mm } else { 0. };
    let (sin, cos) = spec.pitch_cone_angle_rad.portable_sin_cos();
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
                let direct = (a[0] - p[0]).portable_hypot(a[1] - p[1])
                    + (b[0] - q[0]).portable_hypot(b[1] - q[1]);
                let turned = (a[0] - q[0]).portable_hypot(a[1] - q[1])
                    + (b[0] - p[0]).portable_hypot(b[1] - p[1]);
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
        .map(|[p, q]| (p[0] - q[0]).portable_hypot(p[1] - q[1]))
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
            spec.outer_cone_distance_mm / spec.pitch_cone_angle_rad.portable_cos()
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
                assert!(center[0].portable_hypot(center[1]) <= 1e-9, "{center:?}");
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
