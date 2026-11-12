//! Actual five-band geometry and adversarial controls through public stages.
use super::*;

fn graph(profile: &[[f64; 2]], root_arc: bool) -> BoundaryInput {
    let bands = if root_arc {
        profile.len() - 1
    } else {
        profile.len()
    };
    let mut ring = profile.to_vec();
    if root_arc {
        ring.push([-5., 0.]);
    }
    let n = ring.len();
    let rails = profile.len();
    let axis = |z| Axis {
        origin: [0., 0., z],
        normal: [0., 0., 1.],
        reference: [1., 0., 0.],
    };
    let vertices = (0..2)
        .flat_map(|end| ring.iter().map(move |p| [p[0], p[1], 4. * end as f64]))
        .collect::<Vec<_>>();
    let mut edges = Vec::new();
    let mut curves = Vec::new();
    let line = |a: Point, b: Point| Curve::Spline {
        spans: 1,
        poles: (0..=9)
            .map(|i| {
                if i == 0 {
                    a
                } else if i == 9 {
                    b
                } else {
                    std::array::from_fn(|k| a[k] + (b[k] - a[k]) * i as f64 / 9.)
                }
            })
            .collect(),
    };
    for end in 0..2 {
        for i in 0..n {
            let [a, b] = [end * n + i, end * n + (i + 1) % n];
            let curve = curves.len();
            curves.push(if i < bands {
                line(vertices[a], vertices[b])
            } else {
                Curve::Circle {
                    axis: axis(4. * end as f64),
                    radius_mm: 5.,
                }
            });
            edges.push(Edge {
                vertices: [a, b],
                curve,
                forward: true,
            });
        }
    }
    for i in 0..rails {
        let curve = curves.len();
        curves.push(line(vertices[i], vertices[n + i]));
        edges.push(Edge {
            vertices: [i, n + i],
            curve,
            forward: true,
        });
    }
    let mut surfaces = vec![
        Surface::Plane(Axis {
            normal: [0., 0., -1.],
            ..axis(0.)
        }),
        Surface::Plane(axis(4.)),
    ];
    let mut faces = vec![
        Face {
            surface: 0,
            outer: (0..n).rev().map(|i| (i, false)).collect(),
            forward: true,
        },
        Face {
            surface: 1,
            outer: (n..2 * n).map(|i| (i, true)).collect(),
            forward: true,
        },
    ];
    for (i, curve) in curves.iter().enumerate().take(bands) {
        let Curve::Spline { poles: base, .. } = curve else {
            panic!("actual spline rim");
        };
        let poles = base
            .iter()
            .flat_map(|p| (0..=9).map(move |j| [p[0], p[1], 4. * j as f64 / 9.]))
            .collect();
        let surface = surfaces.len();
        surfaces.push(Surface::Tensor(RetainedTensor::new([1, 1], poles).unwrap()));
        faces.push(Face {
            surface,
            outer: vec![
                (i, true),
                (2 * n + (i + 1) % rails, true),
                (n + i, false),
                (2 * n + i, false),
            ],
            forward: true,
        });
    }
    if root_arc {
        let surface = surfaces.len();
        surfaces.push(Surface::Cone {
            axis: axis(0.),
            radius_mm: 5.,
            semi_angle_rad: 0.,
        });
        faces.push(Face {
            surface,
            outer: vec![
                (bands, true),
                (bands + 1, true),
                (2 * n, true),
                (n + bands + 1, false),
                (n + bands, false),
                (2 * n + rails - 1, false),
            ],
            forward: true,
        });
    }
    BoundaryInput {
        vertices,
        curves,
        edges,
        surfaces,
        faces,
    }
}
fn strict_source_proofs(input: &BoundaryInput, core: f64) {
    let mut count = 0;
    for surface in &input.surfaces {
        if let Surface::Tensor(t) = surface {
            count += 1;
            t.certify_radial_chart(true, None).unwrap();
            // Every source xy coefficient pair has dot product above core².
            // Positive Bernstein weights sum to one, proving the entire actual
            // polynomial wall is outside this core, rather than testing samples.
            for p in t.poles() {
                for q in t.poles() {
                    let square = I::point(p[0]) * I::point(q[0]) + I::point(p[1]) * I::point(q[1]);
                    assert!(square.lo > (I::point(core).square()).hi);
                }
            }
        }
    }
    assert_eq!(count, 5);
    let n = input.vertices.len() / 2;
    let mut angle = I::point(0.);
    for i in 0..n {
        let xy = |j: usize| {
            [
                I::point(input.vertices[j][0]),
                I::point(input.vertices[j][1]),
            ]
        };
        angle = angle
            + crate::radial_tensor::ordering::angle(xy(i), xy((i + 1) % n))
                .unwrap()
                .unwrap();
    }
    assert!(angle.lo > std::f64::consts::PI && angle.hi < 3. * std::f64::consts::PI);
    assert!(angle.lo <= std::f64::consts::TAU && angle.hi >= std::f64::consts::TAU);
}
fn owned_and_bound(input: BoundaryInput) -> TrimmedBoundary {
    let owned = BoundaryOwnership::new(input, 1e-7, None).unwrap();
    assert!(owned.topology.edge_faces.iter().all(|e| e[0] != e[1]));
    owned.bind_surfaces(1e-7, None).unwrap()
}
fn source_candidate(text: &str) -> BoundaryOwnership {
    let members = crate::step_read::read_source_compound(text, 1e-7).unwrap();
    assert_eq!(members.len(), 1);
    members
        .into_iter()
        .find_map(|m| {
            if let crate::step_read::SourceImport::Radial(body) = m {
                Some(*body)
            } else {
                None
            }
        })
        .unwrap()
}
#[test]
fn no_cone_axis_proposals_precision_and_cancelled_stages_leave_source_immutable() {
    let input = graph(
        &[[5., 0.], [2., 5.], [-4., 3.], [-4., -3.], [2., -5.]],
        false,
    );
    let bound = owned_and_bound(input.clone());
    assert!(bound.clone().certify_axial(1e-30, 1e-6, None).is_err());
    assert!(bound.clone().certify_axial(1e-7, 10., None).is_err());
    let calls = std::cell::Cell::new(0);
    let observe = || {
        calls.set(calls.get() + 1);
        false
    };
    bound
        .clone()
        .certify_axial(1e-7, 1e-6, Some(&observe))
        .unwrap();
    let total = calls.get();
    assert!(total > 10);
    for at in [1, total / 2, total] {
        calls.set(0);
        let cancel = || {
            calls.set(calls.get() + 1);
            calls.get() >= at
        };
        assert_eq!(
            bound.clone().certify_axial(1e-7, 1e-6, Some(&cancel)).err(),
            Some(Error::Cancelled)
        );
        assert_eq!(bound.input(), &input);
    }
    let axial = bound.clone().certify_axial(1e-7, 1e-6, None).unwrap();
    assert_eq!(
        axial.clone().certify_caps(Some(&|| true)).err(),
        Some(Error::Cancelled)
    );
    let caps = axial.certify_caps(None).unwrap();
    assert_eq!(
        caps.clone().certify_ordering(Some(&|| true)).err(),
        Some(Error::Cancelled)
    );
    let ordered = caps.certify_ordering(None).unwrap();
    assert_eq!(
        ordered.clone().certify_material(1e-7, Some(&|| true)).err(),
        Some(Error::Cancelled)
    );
    assert_eq!(
        ordered.certify_material(1e-7, None).unwrap().input(),
        &input
    );

    // A plane placement only proposes an axis. Moving that placement onto an
    // actual boundary vertex does not prove any positive core or new material.
    let mut bad_axis = input.clone();
    if let Surface::Plane(axis) = &mut bad_axis.surfaces[0] {
        axis.origin = input.vertices[0];
    }
    let bad = owned_and_bound(bad_axis.clone());
    assert!(bad.certify_axial(1e-7, 1e-6, None).is_err());
    assert_eq!(bound.input(), &input);
}
#[test]
fn no_exposed_cone_material_requalifies_actual_proper_and_distant_source_frames() {
    let input = graph(
        &[[5., 0.], [2., 5.], [-4., 3.], [-4., -3.], [2., -5.]],
        false,
    );
    let body = owned_and_bound(input.clone())
        .certify_axial(1e-7, 1e-6, None)
        .unwrap()
        .certify_caps(None)
        .unwrap()
        .certify_ordering(None)
        .unwrap()
        .certify_material(1e-7, None)
        .unwrap();
    for origin in [[100., -30., 50.], [1e6, -2e6, 3e6]] {
        let placed = body
            .transformed_with_interrupt(
                [origin, [0.8, 0.6, 0.], [-0.36, 0.48, 0.8]],
                1e-5,
                1e-4,
                None,
            )
            .unwrap();
        assert_eq!(placed.input().faces.len(), 7);
        let volume = placed.volume_with_interrupt(1e-4, None).unwrap();
        assert!(volume.volume_mm3[0] <= 252. && 252. <= volume.volume_mm3[1]);
        let point = std::array::from_fn(|k| origin[k] + 2. * [0.48, -0.64, 0.6][k]);
        assert_eq!(placed.classify(point, None).unwrap(), Location::Inside);
        assert_eq!(body.input(), &input);
    }
}
fn references(line: &str) -> Vec<usize> {
    line.split('#')
        .skip(1)
        .map(|s| {
            s.chars()
                .take_while(char::is_ascii_digit)
                .collect::<String>()
                .parse()
                .unwrap()
        })
        .collect()
}
#[test]
fn foreign_seven_face_step_needs_cold_physical_proofs_and_rejects_same_count_tampering() {
    let input = graph(
        &[[5., 0.], [2., 5.], [-4., 3.], [-4., -3.], [2., -5.]],
        false,
    );
    let target = input
        .surfaces
        .iter()
        .find_map(|s| {
            if let Surface::Tensor(t) = s {
                Some(t.poles()[55])
            } else {
                None
            }
        })
        .unwrap();
    let body = owned_and_bound(input)
        .certify_axial(1e-7, 1e-6, None)
        .unwrap()
        .certify_caps(None)
        .unwrap()
        .certify_ordering(None)
        .unwrap()
        .certify_material(1e-7, None)
        .unwrap();
    let text = crate::step::write_radial_boundary(&body, 1e-7)
        .unwrap()
        .replace("CARTESIAN_POINT('',", "CARTESIAN_POINT('supplier point',")
        .replace("MANIFOLD_SOLID_BREP('", "MANIFOLD_SOLID_BREP('supplier ");
    assert_eq!(text.matches("ADVANCED_FACE(").count(), 7);
    let candidate = source_candidate(&text);
    assert_eq!(candidate.input().faces.len(), 7);
    let decoded = candidate
        .bind_surfaces(1e-7, None)
        .unwrap()
        .certify_axial(1e-7, 1e-6, None)
        .unwrap();
    assert_eq!(decoded.start_radius_bounds_mm(), None);
    let decoded = decoded
        .certify_caps(None)
        .unwrap()
        .certify_ordering(None)
        .unwrap()
        .certify_material(1e-7, None)
        .unwrap();
    let volume = decoded.volume_with_interrupt(1e-5, None).unwrap();
    assert!(volume.volume_mm3[0] <= 252. && 252. <= volume.volume_mm3[1]);

    let mesh = decoded
        .tessellate_with_interrupt(
            crate::display::Tessellation {
                deflection: 0.1,
                angle: 0.3,
                triangle_limit: 100_000,
            },
            None,
        )
        .unwrap();
    let diagnostic = crate::manufacturing::certify(&mesh).unwrap();
    assert!((diagnostic.signed_volume_mm3 - 252.).abs() < 1.);
    assert_eq!(
        mesh.face_indices
            .iter()
            .copied()
            .collect::<std::collections::BTreeSet<_>>(),
        (0..7).collect()
    );

    let control = text
        .lines()
        .find(|line| {
            if !line.contains("CARTESIAN_POINT(") {
                return false;
            }
            let coords = line
                .rsplit('(')
                .next()
                .unwrap()
                .split(')')
                .next()
                .unwrap()
                .split(',')
                .map(str::parse::<f64>)
                .collect::<Result<Vec<_>, _>>();
            matches!(coords,Ok(p) if p.as_slice()==target)
        })
        .unwrap();
    let id = references(control)[0];
    let folded = text.replacen(
        control,
        &format!(
            "#{id}=CARTESIAN_POINT('supplier point',(-20.,-20.,{}));",
            target[2]
        ),
        1,
    );
    let folded = source_candidate(&folded);
    assert_eq!(folded.input().faces.len(), 7);
    let folded = folded.bind_surfaces(1e-7, None).unwrap();
    assert!(matches!(
        folded.certify_axial(1e-7, 1e-6, None),
        Err(Error::Unresolved)
    ));

    let plane = text.lines().find(|line| line.contains("=PLANE(")).unwrap();
    let axis = references(plane)[1];
    let axis = text
        .lines()
        .find(|line| line.starts_with(&format!("#{axis}=")))
        .unwrap();
    let point = references(axis)[1];
    let original = text
        .lines()
        .find(|line| line.starts_with(&format!("#{point}=")))
        .unwrap();
    let unsupported = text.replacen(
        original,
        &format!("#{point}=CARTESIAN_POINT('supplier point',(0.,0.,0.25));"),
        1,
    );
    let unsupported = source_candidate(&unsupported);
    assert_eq!(unsupported.input().faces.len(), 7);
    assert!(unsupported.bind_surfaces(1e-7, None).is_err());
}
#[test]
fn five_positive_strips_keep_actual_exposed_cone_and_incidence_collars() {
    let input = graph(
        &[
            [3., -4.],
            [6., -3.],
            [7., -1.],
            [7., 1.],
            [6., 3.],
            [3., 4.],
        ],
        true,
    );
    let original = input.clone();
    assert_eq!(
        (input.vertices.len(), input.edges.len(), input.faces.len()),
        (14, 20, 8)
    );
    strict_source_proofs(&input, 4.9);
    let bound = owned_and_bound(input.clone());
    assert_eq!(bound.input(), &input);
    let axial = bound.certify_axial(1e-7, 1e-6, None).unwrap();
    assert!(axial.start_radius_bounds_mm().is_some());
    let caps = axial.certify_caps(None).unwrap();
    assert_eq!(caps.tensor_rail_incidences, 10);
    let body = caps
        .certify_ordering(None)
        .unwrap()
        .certify_material(1e-7, None)
        .unwrap();
    let expected = 4. * (37. + 25. / 2. * (std::f64::consts::TAU - 2. * 4_f64.atan2(3.)));
    let volume = body.volume_with_interrupt(1e-5, None).unwrap();
    assert!(volume.volume_mm3[0] <= expected && expected <= volume.volume_mm3[1]);
    assert_eq!(input, original);
}
#[test]
fn five_positive_strips_without_cone_qualify_only_with_actual_axis_and_core_proof() {
    let input = graph(
        &[[5., 0.], [2., 5.], [-4., 3.], [-4., -3.], [2., -5.]],
        false,
    );
    let original = input.clone();
    assert_eq!(
        (input.vertices.len(), input.edges.len(), input.faces.len()),
        (10, 15, 7)
    );
    assert!(
        !input
            .surfaces
            .iter()
            .any(|s| matches!(s, Surface::Cone { .. }))
    );
    strict_source_proofs(&input, 2.);
    let bound = owned_and_bound(input.clone());
    assert_eq!(bound.input(), &input);
    let axial = bound.certify_axial(1e-7, 1e-6, None).unwrap();
    assert_eq!(axial.start_radius_bounds_mm(), None);
    let caps = axial.certify_caps(None).unwrap();
    assert_eq!(caps.tensor_rail_incidences, 10);
    let body = caps
        .certify_ordering(None)
        .unwrap()
        .certify_material(1e-7, None)
        .unwrap();
    let volume = body.volume_with_interrupt(1e-5, None).unwrap();
    assert!(volume.volume_mm3[0] <= 252. && 252. <= volume.volume_mm3[1]);
    let text = crate::step::write_radial_boundary(&body, 1e-7).unwrap();
    assert!(!text.contains("CONICAL_SURFACE"));
    let decoded = crate::step_read::read_radial_boundary(&text, 1e-7)
        .unwrap()
        .bind_surfaces(1e-7, None)
        .unwrap()
        .certify_axial(1e-7, 1e-6, None)
        .unwrap();
    assert_eq!(decoded.start_radius_bounds_mm(), None);
    let decoded = decoded
        .certify_caps(None)
        .unwrap()
        .certify_ordering(None)
        .unwrap()
        .certify_material(1e-7, None)
        .unwrap();
    let volume = decoded.volume_with_interrupt(1e-5, None).unwrap();
    assert!(volume.volume_mm3[0] <= 252. && 252. <= volume.volume_mm3[1]);
    assert_eq!(input, original);
}

#[test]
fn positive_cap_cycle_and_valid_incidence_do_not_admit_an_interior_fold() {
    let mut input = graph(
        &[[5., 0.], [2., 5.], [-4., 3.], [-4., -3.], [2., -5.]],
        false,
    );
    let Surface::Tensor(t) = &input.surfaces[2] else {
        panic!("actual face");
    };
    let mut poles = t.poles().to_vec();
    poles[5 * 10 + 5][0] = -20.;
    poles[5 * 10 + 5][1] = -20.;
    let folded = RetainedTensor::new([1, 1], poles).unwrap();
    assert_eq!(
        folded.certify_radial_chart(true, None),
        Err(Error::Unresolved)
    );
    input.surfaces[2] = Surface::Tensor(folded);
    let bound = owned_and_bound(input.clone());
    assert_eq!(bound.input(), &input);
    assert!(matches!(
        bound.certify_axial(1e-7, 1e-6, None),
        Err(Error::Unresolved)
    ));
}
#[test]
fn locally_positive_five_bands_cannot_publish_a_double_winding_cap() {
    let profile = (0..5)
        .map(|i| {
            let angle = std::f64::consts::TAU * (2 * i) as f64 / 5.;
            [5. * angle.cos(), 5. * angle.sin()]
        })
        .collect::<Vec<_>>();
    let input = graph(&profile, false);
    let axial = owned_and_bound(input)
        .certify_axial(1e-7, 1e-6, None)
        .unwrap();
    assert_eq!(axial.certify_caps(None).err(), Some(Error::Unresolved));
}
