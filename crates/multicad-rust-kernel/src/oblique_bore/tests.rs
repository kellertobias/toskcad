use super::*;
use crate::manufacturing;
use std::{cell::Cell, f64::consts::PI};

const OPTIONS: Tessellation = Tessellation {
    deflection: 0.002,
    angle: 0.2,
    triangle_limit: 2_000_000,
};

fn normalized(v: V) -> V {
    let length = dot(v, v).sqrt();
    v.map(|x| x / length)
}

/// A right-handed tool frame around a travel direction.
fn tool(base_mm: V, travel: V, radius_mm: f64, height_mm: f64) -> Tool {
    let w = normalized(travel);
    let seed = if w[2].abs() < 0.9 { unit(2) } else { unit(0) };
    let e1 = normalized(cross(seed, w));
    Tool {
        base_mm,
        axes: [e1, cross(w, e1), w],
        radius_mm,
        height_mm,
        literal: true,
    }
}

fn body(size: V, tool: Tool) -> ObliqueBore {
    match subtract(size, tool, None) {
        Ok(Cut::Body(body)) => *body,
        other => panic!("expected a cut body, got {other:?}"),
    }
}

/// Midpoint quadrature of the removed volume over the tool cross-section.
fn removed_volume(size: V, tool: Tool) -> f64 {
    let [e1, e2, w] = tool.axes;
    let (rings, turns) = (1200, 2400);
    let mut volume = 0.;
    for i in 0..rings {
        let radius = tool.radius_mm * (i as f64 + 0.5) / rings as f64;
        for j in 0..turns {
            let (sin, cos) = (TAU * (j as f64 + 0.5) / turns as f64).portable_sin_cos();
            let start = mix(mix(tool.base_mm, radius * cos, e1), radius * sin, e2);
            let mut range = [0., tool.height_mm];
            let mut inside = true;
            for k in 0..3 {
                let bounds = [(0. - start[k]) / w[k], (size[k] - start[k]) / w[k]];
                if w[k] == 0. {
                    inside &= (0. ..=size[k]).contains(&start[k]);
                } else {
                    range[0] = range[0].max(bounds[0].min(bounds[1]));
                    range[1] = range[1].min(bounds[0].max(bounds[1]));
                }
            }
            if inside && range[1] > range[0] {
                volume += (range[1] - range[0]) * radius;
            }
        }
    }
    volume * tool.radius_mm / rings as f64 * TAU / turns as f64
}

/// Check every invariant shared by admitted bodies and return the mesh volume.
fn certify(cut: &ObliqueBore, exact_removed: Option<f64>) -> f64 {
    let size = cut.size_mm();
    let tool = cut.tool();
    let stock = size[0] * size[1] * size[2];
    let (removed, slack) = match exact_removed {
        Some(removed) => (removed, 0.),
        None => (removed_volume(size, tool), 4e-5),
    };
    assert!(
        (stock - cut.volume_mm3() - removed).abs() < slack * removed + cut.volume_error_mm3(),
        "analytic {} removed {}",
        stock - cut.volume_mm3(),
        removed
    );
    for (index, edge) in cut.edges().iter().enumerate() {
        for (side, &face) in edge.faces.iter().enumerate() {
            let face = &cut.faces()[face];
            assert_eq!(face.support, edge.supports[side]);
            let uses = face.loops.iter().flatten().filter(|u| u.0 == index);
            assert_eq!(uses.count(), 1);
        }
        let [first, second] = edge.faces.map(|face| {
            cut.faces()[face]
                .loops
                .iter()
                .flatten()
                .find(|u| u.0 == index)
                .unwrap()
                .1
        });
        assert_ne!(first, second);
        if let Some(arc) = edge.curve {
            let ends = [arc.start_rad, arc.start_rad + arc.sweep_rad].map(|a| arc.point(a));
            let principal = arc.principal();
            for (end, angle) in ends.into_iter().zip([0., arc.sweep_rad]) {
                let (sin, cos) = (principal.start_rad + angle).portable_sin_cos();
                let again = mix(
                    mix(
                        principal.center_mm,
                        principal.semiaxes_mm[0] * cos,
                        principal.u,
                    ),
                    principal.semiaxes_mm[1] * sin,
                    principal.v,
                );
                assert!(dot(sub(end, again), sub(end, again)).sqrt() < 1e-9);
            }
            assert!(dot(principal.u, principal.v).abs() < 1e-9);
            assert!(principal.semiaxes_mm[0] >= principal.semiaxes_mm[1]);
            for (end, vertex) in ends.into_iter().zip(edge.ends) {
                let gap = sub(end, cut.vertices()[vertex]);
                assert!(dot(gap, gap).sqrt() < 1e-9, "arc end {gap:?}");
            }
        }
    }
    let area: f64 = cut.faces().iter().map(|face| face.area_mm2).sum();
    let wall: f64 = cut
        .faces()
        .iter()
        .filter(|face| face.support == Support::ToolSide)
        .map(|face| face.area_mm2)
        .sum();
    assert!(cut.plane(Support::ToolSide).is_none());
    let mesh = cut.mesh(OPTIONS).unwrap();
    let proof = manufacturing::certify(&mesh).unwrap();
    assert!(
        (proof.signed_volume_mm3 - cut.volume_mm3()).abs()
            <= OPTIONS.deflection * wall + 1e-9 * area,
        "mesh {} analytic {}",
        proof.signed_volume_mm3,
        cut.volume_mm3()
    );
    assert!(proof.signed_volume_mm3 >= cut.volume_mm3() - 1e-9 * area);
    let vertices = mesh.vertices.as_chunks::<3>().0;
    let mut areas = vec![0.; cut.faces().len()];
    let mut centers = vec![[0.; 3]; cut.faces().len()];
    for (triangle, &face) in mesh
        .indices
        .as_chunks::<3>()
        .0
        .iter()
        .zip(&mesh.face_indices)
    {
        let [a, b, c] = triangle.map(|i| vertices[i as usize]);
        let normal = cross(sub(b, a), sub(c, a));
        let piece = dot(normal, normal).sqrt() / 2.;
        areas[face as usize] += piece;
        centers[face as usize] = mix(
            centers[face as usize],
            piece / 3.,
            mix(mix(a, 1., b), 1., c),
        );
        let face = &cut.faces()[face as usize];
        if let Some((point, outward)) = cut.plane(face.support) {
            for corner in [a, b, c] {
                assert!(dot(sub(corner, point), outward).abs() < 1e-9);
            }
            assert!(dot(normal, outward) > 0.);
        } else {
            let middle = mix(mix(a, 1., b), 1., c).map(|x| x / 3.);
            let [x, y, _] = cut.chart.local(middle);
            assert!((x.portable_hypot(y) - tool.radius_mm).abs() <= OPTIONS.deflection);
            assert!(dot(normal, cut.chart.normal(WALL, middle)) > 0.);
        }
    }
    for (face, (area, center)) in cut.faces().iter().zip(areas.iter().zip(centers)) {
        assert!(
            (face.area_mm2 - area).abs() < 3e-3 * face.area_mm2 + 1e-9,
            "{:?} area {} mesh {}",
            face.support,
            face.area_mm2,
            area
        );
        if face.support != Support::ToolSide {
            let gap = sub(face.center_mm, center.map(|c| c / area));
            assert!(dot(gap, gap).sqrt() < 2e-3 * size.into_iter().fold(0., f64::max));
        }
    }
    for direction in [unit(0), unit(1), unit(2), normalized([0.3, -0.5, 0.8])] {
        let [low, high] = cut.extent(direction);
        let sampled = vertices.iter().map(|p| dot(*p, direction));
        let reached = [
            sampled.clone().fold(f64::INFINITY, f64::min),
            sampled.fold(f64::NEG_INFINITY, f64::max),
        ];
        assert!(low <= reached[0] + 1e-12 && reached[0] - low < 0.01);
        assert!(high >= reached[1] - 1e-12 && high - reached[1] < 0.01);
    }
    proof.signed_volume_mm3
}

fn supports(cut: &ObliqueBore) -> Vec<(Support, usize)> {
    cut.faces()
        .iter()
        .map(|face| (face.support, face.loops.len()))
        .collect()
}

const SIZE: V = [10., 8., 6.];

fn audit_tool() -> Tool {
    let s = 2f64.sqrt().recip();
    let t = 6f64.sqrt().recip();
    let n = 3f64.sqrt().recip();
    let (u, v) = ([s, -s, 0.], [t, t, -2. * t]);
    Tool {
        base_mm: [1. / 3. - n; 3],
        axes: [u, v, cross(u, v)],
        radius_mm: 0.1,
        height_mm: 20.,
        literal: true,
    }
}

#[test]
fn audit_fixture_keeps_corner_trims_and_the_cylinder_support() {
    let cut = body([10.; 3], audit_tool());
    // Both corners on the diagonal are removed. Each stock face keeps one
    // loop with an ellipse trim and the bore is one face with two rims.
    let stock_faces =
        (0..3).flat_map(|axis| [false, true].map(|high| Support::Stock { axis, high }));
    let expected: Vec<_> = stock_faces
        .map(|face| (face, 1))
        .chain([(Support::ToolSide, 2)])
        .collect();
    assert_eq!(supports(&cut), expected);
    assert_eq!(
        [cut.vertices().len(), cut.edges().len(), cut.faces().len()],
        [12, 18, 7]
    );
    assert_eq!(cut.euler_characteristic(), 0);
    let arcs: Vec<_> = cut.edges().iter().filter_map(|edge| edge.curve).collect();
    assert_eq!(arcs.len(), 6);
    for arc in arcs {
        let principal = arc.principal();
        assert!((principal.semiaxes_mm[0] - 0.1 * 3f64.sqrt()).abs() < 1e-12);
        assert!((principal.semiaxes_mm[1] - 0.1).abs() < 1e-12);
        assert!(arc.sweep_rad > 2. && arc.sweep_rad < 2.2);
    }
    for direction in 0..3 {
        let [low, high] = cut.extent(unit(direction));
        assert!(low.abs() < 1e-12 && (high - 10.).abs() < 1e-12);
    }
    assert!(cut.resolution_mm() < 2e-8 && cut.volume_error_mm3() < 1e-5);
    let wall = cut.faces().last().unwrap();
    // The bore runs corner to corner; its centroid stays on the diagonal.
    assert!(wall.center_mm.iter().all(|c| (c - 5.).abs() < 1e-9));
    // Each rim is the upper envelope of three equal sinusoids in the angle.
    let mean = 10. * 3f64.sqrt() - 0.1 * 3. * 6f64.sqrt() / PI;
    assert!((wall.area_mm2 - TAU * 0.1 * mean).abs() < 1e-9);
    certify(&cut, None);
}

#[test]
fn through_and_blind_cuts_match_exact_oblique_cylinder_volumes() {
    let through = tool([3., 4., -2.], [0.2, 0.1, 1.], 1., 20.);
    let cut = body(SIZE, through);
    assert_eq!(cut.euler_characteristic(), 0);
    assert_eq!(
        supports(&cut)[4..],
        [
            (
                Support::Stock {
                    axis: 2,
                    high: false
                },
                2
            ),
            (
                Support::Stock {
                    axis: 2,
                    high: true
                },
                2
            ),
            (Support::ToolSide, 2),
        ]
    );
    // A full oblique section removes its axis length times the disk area.
    certify(&cut, Some(PI * 6. / through.axes[2][2]));

    // Reversing the travel swaps only the tool cap identities.
    let blind = tool([5., 4., 8.], [0.2, -0.1, -1.], 1., 5.);
    let end = mix(blind.base_mm, blind.height_mm, blind.axes[2]);
    let reversed = tool(end, blind.axes[2].map(|c| -c), 1., 5.);
    let removed = PI * (5. + 2. / blind.axes[2][2]);
    let mut volumes = Vec::new();
    for (tool, cap) in [(blind, Support::ToolEnd), (reversed, Support::ToolStart)] {
        let cut = body(SIZE, tool);
        assert_eq!(cut.euler_characteristic(), 2);
        // The extended cap plane adds no vertex where it leaves the disk.
        assert_eq!([cut.vertices().len(), cut.edges().len()], [10, 14]);
        let floor = cut.faces().iter().find(|face| face.support == cap).unwrap();
        assert!((floor.area_mm2 - PI).abs() < 1e-12);
        let gap = sub(floor.center_mm, end);
        assert!(dot(gap, gap).sqrt() < 1e-12);
        let (point, normal) = cut.plane(cap).unwrap();
        assert!(dot(sub(point, end), sub(point, end)).sqrt() < 1e-12);
        // The floor faces back along the removed tool.
        assert!(dot(normal, sub(blind.base_mm, end)) > 0.);
        assert_eq!(cut.faces().len(), 8);
        volumes.push(certify(&cut, Some(removed)));
    }
    assert!((volumes[0] - volumes[1]).abs() < 1e-6);
}

#[test]
fn partial_caps_generators_and_planar_cuts_keep_exact_trims() {
    // The floor of this blind cut is clipped by a stock face in a chord.
    let cut = body(SIZE, tool([1., 4., 8.], [-0.3, 0., -1.], 1., 5.));
    let floor = cut
        .faces()
        .iter()
        .find(|face| face.support == Support::ToolEnd)
        .unwrap();
    let kinds: Vec<bool> = floor.loops[0]
        .iter()
        .map(|u| cut.edges()[u.0].curve.is_some())
        .collect();
    assert_eq!(kinds.iter().filter(|arc| **arc).count(), 1);
    assert_eq!(kinds.len(), 2);
    assert_eq!(cut.euler_characteristic(), 2);
    certify(&cut, None);

    // An axis parallel to a face pair opens a slot bounded by generators.
    let slot = body(SIZE, tool([-2., -1., 5.5], [1., 1., 0.], 1., 30.));
    let top = Support::Stock {
        axis: 2,
        high: true,
    };
    assert_eq!(
        slot.faces()
            .iter()
            .filter(|face| face.support == top)
            .count(),
        2
    );
    let generators = slot
        .edges()
        .iter()
        .filter(|edge| edge.curve.is_none() && edge.supports == [top, Support::ToolSide]);
    assert_eq!(generators.count(), 2);
    assert_eq!(slot.euler_characteristic(), 2);
    certify(&slot, None);

    // An axis lying in the bottom face cuts a half-round groove.
    let groove = body(SIZE, tool([-2., -1., 0.], [1., 1., 0.], 1., 30.));
    let bottom = Support::Stock {
        axis: 2,
        high: false,
    };
    let generators = groove
        .edges()
        .iter()
        .filter(|edge| edge.curve.is_none() && edge.supports == [bottom, Support::ToolSide]);
    assert_eq!(generators.count(), 2);
    certify(&groove, None);

    // The same axis inside the slab leaves a hole through two side faces,
    // one of them across a stock edge.
    let hole = body(SIZE, tool([-2., -1., 3.], [1., 1., 0.], 1., 30.));
    assert_eq!(hole.euler_characteristic(), 0);
    certify(&hole, None);
    let notch = body(SIZE, tool([-2., 6., 3.], [1., 1., 0.], 1., 30.));
    assert_eq!(notch.euler_characteristic(), 2);
    // The axis meets a stock edge at a right angle: twice the integral of
    // the offset over a half disk.
    certify(&notch, Some(4. / 3.));

    // A wide tool cap alone cuts the corner tetrahedron x + y + z >= 21.
    let chamfer = body(SIZE, tool([9., 7., 5.], [1., 1., 1.], 50., 30.));
    assert!(chamfer.edges().iter().all(|edge| edge.curve.is_none()));
    let cap = chamfer.faces().last().unwrap();
    assert_eq!((cap.support, cap.loops[0].len()), (Support::ToolStart, 3));
    assert!((cap.area_mm2 - 4.5 * 3f64.sqrt()).abs() < 1e-12);
    assert!(
        cap.center_mm
            .iter()
            .zip([8., 6., 4.])
            .all(|(c, e)| (c - e - 1.).abs() < 1e-12)
    );
    certify(&chamfer, Some(4.5));

    // Larger coordinates keep the same relative resolution.
    let large = body(
        [5000., 4000., 3000.],
        tool([1500., 2000., -700.], [0.2, 0.1, 1.], 400., 9000.),
    );
    assert!((large.resolution_mm() - 9000. * RESOLUTION).abs() < 1e-12);
    let removed = PI * 160_000. * 3000. / large.tool().axes[2][2];
    assert!(
        (5000. * 4000. * 3000. - large.volume_mm3() - removed).abs() < large.volume_error_mm3()
    );
}

fn planar_tool(base_mm: V) -> Tool {
    // Exact binary columns: the axis is parallel to the z faces.
    let s = 2f64.sqrt().recip();
    let (u, v) = ([s, -s, 0.], [0., 0., 1.]);
    Tool {
        base_mm,
        axes: [u, v, cross(u, v)],
        radius_mm: 1.,
        height_mm: 30.,
        literal: true,
    }
}

fn rim_tool(x: f64) -> Tool {
    let (u, v) = ([0.6, 0.8, 0.], [0., 0., -1.]);
    Tool {
        base_mm: [x, 4., 3.],
        axes: [u, v, cross(u, v)],
        radius_mm: 1.,
        height_mm: 5.,
        literal: true,
    }
}

#[test]
fn disjoint_and_true_contact_keep_the_stock_but_positive_material_is_removed() {
    let unchanged = |tool| matches!(subtract(SIZE, tool, None), Ok(Cut::Unchanged));
    // Separated by a stock plane, and disjoint beside a stock edge where
    // no stock plane separates.
    assert!(unchanged(tool([30., 4., 3.], [1., 1., 1.], 1., 5.)));
    assert!(unchanged(tool([4., 18., 3.], [1., -1., 0.], 0.5, 12.)));
    // The side touches the top face along one generator.
    assert!(unchanged(planar_tool([12., 9., 7.])));
    // The start rim touches the face x = 0 in one point.
    assert_eq!(rim_tool(-0.6).base_mm[0] + 0.6, 0.);
    assert!(unchanged(rim_tool(-0.6)));

    // A rounded chart of the same operands cannot certify contact, but it
    // still certifies clear separation by margin.
    for touching in [planar_tool([12., 9., 7.]), rim_tool(-0.6)] {
        let rounded = Tool {
            literal: false,
            ..touching
        };
        assert_eq!(subtract(SIZE, rounded, None).err(), Some(Error::Unresolved));
    }
    let rounded = Tool {
        literal: false,
        ..tool([30., 4., 3.], [1., 1., 1.], 1., 5.)
    };
    assert!(unchanged(rounded));
    // The smallest representable penetration is never treated as contact.
    for tool in [
        planar_tool([12., 9., 7f64.next_down()]),
        rim_tool((-0.6f64).next_up()),
    ] {
        assert_eq!(subtract(SIZE, tool, None).err(), Some(Error::Unresolved));
    }
    // A resolvable sliver is removed.
    for tool in [planar_tool([12., 9., 7. - 1e-4]), rim_tool(-0.6 + 1e-4)] {
        let cut = body(SIZE, tool);
        let removed = 480. - cut.volume_mm3();
        assert!(removed > 1e-10 && removed < 1e-4, "{removed}");
        let mesh = cut.mesh(OPTIONS).unwrap();
        let proof = manufacturing::certify(&mesh).unwrap();
        let wall = cut.faces().last().unwrap();
        assert_eq!(wall.support, Support::ToolSide);
        let chord = OPTIONS.deflection * wall.area_mm2;
        assert!((proof.signed_volume_mm3 - cut.volume_mm3()).abs() <= chord + 1e-9);
        assert!(proof.signed_volume_mm3 < 480.);
    }
    // A tool containing the stock leaves nothing.
    let cover = tool([-5.; 3], [1., 1., 1.], 10., 30.);
    assert!(matches!(subtract([1.; 3], cover, None), Ok(Cut::Empty)));
}

#[test]
fn refuses_invalid_unsupported_and_unresolved_sources() {
    let good = tool([3., 4., -2.], [0.2, 0.1, 1.], 1., 20.);
    let refusal = |size, tool| subtract(size, tool, None).err();
    let [u, v, w] = good.axes;
    let invalid = [
        ([0., 8., 6.], good),
        ([f64::NAN, 8., 6.], good),
        ([2e7, 8., 6.], good),
        (
            SIZE,
            Tool {
                radius_mm: -1.,
                ..good
            },
        ),
        (
            SIZE,
            Tool {
                radius_mm: 1e-6,
                ..good
            },
        ),
        (
            SIZE,
            Tool {
                height_mm: f64::INFINITY,
                ..good
            },
        ),
        (
            SIZE,
            Tool {
                base_mm: [f64::NAN, 0., 0.],
                ..good
            },
        ),
        (
            SIZE,
            Tool {
                base_mm: [0., -2e7, 0.],
                ..good
            },
        ),
        (
            SIZE,
            Tool {
                axes: [u, v.map(|c| c * 1.001), w],
                ..good
            },
        ),
        (
            SIZE,
            Tool {
                axes: [v, u, w],
                ..good
            },
        ),
        (
            SIZE,
            Tool {
                axes: [u, mix(v, 1e-6, u), w],
                ..good
            },
        ),
        (
            SIZE,
            Tool {
                axes: [u, [f64::NAN; 3], w],
                ..good
            },
        ),
    ];
    for (size, tool) in invalid {
        assert_eq!(refusal(size, tool), Some(Error::Invalid), "{tool:?}");
    }
    let unsupported = [
        // Axis-aligned tools belong to the orthogonal provider.
        tool([3., 4., -2.], [0., 0., 1.], 1., 20.),
        // An enclosed void and a cut that splits the stock.
        tool([3., 3., 3.], [1., 1., 1.], 0.5, 2.),
        tool([-5., -6., 3.], [1., 1., 0.], 3.5, 30.),
    ];
    for tool in unsupported {
        assert_eq!(refusal(SIZE, tool), Some(Error::Unsupported), "{tool:?}");
    }
    let s = 2f64.sqrt();
    let unresolved = [
        // Nearly parallel to a face pair.
        tool([3., 4., -2.], [1., 1., 1e-6], 1., 20.),
        // The side grazes a stock edge from outside.
        tool([10. + s, 0., 3.], [-1., 1., 0.], 1., 30.),
        // The start cap plane contains a stock edge.
        tool([0., 0., 3.], [1., 1., 0.], 1., 5.),
        // The wall passes through a stock corner.
        tool([-5.; 3], [1., 1., 1.], 10. * (2f64 / 3.).sqrt(), 30.),
    ];
    for tool in unresolved {
        assert_eq!(refusal(SIZE, tool), Some(Error::Unresolved), "{tool:?}");
    }
    // A rim that grazes an extended stock plane clearly outside the stock
    // does not block the cut or a certified separation.
    let grazing = tool([-3., 4., 6.6], [0.6, 0., -0.8], 1., 20.);
    let rim = mix(grazing.base_mm, -1., cross(grazing.axes[2], unit(1)));
    assert!((rim[2] - 6.).abs() < 1e-12 && rim[0] < -3.5);
    assert_eq!(body(SIZE, grazing).euler_characteristic(), 0);
    let outside = Tool {
        base_mm: [-0.6, 20., 3.],
        literal: false,
        ..rim_tool(-0.6)
    };
    assert!(matches!(subtract(SIZE, outside, None), Ok(Cut::Unchanged)));
    // The start cap is centered on a stock corner.
    let corner = tool([0.; 3], [1., 1., 1.], 0.1, 20.);
    assert_eq!(refusal([10.; 3], corner), Some(Error::Unresolved));
    // The same cap plane through a corner far outside its disk is no vertex.
    let remote = tool([1500., 2000., -500.], [0.2, 0.1, 1.], 400., 9000.);
    let large = [5000., 4000., 3000.];
    assert_eq!(dot(remote.base_mm, [0.2, 0.1, 1.]), 0.);
    let cut = body(large, remote);
    assert_eq!([cut.vertices().len(), cut.edges().len()], [10, 14]);
}

#[test]
fn cancellation_and_mesh_budgets_refuse_without_a_partial_result() {
    let calls = Cell::new(0);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    assert!(subtract([10.; 3], audit_tool(), Some(&count)).is_ok());
    let total = calls.get();
    assert!(total > 30);
    for stop in [1, 2, total / 2, total] {
        calls.set(0);
        let interrupt = || {
            calls.set(calls.get() + 1);
            calls.get() == stop
        };
        let result = subtract([10.; 3], audit_tool(), Some(&interrupt));
        assert_eq!(result.err(), Some(Error::Cancelled));
        assert_eq!(calls.get(), stop);
    }

    let cut = body(SIZE, tool([-2., -1., 5.5], [1., 1., 0.], 1., 30.));
    let full = cut.mesh(OPTIONS).unwrap();
    let triangles = full.face_indices.len();
    let limited = |triangle_limit| {
        cut.mesh(Tessellation {
            triangle_limit,
            ..OPTIONS
        })
    };
    assert_eq!(limited(triangles).unwrap(), full);
    for limit in [1, triangles / 8, triangles - 1] {
        assert_eq!(limited(limit), Err(MeshError::BudgetExceeded));
    }
    let invalid = Tessellation {
        deflection: 0.,
        ..OPTIONS
    };
    assert_eq!(cut.mesh(invalid), Err(MeshError::InvalidOptions));
    // The retained body is still exact after every refusal.
    certify(&cut, None);
}

#[test]
fn a_cap_plane_through_a_distant_stock_edge_is_still_resolved() {
    // A blind pocket tilted in XZ: its start-cap plane contains the stock
    // edge x = 0, z = 6, five millimetres from the disk, and its end cap
    // dips back into the top face.
    let pocket = Tool {
        base_mm: [4., 4., 3.],
        axes: [[0.8, 0., -0.6], [0., 1., 0.], [0.6, 0., 0.8]],
        radius_mm: 2.,
        height_mm: 5.,
        literal: true,
    };
    let cut = body(SIZE, pocket);
    assert_eq!(cut.euler_characteristic(), 2);
    let stock = SIZE[0] * SIZE[1] * SIZE[2];
    let removed = removed_volume(SIZE, pocket);
    assert!((stock - cut.volume_mm3() - removed).abs() < 4e-5 * removed);
    // The end cap keeps the circular segment below the top face: the chord
    // lies 5/3 mm from the cap centre of the radius-2 disk.
    let segment = 4. * (5f64 / 6.).portable_acos() - 5. / 3. * (11f64 / 9.).sqrt();
    let floor = cut
        .faces()
        .iter()
        .find(|face| face.support == Support::ToolEnd)
        .unwrap();
    assert!(
        (floor.area_mm2 - segment).abs() < 1e-9,
        "{}",
        floor.area_mm2
    );
    let start = cut
        .faces()
        .iter()
        .find(|face| face.support == Support::ToolStart)
        .unwrap();
    assert!((start.area_mm2 - 4. * PI).abs() < 1e-9);
    // The same coincidence with the edge crossing the disk stays refused.
    let crossing = Tool {
        base_mm: [0.8, 4., 5.4],
        ..pocket
    };
    assert_eq!(
        subtract(SIZE, crossing, None).err(),
        Some(Error::Unresolved)
    );
}
