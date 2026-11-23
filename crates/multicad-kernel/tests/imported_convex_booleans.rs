use multicad_core::{document::EntityId, references::*};
use multicad_kernel::{
    rust_kernel::{
        Rational,
        planar_boolean::{ConvexPiece, Union},
        step,
        step_read::{ImportProfile, read_import_solid},
    },
    *,
};

const TOL: f64 = 1e-7;
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 3e-5, "{a} != {b}");
}
/// Compare analytic geometry, ignoring the direction of a line.
fn same_geometry(a: AnalyticGeometry, b: AnalyticGeometry) {
    match (a, b) {
        (
            AnalyticGeometry::Line { a_mm, b_mm },
            AnalyticGeometry::Line {
                a_mm: start,
                b_mm: end,
            },
        ) => {
            let (start, end) = if (0..3).all(|i| (a_mm[i] - start[i]).abs() < 3e-5) {
                (start, end)
            } else {
                (end, start)
            };
            for i in 0..3 {
                near(a_mm[i], start[i]);
                near(b_mm[i], end[i]);
            }
        }
        (
            AnalyticGeometry::Point { position_mm: a },
            AnalyticGeometry::Point { position_mm: b },
        ) => (0..3).for_each(|i| near(a[i], b[i])),
        (AnalyticGeometry::Plane { frame: a }, AnalyticGeometry::Plane { frame: b }) => {
            (0..3).for_each(|i| near(a.normal()[i], b.normal()[i]));
            near(
                (0..3)
                    .map(|i| a.normal()[i] * (a.origin_mm[i] - b.origin_mm[i]))
                    .sum(),
                0.,
            );
        }
        other => panic!("changed geometry {other:?}"),
    }
}
fn tetra() -> Union {
    let [a, b, c, d] = [[0., 0., 0.], [4., 0., 0.], [0., 4., 0.], [0., 0., 4.]];
    Union::new(
        &[ConvexPiece::new(
            vec![vec![a, c, b], vec![a, b, d], vec![a, d, c], vec![b, c, d]],
            TOL,
        )
        .unwrap()],
        TOL,
    )
    .unwrap()
}
fn placed_text(geometry: &Union, frame: PlaneFrame) -> String {
    step::write_planar_union(geometry, frame.origin, frame.u, frame.v).unwrap()
}
fn transform(frame: PlaneFrame) -> RigidTransform {
    RigidTransform {
        rotation: std::array::from_fn(|i| [frame.u[i], frame.v[i], frame.normal()[i]]),
        translation_mm: frame.origin,
    }
}
fn compare(a: &Solid, b: &Solid) {
    near(a.volume_mm3().unwrap(), b.volume_mm3().unwrap());
    near(
        a.faces().unwrap().iter().map(|f| f.area_mm2).sum(),
        b.faces().unwrap().iter().map(|f| f.area_mm2).sum(),
    );
    let (a, b) = (a.bounds().unwrap(), b.bounds().unwrap());
    for (a, b) in a
        .min
        .into_iter()
        .chain(a.max)
        .zip(b.min.into_iter().chain(b.max))
    {
        near(a, b);
    }
}

#[test]
fn imported_tetrahedron_proves_occupied_convex_material_for_cut_join_and_queries() {
    let policy = GeometryPolicy::default();
    let tilted = PlaneFrame {
        origin: [11., -7., 23.],
        u: [0.8, 0.6, 0.],
        v: [-0.36, 0.48, 0.8],
    };
    for frame in [PlaneFrame::default(), tilted] {
        let text = placed_text(&tetra(), frame);
        assert!(matches!(
            read_import_solid(&text, TOL).unwrap().1.profile,
            ImportProfile::PlanarBoundary { .. }
        ));
        for additional in [
            RigidTransform::default(),
            RigidTransform {
                rotation: [[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]],
                translation_mm: [7., -2., 5.],
            },
        ] {
            let owner = EntityId::new();
            let source = Solid::import_rust_step(&text, owner, policy)
                .unwrap()
                .transformed(additional)
                .unwrap();
            let tool = Solid::rust_box(
                EntityId::new().uuid_bits(),
                [4, 1, 1].map(|n| Rational::new(n, 1).unwrap()),
                PlaneFrame {
                    origin: [1., 0., 0.],
                    ..Default::default()
                },
                policy,
            )
            .unwrap()
            .transformed(transform(frame))
            .unwrap()
            .transformed(additional)
            .unwrap();
            let before = source.topology_snapshot(owner, 0).unwrap();
            let tool_before = tool.topology_snapshot(owner, 0).unwrap();
            let cut = source
                .rust_subtract_planar_tools(&[&tool], EntityId::new(), policy)
                .unwrap()
                .unwrap();
            let join = Solid::union_planar_pieces(&[&source, &tool], EntityId::new(), policy)
                .unwrap()
                .unwrap();
            near(cut.volume_mm3().unwrap(), 26. / 3.);
            near(join.volume_mm3().unwrap(), 38. / 3.);
            let default_cut = source
                .boolean(&tool, BooleanOperation::Subtract, policy)
                .unwrap()
                .unwrap();
            let default_join = source
                .boolean(&tool, BooleanOperation::Union, policy)
                .unwrap()
                .unwrap();
            let default_common = source
                .boolean(&tool, BooleanOperation::Intersect, policy)
                .unwrap()
                .unwrap();
            near(default_cut.volume_mm3().unwrap(), 26. / 3.);
            near(default_join.volume_mm3().unwrap(), 38. / 3.);
            near(default_common.volume_mm3().unwrap(), 2.);
            let common = source
                .rust_planar_intersection_metrics(&tool, policy)
                .unwrap()
                .unwrap();
            near(common.0, 2.);
            let left_normal = additional.rotation.map(|row| {
                row.into_iter()
                    .zip(frame.u)
                    .map(|(a, b)| -a * b)
                    .sum::<f64>()
            });
            let left=before.elements.iter().find(|e|matches!(e.geometry,Some(AnalyticGeometry::Plane {frame})if frame.normal().into_iter().zip(left_normal).all(|(a,b)|(a-b).abs()<1e-9))).unwrap();
            let Some(AnalyticGeometry::Plane { frame: left_plane }) = left.geometry else {
                unreachable!()
            };
            let on_left = |p: [f64; 3]| {
                (0..3)
                    .map(|i| left_plane.normal()[i] * (p[i] - left_plane.origin_mm[i]))
                    .sum::<f64>()
                    .abs()
                    < 1e-8
            };
            let retained: Vec<_> = before
                .elements
                .iter()
                .filter(|element| match element.geometry {
                    Some(AnalyticGeometry::Plane { .. }) => element.key == left.key,
                    Some(AnalyticGeometry::Line { a_mm, b_mm }) => on_left(a_mm) && on_left(b_mm),
                    Some(AnalyticGeometry::Point { position_mm }) => on_left(position_mm),
                    _ => false,
                })
                .map(|element| {
                    (
                        GeometryIntent::capture(&before, element.key).unwrap(),
                        element.geometry,
                    )
                })
                .collect();
            assert_eq!(retained.len(), 7);
            for result in [&cut, &join, &default_cut, &default_join] {
                assert_eq!(result.solid_count().unwrap(), 1);
                let output = result.topology_snapshot(owner, 1).unwrap();
                for (reference, original) in &retained {
                    // The clipping box never reaches x=0. Every boundary of
                    // that face, including the whole A-C edge, stays one
                    // element with its original extent.
                    let Resolution::Resolved { key, .. } = resolve(reference, &output, 1) else {
                        panic!("lost tetrahedron boundary {reference:?}");
                    };
                    let current = output.elements.iter().find(|e| e.key == key).unwrap();
                    same_geometry(current.geometry.unwrap(), original.unwrap());
                }
                result
                    .manufacturing_mesh_with_receipt(Default::default())
                    .unwrap();
                let directory = tempfile::tempdir().unwrap();
                let path = directory.path().join("result.step");
                result.write_step(&path).unwrap();
                let read = Solid::import_rust_step(
                    &std::fs::read_to_string(&path).unwrap(),
                    EntityId::new(),
                    policy,
                )
                .unwrap();
                compare(result, &read);
                if OCCT_LINKED {
                    compare(result, &Solid::read_step(&path, policy).unwrap());
                }
            }
            if OCCT_LINKED {
                let directory = tempfile::tempdir().unwrap();
                let a = directory.path().join("source.step");
                let b = directory.path().join("tool.step");
                source.write_step(&a).unwrap();
                tool.write_step(&b).unwrap();
                let a = Solid::read_step(&a, policy).unwrap();
                let b = Solid::read_step(&b, policy).unwrap();
                compare(
                    &cut,
                    &a.boolean(&b, BooleanOperation::Subtract, policy)
                        .unwrap()
                        .unwrap(),
                );
                compare(
                    &join,
                    &a.boolean(&b, BooleanOperation::Union, policy)
                        .unwrap()
                        .unwrap(),
                );
                let native = a
                    .boolean(&b, BooleanOperation::Intersect, policy)
                    .unwrap()
                    .unwrap();
                near(common.0, native.volume_mm3().unwrap());
            }
            assert_eq!(source.topology_snapshot(owner, 0).unwrap(), before);
            assert_eq!(tool.topology_snapshot(owner, 0).unwrap(), tool_before);
        }
    }
}

fn lines(solid: &Solid) -> Vec<[[f64; 3]; 2]> {
    solid
        .topology_snapshot(EntityId::new(), 0)
        .unwrap()
        .elements
        .into_iter()
        .filter_map(|element| match element.geometry {
            Some(AnalyticGeometry::Line { a_mm, b_mm }) => Some([a_mm, b_mm]),
            _ => None,
        })
        .collect()
}
fn has_line(lines: &[[[f64; 3]; 2]], a: [f64; 3], b: [f64; 3]) -> bool {
    let same = |p: [f64; 3], q: [f64; 3]| (0..3).all(|i| (p[i] - q[i]).abs() < 1e-9);
    lines
        .iter()
        .any(|&[p, q]| (same(p, a) && same(q, b)) || (same(p, b) && same(q, a)))
}

#[test]
fn genuine_splits_and_authored_collinear_vertices_stay_explicit() {
    let policy = GeometryPolicy::default();
    let box_tool = |origin: [f64; 3], size: [[i64; 2]; 3]| {
        Solid::rust_box(
            EntityId::new().uuid_bits(),
            size.map(|[n, d]| Rational::new(n, d).unwrap()),
            PlaneFrame {
                origin,
                ..Default::default()
            },
            policy,
        )
        .unwrap()
    };
    let import = |geometry: &Union, owner| {
        Solid::import_rust_step(&placed_text(geometry, PlaneFrame::default()), owner, policy)
            .unwrap()
    };
    let [a, c] = [[0., 0., 0.], [0., 4., 0.]];

    // A notch through the middle of A-C removes part of that edge. The old
    // whole-edge reference must name both survivors, never pick one.
    let owner = EntityId::new();
    let source = import(&tetra(), owner);
    let before = source.topology_snapshot(owner, 0).unwrap();
    let edge = before
        .elements
        .iter()
        .find(|e| {
            matches!(e.geometry, Some(AnalyticGeometry::Line { a_mm, b_mm })
                if has_line(&[[a_mm, b_mm]], a, c))
        })
        .unwrap();
    let reference = GeometryIntent::capture(&before, edge.key).unwrap();
    let notched = source
        .rust_subtract_planar_tools(
            &[&box_tool([-1., 1., -1.], [[3, 2], [1, 1], [3, 2]])],
            EntityId::new(),
            policy,
        )
        .unwrap()
        .unwrap();
    near(notched.volume_mm3().unwrap(), 32. / 3. - 0.25);
    let output = notched.topology_snapshot(owner, 1).unwrap();
    let Resolution::Ambiguous { candidates } = resolve(&reference, &output, 1) else {
        panic!("a genuinely split edge must stay ambiguous");
    };
    assert_eq!(candidates.len(), 2);
    let segments = lines(&notched);
    assert!(has_line(&segments, a, [0., 1., 0.]));
    assert!(has_line(&segments, [0., 2., 0.], c));
    assert!(!has_line(&segments, a, c));

    // An authored midpoint on A-C is part of the imported boundary and
    // survives the same remote cut that no longer subdivides A-C itself.
    let m = [0., 2., 0.];
    let [b, d] = [[4., 0., 0.], [0., 0., 4.]];
    let authored = Union::new(
        &[ConvexPiece::new(
            vec![
                vec![a, m, c, b],
                vec![a, b, d],
                vec![a, d, c, m],
                vec![b, c, d],
            ],
            TOL,
        )
        .unwrap()],
        TOL,
    )
    .unwrap();
    let remote = box_tool([1., 0., 0.], [[4, 1], [1, 1], [1, 1]]);
    for (source, authored) in [
        (import(&authored, EntityId::new()), true),
        (import(&tetra(), owner), false),
    ] {
        let cut = source
            .rust_subtract_planar_tools(&[&remote], EntityId::new(), policy)
            .unwrap()
            .unwrap();
        let join = Solid::union_planar_pieces(&[&source, &remote], EntityId::new(), policy)
            .unwrap()
            .unwrap();
        for result in [&cut, &join] {
            let segments = lines(result);
            let midpoint = source
                .topology_snapshot(EntityId::new(), 0)
                .unwrap()
                .elements
                .iter()
                .any(|e| {
                    matches!(e.geometry, Some(AnalyticGeometry::Point { position_mm })
                        if position_mm == m)
                });
            assert_eq!(midpoint, authored);
            assert_eq!(has_line(&segments, a, m), midpoint);
            assert_eq!(has_line(&segments, m, c), midpoint);
            assert_eq!(has_line(&segments, a, c), !midpoint);
            assert!(!segments.iter().flatten().any(|&p| p == [0., 1., 0.]));
        }
    }
}
