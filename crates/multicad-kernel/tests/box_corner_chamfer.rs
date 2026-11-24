use multicad_core::{document::EntityId, references::*};
use multicad_kernel::*;

fn selected(stock: &Solid, owner: EntityId, feature: EntityId) -> Vec<TreatmentEdge> {
    stock
        .topology_snapshot(owner, 0)
        .unwrap()
        .elements
        .into_iter()
        .filter(|e| e.key.kind == TopologyKind::Edge)
        .map(|e| TreatmentEdge {
            key: e.key,
            origin: Origin {
                source: feature,
                role: format!("bevel/{}", e.key.ordinal),
            },
        })
        .collect()
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-6, "{a} != {b}");
}

#[test]
fn all_twelve_edges_have_true_planar_corner_closures_and_exact_material() {
    let policy = GeometryPolicy::default();
    let owner = EntityId::new();
    let feature = EntityId::new();
    let stock = Solid::box_named([10.; 3], owner, policy).unwrap();
    let before = stock.topology_snapshot(owner, 0).unwrap();
    let edges = selected(&stock, owner, feature);
    assert_eq!(edges.len(), 12);
    let output = stock
        .treat_edges_named(
            &edges,
            EdgeTreatment::Chamfer { distance_mm: 1. },
            feature,
            policy,
        )
        .unwrap();
    near(output.volume_mm3().unwrap(), 945. + 1. / 3.);
    let after = output.topology_snapshot(feature, 1).unwrap();
    assert_eq!(output.faces().unwrap().len(), 26);
    for edge in &edges {
        let faces: Vec<_> = after
            .elements
            .iter()
            .filter(|e| e.key.kind == TopologyKind::Face && e.origins.contains(&edge.origin))
            .collect();
        assert_eq!(faces.len(), 1, "{:?}", edge.origin);
    }
    let mesh = output.manufacturing_mesh(MeshOptions::default()).unwrap();
    let proof = rust_kernel::manufacturing::certify(&rust_kernel::display::DisplayMesh {
        vertices: mesh.vertices,
        indices: mesh.indices,
        face_indices: mesh.face_indices,
    })
    .unwrap();
    near(proof.signed_volume_mm3, 945. + 1. / 3.);
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("corners.step");
    output.write_step(&path).unwrap();
    let read = Solid::read_step(&path, policy).unwrap();
    near(read.volume_mm3().unwrap(), 945. + 1. / 3.);
    for (origin, occupied) in [([0.49; 3], false), ([0.7, 0.8, 0.8], true), ([4.; 3], true)] {
        let probe = Solid::box_named([0.01; 3], EntityId::new(), policy)
            .unwrap()
            .transformed(RigidTransform {
                translation_mm: origin,
                ..Default::default()
            })
            .unwrap();
        let overlap = output.rust_occupied_overlap_volume(&probe, policy).unwrap();
        let expected = if occupied { 1e-6 } else { 0. };
        assert!(
            (overlap - expected).abs() < 1e-11,
            "{overlap} != {expected}"
        );
    }
    assert_eq!(stock.topology_snapshot(owner, 0).unwrap(), before);
}

#[test]
fn subsets_pose_authored_cap_frames_cancel_and_refuse_atomically() {
    use std::cell::Cell;
    let p = GeometryPolicy::default();
    let owner = EntityId::new();
    let feature = EntityId::new();
    let stock = Solid::box_named([10., 16., 20.], owner, p)
        .unwrap()
        .with_feature_cap_roles(owner, "base.start", "base.end")
        .unwrap();
    let before = stock.topology_snapshot(owner, 0).unwrap();
    let all = selected(&stock, owner, feature);
    for count in 2..=12 {
        let edges = &all[..count];
        let d = 1.;
        let mut length = 0.;
        let mut ends = std::collections::BTreeMap::new();
        for edge in edges {
            let Some(AnalyticGeometry::Line { a_mm: a, b_mm: b }) = before
                .elements
                .iter()
                .find(|e| e.key == edge.key)
                .unwrap()
                .geometry
            else {
                panic!()
            };
            length += (0..3).map(|i| (a[i] - b[i]).abs()).sum::<f64>();
            for point in [a, b] {
                *ends.entry(point.map(f64::to_bits)).or_insert(0) += 1;
            }
        }
        let correction: f64 = ends
            .values()
            .map(|n| match n {
                2 => d * d * d / 3.,
                3 => 2. * d * d * d / 3.,
                _ => 0.,
            })
            .sum();
        let result = stock
            .treat_edges_named(edges, EdgeTreatment::Chamfer { distance_mm: d }, feature, p)
            .unwrap();
        near(
            result.volume_mm3().unwrap(),
            3200. - 0.5 * d * d * length + correction,
        );
        assert_eq!(
            result.faces().unwrap().len(),
            6 + count + ends.values().filter(|&&n| n == 3).count()
        );
        let mut reverse = edges.to_vec();
        reverse.reverse();
        let reordered = stock
            .treat_edges_named(
                &reverse,
                EdgeTreatment::Chamfer { distance_mm: d },
                feature,
                p,
            )
            .unwrap();
        assert_eq!(
            result.topology_snapshot(feature, 1).unwrap(),
            reordered.topology_snapshot(feature, 1).unwrap()
        );
    }
    let pose = RigidTransform {
        rotation: [[0.8, -0.36, 0.48], [0.6, 0.48, -0.64], [0., 0.8, 0.6]],
        translation_mm: [80., 60., 20.],
    };
    let moved = stock.transformed(pose).unwrap();
    let moved_before = moved.topology_snapshot(owner, 0).unwrap();
    let result = moved
        .treat_edges_named(&all, EdgeTreatment::Chamfer { distance_mm: 1. }, feature, p)
        .unwrap();
    let after = result.topology_snapshot(feature, 1).unwrap();
    for face in moved_before
        .elements
        .iter()
        .filter(|e| e.key.kind == TopologyKind::Face)
    {
        let mut intent = GeometryIntent::capture(&moved_before, face.key).unwrap();
        if let Some(anchor) = face.frames.first() {
            intent.attachment = Some(Attachment {
                outward: false,
                anchor: anchor.origin.clone(),
                offset_mm: [0.; 3],
                angle_rad: 0.,
            });
        }
        let Resolution::Resolved { key, frame, .. } = resolve(&intent, &after, 1) else {
            panic!("retained stock face needs no repair")
        };
        if let Some(frame) = frame {
            let old = face
                .frames
                .iter()
                .find(|f| f.origin == intent.attachment.as_ref().unwrap().anchor)
                .unwrap()
                .frame;
            for (a, b) in frame
                .origin_mm
                .into_iter()
                .chain(frame.u)
                .chain(frame.v)
                .zip(old.origin_mm.into_iter().chain(old.u).chain(old.v))
            {
                assert!((a - b).abs() < 1e-10);
            }
            let actual = result
                .faces()
                .unwrap()
                .into_iter()
                .find(|f| f.index == key.ordinal)
                .unwrap();
            let (o, n) = actual.plane.unwrap();
            assert!(
                (0..3)
                    .map(|i| (frame.origin_mm[i] - o[i]) * n[i])
                    .sum::<f64>()
                    .abs()
                    < p.linear_tolerance_mm
            );
        }
    }
    let calls = Cell::new(0);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    stock
        .treat_edges_named_with_interrupt(
            &all,
            EdgeTreatment::Chamfer { distance_mm: 1. },
            feature,
            p,
            Some(&count),
        )
        .unwrap();
    assert!(calls.get() > 20);
    for at in [1, calls.get() / 2, calls.get()] {
        let seen = Cell::new(0);
        let cancel = || {
            seen.set(seen.get() + 1);
            seen.get() == at
        };
        assert_eq!(
            stock
                .treat_edges_named_with_interrupt(
                    &all,
                    EdgeTreatment::Chamfer { distance_mm: 1. },
                    feature,
                    p,
                    Some(&cancel)
                )
                .err()
                .unwrap()
                .kind,
            GeometryErrorKind::Cancelled
        );
        assert_eq!(stock.topology_snapshot(owner, 0).unwrap(), before);
    }
    for d in [5., 6., f64::NAN] {
        assert!(
            stock
                .treat_edges_named(&all, EdgeTreatment::Chamfer { distance_mm: d }, feature, p)
                .is_err()
        );
    }
    assert!(
        stock
            .treat_edges_named(
                &[all[0].clone(), all[0].clone()],
                EdgeTreatment::Chamfer { distance_mm: 1. },
                feature,
                p
            )
            .is_err()
    );
    stock
        .treat_edges_named(&all, EdgeTreatment::Chamfer { distance_mm: 1. }, feature, p)
        .unwrap();
    assert_eq!(stock.topology_snapshot(owner, 0).unwrap(), before);
    assert_eq!(moved.topology_snapshot(owner, 0).unwrap(), moved_before);
}

#[test]
fn interruption_does_not_narrow_the_qualified_parallel_edge_envelope() {
    let p = GeometryPolicy::default();
    let owner = EntityId::new();
    let feature = EntityId::new();
    let stock = Solid::box_named([10., 20., 30.], owner, p).unwrap();
    let topology = stock.topology_snapshot(owner, 0).unwrap();
    let edges:Vec<_>=selected(&stock,owner,feature).into_iter().filter(|edge|matches!(topology.elements.iter().find(|e|e.key==edge.key).unwrap().geometry,Some(AnalyticGeometry::Line{a_mm:a,b_mm:b}) if a[2]!=b[2] && ((a[0]==0.&&a[1]==0.)||(a[0]==10.&&a[1]==20.)))).collect();
    assert_eq!(edges.len(), 2);
    let ordinary = stock
        .treat_edges_named(
            &edges,
            EdgeTreatment::Chamfer { distance_mm: 6. },
            feature,
            p,
        )
        .unwrap();
    let controlled = stock
        .treat_edges_named_with_interrupt(
            &edges,
            EdgeTreatment::Chamfer { distance_mm: 6. },
            feature,
            p,
            Some(&|| false),
        )
        .unwrap();
    near(controlled.volume_mm3().unwrap(), 4920.);
    assert_eq!(
        ordinary.topology_snapshot(feature, 1).unwrap(),
        controlled.topology_snapshot(feature, 1).unwrap()
    );
}

#[test]
fn five_authored_edges_with_nuls_keep_both_corner_identities_distinct() {
    let p = GeometryPolicy::default();
    let owner = EntityId::new();
    let feature = EntityId::new();
    let stock = Solid::box_named([10.; 3], owner, p).unwrap();
    let before = stock.topology_snapshot(owner, 0).unwrap();
    let edges: Vec<_> = selected(&stock, owner, feature)
        .into_iter()
        .filter_map(|mut edge| {
            let Some(AnalyticGeometry::Line { a_mm: a, b_mm: b }) = before
                .elements
                .iter()
                .find(|e| e.key == edge.key)
                .unwrap()
                .geometry
            else {
                panic!()
            };
            let role = if a[0] != b[0] && a[1] == 0. && a[2] == 0. {
                "d"
            } else if a[1] != b[1] && a[0] == 0. && a[2] == 0. {
                "a"
            } else if a[2] != b[2] && a[0] == 0. && a[1] == 0. {
                "b\0c"
            } else if a[1] != b[1] && a[0] == 10. && a[2] == 0. {
                "a\0b"
            } else if a[2] != b[2] && a[0] == 10. && a[1] == 0. {
                "c"
            } else {
                return None;
            };
            edge.origin.role = role.into();
            Some(edge)
        })
        .collect();
    assert_eq!(edges.len(), 5);
    let result = stock
        .treat_edges_named(
            &edges,
            EdgeTreatment::Chamfer { distance_mm: 1. },
            feature,
            p,
        )
        .unwrap();
    let topology = result.topology_snapshot(feature, 1).unwrap();
    near(result.volume_mm3().unwrap(), 976. + 1. / 3.);
    let corners: Vec<_> = topology
        .elements
        .iter()
        .filter(|e| {
            e.key.kind == TopologyKind::Face
                && e.origins
                    .iter()
                    .any(|o| o.role.starts_with("chamfer.corner/"))
        })
        .collect();
    assert_eq!(corners.len(), 2);
    assert_ne!(corners[0].origins, corners[1].origins);
    let captures: Vec<_> = corners
        .iter()
        .map(|e| GeometryIntent::capture(&topology, e.key).unwrap())
        .collect();
    let mut reverse = edges.clone();
    reverse.reverse();
    let edited = stock
        .treat_edges_named(
            &reverse,
            EdgeTreatment::Chamfer { distance_mm: 1.5 },
            feature,
            p,
        )
        .unwrap();
    let after = edited.topology_snapshot(feature, 2).unwrap();
    let resolved: Vec<_> = captures
        .iter()
        .map(|intent| {
            let Resolution::Resolved { key, .. } = resolve(intent, &after, 2) else {
                panic!("corner identity needs no repair")
            };
            key
        })
        .collect();
    assert_ne!(resolved[0], resolved[1]);
    assert_eq!(stock.topology_snapshot(owner, 0).unwrap(), before);
}

#[test]
fn trihedral_cut_retains_far_sharp_geometry_and_deletes_selected_sharp_geometry() {
    let p = GeometryPolicy::default();
    let owner = EntityId::new();
    let feature = EntityId::new();
    let stock = Solid::box_named([10.; 3], owner, p).unwrap();
    let before = stock.topology_snapshot(owner, 0).unwrap();
    let edges:Vec<_>=selected(&stock,owner,feature).into_iter().filter(|edge|matches!(before.elements.iter().find(|e|e.key==edge.key).unwrap().geometry,Some(AnalyticGeometry::Line{a_mm:a,b_mm:b}) if (0..3).filter(|&i|a[i]==0.&&b[i]==0.).count()==2)).collect();
    assert_eq!(edges.len(), 3);
    let far: Vec<_> = before
        .elements
        .iter()
        .filter(|e| match e.geometry {
            Some(AnalyticGeometry::Point { position_mm }) => position_mm == [10.; 3],
            Some(AnalyticGeometry::Line { a_mm: a, b_mm: b }) => {
                e.key.kind == TopologyKind::Edge
                    && (0..3).filter(|&i| a[i] == 10. && b[i] == 10.).count() == 2
            }
            _ => false,
        })
        .collect();
    assert_eq!(far.len(), 4);
    let removed: Vec<_> = before
        .elements
        .iter()
        .filter(|e| {
            matches!(e.geometry,Some(AnalyticGeometry::Point{position_mm}) if position_mm==[0.;3])
                || edges.iter().any(|edge| e.key == edge.key)
        })
        .collect();
    assert_eq!(removed.len(), 4);
    for controlled in [false, true] {
        let never = || false;
        let interrupt = controlled.then_some(&never as &dyn Fn() -> bool);
        let output = stock
            .treat_edges_named_with_interrupt(
                &edges,
                EdgeTreatment::Chamfer { distance_mm: 1. },
                feature,
                p,
                interrupt,
            )
            .unwrap();
        let after = output.topology_snapshot(feature, 1).unwrap();
        for old in &far {
            let intent = GeometryIntent::capture(&before, old.key).unwrap();
            let Resolution::Resolved { geometry, .. } = resolve(&intent, &after, 1) else {
                panic!("untouched sharp geometry needs no repair")
            };
            match (old.geometry, geometry) {
                (
                    Some(AnalyticGeometry::Point { position_mm: a }),
                    Some(AnalyticGeometry::Point { position_mm: b }),
                ) => assert_eq!(a, b),
                (
                    Some(AnalyticGeometry::Line { a_mm: a, b_mm: b }),
                    Some(AnalyticGeometry::Line { a_mm: c, b_mm: d }),
                ) => assert!((a == c && b == d) || (a == d && b == c)),
                _ => panic!("retained finite geometry must agree"),
            }
        }
        for old in &removed {
            let intent = GeometryIntent::capture(&before, old.key).unwrap();
            assert!(matches!(resolve(&intent, &after, 1), Resolution::Missing));
        }
    }
    assert_eq!(stock.topology_snapshot(owner, 0).unwrap(), before);
}

#[test]
fn adding_third_incident_edge_keeps_bevel_identity_distinct_from_corner() {
    let p = GeometryPolicy::default();
    let owner = EntityId::new();
    let feature = EntityId::new();
    let stock = Solid::box_named([10.; 3], owner, p).unwrap();
    let snapshot = stock.topology_snapshot(owner, 0).unwrap();
    let edges: Vec<_> = selected(&stock, owner, feature)
        .into_iter()
        .filter(|e| {
            let Some(AnalyticGeometry::Line { a_mm: a, b_mm: b }) = snapshot
                .elements
                .iter()
                .find(|f| f.key == e.key)
                .unwrap()
                .geometry
            else {
                panic!()
            };
            (0..3).all(|i| a[i] == 0. || b[i] == 0.)
                && (0..3).filter(|&i| a[i] == 0. && b[i] == 0.).count() == 2
        })
        .collect();
    assert_eq!(edges.len(), 3);
    let first = stock
        .treat_edges_named(
            &edges[..2],
            EdgeTreatment::Chamfer { distance_mm: 1. },
            feature,
            p,
        )
        .unwrap();
    let before = first.topology_snapshot(feature, 1).unwrap();
    let captures: Vec<_> = edges[..2]
        .iter()
        .map(|edge| {
            let face = before
                .elements
                .iter()
                .find(|e| e.key.kind == TopologyKind::Face && e.origins.contains(&edge.origin))
                .unwrap();
            GeometryIntent::capture(&before, face.key).unwrap()
        })
        .collect();
    let second = stock
        .treat_edges_named(
            &edges,
            EdgeTreatment::Chamfer { distance_mm: 1. },
            feature,
            p,
        )
        .unwrap();
    let after = second.topology_snapshot(feature, 2).unwrap();
    for intent in &captures {
        assert!(
            matches!(resolve(intent, &after, 2), Resolution::Resolved { .. }),
            "{:?}",
            resolve(intent, &after, 2)
        );
    }
    let corner: Vec<_> = after
        .elements
        .iter()
        .filter(|e| {
            e.key.kind == TopologyKind::Face
                && e.origins
                    .iter()
                    .any(|o| o.role.starts_with("chamfer.corner/"))
        })
        .collect();
    assert_eq!(corner.len(), 1);
    assert!(
        edges
            .iter()
            .all(|edge| !corner[0].origins.contains(&edge.origin))
    );
    assert!(
        edges
            .iter()
            .all(|edge| corner[0].boundary.contains(&edge.origin))
    );
}
