#![cfg(not(feature = "occt"))]
//! Oblique cylindrical subtraction from box stock through the normal Boolean.

use multicad_core::{
    document::EntityId,
    references::{AnalyticGeometry, Origin, TopologyElement, TopologyKind, TopologySnapshot},
};
use multicad_kernel::{
    BRepSelection, BooleanOperation, EdgeTreatment, GeometryErrorKind, GeometryPolicy, MeshOptions,
    PlaneFrame, RigidTransform, Solid, TreatmentEdge,
    rust_kernel::{Circle, Rational, cylinder::AxialCylinder, polygon::PolygonPrism},
};
use std::{
    cell::Cell,
    collections::BTreeSet,
    f64::consts::{PI, TAU},
    time::Instant,
};

type V = [f64; 3];

const OPTIONS: MeshOptions = MeshOptions {
    linear_deflection_mm: 0.002,
    angular_deflection_rad: 0.2,
    maximum_triangles: 2_000_000,
};

fn dot(a: V, b: V) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: V, b: V) -> V {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn normalized(v: V) -> V {
    let length = dot(v, v).sqrt();
    v.map(|x| x / length)
}

/// A tool frame whose normal is the requested travel direction.
fn frame(origin: V, travel: V) -> PlaneFrame {
    let w = normalized(travel);
    let seed = if w[2].abs() < 0.9 {
        [0., 0., 1.]
    } else {
        [1., 0., 0.]
    };
    let u = normalized(cross(seed, w));
    PlaneFrame {
        origin,
        u,
        v: cross(w, u),
    }
}

fn policy() -> GeometryPolicy {
    GeometryPolicy::default()
}

fn stock() -> Solid {
    Solid::box_solid([10., 8., 6.], policy()).unwrap()
}

fn cylinder(origin: V, travel: V, radius: f64, height: f64) -> Solid {
    Solid::cylinder(frame(origin, travel), radius, height, policy()).unwrap()
}

fn subtract(stock: &Solid, tool: &Solid) -> Solid {
    stock
        .boolean(tool, BooleanOperation::Subtract, policy())
        .unwrap()
        .unwrap()
}

fn audit_fixture() -> (Solid, Solid) {
    let body = Solid::box_solid([10.; 3], policy()).unwrap();
    let s = 2f64.sqrt().recip();
    let t = 6f64.sqrt().recip();
    let n = 3f64.sqrt().recip();
    let bore = Solid::cylinder(
        PlaneFrame {
            origin: [1. / 3. - n; 3],
            u: [s, -s, 0.],
            v: [t, t, -2. * t],
        },
        0.1,
        20.,
        policy(),
    )
    .unwrap();
    (body, bore)
}

fn elements(snapshot: &TopologySnapshot, kind: TopologyKind) -> Vec<&TopologyElement> {
    let found = snapshot.elements.iter().filter(|e| e.key.kind == kind);
    found.collect()
}

/// Certified manufacturing mesh volume, checked against the analytic solid,
/// its bounds and its face labels.
fn certified_volume(solid: &Solid) -> f64 {
    let receipt = solid.manufacturing_mesh_with_receipt(OPTIONS).unwrap();
    let mesh = &receipt.mesh;
    assert_eq!(receipt.maximum_displacement_mm, 0.);
    assert_eq!(receipt.edge_pairs as usize * 2, mesh.indices.len());
    let faces = solid.faces().unwrap();
    let vertices = mesh.vertices.as_chunks::<3>().0;
    let origin = vertices[0];
    let mut volume = 0.;
    let mut areas = vec![0.; faces.len()];
    for (triangle, &face) in mesh
        .indices
        .as_chunks::<3>()
        .0
        .iter()
        .zip(&mesh.face_indices)
    {
        let [a, b, c] = triangle.map(|i| {
            let p = vertices[i as usize];
            [p[0] - origin[0], p[1] - origin[1], p[2] - origin[2]]
        });
        volume += dot(a, cross(b, c)) / 6.;
        let normal = cross(
            [b[0] - a[0], b[1] - a[1], b[2] - a[2]],
            [c[0] - a[0], c[1] - a[1], c[2] - a[2]],
        );
        areas[face as usize] += dot(normal, normal).sqrt() / 2.;
        if let Some((_, outward)) = faces[face as usize].plane {
            assert!(dot(normal, outward) > 0.);
        }
    }
    let curved: f64 = faces
        .iter()
        .filter(|f| f.plane.is_none())
        .map(|f| f.area_mm2)
        .sum();
    let exact = solid.volume_mm3().unwrap();
    assert!(
        (volume - exact).abs() <= OPTIONS.linear_deflection_mm * curved + 1e-7,
        "mesh {volume} analytic {exact}"
    );
    for (face, area) in faces.iter().zip(areas) {
        assert!((face.area_mm2 - area).abs() < 3e-3 * face.area_mm2 + 1e-9);
    }
    let bounds = solid.bounds().unwrap();
    for axis in 0..3 {
        let values = vertices.iter().map(|p| p[axis]);
        let low = values.clone().fold(f64::INFINITY, f64::min);
        let high = values.fold(f64::NEG_INFINITY, f64::max);
        assert!(bounds.min[axis] <= low + 1e-9 && low - bounds.min[axis] < 0.01);
        assert!(bounds.max[axis] >= high - 1e-9 && bounds.max[axis] - high < 0.01);
    }
    let display = solid.mesh(OPTIONS).unwrap();
    assert_eq!(display.indices.len(), mesh.indices.len());
    volume
}

/// Independent midpoint quadrature of the tool volume inside a box at the
/// world origin.
fn removed_volume(size: V, tool: PlaneFrame, radius: f64, height: f64) -> f64 {
    let w = tool.normal();
    let (rings, turns) = (1000, 2000);
    let mut volume = 0.;
    for i in 0..rings {
        let r = radius * (i as f64 + 0.5) / rings as f64;
        for j in 0..turns {
            let (sin, cos) = (TAU * (j as f64 + 0.5) / turns as f64).sin_cos();
            let start: V =
                std::array::from_fn(|k| tool.origin[k] + r * (cos * tool.u[k] + sin * tool.v[k]));
            let mut range = [0., height];
            for k in 0..3 {
                let ends = [-start[k] / w[k], (size[k] - start[k]) / w[k]];
                range[0] = range[0].max(ends[0].min(ends[1]));
                range[1] = range[1].min(ends[0].max(ends[1]));
            }
            volume += (range[1] - range[0]).max(0.) * r;
        }
    }
    volume * radius / rings as f64 * TAU / turns as f64
}

#[test]
fn audit_fixture_keeps_the_cylinder_support_and_finite_trims() {
    let (body, bore) = audit_fixture();
    let started = Instant::now();
    let hollow = subtract(&body, &bore);
    let construction = started.elapsed();

    let snapshot = hollow.topology_snapshot(EntityId::new(), 3).unwrap();
    snapshot.validate().unwrap();
    let faces = elements(&snapshot, TopologyKind::Face);
    let edges = elements(&snapshot, TopologyKind::Edge);
    assert_eq!(
        [
            faces.len(),
            edges.len(),
            elements(&snapshot, TopologyKind::Vertex).len()
        ],
        [7, 18, 12]
    );
    // Six trimmed stock planes and one bore face without a planar support.
    let planar = faces
        .iter()
        .filter(|f| matches!(f.geometry, Some(AnalyticGeometry::Plane { .. })));
    assert_eq!(planar.count(), 6);
    assert_eq!(faces.iter().filter(|f| f.geometry.is_none()).count(), 1);
    // Every trim between the bore and a stock plane is a finite arc of the
    // actual section ellipse of the radius 0.1 cylinder.
    let mut arcs = 0;
    for edge in &edges {
        match edge.geometry {
            Some(AnalyticGeometry::EllipticArc {
                frame,
                semiaxes_mm,
                sweep_rad,
                ..
            }) => {
                arcs += 1;
                assert!((semiaxes_mm[0] - 0.1 * 3f64.sqrt()).abs() < 1e-12);
                assert!((semiaxes_mm[1] - 0.1).abs() < 1e-12);
                assert!(sweep_rad > 2. && sweep_rad < 2.2);
                // The center lies on the bore axis x = y = z.
                let c = frame.origin_mm;
                assert!((c[0] - c[1]).abs().max((c[1] - c[2]).abs()) < 0.15);
                let diagonal = [3f64.sqrt().recip(); 3];
                let radial = cross(diagonal, cross(c, diagonal));
                assert!(dot(radial, radial).sqrt() < 1e-12);
            }
            Some(AnalyticGeometry::Line { .. }) => {}
            other => panic!("unexpected edge support {other:?}"),
        }
    }
    assert_eq!(arcs, 6);
    // Anonymous constructors publish no persistent names.
    assert!(snapshot.elements.iter().all(|e| e.origins.is_empty()));

    let frame = PlaneFrame {
        origin: [1. / 3. - 3f64.sqrt().recip(); 3],
        u: [2f64.sqrt().recip(), -(2f64.sqrt().recip()), 0.],
        v: [
            6f64.sqrt().recip(),
            6f64.sqrt().recip(),
            -2. * 6f64.sqrt().recip(),
        ],
    };
    let removed = removed_volume([10.; 3], frame, 0.1, 20.);
    let volume = hollow.volume_mm3().unwrap();
    assert!(
        (1000. - volume - removed).abs() < 4e-5 * removed,
        "{volume} {removed}"
    );
    // The bore is shorter than the full diagonal by its two corner trims.
    assert!(removed < PI * 0.01 * 10. * 3f64.sqrt() && removed > 0.53);
    let bounds = hollow.bounds().unwrap();
    for axis in 0..3 {
        assert!(bounds.min[axis].abs() < 1e-12 && (bounds.max[axis] - 10.).abs() < 1e-12);
    }
    assert_eq!(hollow.solid_count().unwrap(), 1);
    let started = Instant::now();
    let meshed = certified_volume(&hollow);
    println!(
        "oblique audit fixture: construction {construction:?}, certified mesh checks {:?}, \
         analytic {volume:.9}, mesh {meshed:.9}",
        started.elapsed()
    );
    // The sources are untouched.
    assert_eq!(body.volume_mm3().unwrap(), 1000.);
    assert_eq!(body.faces().unwrap().len(), 6);
    assert_eq!(bore.faces().unwrap().len(), 3);
}

/// Unit normal and offset of a plane.
type Plane = (V, f64);

fn face_origins(snapshot: &TopologySnapshot) -> Vec<(Option<Plane>, BTreeSet<Origin>)> {
    elements(snapshot, TopologyKind::Face)
        .into_iter()
        .map(|face| {
            let plane = match face.geometry {
                Some(AnalyticGeometry::Plane { frame }) => {
                    let normal = frame.normal();
                    Some((normal, dot(normal, frame.origin_mm)))
                }
                _ => None,
            };
            (plane, face.origins.clone())
        })
        .collect()
}

fn named_cylinder(feature: EntityId, frame: PlaneFrame, radius: (i64, i64), height: i64) -> Solid {
    let zero = Rational::new(0, 1).unwrap();
    let circle = Circle::new([zero; 2], Rational::new(radius.0, radius.1).unwrap()).unwrap();
    let cylinder = AxialCylinder::new(
        feature.uuid_bits(),
        circle,
        zero,
        Rational::new(height, 1).unwrap(),
    )
    .unwrap();
    let origin = Origin {
        source: feature,
        role: "profile.circle".into(),
    };
    Solid::from_axial_cylinder_sweep(cylinder, frame, [0.; 2], false, origin, policy()).unwrap()
}

fn same_plane(a: Plane, b: Plane) -> bool {
    (0..3).all(|i| (a.0[i] - b.0[i]).abs() < 1e-9) && (a.1 - b.1).abs() < 1e-9
}

#[test]
fn surviving_stock_faces_keep_their_origins_and_the_bore_names_the_tool() {
    let (box_id, tool_id) = (EntityId::new(), EntityId::new());
    let sizes = [10, 8, 6].map(|n| Rational::new(n, 1).unwrap());
    let zero = Rational::new(0, 1).unwrap();
    let ring = vec![
        [zero, zero],
        [sizes[0], zero],
        [sizes[0], sizes[1]],
        [zero, sizes[1]],
    ];
    let curves: Vec<Origin> = ["y.min", "x.max", "y.max", "x.min"]
        .into_iter()
        .map(|role| Origin {
            source: box_id,
            role: role.into(),
        })
        .collect();
    let polygon = |reverse| {
        Solid::from_polygon_prism_sweep(
            PolygonPrism::new(vec![ring.clone()], sizes[2]).unwrap(),
            PlaneFrame {
                origin: [0., 0., if reverse { 6. } else { 0. }],
                ..Default::default()
            },
            reverse,
            box_id,
            &curves,
            &[false; 4],
            policy(),
        )
        .unwrap()
    };
    let stocks = [
        Solid::rust_box(box_id.uuid_bits(), sizes, PlaneFrame::default(), policy()).unwrap(),
        Solid::box_named([10., 8., 6.], box_id, policy()).unwrap(),
        polygon(false),
        polygon(true),
    ];
    // A blind pocket entering the top face keeps its floor inside.
    let tool_frame = frame([5., 4., 8.], [0.2, -0.1, -1.]);
    let tool = named_cylinder(tool_id, tool_frame, (1, 1), 5);
    let owner = EntityId::new();
    let tool_faces = face_origins(&tool.topology_snapshot(owner, 0).unwrap());
    let removed = PI * (5. + 2. / tool_frame.normal()[2]);
    for stock in &stocks {
        let before = face_origins(&stock.topology_snapshot(owner, 0).unwrap());
        assert!(before.iter().all(|(_, origins)| !origins.is_empty()));
        let cut = subtract(stock, &tool);
        assert!((480. - cut.volume_mm3().unwrap() - removed).abs() < 1e-9);
        let snapshot = cut.topology_snapshot(owner, 0).unwrap();
        snapshot.validate().unwrap();
        let after = face_origins(&snapshot);
        assert_eq!(after.len(), 8);
        let mut stock_faces = 0;
        for (plane, origins) in &after {
            let source = plane.and_then(|plane| {
                let same = before.iter().find(|(p, _)| same_plane(p.unwrap(), plane));
                same.map(|(_, origins)| origins)
            });
            if let Some(expected) = source {
                // Each surviving stock plane keeps exactly its source names.
                stock_faces += 1;
                assert_eq!(origins, expected);
                assert!(origins.iter().all(|origin| origin.source == box_id));
            } else if plane.is_some() {
                // The pocket floor is the end cap of the tool.
                assert_eq!(origins, &tool_faces[1].1);
            } else {
                assert_eq!(origins, &tool_faces[2].1);
                assert!(origins.iter().any(|o| o.role == "profile.circle"));
            }
            assert!(!origins.is_empty());
        }
        assert_eq!(stock_faces, 6);
        // Edges and vertices are named by their incident faces.
        for element in &snapshot.elements {
            assert!(!element.origins.is_empty());
        }
        let mouth = elements(&snapshot, TopologyKind::Edge)
            .into_iter()
            .filter(|e| matches!(e.geometry, Some(AnalyticGeometry::EllipticArc { .. })));
        let floor = elements(&snapshot, TopologyKind::Edge)
            .into_iter()
            .filter(|e| matches!(e.geometry, Some(AnalyticGeometry::CircularArc { .. })));
        for rim in mouth.chain(floor) {
            assert!(rim.origins.iter().any(|o| o.source == tool_id));
            assert!(
                rim.origins
                    .iter()
                    .any(|o| o.source == box_id || o.role.starts_with("cap"))
            );
        }
        certified_volume(&cut);
        // Source frames stay attached to their supports: literally at the
        // cut placement and on the moved support after a rigid motion.
        let pose = RigidTransform {
            rotation: rotation([1., 2., 3.], 0.7),
            translation_mm: [4., -3., 2.],
        };
        let moved = cut.transformed(pose).unwrap();
        let shared = cut.share().unwrap();
        assert_eq!(
            face_origins(&shared.topology_snapshot(owner, 0).unwrap()),
            after
        );
        for (solid, source) in [(&cut, Some(stock)), (&moved, None)] {
            let snapshot = solid.topology_snapshot(owner, 0).unwrap();
            snapshot.validate().unwrap();
            let mut anchored = 0;
            for face in elements(&snapshot, TopologyKind::Face) {
                let Some(AnalyticGeometry::Plane { frame: support }) = face.geometry else {
                    assert!(face.frames.is_empty());
                    continue;
                };
                for anchor in &face.frames {
                    anchored += 1;
                    assert!(face.origins.contains(&anchor.origin));
                    // Each source frame stays in its support plane.
                    let normal = support.normal();
                    assert!(dot(anchor.frame.normal(), normal).abs() > 1. - 1e-9);
                    let offset: V =
                        std::array::from_fn(|i| anchor.frame.origin_mm[i] - support.origin_mm[i]);
                    assert!(dot(offset, normal).abs() < 1e-9);
                }
                if let Some(stock) = source {
                    // A surviving stock face keeps its literal source frames.
                    let before = stock.topology_snapshot(owner, 0).unwrap();
                    let same = elements(&before, TopologyKind::Face)
                        .into_iter()
                        .find(|f| f.origins == face.origins);
                    if let Some(same) = same {
                        assert_eq!(same.frames, face.frames);
                    }
                }
            }
            assert!(anchored >= 7, "{anchored}");
        }
        // Moving the cut moves its source frames exactly like moving the
        // sources first.
        let moved_stock = stock.transformed(pose).unwrap();
        let expected = moved_stock.topology_snapshot(owner, 0).unwrap();
        let actual = moved.topology_snapshot(owner, 0).unwrap();
        let mut compared = 0;
        for face in elements(&actual, TopologyKind::Face) {
            let same = elements(&expected, TopologyKind::Face)
                .into_iter()
                .find(|f| f.origins == face.origins);
            if let Some(same) = same {
                compared += 1;
                assert_eq!(same.frames.len(), face.frames.len());
                for (a, b) in same.frames.iter().zip(&face.frames) {
                    let close = |x: V, y: V| (0..3).all(|i| (x[i] - y[i]).abs() < 1e-9);
                    assert_eq!(a.origin, b.origin);
                    assert!(close(a.frame.origin_mm, b.frame.origin_mm));
                    assert!(close(a.frame.u, b.frame.u) && close(a.frame.v, b.frame.v));
                }
            }
        }
        assert_eq!(compared, 6);
        assert_eq!(moved.display_edges(0.01).unwrap().len(), 14);
        // A world-vertical drill is oblique to the same stock once the
        // stock itself is rotated.
        let tilted = stock.transformed(pose).unwrap();
        let center = tilted.bounds().unwrap();
        let center: V = std::array::from_fn(|i| (center.min[i] + center.max[i]) / 2.);
        let drill = Solid::cylinder(
            PlaneFrame {
                origin: [center[0], center[1], center[2] - 20.],
                ..Default::default()
            },
            1.,
            40.,
            policy(),
        )
        .unwrap();
        let drilled = subtract(&tilted, &drill);
        let through = drilled.topology_snapshot(owner, 0).unwrap();
        let bore = elements(&through, TopologyKind::Face)
            .into_iter()
            .find(|face| face.geometry.is_none())
            .unwrap();
        assert!(bore.origins.iter().all(|origin| origin.source != box_id));
        assert!(480. - drilled.volume_mm3().unwrap() > 6. * PI);
        certified_volume(&drilled);
        // The source stock is unchanged.
        assert_eq!(stock.volume_mm3().unwrap(), 480.);
        assert_eq!(
            face_origins(&stock.topology_snapshot(owner, 0).unwrap()),
            before
        );
    }
}

fn rotation(axis: V, angle: f64) -> [[f64; 3]; 3] {
    let [x, y, z] = normalized(axis);
    let (s, c) = angle.sin_cos();
    let t = 1. - c;
    [
        [t * x * x + c, t * x * y - s * z, t * x * z + s * y],
        [t * x * y + s * z, t * y * y + c, t * y * z - s * x],
        [t * x * z - s * y, t * y * z + s * x, t * z * z + c],
    ]
}

#[test]
fn ordinary_rotated_far_and_reversed_placements_agree() {
    let stock = stock();
    // A through hole across the two z faces.
    let through = cylinder([3., 4., -2.], [0.2, 0.1, 1.], 1., 20.);
    let w = normalized([0.2, 0.1, 1.]);
    let removed = PI * 6. / w[2];
    let ordinary = subtract(&stock, &through);
    assert!((480. - ordinary.volume_mm3().unwrap() - removed).abs() < 1e-9);
    assert_eq!(ordinary.faces().unwrap().len(), 7);
    certified_volume(&ordinary);

    // The same travel entered from the far end.
    let end = [3. + 20. * w[0], 4. + 20. * w[1], -2. + 20. * w[2]];
    let reversed = subtract(&stock, &cylinder(end, w.map(|c| -c), 1., 20.));
    assert!((reversed.volume_mm3().unwrap() - ordinary.volume_mm3().unwrap()).abs() < 1e-9);
    certified_volume(&reversed);

    let poses = [
        RigidTransform {
            rotation: rotation([1., 2., 3.], 0.7),
            translation_mm: [4., -3., 2.],
        },
        RigidTransform {
            rotation: rotation([-2., 1., 0.5], 2.1),
            translation_mm: [1e4, -2e4, 3e4],
        },
        RigidTransform {
            translation_mm: [1e4, 1e4, -1e4],
            ..Default::default()
        },
    ];
    for pose in poses {
        let placed_stock = stock.transformed(pose).unwrap();
        let placed_tool = through.transformed(pose).unwrap();
        let cut = subtract(&placed_stock, &placed_tool);
        let volume = cut.volume_mm3().unwrap();
        assert!(
            (volume - ordinary.volume_mm3().unwrap()).abs() < 1e-7,
            "{volume}"
        );
        assert_eq!(cut.faces().unwrap().len(), 7);
        let meshed = certified_volume(&cut);
        // Moving the finished cut gives the same solid as cutting in place.
        let moved = ordinary.transformed(pose).unwrap();
        assert_eq!(moved.volume_mm3().unwrap(), ordinary.volume_mm3().unwrap());
        let (a, b) = (moved.bounds().unwrap(), cut.bounds().unwrap());
        for axis in 0..3 {
            assert!((a.min[axis] - b.min[axis]).abs() < 1e-7);
            assert!((a.max[axis] - b.max[axis]).abs() < 1e-7);
        }
        assert!((certified_volume(&moved) - meshed).abs() < 1e-6);
        let moved_topology = moved.topology_snapshot(EntityId::new(), 0).unwrap();
        moved_topology.validate().unwrap();
        assert_eq!(moved_topology.elements.len(), 31);
        for (x, y) in moved.faces().unwrap().iter().zip(cut.faces().unwrap()) {
            assert!((x.area_mm2 - y.area_mm2).abs() < 1e-7);
            assert!((0..3).all(|i| (x.center_mm[i] - y.center_mm[i]).abs() < 1e-6));
        }
    }
}

#[test]
fn partial_penetration_generators_and_corner_exits_are_exact() {
    let stock = stock();
    // Finite partial penetration: the pocket floor is the tool end cap.
    let blind = subtract(&stock, &cylinder([5., 4., 8.], [0.2, -0.1, -1.], 1., 5.));
    let depth = 5. - 2. / normalized([0.2, -0.1, -1.])[2].abs();
    assert!((480. - blind.volume_mm3().unwrap() - PI * depth).abs() < 1e-9);
    let floor = blind.faces().unwrap();
    let floor = floor
        .iter()
        .find(|f| f.plane.is_some() && (f.area_mm2 - PI).abs() < 1e-9)
        .unwrap();
    // The floor faces back along the removed tool.
    assert!(dot(floor.plane.unwrap().1, [-0.2, 0.1, 1.]) > 0.9);
    certified_volume(&blind);

    // Entered from the far end, the same floor is the tool start cap.
    let end = normalized([0.2, -0.1, -1.]).map(|c| 5. * c);
    let end = [5. + end[0], 4. + end[1], 8. + end[2]];
    let reversed = subtract(&stock, &cylinder(end, [-0.2, 0.1, 1.], 1., 5.));
    assert!((reversed.volume_mm3().unwrap() - blind.volume_mm3().unwrap()).abs() < 1e-9);
    let floors = reversed.faces().unwrap();
    let floor = floors
        .iter()
        .find(|f| f.plane.is_some() && (f.area_mm2 - PI).abs() < 1e-9)
        .unwrap();
    assert!(dot(floor.plane.unwrap().1, [-0.2, 0.1, 1.]) > 0.9);
    certified_volume(&reversed);

    // An axis parallel to the z faces that opens the top face in a slot.
    let slot = subtract(&stock, &cylinder([-2., -1., 5.5], [1., 1., 0.], 1., 30.));
    let snapshot = slot.topology_snapshot(EntityId::new(), 0).unwrap();
    assert_eq!(elements(&snapshot, TopologyKind::Face).len(), 8);
    let top = |p: V| (p[2] - 6.).abs() < 1e-12;
    let generators = elements(&snapshot, TopologyKind::Edge)
        .into_iter()
        .filter(|edge| match edge.geometry {
            Some(AnalyticGeometry::Line { a_mm, b_mm }) => {
                let d = [b_mm[0] - a_mm[0], b_mm[1] - a_mm[1], b_mm[2] - a_mm[2]];
                top(a_mm) && top(b_mm) && (d[0] - d[1]).abs() < 1e-9 && d[0].abs() > 1.
            }
            _ => false,
        });
    assert_eq!(generators.count(), 2);
    certified_volume(&slot);
    // Travel within the chart noise of that face pair is the same cut.
    let noisy = subtract(&stock, &cylinder([-2., -1., 5.5], [1., 1., 1e-17], 1., 30.));
    assert!((noisy.volume_mm3().unwrap() - slot.volume_mm3().unwrap()).abs() < 1e-9);
    assert_eq!(noisy.faces().unwrap().len(), 8);

    // The axis meets a vertical stock edge at a right angle: the notch
    // removes twice the integral of the offset over a half disk.
    let notch = subtract(&stock, &cylinder([-2., 6., 3.], [1., 1., 0.], 1., 30.));
    assert!((480. - notch.volume_mm3().unwrap() - 4. / 3.).abs() < 1e-9);
    certified_volume(&notch);
}

#[test]
fn disjoint_contact_small_penetration_and_containment_are_distinguished() {
    let stock = stock();
    let s = 2f64.sqrt().recip();
    let planar = |z: f64| {
        let frame = PlaneFrame {
            origin: [12., 9., z],
            u: [s, -s, 0.],
            v: [0., 0., 1.],
        };
        Solid::cylinder(frame, 1., 30., policy()).unwrap()
    };
    // Disjoint, and tangent to the top face along one generator.
    for tool in [cylinder([30., 4., 3.], [1., 1., 1.], 1., 5.), planar(7.)] {
        let kept = subtract(&stock, &tool);
        assert_eq!(kept.volume_mm3().unwrap(), 480.);
        assert_eq!(kept.faces().unwrap().len(), 6);
        assert!(kept.rust_prism().is_some());
    }
    // Exact contact needs an exact chart: an integer translation keeps it,
    // a rotated chart or a decimal radius must refuse instead of guessing.
    let shift = RigidTransform {
        translation_mm: [4096., -512., 64.],
        ..Default::default()
    };
    let kept = subtract(
        &stock.transformed(shift).unwrap(),
        &planar(7.).transformed(shift).unwrap(),
    );
    assert_eq!(kept.volume_mm3().unwrap(), 480.);
    let turn = RigidTransform {
        rotation: rotation([1., 2., 3.], 0.7),
        translation_mm: [4., -3., 2.],
    };
    let thin = PlaneFrame {
        origin: [12., 9., 6.1],
        u: [s, -s, 0.],
        v: [0., 0., 1.],
    };
    let thin = Solid::cylinder(thin, 0.1, 30., policy()).unwrap();
    for (target, tool) in [
        (
            stock.transformed(turn).unwrap(),
            planar(7.).transformed(turn).unwrap(),
        ),
        (stock.share().unwrap(), thin),
    ] {
        let error = target
            .boolean(&tool, BooleanOperation::Subtract, policy())
            .err()
            .unwrap();
        assert!(error.message.contains("unresolved"), "{error}");
    }
    // One unit in the last place of penetration is never contact.
    let error = stock
        .boolean(
            &planar(7f64.next_down()),
            BooleanOperation::Subtract,
            policy(),
        )
        .err()
        .unwrap();
    assert_eq!(error.kind, GeometryErrorKind::Unsupported);
    assert!(error.message.contains("unresolved"), "{error}");
    // Resolvable positive material is removed.
    let grazed = subtract(&stock, &planar(7. - 1e-4));
    let removed = 480. - grazed.volume_mm3().unwrap();
    assert!(removed > 1e-6 && removed < 1e-4, "{removed}");
    assert_eq!(grazed.faces().unwrap().len(), 8);
    assert!(certified_volume(&grazed) < 480.);
    // A tool containing the stock leaves no solid.
    let small = Solid::box_solid([1.; 3], policy()).unwrap();
    let cover = cylinder([-5.; 3], [1., 1., 1.], 10., 30.);
    assert!(
        small
            .boolean(&cover, BooleanOperation::Subtract, policy())
            .unwrap()
            .is_none()
    );
    assert_eq!(stock.volume_mm3().unwrap(), 480.);
}

#[test]
fn invalid_unresolved_cancelled_and_over_budget_requests_preserve_the_source() {
    let stock = stock();
    let good = cylinder([3., 4., -2.], [0.2, 0.1, 1.], 1., 20.);
    let check = |stock: &Solid| {
        assert_eq!(stock.volume_mm3().unwrap(), 480.);
        assert_eq!(stock.faces().unwrap().len(), 6);
        assert!(stock.topology_snapshot(EntityId::new(), 0).is_ok());
        // The source still accepts the admitted cut.
        assert_eq!(subtract(stock, &good).faces().unwrap().len(), 7);
    };
    let refused = |tool: &Solid| {
        let error = stock
            .boolean(tool, BooleanOperation::Subtract, policy())
            .err()
            .unwrap();
        check(&stock);
        error
    };
    // Near-tangent: the side grazes a vertical stock edge from outside.
    let grazing = cylinder([10. + 2f64.sqrt(), 0., 3.], [-1., 1., 0.], 1., 30.);
    let error = refused(&grazing);
    assert_eq!(error.kind, GeometryErrorKind::Unsupported);
    assert!(error.message.contains("unresolved"));
    // Nearly parallel to a face pair, beyond the chart noise.
    let error = refused(&cylinder([3., 4., -2.], [1., 1., 1e-6], 1., 20.));
    assert!(error.message.contains("unresolved"));
    // A cut that would split the stock or leave an enclosed void.
    for tool in [
        cylinder([-5., -6., 3.], [1., 1., 0.], 3.5, 30.),
        cylinder([3., 3., 3.], [1., 1., 1.], 0.5, 2.),
    ] {
        let error = refused(&tool);
        assert_eq!(error.kind, GeometryErrorKind::Unsupported);
        assert!(error.message.contains("one connected result shell"));
    }
    // Stock that already carries a round cut is outside this provider.
    let drilled = subtract(&stock, &cylinder([5., 4., -1.], [0., 0., 1.], 1., 8.));
    let error = drilled
        .boolean(&good, BooleanOperation::Subtract, policy())
        .err()
        .unwrap();
    assert_eq!(error.kind, GeometryErrorKind::Unsupported);

    // A tool frame inside the constructor tolerance but outside the
    // certified chart is invalid input for this provider.
    let mut sloppy = frame([3., 4., -2.], [0.2, 0.1, 1.]);
    sloppy.u = sloppy.u.map(|c| c * (1. + 2e-11));
    let sloppy = Solid::cylinder(sloppy, 1., 20., policy()).unwrap();
    assert_eq!(refused(&sloppy).kind, GeometryErrorKind::Input);
    // Polygonal stock that is not a plain box stays with its own provider.
    let lengths = |values: [i64; 2]| values.map(|n| Rational::new(n, 1).unwrap());
    let shapes = [
        vec![vec![lengths([0, 0]), lengths([10, 0]), lengths([0, 8])]],
        vec![
            vec![
                lengths([0, 0]),
                lengths([10, 0]),
                lengths([10, 8]),
                lengths([0, 8]),
            ],
            vec![
                lengths([2, 2]),
                lengths([2, 4]),
                lengths([4, 4]),
                lengths([4, 2]),
            ],
        ],
        vec![vec![
            lengths([1, 0]),
            lengths([10, 0]),
            lengths([10, 8]),
            lengths([1, 8]),
        ]],
    ];
    for rings in shapes {
        let feature = EntityId::new();
        let count = rings.iter().map(Vec::len).sum::<usize>();
        let curves: Vec<Origin> = (0..count)
            .map(|i| Origin {
                source: feature,
                role: format!("edge.{i}"),
            })
            .collect();
        let shape = Solid::from_polygon_prism(
            PolygonPrism::new(rings, Rational::new(6, 1).unwrap()).unwrap(),
            PlaneFrame::default(),
            feature,
            &curves,
            &vec![false; count],
            policy(),
        )
        .unwrap();
        let error = shape
            .boolean(&good, BooleanOperation::Subtract, policy())
            .err()
            .unwrap();
        assert_eq!(error.kind, GeometryErrorKind::Unsupported, "{error}");
        let direct = shape.subtract_oblique_cylinder(&good, policy(), None);
        assert_eq!(direct.err().unwrap().kind, GeometryErrorKind::Unsupported);
    }

    // Invalid requests of the interruptible entry point.
    let fuzzy = GeometryPolicy {
        boolean_fuzzy_mm: 1e-8,
        ..policy()
    };
    let other = GeometryPolicy {
        minimum_feature_mm: 2e-5,
        ..policy()
    };
    let invalid = GeometryPolicy {
        linear_tolerance_mm: f64::NAN,
        ..policy()
    };
    for policy in [fuzzy, other, invalid] {
        let error = stock
            .subtract_oblique_cylinder(&good, policy, None)
            .err()
            .unwrap();
        assert_eq!(error.kind, GeometryErrorKind::Input);
    }
    let aligned = cylinder([5., 4., -1.], [0., 0., 1.], 1., 8.);
    for (target, tool) in [(&stock, &aligned), (&stock, &stock), (&drilled, &good)] {
        let error = target
            .subtract_oblique_cylinder(tool, policy(), None)
            .err()
            .unwrap();
        assert_eq!(error.kind, GeometryErrorKind::Unsupported);
    }
    check(&stock);

    // Cancellation at the first and at a later checkpoint.
    for stop in [1, 20] {
        let calls = Cell::new(0);
        let interrupt = || {
            calls.set(calls.get() + 1);
            calls.get() >= stop
        };
        let error = stock
            .subtract_oblique_cylinder(&good, policy(), Some(&interrupt))
            .err()
            .unwrap();
        assert_eq!(error.kind, GeometryErrorKind::Cancelled);
        assert_eq!(calls.get(), stop);
        check(&stock);
    }
    let never = || false;
    let cut = stock
        .subtract_oblique_cylinder(&good, policy(), Some(&never))
        .unwrap()
        .unwrap();
    assert_eq!(
        cut.volume_mm3().unwrap(),
        subtract(&stock, &good).volume_mm3().unwrap()
    );

    // Tessellation budget exhaustion leaves the cut body usable.
    let tight = MeshOptions {
        maximum_triangles: 64,
        ..OPTIONS
    };
    for error in [cut.mesh(tight).err(), cut.manufacturing_mesh(tight).err()] {
        assert_eq!(error.unwrap().kind, GeometryErrorKind::Limit);
    }
    let invalid = MeshOptions {
        linear_deflection_mm: 0.,
        ..OPTIONS
    };
    assert_eq!(
        cut.mesh(invalid).err().unwrap().kind,
        GeometryErrorKind::Input
    );
    certified_volume(&cut);
    check(&stock);
}

#[test]
fn unsupported_downstream_operations_refuse_explicitly() {
    let stock = stock();
    let cut = subtract(&stock, &cylinder([3., 4., -2.], [0.2, 0.1, 1.], 1., 20.));
    let volume = cut.volume_mm3().unwrap();

    // Exchange: no STEP encoding is claimed for this boundary.
    assert!(!cut.supports_step_export());
    let scratch = tempfile::tempdir().unwrap();
    let path = scratch.path().join("oblique.step");
    let error = cut.write_step(&path).err().unwrap();
    assert_eq!(error.kind, GeometryErrorKind::Unsupported);
    assert!(!path.exists());

    // Further Booleans in either role.
    let block = Solid::box_solid([2.; 3], policy()).unwrap();
    let second = cylinder([7., 4., -2.], [0.1, 0.2, 1.], 0.5, 20.);
    let upright = cylinder([8., 2., -1.], [0., 0., 1.], 0.5, 8.);
    for operation in [
        BooleanOperation::Union,
        BooleanOperation::Subtract,
        BooleanOperation::Intersect,
    ] {
        for other in [&block, &second, &upright] {
            let error = cut.boolean(other, operation, policy()).err().unwrap();
            assert_eq!(error.kind, GeometryErrorKind::Unsupported, "{error}");
            let error = other.boolean(&cut, operation, policy()).err().unwrap();
            assert_eq!(error.kind, GeometryErrorKind::Unsupported, "{error}");
        }
    }

    // Fillets and chamfers on a retained edge.
    let owner = EntityId::new();
    let named = Solid::box_named([10., 8., 6.], owner, policy()).unwrap();
    let named = subtract(&named, &cylinder([3., 4., -2.], [0.2, 0.1, 1.], 1., 20.));
    let snapshot = named.topology_snapshot(owner, 0).unwrap();
    let edge = elements(&snapshot, TopologyKind::Edge)
        .into_iter()
        .find(|e| {
            matches!(e.geometry, Some(AnalyticGeometry::Line { .. })) && !e.origins.is_empty()
        })
        .unwrap();
    let feature = EntityId::new();
    let selected = [TreatmentEdge {
        key: edge.key,
        origin: Origin {
            source: feature,
            role: "edge.0".into(),
        },
    }];
    for treatment in [
        EdgeTreatment::Fillet { radius_mm: 0.5 },
        EdgeTreatment::Chamfer { distance_mm: 0.5 },
    ] {
        let error = named
            .treat_edges_named(&selected, treatment, feature, policy())
            .err()
            .unwrap();
        println!("{error}");
        assert_eq!(error.kind, GeometryErrorKind::Unsupported, "{error}");
    }
    // Queries without a certified provider for this boundary refuse.
    let far = block
        .transformed(RigidTransform {
            translation_mm: [50., 0., 0.],
            ..Default::default()
        })
        .unwrap();
    let whole = BRepSelection::Whole;
    let refusals = [
        cut.inspect_pair(&block, 1e-6, policy()).err(),
        cut.rust_occupied_overlap_volume(&block, policy()).err(),
        cut.minimum_distance(whole, &far, whole).err(),
        cut.split_by_plane(frame([5., 0., 0.], [1., 0., 0.]), policy())
            .err(),
        cut.mirrored_named(PlaneFrame::default(), EntityId::new(), policy())
            .err(),
        cut.to_brep().err(),
    ];
    for refusal in refusals {
        assert_eq!(refusal.unwrap().kind, GeometryErrorKind::Unsupported);
    }
    for face in 0..7 {
        let plane = frame([0.; 3], [0., 0., 1.]);
        assert!(cut.planar_face_contours(face, plane).is_err());
    }
    // Separated bounds still certify zero overlap.
    assert_eq!(
        cut.rust_occupied_overlap_volume(&far, policy()).unwrap(),
        0.
    );
    // Exact display edges come from the retained lines and ellipses.
    let display = cut.topology_with_display(owner, 0, Some(0.01)).unwrap();
    assert_eq!([display.edges.len(), display.vertices.len()], [14, 10]);
    for edge in &display.edges {
        let on_bore = |p: &[f64; 3]| {
            let w = normalized([0.2, 0.1, 1.]);
            let d = [p[0] - 3., p[1] - 4., p[2] + 2.];
            let radial = cross(d, w);
            (dot(radial, radial).sqrt() - 1.).abs() < 1e-9
        };
        if edge.points_mm.len() > 2 {
            assert!(edge.points_mm.iter().all(on_bore));
            let on_cap = |p: &[f64; 3]| p[2].abs() < 1e-9 || (p[2] - 6.).abs() < 1e-9;
            assert!(edge.points_mm.iter().all(on_cap));
        }
    }
    // A compound keeps the cut as an exact member and inherits its limits.
    let members = [&cut, &far];
    let compound = Solid::compound_named(&members, EntityId::new(), policy()).unwrap();
    assert!((compound.volume_mm3().unwrap() - volume - 8.).abs() < 1e-9);
    assert!(!compound.supports_step_export());
    assert!(compound.manufacturing_mesh(OPTIONS).is_ok());
    // Nothing above changed the cut.
    assert_eq!(cut.volume_mm3().unwrap(), volume);
    certified_volume(&cut);
}
