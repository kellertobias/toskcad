use super::*;
fn fixture() -> (EntityId, Solid, Solid, [OrthogonalCountersink; 2]) {
    let policy = GeometryPolicy::default();
    let owner = EntityId::new();
    let geometry = orthogonal::boolean(
        &[
            Box3 {
                min: [0.; 3],
                max: [40., 40., 5.],
            },
            Box3 {
                min: [0.; 3],
                max: [40., 5., 40.],
            },
        ],
        &[],
        Operation::Union,
        policy.linear_tolerance_mm,
    )
    .unwrap()
    .unwrap();
    let stock = Solid::rust_orthogonal_stock_named(owner, geometry, policy).unwrap();
    let tools = [
        [[20., 20., -0.1], [0., 0., 5.2], [1., 0., 0.], [0., 1., 0.]],
        [[20., -0.1, 20.], [0., 5.2, 0.], [0., 0., 1.], [1., 0., 0.]],
    ]
    .into_iter()
    .enumerate()
    .map(
        |(i, [origin, displacement_mm, u, v])| crate::CircularPocketTool {
            frame: PlaneFrame { origin, u, v },
            center_uv: [0.; 2],
            radius_mm: 4.5,
            displacement_mm,
            origin: Origin {
                source: owner,
                role: format!("angle.bore.{i}"),
            },
        },
    )
    .collect::<Vec<_>>();
    let drilled = stock
        .rust_drill_circular_machining(&tools, &[true, true], policy)
        .unwrap();
    let seats = [0, 1].map(|bore| OrthogonalCountersink {
        geometry: Seat {
            bore,
            end: false,
            radius_mm: 8.05,
            depth_mm: 3.55,
        },
        origin: Origin {
            source: owner,
            role: format!("angle.countersink.{bore}"),
        },
    });
    let body = drilled
        .rust_countersink_orthogonal(&seats, policy, None)
        .unwrap();
    (owner, drilled, body, seats)
}
#[test]
fn actual_connected_perpendicular_countersinks_have_one_native_manufacturing_shell() {
    let (owner, _drilled, body, _seats) = fixture();
    let (source, _) = body.rust_orthogonal_countersink_source().unwrap();
    assert_eq!(source.bores().len(), 2);
    assert_eq!(source.seats().len(), 2);
    assert_eq!(body.solid_count().unwrap(), 1);
    let snapshot = body.topology_snapshot(owner, 7).unwrap();
    snapshot.validate().unwrap();
    let faces = body.faces().unwrap();
    assert_eq!(faces.len(), source.boundary().faces.len() + 4);
    assert!(
        faces
            .iter()
            .all(|f| f.area_mm2.is_finite() && f.area_mm2 > 0.)
    );
    let receipt = body
        .manufacturing_mesh_with_receipt(Default::default())
        .unwrap();
    assert_eq!(receipt.edge_pairs as usize * 2, receipt.mesh.indices.len());
    assert_eq!(
        receipt.vertex_links as usize,
        receipt.mesh.vertices.len() / 3
    );
    assert_eq!(receipt.maximum_displacement_mm, 0.);
}
#[test]
fn actual_nonisometric_stored_columns_transport_mass_area_and_plane_normals() {
    let (_, _, body, _) = fixture();
    let original = body.rust_orthogonal_countersink_measures().unwrap();
    let local_faces = body.faces().unwrap();
    let cap = local_faces
        .iter()
        .find(|f| {
            f.plane
                .is_some_and(|(_, n)| n[0].abs() == 1. && n[1] == 0. && n[2] == 0.)
        })
        .unwrap()
        .index as usize;
    let s = 1. + 2_f64.powi(-42);
    let scaled = body
        .transformed(RigidTransform {
            rotation: [[s, 0., 0.], [0., 1., 0.], [0., 0., 1.]],
            translation_mm: [11., -7., 5.],
        })
        .unwrap();
    let measured = scaled.rust_orthogonal_countersink_measures().unwrap();
    assert!(measured.determinant[0] > 1.);
    assert!(measured.volume_mm3[0] > original.volume_mm3[1]);
    assert!(measured.face_area_mm2[cap][0] > original.face_area_mm2[cap][1]);
    assert!(scaled.volume_mm3().unwrap() >= measured.volume_mm3[0]);
    assert!(std::ptr::eq(
        scaled.rust_orthogonal_countersink_source().unwrap().0,
        body.rust_orthogonal_countersink_source().unwrap().0
    ));
    let sheared = body
        .transformed(RigidTransform {
            rotation: [[1., 0., 0.], [2_f64.powi(-42), 1., 0.], [0., 0., 1.]],
            ..Default::default()
        })
        .unwrap();
    let n = sheared.faces().unwrap()[cap].plane.unwrap().1;
    assert!(n[1].abs() < 1e-15);
    assert!(n[2].abs() < 1e-14);
    let receipt = scaled
        .manufacturing_mesh_with_receipt(Default::default())
        .unwrap();
    assert_eq!(receipt.edge_pairs as usize * 2, receipt.mesh.indices.len());
    assert_eq!(
        receipt.vertex_links as usize,
        receipt.mesh.vertices.len() / 3
    );
    let bad = body.transformed(RigidTransform {
        rotation: [[1. + 1e-8, 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        ..Default::default()
    });
    assert_eq!(bad.err().unwrap().kind, GeometryErrorKind::Input);
}
#[test]
fn original_native_stock_and_actual_seats_remain_unchanged_on_refusal_or_cancel() {
    let (owner, drilled, body, seats) = fixture();
    let before = drilled.topology_snapshot(owner, 0).unwrap();
    let mut changed = seats.clone();
    changed[0].geometry.depth_mm = 5.;
    assert_eq!(
        drilled
            .rust_countersink_orthogonal(&changed, GeometryPolicy::default(), None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Unsupported
    );
    let polls = std::cell::Cell::new(0);
    assert_eq!(
        drilled
            .rust_countersink_orthogonal(
                &seats,
                GeometryPolicy::default(),
                Some(&|| {
                    let n = polls.get() + 1;
                    polls.set(n);
                    n > 4
                })
            )
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Cancelled
    );
    assert!(polls.get() > 4);
    changed = seats.clone();
    changed[1].origin = changed[0].origin.clone();
    assert_eq!(
        drilled
            .rust_countersink_orthogonal(&changed, GeometryPolicy::default(), None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Input
    );
    assert_eq!(
        body.rust_countersink_orthogonal(&seats, GeometryPolicy::default(), None)
            .err()
            .unwrap()
            .kind,
        GeometryErrorKind::Unsupported
    );
    assert_eq!(before, drilled.topology_snapshot(owner, 0).unwrap());
}
