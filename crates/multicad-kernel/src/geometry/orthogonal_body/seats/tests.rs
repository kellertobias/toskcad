use super::*;
fn raw_fixture(offset: Point) -> (EntityId, Solid, [OrthogonalCountersink; 2]) {
    let policy = GeometryPolicy::default();
    let owner = EntityId::new();
    let geometry = orthogonal::boolean(
        &[
            Box3 {
                min: offset,
                max: std::array::from_fn(|i| offset[i] + [40., 40., 5.][i]),
            },
            Box3 {
                min: offset,
                max: std::array::from_fn(|i| offset[i] + [40., 5., 40.][i]),
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
            frame: PlaneFrame {
                origin: std::array::from_fn(|i| origin[i] + offset[i]),
                u,
                v,
            },
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
    (owner, drilled, seats)
}
fn fixture() -> (EntityId, Solid, Solid, [OrthogonalCountersink; 2]) {
    let (owner, drilled, seats) = raw_fixture([0.; 3]);
    let body = drilled
        .rust_countersink_orthogonal(&seats, GeometryPolicy::default(), None)
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
#[test]
fn actual_conical_azimuth_weights_have_policy_bounded_centroid_covariance() {
    let (_, _, body, _) = fixture();
    let shear = 2_f64.powi(-42);
    let actual = body
        .transformed(RigidTransform {
            rotation: [[1., 0., 0.], [shear, 1., 0.], [0., 0., 1.]],
            translation_mm: [11., -7., 5.],
        })
        .unwrap();
    let measured = actual.rust_orthogonal_countersink_measures().unwrap();
    assert!(
        measured.maximum_centroid_error_mm <= GeometryPolicy::default().linear_tolerance_mm / 4.
    );
    let (source, frame) = actual.rust_orthogonal_countersink_source().unwrap();
    let seat = source.seats()[1];
    let bore = source.bores()[seat.bore];
    let end = usize::from(seat.end);
    let center = bore.centers_mm[end];
    let normal = source.boundary().faces[bore.faces[end]].normal;
    let axis = (0..3).find(|&i| normal[i] != 0.).unwrap();
    let u = (axis + 1) % 3;
    let v = (axis + 2) % 3;
    let axes = [frame.u, frame.v, frame.normal()];
    let cof = [
        cross(axes[1], axes[2]),
        cross(axes[2], axes[0]),
        cross(axes[0], axes[1]),
    ];
    let slope = (seat.radius_mm - bore.radius_mm) / seat.depth_mm;
    let mut weight = 0.;
    let mut cosine = 0.;
    let mut sine = 0.;
    let mut min = f64::INFINITY;
    let mut max = 0_f64;
    // Numerical regression witness only. The production enclosure above uses
    // whole-support Jacobian bounds, independent of these quadrature samples.
    for i in 0..4096 {
        let theta = TAU * (i as f64 + 0.5) / 4096.;
        let (s, c) = theta.sin_cos();
        let mut n = normal.map(|x| x * slope);
        n[u] -= c;
        n[v] -= s;
        let world =
            std::array::from_fn::<_, 3, _>(|j| (0..3).map(|k| cof[k][j] * n[k]).sum::<f64>());
        let w = world[0].hypot(world[1]).hypot(world[2]);
        weight += w;
        cosine += w * c;
        sine += w * s;
        min = min.min(w);
        max = max.max(w);
    }
    assert!(max - min > 1e-14);
    let r = bore.radius_mm;
    let outer = seat.radius_mm;
    let radial = 2. * (outer * outer + outer * r + r * r) / (3. * (outer + r));
    let depth = seat.depth_mm * (outer + 2. * r) / (3. * (outer + r));
    let mut local = std::array::from_fn(|j| center[j] - normal[j] * depth);
    local[u] += radial * cosine / weight;
    local[v] += radial * sine / weight;
    let numerical = transform_point(frame, local);
    let index = source.boundary().faces.len() + source.bores().len() + 1;
    let bound = measured.face_centroid_mm[index];
    for j in 0..3 {
        assert!(
            bound[j][0] <= numerical[j] && numerical[j] <= bound[j][1],
            "{bound:?} vs {numerical:?}"
        );
    }
    let face = actual.faces().unwrap()[index].clone();
    for (j, b) in bound.iter().enumerate() {
        assert!(b[0] <= face.center_mm[j] && face.center_mm[j] <= b[1]);
    }
}

#[test]
fn ill_conditioned_original_moments_refuse_without_widening_policy_or_mutating_stock() {
    let (owner, drilled, seats) = raw_fixture([1e6; 3]);
    let before = drilled.topology_snapshot(owner, 0).unwrap();
    let error = drilled
        .rust_countersink_orthogonal(&seats, GeometryPolicy::default(), None)
        .err()
        .unwrap();
    assert_eq!(error.kind, GeometryErrorKind::Unsupported);
    assert!(error.message.contains("centroid"), "{error}");
    assert_eq!(before, drilled.topology_snapshot(owner, 0).unwrap());
    // Authoring the same local stock canonically and placing the admitted body
    // keeps finite moment conditioning; no existing source is recentered.
    let (_, _, canonical, _) = fixture();
    let placed = canonical
        .transformed(RigidTransform {
            translation_mm: [1e6; 3],
            ..Default::default()
        })
        .unwrap();
    let bound = placed.rust_orthogonal_countersink_measures().unwrap();
    assert!(bound.maximum_centroid_error_mm <= GeometryPolicy::default().linear_tolerance_mm / 4.);
}
