use multicad_core::{
    document::EntityId,
    references::{GeometryIntent, Resolution, TopologyKind, resolve},
};
use multicad_kernel::{
    BlindThreadCut, BlindThreadSeat, GeometryPolicy, HelicalThread, PlaneFrame, Solid,
    rust_kernel::Rational,
};
fn n(v: i64) -> Rational {
    Rational::new(v, 1).unwrap()
}
fn point(f: PlaneFrame, p: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| f.origin[i] + f.u[i] * p[0] + f.v[i] * p[1] + f.normal()[i] * p[2])
}
fn spec(left_hand: bool) -> HelicalThread {
    HelicalThread {
        major_radius_mm: 1.5,
        pitch_mm: 0.5,
        length_mm: 3.5,
        depth_mm: 0.27,
        crest_width_mm: 0.08,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        radial_taper: 0.,
        starts: 1,
        left_hand,
    }
}
fn cut(
    stock: &Solid,
    feature: EntityId,
    entry: PlaneFrame,
    left: bool,
    depth: f64,
) -> Result<Solid, multicad_kernel::GeometryError> {
    stock.rust_cut_counterbored_blind_internal_thread(
        BlindThreadCut {
            profile: spec(left),
            entry,
            receiving_radius_mm: 1.5,
            minimum_wall_mm: 1.,
            minimum_floor_mm: 1.,
        },
        BlindThreadSeat {
            radius_mm: 3.5,
            depth_mm: depth,
        },
        feature,
        GeometryPolicy::default(),
    )
}
#[test]
fn actual_counterbored_blind_source_preserves_stock_and_finite_reference_geometry() {
    for left in [false, true] {
        for top in [false, true] {
            for frame in [
                PlaneFrame::default(),
                PlaneFrame {
                    origin: [100., -30., 50.],
                    u: [0.8, 0.6, 0.],
                    v: [-0.36, 0.48, 0.8],
                },
            ] {
                let stock_id = EntityId::new();
                let feature = EntityId::new();
                let policy = GeometryPolicy::default();
                let stock =
                    Solid::rust_box(stock_id.uuid_bits(), [n(30), n(20), n(20)], frame, policy)
                        .unwrap();
                let entry = PlaneFrame {
                    origin: point(frame, [10., 10., if top { 20. } else { 0. }]),
                    u: frame.u,
                    v: frame.v.map(|v| if top { -v } else { v }),
                };
                let original = stock.topology_snapshot(stock_id, 0).unwrap();
                let seated = cut(&stock, feature, entry, left, 0.8).unwrap();
                let snapshot = seated.topology_snapshot(stock_id, 1).unwrap();
                assert_eq!(seated.faces().unwrap().len(), 13);
                assert_eq!(
                    snapshot
                        .elements
                        .iter()
                        .filter(|e| e.key.kind == TopologyKind::Face)
                        .count(),
                    13
                );
                let shoulder = snapshot
                    .elements
                    .iter()
                    .find(|e| {
                        e.key.kind == TopologyKind::Face
                            && e.origins
                                .iter()
                                .any(|o| o.role == "thread.counterbore/shoulder")
                    })
                    .unwrap();
                let intent = GeometryIntent::capture(&snapshot, shoulder.key).unwrap();
                let edited = cut(&stock, feature, entry, left, 1.1).unwrap();
                assert!(matches!(
                    resolve(&intent, &edited.topology_snapshot(stock_id, 2).unwrap(), 2),
                    Resolution::Resolved { .. }
                ));
                let receipt = seated
                    .manufacturing_mesh_with_receipt(Default::default())
                    .unwrap();
                assert!(receipt.edge_pairs > 0 && receipt.vertex_links > 0);
                assert!(
                    (multicad_kernel::rust_kernel::manufacturing::certify(
                        &multicad_kernel::rust_kernel::display::DisplayMesh {
                            vertices: receipt.mesh.vertices.clone(),
                            indices: receipt.mesh.indices.clone(),
                            face_indices: receipt.mesh.face_indices.clone()
                        }
                    )
                    .unwrap()
                    .signed_volume_mm3
                        - seated.volume_mm3().unwrap())
                    .abs()
                        < 1.
                );
                assert_eq!(stock.topology_snapshot(stock_id, 0).unwrap(), original);
                assert!(cut(&stock, feature, entry, left, 3.5).is_err());
                assert_eq!(stock.topology_snapshot(stock_id, 0).unwrap(), original);
                // No legacy blind STEP support is inferred from a new valid mesh.
                assert!(!seated.supports_step_export());
            }
        }
    }
}
