use super::*;
#[test]
fn actual_source_solid_publishes_annular_topology_faces_and_closed_manufacturing() {
    for frame in [
        PlaneFrame::default(),
        PlaneFrame {
            origin: [11., -7., 5.],
            u: [0.8, 0.6, 0.],
            v: [-0.36, 0.48, 0.8],
        },
    ] {
        let source = crate::geometry::wood_insert_world_tests::source(false);
        let feature = EntityId::new();
        let timer = std::time::Instant::now();
        let solid =
            Solid::rust_wood_insert(source, feature, frame, GeometryPolicy::default(), None)
                .unwrap();
        assert_eq!(solid.solid_count().unwrap(), 1);
        assert!(solid.volume_mm3().unwrap() > 0.);
        assert_eq!(
            solid.rust_wood_insert_data(),
            Some(WoodInsertData { source, frame })
        );
        let topology = solid
            .topology_with_display(EntityId::new(), 7, Some(0.1))
            .unwrap();
        assert_eq!(
            topology
                .topology
                .elements
                .iter()
                .filter(|e| e.key.kind == TopologyKind::Face)
                .count(),
            source.face_count() as usize
        );
        let faces = solid.faces().unwrap();
        assert_eq!(faces.len(), source.face_count() as usize);
        assert!(faces.iter().all(|f| f.area_mm2 > 0.));
        assert!(faces[0].plane.is_some() && faces[1].plane.is_some());
        assert!(!topology.edges.is_empty() && !topology.vertices.is_empty());
        let proof = solid
            .manufacturing_mesh_with_receipt(MeshOptions::default())
            .unwrap();
        assert!(proof.edge_pairs > 0 && proof.vertex_links > 0);
        assert_eq!(proof.removed_exact_zero_triangles, 0);
        assert_eq!(solid.rust_wood_insert_data().unwrap().source, source);
        eprintln!(
            "actual full insert Solid {}faces {}edges {}triangles {:?}",
            faces.len(),
            topology.edges.len(),
            proof.mesh.indices.len() / 3,
            timer.elapsed()
        );
    }
}
