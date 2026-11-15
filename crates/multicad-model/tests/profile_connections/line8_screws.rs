use super::*;
use multicad_kernel::RigidTransform;
#[test]
fn rust_line8_entire_screws_clear_both_original_members_after_project_reopen() {
    for case in 0..3 {
        let (session, members, request, outputs) =
            super::line8_nuts::fixture(case, ProfileConnectionNut::ZincPlated);
        let before = session.document().clone();
        let mut bytes = Vec::new();
        project_file::write(
            &ProjectFile::new(before.clone()).unwrap(),
            &mut bytes,
            Default::default(),
        )
        .unwrap();
        let reopened = project_file::read(
            &mut bytes.as_slice(),
            Default::default(),
            &project_file::supported_capabilities(),
        )
        .unwrap();
        let document = &reopened.project.document;
        assert_eq!(document, &before);
        let revision = session.revision();
        let compiled = super::line8_nuts::compiler()
            .compile_with_backend(document, revision, GeometryBackend::Rust, None)
            .unwrap_or_else(|e| panic!("case {case}: {e}"));
        let parameters = multicad_core::parameters::evaluate(document).unwrap();
        let snapshots = compiled
            .snapshot()
            .shapes
            .iter()
            .map(|(id, shape)| (*id, shape.topology.clone()))
            .collect();
        for index in [1, 3] {
            let screw = compiled.solid(outputs[index].feature).unwrap();
            let (source, _) = screw.rust_profiled_screw_geometry().unwrap();
            assert_eq!(source.spec().shank_length_mm, 1.25);
            assert!(source.thread().spec().length_mm > 8.);
            assert_eq!(screw.solid_count().unwrap(), 1);
            for (site, member) in [
                (&request.first, &members[0]),
                (&request.second, &members[1]),
            ] {
                let pose = model::resolve_occurrence_pose(
                    document,
                    &site.occurrence,
                    site.member.body,
                    &parameters,
                    revision,
                    &snapshots,
                )
                .unwrap();
                let receipt = screw
                    .rust_profiled_screw_receiver_clearance(
                        compiled.solid(member.feature).unwrap(),
                        RigidTransform {
                            rotation: pose.rotation,
                            translation_mm: pose.translation_mm,
                        },
                        0.001,
                        None,
                    )
                    .unwrap();
                assert!(receipt.overlap_volume_upper_mm3 <= 0.001);
                assert!(std::ptr::eq(receipt.screw().source, source));
                assert_eq!(
                    receipt.boundary_owners.len(),
                    source.profile().segments().count() + 1
                );
            }
            let mesh = screw
                .manufacturing_mesh_with_receipt(Default::default())
                .unwrap();
            assert_eq!(mesh.edge_pairs as usize * 2, mesh.mesh.indices.len());
            assert_eq!(mesh.vertex_links as usize, mesh.mesh.vertices.len() / 3);
        }
        assert_eq!(session.document(), &before);
    }
}
