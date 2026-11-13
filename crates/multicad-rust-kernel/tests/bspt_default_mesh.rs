use multicad_rust_kernel::{
    axial_pipe_thread::authored::{AuthoredSection, AuthoredSpec, AuthoredThread},
    display::Tessellation,
    manufacturing,
};
#[test]
fn canonical_eight_mm_bspt_has_a_closed_source_certified_manufacturing_mesh() {
    let pitch = 25.4 / 28.;
    let source = AuthoredThread::new(
        AuthoredSpec {
            major_radius_mm: 4.739,
            pitch_mm: pitch,
            length_mm: 8.,
            radial_taper: 1. / 32.,
            section: AuthoredSection::BsptRounded {
                depth_mm: 0.5808684050603288,
            },
            crest_radius_limit_mm: None,
        },
        1e-7,
        1e-5,
    )
    .unwrap();
    let body = source
        .raw_boundary_ownership(1e-7, None)
        .unwrap()
        .ownership
        .bind_surfaces(1e-7, None)
        .unwrap()
        .certify_axial(1e-7, 1e-5, None)
        .unwrap()
        .certify_caps(None)
        .unwrap()
        .certify_ordering(None)
        .unwrap()
        .certify_material(1e-7, None)
        .unwrap();
    let before = body.input().clone();
    let source_volume = body.volume_with_interrupt(1e-5, None).unwrap().volume_mm3;
    let started = std::time::Instant::now();
    let mesh = body
        .tessellate_with_interrupt(
            Tessellation {
                deflection: 0.1,
                angle: 0.3,
                triangle_limit: 1_000_000,
            },
            None,
        )
        .unwrap();
    let receipt = manufacturing::certify(&mesh).unwrap();
    eprintln!(
        "public canonical BSPT8 mesh {:?}, {} triangles, volume {}",
        started.elapsed(),
        mesh.indices.len() / 3,
        receipt.signed_volume_mm3
    );
    assert!(receipt.signed_volume_mm3 > 0.);
    // Independent retained-source volume and literal oriented triangle flux
    // must agree at the requested coarse display/manufacturing chord scale.
    // This regression comparison does not replace either C0/C1 source proof.
    assert!(
        source_volume
            .iter()
            .all(|v| (receipt.signed_volume_mm3 - v).abs() < 0.01 * v)
    );
    eprintln!("independent source volume {source_volume:?}");
    assert_eq!(body.input(), &before);
}
