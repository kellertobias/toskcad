//! Actual source candidates into the unchanged whole literal triangle gate.
use super::*;
use multicad_rust_kernel::{
    axial_pipe_thread::authored::{
        AuthoredSection, AuthoredSpec, AuthoredThread, StockClipOptions,
    },
    display::Tessellation,
};
fn actual_source(quarter: bool) -> AuthoredThread {
    let (tpi, e1, l1, l2, od, max) = if quarter {
        (18., 0.49163, 0.2278, 0.4018, 0.540, 0.088)
    } else {
        (27., 0.37360, 0.1615, 0.2639, 0.405, 0.096)
    };
    let p = 25.4 / tpi;
    let k = 1. / 32.;
    let m = 3_f64.sqrt();
    let sides = 1. / (m - k) + 1. / (m + k);
    let depth = p / sides - (0.033 + max) * p;
    AuthoredThread::new(
        AuthoredSpec {
            major_radius_mm: (e1 * 25.4 - 2. * k * l1 * 25.4 + depth) / 2.,
            pitch_mm: p,
            length_mm: l2 * 25.4,
            radial_taper: k,
            section: AuthoredSection::NptFlat {
                depth_mm: depth,
                flat_width_mm: (p - depth * sides) / 2.,
            },
            crest_radius_limit_mm: Some(od * 25.4 / 2.),
        },
        1e-7,
        1e-5,
    )
    .unwrap()
}
#[test]
fn actual_clipped_npt_source_candidate_has_whole_literal_embedding() {
    for quarter in [false, true] {
        let source = actual_source(quarter);
        let atlas = source
            .stock_clip_plan(StockClipOptions::default(), None)
            .unwrap();
        for (deflection, angle, n) in [(0.1, 0.3, 128), (0.02, 0.5, 256)] {
            let plan = atlas
                .cell_plan(n, 128, 512_000_000, 64 * 1024 * 1024, None)
                .unwrap();
            let candidate = plan
                .mesh_candidate(
                    Tessellation {
                        deflection,
                        angle,
                        triangle_limit: 1_000_000,
                    },
                    512_000_000,
                    64 * 1024 * 1024,
                    None,
                )
                .unwrap();
            let data = candidate.mesh();
            let mesh = Mesh {
                vertices: data.vertices.clone(),
                indices: data.indices.clone(),
                face_indices: data.face_indices.clone(),
            };
            let timer = std::time::Instant::now();
            let pairs = certify(&mesh).unwrap();
            eprintln!(
                "NPT whole literal {quarter}/{deflection}/{angle}:{} facets,{pairs} pairs,{:?}",
                mesh.indices.len() / 3,
                timer.elapsed()
            );
            assert!(pairs > 0);
            assert_eq!(
                candidate.connectivity().vertex_links as usize,
                mesh.vertices.len() / 3
            );
            let mut bad = mesh.clone();
            let first: [u32; 3] = bad.indices[..3].try_into().unwrap();
            bad.indices[3..6].copy_from_slice(&first);
            assert!(certify(&bad).is_err());
        }
        assert_eq!(*atlas.source(), source);
    }
}
