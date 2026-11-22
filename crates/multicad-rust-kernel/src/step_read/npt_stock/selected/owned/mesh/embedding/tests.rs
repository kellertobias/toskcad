use super::*;
use multicad_precision::PortableMath as _;
fn check(mesh: &DisplayMesh, caps: &[Cap; 2], bytes: usize) -> Result<LiteralEmbeddingProof> {
    let mut budget = Budget {
        used: 0,
        remaining: WORK,
        bytes,
        peak: bytes,
        limit: BYTES,
        interrupt: None,
    };
    certify(mesh, caps, &mut budget)
}
#[test]
fn actual_literal_embedding_refuses_swaps_folds_double_winding_and_inward_apex() {
    let source = super::super::tests::source(false, false);
    let plan = source.geometry_plan(Default::default(), None).unwrap();
    let candidate = plan.mesh_candidate(Default::default(), None).unwrap();
    let mesh = candidate.mesh();
    let caps = &candidate.caps;
    let bytes = candidate.retained_phase_bytes;
    let proof = check(mesh, caps, bytes).unwrap();
    assert!(
        proof.minimum_radial_determinant() > 0. && proof.minimum_axis_height_clearance_mm() > 0.
    );
    let mut swapped = mesh.clone();
    let a = caps[0].rim[0] as usize;
    let b = caps[0].rim[1] as usize;
    for k in 0..3 {
        swapped.vertices.swap(a * 3 + k, b * 3 + k);
    }
    assert!(crate::manufacturing::certify(&swapped).is_ok());
    assert_eq!(check(&swapped, caps, bytes).err(), Some(Error::Unsupported));
    let mut folded = mesh.clone();
    let index = (0..mesh.vertices.len() / 3)
        .find(|i| {
            caps.iter()
                .all(|cap| cap.apex as usize != *i && !cap.rim.contains(&(*i as u32)))
        })
        .unwrap();
    for k in 0..2 {
        folded.vertices[index * 3 + k] = -folded.vertices[index * 3 + k];
    }
    assert!(crate::manufacturing::certify(&folded).is_ok());
    assert_eq!(check(&folded, caps, bytes).err(), Some(Error::Unsupported));
    let mut double = mesh.clone();
    for &index in &caps[0].rim {
        let at = index as usize * 3;
        let x = mesh.vertices[at];
        let y = mesh.vertices[at + 1];
        let radius = x.portable_hypot(y);
        let angle = 2. * y.portable_atan2(x);
        double.vertices[at] = radius * angle.portable_cos();
        double.vertices[at + 1] = radius * angle.portable_sin();
    }
    assert!(crate::manufacturing::certify(&double).is_ok());
    assert_eq!(check(&double, caps, bytes).err(), Some(Error::Unsupported));
    let mut inward = mesh.clone();
    inward.vertices[caps[0].apex as usize * 3 + 2] =
        mesh.vertices[caps[0].rim[0] as usize * 3 + 2] + 0.001;
    assert!(crate::manufacturing::certify(&inward).is_ok());
    assert_eq!(check(&inward, caps, bytes).err(), Some(Error::Unsupported));
    let interrupted = Budget {
        used: 0,
        remaining: WORK,
        bytes,
        peak: bytes,
        limit: BYTES,
        interrupt: Some(&|| true),
    };
    let mut budget = interrupted;
    assert_eq!(
        certify(mesh, caps, &mut budget).err(),
        Some(Error::Cancelled)
    );
    let mut budget = Budget {
        used: 0,
        remaining: 1,
        bytes,
        peak: bytes,
        limit: BYTES,
        interrupt: None,
    };
    assert_eq!(certify(mesh, caps, &mut budget).err(), Some(Error::Limit));
    let required = candidate.work_used();
    let constrained = plan
        .mesh_candidate(
            LiteralMeshOptions {
                maximum_work: required + 1,
                ..Default::default()
            },
            None,
        )
        .unwrap();
    assert_eq!(
        constrained.certify_embedding(None).err(),
        Some(Error::Limit)
    );
    let verified = candidate.certify_embedding(None).unwrap();
    assert!(verified.candidate().work_used() > required);
    assert_eq!(
        verified.candidate().work_used() + verified.candidate().remaining_work(),
        WORK
    );
    assert!(std::sync::Arc::ptr_eq(&source.0, &verified.source().0));
}
