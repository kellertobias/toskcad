use multicad_core::{
    commands::Session,
    document::Document,
    model::{self, ThreadFamily, ThreadFit, ThreadHand, ThreadMode, ThreadParameters, length},
};
use multicad_kernel::rust_kernel::{
    display::Tessellation,
    rounded_thread::{Spec, Thread, ThroughBore},
    rounded_thread_mesh,
};
fn canonical(family: ThreadFamily, left: bool) -> ThroughBore {
    let mut p = if family == ThreadFamily::Round {
        model::round_thread_presets()[1].parameters()
    } else {
        let row = model::whitworth_thread_presets()
            .into_iter()
            .find(|p| p.family == family)
            .unwrap();
        ThreadParameters {
            family,
            major_diameter: length(row.diameter_mm),
            pitch: length(row.pitch_mm),
            fit: if family == ThreadFamily::Bspp {
                ThreadFit::BsppInternal
            } else {
                ThreadFit::WhitworthMediumInternal
            },
            ..Default::default()
        }
    };
    p.hand = if left {
        ThreadHand::Left
    } else {
        ThreadHand::Right
    };
    let s = Session::new(Document::new("Literal rounded receiver").unwrap()).unwrap();
    let r = p.resolve(s.parameters(), ThreadMode::CutInternal).unwrap();
    p.length = length(2.37 * r.pitch_mm);
    let r = p.resolve(s.parameters(), ThreadMode::CutInternal).unwrap();
    let (crest, root) = if let Some(round) = r.round {
        (round.cutter_crest_radius_mm, round.cutter_root_radius_mm)
    } else {
        let q = r.profile.rounding_radius_mm.unwrap();
        (q, q)
    };
    let t = Thread::new(Spec {
        major_radius_mm: r.major_radius_mm,
        pitch_mm: r.pitch_mm,
        length_mm: r.length_mm,
        depth_mm: r.profile.depth_mm,
        crest_radius_mm: crest,
        root_radius_mm: root,
        flank_angle_deg: r.profile.lower_flank_deg,
        starts: r.starts,
        left_hand: left,
    })
    .unwrap();
    ThroughBore::new(t, r.major_radius_mm + 2.).unwrap()
}
#[test]
fn canonical_rounded_receiver_keeps_literal_internal_cutter_and_closed_mesh() {
    for family in [
        ThreadFamily::Bsw,
        ThreadFamily::Bsf,
        ThreadFamily::Bspp,
        ThreadFamily::Round,
    ] {
        for left in [false, true] {
            let bore = canonical(family, left);
            let t = bore.thread();
            let s = t.spec();
            assert_eq!(s.left_hand, left);
            assert!(bore.volume_mm3() > 0.);
            let mesh = rounded_thread_mesh::tessellate_bore(
                &bore,
                Tessellation {
                    deflection: 0.05,
                    angle: 0.3,
                    triangle_limit: 500_000,
                },
            )
            .unwrap();
            assert_eq!(
                mesh.face_indices
                    .iter()
                    .copied()
                    .collect::<std::collections::BTreeSet<_>>(),
                (0..8).collect()
            );
            if family == ThreadFamily::Round {
                assert_ne!(s.crest_radius_mm, s.root_radius_mm);
            };
        }
    }
}
