use multicad_core::{
    commands::Session,
    document::Document,
    model::{self, ThreadFamily, ThreadFit, ThreadHand, ThreadMode, ThreadParameters, length},
};
use multicad_kernel::rust_kernel::rounded_thread::{Spec, Thread};
fn canonical(family: ThreadFamily, internal: bool, left: bool) -> model::ResolvedThread {
    let params = if family == ThreadFamily::Round {
        let mut p = model::round_thread_presets()[1].parameters();
        p.length = length(2.37 * 3.175);
        p
    } else {
        let row = model::whitworth_thread_presets()
            .into_iter()
            .find(|p| p.family == family)
            .unwrap();
        ThreadParameters {
            family,
            major_diameter: length(row.diameter_mm),
            pitch: length(row.pitch_mm),
            length: length(2.37 * row.pitch_mm),
            fit: if family == ThreadFamily::Bspp {
                if internal {
                    ThreadFit::BsppInternal
                } else {
                    ThreadFit::BsppClassBExternal
                }
            } else if internal {
                ThreadFit::WhitworthMediumInternal
            } else {
                ThreadFit::WhitworthMediumExternal
            },
            ..Default::default()
        }
    };
    let params = ThreadParameters {
        hand: if left {
            ThreadHand::Left
        } else {
            ThreadHand::Right
        },
        ..params
    };
    params
        .resolve(
            Session::new(Document::new("Canonical rounded profile").unwrap())
                .unwrap()
                .parameters(),
            if internal {
                ThreadMode::CutInternal
            } else {
                ThreadMode::NewExternal
            },
        )
        .unwrap()
}
fn descriptor(resolved: &model::ResolvedThread) -> Thread {
    let (crest, root) = if let Some(round) = &resolved.round {
        (round.cutter_crest_radius_mm, round.cutter_root_radius_mm)
    } else {
        let r = resolved.profile.rounding_radius_mm.unwrap();
        (r, r)
    };
    Thread::new(Spec {
        major_radius_mm: resolved.major_radius_mm,
        pitch_mm: resolved.pitch_mm,
        length_mm: resolved.length_mm,
        depth_mm: resolved.profile.depth_mm,
        crest_radius_mm: crest,
        root_radius_mm: root,
        flank_angle_deg: resolved.profile.lower_flank_deg,
        starts: resolved.starts,
        left_hand: resolved.hand == ThreadHand::Left,
    })
    .unwrap()
}
#[test]
fn canonical_whitworth_pipe_and_din405_resolve_without_flattening_or_rewriting() {
    for family in [
        ThreadFamily::Bsw,
        ThreadFamily::Bsf,
        ThreadFamily::Bspp,
        ThreadFamily::Round,
    ] {
        for internal in [false, true] {
            for left in [false, true] {
                let resolved = canonical(family, internal, left);
                let t = descriptor(&resolved);
                let s = t.spec();
                assert_eq!(s.major_radius_mm, resolved.major_radius_mm);
                assert_eq!(s.depth_mm, resolved.profile.depth_mm);
                assert_eq!(s.length_mm, resolved.length_mm);
                assert_eq!(s.left_hand, left);
                assert!(t.volume_mm3() > 0.);
                if family == ThreadFamily::Round && internal {
                    assert_ne!(s.crest_radius_mm, s.root_radius_mm);
                    assert_eq!(
                        s.crest_radius_mm,
                        resolved
                            .round
                            .as_ref()
                            .unwrap()
                            .actual_cutter_crest_circle
                            .radius_mm
                    );
                    assert_eq!(
                        s.root_radius_mm,
                        resolved
                            .round
                            .as_ref()
                            .unwrap()
                            .actual_cutter_root_circle
                            .radius_mm
                    );
                }
            }
        }
    }
}
