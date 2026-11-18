//! A root-rounded metric thread built for one owner is the same validated
//! geometry for every other owner; only the topology names are rebound.
//! The former OCCT owner-thread factory and its bounded root retry have no
//! Rust counterpart, so its work counters are not part of this proof.
use multicad_core::{
    commands::Session,
    document::{Document, EntityId},
    model::{ThreadFit, ThreadHand, ThreadMode, ThreadParameters, length},
    references::Origin,
};
use multicad_kernel::{GeometryPolicy, HelicalThread, PlaneFrame, Solid};

#[test]
fn symmetric_root_retry_retains_only_the_validated_result_and_rebinds_cache_owners() {
    let session = Session::new(Document::new("Root retry cache proof").unwrap()).unwrap();
    let resolved = ThreadParameters {
        major_diameter: length(3.),
        pitch: length(0.5),
        length: length(4.),
        fit: ThreadFit::Metric6gExternal,
        hand: ThreadHand::Left,
        ..Default::default()
    }
    .resolve(session.parameters(), ThreadMode::NewExternal)
    .unwrap();
    let spec = HelicalThread {
        major_radius_mm: resolved.major_radius_mm,
        pitch_mm: resolved.pitch_mm,
        length_mm: resolved.length_mm,
        depth_mm: resolved.profile.depth_mm,
        crest_width_mm: resolved.profile.crest_width_mm,
        lower_flank_deg: resolved.profile.lower_flank_deg,
        upper_flank_deg: resolved.profile.upper_flank_deg,
        radial_taper: resolved.profile.radial_taper,
        starts: resolved.starts,
        left_hand: true,
    };
    let root_radius_mm = resolved
        .metric_fit
        .as_ref()
        .unwrap()
        .external_root_radius_mm
        .unwrap();
    let build = |owner| {
        Solid::rust_helical_thread(
            spec,
            Some(root_radius_mm),
            owner,
            PlaneFrame::default(),
            GeometryPolicy::default(),
        )
        .unwrap()
    };
    let a = EntityId::new();
    let b = EntityId::new();
    let first = build(a);
    assert!(first.uses_rust_geometry());
    assert_eq!(first.solid_count().unwrap(), 1);
    let second = build(b);
    assert_eq!(second.solid_count().unwrap(), 1);
    assert_eq!(first.volume_mm3().unwrap(), second.volume_mm3().unwrap());
    assert_eq!(first.bounds().unwrap(), second.bounds().unwrap());
    let mut expected = first.topology_snapshot(a, 0).unwrap();
    let actual = second.topology_snapshot(b, 0).unwrap();
    expected.owner = b;
    let renamed = |origin: &Origin| {
        assert_eq!(origin.source, a);
        Origin {
            source: b,
            role: origin.role.clone(),
        }
    };
    for element in &mut expected.elements {
        element.origins = element.origins.iter().map(renamed).collect();
        element.boundary = element.boundary.iter().map(renamed).collect();
        for frame in &mut element.frames {
            frame.origin = renamed(&frame.origin);
        }
    }
    assert_eq!(expected, actual);
}
