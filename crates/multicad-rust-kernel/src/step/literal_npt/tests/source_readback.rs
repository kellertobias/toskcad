//! Public source admission: these tests never substitute Selected/retain for
//! the same discovery/fallback entry point consumed by the application.
use super::*;
use crate::step_read::{ReadError, SourceImport, read_source_compound};

#[test]
fn actual_public_source_import_retains_both_npt_sizes_and_original_mapping() {
    for (quarter, posed) in [(false, false), (true, false), (false, true)] {
        let original = text(quarter);
        let input = if posed { mapped(&original) } else { original };
        let mut members = read_source_compound(&input, 1e-7).unwrap();
        assert_eq!(members.len(), 1);
        let SourceImport::LiteralNpt(source) = members.pop().unwrap() else {
            panic!("actual original C11 NPT must retain its source");
        };
        assert!(source.world_volume_mm3()[0] > 0.);
        assert_eq!(source.member().mapping_link().is_some(), posed);
        let replay = verify_exact_replay(&source, Default::default(), None).unwrap();
        assert_eq!(replay.receipt().position_error_mm(), 0.);
        let malformed = input.replacen(
            "B_SPLINE_SURFACE_WITH_KNOTS('',11,11",
            "B_SPLINE_SURFACE_WITH_KNOTS('',11,10",
            1,
        );
        assert_ne!(malformed, input);
        assert!(read_source_compound(&malformed, 1e-7).is_err());
    }
}

fn torx() -> crate::profiled_screw::Screw {
    use crate::{
        head_profile::Head,
        helical_thread::{Spec, Thread},
        hexalobular_drive::{Cone, Drive, TangentSpec},
        profiled_screw::Screw,
        screw_recess::Recess,
    };
    let drive = Drive::from_authored_tangent(
        TangentSpec {
            major_diameter_mm: 2.8,
            minor_diameter_mm: 2.05,
            outer_lobe_radius_mm: 0.28,
        },
        Cone {
            tip_depth_mm: 1.1 + 1.4 / 3f64.sqrt(),
            depth_per_radius: 1. / 3f64.sqrt(),
        },
    )
    .unwrap();
    let recess = Recess::conical_hexalobular(drive, 1e-7, 100_000).unwrap();
    let thread = Thread::new(Spec {
        major_radius_mm: 1.5,
        pitch_mm: 0.5,
        length_mm: 3.,
        depth_mm: 0.2,
        crest_width_mm: 0.0625,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: false,
        root_radius_mm: None,
    })
    .unwrap();
    Screw::new_with_recess(
        thread,
        1.5,
        Head::Crowned {
            radius_mm: 2.8,
            height_mm: 2.4,
            crown_radius_mm: 5.,
            top_land_radius_mm: 2.05,
            neck_radius_mm: 0.125,
        },
        0.,
        0.1,
        recess,
    )
    .unwrap()
}

#[test]
fn actual_public_source_import_preserves_torx_and_npt_in_both_orders() {
    let screw = torx();
    let source = source(false, false);
    let verified = verify_exact_replay(&source, Default::default(), None).unwrap();
    for torx_first in [true, false] {
        let mut writer = Writer::new();
        for add_torx in [torx_first, !torx_first] {
            if add_torx {
                writer
                    .append_profiled_screw(&screw, [-20., 0., 0.], [1., 0., 0.], [0., 1., 0.], 1e-7)
                    .unwrap();
            } else {
                writer
                    .append_verified_literal_npt(&verified, Default::default(), None)
                    .unwrap();
            }
        }
        let text = writer.finish().unwrap();
        let members = read_source_compound(&text, 1e-7).unwrap();
        assert_eq!(members.len(), 2);
        assert_eq!(
            members
                .iter()
                .filter(|m| matches!(m, SourceImport::LiteralNpt(_)))
                .count(),
            1
        );
        let analytic = members
            .iter()
            .find_map(|m| match m {
                SourceImport::Analytic(a) => Some(a),
                _ => None,
            })
            .unwrap();
        assert!(analytic.0.volume_mm3 > 0.);
        assert!(matches!(
            members[usize::from(torx_first)],
            SourceImport::LiteralNpt(_)
        ));
    }
}

#[test]
fn actual_public_source_import_refuses_orphan_topology_and_false_material_orientation() {
    let input = text(false);
    let document = CompoundDocument::parse(&input, Default::default(), None).unwrap();
    let selected = Selected::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap();
    let point = selected.parsed().points()[0].entity();
    let orphan = input.replacen(
        "ENDSEC;\nEND-ISO-10303-21;",
        &format!("#999999=VERTEX_POINT('',#{point});\nENDSEC;\nEND-ISO-10303-21;"),
        1,
    );
    assert!(matches!(
        read_source_compound(&orphan, 1e-7),
        Err(ReadError::Invalid)
    ));
    let face = input
        .lines()
        .find(|line| line.contains("=ADVANCED_FACE("))
        .unwrap();
    let reversed = if face.contains(",.T.)") {
        face.replace(",.T.)", ",.F.)")
    } else {
        face.replace(",.F.)", ",.T.)")
    };
    assert_ne!(face, reversed);
    let invalid = input.replacen(face, &reversed, 1);
    assert!(read_source_compound(&invalid, 1e-7).is_err());
    assert!(read_source_compound(&input, 0.).is_err());
    assert!(read_source_compound("not STEP", 1e-7).is_err());
}

#[test]
fn actual_public_source_import_keeps_one_packet_qualification_work_ceiling() {
    let source = source(false, false);
    let verified = verify_exact_replay(&source, Default::default(), None).unwrap();
    let mut writer = Writer::new();
    // Every member is independently valid literal geometry. Their combined
    // original-source qualification exceeds the public fallback's512M ledger;
    // retaining a new member must not reset that packet-wide allowance.
    for _ in 0..13 {
        writer
            .append_verified_literal_npt(&verified, Default::default(), None)
            .unwrap();
    }
    let text = writer.finish().unwrap();
    assert!(matches!(
        read_source_compound(&text, 1e-7),
        Err(ReadError::Limit)
    ));
    assert!(matches!(
        read_source_compound(verified.standalone_step(), 1e-7)
            .unwrap()
            .as_slice(),
        [SourceImport::LiteralNpt(_)]
    ));
}
