//! Semantic definition of the OWNED repaired material. No original record ID,
//! file text/header, volume-only proxy or declared authoring recipe is a key.
use super::*;
pub(super) fn signature(source: &step_read::PhillipsOwnedPeriodicMaterial) -> ContentHash {
    let mut bytes = b"multicad.periodic-phillips-definition/1".to_vec();
    fn numbers(bytes: &mut Vec<u8>, values: impl IntoIterator<Item = f64>) {
        for value in values {
            bytes.extend(if value == 0. { 0. } else { value }.to_le_bytes());
        }
    }
    // These immutable parameters define the qualified corrected meridian's
    // evaluator, not the unmodified literal STEP material. Its construction
    // includes shared endpoints and the circle correction for each arc.
    for segment in source.head().view().profile_proposal().segments() {
        let role = segment.role.source_role();
        bytes.extend((role.len() as u64).to_le_bytes());
        bytes.extend(role.as_bytes());
        numbers(&mut bytes, segment.start.into_iter().chain(segment.end));
        match segment.curve {
            multicad_rust_kernel::head_profile::Curve::Line => bytes.push(0),
            multicad_rust_kernel::head_profile::Curve::Arc {
                center,
                radius_mm,
                sweep,
            } => {
                bytes.push(1);
                numbers(&mut bytes, center.into_iter().chain([radius_mm, sweep]));
            }
        }
    }
    let thread = source.thread();
    numbers(
        &mut bytes,
        thread.phase_node_expressions_mm().into_iter().flatten(),
    );
    numbers(
        &mut bytes,
        thread
            .radius_nodes_mm()
            .into_iter()
            .chain([thread.pitch_mm(), thread.length_mm()]),
    );
    bytes.push(u8::from(thread.left_hand()));
    numbers(
        &mut bytes,
        source.cap_expressions_mm().into_iter().flatten(),
    );
    numbers(&mut bytes, [source.head_z_shift_mm()]);
    numbers(&mut bytes, source.cavity_mouth_expression_mm());
    let chart = source.chart();
    numbers(&mut bytes, chart.columns.into_iter().flatten());
    numbers(&mut bytes, [chart.outer_top_z]);
    // The independently closed cavity proves its finite authoring-domain caps
    // inactive. Its complete active implicit field is defined by this spec;
    // changing an inactive reach bound does not change physical material.
    let s = source.cavity().spec();
    numbers(
        &mut bytes,
        [
            s.mouth_diameter_mm,
            s.mouth_radius_mm,
            s.b_width_mm,
            s.g_diameter_mm,
            s.tilt_angle_deg,
            s.outer_flank_angle_deg,
            s.floor_angle_deg,
            s.outer_corner_angle_deg,
        ],
    );
    match s.inner {
        multicad_rust_kernel::phillips_drive::Inner::ShallowV {
            angle_deg,
            e_width_mm,
        } => {
            bytes.push(0);
            numbers(&mut bytes, [angle_deg, e_width_mm]);
        }
        multicad_rust_kernel::phillips_drive::Inner::Rounded { radius_mm } => {
            bytes.push(1);
            numbers(&mut bytes, [radius_mm]);
        }
    }
    ContentHash::of(&bytes)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeSet;
    fn file(reverse: bool) -> String {
        let a = crate::geometry::repaired_phillips_body::step::tests::part(
            false,
            PlaneFrame::default(),
        );
        let b = crate::geometry::repaired_phillips_body::step::tests::part(
            true,
            PlaneFrame {
                origin: [5., 10., 18.],
                u: [0., 1., 0.],
                v: [1., 0., 0.],
            },
        );
        let mut w = multicad_rust_kernel::step::Writer::new();
        if reverse {
            b.append_periodic_phillips_step(&mut w, None)
                .unwrap()
                .unwrap();
        }
        w.append_sphere(
            multicad_rust_kernel::sphere::Sphere::new(1.).unwrap(),
            [20., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
        )
        .unwrap();
        a.append_periodic_phillips_step(&mut w, None)
            .unwrap()
            .unwrap();
        if !reverse {
            b.append_periodic_phillips_step(&mut w, None)
                .unwrap()
                .unwrap();
        }
        w.finish().unwrap()
    }
    #[test]
    fn public_mixed_import_retains_real_periodic_sources_and_semantic_member_identity() {
        let feature = EntityId::from_uuid_bits(0xcad352a);
        let first = file(false);
        let second = file(true);
        let imported = Solid::import_rust_step(&first, feature, Default::default()).unwrap();
        let reordered = Solid::import_rust_step(&second, feature, Default::default()).unwrap();
        let a = imported.topology_snapshot(feature, 0).unwrap();
        let b = reordered.topology_snapshot(feature, 0).unwrap();
        let origins = |snapshot: &TopologySnapshot| {
            snapshot
                .elements
                .iter()
                .flat_map(|e| e.origins.iter().cloned())
                .collect::<BTreeSet<_>>()
        };
        assert_eq!(origins(&a), origins(&b));
        assert!(imported.volume_mm3().unwrap() > 90.);
        assert!((imported.volume_mm3().unwrap() - reordered.volume_mm3().unwrap()).abs() < 1e-7);
        let decoded = step_read::read_source_compound(&first, 1e-7).unwrap();
        assert_eq!(decoded.len(), 3);
        assert_eq!(
            decoded
                .iter()
                .filter(|v| matches!(v, step_read::SourceImport::PeriodicPhillips(_)))
                .count(),
            2
        );
        assert_eq!(
            decoded
                .iter()
                .filter(|v| matches!(v, step_read::SourceImport::Analytic(_)))
                .count(),
            1
        );
        let insert = first.rfind("ENDSEC;").unwrap();
        let mut orphan = first.clone();
        orphan.insert_str(insert, "#999999=VERTEX_POINT('',#1);\n");
        assert_eq!(
            Solid::import_rust_step(&orphan, feature, Default::default())
                .err()
                .unwrap()
                .kind,
            GeometryErrorKind::InvalidResult
        );
    }
}
