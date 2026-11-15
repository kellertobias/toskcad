use super::*;
fn stock(quarter: bool) -> crate::axial_pipe_thread::authored::StockClipPlan {
    use crate::axial_pipe_thread::authored::{AuthoredSection, AuthoredSpec, AuthoredThread};
    let (tpi, e1, l1, l2, od, max) = if quarter {
        (18., 0.49163, 0.2278, 0.4018, 0.540, 0.088)
    } else {
        (27., 0.37360, 0.1615, 0.2639, 0.405, 0.096)
    };
    let pitch = 25.4 / tpi;
    let k = 1. / 32.;
    let m = 3f64.sqrt();
    let s = 1. / (m - k) + 1. / (m + k);
    let depth = pitch / s - (0.033 + max) * pitch;
    let spec = AuthoredSpec {
        major_radius_mm: (e1 * 25.4 - 2. * k * l1 * 25.4 + depth) / 2.,
        pitch_mm: pitch,
        length_mm: l2 * 25.4,
        radial_taper: k,
        section: AuthoredSection::NptFlat {
            depth_mm: depth,
            flat_width_mm: (pitch - depth * s) / 2.,
        },
        crest_radius_limit_mm: Some(od * 25.4 / 2.),
    };
    AuthoredThread::new(spec, 1e-7, 1e-6)
        .unwrap()
        .stock_clip_plan(Default::default(), None)
        .unwrap()
}
fn text(quarter: bool) -> String {
    use crate::step::npt_stock as step;
    let source = stock(quarter);
    let curves = step::CurvePlan::new(&source, Default::default(), None).unwrap();
    let surfaces = step::SurfacePlan::new(curves, Default::default(), None).unwrap();
    step::write_source_candidate(&surfaces, 16_000_000, 64 * 1024 * 1024, None)
        .unwrap()
        .into_string()
}
fn mapped(text: &str) -> String {
    let doc =
        crate::step_read::owned_insert::CompoundDocument::parse(text, Default::default(), None)
            .unwrap();
    let solid = doc.represented_members()[0].item_entity();
    let old = text
        .lines()
        .find(|s| s.contains("=ADVANCED_BREP_SHAPE_REPRESENTATION("))
        .unwrap();
    let context = old.rsplit('#').next().unwrap().split(')').next().unwrap();
    let changed = old.replacen(&format!("#{solid}"), "#900020", 1);
    let extra = format!(
        "#900010=CARTESIAN_POINT('',(0.,0.,0.));\n#900011=AXIS2_PLACEMENT_3D('',#900010,$,$);\n#900012=SHAPE_REPRESENTATION('',(#{solid},#900011),#{context});\n#900013=REPRESENTATION_MAP(#900011,#900012);\n#900014=CARTESIAN_POINT('',(11.,-5.,3.));\n#900015=DIRECTION('',(0.,2.,0.));\n#900016=DIRECTION('',(-3.,0.,0.));\n#900017=DIRECTION('',(0.,-4.,3.));\n#900018=(CARTESIAN_TRANSFORMATION_OPERATOR(#900015,#900016,#900014,1.)CARTESIAN_TRANSFORMATION_OPERATOR_3D(#900017)FUNCTIONALLY_DEFINED_TRANSFORMATION('','')GEOMETRIC_REPRESENTATION_ITEM()REPRESENTATION_ITEM(''));\n#900020=MAPPED_ITEM('',#900013,#900018);\n"
    );
    text.replacen(old, &changed, 1).replacen(
        "ENDSEC;\nEND-ISO-10303-21;",
        &format!("{extra}ENDSEC;\nEND-ISO-10303-21;"),
        1,
    )
}
fn source(quarter: bool, placed: bool) -> OwnedMaterial {
    let original = text(quarter);
    let text = if placed { mapped(&original) } else { original };
    let document = CompoundDocument::parse(&text, Default::default(), None).unwrap();
    Selected::parse(
        &document,
        document.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap()
    .retain(Default::default(), None)
    .unwrap()
}
#[test]
fn actual_retained_literal_npt_replays_exact_original_geometry_and_mapping() {
    for (quarter, mapped) in [(false, false), (true, false), (false, true)] {
        let source = source(quarter, mapped);
        eprintln!(
            "REPLAY input retained={} docbytes={} records={} work={}",
            source.retained_budget_bytes(),
            source.original_document_budget_bytes(),
            source.source_entity_count(),
            source.work_used()
        );
        let started = std::time::Instant::now();
        let result = verify_exact_replay(&source, Default::default(), None);
        eprintln!(
            "literal replay quarter={quarter} mapped={mapped} elapsed={:?} result={:?}",
            started.elapsed(),
            result.as_ref().map(|r| r.receipt())
        );
        let verified = result.unwrap();
        assert_eq!(verified.receipt().position_error_mm(), 0.);
        assert_eq!(verified.receipt().original_derivative_error_mm(), 0.);
        assert_eq!(verified.receipt().normal_chord_error(), 0.);
        let mut writer = Writer::new();
        writer
            .append_box(
                [Rational::new(2, 1).unwrap(); 3],
                [-20., 0., 0.],
                [1., 0., 0.],
                [0., 1., 0.],
            )
            .unwrap();
        let receipt = writer
            .append_verified_literal_npt(&verified, Default::default(), None)
            .unwrap();
        assert!(receipt.aggregate_work_used() <= WORK && receipt.accounted_peak_bytes() <= BYTES);
        let text = writer.finish().unwrap();
        let doc = CompoundDocument::parse(&text, Default::default(), None).unwrap();
        let selected = Selected::parse(
            &doc,
            doc.represented_members()[1].item_entity(),
            Default::default(),
            None,
        )
        .unwrap();
        let decoded = selected.retain(Default::default(), None).unwrap();
        assert_eq!(
            source.original_mapping_identity(),
            decoded.original_mapping_identity()
        );
        assert_eq!(source.world_volume_mm3(), decoded.world_volume_mm3());
    }
}
#[test]
fn actual_retained_replay_refuses_changed_exact_decimal_despite_same_float_cache() {
    let source = source(false, false);
    let verified = verify_exact_replay(&source, Default::default(), None).unwrap();
    let tensor = &source.tensors()[0];
    let point = source
        .point(tensor.point_indices()[11 * tensor.spans()[1] + 2])
        .unwrap();
    let index = source.record_index(point.entity()).unwrap();
    let actual_id = verified.trace.nodes[index].id;
    let old = point.coordinates()[0].spelling();
    let (mantissa, exponent) = old.split_once('E').unwrap();
    let changed_number = format!("{mantissa}00000000000000001E{exponent}");
    let key = format!("#{actual_id}=");
    let line = verified
        .packet
        .lines()
        .find(|l| l.starts_with(&key))
        .unwrap();
    let changed = verified
        .packet
        .replacen(line, &line.replacen(old, &changed_number, 1), 1);
    let document = CompoundDocument::parse(&changed, Default::default(), None).unwrap();
    let item = document.represented_members()[0].item_entity();
    let selected = Selected::parse(&document, item, Default::default(), None).unwrap();
    let decoded = selected
        .parsed()
        .points()
        .iter()
        .find(|p| p.entity() == actual_id)
        .unwrap();
    assert_eq!(
        decoded.coordinates()[0].value(),
        point.coordinates()[0].value()
    );
    selected
        .parsed()
        .projection(Default::default(), None)
        .unwrap()
        .occupied(Default::default(), None)
        .unwrap();
    assert_eq!(
        binding::bind_owned_document(
            &source,
            source.member().item_entity(),
            &document,
            selected.entity_ids(),
            item,
            binding::Options {
                maximum_work: 96_000_000,
                maximum_entities: 1_048_576,
                ..Default::default()
            },
            None
        )
        .err(),
        Some(binding::Error::Different)
    );
}
fn prefix() -> Writer {
    let mut writer = Writer::new();
    writer
        .append_box(
            [Rational::new(2, 1).unwrap(); 3],
            [-20., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
        )
        .unwrap();
    writer
}
#[test]
fn actual_retained_replay_uses_one_new_query_ledger_and_atomic_append_frontiers() {
    let source = source(false, false);
    assert_eq!(
        verify_exact_replay(
            &source,
            Options {
                maximum_work: 1,
                ..Default::default()
            },
            None
        )
        .err()
        .unwrap(),
        Error::Limit
    );
    assert_eq!(
        verify_exact_replay(&source, Default::default(), Some(&|| true))
            .err()
            .unwrap(),
        Error::Cancelled
    );
    let verified = verify_exact_replay(&source, Default::default(), None).unwrap();
    let mut good = prefix();
    let receipt = good
        .append_verified_literal_npt(&verified, Default::default(), None)
        .unwrap();
    assert_eq!(
        receipt.end_to_end_work_used(),
        receipt.source_document_work_used()
            + receipt.source_qualification_work_used()
            + receipt.aggregate_work_used()
    );
    assert!(
        receipt.append_work_used() > 0
            && receipt.binding_work_used() > 0
            && receipt.material_work_used() > 0
            && receipt.encoding_work_used() > 0
            && receipt.document_work_used() > 0
    );
    let expected = prefix().finish().unwrap();
    for options in [
        Options {
            maximum_work: 1,
            ..Default::default()
        },
        Options {
            maximum_work: receipt.aggregate_work_used() - 1,
            ..Default::default()
        },
        Options {
            maximum_bytes: 1,
            ..Default::default()
        },
    ] {
        let mut writer = prefix();
        assert_eq!(
            writer
                .append_verified_literal_npt(&verified, options, None)
                .err(),
            Some(Error::Limit)
        );
        assert_eq!(writer.finish().unwrap(), expected);
    }
    let mut writer = prefix();
    assert_eq!(
        writer
            .append_verified_literal_npt(&verified, Default::default(), Some(&|| true))
            .err(),
        Some(Error::Cancelled)
    );
    assert_eq!(writer.finish().unwrap(), expected);
    assert!(verified.retained_budget_bytes() < BYTES);
    assert_eq!(verified.source().point_count(), source.point_count());
    assert!(verified.standalone_step().starts_with("ISO-10303-21"));
}
#[test]
fn actual_native_chain_step_replay_binds_local_geometry_and_earns_whole_operator_correspondence() {
    for (quarter, mapped) in [(false, false), (true, false), (false, true)] {
        let source = source(quarter, mapped);
        let pose = read::LiteralOccurrencePlacement::new(
            [[0.8, -0.36, 0.48], [0.6, 0.48, -0.64], [0., 0.8, 0.6]],
            [17., -8., 13.],
            Default::default(),
            None,
        )
        .unwrap();
        let pose = pose
            .transformed(
                [[0., 0., 1.], [0., 1., 0.], [-1., 0., 0.]],
                [2., 3., 4.],
                Default::default(),
                None,
            )
            .unwrap();
        let started = std::time::Instant::now();
        let candidate = verify_placed_replay(&source, &pose, Default::default(), None);
        eprintln!(
            "placed literal replay quarter={quarter} original_mapped={mapped} elapsed={:?} result={:?}",
            started.elapsed(),
            candidate.as_ref().map(|c| c.receipt())
        );
        let candidate = candidate.unwrap();
        assert_eq!(candidate.occurrence().unwrap().transforms().len(), 2);
        let receipt = candidate.receipt();
        assert!(
            receipt.position_error_mm() <= 1e-7
                && receipt.original_derivative_error_mm() <= 1e-6
                && receipt.normal_chord_error() <= 1e-6
        );
        assert!(receipt.literal_volume_mm3()[0] > 0.);
        let mut writer = prefix();
        writer
            .append_verified_literal_npt(&candidate, Default::default(), None)
            .unwrap();
        let text = writer.finish().unwrap();
        let document = CompoundDocument::parse(&text, Default::default(), None).unwrap();
        let selected = Selected::parse(
            &document,
            document.represented_members()[1].item_entity(),
            Default::default(),
            None,
        )
        .unwrap();
        let decoded = selected.retain(Default::default(), None).unwrap();
        let local = [[0., 0.], [0., 0.], [2., 2.]];
        let wanted = pose.map_box_mm(source.map_box_mm(local).unwrap()).unwrap();
        let actual = decoded.map_box_mm(local).unwrap();
        assert!(
            wanted
                .into_iter()
                .zip(actual)
                .all(|(a, b)| a[0] <= b[1] + 1e-7 && b[0] <= a[1] + 1e-7)
        );
        assert_eq!(source.member().mapping_link().is_some(), mapped);
    }
}
