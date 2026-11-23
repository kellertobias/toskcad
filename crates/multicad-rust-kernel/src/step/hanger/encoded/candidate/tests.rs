use super::*;
use crate::authored_hanger_bolt::{BoundaryOptions, Hanger};
#[test]
fn actual_closed_hanger_candidate_emits_all_trimmed_literal_geometry() {
    for (left, nb) in [(false, false), (true, true)] {
        let source =
            Hanger::new(crate::authored_hanger_bolt::tests::fixture(left, nb), None).unwrap();
        let boundary = source.boundary(BoundaryOptions::default(), None).unwrap();
        let field = EncodedField::prepare(&boundary, CurveOptions::default(), None).unwrap();
        let domains = FiniteDomains::prepare(&field, DomainOptions::default(), None).unwrap();
        let closed = ClosedField::prepare(&domains, DomainOptions::default(), None).unwrap();
        let surfaces = CompactSurfaces::prepare(&field, DomainOptions::default(), None).unwrap();
        let compact =
            CompactBoundary::prepare(&closed, &surfaces, DomainOptions::default(), None).unwrap();
        let vertices = closed
            .literal_vertices(CurveOptions::default(), None)
            .unwrap();
        let events = vertices.cap_events(DomainOptions::default(), None).unwrap();
        let topology = compact
            .trim_caps(&events, DomainOptions::default(), None)
            .unwrap();
        let mut affine = Vec::new();
        let mut roots = Vec::new();
        for segment in 0..events.segments().len() {
            if field.meridians[field.cells[events.segments()[segment].source_cell()].meridian]
                .degree()
                == 2
            {
                roots.push(
                    events
                        .bound_root_cap(segment, CurveOptions::default(), None)
                        .unwrap()
                        .unwrap(),
                );
            } else {
                affine.push(
                    events
                        .bound_affine_cap(segment, CurveOptions::default(), None)
                        .unwrap()
                        .unwrap(),
                );
            }
        }
        let caps = affine
            .iter()
            .map(|c| c as &dyn CapLiteral)
            .chain(roots.iter().map(|c| c as &dyn CapLiteral))
            .collect::<Vec<_>>();
        let helices = (0..compact.edges().len())
            .filter_map(|e| {
                compact
                    .helix(e, &vertices, DomainOptions::default(), None)
                    .unwrap()
            })
            .collect::<Vec<_>>();
        let plains = (0..compact.faces().len())
            .filter_map(|f| {
                compact
                    .plain_wall(f, DomainOptions::default(), None)
                    .unwrap()
            })
            .collect::<Vec<_>>();
        let hr = helices.iter().collect::<Vec<_>>();
        let pr = plains.iter().collect::<Vec<_>>();
        let started = std::time::Instant::now();
        let file =
            encode_candidate(&topology, &caps, &hr, &pr, FileOptions::default(), None).unwrap();
        let used_vertices = topology
            .edges()
            .iter()
            .flat_map(|e| e.vertices)
            .collect::<std::collections::BTreeSet<_>>();
        assert_eq!(
            file.as_str().matches("VERTEX_POINT('").count(),
            used_vertices.len()
        );
        assert!(file.peak_accounted_bytes() <= FileOptions::default().maximum_bytes);
        assert!(file.encoding_work_used() <= FileOptions::default().maximum_encoding_work);
        assert!(file.maximum_stage_work_used() <= FileOptions::default().maximum_stage_work);
        assert!(file.encoding_stages() > topology.faces().len() + topology.edges().len());
        let mut duplicate_caps = caps.clone();
        duplicate_caps[0] = duplicate_caps[1];
        assert!(matches!(
            encode_candidate(
                &topology,
                &duplicate_caps,
                &hr,
                &pr,
                FileOptions::default(),
                None
            ),
            Err(Error::Invalid)
        ));
        assert!(matches!(
            encode_candidate(&topology, &caps, &[], &pr, FileOptions::default(), None),
            Err(Error::Invalid)
        ));
        assert!(matches!(
            encode_candidate(&topology, &caps, &hr, &[], FileOptions::default(), None),
            Err(Error::Invalid)
        ));
        let other_vertices = closed
            .literal_vertices(CurveOptions::default(), None)
            .unwrap();
        let other_events = other_vertices
            .cap_events(DomainOptions::default(), None)
            .unwrap();
        let other_cap = other_events
            .bound_affine_cap(caps[0].segment(), CurveOptions::default(), None)
            .unwrap()
            .unwrap();
        let mut wrong_caps = caps.clone();
        wrong_caps[0] = &other_cap;
        assert!(matches!(
            encode_candidate(
                &topology,
                &wrong_caps,
                &hr,
                &pr,
                FileOptions::default(),
                None
            ),
            Err(Error::Invalid)
        ));
        let mut duplicate_helices = hr.clone();
        duplicate_helices.push(hr[0]);
        assert!(matches!(
            encode_candidate(
                &topology,
                &caps,
                &duplicate_helices,
                &pr,
                FileOptions::default(),
                None
            ),
            Err(Error::Invalid)
        ));
        let mut duplicate_plains = pr.clone();
        duplicate_plains.push(pr[0]);
        assert!(matches!(
            encode_candidate(
                &topology,
                &caps,
                &hr,
                &duplicate_plains,
                FileOptions::default(),
                None
            ),
            Err(Error::Invalid)
        ));
        let wrong_helices = (0..compact.edges().len())
            .filter_map(|e| {
                compact
                    .helix(e, &other_vertices, DomainOptions::default(), None)
                    .unwrap()
            })
            .collect::<Vec<_>>();
        let wrong_refs = wrong_helices.iter().collect::<Vec<_>>();
        assert!(matches!(
            encode_candidate(
                &topology,
                &caps,
                &wrong_refs,
                &pr,
                FileOptions::default(),
                None
            ),
            Err(Error::Invalid)
        ));
        for invalid in [
            FileOptions {
                maximum_encoding_work: 0,
                ..FileOptions::default()
            },
            FileOptions {
                maximum_encoding_work: 256_000_001,
                ..FileOptions::default()
            },
            FileOptions {
                maximum_stage_work: 0,
                ..FileOptions::default()
            },
            FileOptions {
                maximum_stage_work: 16_000_001,
                ..FileOptions::default()
            },
            FileOptions {
                maximum_bytes: 0,
                ..FileOptions::default()
            },
            FileOptions {
                maximum_bytes: 64 * 1024 * 1024 + 1,
                ..FileOptions::default()
            },
        ] {
            assert!(matches!(
                encode_candidate(&topology, &caps, &hr, &pr, invalid, None),
                Err(Error::Invalid)
            ));
        }
        for limited in [
            FileOptions {
                maximum_encoding_work: 1,
                ..FileOptions::default()
            },
            FileOptions {
                maximum_stage_work: 1,
                ..FileOptions::default()
            },
            FileOptions {
                maximum_bytes: 1,
                ..FileOptions::default()
            },
        ] {
            assert!(matches!(
                encode_candidate(&topology, &caps, &hr, &pr, limited, None),
                Err(Error::Limit)
            ));
        }
        assert!(file.as_str().starts_with("ISO-10303-21;"));
        assert!(file.as_str().ends_with("END-ISO-10303-21;\n"));
        assert_eq!(
            file.as_str().matches("ADVANCED_FACE('").count(),
            topology.faces().len()
        );
        assert_eq!(
            file.as_str().matches("EDGE_CURVE('").count(),
            topology.edges().len()
        );
        assert_eq!(
            file.as_str().matches("PCURVE('").count(),
            2 * topology.edges().len()
        );
        assert_eq!(file.as_str().matches("MANIFOLD_SOLID_BREP('").count(), 1);
        assert!(file.retained_bytes() >= file.as_str().len());
        assert!(file.aggregate_known_work_used() > file.encoding_work_used());
        eprintln!(
            "candidate left={left} bytes={} encodingwork={} aggregateknownwork={} maxstage={} stages={} peak={} encoding_seconds={}",
            file.as_str().len(),
            file.encoding_work_used(),
            file.aggregate_known_work_used(),
            file.maximum_stage_work_used(),
            file.encoding_stages(),
            file.peak_accounted_bytes(),
            started.elapsed().as_secs_f64()
        );
        std::fs::write(
            std::env::temp_dir().join(format!("hanger-actual-trimmed-{left}.step")),
            file.as_str(),
        )
        .unwrap();
        assert!(matches!(
            encode_candidate(
                &topology,
                &caps,
                &hr,
                &pr,
                FileOptions {
                    maximum_encoding_work: file.encoding_work_used() - 1,
                    ..FileOptions::default()
                },
                None
            ),
            Err(Error::Limit)
        ));
        assert!(matches!(
            encode_candidate(
                &topology,
                &caps,
                &hr,
                &pr,
                FileOptions::default(),
                Some(&|| true)
            ),
            Err(Error::Cancelled)
        ));
        assert!(matches!(
            encode_candidate(&topology, &[], &hr, &pr, FileOptions::default(), None),
            Err(Error::Invalid)
        ));
    }
}
