use super::*;
#[test]
fn shared_chart_complex_preserves_every_source_boundary_and_collapsed_identity() {
    for fixture in 0..=5 {
        let (s, d) = super::super::super::super::super::tests::fixture(fixture);
        let source = Cutter::new(s, d).unwrap();
        let atlas = source.surface_patches().unwrap();
        let complex = Complex::build(&atlas, &|| false).unwrap();
        assert_eq!(complex.points.len(), 66);
        assert_eq!(complex.edges.len(), 152);
        assert_eq!(complex.sides.len(), 88);
        for (edge, e) in complex.edges.iter().enumerate() {
            assert_eq!(
                complex.edge_point(&atlas, edge, 0.).unwrap().map(bounds),
                complex.points[e.ends[0]].map(bounds)
            );
            assert_eq!(
                complex.edge_point(&atlas, edge, 1.).unwrap().map(bounds),
                complex.points[e.ends[1]].map(bounds)
            );
            for t in [0.25, 0.5, 0.75] {
                assert!(
                    complex
                        .edge_point(&atlas, edge, t)
                        .unwrap()
                        .into_iter()
                        .all(I::valid)
                );
            }
        }
    }
}
#[test]
fn refinement_source_rejects_missing_incidence_unknown_anchors_and_bad_parameters() {
    let (s, d) = super::super::super::super::super::tests::fixture(0);
    let source = Cutter::new(s, d).unwrap();
    let atlas = source.surface_patches().unwrap();
    assert_eq!(Row::from_event(EventRole::FloorApex), Err(Error::Invalid));
    assert_eq!(
        original(
            &atlas,
            VertexKey::Event {
                role: EventRole::FloorApex,
                quarter: 1,
                lower: false
            }
        ),
        Err(Error::Unresolved)
    );
    assert_eq!(
        node_point(&atlas, Node::Original(usize::MAX)).unwrap_err(),
        Error::Invalid
    );
    assert!(matches!(
        Complex::build(&atlas, &|| true),
        Err(Error::Cancelled)
    ));
    let complex = Complex::build(&atlas, &|| false).unwrap();
    for t in [f64::NAN, -0.1, 1.1] {
        assert_eq!(
            complex.edge_point(&atlas, 0, t).unwrap_err(),
            Error::Invalid
        );
    }
    assert_eq!(
        complex.edge_point(&atlas, usize::MAX, 0.5).unwrap_err(),
        Error::Invalid
    );
    let mut broken = atlas.clone();
    broken.patches.pop();
    assert!(matches!(
        Complex::build(&broken, &|| false),
        Err(Error::Unresolved)
    ));
}
