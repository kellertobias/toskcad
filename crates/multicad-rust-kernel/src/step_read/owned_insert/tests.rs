use super::*;
fn file(left: bool) -> String {
    use crate::{
        authored_flat_thread::{AuthoredThread, SourceOptions},
        helical_thread::Spec,
        step::insert::{
            CurveOptions, DomainOptions, DomainPlan, SurfacePlan, TopologyPlan,
            write_local_candidate,
        },
        wood_insert::{BoundaryOptions, Insert},
    };
    let spec = Spec {
        major_radius_mm: 5.,
        pitch_mm: 2.,
        length_mm: 4.,
        depth_mm: 1.,
        crest_width_mm: 0.4,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: left,
        root_radius_mm: None,
    };
    let options = SourceOptions {
        position_limit_mm: 1e-7,
        normal_angle_limit_rad: 1e-7,
    };
    let outer = AuthoredThread::new(spec, options).unwrap();
    let inner = AuthoredThread::new(
        Spec {
            major_radius_mm: 1.5,
            pitch_mm: 0.5,
            depth_mm: (3. - (2.459 + 2.599) / 2.) / 2.,
            crest_width_mm: 0.25 - (3. - (2.675 + 2.775) / 2.) / 3_f64.sqrt(),
            ..spec
        },
        options,
    )
    .unwrap();
    let body = Insert::new(outer, inner, 0.5, None)
        .unwrap()
        .boundary_with_interrupt(BoundaryOptions::default(), None)
        .unwrap();
    let d = DomainPlan::prepare(&body, DomainOptions::default(), None).unwrap();
    let s = SurfacePlan::prepare(&d, CurveOptions::default(), None).unwrap();
    let t = TopologyPlan::prepare(&s, 16_000_000, 64 * 1024 * 1024, None).unwrap();
    write_local_candidate(&t, 16_000_000, 64 * 1024 * 1024, None).unwrap()
}
#[test]
fn actual_file_retains_original_scalar_carriers_and_complete_represented_graph() {
    for left in [false, true] {
        let text = file(left);
        let p = parse(&text, Options::default(), None).unwrap();
        assert_eq!(p.literal_input(), text);
        assert_eq!(p.planes().len(), 2);
        assert_eq!(p.faces().len(), p.tensors().len() + 2);
        assert!(p.work_used() <= 16_000_000);
        assert!(p.retained_budget_bytes() <= 64 * 1024 * 1024);
        assert!(p.literal_entity_count() > p.points().len());
        for point in p.points() {
            assert_eq!(p.point(point.entity()).unwrap().entity(), point.entity());
            for v in point.coordinates() {
                assert_eq!(v.spelling().parse::<f64>().unwrap(), v.value());
                let [a, b] = v.bounds();
                assert!(a <= v.value() && b >= v.value());
            }
        }
        eprintln!(
            "ACTUAL_PARSE left{left} faces{} edges{} points{} work{} bytes{}",
            p.faces().len(),
            p.edges().len(),
            p.points().len(),
            p.work_used(),
            p.retained_budget_bytes()
        );
    }
}

#[test]
fn actual_local_file_has_literal_dual_support_and_full_uv_domain_binding() {
    for left in [false, true] {
        let text = file(left);
        let parsed = parse(&text, Options::default(), None).unwrap();
        let proof = parsed.support_bindings(Options::default(), None).unwrap();
        assert!(std::ptr::eq(proof.source(), &parsed));
        assert_eq!(proof.tensors().len(), parsed.tensors().len());
        assert_eq!(proof.edge_face_indices().len(), parsed.edges().len());
        eprintln!(
            "ACTUAL_BIND left{left} work{} bytes{}",
            proof.work_used(),
            proof.retained_budget_bytes()
        );
    }
}
