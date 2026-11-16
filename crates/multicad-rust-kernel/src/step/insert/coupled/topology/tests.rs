use super::*;
#[test]
fn actual_nextafter_and_standard_coupled_rows_form_closed_oriented_annular_incidence() {
    for (left, length) in [
        (false, 4.),
        (false, 4_f64.next_up()),
        (true, 4_f64.next_up()),
        (false, 4_f64.next_down()),
        (true, 8.),
    ] {
        let source = crate::step_read::owned_insert::tests::graph(left, length);
        let plan = CoupledPlan::prepare(&source, CurveOptions::default(), None).unwrap();
        let top = CoupledTopology::prepare(&plan, 16_000_000, 64 * 1024 * 1024, None).unwrap();
        assert!(std::ptr::eq(top.plan(), &plan));
        assert_eq!(top.patches().len(), plan.domain().slabs().len());
        assert_eq!(top.faces().len(), top.patches().len() + 2);
        assert_eq!(top.vertices().len(), top.used_vertices());
        for vertex in top.vertices() {
            let [row, point] = vertex.coefficient_reference();
            assert_eq!(
                vertex.owner(),
                plan.rows()[row].endpoint_owners()[usize::from(point != 0)]
            );
        }
        let mut uses = vec![[0; 2]; top.edges().len()];
        for face in top.faces() {
            for wire in &face.loops {
                for coedge in wire {
                    uses[coedge.edge][usize::from(coedge.forward)] += 1;
                }
            }
        }
        assert!(uses.iter().all(|u| *u == [1, 1]));
        // Two annuli each contribute0 to Euler characteristic; tensor disks1.
        assert_eq!(top.used_vertices() + top.patches().len(), top.edges().len());
        eprintln!(
            "coupled incidence {left}/{length:?}: V{} E{} tensor{} {}work {}bytes",
            top.used_vertices(),
            top.edges().len(),
            top.patches().len(),
            top.work_used(),
            top.retained_bytes()
        );
    }
}
